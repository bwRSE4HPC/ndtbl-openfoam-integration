#!/bin/bash

set -euo pipefail

SCRIPT_DIR="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
RESULTS_DIR="$SCRIPT_DIR/benchmark-results"
RESULTS_CSV="$RESULTS_DIR/results.csv"
TABLE_SUFFIX="f"
RESOLUTIONS=()
MPI_COUNTS=()

usage() {
    cat <<'EOF'
Usage:
  run_ndtbl_benchmarks.sh [options]

Options:
  --resolution N     Add one table resolution to sweep. Repeatable.
  --mpi N            Add one MPI process count. Use 0 or 'none' for serial.
                     Repeatable. Default sweep: none, 2, 4, 8.
  --table-suffix S   ndtbl precision suffix passed to tableGenerator.py.
                     Use f for float32 or d for float64. Default: f
  --results-dir DIR  Directory for CSV results. Default: benchmark-results
  --help             Show this help

Examples:
  ./run_ndtbl_benchmarks.sh --resolution 21 --resolution 31
  ./run_ndtbl_benchmarks.sh --resolution 21 --mpi 0 --mpi 2 --mpi 4 --mpi 8
EOF
}

normalize_mpi_count() {
    local value="$1"
    if [ "$value" = "none" ] || [ "$value" = "None" ]; then
        printf "0\n"
    else
        printf "%s\n" "$value"
    fi
}

while [ $# -gt 0 ]; do
    case "$1" in
        --resolution)
            RESOLUTIONS+=("$2")
            shift 2
            ;;
        --mpi)
            MPI_COUNTS+=("$(normalize_mpi_count "$2")")
            shift 2
            ;;
        --table-suffix)
            TABLE_SUFFIX="$2"
            shift 2
            ;;
        --results-dir)
            RESULTS_DIR="$2"
            RESULTS_CSV="$RESULTS_DIR/results.csv"
            shift 2
            ;;
        --help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            usage >&2
            exit 1
            ;;
    esac
done

if [ ${#RESOLUTIONS[@]} -eq 0 ]; then
    echo "Please provide at least one --resolution value." >&2
    exit 1
fi

if [ ${#MPI_COUNTS[@]} -eq 0 ]; then
    MPI_COUNTS=(0 2 4 8)
fi

# Deduplicate while preserving order.
deduped_mpi=()
for count in "${MPI_COUNTS[@]}"; do
    seen=false
    for existing in "${deduped_mpi[@]}"; do
        if [ "$existing" = "$count" ]; then
            seen=true
            break
        fi
    done
    if [ "$seen" = false ]; then
        deduped_mpi+=("$count")
    fi
done
MPI_COUNTS=("${deduped_mpi[@]}")

mkdir -p "$RESULTS_DIR"
printf "resolution,table_suffix,mode,mpi_ranks,log_file,table_init_s,time_loop_s,total_memory_mb\n" > "$RESULTS_CSV"

extract_metric() {
    local log_file="$1"
    local label="$2"

    awk -F'= ' -v label="$label" '
        index($0, label) {
            split($2, parts, " ");
            print parts[1];
            exit;
        }
    ' "$log_file"
}

extract_memory() {
    local log_file="$1"

    awk -F': ' '
        /Total memory across all processes:/ {
            split($2, parts, " ");
            print parts[1];
            exit;
        }
    ' "$log_file"
}

check_log_for_failure() {
    local log_file="$1"

    if grep -Eq "FOAM FATAL ERROR|MPI_ABORT|exited due to process rank|Segmentation fault|ERROR:" "$log_file"; then
        echo "Benchmark run failed according to log: $log_file" >&2
        exit 1
    fi
}

run_case() {
    local resolution="$1"
    local mpi_ranks="$2"
    local parallel=false
    local mode="serial"
    local log_file="$SCRIPT_DIR/log.bwRSE4HPCFoam"

    if [ "$mpi_ranks" -gt 0 ]; then
        parallel=true
        mode="mpi"
        log_file="$SCRIPT_DIR/log.mpirun"
    fi

    echo "=== Resolution $resolution | mode=$mode | ranks=$mpi_ranks ==="

    (
        cd "$SCRIPT_DIR"
        ./Allclean
        TABLERESOLUTION="$resolution" \
        TABLE_SUFFIX="$TABLE_SUFFIX" \
        PARALLEL="$parallel" \
        NPROCS="$mpi_ranks" \
        ./Allrun
    )

    if [ ! -f "$log_file" ]; then
        echo "Expected log file not found: $log_file" >&2
        exit 1
    fi

    check_log_for_failure "$log_file"

    local table_init_s
    local time_loop_s
    local total_memory_mb

    table_init_s="$(extract_metric "$log_file" "Total table loading/initialization time")"
    time_loop_s="$(extract_metric "$log_file" "Total time loop runtime")"
    total_memory_mb="$(extract_memory "$log_file")"

    if [ -z "$table_init_s" ] || [ -z "$time_loop_s" ] || [ -z "$total_memory_mb" ]; then
        echo "Missing expected benchmark metrics in log: $log_file" >&2
        exit 1
    fi

    printf "%s,%s,%s,%s,%s,%s,%s,%s\n" \
        "$resolution" \
        "$TABLE_SUFFIX" \
        "$mode" \
        "$mpi_ranks" \
        "$log_file" \
        "$table_init_s" \
        "$time_loop_s" \
        "$total_memory_mb" >> "$RESULTS_CSV"
}

for resolution in "${RESOLUTIONS[@]}"; do
    for mpi_ranks in "${MPI_COUNTS[@]}"; do
        run_case "$resolution" "$mpi_ranks"
    done
done

echo
echo "Benchmark sweep complete."
echo "Results written to: $RESULTS_CSV"
