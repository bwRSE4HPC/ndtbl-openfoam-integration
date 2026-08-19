#!/bin/bash

set -euo pipefail

SCRIPT_DIR="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
RESULTS_DIR="$SCRIPT_DIR/benchmark-results"
RESULTS_CSV="$RESULTS_DIR/results.csv"
TABLE_SUFFIX="f"
RESOLUTIONS=()
MPI_COUNTS=()
TABLE_LAYOUTS=()

usage() {
    cat <<'EOF'
Usage:
  run_ndtbl_benchmarks.sh [options]

Options:
  --resolution N     Add one table resolution to sweep. Repeatable.
  --mpi N            Add one MPI process count. Use 0 or 'none' for serial.
                     Repeatable. Default: 0 (serial).
  --table-suffix S   ndtbl filename suffix to use when TABLE_FILE is not set
                     via Allrun. Default: f
  --table-layout L   Add a table layout to sweep: combined or split.
                     Repeatable. Default: combined
  --results-dir DIR  Directory for CSV results. Default: benchmark-results
  --help             Show this help

Examples:
  ./run_ndtbl_benchmarks.sh --resolution 21 --resolution 31
  ./run_ndtbl_benchmarks.sh --resolution 21 --mpi 0 --mpi 2 --mpi 4 --mpi 8
  ./run_ndtbl_benchmarks.sh --resolution 101 \
      --table-layout combined --table-layout split
  ./run_ndtbl_benchmarks.sh --resolution 501 --table-layout split --mpi 5
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
        --table-layout)
            if [ "$2" != combined ] && [ "$2" != split ]; then
                echo "Invalid table layout: $2 (expected combined or split)" >&2
                exit 1
            fi
            TABLE_LAYOUTS+=("$2")
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
    MPI_COUNTS=(0)
fi

if [ ${#TABLE_LAYOUTS[@]} -eq 0 ]; then
    TABLE_LAYOUTS=(combined)
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

deduped_layouts=()
for layout in "${TABLE_LAYOUTS[@]}"; do
    seen=false
    for existing in "${deduped_layouts[@]}"; do
        if [ "$existing" = "$layout" ]; then
            seen=true
            break
        fi
    done
    if [ "$seen" = false ]; then
        deduped_layouts+=("$layout")
    fi
done
TABLE_LAYOUTS=("${deduped_layouts[@]}")

mkdir -p "$RESULTS_DIR"
printf "resolution,table_suffix,table_layout,table_file_count,total_table_bytes,mode,mpi_ranks,log_file,table_init_s,time_loop_s,total_memory_mb\n" > "$RESULTS_CSV"

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
    local table_layout="$2"
    local mpi_ranks="$3"
    local parallel=false
    local mode="serial"
    local log_file="$SCRIPT_DIR/log.bwRSE4HPCFoam"

    if [ "$mpi_ranks" -gt 0 ]; then
        parallel=true
        mode="mpi"
    fi

    echo "=== Resolution $resolution | layout=$table_layout | mode=$mode | ranks=$mpi_ranks ==="

    (
        cd "$SCRIPT_DIR"
        ./Allclean
        TABLERESOLUTION="$resolution" \
        TABLE_SUFFIX="$TABLE_SUFFIX" \
        TABLE_LAYOUT="$table_layout" \
        PARALLEL="$parallel" \
        NPROCS="$mpi_ranks" \
        ./Allrun
    )

    if [ ! -f "$log_file" ]; then
        echo "Expected log file not found: $log_file" >&2
        exit 1
    fi

    check_log_for_failure "$log_file"

    local saved_log_file
    saved_log_file="$RESULTS_DIR/log-r${resolution}-${TABLE_SUFFIX}-${table_layout}-${mode}-${mpi_ranks}.txt"
    cp "$log_file" "$saved_log_file"

    local table_init_s
    local time_loop_s
    local total_memory_mb
    local table_file_count=1
    local total_table_bytes=0

    local table_files=("$SCRIPT_DIR/constant/Tables.ndtbl")
    if [ "$table_layout" = split ]; then
        table_file_count=3
        table_files=(
            "$SCRIPT_DIR/constant/Tables.Table1.ndtbl"
            "$SCRIPT_DIR/constant/Tables.Table2.ndtbl"
            "$SCRIPT_DIR/constant/Tables.Table3.ndtbl"
        )
    fi

    local table_file
    for table_file in "${table_files[@]}"; do
        if [ ! -f "$table_file" ]; then
            echo "Expected table file not found: $table_file" >&2
            exit 1
        fi
        local table_bytes
        table_bytes="$(wc -c < "$table_file")"
        total_table_bytes=$((total_table_bytes + table_bytes))
    done

    table_init_s="$(extract_metric "$log_file" "Total table loading/initialization time")"
    time_loop_s="$(extract_metric "$log_file" "Total time loop runtime")"
    total_memory_mb="$(extract_memory "$log_file")"

    if [ -z "$table_init_s" ] || [ -z "$time_loop_s" ] || [ -z "$total_memory_mb" ]; then
        echo "Missing expected benchmark metrics in log: $log_file" >&2
        exit 1
    fi

    printf "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" \
        "$resolution" \
        "$TABLE_SUFFIX" \
        "$table_layout" \
        "$table_file_count" \
        "$total_table_bytes" \
        "$mode" \
        "$mpi_ranks" \
        "$saved_log_file" \
        "$table_init_s" \
        "$time_loop_s" \
        "$total_memory_mb" >> "$RESULTS_CSV"
}

for resolution in "${RESOLUTIONS[@]}"; do
    for table_layout in "${TABLE_LAYOUTS[@]}"; do
        for mpi_ranks in "${MPI_COUNTS[@]}"; do
            run_case "$resolution" "$table_layout" "$mpi_ranks"
        done
    done
done

echo
echo "Benchmark sweep complete."
echo "Results written to: $RESULTS_CSV"
