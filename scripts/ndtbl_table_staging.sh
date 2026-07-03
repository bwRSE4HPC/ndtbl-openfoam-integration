#!/bin/bash

ndtbl_staging_usage()
{
    cat <<'EOF'
Usage:
  ndtbl_table_staging.sh --table SOURCE CASE_PATH [TMP_NAME] [options]
  ndtbl_table_staging.sh --from-combustion-properties FILE [options]

Options:
  --table SOURCE CASE_PATH [TMP_NAME]
      Stage one pre-generated ndtbl file. Repeat for multiple files.
      SOURCE is the shared-filesystem source file.
      CASE_PATH is the path used by the OpenFOAM case dictionary.
      TMP_NAME is optional and defaults to basename(CASE_PATH).

  --from-combustion-properties FILE
      Read unique ndtblTables file entries from an OpenFOAM
      combustionProperties dictionary and stage them.

  --case-dir DIR
      Base directory for relative case table paths found in
      combustionProperties. Default: inferred from FILE when possible,
      otherwise the current directory.

  --source-dir DIR
      Base directory for relative source table paths found in
      combustionProperties. Default: case-dir.

  --tmp-subdir DIR
      Subdirectory below TMPDIR for staged table files. Default: ndtblTables

  --evict true|false
      Call posix_fadvise(..., DONTNEED) for staged files. Default: true

  --help
      Show this help.

Examples:
  scripts/ndtbl_table_staging.sh \
      --tmp-subdir myCase \
      --table /shared/tables/Tables.ndtbl constant/Tables.ndtbl

  scripts/ndtbl_table_staging.sh \
      --table /shared/tables/fuel.ndtbl constant/fuel.ndtbl \
      --table /shared/tables/oxidizer.ndtbl constant/oxidizer.ndtbl

  scripts/ndtbl_table_staging.sh \
      --tmp-subdir myCase \
      --from-combustion-properties constant/combustionProperties

Sourceable functions:
  ndtbl_stage_table SOURCE CASE_PATH [TMP_SUBDIR] [TMP_NAME] [EVICT]
  ndtbl_stage_combustion_tables COMBUSTION_PROPERTIES [CASE_DIR] [SOURCE_DIR] [TMP_SUBDIR] [EVICT]
  ndtbl_evict_file_cache PATH
EOF
}

ndtbl_evict_file_cache()
{
    local table_file="$1"

    python3 -c '
import os
import sys

path = sys.argv[1]
fd = os.open(path, os.O_RDONLY)
try:
    if hasattr(os, "posix_fadvise"):
        os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_DONTNEED)
finally:
    os.close(fd)
' "$table_file" || {
        echo "Warning: failed to evict table cache for $table_file" >&2
    }
}

ndtbl_stage_table()
{
    local table_source="$1"
    local case_table_file="$2"
    local table_tmp_subdir="${3:-${NDTBL_TABLE_TMP_SUBDIR:-ndtblTables}}"
    local table_tmp_name="${4:-$(basename -- "$case_table_file")}"
    local evict_cache="${5:-${NDTBL_TABLE_EVICT_CACHE:-true}}"

    if [ ! -f "$table_source" ]; then
        echo "ndtbl table source does not exist: $table_source" >&2
        return 1
    fi

    mkdir -p "$(dirname -- "$case_table_file")" || return 1

    if [ -n "${TMPDIR:-}" ]; then
        local table_tmp_dir="$TMPDIR/$table_tmp_subdir"
        local table_tmp_file="$table_tmp_dir/$table_tmp_name"

        if [ -n "${SLURM_JOB_ID:-}" ]; then
            if ! command -v srun >/dev/null 2>&1; then
                echo "TMPDIR staging in a Slurm job requires srun for per-node copies." >&2
                return 1
            fi

            echo "Copying $table_source to $table_tmp_file once per allocated node"
            local srun_args=(--ntasks-per-node=1)
            if [ -n "${SLURM_JOB_NUM_NODES:-}" ]; then
                srun_args+=(--nodes="$SLURM_JOB_NUM_NODES" --ntasks="$SLURM_JOB_NUM_NODES")
            fi

            srun "${srun_args[@]}" bash -c '
                set -e
                table_tmp_dir="$1"
                table_tmp_file="$2"
                table_source="$3"
                evict_cache="$4"

                mkdir -p "$table_tmp_dir"
                cp "$table_source" "$table_tmp_file"

                if [ "$evict_cache" = true ]; then
                    python3 -c '"'"'
import os
import sys

path = sys.argv[1]
fd = os.open(path, os.O_RDONLY)
try:
    if hasattr(os, "posix_fadvise"):
        os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_DONTNEED)
finally:
    os.close(fd)
'"'"' "$table_tmp_file"
                fi
            ' sh "$table_tmp_dir" "$table_tmp_file" "$table_source" "$evict_cache" || return 1
        else
            echo "Copying $table_source to local TMPDIR: $table_tmp_file"
            mkdir -p "$table_tmp_dir" || return 1
            cp "$table_source" "$table_tmp_file" || return 1
            if [ "$evict_cache" = true ]; then
                ndtbl_evict_file_cache "$table_tmp_file"
            fi
        fi

        rm -f "$case_table_file" || return 1
        echo "Linking $case_table_file -> $table_tmp_file"
        ln -s "$table_tmp_file" "$case_table_file" || return 1
    else
        if [ "$table_source" != "$case_table_file" ]; then
            echo "TMPDIR is not set; copying $table_source to $case_table_file"
            cp "$table_source" "$case_table_file" || return 1
        else
            echo "TMPDIR is not set; using existing table $case_table_file"
        fi
        if [ "$evict_cache" = true ]; then
            ndtbl_evict_file_cache "$case_table_file"
        fi
    fi
}

ndtbl_extract_combustion_table_files()
{
    local combustion_properties="$1"

    python3 -c '
import re
import sys

path = sys.argv[1]
text = open(path, encoding="utf-8").read()

def strip_comments(source):
    result = []
    index = 0
    in_string = False
    while index < len(source):
        char = source[index]
        next_char = source[index + 1] if index + 1 < len(source) else ""

        if in_string:
            result.append(char)
            if char == "\\" and next_char:
                result.append(next_char)
                index += 2
                continue
            if char == "\"":
                in_string = False
            index += 1
            continue

        if char == "\"":
            in_string = True
            result.append(char)
            index += 1
            continue

        if char == "/" and next_char == "/":
            index = source.find("\n", index)
            if index == -1:
                break
            result.append("\n")
            index += 1
            continue

        if char == "/" and next_char == "*":
            end = source.find("*/", index + 2)
            if end == -1:
                break
            result.append("\n" * source[index:end + 2].count("\n"))
            index = end + 2
            continue

        result.append(char)
        index += 1

    return "".join(result)

def matching_brace(source, open_index):
    depth = 0
    index = open_index
    in_string = False
    while index < len(source):
        char = source[index]
        next_char = source[index + 1] if index + 1 < len(source) else ""

        if in_string:
            if char == "\\" and next_char:
                index += 2
                continue
            if char == "\"":
                in_string = False
            index += 1
            continue

        if char == "\"":
            in_string = True
        elif char == "{":
            depth += 1
        elif char == "}":
            depth -= 1
            if depth == 0:
                return index
        index += 1
    return -1

clean = strip_comments(text)
seen = set()
tables = []

for match in re.finditer(r"\bndtblTables\b", clean):
    open_index = clean.find("{", match.end())
    if open_index == -1:
        continue
    close_index = matching_brace(clean, open_index)
    if close_index == -1:
        continue

    block = clean[open_index + 1:close_index]
    for file_match in re.finditer(
        r"(?:^|[\s{};])file\s+(?:\"([^\"]+)\"|([^;\s{}]+))\s*;",
        block,
    ):
        table = file_match.group(1) or file_match.group(2)
        if table not in seen:
            seen.add(table)
            tables.append(table)

for table in tables:
    print(table)
' "$combustion_properties"
}

ndtbl_infer_case_dir()
{
    local combustion_properties="$1"
    local properties_dir
    properties_dir="$(CDPATH= cd -- "$(dirname -- "$combustion_properties")" && pwd)"

    if [ "$(basename -- "$properties_dir")" = "constant" ]; then
        dirname -- "$properties_dir"
    else
        pwd
    fi
}

ndtbl_absolute_or_join()
{
    local base_dir="$1"
    local path="$2"

    case "$path" in
        /*)
            printf "%s\n" "$path"
            ;;
        *)
            printf "%s/%s\n" "$base_dir" "$path"
            ;;
    esac
}

ndtbl_tmp_name_from_path()
{
    local path="$1"
    path="${path#/}"
    path="${path//\//__}"
    if [ -n "$path" ]; then
        printf "%s\n" "$path"
    else
        basename -- "$1"
    fi
}

ndtbl_stage_combustion_tables()
{
    local combustion_properties="$1"
    local case_dir="${2:-}"
    local source_dir="${3:-}"
    local table_tmp_subdir="${4:-${NDTBL_TABLE_TMP_SUBDIR:-ndtblTables}}"
    local evict_cache="${5:-${NDTBL_TABLE_EVICT_CACHE:-true}}"

    if [ ! -f "$combustion_properties" ]; then
        echo "combustionProperties file does not exist: $combustion_properties" >&2
        return 1
    fi

    if [ -z "$case_dir" ]; then
        case_dir="$(ndtbl_infer_case_dir "$combustion_properties")"
    fi
    if [ -z "$source_dir" ]; then
        source_dir="$case_dir"
    fi

    local staged_count=0
    local table_path
    while IFS= read -r table_path; do
        [ -n "$table_path" ] || continue

        local table_source
        local case_table_file
        local table_tmp_name
        table_source="$(ndtbl_absolute_or_join "$source_dir" "$table_path")"
        case_table_file="$(ndtbl_absolute_or_join "$case_dir" "$table_path")"
        table_tmp_name="$(ndtbl_tmp_name_from_path "$table_path")"

        ndtbl_stage_table \
            "$table_source" \
            "$case_table_file" \
            "$table_tmp_subdir" \
            "$table_tmp_name" \
            "$evict_cache" || return 1
        staged_count=$((staged_count + 1))
    done < <(ndtbl_extract_combustion_table_files "$combustion_properties")

    if [ "$staged_count" -eq 0 ]; then
        echo "No ndtblTables file entries found in $combustion_properties" >&2
        return 1
    fi
}

ndtbl_stage_tables_main()
{
    local table_tmp_subdir="${NDTBL_TABLE_TMP_SUBDIR:-ndtblTables}"
    local evict_cache="${NDTBL_TABLE_EVICT_CACHE:-true}"
    local combustion_properties=""
    local case_dir=""
    local source_dir=""
    local table_sources=()
    local case_table_files=()
    local table_tmp_names=()

    while [ "$#" -gt 0 ]; do
        case "$1" in
            --table)
                if [ "$#" -lt 3 ]; then
                    echo "--table requires SOURCE and CASE_PATH" >&2
                    return 1
                fi

                local source="$2"
                local case_path="$3"
                local tmp_name=""
                shift 3

                if [ "$#" -gt 0 ] && [ "${1#--}" = "$1" ]; then
                    tmp_name="$1"
                    shift
                fi

                table_sources+=("$source")
                case_table_files+=("$case_path")
                table_tmp_names+=("$tmp_name")
                ;;
            --from-combustion-properties)
                if [ "$#" -lt 2 ]; then
                    echo "--from-combustion-properties requires a file path" >&2
                    return 1
                fi
                combustion_properties="$2"
                shift 2
                ;;
            --case-dir)
                if [ "$#" -lt 2 ]; then
                    echo "--case-dir requires a value" >&2
                    return 1
                fi
                case_dir="$2"
                shift 2
                ;;
            --source-dir)
                if [ "$#" -lt 2 ]; then
                    echo "--source-dir requires a value" >&2
                    return 1
                fi
                source_dir="$2"
                shift 2
                ;;
            --tmp-subdir)
                if [ "$#" -lt 2 ]; then
                    echo "--tmp-subdir requires a value" >&2
                    return 1
                fi
                table_tmp_subdir="$2"
                shift 2
                ;;
            --evict)
                if [ "$#" -lt 2 ]; then
                    echo "--evict requires true or false" >&2
                    return 1
                fi
                evict_cache="$2"
                shift 2
                ;;
            --help)
                ndtbl_staging_usage
                return 0
                ;;
            *)
                echo "Unknown option: $1" >&2
                ndtbl_staging_usage >&2
                return 1
                ;;
        esac
    done

    if [ "${#table_sources[@]}" -eq 0 ] && [ -z "$combustion_properties" ]; then
        echo "Please provide --from-combustion-properties or at least one --table SOURCE CASE_PATH entry." >&2
        ndtbl_staging_usage >&2
        return 1
    fi

    if [ -n "$combustion_properties" ]; then
        ndtbl_stage_combustion_tables \
            "$combustion_properties" \
            "$case_dir" \
            "$source_dir" \
            "$table_tmp_subdir" \
            "$evict_cache" || return 1
    fi

    local index
    for index in "${!table_sources[@]}"; do
        ndtbl_stage_table \
            "${table_sources[$index]}" \
            "${case_table_files[$index]}" \
            "$table_tmp_subdir" \
            "${table_tmp_names[$index]}" \
            "$evict_cache" || return 1
    done
}

if [ "${BASH_SOURCE[0]}" = "$0" ]; then
    set -euo pipefail
    ndtbl_stage_tables_main "$@"
fi
