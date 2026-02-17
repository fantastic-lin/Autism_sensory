#!/bin/bash
set -euo pipefail

# ============================================================
# 0) Load config and set paths
# ============================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.env"

: "${operation_dir:?operation_dir is not set in config.env}"

input_dir="${input_dir:-$operation_dir/input}"
output_dir="${output_dir:-$operation_dir/output}"

echo "operation_dir=$operation_dir"
echo "input_dir=$input_dir"
echo "output_dir=$output_dir"

# Where per-subject stats files are stored
stats_store_dir="${output_dir}/S2345_output_stats_result"

# JSON file containing a list of subject directories (array of strings)
json_file="${input_dir}/S2345_input_ABIDE1_qc_OK.json"
[ -f "$json_file" ] || { echo "ERROR: missing json_file: $json_file"; exit 1; }

# Read subject list safely (no word-splitting issues)
command -v jq >/dev/null 2>&1 || { echo "ERROR: jq not found in PATH"; exit 1; }
mapfile -t subject_list < <(jq -r '.[]' "$json_file")

# Record missing stats files here
missing_files="${stats_store_dir}/ABIDE1_without_lh_rh_subj.csv"
mkdir -p "$stats_store_dir"
: > "$missing_files"   # truncate existing file (start fresh)

# Hemispheres
hemis=(lh rh)

# ============================================================
# 1) Loop over subjects and extract lines 61–88 from stats files
# ============================================================
for subj in "${subject_list[@]}"; do
  # site_name = parent folder name of subj
  site_name="$(basename "$(dirname "$subj")")"

  # subject_name = subj basename without trailing "_total"
  subject_name="$(basename "$subj" | sed 's/_total$//')"

  for he in "${hemis[@]}"; do
    # Stats directory for this subject
    stats_path="${stats_store_dir}/${site_name}/${subject_name}/stats"
    mkdir -p -- "$stats_path"

    # Input/Output file names
    src="${stats_path}/${subject_name}_${he}.txt"
    dst="${stats_path}/${subject_name}_${he}_1.txt"

    if [ -f "$src" ]; then
      # Extract lines 61–88 and keep first 10 whitespace-separated columns
      rm -f -- "$dst"
      sed -n '61,88p' "$src" | awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' > "$dst"
      echo "Wrote: $dst"
    else
      # Record missing file path
      echo "$src" >> "$missing_files"
      echo "Missing: $src (skipped)"
    fi
  done
done

echo "Done."
echo "Missing list saved to: $missing_files"
