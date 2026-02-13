#!/usr/bin/env bash
set -euo pipefail

operation_dir="/histor/sun/linlin/4_olfactory/ABIDE/abide-master/sMRI/sMRI_git"
subject_dir="FS_git"
export SUBJECTS_DIR="${operation_dir}/${subject_dir}"

json_file="${operation_dir}/FS_successful_download.json"

command -v jq >/dev/null 2>&1 || { echo "ERROR: jq not found in PATH"; exit 1; }
[[ -f "$json_file" ]] || { echo "ERROR: JSON file not found: $json_file"; exit 1; }
[[ -d "$SUBJECTS_DIR" ]] || { echo "ERROR: SUBJECTS_DIR not found: $SUBJECTS_DIR"; exit 1; }

# Safer reading of JSON array into bash array
mapfile -t subject_list < <(jq -r '.[]' "$json_file")
echo "Number of subjects: ${#subject_list[@]}"

missing_files="${operation_dir}/FS_git_without_rh_lh_subj.csv"
: > "$missing_files"  # truncate/create the file

hemis=(lh rh)

for subj in "${subject_list[@]}"; do
  echo "$subj"
  for he in "${hemis[@]}"; do
    src_txt="${SUBJECTS_DIR}/${subj}/stats/${subj}_${he}.txt"
    out_txt="${SUBJECTS_DIR}/${subj}/stats/${subj}_${he}_1.txt"

    if [[ -f "$src_txt" ]]; then
      # Always (re)generate _1.txt from the desired lines
      # NOTE: awk without -F splits on any whitespace robustly
      sed -n '61,66p' "$src_txt" | awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' > "$out_txt"
      echo "Wrote: $out_txt"
    else
      echo "${subj}_${he}.txt" >> "$missing_files"
      echo "Missing: $src_txt (logged to $missing_files)"
    fi
  done
done

echo "Done. Missing list saved to: $missing_files"
