#!/usr/bin/env bash
set -euo pipefail

operation_dir="/histor/sun/linlin/4_olfactory/ABIDE/abide-master/sMRI/sMRI_git"
subject_dir="FS_git"
export SUBJECTS_DIR="${operation_dir}/${subject_dir}"

# -----------------------------
# Output files (overwrite each run)
# -----------------------------
area_out="${operation_dir}/FS_git_surface_area.csv"
vol_out="${operation_dir}/FS_git_gray_matter_volume.csv"

# Write headers
echo "subject,total_surface_area" > "$area_out"
echo "subject,TotalGrayVol,Left_Hippocampus,Right_Hippocampus,Left_Amygdala,Right_Amygdala" > "$vol_out"

# -----------------------------
# Iterate subjects (exclude fsaverage)
# -----------------------------
while IFS= read -r subj_path; do
  subj="$(basename "$subj_path")"
  echo "Processing: $subj"

  # ===== Surface Area =====
  lh_stats="${SUBJECTS_DIR}/${subj}/stats/lh.aparc.stats"
  rh_stats="${SUBJECTS_DIR}/${subj}/stats/rh.aparc.stats"

  if [[ -f "$lh_stats" && -f "$rh_stats" ]]; then
    # Sum column 3 (surf area) from lines 54-87 (your original logic)
    lh_sum="$(sed -n '54,87p' "$lh_stats" | awk '{sum+=$3} END{printf "%.6f", sum}')"
    rh_sum="$(sed -n '54,87p' "$rh_stats" | awk '{sum+=$3} END{printf "%.6f", sum}')"

    total_sum="$(awk -v a="$lh_sum" -v b="$rh_sum" 'BEGIN{printf "%.6f", a+b}')"
    echo "${subj},${total_sum}" >> "$area_out"
  else
    echo "WARNING: missing aparc.stats for $subj (lh/rh). Skipping surface area."
  fi

  # ===== Volume =====
  aseg_stats="${SUBJECTS_DIR}/${subj}/stats/aseg.stats"
