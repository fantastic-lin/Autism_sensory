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

stats_store_dir="${output_dir}/S2345_output_stats_result"
json_file="${input_dir}/S2345_input_ABIDE1_qc_OK.json"

command -v jq >/dev/null 2>&1 || { echo "ERROR: jq not found in PATH"; exit 1; }
[ -f "$json_file" ] || { echo "ERROR: missing json file: $json_file"; exit 1; }

mkdir -p "$stats_store_dir"

mapfile -t subject_list < <(jq -r '.[]' "$json_file")
echo "Loaded subjects: ${#subject_list[@]}"

version="ABIDE1"

# ============================================================
# 1) Find a reference aparc.stats to build a unified ROI list
# ============================================================
first_lh=""
first_rh=""
for subj in "${subject_list[@]}"; do
  site_name="$(basename "$(dirname "$subj")")"
  subject_name="$(basename "$subj" | sed 's/_total$//')"
  subj_dir="${subj}/${subject_name}_freesurfer"
  stats_dir="${subj_dir}/${subject_name}/stats"
  if [[ -f "$stats_dir/lh.aparc.stats" && -f "$stats_dir/rh.aparc.stats" ]]; then
    first_lh="$stats_dir/lh.aparc.stats"
    first_rh="$stats_dir/rh.aparc.stats"
    break
  fi
done

if [[ -z "$first_lh" ]]; then
  echo "ERROR: cannot find any lh.aparc.stats to build ROI header."
  exit 1
fi

# ROI table rows: keep your original logic (54–87)
# ROI name assumed in col1.
# Filter out Hippocampus/Amygdala/Thalamus if they appear (safety).
roi_names=$(sed -n '54,87p' "$first_lh" | awk '{print $1}' \
  | grep -Ev '(Hippocampus|Amygdala|Thalamus)')

# ============================================================
# 2) Surface Area CSV (same ROI set as volume)
#    - SurfArea assumed in col3
# ============================================================
surface_area_csv="${stats_store_dir}/${version}_surface_area.csv"
: > "$surface_area_csv"

# Header: SubjectID + lh_ + rh_ ROI SurfArea columns
echo -n "SubjectID" >> "$surface_area_csv"
for roi in $roi_names; do
  echo -n ",lh_${roi}_SurfArea" >> "$surface_area_csv"
done
for roi in $roi_names; do
  echo -n ",rh_${roi}_SurfArea" >> "$surface_area_csv"
done
echo -n ",TotalSurfaceArea" >> "$surface_area_csv"
echo >> "$surface_area_csv"

for subj in "${subject_list[@]}"; do
  site_name="$(basename "$(dirname "$subj")")"
  subject_name="$(basename "$subj" | sed 's/_total$//')"

  subj_dir="${subj}/${subject_name}_freesurfer"
  stats_dir="${subj_dir}/${subject_name}/stats"
  f_lh="${stats_dir}/lh.aparc.stats"
  f_rh="${stats_dir}/rh.aparc.stats"

  subject_whole_name="${version}_${site_name}_${subject_name}"

  if [[ ! -f "$f_lh" || ! -f "$f_rh" ]]; then
    echo "Missing aparc.stats for ${site_name}/${subject_name}:"
    [[ -f "$f_lh" ]] || echo "  missing $f_lh"
    [[ -f "$f_rh" ]] || echo "  missing $f_rh"
    continue
  fi

  # Build values in SAME ROI order as header
  # We read ROI rows (54–87), filter same way, then take col3 (SurfArea).
  lh_vals=$(sed -n '54,87p' "$f_lh" | awk '{print $1,$3}' \
    | grep -Ev '(Hippocampus|Amygdala|Thalamus)' \
    | awk '{printf "%s%s", (NR==1?"":","), $2}')
  rh_vals=$(sed -n '54,87p' "$f_rh" | awk '{print $1,$3}' \
    | grep -Ev '(Hippocampus|Amygdala|Thalamus)' \
    | awk '{printf "%s%s", (NR==1?"":","), $2}')

  total_area=$(awk -v a="$lh_vals" -v b="$rh_vals" 'BEGIN{print ""}')
  echo "${subject_whole_name},${lh_vals},${rh_vals}" >> "$surface_area_csv"
done

echo "Surface area saved to: $surface_area_csv"

# ============================================================
# 3) Gray Matter Volume CSV (same ROI columns as surface area)
#    - GrayVol assumed in col4
#    - NO Hippocampus/Amygdala/Thalamus
# ============================================================
volume_csv="${stats_store_dir}/${version}_gray_matter_volume.csv"
: > "$volume_csv"

# Header must match ROI set used in surface area
echo -n "SubjectID" >> "$volume_csv"
for roi in $roi_names; do
  echo -n ",lh_${roi}_GrayVol" >> "$volume_csv"
done
for roi in $roi_names; do
  echo -n ",rh_${roi}_GrayVol" >> "$volume_csv"
done
echo >> "$volume_csv"

for subj in "${subject_list[@]}"; do
  site_name="$(basename "$(dirname "$subj")")"
  subject_name="$(basename "$subj" | sed 's/_total$//')"

  subj_dir="${subj}/${subject_name}_freesurfer"
  stats_dir="${subj_dir}/${subject_name}/stats"
  f_lh="${stats_dir}/lh.aparc.stats"
  f_rh="${stats_dir}/rh.aparc.stats"

  subject_whole_name="${version}_${site_name}_${subject_name}"

  if [[ ! -f "$f_lh" || ! -f "$f_rh" ]]; then
    echo "Missing aparc.stats for ${site_name}/${subject_name}:"
    [[ -f "$f_lh" ]] || echo "  missing $f_lh"
    [[ -f "$f_rh" ]] || echo "  missing $f_rh"
    continue
  fi

  # Same ROI order, take col4 (GrayVol)
  lh_vals=$(sed -n '54,87p' "$f_lh" | awk '{print $1,$4}' \
    | grep -Ev '(Hippocampus|Amygdala|Thalamus)' \
    | awk '{printf "%s%s", (NR==1?"":","), $2}')
  rh_vals=$(sed -n '54,87p' "$f_rh" | awk '{print $1,$4}' \
    | grep -Ev '(Hippocampus|Amygdala|Thalamus)' \
    | awk '{printf "%s%s", (NR==1?"":","), $2}')

  echo "${subject_whole_name},${lh_vals},${rh_vals}" >> "$volume_csv"
done

echo "Gray matter volume saved to: $volume_csv"
echo "Done."

