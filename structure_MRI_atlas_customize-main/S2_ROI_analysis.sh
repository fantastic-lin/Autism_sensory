#!/usr/bin/env bash
#PBS -N ROI_analysis
#PBS -l nodes=1:ppn=6
#PBS -j oe
#PBS -q batch

set -euo pipefail

operation_dir="/histor/sun/linlin/4_olfactory/ABIDE/abide-master/sMRI/sMRI_git"
subject_dir="FS_git"
export SUBJECTS_DIR="${operation_dir}/${subject_dir}"

json_file="${operation_dir}/FS_successful_download.json"

# ------------------------------------------------------------
# Safety checks
# ------------------------------------------------------------
command -v jq >/dev/null 2>&1 || { echo "ERROR: jq not found in PATH"; exit 1; }
for cmd in mri_surf2surf mris_anatomical_stats; do
  command -v "$cmd" >/dev/null 2>&1 || { echo "ERROR: $cmd not found in PATH (is FreeSurfer sourced?)"; exit 1; }
done

[[ -d "$operation_dir" ]] || { echo "ERROR: operation_dir not found: $operation_dir"; exit 1; }
[[ -d "$SUBJECTS_DIR"   ]] || { echo "ERROR: SUBJECTS_DIR not found: $SUBJECTS_DIR"; exit 1; }
[[ -f "$json_file"      ]] || { echo "ERROR: JSON file not found: $json_file"; exit 1; }

# Load subject list safely (handles spaces/newlines)
mapfile -t subject_list < <(jq -r '.[]' "$json_file")
echo "Number of subjects: ${#subject_list[@]}"
echo "First few subjects: ${subject_list[@]:0:5}"

# Annotation must exist for both hemispheres
for he in lh rh; do
  annot_path="${SUBJECTS_DIR}/fsaverage/label/${he}.aal_olfactory.annot"
  [[ -f "$annot_path" ]] || { echo "ERROR: Missing annotation: $annot_path"; exit 1; }
done

# ------------------------------------------------------------
# Generate structural measurement file for each subject
# ------------------------------------------------------------
declare -a error_subjects_file=()
declare -a fail_extract_subj=()

counter_success_extract=0
hemis=(lh rh)

for subj in "${subject_list[@]}"; do
  he_ok=0

  # Basic subject directory checks
  subj_dir="${SUBJECTS_DIR}/${subj}"
  if [[ ! -d "$subj_dir" ]]; then
    echo "WARNING: Subject directory not found: $subj_dir"
    fail_extract_subj+=("$subj")
    continue
  fi

  mkdir -p "${subj_dir}/label" "${subj_dir}/stats"

  for he in "${hemis[@]}"; do
    out_stats="${subj_dir}/stats/${subj}_${he}.txt"
    out_annot="${subj_dir}/label/${he}.aal_olf_${subj}.annot"
    src_annot="${SUBJECTS_DIR}/fsaverage/label/${he}.aal_olfactory.annot"

    # If stats already exist and are non-empty, skip recomputation
    if [[ -s "$out_stats" ]]; then
      echo "[SKIP] ${subj} ${he}: stats exists: $out_stats"
      ((he_ok++))
      continue
    fi

    # Remove old empty/partial outputs (safe)
    rm -f "$out_stats" "$out_annot"

    echo "========== ${subj} ${he} =========="
    echo "Step 1: mri_surf2surf"
    if mri_surf2surf \
      --hemi "$he" \
      --srcsubject fsaverage \
      --sval-annot "$src_annot" \
      --trgsubject "$subj" \
      --trgsurfval "$out_annot"; then

      echo "Step 2: mris_anatomical_stats"
      if mris_anatomical_stats \
        -a "$out_annot" \
        -f "$out_stats" \
        -b "$subj" \
        "$subj" "$he"; then

        echo "[OK] Anatomical stats generated: ${subj} ${he}"
        ((he_ok++))
      else
        echo "[ERR] mris_anatomical_stats failed: ${subj} ${he}"
        error_subjects_file+=("${subj}_${he}.txt")
        # keep going to the next hemi/subject
      fi
    else
      echo "[ERR] mri_surf2surf failed: ${subj} ${he}"
      error_subjects_file+=("${subj}_${he}.surf2surf_failed")
      # continue
    fi
  done

  if [[ "$he_ok" -eq 2 ]]; then
    ((counter_success_extract++))
    echo "[DONE] ${subj} generated successfully"
  else
    fail_extract_subj+=("$subj")
    echo "[FAIL] ${subj} incomplete (he_ok=${he_ok})"
  fi
done

echo "============================================================"
echo "Error files (${#error_subjects_file[@]}): ${error_subjects_file[*]}"
echo "Failed subjects (${#fail_extract_subj[@]}): ${fail_extract_subj[*]}"
echo "Successful subjects: ${counter_success_extract} / ${#subject_list[@]}"