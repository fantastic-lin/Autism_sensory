#!/bin/bash
set -euo pipefail

# ============================================================
# Environment
# ============================================================
set +u
source "/histor/sun/linlin/mambaforge/etc/profile.d/conda.sh"
conda activate fmri_env
set -u

# ============================================================
# Load config
#   config.env must define: operation_dir
#   (optional) input_dir, output_dir
# ============================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.env"

: "${operation_dir:?operation_dir is not set in config.env}"

# Default input/output under operation_dir (allow override from config.env)
input_dir="${input_dir:-$operation_dir/input}"
output_dir="${output_dir:-$operation_dir/output}"

echo "operation_dir=$operation_dir"
echo "input_dir=$input_dir"
echo "output_dir=$output_dir"

# ============================================================
# FreeSurfer setup
#   SUBJECTS_DIR points to where fsaverage lives for srcsubject.
# ============================================================
subject_dir="input"  # Change to your subject directory if needed
export SUBJECTS_DIR="${operation_dir}/${subject_dir}"

# ============================================================
# Inputs/Outputs
# ============================================================
stats_store_file="${operation_dir}/output/S2345_output_stats_result"
json_file="${operation_dir}/input/S2345_input_ABIDE1_qc_OK.json"

# ============================================================
# Basic checks
# ============================================================
[ -d "$SUBJECTS_DIR/fsaverage" ] || {
  echo "ERROR: missing fsaverage under SUBJECTS_DIR: $SUBJECTS_DIR/fsaverage"
  exit 1
}
[ -f "$json_file" ] || { echo "ERROR: missing json_file: $json_file"; exit 1; }
command -v jq >/dev/null 2>&1 || { echo "ERROR: jq not found in PATH"; exit 1; }

mkdir -p "$stats_store_file"

# ============================================================
# Read subject list from JSON safely
#   JSON format is expected to be an array of strings (paths).
# ============================================================
mapfile -t subject_list < <(jq -r '.[]' "$json_file")
echo "subject_list: ${subject_list[*]}"

# ============================================================
# Collect failures
# ============================================================
declare -a error_subjects_file=()   # store error tags
declare -a fail_extract_subj=()     # store failed subjects
counter_success_extract=0

hemis=(lh rh)

# ============================================================
# Main loop over subjects
# ============================================================
for subj in "${subject_list[@]}"; do
  echo "----------------------------------------"
  echo "subj_path=$subj"

  he_number=0

  # site folder name = parent directory name of $subj
  site_name="$(basename "$(dirname "$subj")")"

  # subject_name = basename of $subj with suffix "_total" removed
  subject_name="$(basename "$subj" | sed 's/_total$//')"

  echo "site_name=$site_name"
  echo "subject_name=$subject_name"

  # Output folders for this subject
  stats_path="$stats_store_file/${site_name}/${subject_name}/stats"
  label_path="$stats_store_file/${site_name}/${subject_name}/label"
  mkdir -p "$stats_path" "$label_path"

  # Each subject's FreeSurfer directory (where recon-all output is)
  subj_dir="${subj}/${subject_name}_freesurfer"
  echo "subj_dir=$subj_dir"

  # If the subject FreeSurfer directory is missing, mark failure and skip
  if [ ! -d "$subj_dir" ]; then
    echo "!!!!ERROR: subject freesurfer dir not found: $subj_dir"
    fail_extract_subj+=("${site_name}_${subject_name}")
    continue
  fi

  # ----------------------------------------------------------
  # Loop over hemispheres
  # ----------------------------------------------------------
  for he in "${hemis[@]}"; do
    out_stats="$stats_path/${subject_name}_${he}.txt"
    out_annot="$label_path/${he}.aal_sens_${subject_name}.annot"

    # If stats already exist, skip this hemi
    # (If you also want to require out_annot exists, change the condition.)
    if [ -f "$out_stats" ]; then
      echo "Skip existing: $out_stats"
      ((++he_number))
      continue
    fi

    echo "************ step 1: mri_surf2surf *******************"

    # Source annotation on fsaverage (built earlier)
    fsavg_annot="$output_dir/fsaverage/label/${he}.aal_sensory.annot"
    if [ ! -r "$fsavg_annot" ]; then
      echo "!!!!ERROR: missing fsaverage annot: $fsavg_annot"
      error_subjects_file+=("${site_name}_${subject_name}_${he}_missing_fsavg_annot")
      continue
    fi

    # Use the subject's FreeSurfer directory as SUBJECTS_DIR for the target subject
    if SUBJECTS_DIR="$subj_dir" mri_surf2surf \
        --hemi "$he" \
        --srcsubject fsaverage \
        --sval-annot "$fsavg_annot" \
        --trgsubject "$subject_name" \
        --trgsurfval "$out_annot"; then

      echo "************ step 2: mris_anatomical_stats ************"

      if SUBJECTS_DIR="$subj_dir" mris_anatomical_stats \
          -a "$out_annot" \
          -f "$out_stats" \
          -b "$subject_name" \
          "$he"; then

        echo "Anatomical stats generation successful: $out_stats"
        ((++he_number))

      else
        echo "!!!!ERROR: mris_anatomical_stats failed for ${subject_name} ${he}"
        error_subjects_file+=("${site_name}_${subject_name}_${he}_stats_fail")
      fi

    else
      echo "!!!!ERROR: mri_surf2surf failed for ${subject_name} ${he}"
      error_subjects_file+=("${site_name}_${subject_name}_${he}_annot_fail")
    fi
  done

  # Subject summary: both hemispheres succeeded
  if [ "$he_number" -eq 2 ]; then
    ((++counter_success_extract))
    echo "${subject_name} generated successfully"
  else
    fail_extract_subj+=("${site_name}_${subject_name}")
  fi
done

# ============================================================
# Final summary
# ============================================================
echo "========================================"
echo "error files: ${error_subjects_file[*]-}"
echo "fail generated subjects: ${fail_extract_subj[*]-}"
echo "successful subjects: ${counter_success_extract}"
