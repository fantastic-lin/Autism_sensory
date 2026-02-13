#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Config
# ============================================================
operation_dir="/histor/sun/linlin/4_olfactory/ABIDE/abide-master/sMRI/sMRI_git"  # change if needed
subject_dir="FS_git"                                                           # change if needed

export SUBJECTS_DIR="${operation_dir}/${subject_dir}"

# ============================================================
# Safety checks
# ============================================================
for cmd in mri_vol2surf mri_vol2label mris_label2annot; do
  command -v "$cmd" >/dev/null 2>&1 || { echo "ERROR: '$cmd' not found in PATH"; exit 1; }
done

[[ -d "$operation_dir" ]] || { echo "ERROR: operation_dir not found: $operation_dir"; exit 1; }
[[ -d "$SUBJECTS_DIR" ]]   || { echo "ERROR: SUBJECTS_DIR not found: $SUBJECTS_DIR"; exit 1; }

aal_nii="${operation_dir}/AAL116_1mm.nii"
[[ -f "$aal_nii" ]] || { echo "ERROR: AAL NIfTI not found: $aal_nii"; exit 1; }

# Color tables (ctab) must exist for BOTH hemispheres if you use hemi-specific files
ctab_lh="${operation_dir}/lh_aal_olfactory_RGB_note_1.txt"
ctab_rh="${operation_dir}/rh_aal_olfactory_RGB_note_1.txt"
[[ -f "$ctab_lh" ]] || { echo "ERROR: missing ctab: $ctab_lh"; exit 1; }
[[ -f "$ctab_rh" ]] || { echo "ERROR: missing ctab: $ctab_rh"; exit 1; }

# Where labels should be written (FreeSurfer expects under $SUBJECTS_DIR/fsaverage/label)
label_dir="${SUBJECTS_DIR}/fsaverage/label"
mkdir -p "$label_dir"

# ============================================================
# 1) Project AAL volume to fsaverage surface (creates .mgh per hemi)
# ============================================================
# Keep outputs alongside this script / in CWD; you can change paths if you want.
mri_vol2surf --mov "$aal_nii" --mni152reg --hemi rh --out_type mgh --o "rh.aal3.mgh"
mri_vol2surf --mov "$aal_nii" --mni152reg --hemi lh --out_type mgh --o "lh.aal3.mgh"

# ============================================================
# 2) AAL index lists
# ============================================================
# AAL odd = left hemi, even = right hemi.
olf_ids_l=(5 9 15 21 25 29 37 39 41)
olf_ids_r=(6 10 16 22 26 30 38 40 42)

# NOTE: you stated AAL left corresponds to FreeSurfer right, so we swap mapping below.

# ============================================================
# 3) For each hemisphere: create labels then merge into annot
# ============================================================
for hemi in lh rh; do
  if [[ "$hemi" == "lh" ]]; then
    olf_index=("${olf_ids_r[@]}")   # swap: FS lh uses AAL right indices
    ctab="$ctab_lh"
  else
    olf_index=("${olf_ids_l[@]}")   # swap: FS rh uses AAL left indices
    ctab="$ctab_rh"
  fi

  echo "Running hemisphere: ${hemi}; AAL indices: ${olf_index[*]}"

  # Ensure the projected file exists
  surf_mgh="${hemi}.aal3.mgh"
  [[ -f "$surf_mgh" ]] || { echo "ERROR: missing ${surf_mgh}. Did mri_vol2surf succeed?"; exit 1; }

  # ---- Create label files under $SUBJECTS_DIR/fsaverage/label
  label_paths=()
  for id in "${olf_index[@]}"; do
    out_label="${label_dir}/${hemi}.num${id}.label"
    mri_vol2label --i "$surf_mgh" \
      --id "$id" \
      --l "$out_label" \
      --surf fsaverage "$hemi"
    label_paths+=("$out_label")
  done

  # ---- Merge labels into annotation
  # Use an array to pass multiple --l arguments safely
  mris_label2annot --s fsaverage --h "$hemi" \
    --ctab "$ctab" \
    --a aal_olfactory \
    $(printf -- '--l %q ' "${label_paths[@]}")

done

echo "Done. Annotation files should be under: ${SUBJECTS_DIR}/fsaverage/label"
