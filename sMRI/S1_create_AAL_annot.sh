#!/bin/bash
set -euo pipefail

# ============================================================
# 0) Load config
# ============================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.env"

# Required: operation_dir
: "${operation_dir:?operation_dir is not set in config.env}"

# Default input/output under operation_dir (allow override from config.env)
input_dir="${input_dir:-$operation_dir/input}"
output_dir="${output_dir:-$operation_dir/output}"

echo "operation_dir=$operation_dir"
echo "input_dir=$input_dir"
echo "output_dir=$output_dir"

# Check required inputs
[ -f "$input_dir/AAL116_1mm.nii" ] || { echo "ERROR: missing $input_dir/AAL116_1mm.nii"; exit 1; }
[ -f "$input_dir/lh_aal_sensory_RGB_note.txt" ] || { echo "ERROR: missing $input_dir/lh_aal_sensory_RGB_note.txt"; exit 1; }
[ -f "$input_dir/rh_aal_sensory_RGB_note.txt" ] || { echo "ERROR: missing $input_dir/rh_aal_sensory_RGB_note.txt"; exit 1; }

# ============================================================
# 1) Use input as the FreeSurfer SUBJECTS_DIR (generation dir)
#    mris_label2annot outputs to: $SUBJECTS_DIR/fsaverage/label/
# ============================================================
export SUBJECTS_DIR="$operation_dir/input"
[ -d "$SUBJECTS_DIR/fsaverage" ] || { echo "ERROR: missing $SUBJECTS_DIR/fsaverage"; exit 1; }
mkdir -p "$SUBJECTS_DIR/fsaverage/label"

# ============================================================
# 2) Prepare output folders (store intermediate + copied outputs)
# ============================================================
mkdir -p "$output_dir/fsaverage/mgh" "$output_dir/fsaverage/label"

# ============================================================
# 3) Project AAL volume to fsaverage surface (intermediate to output_dir)
# ============================================================
mri_vol2surf --mov "$input_dir/AAL116_1mm.nii" --mni152reg --hemi rh --out_type mgh \
  --o "$output_dir/fsaverage/mgh/rh.aal3.mgh"

mri_vol2surf --mov "$input_dir/AAL116_1mm.nii" --mni152reg --hemi lh --out_type mgh \
  --o "$output_dir/fsaverage/mgh/lh.aal3.mgh"

# ============================================================
# 4) AAL sensory IDs
#    Odd -> left in AAL; even -> right in AAL
#    Note: AAL left corresponds to FreeSurfer right (as you observed)
# ============================================================
sens_ids_l=(41 43 31 33 35 45 15 25 9 5 55 79 37 29 47 53 51 49 21 39 61 59 57 67 17 83 81 77)
sens_ids_R=(42 44 32 34 36 46 16 26 10 6 56 80 38 30 48 54 52 50 22 40 62 60 58 68 18 84 82 78)

# ============================================================
# 5) Generate labels (to output_dir) + build annot (to SUBJECTS_DIR),
#    then copy annot/ctab to output_dir
# ============================================================
hemis=(lh rh)
for he in "${hemis[@]}"; do
  if [ "$he" == "lh" ]; then
    sens_index=("${sens_ids_R[@]}")
  else
    sens_index=("${sens_ids_l[@]}")
  fi

  echo "Processing hemisphere: $he"
  echo "AAL ids: ${sens_index[*]}"

  # 5.1 Create labels (store labels in output_dir)
  for id in "${sens_index[@]}"; do
    out_label="$output_dir/fsaverage/label/${he}.num${id}.label"
    mri_vol2label --i "$output_dir/fsaverage/mgh/${he}.aal3.mgh" \
      --id "$id" \
      --l "$out_label" \
      --surf fsaverage "$he"
  done

  # 5.2 Prepare label args as an array
  label_args=()
  for id in "${sens_index[@]}"; do
    label_args+=(--l "$output_dir/fsaverage/label/${he}.num${id}.label")
  done

  # 5.3 Remove existing annot in SUBJECTS_DIR to avoid "already exists"
  rm -f "$SUBJECTS_DIR/fsaverage/label/${he}.aal_sensory.annot" \
        "$SUBJECTS_DIR/fsaverage/label/${he}.aal_sensory.ctab" 2>/dev/null || true

  # 5.4 Build annot (writes into $SUBJECTS_DIR/fsaverage/label/)
  mris_label2annot --s fsaverage --h "$he" \
    --ctab "$input_dir/${he}_aal_sensory_RGB_note.txt" \
    --a aal_sensory \
    "${label_args[@]}"

  # 5.5 Copy annot (and ctab if exists) to output_dir
  src_annot="$SUBJECTS_DIR/fsaverage/label/${he}.aal_sensory.annot"
  dst_dir="$output_dir/fsaverage/label"

  [ -r "$src_annot" ] || { echo "ERROR: annot not generated/readable: $src_annot"; exit 1; }
  cp -f "$src_annot" "$dst_dir/"

  src_ctab="$SUBJECTS_DIR/fsaverage/label/${he}.aal_sensory.ctab"
  if [ -f "$src_ctab" ]; then
    cp -f "$src_ctab" "$dst_dir/"
  fi
done

echo "Done."
echo "Copied annot outputs:"
echo "  $output_dir/fsaverage/label/lh.aal_sensory.annot"
echo "  $output_dir/fsaverage/label/rh.aal_sensory.annot"

