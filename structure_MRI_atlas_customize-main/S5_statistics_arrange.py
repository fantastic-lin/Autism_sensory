#!/usr/bin/env python3
# coding: utf-8

import os
import pandas as pd
import numpy as np

aspect_columns_surface_area = [
    "Frontal_Sup_Orb_L_5","Frontal_Sup_Orb_R_6","Frontal_Mid_Orb_L_9","Frontal_Mid_Orb_R_10",
    "Frontal_Inf_Orb_L_15","Frontal_Inf_Orb_R_16","Olfactory_L_21","Olfactory_R_22",
    "Frontal_Mid_Orb_L_25","Frontal_Mid_Orb_R_26","Insula_L_29","Insula_R_30"
]

aspect_columns_gray_matter_volume = [
    "Frontal_Sup_Orb_L_5","Frontal_Sup_Orb_R_6","Frontal_Mid_Orb_L_9","Frontal_Mid_Orb_R_10",
    "Frontal_Inf_Orb_L_15","Frontal_Inf_Orb_R_16","Olfactory_L_21","Olfactory_R_22",
    "Frontal_Mid_Orb_L_25","Frontal_Mid_Orb_R_26","Insula_L_29","Insula_R_30",
    "TotalGrayVol","Left_Hippocampus_vol","Right_Hippocampus_vol","Left_Amygdala_vol","Right_Amygdala_vol"
]

SUBJECTS_DIR = "/histor/sun/linlin/4_olfactory/ABIDE/abide-master/sMRI/sMRI_git/FS_git"

# ------------------------------------------------------------
# Subjects
# ------------------------------------------------------------
subdirs = [
    name for name in os.listdir(SUBJECTS_DIR)
    if os.path.isdir(os.path.join(SUBJECTS_DIR, name)) and name != "fsaverage"
]

# ------------------------------------------------------------
# Robust readers for the summary CSV files
# ------------------------------------------------------------
def load_surface_area_dict(csv_path: str) -> dict:
    """
    Accepts either:
      - with header: subject,total_surface_area
      - without header: subj_ID,aspect
    Returns {subj_ID: total_surface_area_float}
    """
    df = pd.read_csv(csv_path)
    if set(df.columns) >= {"subject", "total_surface_area"}:
        return dict(zip(df["subject"].astype(str), pd.to_numeric(df["total_surface_area"], errors="coerce")))
    # fallback to headerless format
    df = pd.read_csv(csv_path, header=None)
    df.columns = ["subj_ID", "aspect"]
    return dict(zip(df["subj_ID"].astype(str), pd.to_numeric(df["aspect"], errors="coerce")))

def load_volume_dict(csv_path: str) -> dict:
    """
    Accepts either:
      - with header: subject,TotalGrayVol,Left_Hippocampus,Right_Hippocampus,Left_Amygdala,Right_Amygdala
      - without header (your older format)
    Returns {subj_ID: {TotalGrayVol:..., Left_Hippocampus_vol:..., ...}}
    """
    df = pd.read_csv(csv_path)
    if "subject" in df.columns:
        # Map possible header names to your desired names
        rename_map = {
            "Left_Hippocampus": "Left_Hippocampus_vol",
            "Right_Hippocampus": "Right_Hippocampus_vol",
            "Left_Amygdala": "Left_Amygdala_vol",
            "Right_Amygdala": "Right_Amygdala_vol",
        }
        df = df.rename(columns=rename_map)
        keep_cols = ["subject", "TotalGrayVol", "Left_Hippocampus_vol", "Right_Hippocampus_vol", "Left_Amygdala_vol", "Right_Amygdala_vol"]
        for c in keep_cols:
            if c not in df.columns:
                df[c] = np.nan
        df[keep_cols[1:]] = df[keep_cols[1:]].apply(pd.to_numeric, errors="coerce")
        return df.set_index("subject")[keep_cols[1:]].to_dict(orient="index")

    # fallback to headerless format
    df = pd.read_csv(csv_path, header=None)
    df.columns = ["subj_ID","TotalGrayVol","Left_Hippocampus_vol","Right_Hippocampus_vol","Left_Amygdala_vol","Right_Amygdala_vol"]
    df[["TotalGrayVol","Left_Hippocampus_vol","Right_Hippocampus_vol","Left_Amygdala_vol","Right_Amygdala_vol"]] = \
        df[["TotalGrayVol","Left_Hippocampus_vol","Right_Hippocampus_vol","Left_Amygdala_vol","Right_Amygdala_vol"]].apply(pd.to_numeric, errors="coerce")
    return df.set_index("subj_ID").to_dict(orient="index")

# ------------------------------------------------------------
# Read subject ROI stats file -> dict
# ------------------------------------------------------------
def read_roi_file(file_path: str, value_col_index: int) -> dict:
    """
    file looks like: key col1 col2 col3...
    We take col0 as key, and column `value_col_index` as value.
    Uses whitespace splitting robustly.
    """
    df = pd.read_csv(file_path, sep=r"\s+", engine="python", header=None)
    if df.shape[1] <= value_col_index:
        raise ValueError(f"{file_path} has only {df.shape[1]} columns; cannot take col {value_col_index}")
    keys = df.iloc[:, 0].astype(str)
    vals = pd.to_numeric(df.iloc[:, value_col_index], errors="coerce")
    return dict(zip(keys, vals))

def build_subject_table(value_col_index: int, expected_cols: list, extra_dict: dict, extra_key_name: str) -> tuple[pd.DataFrame, list]:
    """
    value_col_index: 2 (surface area) or 3 (volume) per your assumption
    expected_cols: columns you want in output
    extra_dict:
      - surface: {subj: total_surface_area}
      - volume: {subj: {TotalGrayVol:..., ...}}
    extra_key_name:
      - "Total_surface_area" or "" (volume dict already has keys)
    Returns: (df, missing_subjects)
    """
    rows = []
    missing = []

    for subj in subdirs:
        lh_path = os.path.join(SUBJECTS_DIR, subj, "stats", f"{subj}_lh_1.txt")
        rh_path = os.path.join(SUBJECTS_DIR, subj, "stats", f"{subj}_rh_1.txt")

        if not (os.path.exists(lh_path) and os.path.exists(rh_path)):
            missing.append(subj)
            continue

        lh = read_roi_file(lh_path, value_col_index)
        rh = read_roi_file(rh_path, value_col_index)

        merged = {}
        merged.update(lh)
        merged.update(rh)

        # add extra measures
        if extra_key_name:  # surface
            merged[extra_key_name] = extra_dict.get(subj, np.nan)
        else:  # volume extras are dict-of-dicts
            extras = extra_dict.get(subj, {})
            if isinstance(extras, dict):
                merged.update(extras)

        # enforce expected columns (others ignored)
        row = {"Subject_ID": subj}
        for c in expected_cols:
            row[c] = merged.get(c, np.nan)
        rows.append(row)

    df = pd.DataFrame(rows)
    return df, missing

# ------------------------------------------------------------
# (1) Surface area DF
# ------------------------------------------------------------
total_surface_area_dict = load_surface_area_dict(
    "/histor/sun/linlin/4_olfactory/ABIDE/abide-master/sMRI/sMRI_git/FS_git_surface_area.csv"
)

surface_area_df, missing_surface = build_subject_table(
    value_col_index=2,
    expected_cols=aspect_columns_surface_area + ["Total_surface_area"],
    extra_dict=total_surface_area_dict,
    extra_key_name="Total_surface_area",
)

surface_area_df.to_csv("FS_git_surface_area_df.csv", index=False)
print(f"[Surface] saved FS_git_surface_area_df.csv, n={len(surface_area_df)}, missing={len(missing_surface)}")

# ------------------------------------------------------------
# (2) Gray matter volume DF
# ------------------------------------------------------------
gray_matter_volume_dict = load_volume_dict(
    "/histor/sun/linlin/4_olfactory/ABIDE/abide-master/sMRI/sMRI_git/FS_git_gray_matter_volume.csv"
)

gray_matter_volume_df, missing_vol = build_subject_table(
    value_col_index=3,
    expected_cols=aspect_columns_gray_matter_volume,
    extra_dict=gray_matter_volume_dict,
    extra_key_name="",  # extras already include the keys
)

gray_matter_volume_df.to_csv("FS_git_gray_matter_volume_df.csv", index=False)
print(f"[Volume] saved FS_git_gray_matter_volume_df.csv, n={len(gray_matter_volume_df)}, missing={len(missing_vol)}")