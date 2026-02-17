import json
from pathlib import Path

import pandas as pd
from dotenv import dotenv_values

# ============================================================
# 0) Paths from config.env
# ============================================================
SCRIPT_DIR = Path(__file__).resolve().parent
cfg = dotenv_values(SCRIPT_DIR / "config.env")

operation_dir = cfg.get("operation_dir")
if not operation_dir:
    raise ValueError("operation_dir is not set in config.env")

operation_dir = operation_dir.rstrip("/")
input_dir = Path(operation_dir) / "input"
output_dir = Path(operation_dir) / "output"

stats_store_dir = output_dir / "S2345_output_stats_result"
json_file = input_dir / "S2345_input_ABIDE1_qc_OK.json"

print(f"operation_dir={operation_dir}")
print(f"input_dir={input_dir}/")
print(f"output_dir={output_dir}/")

if not json_file.exists():
    raise FileNotFoundError(f"Missing JSON: {json_file}")

version = "ABIDE1"

# ============================================================
# 1) Column definitions
# ============================================================
aspect_columns_surface_area = [
    "Total_surface_area",
    "Calcarine_L_43", "Calcarine_R_44",
    "Cingulum_Ant_L_31", "Cingulum_Ant_R_32",
    "Cingulum_Mid_L_33", "Cingulum_Mid_R_34",
    "Cingulum_Post_L_35", "Cingulum_Post_R_36",
    "Cuneus_L_45", "Cuneus_R_46",
    "Frontal_Inf_Orb_L_15", "Frontal_Inf_Orb_R_16",
    "Frontal_Med_Orb_L_25", "Frontal_Med_Orb_R_26",
    "Frontal_Mid_Orb_L_9", "Frontal_Mid_Orb_R_10",
    "Frontal_Sup_Orb_L_5", "Frontal_Sup_Orb_R_6",
    "Fusiform_L_55", "Fusiform_R_56",
    "Heschl_L_79", "Heschl_R_80",
    "Insula_L_29", "Insula_R_30",
    "Lingual_L_47", "Lingual_R_48",
    "Occipital_Inf_L_53", "Occipital_Inf_R_54",
    "Occipital_Mid_L_51", "Occipital_Mid_R_52",
    "Occipital_Sup_L_49", "Occipital_Sup_R_50",
    "Olfactory_L_21", "Olfactory_R_22",
    "Parietal_Inf_L_61", "Parietal_Inf_R_62",
    "Parietal_Sup_L_59", "Parietal_Sup_R_60",
    "Postcentral_L_57", "Postcentral_R_58",
    "Precuneus_L_67", "Precuneus_R_68",
    "Rolandic_Oper_L_17", "Rolandic_Oper_R_18",
    "Temporal_Pole_Sup_L_83", "Temporal_Pole_Sup_R_84",
    "Temporal_Sup_L_81", "Temporal_Sup_R_82",
]

aspect_columns_gray_matter_volume = [
    "Calcarine_L_43", "Calcarine_R_44",
    "Cingulum_Ant_L_31", "Cingulum_Ant_R_32",
    "Cingulum_Mid_L_33", "Cingulum_Mid_R_34",
    "Cingulum_Post_L_35", "Cingulum_Post_R_36",
    "Cuneus_L_45", "Cuneus_R_46",
    "Frontal_Inf_Orb_L_15", "Frontal_Inf_Orb_R_16",
    "Frontal_Med_Orb_L_25", "Frontal_Med_Orb_R_26",
    "Frontal_Mid_Orb_L_9", "Frontal_Mid_Orb_R_10",
    "Frontal_Sup_Orb_L_5", "Frontal_Sup_Orb_R_6",
    "Fusiform_L_55", "Fusiform_R_56",
    "Heschl_L_79", "Heschl_R_80",
    "Insula_L_29", "Insula_R_30",
    "Lingual_L_47", "Lingual_R_48",
    "Occipital_Inf_L_53", "Occipital_Inf_R_54",
    "Occipital_Mid_L_51", "Occipital_Mid_R_52",
    "Occipital_Sup_L_49", "Occipital_Sup_R_50",
    "Olfactory_L_21", "Olfactory_R_22",
    "Parietal_Inf_L_61", "Parietal_Inf_R_62",
    "Parietal_Sup_L_59", "Parietal_Sup_R_60",
    "Postcentral_L_57", "Postcentral_R_58",
    "Precuneus_L_67", "Precuneus_R_68",
    "Rolandic_Oper_L_17", "Rolandic_Oper_R_18",
    "Temporal_Pole_Sup_L_83", "Temporal_Pole_Sup_R_84",
    "Temporal_Sup_L_81", "Temporal_Sup_R_82",
    "TotalGrayVol",
    "Left_Hippocampus_vol", "Right_Hippocampus_vol",
    "Left_Amygdala_vol", "Right_Amygdala_vol",
    "Left_Thalamus_vol", "Right_Thalamus_vol",
]

# ============================================================
# 2) Load subject list
# ============================================================
with json_file.open("r", encoding="utf-8") as f:
    subject_list = json.load(f)

# ============================================================
# 3) Helpers
# ============================================================
def regression_dict_surface_area(filename: Path) -> dict[str, float]:
    df = pd.read_csv(filename, header=None, names=["subj_ID", "aspect"])
    df["aspect"] = pd.to_numeric(df["aspect"], errors="coerce")
    return df.set_index("subj_ID")["aspect"].to_dict()


def regression_dict_vol(filename: Path) -> dict[str, dict]:
    cols = [
        "subj_ID",
        "TotalGrayVol",
        "Left_Hippocampus_vol", "Right_Hippocampus_vol",
        "Left_Amygdala_vol", "Right_Amygdala_vol",
        "Left_Thalamus_vol", "Right_Thalamus_vol",
    ]
    df = pd.read_csv(filename, header=None, names=cols)
    for c in cols[1:]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df.set_index("subj_ID").to_dict(orient="index")


def subject_df_generate(file_path: Path, value_col_index: int, subject_id: str) -> pd.DataFrame:
    df = pd.read_csv(file_path, sep=r"\s+", engine="python", header=None, usecols=[0, value_col_index])
    df.columns = ["key", "value"]
    return pd.DataFrame({subject_id: df.set_index("key")["value"].to_dict()})


subject_without_htxt: list[str] = []


def subj_statistic_df_generated(
    value_col_index: int,
    aspect_columns: list[str],
    aspect_dict: dict,
    mode: str,
) -> pd.DataFrame:
    rows = {}

    for p in map(Path, subject_list):
        site_name = p.parent.name
        subject_name = p.name.replace("_total", "")
        subject_whole_name = f"{version}_{site_name}_{subject_name}"

        file_lh = stats_store_dir / site_name / subject_name / "stats" / f"{subject_name}_lh_1.txt"
        file_rh = stats_store_dir / site_name / subject_name / "stats" / f"{subject_name}_rh_1.txt"

        if file_lh.exists() and file_rh.exists():
            df_lh = subject_df_generate(file_lh, value_col_index, subject_whole_name)
            df_rh = subject_df_generate(file_rh, value_col_index, subject_whole_name)
            merged = pd.concat([df_lh, df_rh], axis=0)

            subj_map = {k: v[subject_whole_name] for k, v in merged.iterrows()}

            if mode == "surface_area":
                subj_map["Total_surface_area"] = float(aspect_dict.get(subject_whole_name, float("nan")))
            elif mode == "gray_volume":
                extra = aspect_dict.get(subject_whole_name, {})
                if isinstance(extra, dict):
                    subj_map.update(extra)
            else:
                raise ValueError(f"Unknown mode: {mode}")

            rows[subject_whole_name] = subj_map
        else:
            subject_without_htxt.append(subject_whole_name)

    df_out = pd.DataFrame.from_dict(rows, orient="index")

    for c in aspect_columns:
        if c not in df_out.columns:
            df_out[c] = pd.NA

    df_out = df_out[aspect_columns].copy()
    df_out.reset_index(inplace=True)
    df_out.rename(columns={"index": "Subject_ID"}, inplace=True)
    return df_out


# ============================================================
# 4) Run and save
# ============================================================
surface_area_csv = stats_store_dir / f"{version}_surface_area.csv"
total_surface_area_dict = regression_dict_surface_area(surface_area_csv)

surface_area_df = subj_statistic_df_generated(
    value_col_index=2,
    aspect_columns=aspect_columns_surface_area,
    aspect_dict=total_surface_area_dict,
    mode="surface_area",
)
surface_area_df.to_csv(stats_store_dir / f"{version}_surface_area_df.csv", index=False)

gray_volume_csv = stats_store_dir / f"{version}_gray_matter_volume.csv"
gray_matter_volume_dict = regression_dict_vol(gray_volume_csv)

gray_matter_volume_df = subj_statistic_df_generated(
    value_col_index=3,
    aspect_columns=aspect_columns_gray_matter_volume,
    aspect_dict=gray_matter_volume_dict,
    mode="gray_volume",
)
gray_matter_volume_df.to_csv(stats_store_dir / f"{version}_gray_matter_volume_df.csv", index=False)

print("Done.")
if subject_without_htxt:
    print(f"Subjects missing lh/rh _1.txt: {len(subject_without_htxt)}")
print(f"Subjects missing lh/rh _1.txt: {len(subject_without_htxt)}")
