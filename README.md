# Autism_sensory

This repository contains code used in research on autism perception and sensory abnormalities, including:

- A custom structural magnetic resonance imaging (sMRI) brain-region segmentation pipeline (FreeSurfer + self-defined AAL atlas)
- A Bayesian analytic framework (TADA) for prioritizing sensory-related autism risk genes

---

## 1. Extracting Structural MRI Metrics from FreeSurfer Using a Self-Defined AAL Brain Atlas

### Part 1: Environment Setup

#### Prerequisites
Ensure **FreeSurfer** is correctly installed and initialized.

#### 1) Set Up FreeSurfer
Add the following to `~/.bashrc`:
```bash
export FREESURFER_HOME=/path/to/your/freesurfer
source $FREESURFER_HOME/SetUpFreeSurfer.sh
```

Then apply:
```bash
source ~/.bashrc
```

#### 2) Set Up Work Directory
Replace operation_dir="/histor/sun/linlin/4_olfactory/ABIDE/abide-master/sMRI_codeocean" with your own work directory path in config.env
---

### Part 2: Analysis Pipeline

#### Run the entire process (S1 → S5)
```bash
bash run_all.sh
```

---

#### Scripts (S1–S5): Functions, Inputs, and Outputs

##### (1) Construct AAL Template
- **Script:** `S1_create_AAL_annot.sh`
- **Inputs:**
  - `./subject_test`
  - `./input/fsaverage`
  - `./input/*h_aal_sensory_RGB_note.txt`
  - `./input/AAL116_1mm.nii`
- **Outputs:**
  - `./output/fsaverage/mgh/*h.aal3.mgh`
  - `./output/fsaverage/label/*h.num{ID}.label`
  - `./output/fsaverage/label/*h_aal_sensory.annot`

---

##### (2) Extract Structural Measurement Values
- **Script:** `S2_ROI_analysis.sh`
- **Inputs:**
  - `./subject_test`
  - `./input/S2345_input_ABIDE1_qc_OK.json`
  - `./output/fsaverage/label/*h_aal_sensory.annot`
- **Outputs:**
  - `./output/S2345_output_stats_result/Caltech/{subjectid}/label/*h.aal_sens_{subjectid}.annot`
  - `./output/S2345_output_stats_result/Caltech/{subjectid}/stats/{subjectid}_*h.txt`

---

##### (3) Simplify Structural Measurement Files
- **Script:** `S3_AAL_olfactory_extract_simplify.sh`
- **Inputs:**
  - `./input/S2345_input_ABIDE1_qc_OK.json`
  - `./output/S2345_output_stats_result/Caltech/{subjectid}/stats/{subjectid}_*h.txt`
- **Output:**
  - `./output/S2345_output_stats_result/Caltech/{subjectid}/stats/{subjectid}_*h_1.txt`

---

##### (4) Extract Default FreeSurfer Measurements
- **Script:** `S4_extract_freesurfer_value.sh`
- **Inputs:**
  - `./subject_test`
  - `./input/S2345_input_ABIDE1_qc_OK.json`
- **Outputs:**
  - `./output/S2345_output_stats_result/ABIDE1_surface_area.csv`
  - `./output/S2345_output_stats_result/ABIDE1_gray_matter_volume.csv`

---

##### (5) Summarize Values into One DataFrame
- **Script:** `S5_statistics_arrange.py`
- **Inputs:**
  - `./output/S2345_output_stats_result/Caltech/{subjectid}/stats/{subjectid}_*h_1.txt`
  - `./output/S2345_output_stats_result/ABIDE1_surface_area.csv`
  - `./output/S2345_output_stats_result/ABIDE1_gray_matter_volume.csv`
- **Outputs:**
  - `./output/S2345_output_stats_result/ABIDE1_surface_are_df.csv`
  - `./output/S2345_output_stats_result/ABIDE1_gray_matter_volume_df.csv`

---

### File Notes
- **AAL Atlas:** `AAL116_1mm.nii`
- **RGB annotation:**
  - `lh_aal_olfactory_RGB_note.txt`
  - `rh_aal_olfactory_RGB_note.txt`
- **subject_test:**
  - The original ABIDE data were preprocessed with **fMRIPrep**.

---

## 2. Bayesian Analytic Framework (TADA) for Autism Risk Prioritization of Sensory-Related Genes

- **Script:** `run_TADA.R`

### Inputs
1) **Prior mutation information**
   - `./TADA_sup_info_mut_hgnc.csv`

2) **Mutation events in the population**
   - `./dn_SNV_info_hgnc.csv`
   - `./in_SNV_info_hgnc.csv`
   - `./cc_SNV_info_hgnc.csv`

3) **Sensory-related gene information**
   - `./genecards_sensory_hgnc_supple.csv`

4) **CNV information**
   - `./CNV_merge_note.csv`

5) **ASD risk genes from SFARI**
   - `./SFARI_ASD_gene_note.csv`

6) **Sensory type annotations for sensory-related genes**
   - `./genecards_hgnc_sensory_info.csv`

### Outputs
1) **Bayes factors (BFs) calculated by TADA**
   - `./SNV_SRG_BF.csv`
   - `./CNV_SRG_BF.csv`
   - `./SNV_CNV_SRG_BF.csv`

2) **Conversion of BF to FDR**
   - `./SNV_CNV_SRG_BF_FDR_pval.csv`
   - `./SNV_CNV_SRG_BF_FDR_pval_ASD.csv`
   - `./SNV_CNV_SRG_BF_FDR_pval_ASD_note.csv`
   - `./BF_sensory.csv`
