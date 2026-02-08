# Autism_sensory
Key codes used in autism perception and sensory abnormalities research include custom sMRI brain region segmentation atlases and TADA analysis.

### 1. Extracting Structural MRI Metrics from FreeSurfer Using a Self-Defined AAL Brain Atlas


## Overview
This repository provides a complete pipeline for processing structural MRI (sMRI) data using FreeSurfer and a custom AAL (Automated Anatomical Labeling) brain atlas.

## Part 1: Environment Setup
### Prerequisites
Ensure FreeSurfer is correctly installed and initialized.

### 1. Set Up FreeSurfer
Add to ~/.bashrc:
```
export FREESURFER_HOME=/path/to/your/freesurfer
source $FREESURFER_HOME/SetUpFreeSurfer.sh
```
Then apply:
```
source ~/.bashrc
```

## Part 2: Analysis Pipeline
### (1) Download Subject Files from ABIDE
- Script: S0_sMRI_1_data_download_git.py
- Inputs: Phenotypic_V1_0b_preprocessed1.csv, each_sub_necessary_files.txt
- Output: FS_successful_download.json

### (2) Construct AAL Template
- Script: S1_create_AAL_annot.sh
- Inputs: Subject dictionary, *h_aal_olfactory_RGB_note.txt
- Outputs: *h.aal3.mgh, *h.num{ID}.label, *h_aal_olfactory.annot

### (3) Extract Structural Measurement Values
- Script: S2_ROI_analysis.sh
- Inputs: FS_successful_download.json, *h.aal_olfactory.annot, *h.aal_olf_{subjectid}.annot
- Output: ${subjectid}*h.txt

### (4) Simplify Structural Measurement Files
- Script: S3_AAL_olfactory_extract_simplify.sh
- Inputs: FS_successful_download.json, ${subjectid}*h.txt
- Output: ${subjectid}_*h_1.txt

### (5) Extract Default FreeSurfer Measurements
- Script: S4_extract_freesurfer_value.sh
- Output: FS_git_surface_area.csv, FS_git_gray_matter_volume.csv

### (6) Summarize Values into One DataFrame
- Script: S5_statistics_arrange.py
- Outputs: FS_git_surface_area_df.csv, FS_git_gray_matter_volume_df.csv

## File Notes
- AAL Atlas: AAL116_1mm.nii
- Required subject files: each_sub_necessary_files.txt
- Subject info: Phenotypic_V1_0b_preprocessed1.csv
- RGB annotation: lh_aal_olfactory_RGB_note.txt, rh_aal_olfactory_RGB_note.txt

### 2. Bayesian analytic framework, TADA, for the Autism risk prioritization of sensory-related genes
- Script: run_TADA.R
- Inputs: 
    1) prior mutation info: ./TADA_sup_info_mut_hgnc.csv
    2) Mutation event in Population:  ./dn_SNV_info_hgnc.csv, ./in_SNV_info_hgnc.csv, ./cc_SNV_info_hgnc.csv
    3) Sensory-related gene information: ./genecards_sensory_hgnc_supple.csv
    4) CNV information:./CNV_merge_note.csv
    5) ASD risk gene in SFARI: ./SFARI_ASD_gene_note.csv
    6) Sensory type note of sensory-related genes: ./genecards_hgnc_sensory_info.csv
- Output: 
    1) BF calculated by TADA:./SNV_ORG_BF.csv, ./CNV_ORG_BF.csv, ./SNV_CNV_ORG_BF.csv
    2) Convert BF to FDR: ./SNV_CNV_ORG_BF_FDR_pval.csv, ./SNV_CNV_ORG_BF_FDR_pval_ASD.csv, ./SNV_CNV_ORG_BF_FDR_pval_ASD_note.csv, ./BF_sensory.csv

- Script: run_TADA.R
- Inputs:
    1) Prior mutation information: ./TADA_sup_info_mut_hgnc.csv
    2) Mutation events in the population:
        • ./dn_SNV_info_hgnc.csv
        • ./in_SNV_info_hgnc.csv
        • ./cc_SNV_info_hgnc.csv
    3) Sensory-related gene information: ./genecards_sensory_hgnc_supple.csv
    4) CNV information: ./CNV_merge_note.csv
    5) ASD risk genes from SFARI: ./SFARI_ASD_gene_note.csv
    6) Sensory type annotations for sensory-related genes: ./genecards_hgnc_sensory_info.csv
- Outputs:
    1) Bayes factors (BFs) calculated by TADA:
        • ./SNV_SRG_BF.csv
        • ./CNV_SRG_BF.csv
        • ./SNV_CNV_SRG_BF.csv
    2) Conversion of BF to FDR:
        • ./SNV_CNV_SRG_BF_FDR_pval.csv
        • ./SNV_CNV_SRG_BF_FDR_pval_ASD.csv
        • ./SNV_CNV_SRG_BF_FDR_pval_ASD_note.csv
        • ./BF_sensory.csv