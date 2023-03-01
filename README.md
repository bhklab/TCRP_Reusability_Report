# Reusability Report: Few-shot learning creates predictive models of drug response that translate from high-throughput screens to individual patients

This is the code to replicate all figures from "Reusability Report: Few-shot learning creates predictive models of drug response that translate from high-throughput screens to individual patients". The goal of this paper is to reproduce the results from the Nature Cancer paper "Few-shot learning creates predictive models of drug response that translate from high-throughput screens to individual patients"[1]

*N.B*: This code is available on Code Ocean for easy, one-click running [here](https://codeocean.com/capsule/8411716/tree/v2)

There are multiple different ways you can interact with this capsule:
# Reproducing results of Reusability Report (within a Code Ocean Capsule)
1. Use the App Panel dropdown on the left sidebar to select a reference-validation combination that you would like to run few-shot TCRP analysis on, as well as baseline comparison against established machine learning methods 
2. Click "Reproducible Run" at the top right of the capsule
3. The results you will receive are the baseline performance results, TCRP performance results and baseline/TCRP performance comparison. (**N.B.** due to the limited compute hours you may have with your Code Ocean compute account, some good examples to run would be `gCSI_PDTC`, `CTRP_PDTC` or `GDSC1_GSE41998` as they are very small in compute runtime)

# Running your own, new clinical context 
1. reference files in `/data/references` and `/data/transfers/` to understand input files that are needed. 
2. Run the correct preprocessing notebook. Substitute dataset names where appropriate Decisions for the appropriate notebook are made as follows:
 a. `/code/preprocessing/extract_features/new_everything_preprocessing.ipynb` if you are substituting **both reference and validation cell line datasets**
 b. `/code/preprocessing/extract_features/new_PDTX_preprocessing.ipynb` if you are susbtituting **validation cell line dataset** 
 c. `/code/preprocessing/extract_features/new_clinical_processing.ipynb` if you are substiting **validation clinical datasets** with R/NR response
3. save drug features under `/data/drug_features` with the corresponding nomenclature of {reference}_{validation}

# Performance 
1. Navigate to `/code/tcrp_model/pipelines/prepare_complete_run.py` and replace lines 33 and 34 wiuth the correct drug feature/tissue name
2. A folder will be created under `/code/tcrp_model/created_models` with the name `created_models_{reference}_{validation}`
3. Navigate to `/code/tcrp_model/created_models/created_models_{reference}_{validation}/baseline_cmd/subcommands` and execute ALL bash scripts in the directory
4. Navigate to `/code/tcrp_model/created_models/created_models_{reference}_{validation}/MAML_cmd/subcommands` and execute ALL bash scripts in the directory
5. For each drug that TCRP was run on, use the notebook `/code/tcrp_model/model/find_max_run.ipynb` to find the most optimal run (substituting tissue/dataset/drug name appropriately). This notebook will tell you which output file has the best TCRP correlation, which you can search for in the corresponding subcommands script to find the full command.
6. For each drug, rerun the optimal TCRP command found in step 5 to replace the current result with the optimized result. 
7. Run `/code/tcrp_model/model_comparisons/plot_results.py` with the corresponding flags of `--reference {reference} --validation {validation}`.
8. Your basleine comparison vs. TCRP PNG file will be saved in `/results/`


Should you have any inquiries or questions, pelase contact [emily.so@mail.utoronto.ca](emily.so@mail.utoronto.ca)
# References
1. Ma, J. et al. Few-shot learning creates predictive models of drug response that translate from high-throughput screens to individual patients. Nat Cancer 2, 233–244 (2021).
