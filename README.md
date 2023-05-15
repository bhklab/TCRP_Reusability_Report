# Reusability Report: Few-shot learning creates predictive models of drug response that translate from high-throughput screens to individual patients (Version 2.3)

This is the code to replicate all figures from "Reusability Report: Few-shot learning creates predictive models of drug response that translate from high-throughput screens to individual patients". The goal of this paper is to reproduce the results from the Nature Cancer paper "Few-shot learning creates predictive models of drug response that translate from high-throughput screens to individual patients"[1]

There are multiple different ways you can interact with this capsule 
# Reproducing results of Reusability Report 
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
1. Navigate to `/code/tcrp_model/pipelines/prepare_complete_run.py` and replace lines 33 and 34 with the correct drug feature/tissue name
2. A folder will be created under `/code/tcrp_model/created_models` with the name `created_models_{reference}_{validation}`

3. Navigate to `/code/tcrp_model/created_models/created_models_{reference}_{validation}/baseline_cmd/subcommands` and execute ALL bash scripts in the directory (the order does not matter)
- Example command: `python3 /root/capsule/code/tcrp_model/model/MAML_DRUG.py --dataset CTRP_PDTC --tissue PDTC --drug Axitinib --K 10 --num_trials 20 --tissue_num 12 --meta_batch_size 10 --meta_lr 0.01 --inner_lr 0.01 --layer 1 --run_name 210803_drug-baseline-models`
4. Navigate to `/code/tcrp_model/created_models/created_models_{reference}_{validation}/MAML_cmd/subcommands` and execute ALL bash scripts in the directory (the order does not matter)
- **N.B** Drugs whose performance can be visualized are drugs that completed both baseline and TCRP execution
5. Use the notebook `/code/tcrp_model/model/find_max_run.ipynb` to find the most optimal run (substituting tissue/dataset/drug name appropriately). This notebook will tell you which output file has the best TCRP correlation, which you can search for in the corresponding subcommands script to find the full command.
6. If, alternatively, you would like to calculate additional performance metrics for measuring TCRP in comparison to baseline results (the default is Pearson's correlation), you can visit `/code/tcrp_model/model/new_score.ipynb` to view how that information can be generated. The new arrays containing the new perforamnce metric should then be saved to replace the original `TCRP_performance.npz` for each dataset-drug-tissue.
6. For each drug, rerun the optimal TCRP command found in step 5 to replace the current result with the optimized result. 
7. Run `/code/tcrp_model/model_comparisons/plot_results.py` with the corresponding flags of `--reference {reference} --validation {validation}`.
8. Your baseline comparison vs. TCRP PNG file will be saved in `/results/`


Should you have any inquiries or questions, pelase contact [emily.so@mail.utoronto.ca](emily.so@mail.utoronto.ca)
# References
1. Ma, J. et al. Few-shot learning creates predictive models of drug response that translate from high-throughput screens to individual patients. Nat Cancer 2, 233–244 (2021).
