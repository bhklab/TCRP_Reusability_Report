#!/bin/bash
set -ex
python /root/capsule/code/tcrp_model/pipelines/generate_fewshot_samples.py --dataset GDSC1_PDTC --tissue PDTC --drug PD-0332991 --K 10 --num_trials 20 --run_name 210803_drug-baseline-models
python /code/tcrp_model/baselines/supp_baseline_DRUG.py --dataset GDSC1_PDTC --tissue PDTC --drug PD-0332991 --K 10 --num_trials 20 --run_name 210803_drug-baseline-models --fewshot_data_path /data/fewshot_data/fewshot_data_GDSC1_PDTC
python /root/capsule/code/tcrp_model/pipelines/generate_fewshot_samples.py --dataset GDSC1_PDTC --tissue PDTC --drug GW_441756 --K 10 --num_trials 20 --run_name 210803_drug-baseline-models
python /code/tcrp_model/baselines/supp_baseline_DRUG.py --dataset GDSC1_PDTC --tissue PDTC --drug GW_441756 --K 10 --num_trials 20 --run_name 210803_drug-baseline-models --fewshot_data_path /data/fewshot_data/fewshot_data_GDSC1_PDTC
