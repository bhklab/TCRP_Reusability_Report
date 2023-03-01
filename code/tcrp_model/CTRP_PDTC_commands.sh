set -ex
python3 /root/capsule/code/tcrp_model/model/MAML_DRUG.py --dataset CTRP_PDTC --tissue PDTC --drug Axitinib --K 10 --num_trials 20 --tissue_num 12 --meta_batch_size 10 --meta_lr 0.01 --inner_lr 0.01 --layer 1 --run_name 210803_drug-baseline-models
python3 /root/capsule/code/tcrp_model/model/MAML_DRUG.py --dataset CTRP_PDTC --tissue PDTC --drug BI-2536 --K 10 --num_trials 20 --tissue_num 20 --meta_batch_size 10 --meta_lr 0.01 --inner_lr 0.01 --layer 2 --run_name 210803_drug-baseline-models
python3 /root/capsule/code/tcrp_model/model/MAML_DRUG.py --dataset CTRP_PDTC --tissue PDTC --drug Bosutinib --K 10 --num_trials 20 --tissue_num 20 --meta_batch_size 10 --meta_lr 0.1 --inner_lr 0.01 --layer 1 --run_name 210803_drug-baseline-models
python3 /root/capsule/code/tcrp_model/model/MAML_DRUG.py --dataset CTRP_PDTC --tissue PDTC --drug KU-55933 --K 10 --num_trials 20 --tissue_num 20 --meta_batch_size 10 --meta_lr 0.01 --inner_lr 0.01 --layer 2 --run_name 210803_drug-baseline-models
python3 /root/capsule/code/tcrp_model/model/MAML_DRUG.py --dataset CTRP_PDTC --tissue PDTC --drug MK-2206 --K 10 --num_trials 20 --tissue_num 20 --meta_batch_size 10 --meta_lr 0.001 --inner_lr 0.001 --layer 2 --run_name 210803_drug-baseline-models
python3 /root/capsule/code/tcrp_model/model/MAML_DRUG.py --dataset CTRP_PDTC --tissue PDTC --drug Paclitaxel --K 10 --num_trials 20 --tissue_num 20 --meta_batch_size 10 --meta_lr 0.01 --inner_lr 0.01 --layer 2 --run_name 210803_drug-baseline-models
python3 /root/capsule/code/tcrp_model/model/MAML_DRUG.py --dataset CTRP_PDTC --tissue PDTC --drug PLX4720 --K 10 --num_trials 20 --tissue_num 6 --meta_batch_size 10 --meta_lr 0.01 --inner_lr 0.1 --layer 2 --run_name 210803_drug-baseline-models


