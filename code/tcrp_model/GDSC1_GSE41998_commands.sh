set -ex
python3 /root/capsule/code/tcrp_model/model/MAML_DRUG.py --dataset GDSC1_GSE41998 --tissue GSE41998_Breast --drug Paclitaxel --K 10 --num_trials 20 --tissue_num 12 --meta_batch_size 10 --meta_lr 0.01 --inner_lr 0.01 --layer 1 --run_name 210803_drug-baseline-models
