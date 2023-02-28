set -ex
python3 /root/capsule/code/tcrp_model/model/MAML_DRUG.py --dataset gCSI_PDTC --tissue PDTC --drug Paclitaxel --K 10 --num_trials 20 --tissue_num 12 --meta_batch_size 10 --meta_lr 0.1 --inner_lr 0.001 --layer 1 --run_name 210803_drug-baseline-models
