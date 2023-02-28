import time
import numpy as np
import random
import os
import sys
import pickle
import copy
import argparse

parser = argparse.ArgumentParser(description="Generate slurm scripts to run experiment")
parser.add_argument('--dataset', default='', help='dataset to run crossvalidation on')
parser.add_argument('--tissue', default='', help='transfer tissue of interest')
parser.add_argument('--drug_list_file', default='selected_drug_list', help='list of drugs to generate code for')
parser.add_argument('--job', default='job', help='Job name for slurm script')
parser.add_argument('--job_id', default='', help='Jod ID to append to create a separate script if needed')
parser.add_argument('--run_name', default='run')

args = parser.parse_args()
drug_list_file, job, job_id = args.drug_list_file, args.job, args.job_id

filepath = os.path.realpath(__file__)
dir_name = os.path.dirname(filepath)
home_dir = os.path.dirname(os.path.dirname(dir_name))
dataset=args.dataset
target_tissue = args.tissue
work_dic = f"/data/drug_features/{dataset}/"
job_directory = f"/code/tcrp_model/created_models/created_models_{dataset}/"

file_handle = open( drug_list_file )

cmd_folder = job_directory + 'MAML_cmd/'
os.system("mkdir -p {}".format(cmd_folder))
cmd_list = []
log_list = []
gene_list = []
for line in file_handle:

	gene = line.rstrip()
	gene_list.append(gene)

	#tissue_list = work_dic + 'sanger_cell_line_data/' + gene + '_tissue_cell_line_list.pkl'
	tissue_list = work_dic + gene + '_tissue_cell_line_list.pkl'

	with open(tissue_list, 'rb') as f:
		tissue_map = pickle.load(f)
	for tissue, tissue_cell_line in tissue_map.items():
			# if len(tissue_cell_line) < 15:
			# 	continue
			if tissue != target_tissue:
				continue
			if len(tissue_cell_line) > 15:
				tissue_num_list = ['6','12','20']
			if len(tissue_cell_line) < 15:
				tissue_num_list = ['6','12','20']
			#cmd_str = '$python ' + '/root/capsule/code/tcrp/pipelines/generate_baseline_job_cv.py' + '--tissue {} --drug {} --K 10 --num_trials 20 --run_name {}'.format(tissue, gene, args.run_name)
			#cmd_list.append(cmd_str)
			
			log_folder = f"/results/{dataset}/run-logs/{gene}/{tissue}/"
			log_list.append(f"mkdir -p {log_folder}")
			# os.system("mkdir -p {}".format(log_folder))

			for meta_lr in ['0.1', '0.01', '0.001']:
				for inner_lr in ['0.1', '0.01', '0.001']:
					for layer in ['1', '2']:
						for tissue_num in  tissue_num_list:
							log_file = log_folder + tissue + '_' + gene + '_' + meta_lr + '_' + inner_lr + '_' + layer + '_' + tissue_num + '.log'
							
							if os.path.exists( log_file ):
								
								normal_finish = 0
								file_handle = open( log_file, 'r')
								for line in file_handle:
									if 'Best loss' in line:
										normal_finish = 1
										break
								if normal_finish == 1:
									continue

							cmd_str = 'python3 ' + home_dir + '/tcrp_model/model/' + 'MAML_DRUG.py --dataset ' + dataset+ ' --tissue ' + tissue + ' --drug ' + gene + ' --K 10 --num_trials 20' + ' --tissue_num ' + tissue_num + ' --meta_batch_size 10 --meta_lr ' + meta_lr + ' --inner_lr ' + inner_lr + ' --layer ' + layer + ' --run_name ' + args.run_name + ' > ' + log_file
							cmd_list.append( cmd_str )						

					cmd_list.append('')

file_handle.close()


subcommand_directory = cmd_folder + "subcommands"
os.system("mkdir -p {}".format(subcommand_directory))
with open(subcommand_directory + '/' + 'subcommands_MAML_{}{}.sh'.format(job, job_id), 'w') as f:
    f.write('#!/bin/bash\n')
    f.write('set -ex\n')
    f.write(f"mkdir -p /results/{dataset}\n")
    f.write(f"mkdir -p /results/{dataset}/run-logs/\n")
    for i in gene_list:
        f.write(f"mkdir -p /results/{dataset}/run-logs/{i}\n")
    f.writelines('\n'.join(log_list) + '\n')
    f.writelines('\n'.join(cmd_list) + '\n')


tempfile = cmd_folder + 'run_MAML_{}{}.sh'.format(job, job_id)

# slurm_output = job_directory + "slurm-logs"
# os.system("mkdir -p {}".format(slurm_output))

# file_handle = open(tempfile,'w')

# file_handle.writelines('#!/bin/bash\n')
# file_handle.writelines("#SBATCH --job-name {}{}\n".format(job, job_id))
# file_handle.writelines("#SBATCH --output={}/{}{}.%j\n".format(slurm_output, job, job_id))
# file_handle.writelines("#SBATCH --cpus-per-task=4\n")
# file_handle.writelines("#SBATCH --mem=64G\n")
# file_handle.writelines("#SBATCH --partition=nrnb-gpu\n")
# file_handle.writelines("#SBATCH --account=nrnb-gpu\n")
# file_handle.writelines("#SBATCH --gres=gpu:1\n\n")

# file_handle.writelines("/usr/bin/bash {}/subcommands_MAML_{}{}.sh\n".format(subcommand_directory, job, job_id))

# file_handle.close()

