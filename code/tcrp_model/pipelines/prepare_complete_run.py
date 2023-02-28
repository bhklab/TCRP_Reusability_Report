from glob import glob
import os
import pickle
def make_chunks(iterable, nchunks): 
    length = len(iterable)
    chunk_size = (length // nchunks) 
    remainder = length % nchunks

    chunks = []
    i = 0
    while i < length:
        if remainder != 0: 
            end = i + chunk_size + 1
            remainder -= 1
        else: 
            end = i + chunk_size

        chunks.append(iterable[i:end])
        i = end

    return chunks


filepath = os.path.realpath(__file__)
dirname = os.path.dirname(filepath)
home_dir = os.path.dirname(os.path.dirname(dirname))


n_gpus = 20

run_mode = "tcrp" # "tcrp" or "baseline"
run_name = "210803_drug-baseline-models"
dataset = "GDSC1_GSE20194"
tissue = "GSE20194_Breast"
out_directory = f"/code/tcrp_model/created_models/created_models_{dataset}/"
if run_mode == 'baseline':
    task_directory = out_directory + "baseline"
else:
    task_directory = out_directory + "MAML"
os.system("mkdir -p {}".format(task_directory))

drugs = glob(f"/data/drug_features/{dataset}/*.pkl")
drugs = [drug.split('/')[-1].split('_tissue_cell_line_list.pkl')[0] for drug in drugs]

with open('/code/tcrp_model/pipelines/priority_drugs') as f: 
    priority = [i.strip() for i in f.readlines()]

missing = set(priority).difference(drugs)


remaining_drugs = list(set(drugs).difference(priority))
priority = list(set(priority).difference(missing))


if dataset=="GDSC1_PDTC":
    remaining_drugs = list(set(drugs).difference(priority))
    priority = list(set(priority).difference(missing))
    remaining_drugs_chunked = make_chunks(remaining_drugs, n_gpus)
    priority_chunked = make_chunks(priority, n_gpus)
else:
    features_drug = []
    for i in drugs:
        drug_tissue_list = f"/data/drug_features/{dataset}/" + i + '_tissue_cell_line_list.pkl'
        with open(drug_tissue_list, 'rb') as f:             
            drug_tissue_map = pickle.load(f)

        # Load data
        drug_tissue_map = {k: v for k, v in drug_tissue_map.items() if len(v) > 0}
        if (len(drug_tissue_map) > 1) and (len(drug_tissue_map[tissue]) > 10):
            features_drug.append(i)
    priority=features_drug
    priority_chunked = make_chunks(priority, n_gpus)
print("Total number of drugs: {}".format(len(priority)))

fewshot_data_path = f"/data/fewshot_data/fewshot_data_{dataset}"


if dataset=="GDSC1_PDTC":
    for (i,b) in enumerate(zip(priority_chunked)):
        _drugs = b[0]
        drug_input_file = task_directory + '/drugs_input_{}'.format(i)
        with open(drug_input_file, 'w') as f: 
            f.write('\n'.join(_drugs))
        if run_mode == "tcrp": 
            cmd = ['python', '/root/capsule/code/tcrp_model/pipelines/generate_MAML_job_cv.py', '--dataset',dataset,'--tissue',tissue,'--run_name', run_name, '--drug_list_file', drug_input_file, '--job_id', str(i), '--job', 'drugs']
            cmd = ' '.join(cmd)
            os.system(cmd)

        elif run_mode == "baseline": 
            cmd = ['python', '/root/capsule/code/tcrp_model/pipelines/generate_baseline_job_cv.py','--dataset',dataset,'--tissue',tissue, '--run_name', run_name, '--drug_list_file', drug_input_file, '--job_id', str(i), '--job', 'drugs', '--fewshot_data_path', fewshot_data_path]
            cmd = ' '.join(cmd)
            os.system(cmd)
        else:
            break
else:
    for (i,b) in enumerate(zip(priority_chunked)):
        _drugs = b[0]
        drug_input_file = task_directory + '/drugs_input_{}'.format(i)
        with open(drug_input_file, 'w') as f: 
            f.write('\n'.join(_drugs))
        if run_mode == "tcrp": 
            cmd = ['python', '/root/capsule/code/tcrp_model/pipelines/generate_MAML_job_cv.py', '--dataset',dataset,'--tissue',tissue,'--run_name', run_name, '--drug_list_file', drug_input_file, '--job_id', str(i), '--job', 'drugs']
            cmd = ' '.join(cmd)
            os.system(cmd)

        elif run_mode == "baseline": 
            cmd = ['python', '/root/capsule/code/tcrp_model/pipelines/generate_baseline_job_cv.py','--dataset',dataset,'--tissue',tissue, '--run_name', run_name, '--drug_list_file', drug_input_file, '--job_id', str(i), '--job', 'drugs', '--fewshot_data_path', fewshot_data_path]
            cmd = ' '.join(cmd)
            os.system(cmd)
        else:
            break
