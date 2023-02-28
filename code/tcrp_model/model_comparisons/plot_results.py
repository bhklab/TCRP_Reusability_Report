
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import os
from collections import OrderedDict
import pandas as pd
import sys
import argparse

import warnings

def fxn():
    warnings.warn("deprecated", DeprecationWarning)

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    fxn()
    
parser = argparse.ArgumentParser()
parser.add_argument('--reference', type=str, default='', help='Reference dataset used in analysis')
parser.add_argument('--tissue', type=str, default='', help='transfer tissue used in analysis')

args = parser.parse_args()
reference = args.reference
analysis_tissue = args.tissue
results_directory  = f"/results/{reference}_{analysis_tissue}/"

#specify baseline performance results path here 
datapath = Path(f"{results_directory}/baseline_performances/")


results = {}
for outer_directory in datapath.glob("*"): 
    drug = outer_directory.stem
    results[drug] = {}
    
    for inner_directory in outer_directory.glob("*"): 
        tissue = inner_directory.stem
        results[drug][tissue] = {}
        if str(inner_directory).split("/")[-1] != ".ipynb_checkpoints":
            data = np.load(inner_directory / "baseline_performance.npz")
            #data = np.load(outer_directory / "baseline_performance.npz")
            for model in ['Linear Regression', 'Nearest Neighbour', 'Random Forest','Neural Network']:
            #for model in ['Linear Regression', 'Nearest Neighbour','Neural Network']:
                zero = data[f"{model}-zero"]
                zero = np.vstack([zero for _ in range(20)]) # There is only 1 possible zero-shot, so expanding for all trials
                performance = np.mean(np.hstack([zero, data[f"{model}-fewshot"]]), axis=0)

                results[drug][tissue][model] = performance    


#specifiy TCRP performance results path here 
datapath = Path(f"{results_directory}/TCRP_performances/")


for i in os.listdir(f"{results_directory}/TCRP_performances/"):
    drug_path = f"{results_directory}/TCRP_performances/{i}/{tissue}"
    for j in os.listdir(drug_path):
        os.system("mv {} {}/TCRP_performance.npz >/dev/null 2>&1".format(drug_path+'/'+j,drug_path))


for outer_directory in datapath.glob("*"): 
    drug = outer_directory.stem
    for inner_directory in outer_directory.glob("*"): 
        tissue = inner_directory.stem
        if str(inner_directory).split("/")[-1] != ".ipynb_checkpoints":
            data = np.load(inner_directory / "TCRP_performance.npz")

            for model in ['TCRP']: 
                zero = data[f"{model}-zero"]
                new_data = data["TCRP-fewshot"] 
                fewshot = np.vstack([new_data for _ in range(10)])
                zero = np.vstack([zero for _ in range(10)]) # There is only 1 possible zero-shot, so expanding for all trials
                performance = np.mean(np.hstack([zero, fewshot]), axis=0)
                results[drug][tissue][model] = performance    


for i in results:
    remove_key = results[i].pop(".ipynb_checkpoints", None)



results[drug][tissue][model] = performance    


TCRP_results = {'TCRP': []}

for drug, d in results.items(): 
    for tissue, d in d.items(): 
        for model, p in d.items(): 
            if model == "TCRP":
                p = np.nan_to_num(p)
                TCRP_results[model].append(p)

for model, ps in TCRP_results.items(): 
    TCRP_results[model] = np.vstack(ps)


results_by_baseline = {'Linear Regression': [], 'Nearest Neighbour': [], 'Random Forest': [],"TCRP":[],"Neural Network":[]}
#results_by_baseline = {'Linear Regression': [], 'Nearest Neighbour': [],"TCRP":[],"Neural Network":[]}
for drug, d in results.items(): 
    for tissue, d in d.items(): 
        for model, p in d.items(): 
            p = np.nan_to_num(p)
            results_by_baseline[model].append(p)
            
for model, ps in results_by_baseline.items(): 
    results_by_baseline[model] = np.vstack(ps)


import scipy.stats as st
std_list = {"Linear Regression": [],"Nearest Neighbour": [],'Random Forest': [], 'TCRP': [], 'Neural Network': []}
for model, ps in results_by_baseline.items(): 
    min_ci = []
    max_ci = []
    for i in range(ps.shape[1]):
        data = ps[:,i]
        ci = st.t.interval(0.5, len(data)-1, loc=np.mean(data), scale=st.sem(data))
        min_ci.append(ci[0])
        max_ci.append(ci[1])
    std_list[model].extend([np.array(min_ci),np.array(max_ci)])


fig, ax = plt.subplots()
fig.set_size_inches(10,10)

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
color_dict={"Linear Regression":"green","Neural Network":"purple","Random Forest":"teal","Nearest Neighbour": "blue","TCRP":"red"}
for model, ps in results_by_baseline.items(): 
    p = ax.plot(np.arange(11), np.mean(ps, axis=0),label=model,marker='o',color=color_dict[model])
    if (reference=="GDSC1") and (analysis_tissue=="PDTC"):
        ax.vlines(np.arange(11),std_list[model][0],std_list[model][1],color=color_dict[model])
ax.legend(loc='lower center', bbox_to_anchor=(0.5, 1.05),
          ncol=3, fancybox=True, shadow=True)
labels = ['Pretrained'] + [str(i) for i in range(1, 11)]
ax.set_xticks(np.arange(11))
ax.set_xticklabels(labels)
if (reference=="GDSC1") and (analysis_tissue=="PDTC"):
    plt.ylim((-0.25,0.4))
    plt.yticks([-0.2,-0.1,0,0.1,0.2,0.3,0.4])
plt.ylabel("Correlation(predicted,actual)")
plt.xlabel("Number of samples")
if (reference=="GDSC1") and (analysis_tissue=="PDTC"):
    plt.xlabel("Number of PDTC models")
else:
    plt.title(f"{reference} {analysis_tissue} results")
plt.savefig(f"/results/{reference}_{analysis_tissue}.png")

if (reference=="GDSC1") and (analysis_tissue=="PDTC"):
    TCRP_drug ={}
    linear_drug = {}
    RF_drug = {}
    KNN_drug = {}
    NN_drug = {}
    for key,value in results.items():
        if False in np.isnan(results[key]['PDTC']["Nearest Neighbour"]):
            p = np.nan_to_num(results[key]["PDTC"]["Nearest Neighbour"])
        TCRP_drug[key] = np.mean(results[key]["PDTC"]["TCRP"])
        linear_drug[key] = np.mean(results[key]["PDTC"]["Linear Regression"])
        NN_drug[key] = np.mean(results[key]["PDTC"]["Neural Network"])
        RF_drug[key] = np.mean(results[key]["PDTC"]["Random Forest"])
        KNN_drug[key]= np.mean(p)
    TCRP_drug = {k: v for k, v in sorted(TCRP_drug.items(), key=lambda item: item[1])}   
    linear_drug = dict(OrderedDict((k, linear_drug[k]) for k in list(TCRP_drug.keys())))
    NN_drug = dict(OrderedDict((k, NN_drug[k]) for k in list(TCRP_drug.keys())))
    RF_drug = dict(OrderedDict((k, RF_drug[k]) for k in list(TCRP_drug.keys())))
    KNN_drug = dict(OrderedDict((k, KNN_drug[k]) for k in list(TCRP_drug.keys())))
    def prepare_points(model_dict):
        items = model_dict.items()
        myList = (items) 
        x, y = zip(*myList) 
        return x,y
    TCRP_x,TCRP_y = prepare_points(TCRP_drug)

    linear_x,linear_y = prepare_points(linear_drug)
    RF_x,RF_y = prepare_points(RF_drug)
    KNN_x,KNN_y = prepare_points(KNN_drug)
    NN_x,NN_y = prepare_points(NN_drug)
    drugs = pd.read_csv("/data/drug_performance.csv")["Drug"].tolist()
    new_RF_y = []
    new_KNN_y = []
    new_NN_y = []
    new_linear_y = []
    new_TCRP_y = []
    for i in drugs:
        if i in RF_x:
            new_RF_y.append(RF_y[RF_x.index(i)])
        else:
            new_RF_y.append(1)
    for i in drugs:
        if i in RF_x:
            new_KNN_y.append(KNN_y[KNN_x.index(i)])
        else:
            new_KNN_y.append(1)
    for i in drugs:
        if i in RF_x:
            new_linear_y.append(linear_y[linear_x.index(i)])
        else:
            new_linear_y.append(1)
    for i in drugs:
        if i in RF_x:
            new_NN_y.append(NN_y[NN_x.index(i)])
        else:
            new_NN_y.append(1)
    for i in drugs:
        if i in RF_x:
            new_TCRP_y.append(TCRP_y[TCRP_x.index(i)])
        else:
            new_TCRP_y.append(1)

    new_RF_y = tuple(new_RF_y)
    new_NN_y = tuple(new_NN_y)
    new_KNN_y = tuple(new_KNN_y)
    new_linear_y = tuple(new_linear_y)
    new_TCRP_y = tuple(new_TCRP_y)
    def reverse_tuple(t):
      #condition checking
        if len(t) == 0:
            return t
        else:
            return(t[-1],)+reverse_tuple(t[:-1])
    drugs = tuple(drugs)
    drugs = reverse_tuple(drugs)
    fig = plt.figure()


    fig, ax = plt.subplots()

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    fig.set_size_inches(8,12)

    RF = plt.scatter(reverse_tuple(new_RF_y),drugs,color="blue",s=100, alpha=0.3)
    RF.set_label("Random Forest")
    linear = plt.scatter(reverse_tuple(new_linear_y),drugs,color="green",s=100, alpha=0.3)
    linear.set_label("Linear Regression")
    KNN = plt.scatter(reverse_tuple(new_KNN_y),drugs,color="navy",s=100, alpha=0.3)
    KNN.set_label("Nearest Neighbours")
    NN = plt.scatter(reverse_tuple(new_NN_y),drugs,color="purple",s=100, alpha=0.3)
    NN.set_label("Neural Network")
    TCRP = plt.scatter(reverse_tuple(new_TCRP_y),drugs,color="red",s=100)
    TCRP.set_label("TCRP")
    plt.xlabel("Correlation (predicted,actual)")
    plt.xlim(-0.5,0.8)
    plt.xticks([-0.5,-0.3,-0.1,0.1,0.3,0.5,0.7])
    #plt.tight_layout()
    plt.legend(loc='lower right',prop={'size': 6})
    plt.savefig(f"/results/{reference}_{analysis_tissue}_dotplot.png")
    ground_truth = pd.read_csv("/data/drug_performance.csv")["corr"].tolist()
    drugs = pd.read_csv("/data/drug_performance.csv")["Drug"].tolist()
    mean_list = []
    for i in drugs:
        if i in results.keys():
            mean_list.append(np.mean((results[i][tissue]["TCRP"]).tolist()))
        else:
            mean_list.append(0)
    diff_list = []
    for i,j in zip(ground_truth,mean_list):
        diff_list.append(abs(i-j))
    fig = plt.figure()
    fig, ax = plt.subplots()

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    fig.set_size_inches(8,12)
    diff_list.reverse()
    drugs.reverse()
    TCRP = plt.scatter(diff_list,drugs,color="red",s=80)
    TCRP.set_label("TCRP")
    plt.xlabel("Correlation (predicted,actual)")
    plt.xlim(-0.2,0.6)
    #plt.xticks([-0.5,-0.3,-0.1,0.1,0.3,0.5,0.7])
    #plt.tight_layout()
    plt.legend(loc='lower right',prop={'size': 6})
    plt.savefig("/results/difference_dotplot.png")