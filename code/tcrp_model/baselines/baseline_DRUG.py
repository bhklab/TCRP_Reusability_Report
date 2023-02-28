import time
import argparse
import numpy as np
import random
import os
import glob
import sys
import pickle
import copy
from data_loading import *
from utils import *
from score import *
from layers import * 
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import RandomForestRegressor
from sklearn.neighbors import KNeighborsRegressor
from copy import deepcopy
import torch
import torch.nn as nn
from torch.autograd import Variable
from torch.nn import Linear, ReLU, CrossEntropyLoss, Sequential, Conv2d, MaxPool2d, Module, Softmax, BatchNorm2d, Dropout
from torch.optim import Adam, SGD

# Training settings	
parser = argparse.ArgumentParser()

parser.add_argument('--dataset', type=str, default='', help='dataset to perform crossvalidation on')
parser.add_argument('--tissue', type=str, default='UPPER_AERODIGESTIVE_TRACT', help='Validation tissue, using the rest tissues for training')
parser.add_argument('--drug', type=str, default='AC220', help='Treated drug')
parser.add_argument('--seed', type=int, default=19, help='Random seed.')
parser.add_argument('--K', type=int, default=5, help='Perform K shot learning')
parser.add_argument('--num_trials', type=int, default=10, help='Number of trials for unseen tissue')
#parser.add_argument('--tissue_list', type=str, default=work_dic + 'cell_line_data/tissue_cell_line_list.pkl', help='Cell line list for different tissues')
parser.add_argument('--log_folder', type=str, default='/results/Log/', help='Log folder')
parser.add_argument('--tissue_num', type=int, default=13, help='Tissue number evolved in the inner update')
parser.add_argument('--run_name', type=str, default='run', help='Run name')
parser.add_argument('--fewshot_data_path', type=str, default=None, help='Path to fewshot data')
parser.add_argument('--layers',type=int,default=1)
parser.add_argument('--lr',type=float,default=0.001)
parser.add_argument('--hidden',type=int,default=20)
parser.add_argument('--RF',type=str,default="True")
class mlp(nn.Module):

	def __init__(self, feature_dim, layer, hidden):
		
		super(mlp, self).__init__()

		self.add_module( 'linear1', nn.Linear(feature_dim, hidden) )
		#self.add_module( 'bn1', nn.BatchNorm1d(hidden) )

		if layer == 2:
			self.add_module( 'linear2', nn.Linear(hidden, hidden) )
			#self.add_module( 'bn2', nn.BatchNorm1d(hidden) )
	
		self.add_module( 'linear3', nn.Linear(hidden, 1) )

		self.layer = layer
		self.loss_fn = nn.MSELoss()

		# self._init_weights() 

	def forward(self, x, weights = None):

		if weights == None:
			
			hidden = F.relu( self._modules['linear1'](x) )
			#out = F.tanh( self._modules['bn1'](out) )
			#hidden = F.relu( self._modules['bn1']( self._modules['linear1'](x) ) )
			
			if self.layer == 2:
				hidden = F.relu( self._modules['linear2']( hidden ) )
				#out = F.tanh( self._modules['bn2'](out) )
				#hidden = F.relu( self._modules['bn2']( self._modules['linear2']( hidden ) ) )

			#out = self._modules['linear3']( hidden )	
			#out = F.tanh( out )
			out = self._modules['linear3']( hidden )

		else:
			
			#hidden = F.tanh( linear(x, weights['linear1.weight'], weights['linear1.bias']) )
			#out = batchnorm(out, weight = weights['bn1.weight'], bias = weights['bn1.bias'], momentum=1)
			#hidden = linear(x, weights['linear1.weight'], weights['linear1.bias'])
			#hidden = relu( batchnorm( hidden, weight = weights['bn1.weight'], bias = weights['bn1.bias']) )
			hidden = relu( linear(x, weights['linear1.weight'], weights['linear1.bias']) )

			if self.layer == 2:
				#hidden = F.tanh( linear(hidden, weights['linear2.weight'], weights['linear2.bias']) )
				#out = batchnorm(out, weight = weights['bn2.weight'], bias = weights['bn2.bias'], momentum=1)
				#hidden = linear( hidden, weights['linear2.weight'], weights['linear2.bias'] )
				#hidden = relu( batchnorm( hidden, weight = weights['bn2.weight'], bias = weights['bn2.bias']) )
				hidden = relu( linear( hidden, weights['linear2.weight'], weights['linear2.bias'] ) )

			#out = F.tanh( linear(hidden, weights['linear3.weight'], weights['linear3.bias']) )
			out = linear(hidden, weights['linear3.weight'], weights['linear3.bias'])

		return out, hidden

	def copy_weights(self, net):
		# Set this module's weights to be the same as those of 'net'
		for m_from, m_to in zip(net.modules(), self.modules()):
			if isinstance(m_to, nn.Linear) or isinstance(m_to, nn.BatchNorm1d):

				m_to.weight.data = m_from.weight.data.clone()

				if m_to.bias is not None:
					m_to.bias.data = m_from.bias.data.clone()

	def net_forward(self, x, weights=None):
		return self.forward(x, weights)

	# def _init_weights(self):
	# 	# Set weights to Gaussian, biases to zero
	# 	torch.manual_seed(1337)
	# 	torch.cuda.manual_seed(1337)
	# 	torch.cuda.manual_seed_all(1337)
		
	# 	#print 'init weights'
		
	# 	for m in self.modules():
	# 		if isinstance(m, nn.BatchNorm1d):
	# 			m.weight.data.fill_(1)
	# 			m.bias.data.zero_()
	# 		elif isinstance(m, nn.Linear):
	# 			n = m.weight.size(1)
	# 			m.weight.data.normal_(0, 0.01)
	# 			#m.bias.data.zero_() + 1
	# 			m.bias.data = torch.ones(m.bias.data.size())

args = parser.parse_args()
dataset = args.dataset
#work_dic = '/share/data/jinbodata/siqi/Cancer_Drug_Xenograft/'
#data_dic = '/share/data/jinbodata/siqi/Cancer_Drug_Xenograft/tissue_test_data/'
#work_dic = '/cellar/users/samsonfong/Projects/tcrp-v2/from-ma/cell_line_lists/'
work_dic = f"/data/drug_features/{dataset}/"
#data_dic = '/cellar/users/samsonfong/Projects/tcrp-v2/from-ma/drug_feature/'
data_dic = f"/data/drug_features/{dataset}/drug_feature/"
filepath = os.path.realpath(__file__)
dir_name = os.path.dirname(filepath)

job_directory = "/data" + '/output/{}/'.format(args.run_name)
if args.fewshot_data_path is None:
	fewshot_data_path = job_directory
else: 
	fewshot_data_path = args.fewshot_data_path

K = args.K
num_trials = args.num_trials

random.seed(args.seed)
np.random.seed(args.seed)
torch.manual_seed(args.seed)

drug_tissue_list = work_dic + args.drug + '_tissue_cell_line_list.pkl'
with open(drug_tissue_list, 'rb') as f:
	drug_tissue_map = pickle.load(f)
drug_tissue_map = {k: v for k, v in drug_tissue_map.items() if len(v) > 0}

# Load data
#train_feature, train_label, tissue_index_list, drug_test_feature, drug_test_label, _ = load_data( drug_tissue_map, args.tissue, args.drug, path=data_dic )
train_feature, train_label, tissue_index_list, drug_test_feature, drug_test_label, _ = load_data_cell_line( drug_tissue_map, args.drug, args.tissue, K, path=data_dic )
feature_dim = train_feature.shape[1]

# Here the training process starts

unseen_train_loader_list = []
unseen_test_loader_list = []

# testing_path_suffix = data_dic + args.drug + '/' + args.tissue + '/'
test_data_path = fewshot_data_path + "/" + args.drug + '/' + args.tissue + '/' 

unseen_train_loader_list, unseen_test_loader_list = [], []

for trial in range(num_trials):
	
	# Sample a few shot learning task. Here we use k training, and use the rest for testing. 
	#unseen_train_loader, unseen_test_loader = get_unseen_data_loader(drug_test_feature, drug_test_label, K, args.inner_batch_size)
	unseen_train, unseen_test = [], []

	for k in range(1,K+1):
		# # Sample a few shot learning task. Here we use k training, and use the rest for testing. 

		train_data = np.load(test_data_path + '{}_{}_{}-shot_{}-trial_train.npz'.format(args.drug, args.tissue, k, trial))
		train_X = train_data['train_X']
		train_y = train_data['train_y']
		unseen_train_loader = [(train_X, train_y)]
		
		test_data = np.load(test_data_path + '{}_{}_{}-trial_test.npz'.format(args.drug, args.tissue, trial))
		test_X = test_data['test_X']
		test_y = test_data['test_y']
		unseen_test_loader = [(test_X, test_y)]

		unseen_train.append( unseen_train_loader )
		unseen_test.append( unseen_test_loader )
	 
	unseen_train_loader_list.append(unseen_train)
	unseen_test_loader_list.append(unseen_test)


def train_linear_baseline(Regressor, train_X, train_y, zero_train, zero_test, 
	unseen_train, unseen_test, **kwargs): 

	model = Regressor(**kwargs)
	model.fit(train_X, train_y.ravel())
	zero_p = make_predictions(model, zero_train, zero_test)

	performances = []
	for nt in range(num_trials): 
		inner_p = []
		for k in range(K):
			tmp1 = unseen_train[nt][k]
			tmp2 = unseen_test[nt][k]

			fs_train_X, fs_train_y = tmp1[0]
			fs_test_X, fs_test_y = tmp2[0]

			X = np.vstack([fs_train_X, train_X])
			y = np.vstack([fs_train_y, train_y])

			model = Regressor(**kwargs)
			model.fit(X, y.ravel())
			out = make_predictions(model, fs_test_X, fs_test_y)
			
			inner_p.append(out)
		performances.append(inner_p)

	performances = np.vstack(performances)

	return zero_p, performances
def torch_train(model,criterion,optimizer,epoch,test_X,test_y):
	tr_loss = 0
	original_test_y = test_y
	# getting the training set
	# getting the validation set
	test_X = test_X.astype(np.float32)
	test_y = test_y.astype(np.float32)
	test_X  = torch.from_numpy(test_X)
	test_y  = torch.from_numpy(test_y)
	x_val, y_val = Variable(test_X), Variable(test_y)
	# converting the data into GPU format
	if torch.cuda.is_available():
		x_val = x_val.cuda()
		y_val = y_val.cuda()

	# clearing the Gradients of the model parameters
	optimizer.zero_grad()

	# prediction for training and validation set
	predictions = model(x_val)

	# computing the training and validation loss
	# loss_val = criterion(predictions, y_val)
	# val_losses.append(loss_val)

	# # computing the updated weights of all the model parameters
	# loss_train.backward()
	# optimizer.step()
	# tr_loss = loss_train.item()
	# if epoch%2 == 0:
	# 	# printing the validation loss
	# 	print('Epoch : ',epoch+1, '\t', 'loss :', loss_val)
	saved_predictions = predictions[1].cpu().detach().numpy()
	saved_predictions = (np.mean(saved_predictions,axis=1)).reshape((test_y.shape[0],1))
	out = np.corrcoef(np.vstack([saved_predictions.ravel(), original_test_y.ravel()]))
	return out[0,1]
def train_cnn(train_X, train_y, zero_train, zero_test, 
	unseen_train, unseen_test, epoch=20,**kwargs): 
	# train_X  = torch.from_numpy(train_x)
	# train_y = train_y.astype(int);
	# train_y = torch.from_numpy(train_y)

	#change parameters here
	model = mlp(zero_train.shape[1],args.layers,args.hidden)
	# defining the optimizer
	optimizer = Adam(model.parameters(), lr=args.lr)
	# defining the loss function
	criterion = CrossEntropyLoss()
	# checking if GPU is available
	if torch.cuda.is_available():
		model = model.cuda()
		criterion = criterion.cuda()
	model.train()
	zero_list = []
	for i in range(epoch):
		zero_list.append(torch_train(model,criterion,optimizer,epoch,zero_train,zero_test))
	zero_p = np.mean(zero_list)
	performances = []
	for nt in range(num_trials): 
		inner_p = []
		for k in range(K):
			tmp1 = unseen_train[nt][k]
			tmp2 = unseen_test[nt][k]

			fs_train_X, fs_train_y = tmp1[0]
			fs_test_X, fs_test_y = tmp2[0]

			X = np.vstack([fs_train_X, train_X])
			y = np.vstack([fs_train_y, train_y])

			# model = Regressor(**kwargs)
			# model.fit(X, y.ravel())
			# out = make_predictions(model, fs_test_X, fs_test_y)
			out_list = []
			for i in range(epoch):
				out_list.append(torch_train(model,criterion,optimizer,epoch,fs_test_X,fs_test_y))
			out = np.mean(out_list)
			inner_p.append(out)
		performances.append(inner_p)
	performances = np.vstack(performances)
	return zero_p, performances

def make_predictions(model, X, y): 
	predictions = model.predict(X)
	out = np.corrcoef(np.vstack([predictions.ravel(), y.ravel()]))

	return out[0,1]


base_line_outpath = f"/results/{dataset}/baseline_performances/" + args.drug + '/' + args.tissue + '/'
os.system("mkdir -p {}".format(base_line_outpath))

if args.RF == "True":
    models = [
        ("Linear Regression", LinearRegression, {}), 
        ("Nearest Neighbour", KNeighborsRegressor, {'n_neighbors':2}), 
        ("Random Forest", RandomForestRegressor, {'n_estimators': 100, 'n_jobs': -1})
    ]
else:
    models = [
        ("Linear Regression", LinearRegression, {}), 
        ("Nearest Neighbour", KNeighborsRegressor, {'n_neighbors':2})
    ]
zero_cnn,cnn_out = train_cnn(train_feature,train_label,drug_test_feature,drug_test_label,unseen_train_loader_list,unseen_test_loader_list)
results = {}
train_label = train_label[:train_feature.shape[0],:]
results["{}-zero".format("Neural Network")] = np.array([zero_cnn])
results["{}-fewshot".format("Neural Network")] = cnn_out
cnn_mean = np.mean(cnn_out)
for name, model, kwargs in models:
	print("Training...", name) 
	train_label[np.where(np.isinf(train_label))] = 0
	zero_p, p = train_linear_baseline(model, train_feature, train_label, drug_test_feature, drug_test_label, 
		unseen_train_loader_list, unseen_test_loader_list, **kwargs)
	print("Done")
	results["{}-zero".format(name)] = np.array([zero_p])
	results["{}-fewshot".format(name)] = np.nan_to_num(p)
mean_dict = {}
for key,value in results.items():
	mean_dict[key] = np.mean(value)
median_dict = {}
for key,value in results.items():
	median_dict[key] = np.median(value)
np.savez(
	base_line_outpath + "baseline_performance", 
	**results
)

