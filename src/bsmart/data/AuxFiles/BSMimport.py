##### Put the path to your BSMArt directory here:

MAINpath=''

####

import argparse
import os
from datetime import datetime
import sys

import shutil
import json
from collections import OrderedDict
import collections

import importlib

import torch
import torch.nn as nn
import torch.optim as optim
from torch.autograd import Variable
from torch.utils.data import Dataset, DataLoader
import pandas as pd
import random

import math


if not os.path.isdir(MAINpath):
    sys.exit('Problem locating BSMArt root directory:\n' + MAINpath)

## Insert root directory to path so that we can import our modules ...
sys.path.insert(0, MAINpath)

bsmartpath=os.path.join(MAINpath,'bsmart')
sys.path.insert(0, bsmartpath)

### Import needed bsmart packages: ml and BSMutil

from bsmart import BSMutil
from bsmart import ml


#########################################################################################
############################   AL MODEL TRAINER & TESTER   ##############################
#                                                                                       #
# usage: python3 BSMimport.py <json_file_name> <csv_file_name> <<model_name>> <<test>>  #
# the <arguments> are mandatory; make sure that json & csv file match                   #
# the <<arguments>> are optional:                                                       #
# <<model_name>> specify name to load or test model under                               #
# <<test>> write anything to indicate that it should test the aforementioned model.     #
# if you want to test, the model_name is mandatory!                                     #
#                                                                                       #
#########################################################################################



if len(sys.argv) < 2:
    raise sys.exit('No json file provided!')

json_file=sys.argv[1]


if not os.path.isfile(json_file):
    raise sys.exit('Json file not found!')

try:
    with open(json_file) as json_data:
        inputs = json.load(json_data, object_pairs_hook=collections.OrderedDict)
except:
    #log.error('Failed to load json file '+file)
    print('Failed to load json file '+json_file)
    raise SystemExit


### now to set up scalers and so on
variables=[x for x in inputs['Variables']]
print(variables)

## Load in any neural networks required
cwd = os.getcwd()
if 'Networks' in inputs:
    localMLpath=os.path.join(cwd,'ML')
    if os.path.isdir(localMLpath):
        sys.path.append(localMLpath)

# width of network
if 'HiddenSize' in inputs['Networks']:
    hidden_size = int(inputs['Networks']['HiddenSize'])
else: hidden_size = 100 

# depth of network
if 'HiddenLayers' in inputs['Networks']:
    hidden_layers = int(inputs['Networks']['HiddenLayers'])
else: hidden_layers = 3

# how much it changes its opinion faced with new data
if 'LearningRate' in inputs['Networks']:
    learning_rate = float(inputs['Networks']['LearningRate'])
else: learning_rate = 1e-3

## decrease learning rate by this over time (handy if network too small for number of points)
#if 'Epsilon' in inputs['Networks']:
#    epsilon = float(inputs['Networks']['Epsilon'])**(K/points_target)
#else: epsilon = 1.

# SGD momentum
if 'SGDmomentum' in inputs['Networks']:
    sgd_momentum = float(inputs['Networks']['SGDmomentum'])
else: sgd_momentum = 0.1

# weight decay
if 'WeightDecay' in inputs['Networks']:
    our_weight_decay = float(inputs['Networks']['WeightDecay'])
else: our_weight_decay = 0.0

# number of steps of discriminator training
if 'DSteps' in inputs['Networks']:
    dsteps = int(inputs['Networks']['DSteps'])
else: dsteps = 5000



csv_file=sys.argv[2]
raw_data=pd.read_csv(csv_file)

#good=raw_data.loc[raw_data['Result']==1.0]
#msmasses=good['ms']

#max_ms=msmasses.max()
#print('Max ms: %f'%max_ms)
#relevant_data=pd.DataFrame(raw_data,columns=['kappa','ms2','mqm2-ms2','mqp2','lams','Result'])

######### ARRRGH the dataset may contain elements that are out of range


#dummyinputs={}
ranges={}


for var in variables:
    tempdf=raw_data[var]
    maxvar=tempdf.max()
    minvar=tempdf.min()
    if minvar > maxvar:
        tmax=maxvar
        maxvar=minvar
        minvar=tmax
    ranges[var]={"RANGE": [minvar, maxvar]}

#dummyinputs={"Variables": ranges}
#print(dummyinputs)

downscalers,upscalers = BSMutil.create_scalers(inputs)
    ################################



goodpoints=raw_data.loc[raw_data['Result']==1.0]
badpoints=raw_data.loc[raw_data['Result']==0.0]


protogoodvars=pd.DataFrame(goodpoints,columns=variables).values.tolist()
protobadvars=pd.DataFrame(badpoints,columns=variables).values.tolist()

lowthreshold=1e-4
highthreshold=1-lowthreshold


def choices(dataset, thismany):
    return [random.choice(dataset) for i in range(thismany)]

def checkthreshold(x):
    for z, var in zip(x,variables):
        Range=inputs['Variables'][var]['RANGE']
        if z < min(Range):
            return False
        if z > max(Range):
            return False
    return True

goodvars=[x for x in protogoodvars if checkthreshold(x)]
badvars=[x for x in protobadvars if checkthreshold(x)]


print('%d good points, %d bad points' %(len(goodvars),len(badvars)))
#### Train on full dataset? IDK, don't care


train_goodvars=goodvars
train_badvars=choices(badvars,len(goodvars))


train_set=ml.ml.ClassifierDataset(downscalers,[],[])
train_set.add_some_points_balance(goodvars,badvars)


print('Training size %d '%len(train_set))


NetworkData=dict()

NetworkData['input_size']=len(variables)
NetworkData['hidden_size']=hidden_size
NetworkData['hidden_layers']=hidden_layers

D=ml.ml.BasicDiscriminator(NetworkData)

#sgd_momentum=0.1
#learning_rate=1e-2
d_optimizer = optim.SGD(D.parameters(), lr=learning_rate, momentum=sgd_momentum,weight_decay=our_weight_decay)

dsteps=400
batch_size=64
criterion = nn.BCELoss()
scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(d_optimizer, 'min',factor=0.5)

train_dataloader = DataLoader(train_set, batch_size, shuffle=True)


# get model info
#def count_parameters( model):
#    return sum(p.numel() for p in model.parameters() if p.requires_grad)

#num_params = count_parameters(D)
num_params = sum(p.numel() for p in D.parameters() if p.requires_grad)
print('Number of parameters: %i' % num_params)

## print # parameters
print('NN geometry: \n'+str(D))


#full_dataloader = DataLoader(full_set,batch_size,shuffle=True)
        ### Basic training
percent_score=100.0

torch.set_num_threads(4)


if len(sys.argv) <= 4: # 4th argument is test model

    D.eval()
    total_error=0
    total_result=0.0
    total_square=0.0
    for x_batch, y_batch in train_dataloader:
    
        d_test_result = D(Variable(x_batch, requires_grad=False))
                
                        
        for my,real in zip(d_test_result, y_batch):
            #with open('training_results.csv', 'a') as f:
            #    print(my[0].item(), file=f)
            total_result+=my
            total_square+=my**2
                    
            if real > 0.5 and my < 0.5:
                total_error=total_error+1
            elif real < 0.5 and my > 0.5:
                total_error=total_error+1    
    percent_score=float(total_error)/float(len(train_set))*100.0
    mean=float(total_result)/float(len(train_set))
    var=float(total_square)/float(len(train_set))-mean**2
    if var > 0.0:
        sigma=math.sqrt(var)
    else:
        sigma=0.0
                        
    print("Initial network before training has percent error %.2f on training set, mean %.3f, sigma %.3f" %(percent_score,mean,sigma))
    
    
    
    for d_index in range(dsteps):
                
                #total_D_error=0
        D.train()
        D.zero_grad()
        total_loss=0
        for x_batch, y_batch in train_dataloader:
                    #D.zero_grad()
            d_real_data = Variable(x_batch, requires_grad=False)
            d_real_labels=Variable(y_batch, requires_grad=False)
                    
            d_real_decision = D(d_real_data)
            #d_decision = D(x_batch)
            d_real_error = criterion(d_real_decision, d_real_labels)
                    #d_real_error = criterion(d_decision, y_batch)
            d_real_error.backward() # compute/store gradients, but don't change params
            total_loss=total_loss+d_real_error
            d_optimizer.step()
        scheduler.step(total_loss)
        if d_index % 50 ==0:  #### compute error on the training set
            D.eval()
            total_error=0
            total_result=0.0
            total_square=0.0
            for x_batch, y_batch in train_dataloader:
                    #d_test_result = D(x_batch)
                d_test_result = D(Variable(x_batch, requires_grad=False))
                for my,real in zip(d_test_result, y_batch):
                    #with open('training_results.csv', 'a') as f:
                    #    print(my[0].item(), file=f)
                    total_result+=my
                    total_square+=my**2
                    #if abs(percent_score-50.0) < 1:
                    #    print("%.3f, %.3f" %(my,real))
                    
                    if real > 0.5 and my < 0.5:
                        total_error=total_error+1
                    elif real < 0.5 and my > 0.5:
                        total_error=total_error+1    
            percent_score=float(total_error)/float(len(train_set))*100.0
            mean=float(total_result)/float(len(train_set))
            var=float(total_square)/float(len(train_set))-mean**2
            if var > 0.0:
                sigma=math.sqrt(var)
            else:
                sigma=0.0
                        
            print("Epoch %d with percent error %.2f on training set, mean %.3f, sigma %.3f, loss %.3f" %(d_index,percent_score,mean,sigma,total_loss))
            if percent_score < 5.0: ## Sufficiently trained on the training set
                print("Training cycle complete")
                break

    if len(sys.argv) > 3:
        dname = str(sys.argv[3])
    else: dname = 'my_state_dict'
    
    torch.save(D.state_dict(), dname + '.pt')
    torch.save(d_optimizer.state_dict(), dname + '_opt.pt')
    print('Saved model as %s' %dname)


if len(sys.argv) > 4:  # argument no. 4 tells you to test
    modelname = sys.argv[3]
    D.load_state_dict(torch.load(modelname + '.pt'))
    d_optimizer.load_state_dict(torch.load(modelname + '_opt.pt'))


D.eval()
total_error=0
total_result=0.0
total_square=0.0
for x_batch, y_batch in train_dataloader:

    d_test_result = D(Variable(x_batch, requires_grad=False))


    for my,real in zip(d_test_result, y_batch):
        #with open('training_results.csv', 'a') as f:
        #    print(my[0].item(), file=f)
        total_result+=my
        total_square+=my**2

        if real > 0.5 and my < 0.5:
            total_error=total_error+1
        elif real < 0.5 and my > 0.5:
            total_error=total_error+1
percent_score=float(total_error)/float(len(train_set))*100.0
mean=float(total_result)/float(len(train_set))
var=float(total_square)/float(len(train_set))-mean**2
if var > 0.0:
    sigma=math.sqrt(var)
else:
    sigma=0.0

print("Final network has percent error %.2f on training set, mean %.3f, sigma %.3f" %(percent_score,mean,sigma))

