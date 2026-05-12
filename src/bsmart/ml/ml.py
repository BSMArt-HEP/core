"""
Machine Learning utility routines for BSMArt; aim is these can either be imported by a scan or by external programs to
load/create networks used by the scan.

"""

import numpy as np
import math

import random as random

import torch
import torch.nn as nn
import torch.optim as optim
from torch.autograd import Variable

from torch.utils.data import Dataset, DataLoader

import sys


def choices(dataset, thismany):
    return [random.choice(dataset) for i in range(thismany)]

class ClassifierDataset(Dataset):
    def __init__(self,scalerfuncs,good_points,bad_points,ratio=1.0):
        self.points=[]
        self.labels=[]
        self.scalerfuncs=scalerfuncs
        self.ngood=0
        self.nbad=0
        self.add_some_points(good_points,bad_points)

    def __len__(self):
        return len(self.points)

    def __getitem__(self, idx):
        return self.points[idx], self.labels[idx]

    def goodratio(self):
        if self.nbad==0:
#            if self.ngood==0:
#                return 1.0
#            else:
            return -1
        else:   
            return float(self.ngood)/float(self.nbad)

    def add_some_points(self, good_points,bad_points):
        
        

        for x in good_points:
            self.points.append(torch.Tensor([sf(y) for sf,y in zip(self.scalerfuncs,x)]))
            #print('Good point '+str(self.points[-1]))
            self.labels.append(torch.Tensor([1]))
        self.ngood+=len(good_points)
        #for x in maxset[:maxlen]:
        for x in bad_points:
            self.points.append(torch.Tensor([sf(y) for sf,y in zip(self.scalerfuncs,x)]))
            self.labels.append(torch.Tensor([0]))
        self.nbad+=len(bad_points)


    def add_some_points_balance(self, good_points,bad_points):
        goodlen=len(good_points)
        badlen=len(bad_points)
        if goodlen==0 or badlen ==0:
            return

        if goodlen <= badlen:
            maxset=bad_points
            
            #maxset=bad_points
            
            minval=0.0
            maxval=1.0
            maxlen=badlen
            minlen=goodlen
            if minlen == 0:
                return
            #if maxlen/minlen > ratio:
            #    maxlen=int(minlen*ratio)
            if sys.version_info[1] < 6: minset=choices(good_points,maxlen)
            else: minset=random.choices(good_points,k=maxlen)
            
        else:
            maxset=good_points
            #maxset=good_points
            minval=1.0
            maxval=0.0
            maxlen=goodlen
            minlen=badlen
            if minlen == 0:
                return
            #if maxlen/minlen > ratio:
            #    maxlen=int(minlen*ratio)
            #maxset=random.sample(good_points,maxlen)
            if sys.version_info[1] < 6: minset=choices(bad_points,maxlen)
            else: minset=random.choices(bad_points,k=maxlen)
        
        
            
        for x in minset:
            self.points.append(torch.Tensor([sf(y) for sf,y in zip(self.scalerfuncs,x)]))
            #print('Good point '+str(self.points[-1]))
            self.labels.append(torch.Tensor([minval]))
        #for x in maxset[:maxlen]:
        for x in maxset:
            self.points.append(torch.Tensor([sf(y) for sf,y in zip(self.scalerfuncs,x)]))
            self.labels.append(torch.Tensor([maxval]))


class BSMNetwork(nn.Module):
    """ 
    Template for user-defined networks. The only restriction is that we always pass inputs 
    in the form of a dictionary to have a common interface.  
    """
    def __init__(self,inputs):
        super().__init__()
        self.inputs=inputs

    def forward(self,x):
        raise NotImplementedError("You need to implement the function 'forward' in the BSMNetwork class")


class BasicDiscriminator(BSMNetwork):
    def __init__(self,inputs):
        super().__init__(inputs)
        self.hiddenmaps = nn.ModuleList()
        self.map1 = nn.Linear(int(inputs['input_size']),int(inputs['hidden_size']))
        
        for _ in range(int(inputs['hidden_layers'])):
            self.hiddenmaps.append(nn.Linear(int(inputs['hidden_size']),int(inputs['hidden_size'])))
        self.lastmap=nn.Linear(int(inputs['hidden_size']),1)
        self.relu=nn.ReLU()
        self.f = torch.sigmoid

    def forward(self,x):
        x = self.f(self.map1(x))
        for hmap in self.hiddenmaps:
            x=self.relu(hmap(x))
        return self.f(self.lastmap(x))


class BasicRegressor(BSMNetwork):
    """
        Neural network mapping to one output layer
    """
    def __init__(self,inputs):
        super().__init__(inputs)
        self.hiddenmaps = nn.ModuleList()
        self.map1 = nn.Linear(int(inputs['input_size']),int(inputs['hidden_size']))
        
        for _ in range(int(inputs['hidden_layers'])):
            self.hiddenmaps.append(nn.Linear(int(inputs['hidden_size']),int(inputs['hidden_size'])))
        self.lastmap=nn.Linear(int(inputs['hidden_size']),1)
        self.relu=nn.ReLU()
        self.f = torch.sigmoid

    def forward(self,x):
        x = self.f(self.map1(x))
        for hmap in self.hiddenmaps:
            x=self.relu(hmap(x))
        return self.lastmap(x)


"""
class Discriminator(nn.Module):
    def __init__(self, input_size, hidden_size, output_size, f):
        super().__init__()
        self.map1 = nn.Linear(input_size, hidden_size)
        self.map2 = nn.Linear(hidden_size, hidden_size)
        self.map3 = nn.Linear(hidden_size, hidden_size)
        #self.map4 = nn.Linear(hidden_size, hidden_size)
        self.map5= nn.Linear(hidden_size, output_size)
       # self.leak=nn.LeakyReLU(0.02)
        self.relu=nn.ReLU()
        self.f = f
    
    def forward(self, x):
        x = self.f(self.map1(x))
    
        #x = self.f(self.map2(x))
        #x = nn.functional.dropout(x,p=0.1)
        #x = self.leak(self.map2(x))
        x=self.relu(self.map2(x))
        x=self.relu(self.map3(x))
        #x=self.relu(self.map4(x))
        return self.f(self.map5(x))
"""

