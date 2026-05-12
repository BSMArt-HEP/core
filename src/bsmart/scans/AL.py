"""
Active Learning scan using neural networks


Active learning scan with a customizable neural network (see settings in JSON file), a smart point selection algorithm which encourages searches around existing good points while maintaining point diversity, and some safeguards against premature training failures.

You can find more information about this type of scan in:
Goodsell, Mark D. and Joury, Ari: "Active learning BSM parameter spaces", `arXiv:2204.13950 <https://arxiv.org/abs/2204.13950>`_, April 2022.
Please cite this paper if you are using this scan.

"""

__meta__ = {
    "name": "AL",
    "requires": ["torch", "pandas", "scipy", "numpy"],
    "settings": {
        "Cores": "Number of cores",
        "Points": "Number of points",
        "InitCSV": "Initial CSV file",
        "TestCSV": "Test CSV file",
        "Benchmark": "Benchmark CSV file",
        "Networks": "Dict of network settings"
    }
}



import sys
from bsmart.core import Scan as Scan
import itertools
from bsmart import debug
from collections import OrderedDict
import numpy as np
import math

import random as random
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
from scipy.stats import truncnorm

import torch
import torch.nn as nn
import torch.optim as optim
from torch.autograd import Variable

from torch.utils.data import Dataset, DataLoader

import time
import os

# matplotlib_is_available = True
# try:
#   from matplotlib import pyplot as plt
# except ImportError:
#   self.log.debug("Will skip plotting; matplotlib is not available.")
#   matplotlib_is_available = False



def create_edge_penalty(a,b):
    #a0=1.0-(b+0.25*a)/(b+a*0.25)
    return lambda x: (1.0-0.25*a/(b+a*x*(1-x)))/12

def choices(dataset, thismany):
    return [random.choice(dataset) for i in range(thismany)]


"""
def scaler_func(vmin,vmax):
    diff = vmax-vmin
    return lambda x : (x-vmin)/diff
"""
def scaler_func(vmin,vmax):
    diff = vmax-vmin
    
    def tfunc(x):
        res = (x-vmin)/diff
        if res < 0:
            return 0.0
        elif res > 1.0:
            return 1.0
        else:
            return res
        
    return tfunc


def ClassifyFunc(my_type,mean,var):
    """ 
    Generates a lambda function to decide whether a point is excluded or not based on one observable.
    """
    tmean=mean
    tvar=var
    if my_type == 'RANGE':
        tmax=max(tmean+tvar,tmean-tvar)
        tmin=min(tmean+tvar,tmean-tvar)
        def tfunc(x):
            if x > tmax:
                return 0.0
            elif x < tmin:
                return 0.0
            else:
                return 1.0
        return tfunc
    elif my_type == 'UPPER':
        tmax=mean
        def tfunc(x):
            if x > tmax:
                return 0.0
            else:
                return 1.0
        return tfunc
    elif my_type == 'LOWER':
        tmin=mean
        def tfunc(x):
            if x < tmin:
                return 0.0
            else:
                return 1.0
        return tfunc
    elif my_type == 'USER':
        return lambda x : x
    #elif my_type == 'OFF':
    else:    # if not specified don't include it
        return lambda  x : 1.0 # Point is allowed


class MyDataset(Dataset):
    def __init__(self,scalerfuncs,good_points,bad_points,ratio=1.0):
        self.points=[]
        self.labels=[]
        self.add_some_points_balance(scalerfuncs,good_points,bad_points)

    def __len__(self):
        return len(self.points)

    def __getitem__(self, idx):
        return self.points[idx], self.labels[idx]

    def add_some_points_balance(self, scalerfuncs,good_points,bad_points):
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
            else:  minset=random.choices(good_points,k=maxlen)
            
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
            self.points.append(torch.Tensor([sf(y) for sf,y in zip(scalerfuncs,x)]))
            self.labels.append(torch.Tensor([minval]))
        #for x in maxset[:maxlen]:
        for x in maxset:
            self.points.append(torch.Tensor([sf(y) for sf,y in zip(scalerfuncs,x)]))
            self.labels.append(torch.Tensor([maxval]))



class Discriminator(nn.Module):
    def __init__(self, input_size, hidden_size, output_size, hidden_layers, f):
        super().__init__()
        self.hiddenmaps=nn.ModuleList()
        self.map1 = nn.Linear(input_size,hidden_size)

        for _ in range(hidden_layers):
            self.hiddenmaps.append(nn.Linear(hidden_size,hidden_size))
        self.lastmap=nn.Linear(hidden_size,1)
        self.relu=nn.ReLU()
        self.f = torch.sigmoid

    def forward(self,x):
        x = self.f(self.map1(x))
        for hmap in self.hiddenmaps:
            x=self.relu(hmap(x))
        return self.f(self.lastmap(x))



class NewScan(Scan):
    """Scanner class for Random Scans"""

    def __init__(self, inputs, log):
        Scan.__init__(self, inputs, log)

        self.total_loss=0
        self.header=[]
        self.header.append('ID')
        self.variances = []
        self.variables = []
        self.batch_size = 64
        
        self.edge_penalty = create_edge_penalty(10,0.01)
        if "Cores" in self.inputs['Setup']:
            self.ncores=int(self.inputs['Setup']['Cores'])
            ####### This is broken on pytorch < 1.9
            ######## use setenv OMP_NUM_THREADS 1; python3 class_threads.py or whatever!
            #torch.set_num_threads(math.ceil(self.ncores/2))
            ## One thread per worker, plus one
            self.log.info('Setting number of cores to %d' % self.ncores)
            #torch.set_num_threads(self.ncores-1)
            torch.set_num_threads(8)
            #torch.set_num_interop_threads(2)
        if "Points" in self.inputs['Setup']:
            self.points_target=int(self.inputs['Setup']['Points'])
            
        
        log.debug('Creating variables')
        for v in self.inputs['Variables']:
            self.variances.append(self.inputs['Variables'][v]['VARIANCE'])
            self.header.append(v) ## append the variable name ...
            self.variables.append(v)
        #self.scalers=self.createscalers()
        log.debug('Creating classifiers')
        self.classifiers = []
        for o in self.inputs['Observables']:
            self.header.append(o)
            obs=self.inputs['Observables'][o]
            class_type='OFF'
            mean=0.0
            var=0.0
            try:
                if 'MEAN' in obs:
                    mean=float(obs['MEAN'])
                if 'VARIANCE' in obs:
                    var=float(obs['VARIANCE'])
                if 'TYPE' in obs:
                    class_type=obs['TYPE']
                if 'MAX' in obs:
                    class_type='UPPER'
                    mean=float(obs['MAX'])
                elif 'MIN' in obs:
                    class_type='LOWER'
                    mean=float(obs['MIN'])
                elif 'RANGE' in obs:
                    class_type='RANGE'
                    trange=obs['RANGE']
                    var=abs(float(trange[1])-float(trange[0]))/2.0
                    mean=(float(trange[1])+float(trange[0]))/2.0
                
                self.log.info('Creating classifier for %s of type %s with mean %f and variance %f' %(o,class_type, mean, var))
                self.classifiers.append(ClassifyFunc(class_type,mean,var))
            except Exception as e:
                log.error('Failed to set up classifier, %s' %e)
                raise

        if 'Diversity Alpha' in self.inputs['Networks']:
            self.diversity_alpha=float(self.inputs['Networks']['Diversity Alpha'])
        else:
            self.diversity_alpha=0.75
        # number K points to add by AL
        if 'K' in self.inputs['Networks']:
            self.K = int(self.inputs['Networks']['K'])
        else: self.K = 10
        if 'Kinitial' in self.inputs['Networks']:
            self.Kinitial = int(self.inputs['Networks']['Kinitial'])
        else: self.Kinitial = self.K


        var_columns = list(self.inputs['Variables'].keys())
        obs_columns = list(self.inputs['Observables'].keys())
        self.got_csv = False # tells us whether we have a functioning input csv file or not
        self.got_benchmark = False  # tells us whether we have a random points dataset for benchmarking


        ## Load CSV
        if 'InitCSV' not in self.inputs['Setup']:
            #raise NameError('No initial csv points given')
            log.debug('No initial CSV file given.') # it's ok, we can totally run without an initial csv file
        else:
            self.input_csv_file=self.inputs['Setup']['InitCSV']

            if not os.path.isfile(self.input_csv_file):
                self.log.error('CSV file not found')
                raise SystemExit

            self.got_csv=True
            log.debug('Loading csv inputs from %s' %self.input_csv_file)
            input_dataframe = pd.read_csv(self.input_csv_file)
            useful_dataset=pd.DataFrame(input_dataframe, columns=var_columns+obs_columns).dropna()
            log.info('Dataframe has %d points' % useful_dataset.shape[0])
            self.create_scalers_from_data(useful_dataset,var_columns)


        ## prepare initial dataset
        log.debug('Creating dataframe')
        if self.got_csv is False: self.create_scalers()

        self.good_points=[]
        self.bad_points=[]
        self.last_good_points=[]
        self.number_seed_points=int(self.K/2)
        firstobs=len(var_columns)

        if self.got_csv is True:
            for row in useful_dataset.itertuples(index=False):
                robs=list(row[firstobs:])
                rvars=list(row[:firstobs])
                res=self.classifypoint(robs)
                if res!=[] and res > 0.5:
                    self.good_points.append(rvars)
                else:
                    self.bad_points.append(rvars)

        else:
            new_random_vars = self.scale_points([[np.random.uniform(0,1) for x in self.variables] for _ in range(self.K)],down=False)
            results=self.RunManager.run_batch(new_random_vars)
            for rvars,res in zip(new_random_vars,results):
                if res!=[] and res > 0.5:
                    self.good_points.append(rvars)
                else:
                    self.bad_points.append(rvars)

        
        if len(self.good_points) <= self.number_seed_points:
            self.last_good_points = self.good_points
        else:
            self.last_good_points = self.good_points[:self.number_seed_points] 


        self.good_test_points = []
        self.bad_test_points = []

        ## Load test CSV
        self.test_csv = False  # bool to see if we have a test csv
        if 'TestCSV' not in self.inputs['Setup']:
            #raise NameError('No initial csv points given')
            log.debug('No test CSV file given.') # we don't always need a test csv file
        else:
            self.test_csv_file=self.inputs['Setup']['TestCSV']
            if not os.path.isfile(self.test_csv_file):
                self.log.error('Test CSV file not found')
                raise SystemExit
            self.test_csv=True
            log.debug('Loading csv inputs from %s' %self.test_csv_file)
            test_dataframe = pd.read_csv(self.test_csv_file)
            test_dataset=pd.DataFrame(test_dataframe, columns=var_columns+obs_columns).dropna()
            log.info('Dataframe has %d points' % test_dataset.shape[0])
            #self.create_scalers_from_data(test_dataset,var_columns)

            for row in test_dataset.itertuples(index=False):
                robs=list(row[firstobs:])
                rvars=list(row[:firstobs])
                res=self.classifypoint(robs)
                if res!=[] and res > 0.5:
                    self.good_test_points.append(rvars)
                else:
                    self.bad_test_points.append(rvars)

        if self.test_csv == False:

            #split off enough points to test with
            sample_size = min(len(self.good_points),len(self.bad_points))
            test_size = int(0.2*sample_size)
            if test_size > 0:
                to_test = set(random.sample(range(sample_size), test_size))
            else:
                to_test = []
            self.good_points = [x for i,x in enumerate(self.good_points) if not i in to_test]
            self.bad_points = [x for i,x in enumerate(self.bad_points) if not i in to_test]
            self.good_test_points = [x for i,x in enumerate(self.good_points) if i in to_test]
            self.bad_test_points = [x for i,x in enumerate(self.bad_points) if i in to_test]

        self.test_set = MyDataset(self.scalers,self.good_test_points,self.bad_test_points,1.0)
        log.info('Test set has %d points' % len(self.test_set))


        ## create initial train set
        self.train_set=MyDataset(self.scalers,self.good_points,self.bad_points,1.0)

        #train_size=len(self.initial_dataset)
        log.info('Initial training set has %d points' % len(self.train_set))
        

        ## Load Benchmark
        if 'Benchmark' not in self.inputs['Setup']:
            log.debug('No benchmark CSV file given.') # we don't always want a benchmark file
        else:
            self.benchmark_file=self.inputs['Setup']['Benchmark']

            if not os.path.isfile(self.benchmark_file):
                self.log.error('Benchmark CSV file not found')
                raise SystemExit

            self.got_benchmark=True
            log.debug('Loading benchmark csv inputs from %s' %self.benchmark_file)
            benchmark_dataframe = pd.read_csv(self.benchmark_file)
            benchmark_dataset=pd.DataFrame(benchmark_dataframe, columns=var_columns+obs_columns).dropna()
            log.info('Benchmark dataframe has %d points' % benchmark_dataset.shape[0])

            self.good_bpoints=[]
            self.bad_bpoints=[]
            firstobs=len(var_columns)

            for row in benchmark_dataset.itertuples(index=False):
                robs=list(row[firstobs:])
                rvars=list(row[:firstobs])
                res=self.classifypoint(robs)
                if res!=[] and res > 0.5:
                    self.good_bpoints.append(rvars)
                else:
                    self.bad_bpoints.append(rvars)
                    
            self.benchmarkset = MyDataset(self.scalers,self.good_bpoints,self.bad_bpoints,1.0)

            log.info('Benchmark set has %d points' % len(self.benchmarkset))


        
        ## CREATE DISCRIMINATOR
        log.debug('Creating discriminator')

        # depth of network
        if 'HiddenLayers' in self.inputs['Networks']:
            self.hidden_layers = int(self.inputs['Networks']['HiddenLayers'])
        else: self.hidden_layers = 3

        # width of network
        if 'HiddenSize' in self.inputs['Networks']:
            self.hidden_size = int(self.inputs['Networks']['HiddenSize'])
        else: self.hidden_size = 300

        # how much it changes its opinion faced with new data
        if 'LearningRate' in self.inputs['Networks']:
            self.learning_rate = float(self.inputs['Networks']['LearningRate'])
        else: self.learning_rate = 1e-3

        # decrease learning rate by this over time (handy if network too small for number of points)
        if 'Epsilon' in self.inputs['Networks']:
            self.epsilon = float(self.inputs['Networks']['Epsilon'])**(self.K/self.points_target)
        else: self.epsilon = 1.

        # SGD momentum
        if 'SGDmomentum' in self.inputs['Networks']:
            self.sgd_momentum = float(self.inputs['Networks']['SGDmomentum'])
        else: self.sgd_momentum = 0.1

        # weight decay
        if 'WeightDecay' in self.inputs['Networks']:
            self.weight_decay = float(self.inputs['Networks']['WeightDecay'])
        else: self.weight_decay = 0.0

        # number of steps of discriminator training
        if 'DSteps' in self.inputs['Networks']:
            self.d_steps = int(self.inputs['Networks']['DSteps'])
        else: self.d_steps = 5000

        # number L points to draw K from
        if 'L' in self.inputs['Networks']:
            self.L = int(self.inputs['Networks']['L'])
        else: self.L = 500
        
        # specify how many new points to propose from good points
        if 'FromGood' in self.inputs['Networks']:
            self.fromGood=float(self.inputs['Networks']['FromGood'])
        else: self.fromGood=0.

        # get name for model file
        if 'MLmodel' in self.inputs['Networks']:
            self.modelfile = self.inputs['Networks']['MLmodel']
        else:
            self.log.info("WARNING: No modelfile specified in JSON input file. Add 'MLmodel' in section 'Setup' to create one.")
            self.modelfile = None
        
        # specify if we want to train on full set after this many training rounds
        if 'FullTrain' in self.inputs['Networks']:
            self.fulltrain=int(self.inputs['Networks']['FullTrain'])
        else: self.fulltrain = 0
        
        # specify if we want to abort training if it's bad 
        if 'AutoStop' in self.inputs['Networks']:
            if self.inputs['Networks']['AutoStop'] != 'True': self.autostop = False
        else: self.autostop = True


        ## Initialise discriminator
        #self.output_size=len(obs_columns)
        self.output_size=1
        self.criterion = nn.BCELoss()
        
        self.D=Discriminator(len(self.variables),self.hidden_size,self.output_size,self.hidden_layers,f=torch.sigmoid)

        self.log.info('NN geometry: \n'+str(self.D))       
        #for name, param in self.D.named_parameters():
        #    if param.requires_grad:
        #        self.log.info(name, param.data)
        
        self.d_optimizer = optim.SGD(self.D.parameters(), lr=self.learning_rate, momentum=self.sgd_momentum, weight_decay=self.weight_decay)

        ## show # parameters
        num_params = self.count_parameters(self.D)
        self.log.info('Number of parameters: %i' % num_params)

        self.train_score = self.initial_training()

        self.n_points_run=len(self.train_set)
        
        ## After this initialisiation, run() adds points


    def initial_training(self):
        self.test_result=50
        initial_training_counter=0
        while self.test_result > 40:
            ## initial training & testing
            self.log.info('Start initial training')
            ## if too few points there is not much point in training
            if len(self.train_set) < 20:
                self.train_score=100
            else:
                self.train_score = self.do_train(self.train_set, self.d_steps)
            if len(self.test_set) ==0:
                self.test_result=0
            else:
                self.test_result=self.do_test(self.test_set)

            ## catch case where input points were bad & discriminator did bogus: reinitialize dataset & discriminator
            if self.train_score > 40:
                self.log.info('Initial training too bad. Reinitializing model & initial dataset. Try number %i' %initial_training_counter)
                self.test_set = MyDataset(self.scalers,self.good_test_points,self.bad_test_points,1.0)
                self.train_set = MyDataset(self.scalers,self.good_points,self.bad_points,1.0)
                self.D=Discriminator(len(self.variables),self.hidden_size,self.output_size,self.hidden_layers,f=torch.sigmoid)
                self.d_optimizer = optim.SGD(self.D.parameters(), lr=self.learning_rate, momentum=self.sgd_momentum, weight_decay=self.weight_decay)

            initial_training_counter += 1
            if initial_training_counter > 9:
                raise SystemExit('Model or dataset is bad. Consider changing settings , e.g. Network or ranges.')
            self.log.info('Finished initial training!')

        return self.train_score 



    # get model info
    def count_parameters(self, model):
        return sum(p.numel() for p in model.parameters() if p.requires_grad)


    def do_test(self, set_to_test):
        if len(set_to_test)<=0: 
            self.log.info('No points to test. Using test set')
            testset=self.test_set
            
        else:
            testset=set_to_test
        if len(testset) == 0:
            return 100
        test_dataloader = DataLoader(testset, self.batch_size, shuffle=True)
        self.D.eval()
        total_error=0
        total_result=0.0
        total_square=0.0
        for x_batch, y_batch in test_dataloader:
            d_test_result = self.D(x_batch)
                       
                
                #for xx,my,real in zip(x_batch,d_test_result, y_batch):
            for my,real in zip(d_test_result, y_batch):
                total_result+=my
                total_square+=my**2
                #if abs(percent_score-50.0) < 1:
                #    self.log.debug("%.3f, %.3f" %(my,real))
                #self.log.debug('[%.3f, %.3f, %.3f, %.3f, %.3f, %.3f] -> mine: %f, real: %f' %(xx[0],xx[1],xx[2],xx[3],xx[4],xx[5],my,real))
                if real > 0.5 and my < 0.5:
                    total_error=total_error+1
                elif real < 0.5 and my > 0.5:
                    total_error=total_error+1    
        percent_score=float(total_error)/float(len(testset))*100.0
        mean=float(total_result)/float(len(testset))
        var=float(total_square)/float(len(testset))-mean**2
        if var > 0.0:
            sigma=math.sqrt(var)
        else:
            sigma=0.0
        self.log.info("Testing returned percent error %.2f on test set, mean %.3f, sigma %.3f" %(percent_score,mean,sigma))
        return percent_score


    def do_train(self, train_set, dsteps):
        train_dataloader = DataLoader(train_set, self.batch_size, shuffle=True)
        ### Basic training
        percent_score=100.0
        for d_index in range(dsteps):
            
            #total_D_error=0
            self.total_loss=0
            self.D.train()
            self.D.zero_grad()
            for x_batch, y_batch in train_dataloader:
                #self.D.zero_grad()
                d_real_data = Variable(x_batch, requires_grad=False)
                d_real_labels=Variable(y_batch, requires_grad=False)
                
                d_real_decision = self.D(d_real_data)
                #d_decision = self.D(x_batch)
                d_real_error = self.criterion(d_real_decision, d_real_labels)
                #d_real_error = self.criterion(d_decision, y_batch)
                d_real_error.backward() # compute/store gradients, but don't change params
                self.total_loss = self.total_loss+d_real_error  # need for scheduler
                self.d_optimizer.step()
            if d_index % 50 ==0:  #### compute error on the training set
                self.D.eval()
                total_error=0
                total_result=0.0
                total_square=0.0
                for x_batch, y_batch in train_dataloader:
                    #d_test_result = self.D(x_batch)
                    d_test_result = self.D(Variable(x_batch, requires_grad=False))
            
                    #for xx,my,real in zip(x_batch,d_test_result, y_batch):
                    for my,real in zip(d_test_result, y_batch):
                        total_result+=my
                        total_square+=my**2
                        #if abs(percent_score-50.0) < 1:
                        #    self.log.debug("%.3f, %.3f" %(my,real))
                #self.log.debug('[%.3f, %.3f, %.3f, %.3f, %.3f, %.3f] -> mine: %f, real: %f' %(xx[0],xx[1],xx[2],xx[3],xx[4],xx[5],my,real))
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
                    
                self.log.info("Epoch %d with percent error %.2f on training set, mean %.3f, sigma %.3f" %(d_index,percent_score,mean,sigma))
                #if percent_score < 5.0: ## Sufficiently trained on the training set
                if percent_score < 0.5: ## Sufficiently trained on the training set
                    self.log.info("Training cycle complete")
                    return percent_score
        return percent_score

    def classifypoint(self,observables):
        res=1.0
        for v, f in zip(observables,self.classifiers):
            res=res*f(float(v))
        
        return res
    
    
    #def postprocess(self,Point, observables, outslha,temp_dir,log, lock=None):
    def postprocess(self,Point, observables, data_point,temp_dir,log, lock=None):
        """ Run classification """
        res=1.0
        for v, f in zip(observables,self.classifiers):
            res=res*f(float(v))
        return res
        
    
    def create_scalers(self):
        self.scalers=[]
        for varname in self.variables:
            varmax=max(self.inputs['Variables'][varname]['RANGE'][0],self.inputs['Variables'][varname]['RANGE'][1])
            varmin=min(self.inputs['Variables'][varname]['RANGE'][0],self.inputs['Variables'][varname]['RANGE'][1])
            diff=varmax-varmin
            self.log.info('Creating scaler for %s between %.4e and %.4e' %(varname,varmin,varmax))
            
            self.scalers.append(scaler_func(varmin,varmax))

    
    def create_scalers_from_data(self,thedataset,frame_variables):
        self.scalers=[]
        for varname in frame_variables:
            varmax=max(thedataset[varname].max(),self.inputs['Variables'][varname]['RANGE'][0],self.inputs['Variables'][varname]['RANGE'][1])
            varmin=min(thedataset[varname].min(),self.inputs['Variables'][varname]['RANGE'][0],self.inputs['Variables'][varname]['RANGE'][1])
            diff=varmax-varmin
            self.log.info('Creating scaler for %s between %.4e and %.4e' %(varname,varmin,varmax))
            
            self.scalers.append(scaler_func(varmin,varmax))

    def distancesq(self,a,b):
        ## dot product of unit vectors
        res=sum([(x-y)**2 for x,y in zip(a,b)])
        if res > 0.01:
            res = 0
        else:
            res = -0.0001/(0.0001+res)
        return res

    
        
    def rescore_diversity(self,full_set,last_new_point):
        ## don't need to recompute the distances each time, just the distances to each new point
        return [self.distancesq(x,last_new_point) for x in full_set]
        
            
    # propose K most useful from a total of L points
    def propose_KfromL_diverse(self, K, L, points=None):
        self.D.eval()

        if points==None: newpoints = [[np.random.uniform(0,1) for x in self.variables] for _ in range(L)]
        else: newpoints=points

        newpoints = [[float(var) for var in point] for point in newpoints]

        results=[self.D(torch.tensor(pt)).item() for pt in newpoints]

        #position_scores=[ sum([ self.edge_penalty(y) for y in x]) for x in newpoints]
        scores=np.array([x*(1-x) for x in results])

        """
        In order to make the equivalence between the scores and the diversity, we need to 
        know the average _differences_ amongst the scores, discarding outliers
        """
        
        scores_no_outliers=[ x for x in scores if x > 0.15]
        if len(scores_no_outliers) < K:
            score_standard_dev=max(scores)-min(scores)
        else:
            score_standard_dev=np.std(scores_no_outliers)


        #scores=[ x+y for x,y in zip(scores,position_scores)]
        ### penalty for being close to the edge/bonus for being near the middle
        #position_scores=[ sum([y*(1-y)*4.0 for y in x]) for x in newpoints]
        
        
        ## get the first point
        myind=np.argmax(scores)
        scores[myind]=-10000 ## give it an illegally low score, it can only be between 0 and 0.25
        returnset=[newpoints[myind]]
        distance_scores=self.rescore_diversity(newpoints,returnset[0])
        #distance_scores=[x+y for x,y in zip(distance_scores,position_scores)]
        #max_distance=max(distance_scores)
        max_distance=max(list(map(abs,distance_scores)))
        #### we need a factor of 0.25/max distance to put the diversity scores on the same footing as the uncertainty score
        if max_distance > 0:
            combined_scores=[ (1-self.diversity_alpha) * x + self.diversity_alpha* y*0.25*score_standard_dev/max_distance for x,y in zip(scores,distance_scores)]
        else:
            combined_scores=scores
        
        for _ in range(K):
            myind=np.argmax(combined_scores)
            scores[myind]=-10000 ## give it an illegally low score, it can only be between 0 and 0.25
            nextpoint=newpoints[myind]
            returnset.append(nextpoint)
            ## now compute the new scores
            new_distance_scores=self.rescore_diversity(newpoints,nextpoint)
            distance_scores=[x + y for x,y in zip(distance_scores,new_distance_scores)]
            max_distance=max(list(map(abs,distance_scores)))
            if max_distance > 0:
                combined_scores=[ (1-self.diversity_alpha) * x + self.diversity_alpha* y*0.25/max_distance for x,y in zip(scores,distance_scores)]
            else:
                combined_scores=scores

        return returnset

    
    def propose_KfromL(self, K, L, points=None):
        self.D.eval()

        if points==None: newpoints = [[np.random.uniform(0,1) for x in self.variables] for _ in range(L)]
        else: newpoints=points

        results=[self.D(torch.tensor(pt)).item() for pt in newpoints]
        scores=np.array([x*(1-x) for x in results])
        
        returnset=[]
        for _ in range(K):
            myind=np.argmax(scores)
            scores[myind]=-1 ## give it an illegally low score, it can only be between 0 and 0.25
            returnset.append(newpoints[myind])

        return returnset


    # propose points in the vicinity of good points (they're not scaled down)
    def propose_fromGood(self, numpoints):
        # python3.5 and lower doesn't have random.choices
        if sys.version_info[1] < 6:
            ingoodpoints = choices(self.last_good_points, numpoints)
        else:
            ingoodpoints = random.choices(self.last_good_points, k=numpoints)

        # list of the upper and lower range of each variable
        varranges = []
        for varname in self.variables:
            varmax=max(self.inputs['Variables'][varname]['RANGE'][0],self.inputs['Variables'][varname]['RANGE'][1])
            varmin=min(self.inputs['Variables'][varname]['RANGE'][0],self.inputs['Variables'][varname]['RANGE'][1])
            varranges.append((varmin, varmax))

        # get points
        param_points = [[truncnorm.rvs((r[0]-x[w])/v, (r[1]-x[w])/v, loc=x[w], scale=v)
                    for w,(v,r) in enumerate(zip(self.variances, varranges))] for x in ingoodpoints]

        return param_points

        """
        if sys.version_info[1] < 6: 
            param_points = [[np.random.normal(x[w],v)
                         for w,v in enumerate(self.variances)] for x in choices(self.last_good_points, numpoints)]
        else:
            param_points = [[np.random.normal(x[w],v)
                         for w,v in enumerate(self.variances)] for x in random.choices(self.last_good_points, k=numpoints)]                    
            """

    # scale up points that are currently scaled between 0 and 1
    def scale_points(self, zeroonepoints, down=False):
        if down==False:
            vars_minmax=[]
            for varname in self.variables:
                varmax=max(self.inputs['Variables'][varname]['RANGE'][0],self.inputs['Variables'][varname]['RANGE'][1])
                varmin=min(self.inputs['Variables'][varname]['RANGE'][0],self.inputs['Variables'][varname]['RANGE'][1])
                vars_minmax.append([varmin,varmax])
            scaledpoints = [[z*(x[1]-x[0])+x[0] for z,x in zip(point, vars_minmax)] for point in zeroonepoints]
        else:
            scaledpoints = [[sf(y) for sf,y in zip(self.scalers,point)] for point in zeroonepoints]
            
            #scaledpoints = [[(z-x[0])/(x[1]-x[0]) for z,x in zip(point, vars_minmax)] for point in zeroonepoints]

        return scaledpoints


    # generate parameter points from a list of scaled down points
    def generate_parameter_points(self, proposedpoints):
        run_points = self.RunManager.run_batch(proposedpoints)
        return run_points

    # save discriminator
    def save_model(self, pc=''):
        if self.modelfile is not None:
            if pc=='':
                torch.save(self.D.state_dict(), str(self.modelfile) + '.pt')
                torch.save(self.d_optimizer.state_dict(), str(self.modelfile) + '_opt.pt')
            else:
                torch.save(self.D.state_dict(), str(self.modelfile) + pc + '.pt')
                torch.save(self.d_optimizer.state_dict(), str(self.modelfile) + pc + '_opt.pt')

    
    def run(self):
        """ Cycle through generating points and training the network until either we have a lot of (good?) points or some criterion is met, e.g. the network is confident of everything"""
        points_to_generate = self.points_target - self.n_points_run
        self.log.info('Starting run of %i points' %points_to_generate)

        self.finishscan=False
        noneload = 0 ## count of batches that yielded 0 new points after one another
        numbatch = 0 ## number of batches
        #testscores = []
        saved05 = False  # whether we've saved the model when it first hit below 5% on training set


        if self.epsilon != None:
            #lambda_lr = lambda numbatch: self.epsilon ** numbatch # this multiplies the learning ratio with epsion
            #scheduler = optim.lr_scheduler.LambdaLR(self.d_optimizer, lr_lambda=lambda_lr)
            scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(self.d_optimizer, 'min',factor=self.epsilon)

        while not self.finishscan and points_to_generate > 0 :
            numbatch = numbatch + 1
            # generate new points
            # the number of points might end up a little larger than the initial target, but it's better than having to deal with ridiculously small batches of new points
            if self.fromGood!=0:
                ### Choose 100% random points if the percentage error is >=50%, very few completely randoms if the error is small
                ### Also if we only have very few, generate plenty of points around whatever good points we have
                random_proportion=2.0*float(min(self.train_score,50))/100
                ### add penalty if the number of good points is not enough to be reliable
                if len(self.good_points) < 10:
                    random_proportion=max(random_proportion,0.75)

                num_random_points = int(self.K*random_proportion)
                
                num_KfromL_points=self.K-num_random_points
                self.log.info('Generating %d purely random points, %d via K from L' %(num_random_points,num_KfromL_points))

                if len(self.last_good_points) > 0:
                    num_random_from_good = int(num_random_points*0.8)
                    num_pure_random = num_random_points-num_random_from_good
                    randompointslist = [[np.random.uniform(0,1) for x in self.variables] for _ in range(num_random_from_good)]
                    goodpointslist = self.propose_fromGood(num_pure_random) # upscaled
                    new_random_vars = self.scale_points(randompointslist, down=False) + goodpointslist  # upscaled
                else:
                    randompointslist = [[np.random.uniform(0,1) for x in self.variables] for _ in range(num_random_points)]
                    new_random_vars = self.scale_points(randompointslist,down=False)  # upscaled

                if num_KfromL_points > 1:
                    num_fromgood_points = int(0.9*self.L)
                    extra_random_vars=[[np.random.uniform(0,1) for x in self.variables] for _ in range(int(0.1*(self.L)))]
                    new_temp_vars=self.scale_points(self.propose_fromGood(num_fromgood_points), down=True)
                    new_temp_vars = new_temp_vars + extra_random_vars  # downscaled
                    new_vars = self.scale_points(self.propose_KfromL_diverse(num_KfromL_points, self.L, points=new_temp_vars), down=False) # upscaled
                    new_vars= new_vars+new_random_vars  # upscaled
                else:
                    new_vars=new_random_vars
                
            else:
                new_vars = self.scale_points(self.propose_KfromL_diverse(self.K, self.L), down=False)  # upscaled

            new_res = self.generate_parameter_points(new_vars)
            self.n_points_run = self.n_points_run + len(new_vars)
            points_to_generate = points_to_generate - len(new_vars)

            ## classify & prepare training
            new_good_points=[]
            new_bad_points=[]
            for var, res in zip(new_vars, new_res):
                if res!=[] and res > 0.5:
                    new_good_points.append(var)
                elif res!=[]:
                    new_bad_points.append(var)

            self.last_good_points = new_good_points + self.last_good_points
            if len(self.last_good_points) > self.number_seed_points:
                self.last_good_points=self.last_good_points[:self.number_seed_points]

            self.good_points = self.good_points + new_good_points
            self.bad_points = self.bad_points + new_bad_points
            self.log.info('Adding %i good points and %i bad points to dataset before remixing & rebalancing.' %(len(new_good_points), len(new_bad_points)))
            ## if no new points are added repeatedly, warn & abort 
            if len(new_good_points)==0 or len(new_bad_points)==0: 
                noneload = noneload + 1
                self.log.info('Zero good or bad points to be added in try number %d.' % noneload)
                if noneload >= 100: 
                    self.log.info('Not getting more valid good & bad points. Please change your settings, e.g. variable ranges, fromGood value, etc. to change this in your next run.')
                    break

            ## train
            else:
                noneload=0 # even if noneload is different before, set it back to 0 because something is working anyway
                
            
                
            # rebalance initial train set with new points 
#            new_train_set = MyDataset(self.scalers,self.good_points,self.bad_points,1.0)
            new_train_set = MyDataset(self.scalers,new_good_points,new_bad_points,1.0) ## trains only on new points -> quicker, more impact of new points 
                
            # test dataset with new points
            self.log.info('Testing remixed & rebalanced dataset.')
            if len(new_train_set) > 100:
                self.log.debug('Length of train set to test: %d' % len(new_train_set))
                new_score = self.do_test(new_train_set)
                #if new_score < 5.: self.save_model()

                
            # train self.train_set
            self.log.info('Start training %i points.' % len(new_train_set))
            if len(new_train_set) < 10:
                self.train_score =100
            else:   
                self.log.debug('Length of train set: %d' % len(new_train_set))
                self.train_score=self.do_train(new_train_set, self.d_steps)
            self.log.info('Finished training.')

            # train on full set from time to time
            if self.fulltrain != 0 and numbatch % self.fulltrain == 0 and self.train_score > 5.0 and noneload < 21:
                self.log.info('Starting training on full dataset.')
                full_train_set = MyDataset(self.scalers,self.good_points,self.bad_points,1.0)
                temp_train_score = self.do_train(full_train_set, self.d_steps)
                if self.autostop == True:
                    belowfive = False
                    # save model if it's pretty good
                    if temp_train_score < 5.:
                        self.save_model()
                        belowfive = True
                    # abort if we're getting into shitty terrain
                    if temp_train_score > 20 and belowfive is True:
                        self.log.debug('Aborting. Model getting bad.')
                        self.log.debug('Consider tuning your model parameters.')
                        self.log.debug('If it was good before, a copy of this model has been saved.')
                        raise SystemExit('Model getting bad. Consider decreasing Epsilon or increasing ')

                self.log.info('Finished training.')

                
            new_points_added = 2* min(len(new_good_points),len(new_bad_points))  ## points effectively added after balancing 
            #self.n_points_run = self.n_points_run + new_points_added
            #points_to_generate = points_to_generate - new_points_added
            self.log.info('%d points left to be generated.' % points_to_generate)

            # test on benchmarking dataset, if so desired 
            if self.got_benchmark is True: # and numbatch%5:
                self.log.info('Testing benchmark now.')
                score = self.do_test(self.benchmarkset)
                if score < 5. and saved05 is False: 
                    self.save_model('05pc')
                    saved05 = True

                ## if it's not getting better on benchmark set, save & exit
                #worsttestscore = max(testscore[-10:-1])
                #if score > worsttestscore:
                #    self.log.debug('Aborting. Model getting bad.')
                #    self.log.debug('Consider tuning your model parameters.')
                #    self.log.debug('If it was good before, a copy of this model has been saved.')
                #    raise SystemExit('Model getting bad.')
                #testscores.append(score)


            self.d_optimizer.step() # this is needed for epsilon

            if self.epsilon != None: 
                scheduler.step(self.total_loss)
                self.log.debug('Batch # %i : Learning Rate %f' %(numbatch, self.d_optimizer.param_groups[0]['lr']))
                #self.log.debug('Batch # %i : Learning Rate %f' %(numbatch, scheduler.get_last_lr()[0]))


        # save discriminator
        self.save_model()
        #if self.modelfile is not None:
        #    torch.save(self.D.state_dict(), str(self.modelfile) + '.pt')
        #    torch.save(self.d_optimizer.state_dict(), str(self.modelfile) + '_opt.pt')

        # final test
        if self.got_benchmark is True:
            self.log.info('Testing on full test set.')
            new_score = self.do_test(self.benchmarkset)
            self.log.info('Test set: Error %.2f after %i points' %(new_score, self.n_points_run))

        self.log.info('Finished run of %d points' % self.n_points_run)

