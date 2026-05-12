"""
Machine Learning Scan


Ported from xBit by F. Staub [arXiv: 1906.03277](https://arxiv.org/abs/1906.03277)
Based on approach proposed by J. Ren,  L. Wu,  J.M. Yang,  J. Zhao [arXiv: 1708.06615](https://arxiv.org/abs/1708.06615)

This version written by M. D. Goodsell (goodsell@lpthe.jussieu) and F. Ibrahimov (Farid.Ibrahimov@liverpool.ac.uk)

Basic idea:
1. Start with random sample
2. Train NN
3. Use NN to propose new points with a 'good likelihood'

In BSMArt, the variables are specified with names and ranges, e.g.:

.. code-block:: json

    "Variables": { 
    		"m0": { "RANGE": [200,2500]},
    	    "m12": { "RANGE": [200,2500]}
    	 },


The scan will sample values within these ranges (by mapping from a hypercube)

The likelihood is computed from observables. Any observable that contains the key 'SCALING' will contribute, e.g. 
	"Observables": { "lres": {"SLHA": ["DUMMY",[1]], "SCALING": "USER"}}

or 
"mh" : { "SLHA": ["MASS", [25]], "SCALING": "GAUSS" , "MEAN": 125.09, "VARIANCE": 3.0}

The scan will *either* attempt to learn the likelihood based on the observables, 
*or* try to learn the values of the observables and use those to construct the likelihoods of proposed points. 

The settings for the scan are all specified in "Setup":

.. code-block:: text

    "Setup":{
            "RunName": "TestMLS",
            "Type": "MLS",
            "Cores": 4,
            "Output File": "MSSM_Output",
            "Spectrum File": "SPheno.spc.MSSM",
            "LR": 0.01,                 (learning rate)
            "Neurons": [50,50,50],
            "Points": 10,               (points per batch) 
            "Iterations": 30,           (number of batches)
            "Epochs" : 1000,            (number of epochs of training)
            "TrainLH": "True",          (Train Likelihood if True, learn observables if False)
            "Classifier": "False",      (include a classifier to filter out possibly invalid points)
            "LogLike": "True"           (Use log likelihood rather than the likelihood)
            "DensityPenaly": "True"     (Include a penalty for points too close to existing ones)
        },

"""

__meta__ = {
    "name": "MLS",
    "requires": ["torch", "numpy"],
    "settings": {
        "RunName": "String",
        "Cores": "Int",
        "Output File": "Path",
        "Spectrum File": "Path",
        "LR": "Float",
        "Neurons": "List[Int]",
        "Points": "Int",
        "Iterations": "Int",
        "Epochs": "Int",
        "TrainLH": "Boolean",
        "Classifier": "Boolean",
        "LogLike": "Boolean",
        "DensityPenaly": "Boolean"
    }
}



from bsmart.core import Scan as Scan
from bsmart.HEPRun import UpdateProgress, DataPoint
from bsmart.HEPRun import HEPRun
from bsmart.BSMlikelihood import MakeLikelihoods, MakeGlobalLikelihood
import bsmart.BSMlikelihood as BSMlikelihood
#from scans.aux import aux   ## was missing! Now copied into mls_NN.py


import bsmart.BSMutil as BSMutil   ## to have scalars and priors!

import bsmart.ml.mls_NN as aux

from bsmart.ml.mls_NN import NNMLS as nnmls

import numpy as np
import math
from bsmart import debug
import itertools
import queue
import os
import shutil
import multiprocessing as mp

import six
import queue

import random

import torch
import torch.nn as nn
import torch.optim as optim
from torch.autograd import Variable




class NewScan(Scan, nnmls):
    def __init__(self, inputs, log): 
        Scan.__init__(self, inputs, log)
        self.use_loglike=eval(self.inputs['Setup'].get('LogLike',"True")) ## use log likelihood by default
        #self.likelihoodfns=MakeLikelihoods(self.inputs['Observables'],self.use_loglike)

        self.get_likelihood,self.observable_masks=MakeGlobalLikelihood(self.inputs['Observables'],self.use_loglike)
        nnmls.__init__(self, inputs)    ## used in an odd way ...
        if "Points" in self.inputs['Setup']:
            self.points_target=int(self.inputs['Setup']['Points'])
        save_inputs = self.points_target
        
        ## scaling so that we specify points in terms of ranges
        self.downscalers,self.upscalers = BSMutil.create_scalers(self.inputs,log)
        
        self.n_variables=len(self.upscalers)

    def initialise(self):
        self.citations=r'''@article{Staub:2019xhl,
    author = "Staub, Florian",
    title = "{xBIT: an easy to use scanning tool with machine learning abilities}",
    eprint = "1906.03277",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    reportNumber = "KA-TP-10-2019",
    month = "6",
    year = "2019"
}

@article{Ren:2017ymm,
    author = "Ren, Jie and Wu, Lei and Yang, Jin Min and Zhao, Jun",
    title = "{Exploring supersymmetry with machine learning}",
    eprint = "1708.06615",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    doi = "10.1016/j.nuclphysb.2019.114613",
    journal = "Nucl. Phys. B",
    volume = "943",
    pages = "114613",
    year = "2019"
}'''

    """
    def get_likelihood(self,observables):
        lh = 0.
		
        if self.use_loglike:
            for  v, f in zip(observables,self.likelihoodfns):
            
                lh = lh + f(float(v))
        else:
            lh=1.0
            for  v, f in zip(observables,self.likelihoodfns):
                lh = lh * f(float(v))
        return lh
    """
    def postprocess(self,Point, observables, data_point,temp_dir,log, lock=None):
        """ 
            If we are learning the likelihood we want to return just the likelihood
            If we are learning the *observables* then we want to return just them
        
        """

        if self.train_lh:
            return self.get_likelihood(observables)
            #lh = self.get_likelihood(observables)
            #return lh   
        else:
            return observables
    

    def generate_parameter_points(self, num_points):
        """
            Generate random points in a unit hypercube
            Can use torch rand for this!
        """
        return torch.rand([num_points,self.n_variables])
   
    def guess_LH(self, x_in):

        if self.train_lh:     
            y_try = self.predictor(Variable(torch.FloatTensor(x_in))).data.cpu().numpy()
        else:
            y_try = self.rescale_data(self.predictor(Variable(torch.FloatTensor(x_in))).data.cpu().numpy())
        

        # time for some classifiers, we could use classifierdataset from ml.py but for now do what is in xBit
        if self.use_classifier:
            y_class = self.classifier(Variable(torch.FloatTensor(x_in))).data.cpu().numpy()
        else: # Not using classifier, so define all points as passing. This is inefficient.
            y_class = [[0.99, 0.1] for _ in y_try]
            

        #filter valid points and find likelihood
        lh_valid = []
        x_valid = []

        for y, x, yc in zip(y_try, x_in, y_class):
            if (yc[0] > yc[1]):
                # If train_lh learn likelihood
                if self.train_lh:
                    lh_valid.append(y)
                else:
                    # otherwise we need to guess the likelihood from the observables
                    lh_valid.append(aux.make_safe_lh(self.get_likelihood(y)))

                x_valid.append(x)

        # Calculate final likelihood (optionally with distance penalty)
        if not self.distance_penalty:
            likelihood = lh_valid
        else:
            likelihood = self.min_Delta_LH(lh_valid, x_valid)

        if self.use_loglike:
            likelihood = [aux.exp_safe(l) for l in likelihood] # we actually only use this for comparison with target, never for training etc
        return likelihood, x_valid
    
    def min_Delta_LH(self, lh, xval):    # implements the distance penalty

        #all_points = np.array([x[0] for x in self.sample])
        all_points = self.sample[0].cpu().numpy()
        lh_out = []
        norm=1  # have rescaled distances onto the hypercube

        # This will penalise points that are too close to existing ones
        if self.use_loglike:  # if we use log likelihoods need to take log of the penalty too!
            for l,x in zip(lh, xval):
                min_d = np.amax(np.abs((all_points - np.array(x))/norm), axis = 1).min()
                lh_out.append(l+2.0*math.log(min_d))
        else:
            for l,x in zip(lh, xval):
                min_d = np.amax(np.abs((all_points - np.array(x))/norm), axis = 1).min()
                lh_out.append(l*min_d**2)
        return lh_out
    
    def new_points(self, good_points, points):
        # estimate new points with good likelihood
        self.log.info("new points initiated")
        self.predictor.eval()
        
        """
            Assumes that the likelihood is positive and bounded above by 1. Larger values are better.
        """
       
        
        best_lh=[]
        increase=1
        while (len(best_lh)<5):
            x_try = self.generate_parameter_points(max(increase*points,100000))  # silly but what is in xBit
            #x_try = self.generate_parameter_points(min((increase**2)*points,100000))
            likelihood, x_try = self.guess_LH(x_try)
            best_lh = [l for l in likelihood if l < 1 and l > 0]
            increase += 1
            if increase > 100:
                break
        if len(best_lh) > 0:
            
            best_lh = max(best_lh)
            count = 1
            count_all = 0
            self.log.info("max LH proposed by NN: %s; " % (str(best_lh)))
            while count < int(0.9 * good_points) and count_all < 10: ## don't want arbitrary cycles
                count_all = count_all + 1            
                for l, x in zip(likelihood, x_try):
                    if 2. > l > aux.limit_lh(count_all) * best_lh:  
                        # 2 > might be necessary for the LH fit because some good points are 
                        # predicted to have slightly higher LH;

                        if count == 1:
                            new_good_points=x.reshape((1,-1))
                        else:
                            new_good_points = torch.cat((new_good_points,x.reshape((1,-1))),0)
                        count = count + 1
                        if count > int(0.9 * good_points):
                            break        

                if count < int(0.9 * good_points): # don't need to do this if we reached the number
                    self.log.info('Generating more points, currently have %d out of %d' %(count,int(0.9*good_points)))
                    x_try = self.generate_parameter_points(max(increase*points,100000))            
                    likelihood, x_try = self.guess_LH(x_try)

            # add 10% random points
            new_x = torch.cat((self.generate_parameter_points(int(0.1 * good_points)),new_good_points),0)
        else:
            # no valid points found by the NN! That's likely an overconstraint from the classifier. 
            # in order to contue: let's take 100% random points
            self.log.info('No interesting points predicted!')
            new_x = self.generate_parameter_points(good_points) 

        return new_x


    def upscale(self,input_torch):
        """
            Creates a list of inputs to feed to run_batch when given a tensor of dim (*, num_variables)
        """
        
        return [[f(float(x)) for x,f in zip(p,self.upscalers)] for p in input_torch]

    def call_batch(self,input_torch):
        """
            Handles running of the points

            input_torch is a torch tensor
        """
        upscaled=[[f(float(x)) for x,f in zip(p,self.upscalers)] for p in input_torch]

        results=self.RunManager.run_batch(upscaled)
        total_pts=input_torch.size(0)
        new_sample=[None,[]]
        mask=[]
        for rvars,res in zip(input_torch,results):
                if res!=[]:
                    self.all_valid.append(list(map(float,rvars)))
                    mask.append(True)
                    # results are always a list, but returned as strings for genericity purposes
                    # So map to floats

                    #new_sample[1].append(list(map(float,res)))
                    # only keep relevant observables
                    new_sample[1].append([float(r) for r,om in zip(res,self.observable_masks) if om ]) 
                else:
                    mask.append(False)
                    self.all_invalid.append(list(map(float,rvars)))
        mask_tensor = torch.tensor(mask, dtype=torch.bool)
        new_sample[0] = input_torch[mask_tensor]
        self.sample[0] = torch.cat((self.sample[0],new_sample[0]),dim=0)
        self.sample[1].extend(new_sample[1])

    def run(self):
        
        save_points = self.inputs['Setup']['Points']

        self.all_valid=[]
        self.all_invalid=[]
        self.sample= [torch.empty(0,self.input_dim),[]]
        self.log.info("Starting MLS scan with initial batch of %d points" % (5*save_points))
        initial_points = self.generate_parameter_points(5*save_points) 
        self.call_batch(initial_points)


        for i in range(1, self.iterations):
            self.log.info('Running iteration %d of %d' %(i,self.iterations))
            # self.sample should only contain valid points
            x = self.sample[0]
            """
                If we have set train_lh then we learn the likelihood:
                    In this case the results of postprocessing should just be the likelihood.
                    We have "make_safe_lh" to limit its range, curbing very large values via a tanh function

                If train_lh is False then we learn the observables:
                    In this case postprocessing returns a list of observables. 
                    We scale these so that we don't have large variations
            """

            if self.train_lh:
                y_scale = [aux.make_safe_lh(float(xx)) for xx in self.sample[1]]
            else:          
                y_scale = self.scale_data(self.sample[1])

            self.set_predictor()
            self.do_train(x, y_scale, self.predictor)

            # Include a classifier to distinguish valid/invalid points
            if self.use_classifier:
                self.set_classifier()
                x_class = self.all_valid + self.all_invalid[0:min([len(self.all_valid),len(self.all_invalid)-1])] 
                y_class = [[0.99, 0.01] for j in self.all_valid] + \
                          [[0.01, 0.99] for j in self.all_invalid[0:min([len(self.all_valid),len(self.all_invalid)-1])]]
                self.do_train(x_class, y_class,self.classifier)

            proposed_points = self.new_points(self.inputs['Setup']['Points'], len(self.inputs['Variables'])**2 * self.inputs['Setup']['Points'])
            self.call_batch(proposed_points)

