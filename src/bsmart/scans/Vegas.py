"""

Importance Sampling using vegas

Written by M. Goodsell



"""


from bsmart.core import Scan as Scan

import numpy as np
import math
from bsmart import debug


import os
import shutil

from bsmart import BSMutil


#from bsmart.scans import MCMC
import copy

from bsmart.BSMlikelihood import MakeLikelihoods, MakeGlobalLikelihood

import sys, os

import numpy as np

import vegas




class NewScan(Scan):
        """Scanner class for Vegas Scans"""
        
        def __init__(self, inputs, log):
                Scan.__init__(self, inputs, log)
                

                self.downscalers,self.upscalers = BSMutil.create_scalers(self.inputs,log)
                self.n_variables=len(self.upscalers)
                self.use_loglike=eval(self.inputs['Setup'].get('LogLike',"True"))
                self.Points=int(self.inputs['Setup'].get('Points',10000))
                """Number of points to sample. Default: 10000."""
                # batch size for training
                self.batch_size = int(self.inputs['Setup'].get('Batch Size',100))
                if 'Iterations' in self.inputs['Setup']:
                        self.batch_size=int(self.Points/self.inputs['Setup']['Iterations'])
                # Create the likelihood functions based on 'scaling' property of variables, if they have it; otherwise use gaussian likelihood    
                self.get_likelihood,self.observable_masks=MakeGlobalLikelihood(self.inputs['Observables'],self.use_loglike)
                """
                self.likelihoodfns=MakeLikelihoods(self.inputs['Observables'],True) ## use log likelihood
                """
                #self.input_dim=sum([1 for x in self.observable_masks if x])
                #self.flow=RealNVP(dim=self.n_variables, hidden_dim=64, num_layers=4).to(device)


                self.epochs = int(self.inputs["Setup"].get("Epochs",1000))
                self.LR=float(self.inputs["Setup"].get("LR",1e-3))
                self.patience=int(self.inputs["Setup"].get("Patience",30))
                
                self.all_valid=[]
                self.all_invalid=[]
                #self.sample= [torch.empty(0,self.input_dim),[]]
                # Only 
                #self.sample= [torch.empty(0,self.n_variables),torch.empty(0)]


        def postprocess(self,Point, observables, data_point,temp_dir,log, lock=None):
                """ return the likelihood; we won't get this far if the point failed to be generated """
                return self.get_likelihood(observables)
                #lik =self.get_likelihood(observables)
                #return lik
                
        def call_batch(self,input_numpy):
                """
                Handles running of the points

                input_torch is a torch tensor
                """
                upscaled=[[f(float(x)) for x,f in zip(p,self.upscalers)] for p in input_numpy]

                results=self.RunManager.run_batch(upscaled)
               
                mask=[]
                new_results=[]
                for rvars,res in zip(input_numpy,results):
                        if res!=[]:
                                self.all_valid.append(list(map(float,rvars)))
                                mask.append(True)
                                new_results.append(res)
                        else:
                                mask.append(False)
                                new_results.append(0.0)
                                self.all_invalid.append(list(map(float,rvars)))
                                
                return np.array(new_results)
                

        def run(self):



                self.log.info("Starting Vegas scan!")
                
                self.log.info(f"batch size: {self.batch_size}, points: {self.Points}")
                
                @vegas.lbatchintegrand
                def f(x):
                        return self.call_batch(x)
                
                # create a box that is [0,1]**n_variables

                # dls_vegas.neval_batch=100
                # dls_vegas.min_neval_batch=100
                integ = vegas.Integrator(self.n_variables * [[0, 1]])
                warmup = integ(f, nitn=2, neval=200)
                result = integ(f, nitn=5, neval=1000)
                print(result.summary())
                print('result = %s    Q = %.2f' % (result, result.Q))