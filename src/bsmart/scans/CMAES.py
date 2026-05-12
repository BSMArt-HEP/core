"""
Optimisation using CMAES


https://github.com/CyberAgentAILab/cmaes


BSMArt scan written by M. Goodsell


Requires:
        

pip3 install cmaes



"""

__meta__ = {
    "name": "CMAES",
    "requires": ["cmaes", "numpy"],
    "settings": {
        "Points": "Number of points",
        "CMAES MaxLoss": "Float",
        "CMAESEpisodes": "Int",
        "CMAESPopulationSize": "Int",
        "CMAESMaxGenerations": "Int",
        "CMAESMean": "Float",
        "CMAESSigma": "Float",
        "Sigma Tolerance": "Float",
        "Initial Sample": "Path to file"
    }
}



from bsmart.core import Scan as Scan

import numpy as np
import math
from bsmart import debug


import os
import shutil

from bsmart import BSMutil


from bsmart.BSMlikelihood import MakeLikelihoods, MakeGlobalLikelihood, safe_float

import sys, os

import matplotlib.pyplot as plt
import numpy as np

import cmaes

from cmaes import CMA

import math

class NewScan(Scan):
        """Scanner class for CMAES Scans"""

        def initialise(self):
                if self.runsettings.store_points_in_memory and self.runsettings.store_invalid_points:
                        self.naive = False 
                        self.runsettings.invalid_return_value = 0
                else:
                        self.naive = True  # we treat invalid points as bad
                        self.runsettings.invalid_return_value = []


        def __init__(self, inputs, log):
                Scan.__init__(self, inputs, log)
                

                self.downscalers,self.upscalers = BSMutil.create_scalers(self.inputs,log)
                self.n_variables=len(self.upscalers)
                #self.use_loglike=eval(self.inputs['Setup'].get('LogLike',"True"))
                self.Points=self.inputs['Setup'].get('Points',10000)
                
                self.maxloss =  np.log(1 + np.finfo(np.float64).max) + 1

                # Create the likelihood functions based on 'scaling' property of variables, if they have it; otherwise use gaussian likelihood    
                #elf.get_likelihood,self.observable_masks=MakeGlobalLikelihood(self.inputs['Observables'],self.use_loglike)
                #self.get_likelihood,self.observable_masks=MakeGlobalLikelihood(self.inputs['Observables'],return_type="NLL") # need a negative log likelihood
                self.likelihood_fns, self.observable_masks = MakeLikelihoods(self.inputs["Observables"], loglike=True)

                #self.likelihoodfns=MakeLikelihoods(self.inputs['Observables'],True) ## use log likelihood
                if "CMAES MaxLoss" in self.inputs["Setup"]:
                        if self.inputs["Setup"].get("CMAES MaxLoss") is not None:
                                self.maxloss = float(self.inputs["Setup"].get("CMAES MaxLoss"))
                #self.all_valid=[]
                #self.all_invalid=[]
                if 'CMAESEpisodes' in self.inputs['Setup']:
                        self.episodes=int(self.inputs['Setup'].get('CMAESEpisodes'))
                else:
                        #self.population_size=2*self.runsettings.cores
                        self.episodes = 1

                
                if 'CMAESPopulationSize' in self.inputs['Setup']:
                        self.population_size=int(self.inputs['Setup'].get('CMAESPopulationSize'))
                else:
                        #self.population_size=2*self.runsettings.cores
                        self.population_size = None

                if "CMAESMaxGenerations" in self.inputs["Setup"]:
                        self.max_n_generations = int(self.inputs["Setup"]["CMAESMaxGenerations"])
                else:
                        if self.population_size is not None:
                                self.max_n_generations = int(self.Points/self.population_size)
                        else:
                                self.max_n_generations = 100

                if 'CMAESMean' in self.inputs['Setup']:
                        #self.mean=int(self.inputs['Setup'].get('CMAESMean'))
                        self.mean=float(self.inputs['Setup'].get('CMAESMean'))*np.ones(self.n_variables)
                        self.fixed_mean = True
                else: # random mean
                        self.mean=np.random.rand(self.n_variables)
                        self.fixed_mean = False
                if 'CMAESSigma' in self.inputs['Setup']:
                        self.sigma=float(self.inputs['Setup'].get('CMAESSigma'))
                else: # use a reasonable default value
                        self.sigma=1.0/math.sqrt(self.n_variables)

                self.sigma_tolerance = 1e-4
                if 'Sigma Tolerance' in self.inputs['Setup']:  # to check for convergence 
                        self.sigma_tolerance = float(self.inputs['Setup']["Sigma Tolerance"])

                if 'Initial Sample' in self.inputs['Setup']:
                        from bsmart.ml.seed_points import get_previous_sample

                        seed_inputs,_ = get_previous_sample(sample_file=self.inputs['Setup']['Initial Sample'],
                                                inputs=self.inputs,
                                                nsamples = self.population_size,
                                                returntype =   "NLL")
                        """
                        seed_inputs,nlls = get_previous_sample(sample_file=self.inputs['Setup']['Initial Sample'],
                                                inputs=self.inputs,
                                                nsamples = self.population_size,
                                                returntype =   "NLL")
                        
                        print(seed_inputs.shape)
                        print(nlls.shape)
                        print(f"Minimum nll: {np.min(nlls)}")
                        print(np.hstack([seed_inputs, nlls.reshape(-1,1)]))  # nlls reports size as (n,), need to change to (n,1)
                        """

                        # now we have to scale the inputs:
                        #self.seed_inputs = np.array([[f(x) for x,f in zip(point,self.downscalers)] for point in seed_inputs])

                        # NB want to set the mean to be the best point, which is the first one in the list

                        self.mean = np.array([f(x) for x,f in zip(seed_inputs[0],self.downscalers)]) #self.seed_inputs[0]
                        # self.sigma=0.01 # less variation
                        # self.seed_inputs = None
                        self.seed_inputs = np.array([[f(x) for x,f in zip(point,self.downscalers)] for point in seed_inputs])
                else:
                        self.seed_inputs = None
                


        def postprocess(self,Point, observables, data_point,temp_dir,log, lock=None):
                """ return the likelihood; we won't get this far if the point failed to be generated """
                """ calculate likelihood """
                # For this algorithm we are minimising the objective function
                #res=1e100
                if self.naive:
                        return np.sum(self.get_losses(observables))
                else:
                        return 1.0
   
            

        def smooth_cap_loss(self,x):
                """
                Caps the loss by applying a sigmoid.
                This is useful for losses that are unbounded.                   

                """
                return -self.maxloss*np.expm1(x/self.maxloss) # assume x is negative and want to change its sign
   
        def get_losses(self,observables):
                """ Returns a list of losses.
                        
                The C-function loss should be zero if the observable is within allowed bounds, and greater than zero outside it.
                                If we include a validity condition, then other observables, evaluated afterwards, will be returned as NaN. 
                                        -> We should assign these the maximum loss ~ 710. 
                                
                                This is then compatible with the explicit use of 'hierarchical observables', where all other observables are set to the maximum.
                """

                likeit=iter(self.likelihood_fns)

                return [self.smooth_cap_loss((next(likeit))(val)) if mask and not math.isnan(val := safe_float(v)) else float((next(likeit) and False) or self.maxloss) 
                                        for v, mask in zip(observables, self.observable_masks) if mask]
                
        

        def run(self):
                """
                CMAES scan: just runs CMA-ES. It takes a minus-log-likelihood at the function to minimise.
                Excellent algorithm for finding best points

                
                
                """

                

                self.log.info("Starting CMA-ES scan!")
                bounds=np.zeros((self.n_variables,2))
                bounds[:,1] = np.ones(self.n_variables)
                #print(bounds)
                #optimizer=CMA(mean=np.zeros(self.n_variables),sigma=1.0)
                
                                #mean=0.5*np.ones(self.n_variables),sigma=0.3)
                """
                for generation in range(self.Points):
                        solutions = []
                        for _ in range(optimizer.population_size):
                                x = optimizer.ask()
                                value = loglike(x)
                                solutions.append((x, value))
                                upvals=[f(1.0/(1+math.exp(-x))) for x,f in zip(x,self.upscalers)]
                                print(f"#{generation} {value} ({x}) -> ({upvals})")
                        optimizer.tell(solutions)

                """
                for episode in range(self.episodes):
                        if not self.fixed_mean:
                                self.mean = np.random.rand(self.n_variables)
                        optimizer=CMA(bounds=bounds,
                                mean=self.mean,
                                sigma=self.sigma, 
                                population_size=self.population_size,
                                lr_adapt=True)
                        self.log.info('Starting CMAES run for episode '+str(episode+1))         

                        for generation in range(self.max_n_generations):
                                

                                batch=[]
                                cube_batch=[]
                                for _ in range(optimizer.population_size):
                                        x = optimizer.ask()
                                        cube_batch.append(x)
                                        #batch.append([f(1.0/(1+math.exp(-y))) for y,f in zip(x,self.upscalers)])
                                        batch.append([f(y) for y,f in zip(x,self.upscalers)])
                                results=self.RunManager.run_batch(batch)
                                #batch.append([f(1.0/(1+math.exp(-y))) for y,f in zip(x,self.upscalers)])
                                
                                #solutions = [(x,value) for x,value in zip(cube_batch,results)]
                                # check for bad values

                                if self.naive:
                                        solutions = [(x,value) if isinstance(value, (int, float)) else (x,self.maxloss) for x,value in zip(cube_batch,results)]
                                else: # keeping 'invalid' points
                                        batch_observables = [
                                                point["observables"] for point in self.RunManager.all_batch_points
                                                ]
                                        losses = np.array([self.get_losses(obs) for obs in batch_observables])
                                        fitnesses = np.sum(losses, axis=1)
                                        solutions=[(pop, value) for pop, value in zip(cube_batch, fitnesses)]

                                #print(solutions)
                                self.log.info(f"Losses for batch: {solutions}")
                                optimizer.tell(solutions)

                                stds = np.std(np.array(cube_batch),axis=0)
                                print(f"std: {stds} ")
                                current_sigma=max(stds)
                                #current_sigma = optimizer.sigma 
                                if any([optimizer.should_stop(), current_sigma < self.sigma_tolerance]):
                                        self.log.info(f"Converged afer {generation} generations")
                                        break
              
                self.log.info(f"Finished CMAES run")