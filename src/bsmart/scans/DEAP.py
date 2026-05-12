"""
Optimisation using genetic algorithm DEAP





BSMArt scan written by M. Goodsell


Requires:
        

pip3 install deap

https://github.com/DEAP/deap

"""

__meta__ = {
    "name": "DEAP",
    "requires": ["deap", "numpy"],
    "settings": {
        "Points": "Int",
        "DEAP": { 
                "Population Size": "Int", 
                "Generations": "Int", 
                "Global Only": "Bool"

        },
        "n_cores": "Int"
    }
}


# NB from deap/base.py the fitness has:
#     @property
#     def valid(self):
#         return len(self.wvalues) != 0




# from deap import creator, base, tools, algorithms

# creator.create("FitnessMax", base.Fitness, weights=(1.0,))
# creator.create("Individual", list, fitness=creator.FitnessMax)

# toolbox = base.Toolbox()

# toolbox.register("attr_bool", random.randint, 0, 1)
# toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.attr_bool, n=100)
# toolbox.register("population", tools.initRepeat, list, toolbox.individual)

# def evalOneMax(individual):
#     return sum(individual),

# toolbox.register("evaluate", evalOneMax)
# toolbox.register("mate", tools.cxTwoPoint)
# toolbox.register("mutate", tools.mutFlipBit, indpb=0.05)
# toolbox.register("select", tools.selTournament, tournsize=3)

# population = toolbox.population(n=300)

# NGEN=40
# for gen in range(NGEN):
#     offspring = algorithms.varAnd(population, toolbox, cxpb=0.5, mutpb=0.1)
#     fits = toolbox.map(toolbox.evaluate, offspring)
#     for fit, ind in zip(fits, offspring):
#         ind.fitness.values = fit
#     population = toolbox.select(offspring, k=len(population))
# top10 = tools.selBest(population, k=10)




from bsmart.core import Scan as Scan

import numpy as np
import math
from bsmart import debug


import os
import shutil

from bsmart import BSMutil


from bsmart.BSMlikelihood import MakeLikelihoods, MakeGlobalLikelihood, safe_float

import sys, os

#import matplotlib.pyplot as plt
import numpy as np
import random

from deap import creator, base, tools, algorithms



import math

class NewScan(Scan):
        """Scanner class for NF Scans"""
        def initialise(self): # overload initial settings. 
                self.runsettings.store_points_in_memory = True # Force storage of info for hierarchical loss computation
                self.runsettings.store_invalid_points = True  # Need this to keep points where 
                self.runsettings.invalid_return_value = 0


        def __init__(self, inputs, log):
                Scan.__init__(self, inputs, log)
                

                self.downscalers,self.upscalers = BSMutil.create_scalers(self.inputs,log)
                self.n_variables=len(self.upscalers)
                #self.use_loglike=eval(self.inputs['Setup'].get('LogLike',"True"))
                self.Points=self.inputs['Setup'].get('Points',10000)
                
                
                # Create the likelihood functions based on 'scaling' property of variables, if they have it; otherwise use gaussian likelihood    
                #self.get_likelihood,self.observable_masks=MakeGlobalLikelihood(self.inputs['Observables'],self.use_loglike)
                #self.get_likelihood,self.observable_masks=MakeGlobalLikelihood(self.inputs['Observables'],loglike=True)
                self.likelihood_fns,self.observable_masks=MakeLikelihoods(self.inputs['Observables'],loglike=True)
                #self.likelihoodfns=MakeLikelihoods(self.inputs['Observables'],True) ## use log likelihood
                self.live_observables = sum(1 if mask else 0 for mask in self.observable_masks)
                #self.all_valid=[]
                #self.all_invalid=[]
                self.pop_size = max(10,self.inputs['Setup'].get('n_cores',1))
                self.NGEN = 40
                self.global_only = False
                if 'DEAP' in self.inputs['Setup']:
                        self.pop_size = self.inputs['Setup']['DEAP'].get('Population Size',self.pop_size)
                        self.NGEN = self.inputs['Setup']['DEAP'].get('Generations',self.NGEN)
                        self.global_only = self.inputs['Setup']['DEAP'].get('Global Only',self.global_only)
                
                #self.maxloss = 1e100
                self.maxloss = 1e8

        def postprocess(self,Point, observables, data_point,temp_dir,log, lock=None):
                """ return the likelihood; we won't get this far if the point failed to be generated """
                """ calculate NLL -> minimise this """
                # For this algorithm we are minimising the objective function
                # 

                return 1.0

                """
                res=1e100
                #if self.use_loglike:
                if True:
                        try:
                                res = -self.get_likelihood(observables)
                        except:
                                pass
                """
                """
                else:
                        try:
                                res = -math.log(self.get_likelihood(observables))
                        except:
                                pass
                """
                #print(f"Res? {res}")
                # return "weights"; this should be a tuple, but we sort this later
                return res #-math.log(self.get_likelihood(observables))

                #return self.get_likelihood(observables)
   
        def get_losses(self,observables):
                """ Returns a list of losses.
                
                The C-function loss should be zero if the observable is within allowed bounds, and greater than zero outside it.
                        If we include a validity condition, then other observables, evaluated afterwards, will be returned as NaN. 
                                -> We should assign these the maximum loss ~ 710. 
                        
                        This is then compatible with the explicit use of 'hierarchical observables', where all other observables are set to the maximum.
                """

                likeit=iter(self.likelihood_fns)
                #res = [(f := next(likeit))(val) if mask and not math.isnan(val := safe_float(v)) else ((next(likeit) and False) or self.maxloss if mask else self.maxloss)
    		#		for v, mask in zip(observables, self.observable_masks)]
                #print(f"Result: {res}")
                #return res
                #return [(f := next(likeit))(val) if mask and not math.isnan(val := safe_float(v)) else ((next(likeit) and False) or self.maxloss if mask else self.maxloss)
    		#		for v, mask in zip(observables, self.observable_masks)]
                #res = [(f := next(likeit))(val) if mask and not math.isnan(val := safe_float(v)) else (next(likeit) and False) or self.maxloss 
    		#		for v, mask in zip(observables, self.observable_masks) if mask]
                #print(f"These obs: {[obs for obs,mask in zip(observables,self.observable_masks) if mask]}")
                #print(f"Give losses: {res}")
                #return res
                

                return tuple([(f := (-1.0/self.maxloss)*(next(likeit))(val)) if mask and not math.isnan(val := safe_float(v)) else float((next(likeit) and False) or 1.0) 
    				for v, mask in zip(observables, self.observable_masks) if mask])
                #return tuple([(f := (-1.0)*(next(likeit))(val)) if mask and not math.isnan(val := safe_float(v)) else float((next(likeit) and False) or self.maxloss) 
    		#		for v, mask in zip(observables, self.observable_masks) if mask])


        

        def run(self):
                """
                DEAP scan
                
                Version 1: just use it as an optimisation over a single global likelihood
                @todo: allow multiple likelihoods for different observables, with hierarchical evaluation
                
                
                """
                

                self.log.info("Starting DEAP scan!")

                if self.global_only:
                        creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
                else: #multi weight
                        creator.create("FitnessMin", base.Fitness, weights=(-1.0,)*self.live_observables)

                creator.create("Individual", list, fitness=creator.FitnessMin)

                self.toolbox = base.Toolbox()
                

                #self.toolbox.register("attr_bool", random.randint, 0, 1)
                self.toolbox.register("random_float", random.random)
                self.toolbox.register("individual", tools.initRepeat, creator.Individual, self.toolbox.random_float, n=self.n_variables)
                self.toolbox.register("population", tools.initRepeat, list, self.toolbox.individual)

                #def evalOneMax(individual):
                #        return sum(individual),

                #self.toolbox.register("evaluate", evalOneMax)
                #self.toolbox.register("mate", tools.cxTwoPoint)
                self.toolbox.register("mate", tools.cxBlend,alpha=0.2)  # better for floats
                # mutFlipBit is usually for bools
                #self.toolbox.register("mutate", tools.mutFlipBit, indpb=0.05)
                # deap.tools.mutGaussian(individual, mu, sigma, indpb)
                self.toolbox.register("mutate",tools.mutGaussian, mu=0.0, sigma=0.2, indpb=1.0/math.sqrt(self.n_variables))
                
                self.toolbox.register("select", tools.selTournament, tournsize=3)

                population = self.toolbox.population(n=self.pop_size)

                #NGEN=40
                #self.NGEN
                for gen in range(self.NGEN):
                        # cxpb = mating prob, mutpb = mutation prob
                        offspring = algorithms.varAnd(population, self.toolbox, cxpb=0.4, mutpb=0.05)
                        #print(offspring)
                        batch=[]
                        cube_batch=[]
                        for child in offspring:
                                cube_batch.append(child)
                                #batch.append([f(1.0/(1+math.exp(-y))) for y,f in zip(x,self.upscalers)])
                                batch.append([f(y) for y,f in zip(child,self.upscalers)])
                        print(f"Batch: {batch}")
                        
                        results=self.RunManager.run_batch(batch)
                        #print(self.RunManager.all_batch_points)
                        #full_losses = [(x,value) if isinstance(value, (int, float)) else (x,1e100) for x,value in zip(cube_batch,results)]
                        batch_observables=[point['observables'] for point in self.RunManager.all_batch_points]
                        keep_observables=[[obs for obs,mask in zip(point_obs,self.observable_masks) if mask] for point_obs in batch_observables]
                        #print(f"Got observables: {batch_observables}")
                        weights = [self.get_losses(obs) for obs in batch_observables]

                        for fit, ind,result in zip(weights, offspring,results):
        
                                ind.fitness.values = fit
                        print(f"weights: {weights}")
                        #weights = [(value,) if isinstance(value, (int, float)) else (1e100,) for value in results]
                        # """
                        # valids = [True  if isinstance(value, (int, float)) else False for value in results]
                        # weights = [(value,) if valid else (1e100,) for value,valid in zip(results,valids)]
                        # print(f"Results: {results}")
                        
                        # for fit, ind,valid in zip(weights, offspring,valids):
                        #         if valid:
                        #                 ind.fitness.values = fit
                        # """
                        # """
                        #         #The valid property is set when we give fit values
                        #         if ind.fitness.valid:
                        #                 print(f"{fit} ->valid")
                        #         else:
                        #                 print(f"{fit} -> INVALID")
                        # """
                        population = self.toolbox.select(offspring, k=len(population))
                top10 = tools.selBest(population, k=10)
                self.log.info(f"Best points: {top10}")

                
                
                
