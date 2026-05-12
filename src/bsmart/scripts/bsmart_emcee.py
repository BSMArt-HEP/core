#!/usr/bin/env python3

import argparse
import os
from datetime import datetime
import sys

import shutil
import json
import collections

import importlib

# Package imports
from bsmart import debug
from bsmart.BSMArtBlackBox import get_blackbox
from bsmart.ml.seed_points import get_previous_sample
from bsmart import BSMutil
import emcee
import time
import multiprocessing as mp

# Thanks to python 3.14 ...
try:
    if mp.get_start_method() != 'fork':
        mp.set_start_method('fork', force=True)
except (ValueError, RuntimeError): 
    pass

import numpy as np

# Global variables handling inside logic or passing
RunManager = None
upscale_handler = None
likelihood_handler = None

def emcee_init_workers(shared_map):
    pid = os.getpid()

    # Assign a worker index if not already assigned
    if pid not in shared_map:
        shared_map[pid] = len(shared_map)
    #print(shared_map)

def emcee_log_prob(theta,shared_id,shared_map,mylock):
    global RunManager,upscale_handler,likelihood_handler
    pid = os.getpid()

    # first check if any value is out of range: we won't even run it if that is the case

    if any([x > 1.0 for x in list(theta)]) or any ([x < 0.0 for x in list(theta)]):
        return -1e100

    # Assign a worker index if not already assigned
    with mylock:
        if pid not in shared_map:
            shared_map[pid] = len(shared_map)
        #print(shared_map)
        worker_number=shared_map[pid]
    
        shared_id.value+=1
    try:
        scaled_up=upscale_handler.upscale(list(theta))
        #print(f"scaled up: {scaled_up}")
        RunManager.CoreRuns[worker_number].last_point_dict= {}
        res=RunManager.CoreRuns[worker_number].run_point(upscale_handler.upscale(list(theta)),shared_id.value,mylock)
        #print(f"Worker {worker_number} ran correctly and returned {res}")
    except Exception as e:
        #print(f"Worker {worker_number} failed with {e}")
        if likelihood_handler is not None:
            try:
                res = likelihood_handler.get_LL(RunManager.CoreRuns[worker_number].last_point_dict["observables"])
                #print(f"Worker {worker_number} returned {res}")
            except Exception as e:
                #print(f"Worker {worker_number} failed with {e}")
                res = -1e100
        else:
            #print("No likelihood handler!")
            res = -1e100
    #res=0.5
    #print(f"Worker {worker_number} returned {res}")
    return res


class UpScalerHandler:
    """ 
    This class is used to upscale the parameters to the range of the black box, required in order to pickle.
    """
    def __init__(self,inputs,log):
        _,self.upscalers = BSMutil.create_scalers(inputs,log)
    def upscale(self,theta):
        return [f(y) for y,f in zip(list(theta),self.upscalers)]

# class LikelihoodHandler:
#     def __init__(self,inputs,log):
#         from bsmart.BSMlikelihood import MakeLikelihoods, safe_float
#         self.maxloss = sys.float_info.max
#         self.likelihood_fns, self.observable_masks = MakeLikelihoods(inputs["Observables"], loglike=True)

#     def smooth_cap_loss(x):
#             """
#             Caps the loss by applying a sigmoid.
#             This is useful for losses that are unbounded.                   

#             """
#             return self.maxloss*math.expm1(x/self.maxloss) # 


#     def get_NLL(self,observables):
#         """ return the likelihood; we won't get this far if the point failed to be generated """
                    
#         likeit=iter(self.likelihood_fns)

#         return -1.0*math.sum([self.smooth_cap_loss((next(likeit))(val)) if mask and not math.isnan(val := safe_float(v)) else float((next(likeit) and False) or self.maxloss) 
#                                             for v, mask in zip(observables, self.observable_masks) if mask])

from bsmart.BSMlikelihood import MakeLikelihoods, safe_float

class LikelihoodHandler:
    def __init__(self,inputs,log):
        self.maxloss =  -(np.log(1 + np.finfo(np.float64).max) + 1)
        self.likelihood_fns, self.observable_masks = MakeLikelihoods(inputs["Observables"], loglike=True)

    def smooth_cap_loss(self,x):
            """
            Caps the loss by applying a sigmoid.
            This is useful for losses that are unbounded.                   

            """
            return -self.maxloss*np.expm1(-x/self.maxloss) # 


    def get_LL(self,observables):
        """ return the likelihood; we won't get this far if the point failed to be generated """
                    
        likeit=iter(self.likelihood_fns)

        return np.sum([self.smooth_cap_loss((next(likeit))(val)) if mask and not np.isnan(val := safe_float(v)) else float((next(likeit) and False) or self.maxloss) 
                                            for v, mask in zip(observables, self.observable_masks) if mask])


def main():
    global RunManager,upscale_handler,likelihood_handler
    try:
        parser = argparse.ArgumentParser(
            description='Please give the name of the input file.')
        parser.add_argument('inputfile',
                            metavar='File', type=str,
                            help='Input file name')
        # Are --short and --csv really necessary any more? 
        parser.add_argument("--short", help="Store output in short tabbed form",
                            action="store_true")
        parser.add_argument("--csv", help="Store output in csv form",
                            action="store_true")
        parser.add_argument("--debug", help="write debug information",
                            action="store_true")
        parser.add_argument("--NoMPI", help="Do not check for MPI",
                            action="store_true")
        """
            New feature to help users (not yet implemented) 
        """
        parser.add_argument("--Settings", help="List all settings for your selected scan and tools!",
                            action="store_true")
        args = parser.parse_args()
    except:
        print("No input file provided\n")

        raise SystemExit

    ### Now load in the input file
    try:
        RunManager,inputs=get_blackbox(args.inputfile,vars(args),likelihood_return_type='LL')
    except Exception as e:
        print('Failed to load input file!' + str(e))
        raise SystemExit


    #RunManager=bsmart.HEPRun(self.runsettings, self.write_lh_file, self.postprocess,self.log)

    downscalers,upscalers = BSMutil.create_scalers(inputs,RunManager.log)
    upscale_handler=UpScalerHandler(inputs,RunManager.log)


    if RunManager.settings.store_points_in_memory and RunManager.settings.store_invalid_points:
            
            naive = False 
            RunManager.settings.invalid_return_value = 0
            likelihood_handler=LikelihoodHandler(inputs,RunManager.log)
    else:
            
            naive = True  # we treat invalid points as bad
            RunManager.settings.invalid_return_value = []
            likelihood_handler=None


    ### EMCEE related settings

    nwalkers=40
    ndim=len(inputs['Variables'])
    nsteps = 100



    if 'Steps' in inputs['Setup']:
        nsteps = int(inputs['Setup']['Steps'])


    ncores=1
    if 'Cores' in inputs['Setup']:
        ncores=int(inputs['Setup']['Cores'])

    if 'Walkers' in inputs['Setup']:
        nwalkers = int(inputs['Setup']['Walkers'])
    else: # use a large number of walkers
        nwalkers = ncores * 10


    if 'Initial Sample' in inputs['Setup'] and os.path.exists(inputs['Setup']['Initial Sample']):
        unscaled_p0,_=get_previous_sample(inputs['Setup']['Initial Sample'],inputs,nwalkers)

        p0 =np.array([[f(y) for y,f in zip(pt,downscalers)] for pt in unscaled_p0])
        
        #print(f"Starting with {p0}")
    else:
        p0 = np.random.rand(nwalkers, ndim)
    #print(f"p0: {p0.shape}, nwalkers: {nwalkers}")
    #print(f"p0: {p0}")

    with mp.Manager() as manager:
        shared_map = manager.dict()
        shared_id = manager.Value('i', 0)
        lock=manager.Lock()
        
        
        with mp.Pool(processes=ncores, initializer=emcee_init_workers, initargs=(shared_map,)) as pool:
            
            #sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob, pool=pool)
            sampler = emcee.EnsembleSampler(nwalkers, ndim, emcee_log_prob, pool=pool,args=(shared_id,shared_map,lock))
            start = time.time()
            sampler.run_mcmc(p0, nsteps, progress=True)
            end = time.time()
            multi_time = end - start
            print("Multiprocessing took {0:.1f} seconds".format(multi_time))
            #print("{0:.1f} times faster than serial".format(serial_time / multi_time))
            #samples = sampler.get_chain(flat=True)
            samples = sampler.get_chain()  ## this will be of dimension (nsteps, nwalkers, nvars)
            allsamples=sampler.get_chain(flat=True)   # this is of dimension (nsteps*nwalkers, nvars)

            ## Note
            #[int(f"Retained {len(samples)} samples out of {len(allsamples)}")
            # NB if we don't advance, then the parameters in the chain are repeated.
            print(f"Collected {len(allsamples)} samples")
            print("Mean acceptance fraction: {0:.3f}".format(np.mean(sampler.acceptance_fraction)))
            
            # Now we need to write a results file with the chain
            varlist=list(inputs['Variables'].keys())
            #np.savetxt("output.csv", data, delimiter=',', fmt='% s', comments='', header=header)

            resultsdir=os.path.join(inputs['Setup']['cwd'], inputs['Setup']['RunName'],'Results')
            if os.path.exists(resultsdir):
                shutil.rmtree(resultsdir)
            os.makedirs(resultsdir)

            upscaled_samples=np.array([[f(y) for y,f in zip(list(theta),upscalers)] for theta in allsamples])

            np.savetxt(os.path.join(resultsdir,"Results.csv"), upscaled_samples, delimiter=',', comments='', header=','.join(varlist))


            plotdir=os.path.join(inputs['Setup']['cwd'], inputs['Setup']['RunName'], 'Plots')
            if os.path.exists(plotdir):

                os.chdir(plotdir)
                from bsmart import BSMplots
                # Can't do this as don't have a scan object
                #BSMplots.make_plots()

                ## Make an extra corner plot 

                plotstub=os.path.join(resultsdir,"corner_plot")
                BSMplots.make_auto_corner_csv(plotstub,os.path.join(resultsdir,"Results.csv"),varlist,{})
                os.system('python3 '+plotstub+'.py')

                plotstub=os.path.join(resultsdir,'fancy_auto_corner')
                BSMplots.make_fancy_auto_corner_csv(plotstub,os.path.join(resultsdir,"Results.csv"),varlist,{})
            
                os.system('python3 '+plotstub+'.py')

if __name__ == "__main__":
    main()
