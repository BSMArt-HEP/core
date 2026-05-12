"""

Affine MCMC scan using the Goodman & Weare algorithm.

This is an excellent robust MCMC algorithm for larger numbers of dimensions. 

The chain, and automatic plots, are stored in a Results subdirectory. Results are updated after each pass. 

The relevant scan-specific settings are:

.. code-block:: text

    "Setup" : {
        "Steps": "int, Number of steps",
        "Walkers" : "int, Number of walkers, default 10 * number of variables"
    }





To store observable information inside the results file, it is necessary to set:

.. code-block:: json

    "Setup" : {
        "store_points_in_memory": "True",
        "store_invalid_points": "True"
    }

"""

__meta__ = {
    "name": "AffineMC",
    "requires": ["matplotlib","numpy","corner","pandas","seaborn"],
    "settings": {
        "Steps": "int, Number of steps",
        "Walkers" : "int, Number of walkers, default 10 * number of variables"
    }
}



from bsmart.core import Scan as Scan

import numpy as np
import math
from bsmart import debug


import os

from bsmart import BSMutil


from bsmart.BSMlikelihood import MakeGlobalLikelihood

import sys, os

import numpy as np



import math

class NewScan(Scan):
        """Scanner class for Affine MCMC Scans"""

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
                self.Steps=self.inputs['Setup'].get('Steps',100)
                if 'Walkers' in self.inputs['Setup']:
                        self.n_walkers = int(self.inputs['Setup']['Walkers'])
                else:
                        #self.n_walkers = 10 * self.runsettings.cores
                        self.n_walkers = 20 * self.n_variables
                #self.maxloss =  np.log(1 + np.finfo(np.float64).max) + 1
                self.maxloss = math.log(1+sys.float_info.max) + 1 
                # Create the likelihood functions based on 'scaling' property of variables, if they have it; otherwise use gaussian likelihood    
                #elf.get_likelihood,self.observable_masks=MakeGlobalLikelihood(self.inputs['Observables'],self.use_loglike)
                #self.get_likelihood,self.observable_masks=MakeGlobalLikelihood(self.inputs['Observables'],return_type="NLL") # need a negative log likelihood

                if self.naive:
                       self.getNLL,self.masks = MakeGlobalLikelihood(self.inputs['Observables'],return_type="NLL", smooth_losses= False) 
                else:
                       self.getNLL,self.masks = MakeGlobalLikelihood(self.inputs['Observables'],return_type="NLL", smooth_losses= True, max_loss=self.maxloss)


                self.n_observables = len(self.masks) # count *all* observables
                self.n_active_observables = sum(1 if mask else 0 for mask in self.masks)
                if 'Initial Sample' in self.inputs['Setup']:
                        from bsmart.ml.seed_points import get_previous_sample

                        seed_inputs,_ = get_previous_sample(sample_file=self.inputs['Setup']['Initial Sample'],
                                                inputs=self.inputs,
                                                nsamples = self.n_walkers,
                                                returntype =   "NLL")
                        # now we have to scale the inputs:
                        #self.seed_inputs = np.array([[f(x) for x,f in zip(point,self.downscalers)] for point in seed_inputs])

                        # NB want to set the mean to be the best point, which is the first one in the list

                        
                        self.seed_inputs = np.array([[f(x) for x,f in zip(point,self.downscalers)] for point in seed_inputs])
                else:
                        self.seed_inputs = None
                


        def postprocess(self,Point, observables, data_point,temp_dir,log, lock=None):
                """ return the likelihood; we won't get this far if the point failed to be generated """
                """ calculate likelihood """
                # For this algorithm we are minimising the objective function
                #res=1e100
                if self.naive:
                        return self.getNLL(observables)
                else:
                        return 1.0
                 
        def run_batch(self,thetas):
            # Scale thetas, run the batch, return the likelihoods
            
            batch = [ [f(y) for y,f in zip(list(theta),self.upscalers)] for theta in thetas]
            results = self.RunManager.run_batch(batch)
            
            if self.naive: # postprocess returns the NLL, or [] if invalid
                losses = np.array([value if isinstance(value, (int, float)) else self.maxloss for value in results])
                batch_observables = None
            else: # keeping 'invalid' points
                batch_observables = np.array([point["observables"] for point in self.RunManager.all_batch_points])
                losses = np.array([self.getNLL(obs) for obs in batch_observables])
            #print(f"losses: {losses}, of size {len(losses)}")
            return losses, batch_observables # want to keep the values of observables, only possible if not naive

                


        def run(self):
                """
                Affine MCMC: Uses Goodman & Weare Affine Invariant Ensemble Sampler
                
                """
                self.log.info(f"Starting Affine MCMC scan with {self.n_walkers} walkers for {self.Steps} steps.")

                n_dim = self.n_variables
                n_walkers = self.n_walkers

                # Check if n_walkers is even, required for the split update
                if n_walkers % 2 != 0:
                        self.log.warning(f"Number of walkers {n_walkers} is not even! Increasing by 1.")
                        n_walkers += 1
                        self.n_walkers += 1

                # 1. Initialization
                current_positions = np.zeros((n_walkers, n_dim))

                if self.seed_inputs is not None:
                        # Use seed inputs if available
                        n_seeds = len(self.seed_inputs)
                        if n_seeds >= n_walkers:
                                current_positions = self.seed_inputs[:n_walkers]
                        else:
                                current_positions[:n_seeds] = self.seed_inputs
                                # Fill the rest randomly
                                current_positions[n_seeds:] = np.random.rand(n_walkers - n_seeds, n_dim)
                else:
                        # Random initialization in [0, 1] hypercube
                        current_positions = np.random.rand(n_walkers, n_dim)

                # Evaluate initial NLLs
                #batch_phys = [ [f(x) for x, f in zip(pos, self.upscalers)] for pos in current_positions ]
                #current_nlls = np.array(self.RunManager.run_batch(batch_phys))
                current_nlls, current_observables = self.run_batch(current_positions)
                #print(f"current_nlls: {current_nlls}, current_observables: {current_observables}")
                # Handle potential invalid results in initial batch (though unlikely if random in [0,1])
                # If run_batch returns something other than float (e.g. None or object), handle it.
                # Assuming run_batch returns list of floats/numeric based on CMAES.py analysis.
                
                # Setup specific for Affine Invariant stretch move
                a = 2.0 # Scale parameter, usually 2.0
                
                # Arrays to store chain
                # shape: (Steps, Walkers, D)
                chain = np.zeros((self.Steps, n_walkers, n_dim))
                chain_nll = np.zeros((self.Steps, n_walkers))
                chain_observables = np.zeros((self.Steps, n_walkers, self.n_observables))

                accepted = 0
                total_proposals = 0

                varlist = list(self.inputs['Variables'].keys())
                if self.naive:
                    obslist=[]
                else:
                    obslist = list(self.inputs['Observables'].keys())
                header = ','.join(varlist+obslist+['NLL'])
                # Save Results
                run_name = self.inputs['Setup'].get('RunName', 'AffineMC_Run')
                results_dir = os.path.join(self.inputs['Setup']['cwd'], run_name, 'Results')
                if not os.path.exists(results_dir):
                        os.makedirs(results_dir)
                output_file = os.path.join(results_dir, "Results.csv")
                with open(output_file, 'w') as f:
                    f.write(header + '\n')
                for step in range(self.Steps):
                        
                        # clear stored points from memory
                        if not self.naive:
                            self.RunManager.all_points = []
                            self.RunManager.valid_points = []
                            
                        # Split walkers into two subsets: S1 (first half) and S2 (second half)
                        half_k = n_walkers // 2
                        
                        # We update S1 using S2, then S2 using S1
                        # This allows parallel updates within each subset
                        
                        #splits = [ (np.arange(0, half_k), np.arange(half_k, n_walkers)), np.arange(half_k, n_walkers), np.arange(0, half_k)) ]

                        # random split
                        all_indices = np.arange(n_walkers)
                        np.random.shuffle(all_indices)
                        splits = [(all_indices[:half_k], all_indices[half_k:]), 
                                  (all_indices[half_k:], all_indices[:half_k])]

                        for active_indices, passive_indices in splits:
                                n_active = len(active_indices)
                                n_passive = len(passive_indices)
                                
                                # Current active walkers
                                X_active = current_positions[active_indices]
                                NLL_active = current_nlls[active_indices]
                                
                                # Passive walkers to sample from (complementary set)
                                X_passive = current_positions[passive_indices]

                                # Generate proposals for all active walkers
                                # For each active walker, choose a random passive walker
                                # We can just sample n_active indices from 0 to n_passive-1
                                random_passive_indices = np.random.randint(0, n_passive, size=n_active)
                                X_j = X_passive[random_passive_indices]

                                # Stretch factor Z
                                # g(z) ~ 1/sqrt(z) in [1/a, a]
                                # CDF: G(z) = (sqrt(z) - sqrt(1/a)) / (sqrt(a) - sqrt(1/a))
                                # Inverse CDF: z = ( (sqrt(a) - sqrt(1/a))*u + sqrt(1/a) )^2
                                u = np.random.rand(n_active)
                                z = ((a - 1.0) * u + 1.0)**2 / a 
                                # Note: The standard implementation (emcee) uses ((a-1)*u + 1)^2 / a  Check:
                                # if u=0 -> 1/a. if u=1 -> a^2/a = a. Correct.

                                # Proposal: Y = X_j + z * (X_active - X_j)
                                Y = X_j + z[:, None] * (X_active - X_j)

                                # Check bounds [0, 1]
                                # If OOB, NLL = infinity
                                # We'll identify valid proposals
                                in_bounds = np.all((Y >= 0.0) & (Y <= 1.0), axis=1)
                                
                                # Prepare batch evaluation
                                # Only evaluate points in bounds
                                valid_indices_local = np.where(in_bounds)[0]
                                
                                NLL_new_active = np.full(n_active, np.inf)

                                observables_valid = None
                                if len(valid_indices_local) > 0:
                                        # Extract valid proposals
                                        Y_valid = Y[valid_indices_local]
                                        
                                        # Upscale valid proposals
                                        #batch_phys_valid = [ [f(x) for x, f in zip(pos, self.upscalers)] for pos in Y_valid ]
                                        
                                        # Run batch
                                        #results_valid = self.RunManager.run_batch(batch_phys_valid)
                                        results_valid, observables_valid = self.run_batch(Y_valid)
                                        
                                        # Store results
                                        # Make sure to handle potential non-float returns if any
                                        
                                        
                                        for idx, res in enumerate(results_valid):
                                            if isinstance(res, (int, float)):
                                                NLL_new_active[valid_indices_local[idx]] = res
                                            else:
                                                NLL_new_active[valid_indices_local[idx]] = self.maxloss
                                       

                                # 4. Metropolis-Hastings Acceptance
                                # log_ratio = (D - 1) * log(z) + NLL_old - NLL_new
                                # Accept if log(U) < log_ratio
                                
                                log_u = np.log(np.random.rand(n_active))
                                log_ratio = (n_dim - 1) * np.log(z) + NLL_active - NLL_new_active
                                
                                accept = log_u < log_ratio
                                
                                # Update positions and NLLs
                                active_indices_accepted = active_indices[accept]
                                local_indices_accepted = np.where(accept)[0]
                                
                                if len(active_indices_accepted) > 0:
                                        current_positions[active_indices_accepted] = Y[local_indices_accepted]
                                        current_nlls[active_indices_accepted] = NLL_new_active[local_indices_accepted]
                                        accepted += len(active_indices_accepted)
                                        if observables_valid is not None: # only true if not naive ...
                                            #print(observables_valid.shape)
                                            #print(current_observables.shape)
                                            
                                            # Map local accepted indices to indices in observables_valid
                                            # We know valid_indices_local contains all potentially accepted indices
                                            valid_positions = np.searchsorted(valid_indices_local, local_indices_accepted)
                                            current_observables[active_indices_accepted] = observables_valid[valid_positions]
                                
                                total_proposals += n_active

                        # Store in chain
                        chain[step] = current_positions
                        chain_nll[step] = current_nlls

                        flat_step = chain[step].reshape(-1, n_dim)
                        flat_step_phys = np.array([ [f(x) for x, f in zip(pos, self.upscalers)] for pos in flat_step ])
                        if not self.naive:
                            chain_observables[step] = current_observables
                            with open(output_file, 'a') as f:
                                np.savetxt(f, np.hstack((flat_step_phys, chain_observables[step].reshape(-1, self.n_observables), chain_nll[step].reshape(-1, 1))), delimiter=',', comments='')
                            #np.savetxt(output_file, np.hstack((flat_step_phys, chain_observables[step].reshape(-1, self.n_observables), chain_nll[step].reshape(-1, 1))), delimiter=',', comments='', header=header)
                        else:
                            with open(output_file, 'a') as f:
                                np.savetxt(f, np.hstack((flat_step_phys, chain_nll[step].reshape(-1, 1))), delimiter=',', comments='')
                        if (step + 1) % 10 == 0 or step == 0:
                             self.log.info(f"Step {step+1}/{self.Steps} completed. Mean NLL: {np.mean(current_nlls):.4f}. Acceptance rate: {accepted/total_proposals:.4f}")
                             

                self.log.info(f"Finished AffineMC run. Final Acceptance rate: {accepted/total_proposals:.4f}")
                
                
                
                # Flatten chain for CSV export: (Steps * Walkers, D)
                #flat_chain = chain.reshape(-1, n_dim)
                
                # Upscale all points for saving
                #flat_phys = np.array([ [f(x) for x, f in zip(pos, self.upscalers)] for pos in flat_chain ])
                #np.savetxt(output_file, flat_phys, delimiter=',', comments='', header=header)

                #np.savetxt(output_file, np.hstack((flat_phys, chain_observables.reshape(-1, self.n_observables), chain_nll.reshape(-1, 1))), delimiter=',', comments='', header=header)
                self.log.info(f"Results saved to {output_file}")

                if not self.runsettings.debug:
                    print(f"Finished AffineMC run. Final Acceptance rate: {accepted/total_proposals:.4f}")
                    print(f"Results saved to {output_file}")
                from bsmart import BSMplots

                ## Make an extra corner plot 

                plotstub=os.path.join(results_dir,"corner_plot")
                BSMplots.make_auto_corner_csv(plotstub,os.path.join(results_dir,"Results.csv"),varlist,{})
                os.system('python3 '+plotstub+'.py')

                plotstub=os.path.join(results_dir,'fancy_auto_corner')
                BSMplots.make_fancy_auto_corner_csv(plotstub,os.path.join(results_dir,"Results.csv"),varlist,{})
            
                os.system('python3 '+plotstub+'.py')