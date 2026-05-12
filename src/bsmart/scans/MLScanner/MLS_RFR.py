"""
MLScanner RFR method
--------------------


MLScanner method `MLS_RFR` based on code from 

* `<https://github.com/AHamamd150/MLscanner>`_
* `<https://arxiv.org/abs/2207.09959>`_

This scan implements an active learning strategy using a Random Forest Regressor (RFR)
to efficiently find "good" points in a parameter space. A point is considered "good" if its
primary observable (NLL) is below a specified threshold.

The process is as follows:

1.  **Initialization**: The scan begins by evaluating a small set of randomly generated points (`Bootstrap_Points`).
    It can also load an initial dataset from a CSV file (`InitCSV`).

2.  **Initial Training**: A Random Forest Regressor is trained on this initial dataset to predict
    the Negative Log Likelihood (NLL) from the input parameters.

3.  **Active Learning Loop**: The scan enters a loop to iteratively discover new good points until
    a `Target_Points` count is reached. In each iteration:

    a. A large number of `Candidate_Points` are randomly generated.
    b. The trained RFR model predicts the NLL for these candidates.
    c. candidates with the lowest predicted NLL (best quality), plus a small `Random_Fraction`,
       are selected for evaluation by the physics code.
    d. **Retraining**: The RFR is retrained with the newly discovered points, becoming
       progressively better at identifying promising regions (low NLL).

4.  **Data Collection**: All discovered good points (NLL < Threshold) are returned.

This method is particularly effective for high-dimensional parameter spaces where exhaustive
scanning is computationally prohibitive.
"""

__meta__ = {
 "name": "MLS_RFR",
 "requires": ["sklearn", "pandas", "numpy"],
 "settings": {
    "Networks": {
        "Iterations": "Number of active learning iterations (default: 10).",
        "Candidate_Points": "Number of candidate points to generate and score in each iteration (default: 500).",
        "Bootstrap_Points": "Number of initial random points to evaluate (default: 100).",
        "Points_Per_Iteration": "Number of candidate points to evaluate in each iteration (default: 300).",
        "Threshold_Value": "The threshold for the NLL to classify a point as 'good' (default: 1).",
        "Random_Fraction": "Fraction of points per iteration to be selected randomly, for exploration (default: 0.2).",
        "Estimators": "Number of trees in the forest (default: 100).",
        "Max_Depth": "Maximum depth of the tree (default: 30).",
        "Min_Samples_Split": "The minimum number of samples required to split an internal node (default: 2).",
        "Min_Samples_Leaf": "The minimum number of samples required to be at a leaf node (default: 1).",
        "Verbose": "Verbosity level (default: 0)."
    },
    "Setup": {
        "InitCSV": "Path to an optional CSV file with initial points to seed the scan.",
        "Points": "Number of points to generate *in total* before stopping (default: 1000)"
    }
    }
}

from bsmart.core import Scan as Scan
import os
import sys

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
import math
from bsmart.BSMlikelihood import MakeLikelihoods, MakeGlobalLikelihood, safe_float
from bsmart import debug

def generate_param_points(inputs, num_points):
    variables_range = []
    for varname in inputs['Variables']:
        if 'RANGE' in inputs['Variables'][varname]:
            varmin = inputs['Variables'][varname]['RANGE'][0]
            varmax = inputs['Variables'][varname]['RANGE'][1]
            variables_range.append(np.random.uniform(varmin, varmax, num_points))
    return np.array(variables_range).T

class NewScan(Scan):

    def initialise(self):
        """ Need to make sure we override certain settings """
        self.runsettings.store_points_in_memory = True
        self.runsettings.invalid_return_value = [0]

        if self.runsettings.store_invalid_points:
            self.naive = False 
            self.runsettings.invalid_return_value = [0]
        else:
            self.naive = True  # we treat invalid points as bad
            self.runsettings.invalid_return_value = []
        self.citations = """@article{Hammad:2022wpq,
    author = "Hammad, A. and Park, Myeonghun and Ramos, Raymundo and Saha, Pankaj",
    title = "{Exploration of parameter spaces assisted by machine learning}",
    eprint = "2207.09959",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    doi = "10.1016/j.cpc.2023.108902",
    journal = "Comput. Phys. Commun.",
    volume = "293",
    pages = "108902",
    year = "2023"
}        
"""

    def __init__(self, inputs, log):
        Scan.__init__(self, inputs, log)
        print("start scan init")

        # Helper to get setting from Networks or Setup
        def get_setting(key, default, cast_type=int):
            val = self.inputs['Networks'].get(key)
            if val is None:
                val = self.inputs['Setup'].get(key)
            if val is None:
                return default
            return cast_type(val)

        self.iteration = get_setting('Iterations', 10)
        self.candidate_points = get_setting('Candidate_Points', 500)
        self.bootstrap_points = get_setting('Bootstrap_Points', 100)
        self.points_per_iter = get_setting('Points_Per_Iteration', 300)
        self.threshold_value = get_setting('Threshold_Value', 1, float)
        #self.target_points = get_setting('Target_Points', 20000)
        self.function_dim = int(len(self.inputs['Variables']))
        self.random_fraction = get_setting('Random_Fraction', 0.2, float)
        
        # RFR specific settings
        self.n_estimators = get_setting('Estimators', 100)
        self.max_depth = get_setting('Max_Depth', 30)
        self.min_samples_split = get_setting('Min_Samples_Split', 2)
        self.min_samples_leaf = get_setting('Min_Samples_Leaf', 1)
        
        self.verbose = get_setting('Verbose', 0)
        print(f"verbos: {self.verbose}")

        if "Cores" in self.inputs['Setup']:
            self.ncores = int(self.inputs['Setup']['Cores'])
            self.log.info('Setting number of cores to %d' % self.ncores)

        self.target_points = int(self.inputs['Setup'].get('Points', 10000))
        self.n_variables = len(self.inputs['Variables'])

        print("make likelihoods")
        self.maxloss =  np.log(1 + np.finfo(np.float64).max) + 1
        self.likelihood_fns, self.observable_masks = MakeLikelihoods(self.inputs["Observables"], loglike=True)  

        self.primary_observable = 'NLL'

        self.InitCSV = self.inputs['Setup'].get('InitCSV')
        if self.InitCSV:
            if os.path.isfile(self.InitCSV):
                self.log.info(f'Will use InitCSV file: {self.InitCSV}')
            else:
                self.log.error(f'InitCSV file not found: {self.InitCSV}')
                self.InitCSV = None
        
        print("finished scan init")

    def postprocess(self, Point, observables, data_point, temp_dir, log, lock=None):
        """ return the likelihood; we won't get this far if the point failed to be generated """
        if self.naive:
            return np.sum(self.get_losses(observables))
        else:
            return 1.0

    def smooth_cap_loss(self, x):
        """
        Caps the loss by applying a sigmoid.
        This is useful for losses that are unbounded.                   
        """
        return -self.maxloss*np.expm1(x/self.maxloss) # assume x is negative and want to change its sign
    
    def get_losses(self, observables):
        """ Returns a list of losses. """
        likeit = iter(self.likelihood_fns)
        return [self.smooth_cap_loss((next(likeit))(val)) if mask and not math.isnan(val := safe_float(v)) else float((next(likeit) and False) or self.maxloss) 
                for v, mask in zip(observables, self.observable_masks) if mask]

    def extract_from_valid_points(self, valid_points):
        # inputs contains the parameter values
        params = np.array([point['inputs'] for point in valid_points])
        batch_observables = [point["observables"] for point in valid_points]

        # Calculate NLL (fitness) for each point using get_losses logic
        losses = np.array([self.get_losses(obs) for obs in batch_observables])
        fitnesses = np.sum(losses, axis=1)

        return params, fitnesses

    def run(self):
        RFR = RandomForestRegressor(
            n_estimators=self.n_estimators, max_depth=self.max_depth, min_samples_split=self.min_samples_split,
            min_samples_leaf=self.min_samples_leaf, n_jobs=self.ncores, random_state=42
        )

        # 1. Bootstrap
        self.RunManager.run_batch(generate_param_points(self.inputs, self.bootstrap_points).tolist())
        if not self.RunManager.valid_batch_points: sys.exit("No valid points found in initial batch!")
        
        all_params, all_nll = self.extract_from_valid_points(self.RunManager.valid_batch_points)

        # 2. Init CSV
        # 2. Init CSV
        if self.InitCSV:
            self.log.info(f'Loading InitCSV: {self.InitCSV}')
            df = pd.read_csv(self.InitCSV).dropna()
            csv_params = df[list(self.inputs['Variables'].keys())].values
            
            if self.primary_observable in df.columns:
                csv_nll = pd.to_numeric(df[self.primary_observable], errors='coerce').fillna(1e10).values
            else:
                self.log.info(f" '{self.primary_observable}' not found in InitCSV. Recalculating from Observables...")
                obs_keys = list(self.inputs['Observables'].keys())
                missing_obs = [k for k in obs_keys if k not in df.columns]
                if missing_obs: self.log.warning(f"InitCSV missing observables: {missing_obs}")

                csv_obs_data = []
                for k in obs_keys:
                    if k in df.columns:
                        csv_obs_data.append(pd.to_numeric(df[k], errors='coerce').values)
                    else:
                        csv_obs_data.append(np.full(len(df), np.nan))
                
                csv_obs_data = np.column_stack(csv_obs_data)
                csv_nll = np.array([np.sum(self.get_losses(row)) for row in csv_obs_data])

            all_params = np.vstack([all_params, csv_params])
            all_nll = np.concatenate([all_nll, csv_nll])

        # 3. Filter Failures
        mask = all_nll < 1e9
        if not np.any(mask): sys.exit('No valid training data after filtering failures!')
        all_params, all_nll = all_params[mask], all_nll[mask]

        # 4. Initial Training
        n_good = np.sum(all_nll < self.threshold_value)
        initial_good_points = n_good
        self.log.info(f"Initial good points: {initial_good_points}. Target new points: {self.target_points}. Stopping at: {initial_good_points + self.target_points}")

        self.log.info(f'Training: {len(all_nll)} pts, {n_good} good ({n_good/len(all_nll)*100:.1f}%). Stats: Mean={np.mean(all_nll):.2f}, Std={np.std(all_nll):.2f}')
        RFR.fit(all_params, all_nll)

        # Validation
        idx = np.random.choice(len(all_params), min(5, len(all_params)), replace=False)
        for p, a in zip(RFR.predict(all_params[idx]), all_nll[idx]):
            self.log.info(f'Pred: {p:.2f} vs Act: {a:.2f} (Err: {abs(p-a):.2f})')

        # 5. Active Learning Loop
        run_num = 0
        while n_good < (self.target_points + initial_good_points):
            run_num += 1
            
            # Candidates
            cands = generate_param_points(self.inputs, self.candidate_points)
            ml_cands = cands[RFR.predict(cands) < self.threshold_value]
            
            # Batch Construction
            n_ml = int(self.points_per_iter * (1 - self.random_fraction))
            batch = ml_cands[:n_ml]
            n_rnd = self.points_per_iter - len(batch)
            self.log.info(f"ML predicts {len(batch)} good points; adding {n_rnd} random points")
            #n_rnd = int(self.points_per_iter * self.random_fraction)
            
            rnd_pts = generate_param_points(self.inputs, n_rnd)
            
            batch = np.vstack([batch, rnd_pts]) if len(batch) > 0 else rnd_pts
            
            # Execution
            self.RunManager.run_batch(batch.tolist())
            if not self.RunManager.valid_batch_points: continue
            
            # Update Data
            new_params, new_nll = self.extract_from_valid_points(self.RunManager.valid_batch_points)
            all_params = np.vstack([all_params, new_params])
            all_nll = np.concatenate([all_nll, new_nll])
            
            # Retrain (ignoring failures)
            mask = all_nll < 1e9
            RFR.fit(all_params[mask], all_nll[mask])
            
            n_good = np.sum(all_nll < self.threshold_value)
            new_good = np.sum(new_nll < self.threshold_value)
            self.log.info(f'Run {run_num}: +{new_good} good. Total: {n_good}')

        self.log.info(f'Target reached: {n_good} points.')
        mask = all_nll < self.threshold_value
        return all_params[mask], all_nll[mask]
