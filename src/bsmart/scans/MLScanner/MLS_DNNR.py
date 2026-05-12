"""!
@brief MLScanner DNNR (Deep Neural Network Regressor) method
@ingroup scans

This scan implements an active learning strategy using a Deep Neural Network Regressor (DNNR)
to efficiently find "good" points in a parameter space. A point is considered "good" if its
primary observable (NLL) is below a specified threshold.

The process is as follows:
1.  **Initialization**: The scan begins by evaluating a small set of randomly generated points (`Bootstrap_Points`).
    It can also load an initial dataset from a CSV file (`InitCSV`).

2.  **Initial Training**: A Deep Neural Network is trained on this initial dataset to predict
    the Negative Log Likelihood (NLL) from the input parameters. The network uses PyTorch.

3.  **Active Learning Loop**: The scan enters a loop to iteratively discover new good points until
    a `Target_Points` count is reached. In each iteration:
    a. A large number of `Candidate_Points` are randomly generated.
    b. The trained DNN predicts the NLL for these candidates.
    c. candidates with the lowest predicted NLL (best quality), plus a small `Random_Fraction`,
       are selected for evaluation by the physics code.
    d. **Retraining**: The DNN is retrained with the newly discovered points, becoming
       progressively better at identifying promising regions (low NLL).
       Training uses MSE Loss, AdamW optimizer, and early stopping.

4.  **Data Collection**: All discovered good points (NLL < Threshold) are returned.

This method is particularly effective for high-dimensional parameter spaces where exhaustive
scanning is computationally prohibitive.
!"""

__meta__ = {
 "name": "MLS_DNNR",
 "requires": ["torch", "pandas", "numpy", "sklearn"],
 "settings": {
    "Networks": {
        "Iterations": "Number of active learning iterations (default: 10).",
        "Candidate_Points": "Number of candidate points to generate and score in each iteration (default: 500).",
        "Bootstrap_Points": "Number of initial random points to evaluate (default: 100).",
        "Points_Per_Iteration": "Number of candidate points to evaluate in each iteration (default: 300).",
        "Threshold_Value": "The threshold for the NLL to consider a point 'good' (default: 1).",
        "Random_Fraction": "Fraction of points per iteration to be selected randomly (default: 0.2).",
        "Batch Size": "Batch size for training the neural network (default: 500).",
        "Neurons": "Number of neurons in the hidden layers of the DNN (default: 100).",
        "Epochs": "Number of epochs for training the DNN in each cycle (default: 250).",
        "LearningRate": "Learning rate for the AdamW optimizer (default: 1e-2).",
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
import torch
import torch.nn as nn
import sklearn
import random as random
import math
from bsmart.BSMlikelihood import MakeLikelihoods, MakeGlobalLikelihood, safe_float
from torch.utils.data import TensorDataset, DataLoader
from bsmart import debug

def generate_param_points(inputs, num_points):
    variables_range = []
    for varname in inputs['Variables']:
        if 'RANGE' in inputs['Variables'][varname]:
            varmin = inputs['Variables'][varname]['RANGE'][0]
            varmax = inputs['Variables'][varname]['RANGE'][1]
            variables_range.append(np.random.uniform(varmin, varmax, num_points))
    return np.array(variables_range).T

def torchfit(x, y, model, criterion, optimizer, epochs, batch_size, verbose=True, patience=10, min_delta=1e-5):
    x_tensor = torch.FloatTensor(x)
    y_tensor = torch.FloatTensor(y).reshape(-1, 1)
    train_dataset = TensorDataset(x_tensor, y_tensor)
    train_loader = DataLoader(train_dataset, batch_size=int(batch_size), shuffle=True)

    best_loss = float('inf')
    patience_counter = 0
    best_model_state = None


    for epoch in range(epochs):
        model.train()
        epoch_loss = 0.0
        num_batches = 0
        for batch_x, batch_y in train_loader:
            outputs = model(batch_x)
            loss = criterion(outputs, batch_y)
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            epoch_loss += loss.item()
            num_batches += 1
        if verbose and (epoch % 50 == 0 or epoch == epochs - 1):
            avg_loss = epoch_loss / num_batches
            print(f'  Epoch {epoch}/{epochs}, Loss: {avg_loss:.6f}')

        # Early stopping check
        if epoch_loss < best_loss - min_delta:
            best_loss = epoch_loss
            patience_counter = 0
            best_model_state = model.state_dict()
        else:
            patience_counter += 1
            if patience_counter >= patience:
                if verbose:
                    print(f'Early stopping triggered at epoch {epoch}. Best loss: {best_loss:.6f}')
                break
    
    if best_model_state:
        model.load_state_dict(best_model_state)

    return best_loss

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



        self.candidate_points = int(self.inputs['Networks'].get('Candidate_Points', 5))
        self.bootstrap_points = int(self.inputs['Networks'].get('Bootstrap_Points', 10))
        self.points_per_iter = int(self.inputs['Networks'].get('Points_Per_Iteration', 5))

        self.threshold_value = float(self.inputs['Networks'].get('Threshold_Value', 1))

        self.function_dim = int(len(self.inputs['Variables']))
        self.random_fraction = float(self.inputs['Networks'].get('Random_Fraction', 0.2))
        self.batch_size = int(self.inputs['Networks'].get('Batch Size', 5))
        self.neurons = int(self.inputs['Networks'].get('Neurons', 100))
        self.epoch = int(self.inputs['Networks'].get('Epochs', 1000))
        self.learning_rate = float(self.inputs['Networks'].get('LearningRate', 1e-2))
        self.dropout = 0.2

        self.verbose = int(self.inputs['Networks'].get('Verbose', 0))
        print(f"verbos: {self.verbose}")
        if "Cores" in self.inputs['Setup']:
            self.ncores = int(self.inputs['Setup']['Cores'])
            torch.set_num_threads(min(8, self.ncores))

        self.target_points = int(self.inputs['Setup'].get('Points', 10000))
        self.n_variables = len(self.inputs['Variables'])

        self.maxloss = np.log(1 + np.finfo(np.float64).max) + 1 # copied from Mark's code
        self.likelihood_fns, self.observable_masks = MakeLikelihoods(self.inputs["Observables"], loglike=True)  # there is no primary observable of choice. the primary observable is the likelihood. 
        self.primary_observable = 'NLL'

        self.init_csv_file = self.inputs['Setup'].get('InitCSV')
        if self.init_csv_file:
            if os.path.isfile(self.init_csv_file):
                self.log.info(f'Will use InitCSV file: {self.init_csv_file}')
            else:
                self.log.error(f'InitCSV file not found: {self.init_csv_file}')
                self.init_csv_file = None

    def network_architecture_selector(self, n_training_points):
        
        samples_per_dim = n_training_points / self.function_dim

        if n_training_points >= 5000:
            self.log.info(f'Network: SiLU + Dropout(0.2) + Funnel architecture')
            self.log.info(f'  Reason: Large dataset ({n_training_points} points, {samples_per_dim:.1f} samples/dim)')
            self.log.info(f'  Architecture: {self.function_dim}→{self.neurons*2}→{self.neurons}→{self.neurons//2}→1')
            
            dnn_model = nn.Sequential(
                nn.Linear(self.function_dim, self.neurons * 2),
                nn.SiLU(),
                nn.Dropout(0.2),
                nn.Linear(self.neurons * 2, self.neurons),
                nn.SiLU(),
                nn.Dropout(0.2),
                nn.Linear(self.neurons, self.neurons // 2),
                nn.SiLU(),
                nn.Dropout(0.2),
                nn.Linear(self.neurons // 2, 1)
            )
        
        else:
            self.log.info(f'Network: ReLU + Flat architecture')
            self.log.info(f'  Reason: Limited data ({n_training_points} points, {samples_per_dim:.1f} samples/dim)')
            self.log.info(f'  Architecture: {self.function_dim}→{self.neurons}→{self.neurons}→{self.neurons}→{self.neurons}→{self.neurons}→1')
            
            dnn_model = nn.Sequential(
                nn.Linear(self.function_dim, self.neurons),
                nn.ReLU(),
                nn.Linear(self.neurons, self.neurons),
                nn.ReLU(),
                nn.Linear(self.neurons, self.neurons),
                nn.ReLU(),
                nn.Linear(self.neurons, self.neurons),
                nn.ReLU(),
                nn.Linear(self.neurons, self.neurons),
                nn.ReLU(),
                nn.Linear(self.neurons, 1)
            )
        
        return dnn_model
    

    def postprocess(self, Point, observables, data_point, temp_dir, log, lock=None):  # type: ignore
        # copied from Mark's code
        if self.naive:
            return np.sum(self.get_losses(observables))
        else:
            return 1.0
        
    def get_losses(self,observables):
            # copied from Mark's code

        """ Returns a list of losses.
                    
        The C-function loss should be zero if the observable is within allowed bounds, and greater than zero outside it.
        If we include a validity condition, then other observables, evaluated afterwards, will be returned as NaN. 
         -> We should assign these the maximum loss ~ 710.                 
        This is then compatible with the explicit use of 'hierarchical observables', where all other observables are set to the maximum.
        """

        likeit=iter(self.likelihood_fns)

        return [self.smooth_cap_loss((next(likeit))(val)) if mask and not math.isnan(val := safe_float(v)) else float((next(likeit) and False) or self.maxloss)  for v, mask in zip(observables, self.observable_masks) if mask]
    
    def smooth_cap_loss(self, x):
        # copied from Mark's code
        return -self.maxloss*np.expm1(x/self.maxloss) 
    
    def extract_from_valid_points(self,valid_points):
        params = np.array([point['inputs'] for point in valid_points])
        batch_observables = [point["observables"] for point in valid_points]

        losses = np.array([self.get_losses(obs) for obs in batch_observables])
        fitnesses = np.sum(losses, axis=1)

        return params, fitnesses

    def run(self):
        # 1. Bootstrap
        self.RunManager.run_batch(generate_param_points(self.inputs, self.bootstrap_points).tolist())
        if not self.RunManager.valid_batch_points: sys.exit("No valid points found in initial batch!")
        
        all_params, all_nll = self.extract_from_valid_points(self.RunManager.valid_batch_points)

        # 2. Init CSV
        # 2. Init CSV
        if self.init_csv_file:
            self.log.info(f'Loading InitCSV: {self.init_csv_file}')
            df = pd.read_csv(self.init_csv_file).dropna()
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

        # 4. Setup Model
        dnn_model = self.network_architecture_selector(len(all_params))
        optimizer = torch.optim.AdamW(dnn_model.parameters(), lr=self.learning_rate, weight_decay=1e-4)
        criterion = nn.MSELoss()

        # 5. Initial Training
        self.log.info(f'Initial training on {len(all_params)} points...')
        torchfit(all_params, all_nll, dnn_model, criterion, optimizer, self.epoch, self.batch_size, verbose=bool(self.verbose))
        
        # Validation
        dnn_model.eval()
        with torch.no_grad():
            idx = np.random.choice(len(all_params), min(5, len(all_params)), replace=False)
            preds = dnn_model(torch.FloatTensor(all_params[idx])).squeeze().numpy()
            for p, a in zip(preds, all_nll[idx]):
                self.log.info(f'Pred: {p:.2f} vs Act: {a:.2f} (Err: {abs(p-a):.2f})')

        # 6. Active Learning Loop
        n_good = np.sum(all_nll < self.threshold_value)
        initial_good_points = n_good
        self.log.info(f"Initial good points: {initial_good_points}. Target new points: {self.target_points}. Stopping at: {initial_good_points + self.target_points}")
        
        run_num = 0

        while n_good < (self.target_points + initial_good_points):
            run_num += 1
            
            # Candidates
            cands = generate_param_points(self.inputs, self.candidate_points)
            
            dnn_model.eval()
            with torch.no_grad():
                preds = dnn_model(torch.FloatTensor(cands)).squeeze().numpy()
            
            ml_cands = cands[preds < self.threshold_value]
            
            # Batch Construction
            n_ml = int(self.points_per_iter * (1 - self.random_fraction))
            batch = ml_cands[:n_ml]
           
            #n_rnd = int(self.points_per_iter * self.random_fraction)
            n_rnd = self.points_per_iter - len(batch)
            self.log.info(f"ML predicts {len(batch)} good points; adding {n_rnd} random points")
            rnd_pts = generate_param_points(self.inputs, n_rnd)
            
            batch = np.vstack([batch, rnd_pts]) if len(batch) > 0 else rnd_pts
            
            # Execution
            self.RunManager.run_batch(batch.tolist())
            if not self.RunManager.valid_batch_points: continue
            
            # Update Data
            new_params, new_nll = self.extract_from_valid_points(self.RunManager.valid_batch_points)
            all_params = np.vstack([all_params, new_params])
            all_nll = np.concatenate([all_nll, new_nll])
            
            # Retrain
            mask = all_nll < 1e9
            torchfit(all_params[mask], all_nll[mask], dnn_model, criterion, optimizer, self.epoch, self.batch_size, verbose=True)
            
            n_good = np.sum(all_nll < self.threshold_value)
            new_good = np.sum(new_nll < self.threshold_value)
            self.log.info(f'Run {run_num}: +{new_good} good. Total: {n_good}')

        self.log.info(f'Target reached: {n_good} points.')
        mask = all_nll < self.threshold_value
        return all_params[mask], all_nll[mask]
