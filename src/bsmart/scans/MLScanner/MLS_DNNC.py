"""
MLScanner DNNC method
---------------------


MLScanner method `MLS_DNNC` based on code from 

* `<https://github.com/AHamamd150/MLscanner>`_
* `<https://arxiv.org/abs/2207.09959>`_

This scan implements an active learning strategy using a Deep Neural Network Classifier (DNNC)
to efficiently find "good" points in a parameter space. A point is considered "good" if the likelihood is below a specified threshold. This is a generalisation of the original algorithms; in the original package the scan looked for a `primary observable` to compare to a threshold. Since the likelihood can be set as "EXPUSER" for a given observable, the original case can also be accommodated -- but a likelihood is more generally useful. 

The process is as follows:

1.  **Initialization**: The scan begins by evaluating a small set of randomly generated points
    (``Bootstrap_Points``). It can also load an initial dataset from a CSV file (``InitCSV``).

2.  **Initial Training**: A DNN classifier is trained on this initial dataset. The points are
    labeled as "good" (1) or "bad" (0) based on the ``Threshold_Value``. The loss function
    can be weighted to handle the initial class imbalance.

3.  **Active Learning Loop**: The scan enters a loop to iteratively discover new good points until
    a ``Target_Points`` count is reached. In each iteration:
    
    a. A large number of ``Candidate_Points`` are randomly generated.
    b. The trained DNN predicts which of these candidates are likely to be good.
    c. The most promising candidates, plus a small ``Random_Fraction``, are selected for
       evaluation by the physics code.
    d. **SMOTE for Imbalance**: If too few good points are found, the Synthetic Minority
       Over-sampling TEchnique (SMOTE) is used to generate synthetic good points to
       balance the dataset and improve classifier performance in sparse regions.
    e. **Retraining**: The DNN is retrained with the newly discovered points, becoming
       progressively better at identifying promising regions of the parameter space.

4.  **Data Collection**: All discovered good points and their observable values are saved to a CSV file.

This method is particularly effective for high-dimensional parameter spaces where exhaustive
scanning is computationally prohibitive.
"""


__meta__ = {
 "name": "MLS_DNNC",
 "requires": ["imblearn", "sklearn", "torch", "pandas"],
 "settings": {
    "Networks": {
        "Iterations": "Number of active learning iterations (default: 10).",
        "Candidate_Points": "Number of candidate points to generate and score in each iteration (default: 500).",
        "Bootstrap_Points": "Number of initial random points to evaluate (default: 100).",
        "Points_Per_Iteration": "Number of candidate points to evaluate in each iteration (default: 300).",
        "Kinitial": "Number of neighbors for the SMOTE algorithm (default: 5).",
        "Threshold_Value": "The threshold for the primary observable to classify a point as 'good' (default: 100).",
        "Random_Fraction": "Fraction of points per iteration to be selected randomly, for exploration (default: 0.2).",
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
from imblearn.over_sampling import SMOTE
import random as random
import math
from bsmart.BSMlikelihood import MakeLikelihoods, MakeGlobalLikelihood, safe_float



def labeler(obs_values, threshold):
    obs_values = np.array(obs_values).flatten()
    return (obs_values < threshold).astype(float)

def generate_param_points(inputs, num_points):
    variables_range = []
    for varname in inputs['Variables']:
        if 'RANGE' in inputs['Variables'][varname]:
            varmin = inputs['Variables'][varname]['RANGE'][0]
            varmax = inputs['Variables'][varname]['RANGE'][1]
            variables_range.append(np.random.uniform(varmin, varmax, num_points))
    return np.array(variables_range).T

def torchfit(X, y, model, criterion, optimizer, epochs, batch_size, verbose=True, patience=10, min_delta=1e-5):
    X_tensor = torch.FloatTensor(X)
    y_tensor = torch.FloatTensor(y).reshape(-1, 1)
    train_dataset = torch.utils.data.TensorDataset(X_tensor, y_tensor)
    train_loader = torch.utils.data.DataLoader(train_dataset, batch_size=int(batch_size), shuffle=True)

    best_loss = float('inf')
    patience_counter = 0
    best_model_state = None

    for epoch in range(epochs):
        model.train()
        running_loss = 0.0
        num_batches = 0
        for batch_x, batch_y in train_loader:
            optimizer.zero_grad()
            outputs = model(batch_x)
            loss = criterion(outputs, batch_y)
            loss.backward()
            optimizer.step()
            running_loss += loss.item()
            num_batches += 1
        
        epoch_loss = running_loss / num_batches

        if verbose and (epoch % 100 == 0 or epoch == epochs - 1):
            print(f'  Epoch {epoch}/{epochs}, Loss: {epoch_loss:.6f}')

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



"""
def extract_from_valid_points(valid_points):
    params = np.array([point['inputs'] for point in valid_points])
    obs = np.array([float(point['observables'][0]) if 'observables' in point and len(point['observables']) > 0 else 1e10
                    for point in valid_points], dtype=np.float64)
    return params, obs
"""

class NewScan(Scan):

    def initialise(self):
        """ Need to make sure we override certain settings """
        self.runsettings.store_points_in_memory = True
        self.runsettings.invalid_return_value = 0

        if self.runsettings.store_invalid_points:
            self.naive = False 
            self.runsettings.invalid_return_value = 0
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
        
        self.iteration = int(self.inputs['Networks'].get('Iterations', 10))
        
        self.candidate_points = int(self.inputs['Networks'].get('Candidate_Points', 500))
        self.bootstrap_points = int(self.inputs['Networks'].get('Bootstrap_Points', 100))
        self.points_per_iter = int(self.inputs['Networks'].get('Points_Per_Iteration', 300))
        self.smote_neighbours = int(self.inputs['Networks'].get('Kinitial', 5))
        #self.threshold_value = float(self.inputs['Networks'].get('Threshold_Value', 100))
        self.threshold_value = float(self.inputs['Networks'].get('Threshold_Value', 1))
        #self.target_points = int(self.inputs['Networks'].get('Target_Points', 20000))
        self.function_dim = int(len(self.inputs['Variables']))
        self.random_fraction = float(self.inputs['Networks'].get('Random_Fraction', 0.2))
        self.batch_size = int(self.inputs['Networks'].get('Batch Size', 500))
        self.neurons = int(self.inputs['Networks'].get('Neurons', 100))
        self.epoch = int(self.inputs['Networks'].get('Epochs', 250))
        self.learning_rate = float(self.inputs['Networks'].get('LearningRate', 1e-2))
        self.n_estimators = int(self.inputs['Networks'].get('Estimators', 100))
        self.max_depth = int(self.inputs['Networks'].get('Max_Depth', 30))
        self.verbose = int(self.inputs['Networks'].get('Verbose', 0))
        print(f"verbos: {self.verbose}")
        if "Cores" in self.inputs['Setup']:
            self.ncores = int(self.inputs['Setup']['Cores'])
            torch.set_num_threads(min(8, self.ncores))

        #self.num_points = int(self.inputs['Setup'].get('Points', 10000))

        self.target_points = int(self.inputs['Setup'].get('Points', 10000))
        self.n_variables = len(self.inputs['Variables'])

        #self.loglike = eval(self.inputs['Setup'].get('LogLike',"False"))
        # choose primary observable from likelihoods

        print("make likelihoods")
        self.maxloss =  np.log(1 + np.finfo(np.float64).max) + 1

        # Create the likelihood functions based on 'scaling' property of variables, if they have it; otherwise use gaussian likelihood    
        #elf.get_likelihood,self.observable_masks=MakeGlobalLikelihood(self.inputs['Observables'],self.use_loglike)
        #self.get_likelihood,self.observable_masks=MakeGlobalLikelihood(self.inputs['Observables'],return_type="NLL") # need a negative log likelihood
        self.likelihood_fns, self.observable_masks = MakeLikelihoods(self.inputs["Observables"], loglike=True)

        ##
        #self.likelihood, self.observable_masks = MakeGlobalLikelihood(self.inputs["Observables"], self.loglike)
        """
        self.header = []
        if 'Observables' in self.inputs:
            for observable_name in self.inputs['Observables']:
                self.header.append(observable_name)
        if len(self.header) == 0:
            sys.exit("No observables found in inputs!")
        self.primary_observable = self.header[0]
        """

        self.primary_observable = 'NLL'
        self.init_csv_file = None
        self.init_csv_bool = False
        if 'InitCSV' in self.inputs['Setup']:
            self.init_csv_file = self.inputs['Setup']['InitCSV']
            if not os.path.isfile(self.init_csv_file):
                self.log.error('InitCSV file not found: %s' % self.init_csv_file)
                self.init_csv_file = None
            else:
                self.log.info('Will use InitCSV file: %s' % self.init_csv_file)
                self.init_csv_bool = True

        print("finished scan init")

    def loss_fn_selector(self, training_labels):

        criterion = None
        n_good = int(np.sum(training_labels))
        n_total = len(training_labels)
        good_ratio = n_good / n_total if n_total > 0 else 0

        if self.init_csv_bool:
            criterion = nn.BCEWithLogitsLoss()
        elif good_ratio < 0.1:
            pos_weight = torch.tensor([(1.0 - good_ratio) / good_ratio])
            criterion = nn.BCEWithLogitsLoss(pos_weight=pos_weight)
        else:
            criterion = nn.BCEWithLogitsLoss()

        return criterion

    def network_architecture_selector(self, n_training_points):

        if n_training_points >= 5000:

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

    def extract_from_valid_points(self,valid_points):
        params = np.array([point['inputs'] for point in valid_points])
        batch_observables = [point["observables"] for point in valid_points]

        losses = np.array([self.get_losses(obs) for obs in batch_observables])
        fitnesses = np.sum(losses, axis=1)

        return params, fitnesses
    
    def run(self):

        initial_points = []
        for _ in range(self.bootstrap_points):
            point = [np.random.uniform(self.inputs['Variables'][var]['RANGE'][0], self.inputs['Variables'][var]['RANGE'][1]) for var in self.inputs['Variables']]
            initial_points.append(point)

        self.RunManager.run_batch(initial_points)

        if len(self.RunManager.valid_batch_points) == 0:
            sys.exit("No valid points found in initial batch!")

        bootstrap_params, bootstrap_obs = self.extract_from_valid_points(self.RunManager.valid_batch_points)

        print(bootstrap_params)
        print(bootstrap_obs)
        if self.init_csv_file is not None:
            self.log.info('Loading initial points from CSV: %s' % self.init_csv_file)
            var_columns = list(self.inputs['Variables'].keys())
            #obs_columns = [self.primary_observable]
            obs_columns = ['NLL']
            input_dataframe = pd.read_csv(self.init_csv_file)
            useful_dataset = pd.DataFrame(input_dataframe, columns=var_columns+obs_columns).dropna()
            self.log.info('CSV has %d points' % useful_dataset.shape[0])

            csv_vars = useful_dataset[var_columns].values
            csv_obs = pd.to_numeric(useful_dataset[self.primary_observable], errors='coerce').astype(np.float64).values

            training_params = np.append(bootstrap_params, csv_vars, axis=0)
            all_obs_values = np.append(bootstrap_obs, csv_obs).astype(np.float64)
        else:
            training_params = bootstrap_params
            all_obs_values = bootstrap_obs

        all_labels = labeler(all_obs_values, self.threshold_value)
        good_obs_values = all_obs_values[all_labels == 1]

        dnn_model = self.network_architecture_selector(len(training_params))
        criterion = self.loss_fn_selector(all_labels)
        optimizer = torch.optim.AdamW(dnn_model.parameters(), lr=self.learning_rate, weight_decay=1e-4)

        torchfit(training_params, all_labels, dnn_model, criterion, optimizer, self.epoch, self.batch_size, verbose=True)

        good_params = training_params[all_labels == 1]
        good_labels = all_labels[all_labels == 1]
        bad_params = training_params[all_labels == 0][:len(good_labels)]
        bad_labels = all_labels[all_labels == 0][:len(good_labels)]

        if len(good_params) < self.function_dim:
            sys.exit('First training cannot find any good points, please try different size..')

        var_names = list(self.inputs['Variables'].keys())
        accumulated_csv_path = os.path.join(self.inputs['Setup']['cwd'], self.inputs['Setup']['RunName'],
                                            'Spectrum_Files', 'Accumulated_points_Model_dnn_classifier.csv')

        initial_good_data = pd.DataFrame(good_params, columns=var_names)
        initial_good_data[self.primary_observable] = good_obs_values
        initial_good_data.to_csv(accumulated_csv_path, index=False)
        self.log.info(f'Initialized accumulated CSV with {len(good_params)} initial good points')

        run_number = 0

        while len(good_labels) < self.target_points:
            run_number += 1
            previous_good_count = len(good_params)
            candidate_params = generate_param_points(self.inputs, self.candidate_points)

            dnn_model.eval()
            with torch.no_grad():
                candidate_params_tensor = torch.FloatTensor(candidate_params)
                logits = dnn_model(candidate_params_tensor).squeeze()
                predicted_scores = torch.sigmoid(logits).numpy()

            sorted_indices = np.argsort(predicted_scores)[::-1]
            batch_params = candidate_params[sorted_indices][:round(self.points_per_iter * (1 - self.random_fraction))]
            batch_params = np.append(batch_params, candidate_params[:round(self.points_per_iter * self.random_fraction)], axis=0)

            self.RunManager.run_batch(batch_params.tolist())

            if len(self.RunManager.valid_batch_points) == 0:
                continue

            validated_params, new_obs_values = self.extract_from_valid_points(self.RunManager.valid_batch_points)

            new_labels = labeler(new_obs_values, self.threshold_value)

            if len(new_labels[new_labels == 1]) <= self.smote_neighbours:

                ## Stupid garbage in original MLScanner
                #smote_input_labels = np.append(np.zeros(len(new_labels[new_labels == 0])), np.ones(len(good_params[-int(self.function_dim):])), axis=0)
                #smote_input_params = np.append(candidate_params[:len(new_labels[new_labels == 0])], good_params[np.random.randint(len(good_params), size=self.function_dim), :], axis=0)

                ## want to have enough neighbours:
                nbad= len(new_labels[new_labels == 0])
                ngood = len(new_labels[new_labels == 1])
                
                smote_input_labels = np.append(np.zeros(nbad), np.ones(2*self.smote_neighbours), axis=0)
                                    
                #smote_input_params = np.append(validated_params[new_labels == 0], validated_params[new_labels == 1], axis=0)
                #smote_input_params = np.append(smote_input_params,good_params[np.random.randint(len(good_params),size=self.smote_neighbours-ngood)],axis=0)
                smote_input_params = np.concatenate([validated_params[new_labels == 0], validated_params[new_labels == 1], good_params[np.random.randint(len(good_params), size=(2*self.smote_neighbours - ngood))]], axis=0)

                print(f"bad: {nbad}, good: {ngood}, smote size: {smote_input_params.shape}, neighbours: {self.smote_neighbours}")
                #oversample = SMOTE(k_neighbors=self.smote_neighbours, n_jobs=self.ncores)
                oversample = SMOTE(k_neighbors=self.smote_neighbours)
                synthetic_params, _ = oversample.fit_resample(smote_input_params, smote_input_labels)

                self.RunManager.run_batch(synthetic_params.tolist())

                smote_validated_params, smote_obs_values = self.extract_from_valid_points(self.RunManager.valid_batch_points)
                smote_labels = labeler(smote_obs_values, self.threshold_value)


                good_obs_values = np.append(good_obs_values, new_obs_values[new_labels == 1])
                good_params = np.append(good_params, validated_params[new_labels == 1], axis=0)
                good_labels = np.append(good_labels, new_labels[new_labels == 1], axis=0)
                good_obs_values = np.append(good_obs_values, smote_obs_values[smote_labels == 1])
                good_params = np.append(good_params, smote_validated_params[smote_labels == 1], axis=0)
                good_labels = np.append(good_labels, smote_labels[smote_labels == 1], axis=0)
                bad_params = np.append(bad_params, validated_params[new_labels == 0][:len(smote_labels[smote_labels == 1])], axis=0)
                bad_labels = np.append(bad_labels, new_labels[new_labels == 0][:len(smote_labels[smote_labels == 1])], axis=0)
                collected_params = np.concatenate([smote_validated_params[smote_labels == 1], validated_params[new_labels == 0][:len(smote_labels[smote_labels == 1])]], axis=0)
                balanced_training_labels = np.concatenate([smote_labels[smote_labels == 1], new_labels[new_labels == 0][:len(smote_labels[smote_labels == 1])]], axis=0)

            elif (len(new_labels[new_labels == 1]) < 0.5 * round(self.points_per_iter * (1 - self.random_fraction)) and len(new_labels[new_labels == 1]) > self.function_dim):
                # Need to check what this is supposed to do and make sure it is doing it!!

                # more garbage!!
                #smote_input_labels = np.append(np.zeros(len(new_labels[new_labels == 0])), np.ones(len(new_labels[new_labels == 1])), axis=0)
                #smote_input_params = np.append(candidate_params[:len(new_labels[new_labels == 0])], validated_params[new_labels == 1], axis=0)
                smote_input_labels = np.append(np.zeros(len(new_labels[new_labels == 0])), np.ones(len(new_labels[new_labels == 1])), axis=0)
                smote_input_params = np.append(validated_params[new_labels == 0], validated_params[new_labels == 1], axis=0)
                #oversample = SMOTE(k_neighbors=self.smote_neighbours, n_jobs=self.ncores)
                oversample = SMOTE(k_neighbors=self.smote_neighbours)
                synthetic_params, _ = oversample.fit_resample(smote_input_params, smote_input_labels)

                self.RunManager.run_batch(synthetic_params.tolist())
                smote_validated_params, smote_obs_values = self.extract_from_valid_points(self.RunManager.valid_batch_points)
                smote_labels = labeler(smote_obs_values, self.threshold_value)

                good_obs_values = np.append(good_obs_values, new_obs_values[new_labels == 1])
                good_params = np.append(good_params, validated_params[new_labels == 1], axis=0)
                good_labels = np.append(good_labels, new_labels[new_labels == 1], axis=0)
                good_obs_values = np.append(good_obs_values, smote_obs_values[smote_labels == 1])
                good_params = np.append(good_params, smote_validated_params[smote_labels == 1], axis=0)
                good_labels = np.append(good_labels, smote_labels[smote_labels == 1], axis=0)
                bad_params = np.append(bad_params, validated_params[new_labels == 0][:len(smote_labels[smote_labels == 1])], axis=0)
                bad_labels = np.append(bad_labels, new_labels[new_labels == 0][:len(smote_labels[smote_labels == 1])], axis=0)
                collected_params = np.concatenate([smote_validated_params[smote_labels == 1], validated_params[new_labels == 0][:len(smote_labels[smote_labels == 1])]], axis=0)
                balanced_training_labels = np.concatenate([smote_labels[smote_labels == 1], new_labels[new_labels == 0][:len(smote_labels[smote_labels == 1])]], axis=0)

            else:  # shouldn't we check that we include enough bad points?
                good_obs_values = np.append(good_obs_values, new_obs_values[new_labels == 1])
                good_params = np.append(good_params, validated_params[new_labels == 1], axis=0)
                good_labels = np.append(good_labels, new_labels[new_labels == 1], axis=0)
                bad_params = np.append(bad_params, validated_params[new_labels == 0][:len(validated_params[new_labels == 1])], axis=0)
                bad_labels = np.append(bad_labels, new_labels[new_labels == 0][:len(validated_params[new_labels == 1])], axis=0)
                collected_params = np.concatenate([validated_params[new_labels == 1], validated_params[new_labels == 0]], axis=0)
                balanced_training_labels = np.concatenate([new_labels[new_labels == 1], new_labels[new_labels == 0]], axis=0)

            # Check if we've reached target before expensive retraining
            if len(good_labels) >= self.target_points:
                self.log.info(f'Target reached ({len(good_labels)} >= {self.target_points}). Skipping retraining.')
            else:
                collected_params = np.concatenate([good_params, bad_params[:len(good_params)]], axis=0)
                balanced_training_labels = np.concatenate([good_labels, bad_labels[:len(good_labels)]], axis=0)
                shuffled_params, shuffled_labels = sklearn.utils.shuffle(collected_params, balanced_training_labels)

                if shuffled_params.shape[0] == 0:
                    continue

                torchfit(shuffled_params, shuffled_labels, dnn_model, criterion, optimizer, self.epoch, self.batch_size, verbose = True)

            new_good_count = len(good_params) - previous_good_count

            if new_good_count > 0:
                new_good_params = good_params[previous_good_count:]
                new_good_obs = good_obs_values[previous_good_count:]

                new_good_data = pd.DataFrame(new_good_params, columns=var_names)
                new_good_data[self.primary_observable] = new_good_obs
                new_good_data.to_csv(accumulated_csv_path, mode='a', header=False, index=False)

            self.log.info(f'DNN_model - Run Number {run_number} - Found {new_good_count} new good points - Total collected: {len(good_params)}')

        self.log.info(f'Reached target! Total good points collected: {len(good_params)}')
        self.log.info(f'Accumulated points saved to: {accumulated_csv_path}')

        return good_params, good_obs_values
