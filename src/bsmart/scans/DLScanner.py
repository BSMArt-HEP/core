"""
DLScanner 
---------


DLScanner method `DLScanner` based on code from 

* `<https://github.com/raalraan/DLScanner/>`_
* `<https://arxiv.org/abs/2412.19675>`_

Adapted by Mark Goodsell and Farid Ibrahimov. 
We have tried to stay as close to the original as possible. However, important differences are:
* Instead of comparing directly to a threshold value of an observable, we compare to the negative log likelihood of the model (similarly to the MLScanner methods). This is more general and hopefully more useful, since it allows many observables to be compared at once. In order to directly compare to the value of the observable, you can set "SCALING": "EXPUSER". In order to consider the value to be lower than the threshold in that case, you would need to define a function to be -1*observable_value.
* We use pytorch instead of keras. 

"""


__meta__ = {
 "name": "DLScanner",
 "requires": ["vegas","torch"],
 "settings": {
    "Networks": {
        "Iterations": "Number of active learning iterations (default: 100).",
        "Method": "Method to use: 'Classifier' or 'Regressor' (default: 'Classifier').",
        "K": "Number of points to generate in each iteration (default: 10).",
        "L": "Number of candidate points to test with the model in each iteration (default: 500).",
        "Threshold": "Negative Log Likelihood threshold below which a point is classified as 'good' (default: 10).",
        "Hidden_Layers": "Number of hidden layers in the DNN (default: 3).",
        "Neurons": "Number of neurons in the hidden layers (default: 100).",
        "Epochs": "Number of training epochs (default: 250).",
        "Batch Size": "Batch size for training (default: 32).",
        "LearningRate": "Learning rate for the optimizer (default: 1e-2).",
        "WeightDecay": "Weight decay for the optimizer (default: 1e-4).",
        "Use_Vegas_Map": "Whether to use Vegas map for sampling (default: True).",
        "Vegas_Frac": "Fraction of points to sample using Vegas map (default: 0.1).",
        "Vegas_Ninc": "Number of increments for Vegas map (default: 100).",
        "Random_Points": "Number of random points to add in each iteration (default: K/10).",
        "Verbose": "Verbosity level (default: True)."
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


import random as random
import math
from bsmart.BSMlikelihood import MakeGlobalLikelihood

from bsmart.ml.seed_points import get_previous_sample
from bsmart import BSMutil

#import vegas

#print(f"CWD: {os.getcwd()}")
from bsmart.ml.dls_vegas import vegas_map_samples
import bsmart.ml.DLML as DLML
from bsmart.ml.DLML import NumpyScaler

def torchfit(X, y, model, criterion, optimizer, epochs, batch_size, verbose=True, patience=50, min_delta=1e-5):
    X_tensor = torch.FloatTensor(X)
    y_tensor = torch.FloatTensor(y).reshape(-1, 1)
    #print(f"Input X_tensor shape: {X_tensor.shape}")
    #print(f"Loss y_tensor shape: {y_tensor.shape}")
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
        self.citations = """@article{Hammad:2024tzz,
author = "Hammad, A. and Ramos, Raymundo",
title = "{DLScanner: A parameter space scanner package assisted by deep learning methods}",
eprint = "2412.19675",
archivePrefix = "arXiv",
primaryClass = "hep-ph",
doi = "10.1016/j.cpc.2025.109659",
journal = "Comput. Phys. Commun.",
volume = "314",
pages = "109659",
year = "2025"
}
"""


    def __init__(self, inputs, log):
        Scan.__init__(self, inputs, log)
        self.log.info("Start DLScanner init")
        

        self.downscalers,self.upscalers = BSMutil.create_scalers(self.inputs,log)
        self.n_variables = len(self.inputs['Variables'])
        self.ndim = self.n_variables
        hlayers = self.inputs['Networks'].get('Hidden_Layers', 3)
        neurons = self.inputs['Networks'].get('Neurons', 100)
        self.method = self.inputs['Networks'].get('Method', 'Classifier')
        self.K = self.inputs['Networks'].get('K', 10)
        self.L = self.inputs['Networks'].get('L', 500)

        self.outdim = 1 # We only predict the NLL


        self.verbose = eval(self.inputs['Networks'].get('Verbose', "True"))
        self.epoch = int(self.inputs['Networks'].get('Epochs', 250))
        self.learning_rate = float(self.inputs['Networks'].get('LearningRate', 1e-2))
        self.weight_decay = float(self.inputs['Networks'].get('WeightDecay', 1e-4))
        self.threshold_suggest = float(self.inputs['Networks'].get('Threshold', 10))
        self.vegas_ninc = int(self.inputs['Networks'].get('Vegas_Ninc', 100))
        self.use_vegas_map = eval(self.inputs['Networks'].get('Use_Vegas_Map', "True"))
        self.vegas_frac = float(self.inputs['Networks'].get('Vegas_Frac', 0.1))
        self.steps = int(self.inputs['Networks'].get('Iterations', 100))
        self.batch_size = int(self.inputs['Networks'].get('Batch Size', 32))
        self.randpts = int(self.inputs['Networks'].get('Random_Points', self.K/10 + 0.5))

        mlinputs={'input_size':self.ndim,'hidden_layers':hlayers,'hidden_size':neurons}

        self.classify = False
        if self.method == "Classifier":
            self.classify = True
            self.model = DLML.MLP_Classifier(mlinputs)
            
            #self.loss = "binary_crossentropy"
            self.criterion = nn.BCELoss()

        elif self.method == "Regressor":
            self.model = DLML.MLP_Regressor(mlinputs)
  
            #self.loss = "mean_absolute_error"
            self.criterion = nn.MSELoss(reduction='mean')
        #print(f"Initial weights: {list(self.model.parameters())}")
        #print("Got here")
        ###########################################
        # self.iteration = int(self.inputs['Networks'].get('Iterations', 10))
        
        # self.candidate_points = int(self.inputs['Networks'].get('Candidate_Points', 500))
        # self.bootstrap_points = int(self.inputs['Networks'].get('Bootstrap_Points', 100))
        # self.points_per_iter = int(self.inputs['Networks'].get('Points_Per_Iteration', 300))
        # self.smote_neighbours = int(self.inputs['Networks'].get('Kinitial', 5))
        # #self.threshold_value = float(self.inputs['Networks'].get('Threshold_Value', 100))
        
        # #self.target_points = int(self.inputs['Networks'].get('Target_Points', 20000))
        # self.function_dim = int(len(self.inputs['Variables']))
        # self.random_fraction = float(self.inputs['Networks'].get('Random_Fraction', 0.2))
        # self.batch_size = int(self.inputs['Networks'].get('Batch Size', 500))
        # self.neurons = int(self.inputs['Networks'].get('Neurons', 100))
        
        
        # self.n_estimators = int(self.inputs['Networks'].get('Estimators', 100))
        # self.max_depth = int(self.inputs['Networks'].get('Max_Depth', 30))
        
        #print(f"verbose: {self.verbose}")
        if "Cores" in self.inputs['Setup']:
            self.ncores = int(self.inputs['Setup']['Cores'])
            torch.set_num_threads(min(8, self.ncores))

        #self.num_points = int(self.inputs['Setup'].get('Points', 10000))

        self.target_points = int(self.inputs['Setup'].get('Points', 10000))
        

        #self.loglike = eval(self.inputs['Setup'].get('LogLike',"False"))
        # choose primary observable from likelihoods

        #print("make likelihoods")
        self.maxloss =  np.log(1 + np.finfo(np.float64).max) + 1

        # Create the likelihood functions based on 'scaling' property of variables, if they have it; otherwise use gaussian likelihood    
        #elf.get_likelihood,self.observable_masks=MakeGlobalLikelihood(self.inputs['Observables'],self.use_loglike)
        #self.get_likelihood,self.observable_masks=MakeGlobalLikelihood(self.inputs['Observables'],return_type="NLL") # need a negative log likelihood
        #self.likelihood_fns, self.observable_masks = MakeLikelihoods(self.inputs["Observables"], loglike=True)
        if self.naive:
            self.getNLL,self.masks = MakeGlobalLikelihood(self.inputs['Observables'],return_type="NLL", smooth_losses= False) 
        else:
            self.getNLL,self.masks = MakeGlobalLikelihood(self.inputs['Observables'],return_type="NLL", smooth_losses= True, max_loss=self.maxloss)
        ##
        #self.likelihood, self.observable_masks = MakeGlobalLikelihood(self.inputs["Observables"], self.loglike)
        bounds=np.zeros((self.n_variables,2))
        bounds[:,1] = np.ones(self.n_variables)
        self.limits_arr = bounds

        self.rng = np.random.default_rng()
        self.samples0 = [] 
        self.nlls0 = []
        if 'Initial Sample' in self.inputs['Setup'] and os.path.exists(self.inputs['Setup']['Initial Sample']):
            self.variables0, self.nlls0 = get_previous_sample(self.inputs['Setup']['Initial Sample'],returntype="NLL")

            self.samples0 = np.array([[f(y) for y,f in zip(pt,self.downscalers)] for pt in self.variables0])

        self.log.info("finished scan init")


    def postprocess(self,Point, observables, data_point,temp_dir,log, lock=None):
        """ 
        If we are keeping 'invalid' points, the NLL is computed separately, so we just return 1.0 for valid points
        Otherwise the return value is the NLL.
        """

        # For this algorithm we are minimising the objective function
        #res=1e100
        if self.naive:
            return self.getNLL(observables)
        else:
            return 1.0
   
        

    def extract_from_valid_points(self,valid_points):
        params = np.array([point['inputs'] for point in valid_points])
        batch_observables = [point["observables"] for point in valid_points]

        losses = np.array([self.get_losses(obs) for obs in batch_observables])
        fitnesses = np.sum(losses, axis=1)

        return params, fitnesses
    

    def run_batch(self,thetas):
        # Scale thetas, run the batch, return the likelihoods
            
        batch = [ [f(y) for y,f in zip(list(theta),self.upscalers)] for theta in thetas]
        results = self.RunManager.run_batch(batch)
            
        if self.naive: # postprocess returns the NLL, or [] if invalid
            
            losses = np.array([value if isinstance(value, (int, float)) else self.maxloss for value in results])

            #batch_inputs = []
            #for cube_input,value in zip(thetas,results):
            #    if isinstance(value, (int, float)):
            #        batch_inputs.append(cube_input)
                
        else: # keeping 'invalid' points
            batch_observables = np.array([point["observables"] for point in self.RunManager.all_batch_points])
            losses = np.array([self.getNLL(obs) for obs in batch_observables])
            #batch_inputs = thetas
        #print(f"losses: {losses}, of size {len(losses)}")

        if self.classify:
            losses = np.array([1 if loss < self.threshold_suggest else 0 for loss in losses])

        return thetas,losses # want to keep the values of observables, only possible if not naive


    def genrand(self, npts):
        dim = self.ndim
        limits_arr = self.limits_arr
        p0 = self.rng.uniform(
            low=limits_arr[:, 0],
            high=limits_arr[:, 1],
            size=[npts, dim]
        )
        return p0

    
    def run(self):


        if len(self.samples0) == 0:
            self.samples0 = self.genrand(self.K)

            xinit,outinit = self.run_batch(self.samples0)

            if len(outinit) == 0:
                sys.exit("No valid points found in initial batch!")
        else:
            xinit = self.samples0
            outinit = self.nlls0

        #print(f"Initial batch: {self.samples0,xinit,outinit}")
        #self.criterion = self.loss_fn_selector(outinit)
        # sel#
        #self.optimizer = torch.optim.AdamW(self.model.parameters(), lr=self.learning_rate, weight_decay=1e-4)
        self.optimizer= torch.optim.AdamW(self.model.parameters(), lr=self.learning_rate, weight_decay=self.weight_decay)
        

        self.model.train()

        if not self.classify:
            self.scaler = NumpyScaler()
            self.scaler.fit(outinit)
            outinit_scaled = self.scaler.transform(outinit)
        else:
            outinit_scaled = outinit
        initial_loss=torchfit(xinit, outinit_scaled, self.model, self.criterion, self.optimizer, self.epoch, self.batch_size, verbose=True)

        self.samples = self.samples0
        self.samples_out = outinit
        self.histories = []
        for j in range(self.steps):
            xsug, llsug = self.suggestpts(vegas_ninc=self.vegas_ninc)

            if len(xsug) == 0 or len(llsug) == 0:
                self.log.warning("No valid points found in batch!")
                continue
            self.samples = np.append(self.samples, xsug, axis=0)
            #self.samples_list.append(xsug)
            #print(f"samples out: {self.samples_out}, llsug: {llsug}")
            self.samples_out = np.append(self.samples_out, llsug, axis=0)
            #self.samples_out_list.append(llsug)

            if self.classify:
                nintgt = (self.samples_out > 0.5).sum()
                xtrain_here, ytrain_here = self.class_equalize()
                if self.verbose:
                    print(
                        "Step {}, Number of points in-target: {}".format(
                            j + 1, nintgt
                        )
                    )

                
            else:
                xtrain_here = self.samples
                ytrain_here = self.samples_out

            self.model.reset_weights()
            self.optimizer = torch.optim.AdamW(self.model.parameters(), lr=self.learning_rate, weight_decay=self.weight_decay)
            if not self.classify:
                self.scaler = NumpyScaler()
                self.scaler.fit(ytrain_here)
                ytrain_here_scaled = self.scaler.transform(ytrain_here)
            else:
                ytrain_here_scaled = ytrain_here
            
            self.histories.append(
                torchfit(
                    xtrain_here, ytrain_here_scaled,self.model, self.criterion, self.optimizer,
                    epochs=self.epoch, batch_size=self.batch_size, verbose=True
                ))

        self.log.info(f"Scan finished after {j} steps with best loss {min(self.histories)}")


    def suggestpts(
        self,
        npts=None, randpts=None, L=None, limits=None,
        verbose=None, use_vegas_map=None, vegas_frac=None,
        threshold=None, vegas_ninc=None,
    ):
        # randpts: number of points chosen at random that will be
        #     added to the suggested points, randpts < self.K
        # L: number of points that will be used to test the
        #     model in order to obtain suggested points,
        #     L > self.K
        if npts is None:
            npts = self.K
        if randpts is None:
            randpts = self.randpts
        if L is None:
            L = self.L
        if verbose is None:
            verbose = self.verbose
        if use_vegas_map is None:
            use_vegas_map = self.use_vegas_map
        if vegas_frac is None:
            vegas_frac = self.vegas_frac
        if threshold is None:
            threshold = self.threshold_suggest
        if vegas_ninc is None:
            vegas_ninc = self.vegas_ninc
        _randpts = randpts
        # Try to  predict the observable for several points using what the
        # machine learned
        xtry = self.genrand(L)

        if use_vegas_map:
            try:
                vegas_pts = int(vegas_frac*L + 0.5)
                # Train a vegas map
                if verbose:
                    print(f"Training vegas map to obtain {vegas_pts} points using accumulated samples")
                map_vg = vegas_map_samples(
                    self.samples, self.samples_out.flatten(), self.limits_arr,
                    ninc=vegas_ninc
                )
                self.vegas_map_gen = map_vg
                # Discard jacobian
                xtry_vg, _ = map_vg(L)
                #print(f"Vegas suggested points: {xtry_vg}")
                xtry = np.concatenate([
                    xtry[:L - vegas_pts],
                    xtry_vg[:int(vegas_frac*L)]
                ])
            except Exception as e: 
                self.log.error(f"Error in vegas map: {e}, using random points")
                pass

        self.model.eval()
        if not self.classify:
            ptry = self.scaler.inverse_transform(self.model(torch.tensor(xtry, dtype=torch.float32)).detach().numpy())
        else:
            ptry = self.model(torch.tensor(xtry, dtype=torch.float32)).detach().numpy()
        
        #print(f"xtry: {xtry}, ptry: {ptry}")
        
        #xcand = xtry[(ptry > threshold).flatten()]

        if self.classify:
            xcand = xtry[(ptry > 0.5).flatten()]
        else:
            # using NLL so lower values are better
            xcand = xtry[(ptry < threshold).flatten()]
        if xcand.shape[0] < npts - randpts:
            _randpts = npts - xcand.shape[0]
            if verbose:
                print(
                    "Tried {} points, {} points survived, filling with random points".format(
                        L, xcand.shape[0]
                        )
                    )
        else:
            if verbose:
                print(
                    "Tried {} points, {} points survived, choosing {} points".format(
                        L, xcand.shape[0], npts
                        )
                    )

        # Use the points according to the class predicted by the model
        # but pass the correct class. In this way, the points that the
        # model got wrong should be corrected
        xsel = xcand[:npts - _randpts]
        # TODO Check what would be good/interesting stats
        # check1 = (llsel < 100 + 2*10)*(llsel > 100 - 2*10)
        # check2 = ((ptry[fltrcand][:npts - _randpts] - llsel)
        #           / (ptry[fltrcand][:npts - _randpts] + llsel))
        # stats1 = np.sum(check1)/llsel.shape[0]
        # stats2 = np.sum(check2 < 1.e-3)/llsel.shape[0]
        # stats3 = np.sum(check1*(check2 < 1.e-3))/ll1sel.shape[0]

        # Append `_randpts` more points chosen at random
        xnew = self.genrand(_randpts)
        xout = np.append(xsel, xnew, axis=0)
        # llout = np.empty((npts, 1))
        # for k in range(npts):
        #     llout[k] = np.array(self.ll(xout[k]))


        #llout = np.array(self.ll(xout))
        if self.classify:
            ptry_keep = ptry[(ptry > 0.5).flatten()][:len(xsel)]
        else:
            # using NLL so lower values are better
            ptry_keep = ptry[(ptry < threshold).flatten()][:len(xsel)]
        
        xout_keep,llout = self.run_batch(xout)

        # print losses for comparison
        # print(f"Compare losses to predictions {llout[:len(xsel)].shape}, {ptry_keep.shape}")
        # print(np.hstack([llout[:len(xsel)].reshape(-1,1),ptry_keep.reshape(-1,1)]))
        #return xout_keep, llout.reshape(xout_keep.shape[0], self.outdim)

        self.log.info(f"Compare losses to predictions {llout[:len(xsel)].shape}, {ptry_keep.shape}")
        self.log.info(np.hstack([llout[:len(xsel)].reshape(-1,1),ptry_keep.reshape(-1,1)]))
        return xout_keep, llout.reshape(xout_keep.shape[0])  # dimension 1 output

    def class_equalize(self):
        ccount = [
            (self.samples_out == 0.0).sum(),
            (self.samples_out == 1.0).sum()
        ]

        xsam0 = self.samples[self.samples_out.flatten() == 0.0]
        xsam1 = self.samples[self.samples_out.flatten() == 1.0]
        ysam0 = self.samples_out[self.samples_out.flatten() == 0.0]
        ysam1 = self.samples_out[self.samples_out.flatten() == 1.0]

        xtrain_here = np.empty((0, self.ndim))
        #ytrain_here = np.empty((0, self.outdim))
        ytrain_here = np.empty((0))

        if ccount[0] != ccount[1]:
            if ccount[1] < ccount[0]:
                while xtrain_here.shape[0] < ccount[0]:
                    xtrain_here = np.append(
                        xtrain_here,
                        xsam1,
                        axis=0
                    )
                    ytrain_here = np.append(
                        ytrain_here,
                        ysam1,
                        axis=0
                    )
                xtrain_here = np.append(
                    xtrain_here[:ccount[0]],
                    xsam0,
                    axis=0
                )
                ytrain_here = np.append(
                    ytrain_here[:ccount[0]],
                    ysam0,
                    axis=0
                )
            elif ccount[0] < ccount[1]:
                while xtrain_here.shape[0] < ccount[1]:
                    xtrain_here = np.append(
                        xtrain_here,
                        xsam0,
                        axis=0
                    )
                    ytrain_here = np.append(
                        ytrain_here,
                        ysam0,
                        axis=0
                    )
                xtrain_here = np.append(
                    xtrain_here[:ccount[1]],
                    xsam1,
                    axis=0
                )
                ytrain_here = np.append(
                    ytrain_here[:ccount[1]],
                    ysam1,
                    axis=0
                )
            else:
                xtrain_here = self.samples
                ytrain_here = self.samples_out

        return xtrain_here, ytrain_here