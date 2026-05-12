"""
Optimisation using CMAES with Novelty Detection


Requires `pip install cmaes`, see:
https://github.com/CyberAgentAILab/cmaes


BSMArt scan written by ....



This scan uses cmaes to optimise a log-likelihood constructed from selected observables.

All desired observables must use the scaling CFUNCTION with a range of values for which the observable is in an acceptable range (other observables should have "SCALING": "OFF" or at least no MEAN/VARIANCE/RANGE set so that they will be ignored.)

It is recommended to use this also with a validity condition, so that no further codes are evaluated when the observables are outside of this range, and the code chain can focus on first obtaining acceptable values for that/those observables.

For example, here we tell BSMArt to run SPheno and check the Higgs mass lies within 3 GeV of its observed value. We also record the lightest neutralino and selectron masses, construct a derived observable based on the difference of their masses, and insist that it is greater than zero:

.. code-block:: json

    "Codes": {
            "SPheno": {
                "Command": "< >",
                "InputFile": "< >",
                "OutputFile": "< >",
                "Observables":
                {
                            "mh": {
                            "SLHA": [
                                    "MASS",
                                    [25]
                            ],
                            "SCALING": "CFUNCTION",
                            "MEAN": 125.09,
                            "VARIANCE": 3,
                            "VALID": [122.09,128.09]
                            },
                            "mchi1": {
                                    "SLHA" : [ "MASS", [1000022]], "SCALING": "OFF"]
                            },
                            "mse1": {
                                    "SLHA" : [ "MASS", [1000011]], "SCALING": "OFF"]
                            },

                            "dmse1mchi1" : {
                            "FUNCTION": "mse1-mchi1",
                            "SCALING": "CFUNCTION",
                            "RANGE": [0,100000]
                            }
                 },
                 "Conditions": ["mse1>mchi1"],
                 "Run": "True"
            }











"""

__meta__ = {
    "name": "CMAES_ND",
    "requires": ["cmaes", "pyod", "sklearn", "traceback", "pandas"],
    "settings": {
        "CMAESPopulationSize": "Number of points to run per batch (optional, auto-calculated if not set)",
        "CMAESMean": "Initial mean for CMA-ES optimizer (optional, random if not set)",
        "CMAESSigma": "Initial step size for CMA-ES (optional, default: 1/sqrt(n_variables))",
        "CMAESMaxGenerations": "Number of batches to run (default: 100)",
        "CMAES MaxLoss": "Maximum loss value override (optional)",
        "CMAESPathSeeds": "Path to seed file (.csv or .parquet) for initial mean (optional)",
        "CMAESNoveltyDetection": "Enable novelty detection using HBOS (default: True)",
        "CMAESHBOSNBins": "Number of bins for HBOS novelty detection (default: 100)",
        "CMAESFocus": "Dict with Variables and Observables lists for novelty detection focus (optional)",
        "CMAESHierarchy": "List of hierarchical observables to focus on first (optional)",
        "CMAESNoGoodPointsPatience": "Number of generations without good points before stopping (default: 5000)",
        "CMAESGoodPointsEarlyStop": "Target number of good points for early stopping (default: 5000)",
        "CMAESBestLossPatience": "Patience for best loss counter (default: 200)",
        "CMAESBestNValidConstraintsPatience": "Patience for best n valid constraints counter (default: 2000)",
        "Initial Sample": "Path to initial sample file for warm start (optional)",
    },
}

import math
import os
import shutil
import sys
import traceback
from collections import OrderedDict

from bsmart import BSMutil
import cmaes
from bsmart import debug
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# from BSMlikelihood import MakeLikelihoods, MakeGlobalLikelihood
from bsmart.BSMlikelihood import MakeLikelihoods, safe_float
from cmaes import CMA
from bsmart.core import Scan as Scan
from pyod.models.hbos import HBOS
from sklearn.preprocessing import minmax_scale


class HBOSPenaliser:
    def __init__(self, n_bins=100, *args, **kwargs):
        self.n_bins = n_bins
        self.hbos = HBOS(n_bins=self.n_bins)

    def fit(self, X):
        self.hbos = self.hbos.fit(X)
        return self

    def get_penalties(self, x):
        try:
            x = np.float64(x)
            _x = np.nan_to_num(x, -np.inf)
            scores = self.hbos.predict_proba(_x)
        except Exception as e:
            print("-" * 60)
            traceback.print_exc(file=sys.stdout)
            print("-" * 60)
            print(x)
            print(_x)
            exit(1)
        return scores[:, 0]


class Optimiser:
    def __init__(self, direction="min", *args, **kwargs):
        self.direction = direction
        self.args = args
        self.kwargs = kwargs

    def _sigmoid(self, x):
        return 1 / (1 + np.exp(-x))

    def fit(self, X):
        return self

    def get_penalties(self, X):
        factor = 1 if self.direction == "min" else -1
        _X = minmax_scale(X, feature_range=(-1, 1))
        return self._sigmoid(factor * _X).mean(axis=1)


PENALTIES = {"HBOS": HBOSPenaliser, "Optimiser": Optimiser}


class NewScan(Scan):
    """Scanner class for CMAES-ND Scans"""

    def initialise(self):  # overload initial settings.
        self.runsettings.store_points_in_memory = (
            True  # Force storage of info for hierarchical loss computation
        )
        self.runsettings.store_invalid_points = True  # Need this to keep points where
        self.runsettings.invalid_return_value = 0
        self.citations = """@article{deSouza:2025uxb,
    author = "de Souza, Fernando Abreu and Castro, Nuno Filipe and Crispim Rom{\\~a}o, Miguel and Porod, Werner",
    title = "{Exploring scotogenic parameter spaces and mapping uncharted dark matter phenomenology with multi-objective search algorithms}",
    eprint = "2505.08862",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    reportNumber = "IPPP/25/29",
    doi = "10.1007/JHEP10(2025)116",
    journal = "JHEP",
    volume = "10",
    pages = "116",
    year = "2025"
}
"""

    def __init__(self, inputs, log):
        Scan.__init__(self, inputs, log)

        self.downscalers, self.upscalers = BSMutil.create_scalers(self.inputs, log)
        self.n_variables = len(self.upscalers)
        self.use_loglike = eval(self.inputs["Setup"].get("LogLike", "True"))
        self.Points = self.inputs["Setup"].get("Points", 10000)

        self.Variables = self.inputs["Variables"]

        self.maxloss = np.log(1 + np.finfo(np.float64).max) + 1
        if "CMAES MaxLoss" in self.inputs["Setup"]:
            if self.inputs["Setup"].get("CMAES MaxLoss") is not None:
                self.maxloss = float(self.inputs["Setup"].get("CMAES MaxLoss"))

        # Create the likelihood functions based on 'scaling' property of variables, if they have it; otherwise use gaussian likelihood
        # self.get_likelihood,self.observable_masks=MakeGlobalLikelihood(self.inputs['Observables'],self.use_loglike)

        # self.likelihood_fns,self.observable_masks=MakeLikelihoods(self.inputs['Observables'],self.use_loglike)
        self.likelihood_fns, self.observable_masks = MakeLikelihoods(
            self.inputs["Observables"], loglike=True
        )
        self.Observables = list(self.inputs["Observables"].keys())
        self.observable_indices = {obs: idx for idx, obs in enumerate(self.Observables)}
        self.Constraints = [
            obs for obs, mask in zip(self.Observables, self.observable_masks) if mask
        ]
        # Build constraint_likelihoods by matching observables to likelihood functions
        # Note: MakeLikelihoods only creates likelihood functions for non-OFF observables,
        # so we need to match by iterating through observables in the same order
        likelihood_fn_iter = iter(self.likelihood_fns)
        self.constraint_likelihoods = []
        for obs, mask in zip(self.Observables, self.observable_masks):
            if mask:
                # This observable should have a likelihood function
                try:
                    self.constraint_likelihoods.append(next(likelihood_fn_iter))
                except StopIteration:
                    self.log.error(f"Missing likelihood function for constraint: {obs}")
                    raise ValueError(
                        f"Missing likelihood function for constraint: {obs}. "
                        f"This suggests a mismatch between observable masks and likelihood functions."
                    )

        ## Debugging: Log constraint analysis
        self.log.info("=" * 80)
        self.log.info("CMAES-ND: Constraint Analysis")
        self.log.info("=" * 80)
        self.log.info(f"Total observables found: {len(self.Observables)}")
        self.log.info("")
        self.log.info("Observable details:")
        for obs, mask in zip(self.Observables, self.observable_masks):
            scaling = self.inputs["Observables"][obs].get("SCALING", "NOT SET")
            mask_status = "ACTIVE (mask=True)" if mask else "INACTIVE (mask=False)"
            self.log.info(f"  - {obs}: SCALING={scaling}, {mask_status}")
        self.log.info("")
        self.log.info(f"Active constraints (mask=True): {len(self.Constraints)}")
        self.log.info(f"Constraint names: {self.Constraints}")
        self.log.info(
            f"constraint_likelihoods count: {len(self.constraint_likelihoods)}"
        )
        if len(self.Constraints) != len(self.constraint_likelihoods):
            self.log.error(
                f"ERROR: Length mismatch! Constraints={len(self.Constraints)}, "
                f"constraint_likelihoods={len(self.constraint_likelihoods)}"
            )
        self.log.info("")
        self.log.info("Constraint-to-likelihood mapping:")
        for i, (obs_name, likelihood_fn) in enumerate(
            zip(self.Constraints, self.constraint_likelihoods)
        ):
            self.log.info(
                f"  [{i}] {obs_name} -> likelihood_fn type: {type(likelihood_fn).__name__}"
            )
        self.log.info("=" * 80)
        self.log.info("")

        ## Check that the observables used as constraints are defined with the CFUNCTION method.

        obsok = True
        for obs in self.Constraints:
            if self.inputs["Observables"][obs]["SCALING"] != "CFUNCTION":
                self.log.error(
                    "Observable "
                    + obs
                    + " should be defined with the CFUNCTION scaling method. Cannot proceed with the scan!"
                )
        if not obsok:
            raise NameError("Incorrectly defined observables")

        self.population_size = None
        if "CMAESPopulationSize" in self.inputs["Setup"]:
            if self.inputs["Setup"].get("CMAESPopulationSize") is not None:
                self.population_size = int(
                    self.inputs["Setup"].get("CMAESPopulationSize")
                )
        self.mean = np.random.rand(self.n_variables)
        if "CMAESMean" in self.inputs["Setup"]:
            if self.inputs["Setup"].get("CMAESMean") is not None:
                self.mean = float(self.inputs["Setup"].get("CMAESMean")) * np.ones(
                    self.n_variables
                )

        self.sigma = 1.0 / math.sqrt(self.n_variables)
        if "CMAESSigma" in self.inputs["Setup"]:
            if self.inputs["Setup"].get("CMAESSigma") is not None:
                self.sigma = float(self.inputs["Setup"].get("CMAESSigma"))

        if "CMAESPathSeeds" in self.inputs["Setup"]:
            path_seeds = self.inputs["Setup"].get("CMAESPathSeeds")
            try:
                if ".csv" in path_seeds:
                    centroid_seeds = pd.read_csv(path_seeds)
                else:
                    centroid_seeds = pd.read_parquet(path_seeds)
                print(f"Using seeds in {path_seeds}")
            except:
                raise Exception("Seeds file for CMAES must be .csv or .parquet!")
            # Extract values in original space from the seed file
            seed_values_original = np.array(
                centroid_seeds.sample(n=1)[self.inputs["Variables"].keys()]
                .iloc[0]
                .to_list()
            )
            # Apply downscalers to convert from original space to normalized [0,1] space for CMA-ES
            self.mean = np.array(
                [f(x) for x, f in zip(seed_values_original, self.downscalers)]
            )
        self.novelty_detection = True
        if "CMAESNoveltyDetection" in self.inputs["Setup"]:

            self.novelty_detection = eval(
                str(self.inputs["Setup"]["CMAESNoveltyDetection"])
            )

        if self.novelty_detection:
            # self.hbos_n_bins = int(self.inputs["Setup"].get("CMAESHBOSNBins", 100))
            # self.nd_estimator = HBOSPenaliser(n_bins=self.hbos_n_bins)
            self.nd_estimator = PENALTIES[
                self.inputs["Setup"].get("CMAESNoveltyDetection", "HBOS")
            ](**self.inputs["Setup"].get("CMAESNoveltyDetectionParams", {}))
            self.log.info(f"Using novelty detection estimator: {self.nd_estimator}")

        if "CMAESFocus" in self.inputs["Setup"]:
            self.focus = self.inputs["Setup"]["CMAESFocus"]
        else:
            self.focus = {
                "Variables": list(self.Variables.keys()),
                "Observables": list(self.Observables),
            }
        if "CMAESMaxGenerations" in self.inputs["Setup"]:
            self.max_n_generations = self.inputs["Setup"]["CMAESMaxGenerations"]
        else:
            self.max_n_generations = 100
        if (
            "CMAESHierarchy" in self.inputs["Setup"]
            and len(self.inputs["Setup"]["CMAESHierarchy"]) > 0
        ):
            self.use_hierarchy = True
            self.hierarchy_list = self.inputs["Setup"]["CMAESHierarchy"]
            invalid_hierarchy = [
                obs for obs in self.hierarchy_list if obs not in self.Constraints
            ]
            if len(invalid_hierarchy) > 0:
                raise NameError(
                    "Hierarchy observables must be constrained observables. "
                    f"Invalid entries: {invalid_hierarchy}"
                )
            self.hierarchy_indices = sorted(
                [self.Constraints.index(obs) for obs in self.hierarchy_list]
            )
        else:
            self.use_hierarchy = False

        self.resultsdir = os.path.join(
            self.inputs["Setup"]["cwd"], self.inputs["Setup"]["RunName"], "Results"
        )
        self.results_csv = os.path.join(self.resultsdir, "Results.csv")
        self.resume = False

        if os.path.exists(self.resultsdir) and not self.resume:
            shutil.rmtree(self.resultsdir)

        if not os.path.exists(self.resultsdir):
            os.makedirs(self.resultsdir)
        # Create the likelihood functions based on 'scaling' property of variables, if they have it; otherwise use gaussian likelihood
        # self.get_likelihood,self.observable_masks=MakeGlobalLikelihood(self.inputs['Observables'],self.use_loglike)

        if "Initial Sample" in self.inputs["Setup"]:
            from bsmart.ml.seed_points import get_previous_sample

            seed_inputs, _ = get_previous_sample(
                sample_file=self.inputs["Setup"]["Initial Sample"],
                inputs=self.inputs,
                nsamples=self.population_size,
                returntype="NLL",
            )
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
            # self.seed_inputs = np.array([[f(x) for x,f in zip(point,self.downscalers)] for point in seed_inputs])

            # NB want to set the mean to be the best point, which is the first one in the list

            self.mean = np.array(
                [f(x) for x, f in zip(seed_inputs[0], self.downscalers)]
            )  # self.seed_inputs[0]
            # self.sigma=0.01 # less variation
            # self.seed_inputs = None
            self.seed_inputs = np.array(
                [
                    [f(x) for x, f in zip(point, self.downscalers)]
                    for point in seed_inputs
                ]
            )
        else:
            self.seed_inputs = None
        # self.all_valid=[]
        # self.all_invalid=[]

    def postprocess(self, Point, observables, data_point, temp_dir, log, lock=None):

        # Codes have run successfully; return 1. The actual loss is computed elsewhere, taking into account

        return 1.0  # -math.log(self.get_likelihood(observables))

    def get_losses(self, observables):
        """Returns a list of losses.

        The C-function loss should be zero if the observable is within allowed bounds, and greater than zero outside it.
                If we include a validity condition, then other observables, evaluated afterwards, will be returned as NaN.
                        -> We should assign these the maximum loss ~ 710.

                This is then compatible with the explicit use of 'hierarchical observables', where all other observables are set to the maximum.
        """

        losses = []
        for obs_name, likelihood in zip(self.Constraints, self.constraint_likelihoods):
            obs_idx = self.observable_indices[obs_name]
            val = safe_float(observables[obs_idx])
            if math.isnan(val):
                losses.append(self.maxloss)
            else:
                losses.append((-1.0) * likelihood(val))
        return losses

        """ 
                Equivalent to 
                losses=[]
		for v, mask in zip(observables, self.observable_masks):
    		        if mask:
				f=next(likeit)
				val = safe_float(v)
				if not math.isnan(val):
                                        losses.append(val)
                                else:
                                        losses.append(self.maxloss)
                return losses
						

                """

    def run(self):
        """
        Algorithm:


        Our initial samples are in [0,1] since we use scalers.
        Subsequent rounds use the normalizing flow to generate samples
        """

        """
                def C_function(obs, obs_lower, obs_upper):
                        C_max = np.log(1 + np.finfo(np.float64).max) + 1
                        C_value = C_max
                        if isinstance(obs, str):
                                try:
                                        obs = float(obs)
                                except:
                                        raise Exception(f"Observable value {obs} not a string or float!")
                        if isinstance(obs, float):
                                C_value = np.log10(max(0.0, -obs + obs_lower, obs - obs_upper) + 1)
                        return C_value  
                """

        self.log.info("Starting CMA-ES scan!")

        # print(f"Do we store outputs? {self.RunManager.settings.store_outputs}")
        # CMA-ES bounded between [0,1]
        bounds = np.zeros((self.n_variables, 2))
        bounds[:, 1] = np.ones(self.n_variables)

        """
                if self.seed_inputs is not None: # warm start
                        cube_batch = self.seed_inputs
                        batch = [[f(y) for y,f in zip(x,self.upscalers)] for x in self.seed_inputs]
                        results=self.RunManager.run_batch(batch)
                        #print(self.RunManager.all_batch_points)
                        #full_losses = [(x,value) if isinstance(value, (int, float)) else (x,1e100) for x,value in zip(cube_batch,results)]
                        batch_observables=[point['observables'] for point in self.RunManager.all_batch_points]
                        keep_observables=[[obs for obs,mask in zip(point_obs,self.observable_masks) if mask] for point_obs in batch_observables]
                        #print(f"Got observables: {batch_observables}")
                        losses_raw = np.array([self.get_losses(obs) for obs in batch_observables])
                        
                else:
                        initial_cov = None
                """

        initial_cov = None
        optimizer = CMA(
            bounds=bounds,
            mean=self.mean,
            sigma=self.sigma,
            population_size=self.population_size,
            cov=initial_cov,
        )

        gen_idx = 0

        no_good_points_counter = 0
        no_good_points_counter_patience = int(
            self.inputs["Setup"].get("CMAESNoGoodPointsPatience", 5000)
        )
        n_good_points_early_stop = int(
            self.inputs["Setup"].get("CMAESGoodPointsEarlyStop", 5000)
        )

        best_loss = np.inf
        best_loss_counter = 0
        best_loss_patience = int(self.inputs["Setup"].get("CMAESBestLossPatience", 200))

        best_n_valid_constraints = 0.0
        best_n_valid_constraints_counter = 0
        best_n_valid_constraints_patience = int(
            self.inputs["Setup"].get("CMAESBestNValidConstraintsPatience", 2000)
        )

        good_points_counter = 0

        all_good_points = pd.DataFrame()
        for generation in range(self.max_n_generations):

            batch = []
            cube_batch = []

            self.RunManager.all_points = []  # don't need to store these overall
            self.RunManager.valid_points = []

            if generation == 0 and self.seed_inputs is not None:
                cube_batch = self.seed_inputs
                batch = [[f(y) for y, f in zip(x, self.upscalers)] for x in cube_batch]
            else:
                for _ in range(optimizer.population_size):
                    x = optimizer.ask()
                    cube_batch.append(x)
                    # batch.append([f(1.0/(1+math.exp(-y))) for y,f in zip(x,self.upscalers)])
                    batch.append([f(y) for y, f in zip(x, self.upscalers)])

            results = self.RunManager.run_batch(batch)
            # print(self.RunManager.all_batch_points)
            # full_losses = [(x,value) if isinstance(value, (int, float)) else (x,1e100) for x,value in zip(cube_batch,results)]
            batch_observables = [
                point["observables"] for point in self.RunManager.all_batch_points
            ]

            # print(f"Got observables: {batch_observables}")
            losses_raw = np.array([self.get_losses(obs) for obs in batch_observables])

            """
                        #Ordering of points in all_batch_points was not the same as for the inputs!!
                        for point,inpoint in zip(self.RunManager.all_batch_points,batch):
                                print(f" {list(map(float,point['inputs']))},{list(map(float,inpoint))}")
                        """
            # losses_raw = np.array(self.get_losses(point['observables']) for point in self.RunManager.all_batch_points)
            """
                        losses_raw = np.array([[C_function(obs_val, self.Observables[obs]["MEAN"] - self.Observables[obs]["VARIANCE"], self.Observables[obs] ["MEAN"] + self.Observables[obs]["VARIANCE"]) for obs_val, obs in zip(res, self.Observables)] for res in results])

                        """
            # self.log.info(f"Losses shape: {losses_raw.shape}")
            # print(f"Observables: {self.Observables}")
            # self.log.info(f"relevant obs: {batch_observables}")
            # self.log.info(f"Losses: {losses_raw}")

            """
                        
                        The C-function loss should be zero if the observable is within allowed bounds, and greater than zero outside it.
                        If we include a validity condition, then other observables, evaluated afterwards, will be returned as NaN. 
                                -> We should assign these the maximum loss ~ 710. 
                        
                        This is then compatible with the explicit use of 'hierarchical observables', where all other observables are set to the maximum.
                         
                        In other words, if we use validity conditions then we probably don't need to define hierarchical observables at all.

                        """

            if self.use_hierarchy:  # focus first on "hierarchical observables"
                h_losses = np.ones_like(losses_raw) * (
                    np.log(1 + np.finfo(np.float64).max) + 1
                )
                for idx, loss in enumerate(losses_raw):
                    print(idx, loss)
                    if np.all(loss[self.hierarchy_indices] == 0):
                        h_losses[idx] = losses_raw[idx]
                    else:
                        h_losses[idx, self.hierarchy_indices] = loss[
                            self.hierarchy_indices
                        ]
                losses = h_losses.copy()
            else:
                losses = losses_raw.copy()

            # print(f"Losses: {losses}")
            fitnesses_raw = np.sum(losses, axis=1)
            losses_scaled = minmax_scale(losses, axis=0)
            fitnesses = np.sum(losses_scaled, axis=1)

            # Build arrays for DataFrame
            cube_arr = np.array(cube_batch)
            batch_arr = np.array(batch)
            obs_arr = np.array(batch_observables)
            losses_arr = np.array(losses)
            total_loss_arr = np.array(fitnesses_raw).reshape(-1, 1)

            # Ensure all arrays are 2D
            if cube_arr.ndim == 1:
                cube_arr = cube_arr.reshape(-1, 1)
            if batch_arr.ndim == 1:
                batch_arr = batch_arr.reshape(-1, 1)
            if obs_arr.ndim == 1:
                obs_arr = obs_arr.reshape(-1, 1)
            if losses_arr.ndim == 1:
                losses_arr = losses_arr.reshape(-1, 1)

            # Verify array shapes match expected dimensions
            n_points = len(batch)
            n_vars = len(self.Variables)
            n_obs = len(self.Observables)
            n_constraints = len(self.Constraints)

            assert cube_arr.shape == (
                n_points,
                n_vars,
            ), f"cube_arr shape {cube_arr.shape} != expected ({n_points}, {n_vars})"
            assert batch_arr.shape == (
                n_points,
                n_vars,
            ), f"batch_arr shape {batch_arr.shape} != expected ({n_points}, {n_vars})"
            assert obs_arr.shape == (
                n_points,
                n_obs,
            ), f"obs_arr shape {obs_arr.shape} != expected ({n_points}, {n_obs})"
            assert losses_arr.shape == (
                n_points,
                n_constraints,
            ), f"losses_arr shape {losses_arr.shape} != expected ({n_points}, {n_constraints})"

            # Build column names
            box_var_cols = ["box_" + v for v in self.Variables.keys()]
            var_cols = list(self.Variables.keys())
            obs_cols = list(self.Observables)
            constraint_cols = ["C" + c for c in self.Constraints]

            all_columns = (
                box_var_cols + var_cols + obs_cols + constraint_cols + ["total_loss"]
            )

            all_points = pd.DataFrame(
                np.hstack([cube_arr, batch_arr, obs_arr, losses_arr, total_loss_arr]),
                columns=all_columns,
            )
            # all_points = pd.to_numeric(all_points,errors="coerce")

            good_points = all_points.query("total_loss == 0")

            # print(all_points['fitness'].dtype)
            # print(f"These are good: {all_points[all_points['fitness']< 1.0]}")
            penalties = np.ones_like(fitnesses)
            if len(good_points) > 0:
                self.log.info(f"Found good points in this batch:\n{good_points}")
            if (
                (len(all_good_points) > 1)
                and (len(good_points) > 1)
                and (self.novelty_detection)
            ):
                self.log.info("///////////////////////////////////////////////////////")
                self.log.info(all_good_points)
                self.log.info(self.nd_estimator)
                losses_scaled = minmax_scale(losses, axis=0)
                fitnesses = np.sum(losses_scaled, axis=1)
                # For variables, use box space (normalized [0,1]); for observables, use actual values
                box_var_columns = ["box_" + v for v in self.focus.get("Variables", [])]
                obs_columns = self.focus.get("Observables", [])
                penalty_columns = box_var_columns + obs_columns
                penalty_estimator = self.nd_estimator.fit(
                    X=all_good_points[penalty_columns].values.reshape(
                        -1, len(penalty_columns)
                    )
                )
                self.log.info(good_points)
                penalties = np.zeros_like(fitnesses)

                penalties = penalty_estimator.get_penalties(
                    all_points[penalty_columns].values.reshape(-1, len(penalty_columns))
                )
                ##

            # else:
            #     # Don't do novelty detection yet since we are still trying to find a minimum
            #     self.log.info(
            #         f"{len(good_points)} good points found; {len(all_good_points)} found previously"
            #     )
            #     penalties = np.ones_like(fitnesses)
            #     # self.log.info(f"\n{all_points}")

            # Mask penalties: 1 for invalid points, ND penalty for valid points
            penalties = np.where(np.all(np.array(losses) == 0, axis=1), penalties, 1.0)
            fitnesses = fitnesses + penalties

            # fitnesses = np.where(
            #     np.all(np.array(losses) == 0, axis=1), penalties, fitnesses + 1
            # )

            # all_points["total_loss"] = np.sum(losses, axis=1)
            all_points["ratio_n_valid_constraints"] = np.where(
                losses == 0.0, 1, 0
            ).mean(axis=1)
            all_points["good_point"] = (all_points["total_loss"] == 0.0).astype(int)

            all_points["fitness"] = fitnesses - penalties
            all_points["penalty"] = penalties

            # self.log.info(f"\n{all_points}")
            # Select only the columns we want to display
            display_cols = (
                self.Constraints  # constrained observables
                + constraint_cols  # their losses
                + [
                    "total_loss",
                    "ratio_n_valid_constraints",
                    "good_point",
                    "fitness",
                    "penalty",
                ]
            )
            self.log.info(f"\n{all_points[display_cols]}")

            if not os.path.exists(self.results_csv):

                all_points.to_csv(self.results_csv)
            else:
                all_points.to_csv(self.results_csv, mode="a", header=False)
            # solutions = [(pop,value) for pop,value in zip(cube_batch, fitnesses_raw)]
            solutions = [(pop, value) for pop, value in zip(cube_batch, fitnesses)]
            # all_good_points = pd.concat([all_good_points, good_points])
            optimizer.tell(solutions)

            gen_best_loss = losses.sum(axis=1).min()
            if gen_best_loss < best_loss:
                best_loss = gen_best_loss
                best_loss_counter = 0
            else:
                best_loss_counter += 1

            gen_best_n_valid_constraints = all_points["ratio_n_valid_constraints"].max()
            if gen_best_n_valid_constraints > best_n_valid_constraints:
                best_n_valid_constraints = gen_best_n_valid_constraints
                best_n_valid_constraints_counter = 0
            else:
                best_n_valid_constraints_counter += 1

            good_points_counter += all_points["good_point"].sum()
            # good_points_found += results["good_point"].sum()

            if all_points["good_point"].sum() > 0:
                all_good_points = pd.concat(
                    [all_good_points, all_points[all_points.good_point == 1]],
                    ignore_index=True,
                )
                no_good_points_counter = 0
            else:
                no_good_points_counter += 1

            # =====================================================================================

            gen_idx += 1

            if good_points_counter >= n_good_points_early_stop:
                print(
                    f"Stopped! Max number of {n_good_points_early_stop} good points reached."
                )
                break  # Always break when target number of good results have been reached
            if gen_idx >= self.max_n_generations:
                print(f"Stopped! Max number of {gen_idx} generations reached.")
                break  # Always use N_GEN as maximum budget

            if any(
                [
                    optimizer.should_stop(),
                    all(
                        [
                            best_loss_counter >= best_loss_patience,
                            best_n_valid_constraints_counter
                            >= best_n_valid_constraints_patience,
                        ]
                    ),
                ]
            ):

                if no_good_points_counter < no_good_points_counter_patience:
                    continue  # Still finding good results
                else:
                    print(
                        f"Stopped! Patience of {no_good_points_counter_patience} generations with no good points reached"
                    )
                    break  # Else, break as we're not restarting
