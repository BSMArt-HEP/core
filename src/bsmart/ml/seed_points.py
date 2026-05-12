"""
An algorithm for choosing a number of initial points to launch a scan, based on an already scored sample.



"""

import numpy as np
from scipy.spatial.distance import cdist
from sklearn.neighbors import NearestNeighbors
from sklearn.cluster import DBSCAN

from collections import OrderedDict
#from typing import Any, Callable, Dict, List, Literal, Optional, Text, Tuple, Union

"""

using __all__ will explicitly tell what functions to import

By convention I could also put _ in front of the 'private' functions so that they are ignored by from seed_points import *
"""

__all__ = ["select_seeds"]


def find_local_minima(X, nlls, k=10):
    nbrs = NearestNeighbors(n_neighbors=min(k+1, len(X))).fit(X)
    _, indices = nbrs.kneighbors(X)
    local_minima = []
    for i, neigh in enumerate(indices):
        neigh = neigh[neigh != i]
        if np.all(nlls[i] <= nlls[neigh]):
            local_minima.append(i)
    return np.array(local_minima, dtype=int)




def remove_outliers(nlls_all,do_cutoff=True, cutoff_multiplier=100.0):
    if do_cutoff:
        min_val = np.nanmin(nlls_all)
        cutoff_value = cutoff_multiplier * min_val
        mask = np.isfinite(nlls_all) & (nlls_all <= cutoff_value)
    else:
        mask = np.isfinite(nlls_all)

    if not np.any(mask):
        raise ValueError("All points filtered out by cutoff.")
    
    return mask

def select_seeds(coords: np.ndarray,
                    nlls: np.ndarray,
                    nseeds: int =10,
                    CUTOFF_MULTIPLIER: float = 100.0,
                    SCALE_METHOD: str = "minmax",
                    EPSILON: float = 0.1,
                    MIN_SAMPLES: int = 10,
                    ALPHA: float = 0.7,
                    ) -> np.ndarray:

    """
      Returns a list of the original ids of the points
      
        If they come from a dataframe they can be read off using df.iloc[]
        If they come from a numpy array then it's just coords[] and nlls[] that are wanted. 

    """
    try:
        assert len(coords) == len(nlls)
    except:
        raise ValueError("Length of coordinate and nll arrays not equal in select_seeds")
    
    # create a numpy array with a list of ids 
    original_ids = np.arange(len(nlls))

    # First apply basic mask on outliers
    print(f"coords: {coords}")
    print(f"nlls: {nlls}")
    first_mask=remove_outliers(nlls,True,CUTOFF_MULTIPLIER)
    first_ids=original_ids[first_mask]
    first_nlls = nlls[first_mask]
    first_coords = coords[first_mask]


    # Now we rescale

    # --- Step 2: Scale coordinates ---
    if SCALE_METHOD == "zscore":
        mu = first_coords.mean(axis=0)
        sigma = first_coords.std(axis=0)
        sigma[sigma == 0] = 1.0
        X = (coords - mu) / sigma
    elif SCALE_METHOD == "minmax":
        mn = first_coords.min(axis=0)
        mx = first_coords.max(axis=0)
        rng = mx - mn
        rng[rng == 0] = 1.0
        X = (first_coords - mn) / rng
    else:
        X = first_coords.copy()

    # Next look for 'local minima' 
    candidate_idx = find_local_minima(X, first_nlls, k=nseeds)
    #print(f"Initial points: {candidate_idx}")
    # --- Step 4: Merge minima using DBSCAN ---
    if len(candidate_idx) > 0:
        coords_minima = X[candidate_idx]
        clustering = DBSCAN(eps=EPSILON, min_samples=MIN_SAMPLES).fit(coords_minima)
        labels = clustering.labels_
        unique_labels = set(labels)
        #print(f"Unique labels: {unique_labels}")
        chosen_idx_filtered = []
        for lbl in unique_labels:
            cluster_points = candidate_idx[labels == lbl]
            best_point = cluster_points[np.argmin(first_nlls[cluster_points])]
            chosen_idx_filtered.append(best_point)

        chosen_idx_filtered = np.array(chosen_idx_filtered, dtype=int)
    #print(chosen_idx_filtered)

    #print('==================')
    # --- Step 5: If we have more or less than the target number of seeds, need to adjust
    if len(chosen_idx_filtered) != nseeds:
        nll_norm = (first_nlls - first_nlls.min()) / (first_nlls.max() - first_nlls.min() + 1e-12)

        
        if  len(chosen_idx_filtered) > nseeds: # too many! Just select the 'best' ones from the minima we have
            # start with the best one
            selected = [chosen_idx_filtered[np.argmin(first_nlls[chosen_idx_filtered])]]
            all_idx=chosen_idx_filtered
        else: # too few! Supplement with 'best' points
            all_idx = np.arange(len(X))  # nb these are not the original ids!
            selected = list(chosen_idx_filtered) # take the originals

        while len(selected) < min(nseeds, len(X)):
            #all_idx = np.arange(len(X))
            remaining = np.setdiff1d(all_idx, np.array(selected), assume_unique=False)

            dists = cdist(X[remaining], X[selected])
            min_dists = dists.min(axis=1) if len(selected) > 0 else np.ones(len(remaining))

            dist_norm = (min_dists - min_dists.min()) / (min_dists.max() - min_dists.min() + 1e-12)
            scores = (1 - ALPHA) * dist_norm - ALPHA * nll_norm[remaining]

            next_idx = remaining[np.argmax(scores)]
            selected.append(next_idx)

        chosen_idx_filtered = np.array(selected, dtype=int)
    
    # map to original set
    # Now return the indices of the chosen seeds

    selected_ids = first_ids[chosen_idx_filtered]

    return selected_ids






def get_previous_sample(sample_file: str,
                        inputs: OrderedDict,
                        nsamples: int = None,
                        returntype: str = "LogLike",
                        ): 
    """
        Routine to load data from a previous sample and return a couple of numpy arrays of the values and likelihoods for 
    """
    print('Attempting to load prior sample '+sample_file)
    try:

        import pandas as pd
        df = pd.read_csv(sample_file)
        columns=df.columns.tolist()
        sample_inputs= df[inputs['Variables'].keys()].to_numpy()
        if  'Sample LogLike Name' in inputs["Setup"]:
                    negloglikes =-df[inputs["Setup"]["Sample LogLike Name"]].to_numpy()
        elif  'Sample NLL Name' in inputs["Setup"]:
                    negloglikes = df[inputs["Setup"]["Sample NLL Name"]].to_numpy()
        elif  'Sample Likelihood Name' in inputs["Setup"]:
                    negloglikes = -np.log(df[inputs["Setup"]["Sample Likelihood Name"]].to_numpy())
        elif  'Likelihood' in columns:
                    negloglikes = -np.log(df['Likelihood'].to_numpy())
        elif  'Result' in columns:
                    if eval(inputs["Setup"].get("LogLike","True")):
                        negloglikes = -df['Result'].to_numpy()
                    else:
                        negloglikes = -np.log(df['Likelihood'].to_numpy())
        else:
                    raise ValueError("Missing required columns for likelihood calculation")
        print('Loaded Initial sample '+sample_file)


        if returntype == "LogLike":    
            prob_outs = -negloglikes
        elif returntype == "NLL":
            prob_outs = negloglikes
        elif returntype == "Likelihood":
            prob_outs = np.exp(-negloglikes)

        if nsamples is None:
            return sample_inputs,prob_outs

        if "Sample Best Only" in inputs["Setup"] and eval(inputs["Setup"]["Sample Best Only"]): # just sort and keep best
            idx = np.argpartition(negloglikes, nsamples)[:nsamples]
            # sort those n indices by actual fitness value
            idx = idx[np.argsort(negloglikes[idx])]

            return sample_inputs[idx],prob_outs[idx]


        masks =  select_seeds(sample_inputs,negloglikes,nsamples)
        
        return sample_inputs[masks],prob_outs[masks]
       
    except Exception as e:
            raise ValueError('Failed to load initial sample from file: '+str(e))
            










"""





# --- Parameters ---
#CSV_FILE = "samples.csv"  # change to your filename
CSV_FILE = "RCMAES/Spectrum_Files/AU1xSM1_output.csv"  # change to your filename
N_POINTS = 10
ALPHA = 0.7   # trade-off for fallback: closer to 1 favors low NLL, closer to 0 favors diversity

# Scaling (distance space)
SCALE_METHOD = "minmax"   # "zscore" or "minmax" or "none"

# NLL handling
do_cutoff = True
cutoff_multiplier = 100.0   # keep points with NLL <= cutoff_multiplier * min(NLL)

# DBSCAN clustering for merging minima
EPSILON = 0.1  # distance threshold in scaled coordinates
MIN_SAMPLES = 1

# --- Load data ---
df = pd.read_csv(CSV_FILE)

variables= { "ysv" : {"RANGE": [ -0.7,0.7] , "VARIANCE": 0.1},"yse" : {"RANGE": [ -0.7,0.7] , "VARIANCE": 0.1}, "ysl" : {"RANGE": [ -0.7,0.7] , "VARIANCE": 0.1}, "lamHS" : {"RANGE": [ -1.0,1.0] , "VARIANCE": 0.1},"lamH" : {"RANGE": [ 0.01,0.4] , "VARIANCE": 0.1} ,"lamS" : {"RANGE": [ 0.0,0.7] , "VARIANCE": 0.1}, "vS" : {"RANGE": [ 1e4,5.0e5] , "VARIANCE": 1e4}}

feature_cols=variables.keys()
nll_col= 'Result'
#feature_cols = df.columns[:-1]
#nll_col = df.columns[-1]

# Raw arrays
coords_all = df[feature_cols].values
nlls_all = df[nll_col].values
orig_idx_all = np.arange(len(df))

# --- Step 1: Cutoff mask ---
if do_cutoff:
    min_val = np.nanmin(nlls_all)
    cutoff_value = cutoff_multiplier * min_val
    mask = np.isfinite(nlls_all) & (nlls_all <= cutoff_value)
else:
    cutoff_value = None
    min_val = np.nanmin(nlls_all)
    mask = np.isfinite(nlls_all)

if not np.any(mask):
    raise ValueError("All points filtered out by cutoff.")

# Apply mask
df_f = df.loc[mask].copy()
coords = df_f[feature_cols].values
nlls = df_f[nll_col].values
orig_idx = df_f.index.to_numpy()

# --- Step 2: Scale coordinates ---
if SCALE_METHOD == "zscore":
    mu = coords.mean(axis=0)
    sigma = coords.std(axis=0)
    sigma[sigma == 0] = 1.0
    X = (coords - mu) / sigma
elif SCALE_METHOD == "minmax":
    mn = coords.min(axis=0)
    mx = coords.max(axis=0)
    rng = mx - mn
    rng[rng == 0] = 1.0
    X = (coords - mn) / rng
else:
    X = coords.copy()

# --- Step 3: Graph-based local minima detection ---
def find_local_minima(X, nlls, k=10):
    nbrs = NearestNeighbors(n_neighbors=min(k+1, len(X))).fit(X)
    _, indices = nbrs.kneighbors(X)
    local_minima = []
    for i, neigh in enumerate(indices):
        neigh = neigh[neigh != i]
        if np.all(nlls[i] <= nlls[neigh]):
            local_minima.append(i)
    return np.array(local_minima, dtype=int)

candidate_idx = find_local_minima(X, nlls, k=10)
print(f"Initial points: {candidate_idx}")
# --- Step 4: Merge minima using DBSCAN ---
if len(candidate_idx) > 0:
    coords_minima = X[candidate_idx]
    clustering = DBSCAN(eps=EPSILON, min_samples=MIN_SAMPLES).fit(coords_minima)
    labels = clustering.labels_
    unique_labels = set(labels)
    print(f"Unique labels: {unique_labels}")
    chosen_idx_filtered = []
    for lbl in unique_labels:
        cluster_points = candidate_idx[labels == lbl]
        best_point = cluster_points[np.argmin(nlls[cluster_points])]
        chosen_idx_filtered.append(best_point)

    chosen_idx_filtered = np.array(chosen_idx_filtered, dtype=int)
print(chosen_idx_filtered)

print('==================')
# --- Step 5: Fallback if fewer than N_POINTS ---
if len(chosen_idx_filtered) < N_POINTS:
    nll_norm = (nlls - nlls.min()) / (nlls.max() - nlls.min() + 1e-12)
    selected = list(chosen_idx_filtered)

    while len(selected) < min(N_POINTS, len(X)):
        all_idx = np.arange(len(X))
        remaining = np.setdiff1d(all_idx, np.array(selected), assume_unique=False)

        dists = cdist(X[remaining], X[selected])
        min_dists = dists.min(axis=1) if len(selected) > 0 else np.ones(len(remaining))

        dist_norm = (min_dists - min_dists.min()) / (min_dists.max() - min_dists.min() + 1e-12)
        scores = (1 - ALPHA) * dist_norm - ALPHA * nll_norm[remaining]

        next_idx = remaining[np.argmax(scores)]
        selected.append(next_idx)

    chosen_idx_filtered = np.array(selected, dtype=int)

# --- Step 6: Map back to original DataFrame ---
chosen_orig_idx = orig_idx[chosen_idx_filtered]
selected_samples = df.iloc[chosen_orig_idx]

# --- Results ---
print({
    "total_points": int(len(df)),
    "kept_after_cutoff": int(len(df_f)),
    "dropped_by_cutoff": int(len(df) - len(df_f)),
    "min_nll": float(min_val),
    "cutoff_value": (None if cutoff_value is None else float(cutoff_value)),
    "selected_count": int(len(selected_samples)),
    "candidate_local_minima": int(len(candidate_idx)),
})
print("Selected samples:")
print(selected_samples)

selected_samples.to_csv("selected_minima.csv", index=False)
print("Saved to selected_minima.csv")














"""