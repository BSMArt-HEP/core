"""
Active Learning Scan to find a decision boundary using Gaussian Processes



Entropy gain taken from https://github.com/diana-hep/excursion
by Lukas Heinrich,Gilles Louppe, Kyle Cranmer
(unfortunately that code does not work any more)

Adapted by Mark Goodsell to remove integration over fixed grids and replace by integration using \f[p(c | z)\f]
Plus smaller choice of points to sample (both save a lot of time)
Plus some improvement of robustness

Settings
--------

.. code-block:: text

    Setup: {
        "Type": "ContourGP",
    	    "Points": default 60, int -> number of points to check
            "Initial Points": default 5, int -> number of initial random points to use for initialisation
            "Noise" : default None, float -> signal noise to include in the Gaussian Process
            "Grid Size":  default 100, int -> Number of points to include in the integration grids
            "Sample Size":  default 50, int -> Number of Samples to check for the maximum entropy gain
            "Sample Multiplier":    default 5, int -> We select the sample from this number * the sample size of random points, where the
                                    ones to be checked are those closest in value to the decision boundary
            "Function": if present, this specifies a function of variables and/or observables to solve for f(x) = threshold (see below)
                        if not present, we construct a likelihood function from the observables in the usual way. 
                        One very useful option is to set "SCALING": "USER" in the observable if it is e.g. signal strength or signal strength -1.
            "Threshold": default 0.0, float ->  value to subtract from the above function so that the contour we are looking for solves
                                                f(x) = 0
            "Fail Value": default 10.0, float -> value to assign to failed points
            "Run Name": string -> name of the run       

    }


"""

__meta__ = {
    "name": "ContourGP",
    "requires": ["scikit-learn", "numpy", "pandas"],
    "settings": {
        "Points": "Int",
        "Initial Points": "Int",
        "Noise": "Float",
        "Grid Size": "Int",
        "Sample Size": "Int",
        "Sample Multiplier": "Int",
        "Input_CSV_File": "Path",
        "Fail Value": "Float",
        "Threshold": "Float",
        "Function": "String"
    }
}


import os,sys
from bsmart.core import Scan as Scan
import itertools
import numpy as np
import math
from bsmart import debug

from bsmart import BSMutil

from bsmart.BSMlikelihood import MakeLikelihoods



import sklearn
from sklearn.gaussian_process.kernels import ConstantKernel, RBF, WhiteKernel
#import torch
import numpy as np
from bsmart.HEPRun import UpdateProgress, DataPoint



from scipy.linalg import cho_solve
from scipy.stats import norm


## Get sklearn warnings to shut up
def warn(*args, **kwargs):
    pass
import warnings
warnings.warn = warn


def getnewcandidates(GP,ndim=2,npoints=100,multiplier=6):
    """
        get a bunch of test candidates selected from those near to the boundary
    """
    allpoints=[]
    allyfuncs=[]
    for _ in range(npoints*multiplier):
        point=[np.random.uniform(0.0,1.0) for _ in range(ndim)]
        mu =GP.predict([point],return_std = False)
        allpoints.append(point)
        #allyfuncs.append(mu[0]*(1.0-mu[0]))
        allyfuncs.append(abs(mu[0]))
    allpoints=np.array(allpoints)
    allyfuncs=np.array(allyfuncs)
    return allpoints[np.argsort(allyfuncs)[:npoints]]
    #return allpoints[np.argsort(allyfuncs)[-npoints:]]

def getgrids(GP,ndim=2,npoints=50,thresholds=[-np.inf,0.0,np.inf]):
    """
        get grids via shooting algorithm according to p(c | z)
    """
    
    def getpcj(point,j):
        """
        probability of being in each category for a given point.
        Assuming GP evaluation is expensive, we could test each category for each point. 
        However, it is probably better to have separate points for each category, so I do that.
        """
        #cjs=[]
        mu, std =GP.predict([point],return_std = True)
        #print(pred)
        #
        
        #for j in range(len(thresholds) - 1):
        alpha_j = (thresholds[j] - mu[0]) / std[0]
        beta_j = (thresholds[j+1] - mu[0]) / std[0]
        pcj = norm.cdf(beta_j) - norm.cdf(alpha_j)
        return pcj
    
    grids=[]
    for region in range(len(thresholds) - 1):
        thisgrid=[]
        alldata=[]
        for _ in range(npoints):
            point=[np.random.uniform(0.0,1.0) for _ in range(ndim)]
            pcj=getpcj(point,region)
            alldata.append([point,pcj])
        # now get the maximum 
        
        maxpcj=max([x[1] for x in alldata])
        if maxpcj < 1e-3:
            maxpcj = 1.0

        for pt in alldata:
            if np.random.uniform(0.0,1.0) < pt[1]/maxpcj:
                thisgrid.append(pt)
        
        tries=0
        while len(thisgrid) < npoints and tries < 100*npoints:
            tries=tries+1
            point=[np.random.uniform(0.0,1.0) for _ in range(ndim)]
            pcj=getpcj(point,region)

            if pcj > maxpcj: # refactor all points
                newgrid=[]
                for pt in thisgrid:
                    if np.random.uniform(0.0,1.0) < maxpcj/pcj:
                        newgrid.append(pt)
                thisgrid=newgrid
                maxpcj=pcj
                thisgrid.append([point,pcj])
            else:
                if np.random.uniform(0.0,1.0) < pcj/maxpcj:
                    thisgrid.append([point,pcj])
        grids.append(np.array([x[0] for x in thisgrid])) # only keep the points, not the probabilities
    return grids


def approx_mi_vec(mu, cov, threshold_start,threshold_end):
    

    std_sx = []
    #for j in range(len(thresholds) - 1):
    if True:
        mu1 = mu[:, 0]
        std1 = cov[:, 0, 0] ** 0.5
        mu2 = mu[:, 1]
        std2 = cov[:, 1, 1] ** 0.5
        rho = cov[:, 0, 1] / (std1 * std2)
        alpha_j = (threshold_start - mu2) / std2
        beta_j = (threshold_end - mu2) / std2
        c_j = norm.cdf(beta_j) - norm.cdf(alpha_j)

        # \sigma(Y(X)|S(x')=j)
        npdfb=norm.pdf(beta_j)
        npdfa=norm.pdf(alpha_j)
        npdfb[~np.isfinite(beta_j)]=0.0
        npdfa[~np.isfinite(alpha_j)]=0.0
        beta_j[~np.isfinite(beta_j)] = 0.0
        alpha_j[~np.isfinite(alpha_j)] = 0.0
        b_phi_b = beta_j * npdfb
        a_phi_a = alpha_j * npdfa

        deltaab=npdfb-npdfa
        
        c_j[ c_j ==0.0] = 1.0 # always multiplied by 0

        mu_cond = mu1 - std1 * rho / c_j * deltaab
        var_cond = (mu1 ** 2 - 2 * mu1 * std1 * (rho / c_j * deltaab) +
                    std1 ** 2 * (1. - (rho ** 2 / c_j) * (b_phi_b - a_phi_a)) -
                    mu_cond ** 2)
        
        
        std_sx_j = var_cond ** 0.5

        #std_sx.append(std_sx_j)
        """
        b_phi_b = beta_j * norm.pdf(beta_j)
        b_phi_b[~np.isfinite(beta_j)] = 0.0
        a_phi_a = alpha_j * norm.pdf(alpha_j)
        a_phi_a[~np.isfinite(alpha_j)] = 0.0

        mu_cond = mu1 - std1 * rho / c_j * (norm.pdf(beta_j) - norm.pdf(alpha_j))
        var_cond = (mu1 ** 2 - 2 * mu1 * std1 * (rho / c_j * (norm.pdf(beta_j) - norm.pdf(alpha_j))) +
                    std1 ** 2 * (1. - (rho ** 2 / c_j) * (b_phi_b - a_phi_a)) -
                    mu_cond ** 2)
        std_sx_j = var_cond ** 0.5

        std_sx.append(std_sx_j)
        """

        # Entropy
        CONSTANT = (2 * np.e * np.pi) ** 0.5
        logconst=np.log(CONSTANT)
    
    #for j in range(len(thresholds) - 1):
    
        mu1 = mu[:, 0]
        std1 = cov[:, 0, 0] ** 0.5
        mu2 = mu[:, 1]
        std2 = cov[:, 1, 1] ** 0.5
        h = np.log(std1 * CONSTANT)
        p_j = norm(mu2, std2).cdf(threshold_end) - norm(mu2, std2).cdf(threshold_start)
        #dec = p_j * np.log(std_sx[j] * CONSTANT)
        dec = p_j * (np.log(std_sx_j)+logconst)
        h[p_j > 0.0] -= dec[p_j > 0.0]

    return h


def info_gain(x_candidate, gp, thresholds, meanXs,noise=0.0):  ## don't bother including noise in info gain
    
    info_total=0.0
    for j in range(len(thresholds) - 1):
        meanX=meanXs[j]
        n_samples = len(meanX)
        X_all = np.concatenate([np.array([x_candidate]), meanX]).reshape(1 + n_samples, -1)
        K_trans_all = gp.kernel_(X_all, gp.X_train_) ## off-diagonal -> no noise
        y_mean_all = K_trans_all.dot(gp.alpha_) + gp._y_train_mean
        v_all = cho_solve((gp.L_, True), K_trans_all.T)

        mus = np.zeros((n_samples, 2))
        mus[:, 0] = y_mean_all[0]
        mus[:, 1] = y_mean_all[1:]

        covs = np.zeros((n_samples, 2, 2))
        
        c = gp.kernel_(X_all[:1], X_all)  # 0:1 means just the 0th entry, but keep same number of indices
          #-> x_candidate

        covs[:, 0, 0] = c[0, 0]+noise
        covs[:, 1, 1] = c[0, 0]+noise  # all the diagonal entries are identical because they are zero distance from themselves!
        #covs[:, 1, 1] = c2[:]+noise   
        covs[:, 0, 1] = c[0, 1:] 
        covs[:, 1, 0] = c[0, 1:]

        x_train_len = len(gp.X_train_)
        K_trans_all_repack = np.zeros((n_samples, 2, x_train_len))
        K_trans_all_repack[:, 0, :] = K_trans_all[0, :]
        K_trans_all_repack[:, 1, :] = K_trans_all[1:]
        v_all_repack = np.zeros((n_samples, x_train_len, 2))
        v_all_repack[:, :, 0] = v_all[:, 0]
        v_all_repack[:, :, 1] = v_all[:, 1:].T
        covs -= np.einsum('...ij,...jk->...ik', K_trans_all_repack, v_all_repack)

        mi = approx_mi_vec(mus, covs, thresholds[j],thresholds[j+1])
        mi[~np.isfinite(mi)] = 0.0
        info_total-=np.mean(mi)
    return info_total
    #return -np.mean(mi)

class NewScan(Scan):
    """Scanner class for Grid Scans"""
 
    def __init__(self, inputs, log):
        Scan.__init__(self, inputs, log)
        self.downscalers,self.upscalers = BSMutil.create_scalers(self.inputs,log)
       
        self.ndim=len(self.upscalers)
        if 'Points' in self.inputs['Setup']:
            self.maxpoints=self.inputs['Setup']['Points']
        else:
            self.maxpoints=100

        
        self.noise=self.inputs['Setup'].get('Noise',None)
        self.initialpoints=int(self.inputs['Setup'].get('Initial Points',5))
        self.gridsize=int(self.inputs['Setup'].get('Grid Size',100))
        self.samplestocheck=int(self.inputs['Setup'].get('Sample Size',50))
        self.sample_multiplier=int(self.inputs['Setup'].get('Sample Multiplier',5))

        self.bootstrap_file=None
        if "Input_CSV_File" in self.inputs["Setup"]:
            in_file=self.inputs["Setup"]["Input_CSV_File"]
            if os.path.exists(in_file):
                self.bootstrap_file=in_file

        self.fail_value=float(self.inputs['Setup'].get('Fail Value',10.0))
        if len(self.upscalers) == 0:
            log.error('You need to specify at least one variable with a RANGE!')
            raise

        if "Threshold" not in self.inputs["Setup"]:
            self.threshold=0.0
        else:
            self.threshold=float(self.inputs["Setup"]["Threshold"])
        

        self.function_mode=False
        if 'Function' in self.inputs["Setup"]: ### will define the contour via some function
            self.function_mode=True
            self.contour_function=self.inputs["Setup"]["Function"]
        else:
            self.likelihoodfns=MakeLikelihoods(self.inputs['Observables'])


        
    def postprocess(self,Point, observables, data_point,temp_dir,log, lock=None):
        """ Get the function as a product of the likelihoods """
        if self.function_mode:  # use a function of variables and observables to get the contour
            try:
                all_dict=data_point.get_all_dict()
                res=float(eval(self.contour_function,all_dict))
                return res
            except:
                return self.fail_value
        else:  # use likelihood functions instead
            try:
                lh = 1.
                for  v, f in zip(observables,self.likelihoodfns):
                    lh = lh * f(float(v))
                return lh -self.threshold #-self.threshold
            except:
                return self.fail_value


    def get_random_point(self):
            return [np.random.uniform(0,1.0)
                    for _ in range(self.ndim)]
    
    def run(self):
        def optifunc(inpoints):
            
            points=[[self.upscalers[n](x[n]) for n in range(self.ndim)] for x in inpoints]
            
            try:
                ## Update progress bar for batches, not for individual points
                if len(points) > 1:
                    outlh=self.RunManager.run_batch(points)
                    
                else:
                    outlh=[self.RunManager.CoreRuns[0].run_point(points[0])]
                outlh = [ x  if x != [] else self.fail_value for x in outlh ]
                #print(outlh)
                res=np.array(outlh)
                return res 
            except:
                return np.array([self.fail_value]*len(points))

        mykernel=ConstantKernel() * RBF(length_scale_bounds=[0.1,0.5], length_scale=[0.1]*self.ndim)
       
        if self.noise is not None:
            self.myGP=sklearn.gaussian_process.GaussianProcessRegressor(kernel=mykernel,alpha=self.noise)
        else:
            self.myGP=sklearn.gaussian_process.GaussianProcessRegressor(kernel=mykernel)

        ## first get some random points
        self.all_points_x=[]
        self.all_points_y=[]
        start_points=[]
        if self.bootstrap_file is None:
            for _ in range(self.initialpoints):
                start_points.append(self.get_random_point())
            self.all_points_x = np.array(self.all_points_x + start_points)

            self.log.info("Running initial batch of random points")
            self.all_points_y = optifunc(start_points)
        else:       
            import pandas as pd
            input_data=pd.read_csv(self.bootstrap_file)
            invars=[x for x in self.inputs["Variables"].keys()]
            
            initial_points=pd.DataFrame(input_data,columns=invars).to_numpy().tolist()
            
            start_points=[[sf(y) for sf,y in zip(self.downscalers,point)] for point in initial_points]            
            
            
            self.all_points_x = np.array(start_points)
            self.all_points_y=pd.DataFrame(input_data,columns=['Result']).to_numpy()[:,0]


        ## gridsearch for evaluation
        import itertools
        #self.pointsgrid=np.array([np.linspace(-0.5,0.5,10)  for _ in range(self.ndim)]).transpose()
        self.pointsgrid=np.array(list(itertools.product(*[np.linspace(0.0,1.0,26) for _ in range(self.ndim)])))
        
        
        
        self.myGP=self.myGP.fit(start_points,self.all_points_y)

        self.log.info("Running main loop")

        while self.RunManager.nextpointID <= self.maxpoints:
            newpoints=getnewcandidates(self.myGP,self.ndim,self.samplestocheck,self.sample_multiplier)
        
            newgrids=getgrids(self.myGP,self.ndim,self.gridsize)
            for ng, grid in enumerate(newgrids):
                if len(grid) == 0:
                    newgrids[ng] = self.pointsgrid
            infos=np.array([info_gain(x_candidate, self.myGP, [-np.inf,0.0,np.inf], newgrids) for x_candidate in newpoints])
            
            bestpoint_n=np.argmin(infos)

            bestpoint=np.array([list(newpoints[bestpoint_n])])
            
            self.all_points_x = np.vstack((self.all_points_x,bestpoint))
            newy=optifunc(bestpoint)
            self.all_points_y = np.hstack((self.all_points_y,newy))

            self.myGP=self.myGP.fit(self.all_points_x,self.all_points_y)
            UpdateProgress(self.RunManager.nextpointID-1, self.maxpoints)



