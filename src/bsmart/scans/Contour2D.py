"""
Basic scan to find contours in two variables


You can specify the contour either via the option 'Function' (as a function of variables and observables, e.g. "mh" or "mh-125.0"). If not it will create
likelihoods based on the observables. In both cases you need to specify a 'Threshold' for the contour value. 

Note this comes along with an upgrade to the BSMplots to allow contour plots in 'csv' mode. You just add 'Contours':[125.0] etc to the options in the plot and make sure 
three variables are specified.


"""

__meta__ = {
    "name": "Contour2D",
    "requires": ["scikit-image", "scipy", "numpy", "pandas"],
    "settings": {
        "Bad value": "Float",
        "N0": "Initial grid points",
        "Threshold": "Float",
        "Function": "Contour function",
        "Input_CSV_File": "Path"
    }
}



import sys
from bsmart.core import Scan as Scan
import itertools
from bsmart import debug
from collections import OrderedDict
import numpy as np
import math
from bsmart.HEPRun import DataPoint
import random as random

import bsmart.BSMutil as BSMutil

from bsmart.BSMlikelihood import MakeLikelihoods

import time
import os
from skimage import measure
from scipy.interpolate import LinearNDInterpolator





def distance(pt1,pt2):
    return (pt1[0]-pt2[0])**2+(pt1[1]-pt2[1])**2

def removeclose(pts,threshold):
    final=[]
    for pt in pts:
        mindistance=4000.0
        for ptb in final:
            dist=distance(pt,ptb)
            #if dist 
            if dist < mindistance:
                mindistance = dist
            if mindistance < threshold:
                break
        if mindistance > threshold:
            final.append(pt)
    return final





class NewScan(Scan):
    """ """

    def initialise(self):
        """ This is called during init, but allows us to override certain things like the invalid_return_value """

        self.runsettings.invalid_return_value=self.bad_value

    def __init__(self, inputs, log):


        if "Bad value" not in inputs["Setup"]:
            self.bad_value=-1.0
        else:
            self.bad_value=float(inputs["Setup"]["Bad value"])

        Scan.__init__(self, inputs, log)

        ## scale all inputs to be between 0 and 1
        self.downscalers,self.upscalers = BSMutil.create_scalers(self.inputs,log)

        if len(self.upscalers) != 2:
            log.error('You need to specify exactly two variables with a RANGE!')
            raise

            

        ## Will store all the points naively

        self.all_points=[]
        self.allZ=[]

        if "N0" in self.inputs["Setup"]:    
            self.N0=int(self.inputs["Setup"]["N0"])
        else:
            self.N0=5
                        
                        
        self.Ngrid=self.N0
        self.nits=8
        self.threshold_to_grid=0.0001
        self.threshold_contour=0.001

        self.x=[]
        self.y=[]
        self.X=[]
        self.Y=[]
        self.Z=[]
        
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

        self.bootstrap_file=None
        if "Input_CSV_File" in self.inputs["Setup"]:
            in_file=self.inputs["Setup"]["Input_CSV_File"]
            if os.path.exists(in_file):
                self.bootstrap_file=in_file
                        
                
                        
        

    def min_distance(self,pt):
        mindistance = 4000.0
        for ptb in self.all_points:
            dist=distance(pt,ptb)
            if dist < mindistance:
                mindistance = dist
        return mindistance

    


    def FirstGrid(self):
        if self.bootstrap_file is None:
            self.x = np.linspace(0, 1.0, self.N0+1)
            self.y = np.linspace(0, 1.0, self.N0+1)
                        
        
            firstbatch=[[xx,yy] for xx in self.x for yy in self.y]
            proposedpoints=[[sf(y) for sf,y in zip(self.upscalers,point)] for point in firstbatch]
            
            results =self.RunManager.run_batch(proposedpoints)
            Zgrid=np.array(results).reshape(len(self.x),len(self.y))
            contours = measure.find_contours(Zgrid,self.threshold)

        
            self.all_points,self.allZ = map(list,zip(*((xy, z) for xy,z in zip(firstbatch,results) if z != self.bad_value)))
            
            return contours
        else:
            import pandas as pd
            input_data=pd.read_csv(self.bootstrap_file)
            invars=[x for x in self.inputs["Variables"].keys()]
            
            initial_points=pd.DataFrame(input_data,columns=invars).to_numpy().tolist()
            initial_results=pd.DataFrame(input_data,columns=['Result']).to_numpy()[:,0].tolist()
            self.all_points=[[sf(y) for sf,y in zip(self.downscalers,point)] for point in initial_points]            
            self.allZ=initial_results

            interp=LinearNDInterpolator(self.all_points,self.allZ)
            self.Ngrid=1000
            xI = np.linspace(0,1,self.Ngrid+1)
            yI = np.linspace(0,1,self.Ngrid+1)

            XI, YI = np.meshgrid(xI,yI)
            ZI = interp(XI, YI)

            contours = measure.find_contours(ZI,self.threshold)
            return contours
            


    def postprocess(self,Point, observables, data_point,temp_dir,log, lock=None):
        """ Get the function as a product of the likelihoods """
        if self.function_mode:  # use a function of variables and observables to get the contour
            all_dict=data_point.get_all_dict()
            res=float(eval(self.contour_function,all_dict))
            return res
        else:  # use likelihood functions instead
            lh = 1.
            for  v, f in zip(observables,self.likelihoodfns):
                lh = lh * f(float(v))
            return lh #-self.threshold

    

    def run(self):
        """ """
        ## First get the initial grid
        contours=self.FirstGrid()

        
        for _ in range(self.nits):
            newpoints=[]
            for cnts in contours:
                new_points=[]
    
                new_points=[[cnt[1]/(self.Ngrid),cnt[0]/(self.Ngrid)] for cnt in cnts]
        
                newpoints.extend(new_points)

            newpoints=removeclose(newpoints,self.threshold_contour)

            tocheck=[]

            for pt in newpoints:
                if self.min_distance(pt) < self.threshold_to_grid:
                    continue
                else:
                    tocheck.append(pt)
            if len(tocheck) == 0:
                break

            batch=[[sf(y) for sf,y in zip(self.upscalers,point)] for point in tocheck]

            newZ=self.RunManager.run_batch(batch)

            goodbatch,goodZ = zip(*((xy, z) for xy,z in zip(tocheck,newZ) if z != self.bad_value))

            self.all_points.extend(goodbatch)
            self.allZ.extend(goodZ)

            interp=LinearNDInterpolator(self.all_points,self.allZ)
            self.Ngrid=1000
            xI = np.linspace(0,1,self.Ngrid+1)
            yI = np.linspace(0,1,self.Ngrid+1)

            XI, YI = np.meshgrid(xI,yI)
            ZI = interp(XI, YI)

            contours = measure.find_contours(ZI,self.threshold)    


        """
        ## Now prepare a contour plot ... ?
    
        import matplotlib.pyplot as plt
        import matplotlib as mpl

        allscaled=[[sf(y) for sf,y in zip(self.upscalers,point)] for point in self.all_points]

        interp=LinearNDInterpolator(allscaled,self.allZ)
        self.Ngrid=1000
        xI = [self.upscalers[0](x) for x in np.linspace(0,1,self.Ngrid+1)]
        yI = [self.upscalers[1](y) for y in np.linspace(0,1,self.Ngrid+1)]

        

        XI, YI = np.meshgrid(xI,yI)
        ZI = interp(XI, YI)
        fig, ax = plt.subplots()
        CS1 = ax.contour(XI, YI, ZI, levels=[self.threshold])
        plt.grid('on')

        plt.show()

        """
