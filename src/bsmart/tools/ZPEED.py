"""
Tool to run ZPEED from https://github.com/kahlhoefer/ZPEED


Some tweaks to bring up to date

"""
from __future__ import division

__meta__ = {
    "name": "ZPEED",
    "requires": ["numpy", "scipy"],
    "settings": {
        "Path": "Path to ZPEED installation"
    }
}




import os
import shutil

from bsmart import debug
from bsmart.HEPRun import HepTool,DataPoint


from bsmart import zslha
#import fnmatch
import sys
#import math

import importlib
from contextlib import contextmanager

import numpy as np
import scipy.optimize as optimize
from scipy.stats import norm
#import scipy.integrate as integrate

#Makes sure no square root of a negative value is taken
def sqrt_safe(x):
    return np.sqrt(max(x,0))

#For a given function chi2(mu) this function finds the value of mu that minimizes chi2, ensuring mu >= 0
def minimum(chi2f):

    chi2_min = optimize.minimize_scalar(chi2f, bracket=(0.5,1.5), method = 'brent', options = {'tol':.01, 'maxiter':100} )

    if not(chi2_min['success']):
        raise NameError('Warning: Failed to find minimal chi2')

    if chi2_min['x'] > 0:
        mu_min = chi2_min['x']
    else:
        mu_min = 0

    return mu_min

def get_likelihood(chi2, chi2_Asimov):

    muhat = minimum(chi2)
    Delta_chi2 = lambda mu: chi2(mu) - chi2(min(muhat,mu))
    CLs = lambda mu: (1 - norm.cdf(sqrt_safe(Delta_chi2(mu))))/norm.cdf(sqrt_safe(chi2_Asimov(mu)) - sqrt_safe(Delta_chi2(mu)))
    return [chi2(1),Delta_chi2(1), CLs(1)]




@contextmanager
def temporary_path(path):
    
    # Context manager to temporarily add a path to sys.path and change the working directory.
    
    if not path:
        yield
        return
        
    original_cwd = os.getcwd()
    sys.path.insert(0, path)
    try:
        os.chdir(path)
        yield
    finally:
        os.chdir(original_cwd)
        if sys.path and sys.path[0] == path:
            sys.path.pop(0)

 
class NewTool(HepTool):
    """ Runs madgraph and extracts the cross-section plus uncertainty """
    def __init__(self, name, settings,global_settings=None):
        HepTool.__init__(self, name, settings,global_settings)
        
        path = self.settings.get('Path')
        
        with temporary_path(path):
            #self.chi2_CLs = importlib.import_module('chi2_CLs') ## must handle the chi2 here because of python version issue 
            self.ATLAS_13TeV_calibration=importlib.import_module('ATLAS_13TeV_calibration')
            self.ATLAS_13TeV=importlib.import_module('ATLAS_13TeV')
            self.df=importlib.import_module('dileptons_functions')
            self.ATLAS_13TeV.exclude_negative_events=1
    

    def run(self, spc_file, temp_dir, log,data_point):
        

        if data_point is None or data_point.spc is None:
            spc=zslha.read(spc_file)
        else:
            spc=data_point.spc


        bestCLs=1.0
        totallogL=0.0
        totalDeltaLogL=0.0

        for nzprime in range(5):

            start=10*nzprime
            

            if str(1+start) not in spc.blocks['ZPEEDIN']:
                break

            if str(1+start) in spc.blocks['ZPEEDIN']:

                Zp_model={}
               
                pdg=int(spc.Value('ZPEEDIN',[1+start]))
                Zp_model['MZp'] = float(spc.Value('MASS',[pdg]))
                Zp_model['Gamma'] = float(spc.Width(pdg))
                #print('Gamma: '+str(Zp_model['Gamma']))
                Zp_model['gdv']=float(spc.Value('ZPEEDIN',[start+2]))
                Zp_model['gda']=-float(spc.Value('ZPEEDIN',[start+3]))
                Zp_model['guv']=float(spc.Value('ZPEEDIN',[start+4]))
                Zp_model['gua']=-float(spc.Value('ZPEEDIN',[start+5]))
                Zp_model['gev']=float(spc.Value('ZPEEDIN',[start+6]))
                Zp_model['gea']=-float(spc.Value('ZPEEDIN',[start+7]))
                Zp_model['gmuv']=float(spc.Value('ZPEEDIN',[start+8]))
                Zp_model['gmua']=-float(spc.Value('ZPEEDIN',[start+9]))

                #Step 2: Calculate differential cross section (including detector efficiency)

                Zp_model['glv']=Zp_model['gev']
                Zp_model['gla']=Zp_model['gea']
                ee_signal = lambda x : self.ATLAS_13TeV_calibration.xi_function(x, "ee") * self.df.dsigmadmll(x, Zp_model, "ee")
                ee_signal_with_interference = lambda x : self.ATLAS_13TeV_calibration.xi_function(x, "ee") * self.df.dsigmadmll_wint(x, Zp_model, "ee")


                mm_signal = lambda x : self.ATLAS_13TeV_calibration.xi_function(x, "mm") * self.df.dsigmadmll(x, Zp_model, "mm")	
                mm_signal_with_interference = lambda x : self.ATLAS_13TeV_calibration.xi_function(x, "mm") * self.df.dsigmadmll_wint(x, Zp_model, "mm")	

                

                #Step 3: Create likelihood functions

                #window_size=min(3.*Zp_model['Gamma'],0.01*Zp_model['MZp'])
                window_size=3.*Zp_model['Gamma']
                Mlow = Zp_model['MZp'] - window_size
                Mhigh = Zp_model['MZp'] + window_size

                sig_range = [Mlow,Mhigh]

                
                
                
                # chi2, chi2_Asimov = self.ATLAS_13TeV.calculate_chi2(ee_signal, mm_signal, signal_range=sig_range)
                #     chi2_with_interference, chi2_Asimov_with_interference = self.ATLAS_13TeV.calculate_chi2(ee_signal_with_interference, mm_signal_with_interference, signal_range=sig_range)

                

                #     # Step 4: Evaluate test statistic
                #     result = self.chi2_CLs.get_likelihood(chi2, chi2_Asimov)
                #     result_with_interference = self.chi2_CLsget_likelihood(chi2_with_interference, chi2_Asimov_with_interference)
                #     log.info("Without interference")
                # log.info("-2 log L:       ", result[0])
                # log.info("-2 Delta log L: ", result[1])
                # log.info("CLs:            ", result[2])

                # log.info("With interference")
                # log.info("-2 log L:       ", result_with_interference[0])
                # log.info("-2 Delta log L: ", result_with_interference[1])
                # log.info("CLs:            ", result_with_interference[2])


                result=[0.0,0.0,0.0]
                try:
                    
                    chi2_with_interference, chi2_Asimov_with_interference = self.ATLAS_13TeV.calculate_chi2(ee_signal_with_interference, mm_signal_with_interference, signal_range=sig_range)

                    # Step 4: Evaluate test statistic
                    #result = self.chi2_CLs.get_likelihood(chi2_with_interference, chi2_Asimov_with_interference)
                    result = get_likelihood(chi2_with_interference, chi2_Asimov_with_interference)
                    log.info("With interference")
                except Exception as ef:    ## Negative number?
                    log.info("Interference failed: "+str(ef))
                    try:
                        chi2, chi2_Asimov = self.ATLAS_13TeV.calculate_chi2(ee_signal, mm_signal, signal_range=sig_range)
                        result = get_likelihood(chi2, chi2_Asimov)
                        #result = self.chi2_CLs.get_likelihood(chi2, chi2_Asimov)
                        log.info("Without interference")
                    except Exception as eff:
                        log.error('Failed even without interference '+str(eff))
                        result[0]=0.0
                        result[1]=0.0
                        result[2]=1.0
                log.info("-2 log L:       %10.4e"% result[0])
                log.info("-2 Delta log L: %10.4e"% result[1])
                log.info("CLs:            %.5f"% result[2])

                
                
                bestCLs=min(result[2],bestCLs)
                totallogL=totallogL+result[0]
                totalDeltaLogL=totalDeltaLogL+result[1]

        OF = open(spc_file,"a")

        OF.write('BLOCK ZPEED #\n')
        OF.write(' 1 %10.6E # -2 log L\n' % (totallogL))
        OF.write(' 2 %10.6E # -2 Delta log L\n' % (totalDeltaLogL))
        OF.write(' 3 %10.6E # Best CLs\n' % (bestCLs))
        
        OF.close()



    
        if data_point is not None:
            if data_point.spc is None:
                data_point.spc = spc
        
            data_point.spc.blocks['ZPEED']={'1': totallogL, '2': totalDeltaLogL, '3': bestCLs}
            data_point.spc.blockcomments['ZPEED'] = {'1': '-2 log L', '2': '-2 Delta log L', '3': 'Best CLs'}
