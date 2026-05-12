"""

Routines for interpolating cross-sections

Return results in pb

Options are:

"""


import os
import numpy as np
import math
import sys

import re
from collections import OrderedDict

from scipy.interpolate import UnivariateSpline
#from scipy.interpolate import interp1d
from scipy import interpolate



class XS():
    def __init__(self):
        self.xsdata=[]
        self.xsint=None
        self.minmass=0.0
        self.maxmass=0.0
        self.mass_variable=''

        
    def loadxs(self,xsinfodict):

        xsfile=xsinfodict['File']
        try:
            IF=open(xsfile,'r')
        except:
            raise NameError("Couldn't open XS file "+xsfile)
        xsdata=[]
        uncdata=[]

        if 'Columns' in xsinfodict:
            cols_to_take=list(map(int,xsinfodict['Columns']))
        else:
            cols_to_take=[0,1,2]
        maxcol=max(cols_to_take)
        ncol=len(cols_to_take)

        if 'Mass Variable' in xsinfodict:
            self.mass_variable=xsinfodict['Mass Variable']
        else:
            self.mass_variable=''

        
        if 'Units' in xsinfodict:
            if xsinfodict['Units']=='pb':
                self.xsmultiplier=1.0
            elif xsinfodict['Units']=='fb':
                self.xsmultiplier=0.001
            
        else:  ## Assume input in fb by default
            self.xsmultiplier=0.001

        for line in IF:
            if line.startswith('#'):
                continue
            sline=line.strip()

            
            data=sline.split(',')
            if len(data) < ncol: # try whitespace delimiter instead
                data=sline.split()
                if len(data) < ncol:
                    continue

        
            xsdata.append([float(data[x]) for x in cols_to_take])

        IF.close()
        self.xsdata=np.array(xsdata)
        self.xsint=interpolate.UnivariateSpline(self.xsdata[:,0],self.xsdata[:,1],k=3,s=0)
        if ncol > 2:
            self.xserr=interpolate.UnivariateSpline(self.xsdata[:,0],self.xsdata[:,2],k=3,s=0)
        else:
            self.xserr=lambda x: 0
        self.minmass=min(self.xsdata[:,0])
        self.maxmass=max(self.xsdata[:,0])
        
    def __init__(self,xsinfodict):
        self.loadxs(xsinfodict)
        
        
    def getxs(self,mass):
        if mass < self.minmass:
            #print("Mass out of bounds")
            raise NameError('Mass too low')
            #return 0.0
        if mass > self.maxmass:
            #print("Mass out of bounds")
            raise NameError('Mass too high')
            #return 0.0
        res=self.xsint(mass)
        
        return res

    def getxserr(self,mass):
        if mass < self.minmass:
            raise NameError('Mass too low')
            #print("Mass out of bounds")
            #return 0.0
        if mass > self.maxmass:
            raise NameError('Mass too high')
            #print("Mass out of bounds")
            #return 0.0
        res=self.xserr(mass)
        
        return res
    
    def get_xs_info(self,var_dict):
        try:
            massval=float(eval(self.mass_variable,var_dict))
        except:
            raise NameError('Could not determine mass')
        
        try:
            xs=self.xsmultiplier*self.getxs(massval)
            xserr=self.xsmultiplier*self.getxserr(massval)
        except Exception as e:
            raise NameError('Could not extract cross-sections or uncertainties: ' +str(e))
        
        return xs,xserr
