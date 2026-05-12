"""
Code to perform very basic direct detection DM limits using interpolation based on data scraped from plots.

This should ultimately be upgraded.
But in any case it is largely unneeded because MicrOMEGAs can compute its own p-value.

"""


from scipy.interpolate import UnivariateSpline
#from scipy.interpolate import interp1d
from scipy import interpolate

import numpy as np

import os

from math import log, log10

import sys

class DMlimits: 

    def __init__(self):
        ''' Here I'm assuming that sys.path[0] is the root bsmart path '''
        self.dataSI = None
        self.dataSDp = None
        self.dataSDn = None
        
        from importlib import resources
        ref = resources.files('bsmart') / 'dmlimits'
        
        with resources.as_file(ref / 'DD_SI.dat') as path:
            self.dataSI = np.loadtxt(path)
        with resources.as_file(ref / 'DD_SD_P.dat') as path:
            self.dataSDp = np.loadtxt(path)
        with resources.as_file(ref / 'DD_SD_N.dat') as path:
            self.dataSDn = np.loadtxt(path)
        self.PLANCKOH2=0.11933
        self.PLANCKERROR=0.00091
        self.PLANCKLIMIT=0.122

        
        self.fSI = interpolate.UnivariateSpline(self.dataSI[:,0],[log10(x) for x in self.dataSI[:,1]], k=3, s=0)

        self.fSDp = interpolate.UnivariateSpline(self.dataSDp[:,0],self.dataSDp[:,1], k=3, s=0)
        self.fSDn = interpolate.UnivariateSpline(self.dataSDn[:,0],self.dataSDn[:,1], k=3, s=0)


    def DMexclusion(self,Mdm,Oh2,SI,SDp,SDn):
        if(Oh2 > self.PLANCKLIMIT):
            return 0

        ## don't interpolate beyond the range of our scraped data
        if Mdm > 10000.0 or Mdm < 2.0:  
            return 1
#        print("exclusion limit here: %.4e" % self.fSI([Mdm]))
        relativeXS=Oh2/self.PLANCKOH2
        
#        print(log10(relativeXS*SI))
        try:
            if (log10(relativeXS*SI) > self.fSI([Mdm])):
                return 0
        except:
            if relativeXS*SI > 1:
                """ can only fail if the XS is too large or too small """
                return 0
            else:
                return 1

        if (Mdm < 1001.0):
            if(relativeXS*SDp > self.fSDp([Mdm])):
                return 0

            if(relativeXS*SDn > self.fSDn([Mdm])):
                return 0

        return 1

    
    def excludefile(self,filename):
        ''' open file omg.out and read in the relevant lines, then return exclusion (or not) as 1 or 0 '''
        try:
            IF = open(filename,'r+')
        except:
            raise

        DMread={}
        for line in IF:
            li = line.strip().upper()
            if li.startswith('#') or len(li) < 1:
                continue
           
            li=li.replace('\t',' ')
            data=li.split()
            if "#" in data:
                data = data[:data.index("#")]

            if len(data)>1:
                DMread[int(data[0])]=float(data[1])

        ## Now extract Mdm,Oh2,SI,SDp,SDn
        Mdm=DMread[3]
        Oh2=DMread[1]
        SDp=DMread[202]*1e-36
        SDn=DMread[204]*1e-36
        SI=(54.0*DMread[201]+(131.0-54.0)*DMread[203])*1e-36/131.0

        excl=self.DMexclusion(Mdm,Oh2,SI,SDp,SDn)
        #IF.write('501 %d # Direct detection exclusion\n' % excl)
        IF.close()
        return excl
        
             
