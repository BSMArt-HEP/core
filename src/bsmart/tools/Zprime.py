"""
Tool to run very basic Z prime limits, based (for now) on Z' explorer arXiv:2109.13194


Takes as input the cross-section for the Z' prime and the branching ratios, read from the slha file


"""
__meta__ = {
    "name": "Zprime",
    "requires": ["scipy"],
    "settings": {
        "PDGs": "List of Z prime PDGs (optional: default is 31)"
    }
}


import os

import shutil

from bsmart import debug
from bsmart.HEPRun import HepTool,DataPoint


from bsmart import zslha
#import fnmatch
import sys
from importlib import resources
#import math
from scipy.interpolate import UnivariateSpline
#from scipy.interpolate import interp1d
from scipy import interpolate

import numpy as np

class XSBR():
    def __init__(self, xsfile):
        """
        Initialize XSBR with a file path or a file-like object.
        """
        self.xsdata = []
        self.xsint = None
        self.minmass = 0.0
        self.maxmass = 0.0
        self.loadxs(xsfile)

    def loadxs(self, xsfile):
        # Determine if xsfile is a path or a file-like object
        close_file = False
        if hasattr(xsfile, 'read'):
            IF = xsfile
        else:
            try:
                IF = open(xsfile, 'r')
                close_file = True
            except:
                raise NameError("Couldn't open XS file " + str(xsfile))

        xsdata = []
        ncol = 2
        cols_to_take = [0, 1]

        try:
            for line in IF:
                if line.startswith('#'):
                    continue
                sline = line.strip()

                data = sline.split()
                if len(data) < ncol:  # try whitespace delimiter instead
                    data = sline.split(',')
                    if len(data) < ncol:
                        continue
                if len(xsdata) > 1 and float(data[0]) <= xsdata[-1][0]:
                    continue

                xsdata.append([float(data[x]) for x in cols_to_take])
        finally:
            if close_file:
                IF.close()

        self.xsdata = np.array(xsdata)
        
        # Only create spline if data was successfully loaded
        if len(self.xsdata) > 0:
            self.xsint = interpolate.UnivariateSpline(self.xsdata[:, 0], self.xsdata[:, 1], k=3, s=0)
            self.minmass = min(self.xsdata[:, 0])
            self.maxmass = max(self.xsdata[:, 0])

    def getxs(self, mass):
        if mass < self.minmass:
            return 1e10
        if mass > self.maxmass:
            return 1e10
        res = self.xsint(mass)
        return res

 
class NewTool(HepTool):
    """ Runs madgraph and extracts the cross-section plus uncertainty """
    def __init__(self, name, settings,global_settings=None):
        HepTool.__init__(self, name, settings,global_settings)

        #if 'Data Path' not in settings or not os.path.exists(settings['Data Path']):
        #    raise NameError('Require "Data Path" in settings for Z prime limits')
        
        data_root = resources.files('bsmart') / 'data' / 'ZprimeExplorer'
        
        def load_limit(filename):
            try:
                # Open resource wrapper that works for both files and zip-compressed packages
                with (data_root / filename).open('r') as f:
                    return XSBR(f)
            except Exception as e:
                # Fallback or detailed error message
                raise NameError(f"Failed to load data file '{filename}' from {data_root}: {e}")

        self.limjj = load_limit('jj.dat')
        self.limbb = load_limit('bb.dat')
        self.limtt = load_limit('tt.dat')
        self.limee = load_limit('ee.dat')
        self.limmumu = load_limit('mumu.dat')
        self.limtautau = load_limit('tautau.dat')
        self.limWW = load_limit('WW.dat')
        self.limZH = load_limit('ZH.dat')


        if 'PDGs' not in settings:
            self.PDGs=[31]
        else:
            self.PDGs=settings['PDGs']



    def run(self, spc_file, temp_dir, log,data_point):
        

        if data_point is None or data_point.spc is None:
            spc=zslha.read(spc_file)
        else:
            spc=data_point.spc


        bestlim=0.0
        bestlabel=''
        bestpdg=0
        limlabels=['jj','bb','tt','ee','mumu','tautau','WW','ZH']

        for pdg in self.PDGs:
            try:
                pointXS=float(spc.Value('XSECTION',[1.3e4,2212,2212,pdg]))
            except:
                log.info('Could not find 13 TeV cross-section for Z prime %d, assuming that it is out of range' %pdg)
                continue
                #raise NameError('Could not find 13 TeV cross-section for Z prime %d' %pdg)
            
            try:
                zpmass=float(spc.Value('MASS',[pdg]))*1e-3  ### limits in TeV
            except:
                #self.log.info('Could not find mass Z prime %d,' %pdg)
                raise NameError('Could not find mass for Z prime %d' %pdg)

          
            ZHBR=float(spc.safeValue('BR',[pdg,[25,23]],0.0))
            WWBR=float(spc.safeValue('BR',[pdg,[24,-24]],0.0))
            visbr=sum([float(spc.safeValue('BR',[pdg,[x,-x]],0.0)) for x in range(1,7)])

            jjbr=sum([float(spc.safeValue('BR',[pdg,[x,-x]],0.0)) for x in range(1,5)])

            lims=[]

            

            lims.append(pointXS*jjbr/self.limjj.getxs(zpmass))
            lims.append(pointXS*float(spc.safeValue('BR',[pdg,[5,-5]],0.0))/self.limbb.getxs(zpmass))
            lims.append(pointXS*float(spc.safeValue('BR',[pdg,[6,-6]],0.0))/self.limtt.getxs(zpmass))
            lims.append(pointXS*float(spc.safeValue('BR',[pdg,[11,-11]],0.0))/self.limee.getxs(zpmass))
            lims.append(pointXS*float(spc.safeValue('BR',[pdg,[13,-13]],0.0))/self.limmumu.getxs(zpmass))
            lims.append(pointXS*float(spc.safeValue('BR',[pdg,[15,-15]],0.0))/self.limtautau.getxs(zpmass))
            lims.append(pointXS*float(spc.safeValue('BR',[pdg,[24,-24]],0.0))/self.limWW.getxs(zpmass))
            lims.append(pointXS*float(spc.safeValue('BR',[pdg,[25,23]],0.0))/self.limZH.getxs(zpmass))

            maxlim=max(lims)
            maxlabel=limlabels[lims.index(maxlim)]

            if maxlim > bestlim:
                bestlim=maxlim
                bestlabel=maxlabel
                bestpdg=pdg
        
        OF = open(spc_file,"a")

        OF.write('BLOCK ZPRIME #\n')
        OF.write(' 1 %10.6E # Best Limit ratio (>1 excluded), from channel %s\n' % (bestlim,bestlabel))
        OF.write(' 2 %d # pdg for best limit \n' % bestpdg)
        
        OF.close()



    
        if data_point is not None:
            if data_point.spc is None:
                data_point.spc = spc
        
            data_point.spc.blocks['ZPRIME']={'1': bestlim, '2': bestpdg}
            data_point.spc.blockcomments['ZPRIME'] = {'1': 'Best Limit ratio (>1 excluded), from channel '+bestlabel, '2': 'pdg for best limit'}
