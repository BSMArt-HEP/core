"""
Flavio


This will attempt to compute observables using `flavio <https://flav-io.github.io/docs/>`_. 

The name of each observable should correspond to the name in flavio; see `flavio's observables <https://flav-io.github.io/docs/observables.html>`_.

This is still very basic and intended for development usage.




"""

__meta__ = {
    "name": "flavio",
    "requires": ["flavio", "wilson"],
    "settings": {
        "Observables": "List of observables",
        "InputFile": "Path to WC file"
    }
}




import flavio

from wilson import Wilson




import os
import shutil
#import subprocess
from bsmart import debug
from bsmart.HEPRun import HepTool, DataPoint

from wilson import Wilson
from flavio.classes import Observable, Prediction
from math import sqrt,pi
from bsmart import zslha

"""
Definitions of lepton EDMs which for some reason aren't included in flavio
"""


def Ledm(wc_obj, par, l, eft, scale):
    r"""EDM of lepton `l` at scale `scale`"""
    wc = wc_obj.get_wc(sector='dF=0', scale=scale, par=par, eft=eft, basis='flavio')
    GF = par['GF']
    m = par['m_' + l]
    pre= 0.1973269718e-13*0.30286190409413794393*m*GF / 2/sqrt(2) / pi**2
    C7 = wc['C7_' + 2 * l]  # e.g. C7_mumu
    return pre * C7.imag

def eEDM(wc_obj, par):
    r"""EDM of the muon """
    scale = flavio.config['renormalization scale']['e g-2']
    return Ledm(wc_obj, par, l='e', eft='WET-3', scale=scale)

def muEDM(wc_obj, par):
    r"""EDM of the muon """
    scale = flavio.config['renormalization scale']['mudecays']
    return Ledm(wc_obj, par, l='mu', eft='WET-3', scale=scale)



class NewTool(HepTool):
    """ overload the init to initialise the DM limit calculator, but name and settings are already given in HepTool """
    def __init__(self, name, settings,global_settings=None):
        HepTool.__init__(self, name, settings,global_settings)
        self.wanted_observables=[str(obs) for obs in settings['Observables']] ## works for list or OrderedDict
        #self.wanted_observables=list(settings['Observables'].keys())  ## now assume that 
        Observable('mu_EDM')
        Prediction('mu_EDM',muEDM)

        Observable('e_EDM')
        Prediction('e_EDM',eEDM)
        
    def run(self, spc_file, temp_dir, log,data_point):
        try:
            IF =open(self.settings['InputFile'],'r')
            my_w=Wilson.load_wc(IF)
            IF.close()
        except:
            raise
        

        myres=[flavio.np_prediction(myobs,my_w)-flavio.sm_prediction(myobs) for myobs in self.wanted_observables]

        try:
            OF=open(spc_file,"a")
        except:
            log.debug('Failed to open '+str(spc_file)+' for appending')
            raise
        
        if data_point.spc is None:
            data_point.spc=zslha.SLHA()
        data_point.spc.blocks['FLAVIO']={}
        data_point.spc.blockcomments['FLAVIO']={}
        OF.write('BLOCK FLAVIO #\n')
        for i,res in enumerate(myres):
            #data_point.obs_dict[self.wanted_observables[i]] = res   # for future use
            OF.write('%3d  %.8E # %s\n' %(i+1,res,self.wanted_observables[i]))
            data_point.spc.blocks['FLAVIO'][str(i+1)] = res
            data_point.spc.blockcomments['FLAVIO'][str(i+1)]=self.wanted_observables[i]
        OF.close()
