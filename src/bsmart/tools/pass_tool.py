"""
Dummy class just for generating parameter cards


"""

__meta__ = {
    "name": "pass_tool",
    "requires": [],
    "settings": {}
}


import os

from bsmart import debug
from bsmart.HEPRun import HepTool, DataPoint
from bsmart import zslha

#from time import sleep

class NewTool(HepTool):
    """ overload the init to initialise the DM limit calculator, but name and settings are already given in HepTool """
    def __init__(self, name, settings,global_settings=None):
        HepTool.__init__(self, name, settings,global_settings)

        
    def run(self, spc_file, temp_dir, log,data_point):
        #sleep(0.2)
        if not self.no_slha_output and data_point.spc is None:
            data_point.spc=zslha.read(spc_file)
        return
