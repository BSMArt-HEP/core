"""
FlexibleSUSY


"""

__meta__ = {
    "name": "FlexibleSUSY",
    "requires": [],
    "settings": {
        "Command": "Path to FlexibleSUSY executable",
        "InputFile": "Path to input SLHA file"
    }
}

import os
import shutil
from bsmart import debug
from bsmart import zslha
from bsmart.HEPRun import HEPRun, RunSettings, CoreRunner, HepTool, DataPoint

class NewTool(HepTool):
    def __init__(self, name, settings,global_settings=None):
        HepTool.__init__(self, name, settings,global_settings)
        """ Check that we can find the command, otherwise not much point carrying on"""
        self.command=os.path.expanduser(self.settings['Command'])
        if not os.path.isfile(self.command):
            #fdir=os.path.dirname
            raise NameError("Cannot identify SPheno executable. Please give it (without additional command line parameters) under SPheno->settings->Command")
        
    def run(self, spc_file, temp_dir, log,data_point):

        ## default have a timeout of 600s.
        debug.command_line_log(self.command+' --slha-input-file='+self.settings['InputFile']+' --slha-output-file='+spc_file,log,ctimeout=600)
        #debug.command_line_log(self.settings['Command']+' '+self.settings['InputFile']+' '+spc_file,log,ctimeout=60)
        if os.path.exists(spc_file):
            data_point.spc = zslha.read(spc_file)
        
