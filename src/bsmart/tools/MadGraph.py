"""
Tool to just run MadGraph. Pretty much everything is in MGaux.py because it is 
also used by other tools ... 


"""

__meta__ = {
    "name": "MadGraph",
    "requires": [],
    "settings": {
        "Command": "Path to MadGraph"
    }
}

import os
import shutil

from bsmart import debug
from bsmart.HEPRun import HepTool,DataPoint
from bsmart import zslha

from bsmart.collider.MGaux import readxsdata,write_command_file,MadGraphRunner


class NewTool(HepTool):
    """ Runs madgraph and extracts the cross-section plus uncertainty """
    def __init__(self, name, settings,global_settings=None):
        HepTool.__init__(self, name, settings,global_settings)

        if hasattr(global_settings,'cores') and global_settings.cores > 1:
            raise NameError("MadGraph tool cannot be run in multicore mode!")
        """
        Set up MadGraph

        """

        if 'Command' in self.settings:
            self.settings['MadGraph'] = self.settings['Command']
        try:
            self.MGrunner=MadGraphRunner(self.settings)
        except Exception as e:
            raise NameError('Failed to initialise MadGraph: '+str(e))
        

    def run(self, spc_file, temp_dir, log,data_point):
        #Run MadGraph

        MGinfo=self.MGrunner.run(spc_file,temp_dir,log,data_point)