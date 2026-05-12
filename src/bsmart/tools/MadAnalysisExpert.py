"""
Tool to run madanalysis through the PAD. 


Need to specify in tools:

.. code-block:: json

    "MadAnalysisPADSFS" : {
        "MadGraph" : "/path/to/modeldir/bin/generate_events",
        "MA5 Path" : "/path/to/ma5",
        "MA5 Card" : "MA5 card for the PAD",
        "Process" : "process_name_without_spaces" (optional)
    }



### look in MA5 recast_configuration.py for lots of useful info

"""

__meta__ = {
    "name": "MadAnalysisExpert",
    "requires": [],
    "settings": {
        "Command": "Path to ma5",
        "Process": "Process name"
    }
}


import os
import shutil

from bsmart import debug
from bsmart.HEPRun import HepTool,DataPoint

#from collider.pythiacfg import PYTHIACFG

from bsmart.collider.MGaux import readxsdata,write_command_file,MadGraphRunner

#import yaml
from bsmart import zslha
#import fnmatch
import sys
#import math




###############################################################

 
class NewTool(HepTool):
    """ Runs madgraph and extracts the cross-section plus uncertainty """
    def __init__(self, name, settings,global_settings=None):
        HepTool.__init__(self, name, settings,global_settings)
        
        if 'Process' in self.settings:
            self.process = self.settings['Process'] ## to allow labelling the processes
        else:
            self.process='' 

        self.sample_name=self.process+'_sample'
        
        self.ma5_cmd=self.settings['Command']

        

        if hasattr(global_settings,'cores') and global_settings.cores > 1:
            raise NameError("MadGraph/MadAnalysis cannot be run in multicore mode!")
        # Settings
        """
        Set up MadGraph

        """
        try:
            self.MGrunner=MadGraphRunner(settings)
        except Exception as e:
            raise NameError('Failed to initialise MadGraph: '+str(e))



    def run(self, spc_file, temp_dir, log,data_point):
        
        MGinfo=self.MGrunner.run(spc_file,temp_dir,log,data_point)

        hepmc_events=MGinfo['hepmc']
        if hepmc_events is None:
            log.error('No hepmc events produced')
            raise

        ## set up the MA5 run directory
        tempMA5dir=os.path.join(temp_dir,'MA5')
        templaunchdir=os.path.join(tempMA5dir,'launch')
        if os.path.exists(tempMA5dir):
            shutil.rmtree(tempMA5dir)

        os.makedirs(templaunchdir)

        ## create the sample file
        
        ma5proc=os.path.join(templaunchdir,self.sample_name)
        

        with open(ma5proc,'w') as f:
            f.write(hepmc_events)
        
        
        ## Now to run
        #debug.command_line_log(self.ma5_cmd+' '+ma5proc, log)
        tpid=debug.command_line_log_PID(self.ma5_cmd+' '+ma5proc, templaunchdir, log)
        
       ## presumably do something with the results ...

        # Now do we store the cutflows and yoda files too? These can be picked up by copying everything

        ## something here for storing the events ... 
        #if self.store_events:
            

