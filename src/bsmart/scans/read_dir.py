"""
Scan to read input files stored in a given directory

Scan to read input files stored in a given directory:

.. code-block:: json

    "Setup": {
        "Input_Dir": "<path to directory containing input files>"

    }

"""

__meta__ = {
    "name": "read_dir",
    "requires": [],
    "settings": {
        "Input_Dir": "Path to input directory"
    }
}



from bsmart.core import Scan as Scan
import itertools
import numpy as np
import math
from bsmart import debug

import sys
import shutil
import os

class NewScan(Scan):
    """ Class to run codes over a series of points stored in a given directory. The entries corresponding to the variables are stored in the json input as variables with the name given in the first line"""
 
    def __init__(self, inputs, log):
        #print('Hello')
        Scan.__init__(self, inputs, log)
        

    def initialise(self):
        ## set the default to store the full files as output
        #print('Do initialise')
        self.runsettings.store_outputs=True
        self.runsettings.tabbed_output=True
        ## check that the input file is specified and exists
        if 'Input_Dir' not in self.inputs['Setup']:
            self.log.error('No Input_Dir (input directory) specified')
            raise SystemExit
        
        self.input_dir=self.inputs['Setup']['Input_Dir']
        self.log.debug('Looking for input files in directory: '+self.input_dir)
        #print(self.input_dir)
        if not os.path.isdir(self.input_dir):
            self.log.error('Input directory not found')
            raise SystemExit

        self.inputfiles=[[f] for f in os.listdir(self.input_dir)]
        #print(self.inputfiles)
        #self.fullinputfiles=[[os.path.join(self.input_dir,f)] for f in self.inputfiles]
        #print(self.fullinputfiles)
        self.n_input_data=len(self.inputfiles)
        #print('done initialise')

    def write_lh_file(self,point, dir, name,values_dict=None):
        """ overload Les Houches input to just copy the file to the appropriate directory """
        """ point should just be the [full filename] """
        fullpoint=os.path.join(self.input_dir,point[0])
        lhname=os.path.join(dir,name)
        #print('Copying ' +fullpoint+ ' to ' + lhname )
        shutil.copyfile(fullpoint,lhname)
        
    def run(self):
   
        self.next_batch_index=0

        if 'Batch_size' in self.inputs['Setup']:
            batch_size=eval(self.inputs['Setup']['Batch_size'])
        else:
            batch_size=self.n_input_data

            
        while self.next_batch_index < self.n_input_data:
            all_points=self.generate_parameter_points(batch_size)
            self.RunManager.run_batch(all_points)

    def generate_parameter_points(self,npoints):
        if (self.next_batch_index+npoints) > self.n_input_data:
            lastindex=self.n_input_data
        else:
            lastindex=self.next_batch_index+npoints+1
            
        all_points = self.inputfiles[self.next_batch_index:lastindex]
        self.next_batch_index = lastindex        
    
        return all_points       
