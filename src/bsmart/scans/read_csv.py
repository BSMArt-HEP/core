"""
Scan that reads values for the variables from a Comma Separated Values file whose first line contains a list of the variable names (as used by pandas) 




Requires:

.. code-block:: json

    "Setup: {
        "Input_CSV_File" : "<csv filename>"
    }


"""

__meta__ = {
    "name": "Read CSV",
    "requires": ["pandas"],
    "settings": {
        "Input_CSV_File": "Path to input CSV file"
    }
}


from bsmart.core import Scan as Scan
import itertools
import numpy as np
import math
from bsmart import debug

import sys
import pandas as pd
import os

class NewScan(Scan):
    """ Class to run codes over a series of points stored in a csv file. The entries corresponding to the variables are stored in the json input as variables with the name given in the first line"""
 
    def __init__(self, inputs, log):
        Scan.__init__(self, inputs, log)
        

    def initialise(self):
        ## set the default to store the full files as output
        #self.runsettings.store_outputs=True

        ## check that the input file is specified and exists
        if 'Input_CSV_File' not in self.inputs['Setup']:
            self.log.error('No Input_CSV_File specified')
            raise SystemExit
        
        self.input_csv_file=self.inputs['Setup']['Input_CSV_File']
        if not os.path.isfile(self.input_csv_file):
            self.log.error('CSV file not found')
            raise SystemExit

        

    def run(self):
        self.input_data=pd.read_csv(self.input_csv_file)
        #self.input_variables=pd.DataFrame(self.input_data, columns=[ self.inputs['Variables'][var] for var in self.inputs['Variables']]).values.tolist()
        self.input_variables=pd.DataFrame(self.input_data, columns=[ var for var in self.inputs['Variables']]).values.tolist()
        
        self.n_input_data=self.input_data.shape[0]
        self.next_batch_index=0
        #print(self.input_variables)
        if 'Batch_size' in self.inputs['Setup']:
            batch_size=eval(self.inputs['Setup']['Batch_size'])
        else:
            batch_size=self.n_input_data
        if batch_size == 0 or self.n_input_data == 0:
            self.log.error('Zero batch size')
            raise
        
        while self.next_batch_index < self.n_input_data:
            all_points=self.generate_parameter_points(batch_size)
            
            self.RunManager.run_batch(all_points)

    def generate_parameter_points(self,npoints):
        if self.next_batch_index+npoints >= self.n_input_data:
            lastindex=self.n_input_data
        else:
            lastindex=self.next_batch_index+npoints
        #print("last index: "+str(lastindex))
        all_points = self.input_variables[self.next_batch_index:(lastindex+1)]
    
        self.next_batch_index = lastindex
        #print("next batch index: %d, n_input_data: %d" %(self.next_batch_index,self.n_input_data))
    
        return all_points


    def postprocess(self,Point, observables, data_point,temp_dir,log, lock=None):
        """ No postprocessing """
        return ''
