"""
Random Scan


Gives random points within specified ranges for the variables:

.. code-block:: json

    "Variables": {
        "var1" :{ "RANGE": [100,1000]}

    }

etc

The only other setting is the number of points generated, specified in the  Setup block:

.. code-block:: json

    "Setup":{
        "Points": 1000
    }



"""

__meta__ = {
    "name": "Random",
    "requires": ["numpy"],
    "settings": {
        "Points": "Int"
    }
}



from bsmart.core import Scan as Scan
import itertools
import numpy as np
import math
from bsmart import debug



class NewScan(Scan):
    """Scanner class for Random Scans"""
 
    def __init__(self, inputs, log):
        Scan.__init__(self, inputs, log)
        # set number of points
        if 'Points' in self.inputs['Setup']:
            self.maxpoints=self.inputs['Setup']['Points']
        else:
            raise NameError("Random Scan needs a number of points!")
            #self.maxpoints=0
        

    #def initialise(self):
    #    ## set the default to store tabbed output
    #    self.runsettings.tabbed_output=True
        
        
    def run(self):
        all_points=self.generate_parameter_points()        
        self.RunManager.run_batch(all_points)

       
    def get_random_point(self):
            return [np.random.uniform(x['RANGE'][0], x['RANGE'][1])
                    for x in self.inputs['Variables'].values()]
 
    def generate_parameter_points(self):
        all_points = []
        for _ in range(self.maxpoints):
            all_points.append(self.get_random_point())
        return all_points       

    def postprocess(self,Point, observables, data_point,temp_dir,log, lock=None):
        """ Not needed """
        return ''
