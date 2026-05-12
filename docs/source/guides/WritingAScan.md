Writing your own scan
=====================

The inbuilt scans are contained within python scripts in the package in a subdirectory called `scans`. However, if you want to write your own, `BSMArt` will search the current directory for a directory called `Scans` and look for a python script with the name of the scan you have specified. For example, a scan called `AwesomeScan` would look for a file called `AwesomeScan.py` in the `Scans` directory:

```bash
Scans/
    AwesomeScan.py
```








```python

from core import Scan as Scan
import itertools
import numpy as np
import math
import debug

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
        

    def initialise(self):
        ## set the default to store tabbed output
        self.runsettings.tabbed_output=True
        
        
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

```




In addition, for documentation you can include:

- A top-level docstring; this should be of the form """
- A `__meta__` dictionary describing the scan’s metadata. This typically consists of:
    - `name`: The name of the scan.
    - `requires`: A list of required packages.
    - `Settings`: A dictionary of settings for the scan.



For information about general scan settings see {py:mod}`bsmart.core`