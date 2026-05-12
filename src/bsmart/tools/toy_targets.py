"""
Tool to generate example optimisation targets


This tool defines example target functions that are commonly used to test Black Box optimisation algorithms.

The results are stored directly as observables:
* Rosenbrock
* Rastrigin
* Ackley
* Himmelblau
* Easom

So these can be obtained directly as observables using "FUNCTION": "Rosenbrock" etc; alternatively they can be obtained in the block TOYTARGETS:

Block TOYTARGETS 
  1  <>   # Rosenbrock function
  2  <>   # Rastrigin
  3  <>   # Ackley
  4  <>   # Himmelblau
  5  <>   # Easom

But be warned that no SLHA file is read or written for this tool: data is held in memory only. 

Note that, since the functions are positive-definite and supposed to be minimised, this is equivelent to a negative-log-likelihood, so they should be declared with scaling "MINUSEXPUSER" (since this would transform the value into a 'likelihood' between 0 and 1)

The following parameters can be configured in the settings for this tool:
* `Rosenbrock a`: The 'a' parameter for the Rosenbrock function (default: 1.0).
* `Rosenbrock b`: The 'b' parameter for the Rosenbrock function (default: 100.0).
* `Rastrigin A`: The 'A' parameter for the Rastrigin function (default: 10.0).
* `Ackley a`: The 'a' parameter for the Ackley function (default: 20.0).
* `Ackley b`: The 'b' parameter for the Ackley function (default: 0.2).
* `Ackley c`: The 'c' parameter for the Ackley function (default: 2*pi).
* `Himmelblau a`: The 'a' parameter for Himmelblau's function (default: 11.0).
* `Himmelblau b`: The 'b' parameter for Himmelblau's function (default: 7.0).
* `Easom center`: The [x, y] center for the Easom function (default: [pi, pi]).

Example usage in a BSMArt JSON file:
.. code-block:: json

    "Codes": {
            "toy_targets": {
                "Run": "True",
                "Rosenbrock a": 1.5,
                "Observables":
                {
                            "Rosen": {
                            "FUNCTION": "Rosenbrock",
                            "SCALING" : "MINUSEXPUSER"
                            },
                            "Rastrigin" : {
                            "SLHA":  [ "TOYTARGETS", [2]],
                            "SCALING": "MINUSEXPUSER"
                            }

                 }
            }



"""

__meta__ = {
 "name": "toy_targets",
 "requires": [],
 "settings": {
        'Rosenbrock a': "The 'a' parameter for the Rosenbrock function (default: 1.0).",
        'Rosenbrock b': "The 'b' parameter for the Rosenbrock function (default: 100.0).",
        'Rastrigin A': "The 'A' parameter for the Rastrigin function (default: 10.0).",
        'Ackley a': "The 'a' parameter for the Ackley function (default: 20.0).",
        'Ackley b': "The 'b' parameter for the Ackley function (default: 0.2).",
        'Ackley c': "The 'c' parameter for the Ackley function (default: 2*pi).",
        'Himmelblau a': "The 'a' parameter for Himmelblau's function (default: 11.0).",
        'Himmelblau b': "The 'b' parameter for Himmelblau's function (default: 7.0).",
        'Easom center': "The [x, y] center for the Easom function (default: [pi, pi])."
    }
}

import os

from bsmart import debug
from bsmart.HEPRun import HepTool, DataPoint
from bsmart import zslha


import math
class NewTool(HepTool):
    """ overload the init to initialise variables """
    def __init__(self, name, settings,global_settings=None):
        HepTool.__init__(self, name, settings,global_settings)
        self.Rosenbrock_a = self.settings.get('Rosenbrock a',1.0)
        self.Rosenbrock_b = self.settings.get('Rosenbrock b',100.0)
        self.Rastrigin_A = self.settings.get('Rastrigin A', 10.0)
        self.Ackley_a = self.settings.get('Ackley a', 20.0)
        self.Ackley_b = self.settings.get('Ackley b', 0.2)
        self.Ackley_c = self.settings.get('Ackley c', 2.0 * math.pi)
        self.Himmelblau_a = self.settings.get('Himmelblau a', 11.0)
        self.Himmelblau_b = self.settings.get('Himmelblau b', 7.0)
        self.Easom_center = self.settings.get('Easom center', [math.pi, math.pi])
        
    def run(self, spc_file, temp_dir, log,data_point):
        vars = data_point.inputs # don't care about the names!
        n = len(vars)
        x = vars

        blocks = {}
        blockcomments = {}

        if n > 0:
            # Rosenbrock function
            if n == 1:
                rosenbrock = (self.Rosenbrock_a - x[0])**2.0
            else: # n > 1
                rosenbrock = sum(self.Rosenbrock_b*(x[i+1]-x[i]**2.0)**2.0 + (self.Rosenbrock_a-x[i])**2.0 for i in range(n-1))
            blocks['1'] = rosenbrock
            blockcomments['1'] = 'Rosenbrock function'
            data_point.obs_dict['Rosenbrock'] = rosenbrock

            # Rastrigin function
            rastrigin = self.Rastrigin_A*n + sum(xi**2 - self.Rastrigin_A * math.cos(2 * math.pi * xi) for xi in x)
            blocks['2'] = rastrigin
            blockcomments['2'] = 'Rastrigin function'
            data_point.obs_dict['Rastrigin'] = rastrigin

            # Ackley function
            sum_sq_term = 1.0/n * sum(xi**2 for xi in x)
            sum_cos_term = 1.0/n * sum(math.cos(self.Ackley_c*xi) for xi in x)
            ackley = -self.Ackley_a * math.exp(-self.Ackley_b * math.sqrt(sum_sq_term)) - math.exp(sum_cos_term) + self.Ackley_a + math.exp(1)
            blocks['3'] = ackley
            blockcomments['3'] = 'Ackley function'
            data_point.obs_dict['Ackley'] = ackley

            if n > 1:
                # Himmelblau's function
                himmelblau = (x[0]**2 + x[1] - self.Himmelblau_a)**2 + (x[0] + x[1]**2 - self.Himmelblau_b)**2
                blocks['4'] = himmelblau
                blockcomments['4'] = "Himmelblau's function"
                data_point.obs_dict['Himmelblau'] = himmelblau

                # Easom function
                easom = -math.cos(x[0]) * math.cos(x[1]) * math.exp(-((x[0] - self.Easom_center[0])**2 + (x[1] - self.Easom_center[1])**2))
                blocks['5'] = easom
                blockcomments['5'] = "Easom function"
                data_point.obs_dict['Easom'] = easom

        data_point.spc = zslha.SLHA()
        data_point.spc.blocks['TOYTARGETS']=blocks
        data_point.spc.blockcomments['TOYTARGETS']=blockcomments
        return
