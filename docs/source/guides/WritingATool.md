Adding your own tool
====================

The inbuilt tools are contained within python scripts in the package in a subdirectory called [`tools`](../api/bsmart.tools.rst). However, if you want to write your own, `BSMArt` will search the current directory for a directory called `Tools` and look for a python script with the name of the tool you have specified. For example, a tool called `AwesomeTool` would look for a file called `AwesomeTool.py` in the `Tools` directory:

```bash
Tools/
    AwesomeTool.py
```

Then you write your json file accordingly:





For example, this is the code for the `FlexibleSUSY` tool:

```python

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
        

```



In addition, for documentation you can include:

- A top-level docstring; this should be of the form """
- A `__meta__` dictionary describing the scan’s metadata. This typically consists of:
    - `name`: The name of the tool.
    - `requires`: A list of required packages.
    - `Settings`: A dictionary of settings for the tool.



For information about general tool settings see {py:mod}`bsmart.HEPRun`