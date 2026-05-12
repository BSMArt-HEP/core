__meta__ = {
    "name": "VacStab",
    "requires": [],
    "settings": {
        "Command": "Path to executable",
        "InputFile": "Input file"
    }
}

import os
import shutil

from bsmart import debug
from bsmart.HEPRun import HepTool,DataPoint

class NewTool(HepTool):
    """ Do not need to overload init """


    
    def run(self, spc_file, temp_dir, log,data_point):
        # Settings
        

        ## Don't need an output file, will append to the SLHA
        #output_file = self.settings['OutputFile']
        
        

        if spc_file != self.settings['InputFile']:
            shutil.copyfile(spc_file,self.settings['InputFile'])
            command = self.settings['Command']+' '+os.path.join(temp_dir,self.settings['InputFile'])+' '+os.path.join(temp_dir,spc_file)
        else:
            command = self.settings['Command']+' '+os.path.join(temp_dir,self.settings['InputFile']) ## will automatically append
            
        # Run 
        debug.command_line_log(command, log)
        
        

       
