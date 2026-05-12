"""
Tool to run Vevacious++


Requires in the tool settings:

* 'Initialization File' with an *absolute path*. You will have to set that up yourself, unless it is for THDM or the MSSM. See the Vevacious documentation

* 'Command' with an *absolute path* for the vevacious++ executable
 
Optional: if you specify the name of the input xml file we use that, otherwise it will create a file 'VevaciousInput.xml'

NB: some question about thread safety but it may be ok

"""

__meta__ = {
    "name": "VevaciousPlusPlus",
    "requires": [],
    "settings": {
        "Command": "Path to executable",
        "Initialization File": "Path to init file",
        "InputFile": "Path to input XML"
    }
}


import os
from bsmart import debug
from bsmart.HEPRun import HepTool,DataPoint
from bsmart import zslha

class NewTool(HepTool):
    """ overload the init to initialise the DM limit calculator, but name and settings are already given in HepTool """
    def __init__(self, name, settings,global_settings=None):
        HepTool.__init__(self, name, settings,global_settings)        
        
        self.command = self.settings['Command']
        self.initfile=self.settings['Initialization File']
        
        if 'InputFile' in self.settings:
            self.inputfile=self.settings['InputFile']
        else:
            self.inputfile='VevaciousInput.xml'

        self.spcname=''
        
            
    def run(self, spc_file, temp_dir, log,data_point):
        
        if not os.path.exists(spc_file):
            log.debug('No spc; not running Vevacious')
            raise

        ## Create the input file for first use, don't need to overwrite it 
        if not os.path.exists(self.inputfile):
            tpath=os.path.join(temp_dir,spc_file)
            if tpath != self.spcname:
                self.spcname=tpath
                
                with open(self.inputfile,'w') as IF:
                    IF.write('<VevaciousPlusPlusMainInput>\n')
                    IF.write('  <InitializationFile>\n')
                    IF.write('   '+self.initfile+'\n')
                    IF.write('  </InitializationFile>\n\n')
                    IF.write('  <SingleParameterPoint>\n')
                    IF.write('     <RunPointInput>\n')
                    IF.write('   '+tpath+'\n')
                    IF.write('     </RunPointInput>\n')
                    IF.write('     <OutputFilename>\n')
                    IF.write('   '+os.path.join(temp_dir,'Vevacious.vout')+'\n')
                    IF.write('     </OutputFilename>\n')
                    IF.write('     <AppendLhaOutputToLhaInput></AppendLhaOutputToLhaInput>\n')
                    IF.write('  </SingleParameterPoint>\n\n')
                    IF.write('</VevaciousPlusPlusMainInput>\n')
                
        # Run Vevacious, don't allow it to run for more than two minutes
        if debug.psutil_available:  # single core if psutil is available
            debug.onecore_command_log(self.command + ' '+self.inputfile, log,ctimeout=120)
        else:
            debug.command_line_log(self.command + ' '+self.inputfile, log,ctimeout=120)

        ## Output should be appended to the SLHA file: need to update the data point if it is there
        if os.path.exists(spc_file):
            if data_point is None:
                data_point.spc=zslha.read(spc_file)
            else:
                updated_data=zslha.read(spc_file,None,None,'BLOCK VEVACIOUSSTABILITY')
                data_point.spc.add(updated_data)
