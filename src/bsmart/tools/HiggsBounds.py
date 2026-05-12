"""
Tool to run HiggsBounds


"""

__meta__ = {
    "name": "HiggsBounds",
    "requires": [],
    "settings": {
        "Command": "Path to HiggsBounds executable",
        "Options": "Optional command line arguments",
        "Neutral Higgs": "Number of neutral Higgs",
        "Charged Higgs": "Number of charged Higgs",
        "OutputFile": "Path to output file"
    }
}

import os
import shutil
#import subprocess
from bsmart import debug
from bsmart.HEPRun import HepTool, DataPoint
from bsmart import zslha

class NewTool(HepTool):
    """ overload the init, name and settings are already given in HepTool """
    def __init__(self, name, settings,global_settings=None):
        HepTool.__init__(self, name, settings,global_settings)
    
        self.command = self.settings['Command']

        if 'Options' in self.settings: ## e.g. LandH effC 3 1
            self.options=self.settings['Options']
        else:
            self.options=''

        if 'Neutral Higgs' in self.settings and 'Charged Higgs' in self.settings:
            #### Override the options line ...
            self.neutralhiggs=int(self.settings['Neutral Higgs'])
            self.chargedhiggs=int(self.settings['Charged Higgs'])
            if self.options == '':
                self.options='LandH effC %d %d' %(self.neutralhiggs,self.chargedhiggs)
            
        if 'OutputFile' in self.settings:
            self.output_file = self.settings['OutputFile']
        else:
            self.output_file = 'HiggsBounds_results.dat'
        
        #self.possible_input_files=['BR_H_NP.dat','BR_Hplus.dat','BR_t.dat','effC.dat','Key.dat','MH_GammaTot.dat','MHplus_GammaTot.dat']
        
    def run(self, spc_file, temp_dir, log,data_point):
        
        if not os.path.exists(spc_file):
            log.debug('No spc; not running HiggsBounds')
            raise
        
        # Clean up previous output
        if os.path.exists(self.output_file):
            os.remove(self.output_file)
        
                
        # Run HiggsBounds
        debug.command_line_log(self.command + " " + self.options + " " + temp_dir + "/", log,ctimeout=60)
        
        # ### clean up on a succesful run; doesn't matter if we have a failed run!
        # for infile in self.possible_input_files:
        #     #full_infile=os.path.join(temp_dir,infile)
        #     if os.path.exists(infile):
        #         os.remove(infile)
        # Reading the Output file
        if os.path.exists(self.output_file):
            
            labels=[]
            for line in open(self.output_file):
                li = line.strip()
                if not li.startswith("#"):
                    results = list(filter(None, line.rstrip().split(' ')))
                elif li.startswith("#cols:"):
                    labels= list(filter(None, line.rstrip().split(' ')))
            
                # Append output to the SPheno file
            ## There should be only one line of valid output, so keep the last one and transform it into a block

            try:
                OF=open(spc_file,"a")
            except:
                log.debug('Failed to open '+str(spc_file)+' for appending')
                raise
            if data_point.spc is None:
                data_point.spc=zslha.SLHA()
            OF.write('BLOCK HIGGSBOUNDS #\n')
            data_point.spc.blocks['HIGGSBOUNDS']={}
            data_point.spc.blockcomments['HIGGSBOUNDS']={}
            obsratio=-1.0
            for i in range(1, len(results)):
                #OF.write(str(i) + " " + str(results[i]) + " #")
                OF.write('%3d  %.8E # '%(i,float(results[i])))
                data_point.spc.blocks['HIGGSBOUNDS'][str(i)] = float(results[i])
                
                if i < len(labels)-1:
                    OF.write(' '+str(labels[i+1]))
                    data_point.spc.blockcomments['HIGGSBOUNDS'][str(i)] = str(labels[i+1])
                    if str(labels[i+1]) == 'obsratio':
                        obsratio=float(results[i])
                OF.write('\n')
                #print(str(spc_file) + ': ' + str(i) + ' ' +str(results[i]) + ' # '+str(labels[i+1]) )
            #print('Finished writing to spc')
            if obsratio > -1:
                OF.write('101 %.8E # HiggsBounds obs ratio\n' % obsratio) 
                data_point.spc.blocks['HIGGSBOUNDS']['101'] = obsratio
                data_point.spc.blockcomments['HIGGSBOUNDS']['101'] = 'HiggsBounds obs ratio'
            OF.close()
        else:
            log.debug("No HiggsBounds output!")
            raise
        """
        if os.path.exists(output_file):
            debug.command_line_log("echo \"Block HIGGSBOUNDS # \" >> "
                                + os.path.join(temp_dir,spc_file), log)
            for line in open(output_file):
                li = line.strip()
                if not li.startswith("#"):
                    results = list(filter(None, line.rstrip().split(' ')))
                # Append output to the SPheno file
            ## There should be only one line of valid output, so keep the last one and transform it into a block
            for i in range(1, len(results)):
                debug.command_line_log("echo \"" + str(i) + " " + str(results[i])
                                    + " # \" >> " + os.path.join(temp_dir,spc_file), log)
        else:
            log.error("No HiggsBounds output!",
                    self.command + " " + self.options + " " + temp_dir)
            raise
        """
