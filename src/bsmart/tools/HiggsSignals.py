"""
Tool to run HiggsSignals


"""

__meta__ = {
    "name": "HiggsSignals",
    "requires": [],
    "settings": {
        "Command": "Path to HiggsSignals executable",
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
            self.uncertainty_file='MHall_uncertainties.dat'
            self.fake_uncertainties=True
            if self.options == '':
                self.options='latestresults 2 effC %d %d' %(self.neutralhiggs,self.chargedhiggs)
        else:
            self.fake_uncertainties=False
        
            
        if 'OutputFile' in self.settings:
            self.output_file = self.settings['OutputFile']
        else:
            self.output_file = 'HiggsSignals_results.dat'
        #self.Pvalue_index=-1
        
    def run(self, spc_file, temp_dir, log,data_point):
        
        if not os.path.exists(spc_file):
            log.debug('No spc; not running HiggsSignals')
            raise
        
        # Clean up previous output
        if os.path.exists(self.output_file):
            os.remove(self.output_file)
        
        if self.fake_uncertainties:
            ### Create the uncertainty file
            if os.path.exists(self.uncertainty_file):
                os.remove(self.uncertainty_file)

            with open(self.uncertainty_file,'w') as UncF:
                UncF.write('1 ')
                for _ in range(self.neutralhiggs):
                    UncF.write(' 5.0')
                for _ in range(self.chargedhiggs):
                    UncF.write(' 3.0')
                UncF.write('\n')

            
        # Run HiggsSignals
        debug.command_line_log(self.command + " " + self.options + " " + temp_dir + "/", log,ctimeout=60)
        
        ### clean up on a succesful run; doesn't matter if we have a failed run!
        #for infile in self.possible_input_files:
        #    #full_infile=os.path.join(temp_dir,infile)
        #    if os.path.exists(infile):
        #        os.remove(infile)
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
            OF.write('BLOCK HIGGSSIGNALS #\n')
            data_point.spc.blocks['HIGGSSIGNALS']={}
            data_point.spc.blockcomments['HIGGSSIGNALS']={}
            pvalue=-1.0
            for i in range(1, len(results)):
                #OF.write(str(i) + " " + str(results[i]) + " #")
                OF.write('%3d  %.8E # '%(i,float(results[i])))
                data_point.spc.blocks['HIGGSSIGNALS'][str(i)] = float(results[i])
                if i < len(labels)-1:
                    OF.write(' '+str(labels[i+1]))
                    data_point.spc.blockcomments['HIGGSSIGNALS'][str(i)] = str(labels[i+1])
                    if str(labels[i+1]) == "Pvalue":
                        pvalue=float(results[i])
                OF.write('\n')
                #print(str(spc_file) + ': ' + str(i) + ' ' +str(results[i]) + ' # '+str(labels[i+1]) )
            #print('Finished writing to spc')
            if pvalue > -1:
                OF.write('101 %.8e # HiggsSignals P value\n' % pvalue) 
                data_point.spc.blocks['HIGGSSIGNALS']['101'] = pvalue
                data_point.spc.blockcomments['HIGGSSIGNALS']['101'] = 'HiggsSignals pvalue'
            OF.close()
        else:
            log.debug("No HiggsSignals output!")
            raise
