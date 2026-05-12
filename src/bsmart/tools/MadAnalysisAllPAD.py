"""
Tool to run madanalysis through the PAD, either via SFS or the regular one


Need to specify in tools:

.. code-block:: text

    "MadAnalysisAllPAD" : {
        "MadGraph" : "/path/to/modeldir/bin/generate_events",
        "MA5 Path" : "/path/to/ma5",
        "Analyses" [list,of,analyses],  OR "MA5 Card": "/path/to/MA5_recasting_card",
        "Process" : "process_name_without_spaces" (optional)
    }



### look in MA5 recast_configuration.py for lots of useful info

"""

__meta__ = {
    "name": "MadAnalysisAllPAD",
    "requires": [],
    "settings": {
        "MadGraph": "Path to generate_events",
        "MA5 Path": "Path to MadAnalysis5",
        "Analyses": "List of analyses",
        "MA5 Card": "Path to MA5 recasting card",
        "Process": "Process name"
    }
}


import os
import shutil

from bsmart import debug
from bsmart.HEPRun import HepTool,DataPoint

#from collider.pythiacfg import PYTHIACFG

from bsmart.collider.MGaux import readxsdata,write_command_file,MadGraphRunner

#import yaml
from bsmart import zslha
#import fnmatch
import sys
#import math

def ReadMASummary(summaryfile):
    RF=open(summaryfile,'r')
    """
            headers should be:
            # dataset name                analysis name               signal region                                               best?     sig95(exp)          sig95(obs)          1-CLs || ...
            
    """

    bestana='None'
    bestregion='None'
    bestsig95=-1
    bestCLs=0.0
    foundbest=False

    for line in RF:
        sline=line.strip()
        if sline.startswith('#'):
            continue
        data=sline.split()

        if len(data) < 7:
            continue
        if int(data[3]) == 1:


            if not foundbest:
                foundbest=True
                bestana=data[1]
                bestregion=data[2]
                bestsig95=float(data[5])
                bestCLs=float(data[6])
            else:
                newbestsig95=float(data[5])
                if newbestsig95 > 0 and newbestsig95 < bestsig95:
                    bestsig95 = newbestsig95
                    bestana=data[1]
                    bestregion=data[2]
                    bestCLs=float(data[6])
    RF.close()

    return { 'Analysis':bestana, 'Region':bestregion, 'CLs': bestCLs, 'sig95': bestsig95}

###############################################################

 
class NewTool(HepTool):
    """ Runs madgraph and extracts the cross-section plus uncertainty """
    def __init__(self, name, settings,global_settings=None):
        HepTool.__init__(self, name, settings,global_settings)
        
        if 'Process' in self.settings:
            self.process = self.settings['Process'] ## to allow labelling the processes
        else:
            self.process='' 

        self.sample_name=self.process+'_sample'
        
        if 'MA5 Path' not in self.settings or not os.path.exists(self.settings['MA5 Path']):
            raise NameError("Need to specify MA5 Path in the tool settings!")

        self.ma5_cmd=os.path.join(self.settings['MA5 Path'],'bin','ma5')
       

        if hasattr(global_settings,'cores') and global_settings.cores > 1:
            raise NameError("MadGraph/MadAnalysis cannot be run in multicore mode!")
        # Settings

        
        
        #self.madgraph_cmd = self.settings['MadGraph']
        
        """
        Set up MadGraph

        """
        try:
            self.MGrunner=MadGraphRunner(settings)
        except Exception as e:
            raise NameError('Failed to initialise MadGraph: '+str(e))

        """ Read the recast_config.dat to get info about the analyses and the cards they require """
        
        self.analysismap={}
        self.analysistypemap={}
        with open(os.path.join(self.settings['MA5 Path'],'tools','PADForSFS','Input','recast_config.dat'),'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue

                data1=line.split('|')
                if len(data1) != 2:
                    continue
                cardname=data1[0].strip()
                analyses_for_card=data1[1].strip().split()
                for ana in analyses_for_card:
                    self.analysismap[ana]=cardname
                    self.analysistypemap[ana]='vSFS'
        with open(os.path.join(self.settings['MA5 Path'],'tools','PAD','Input','recast_config.dat'),'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue

                data1=line.split('|')
                if len(data1) != 2:
                    continue
                cardname=data1[0].strip()
                analyses_for_card=data1[1].strip().split()
                for ana in analyses_for_card:
                    self.analysismap[ana]=cardname
                    self.analysistypemap[ana]='v1.2'
        

        #if 'MA5 Card' not in self.settings:
        #    raise NameError("Need (for now) to name a MA5 Card!")
        self.ma5_card=None
        if 'MA5 Card' in self.settings:
            tpath=os.path.dirname(self.settings['MA5 Card'])
            if tpath == '':
            
                self.ma5_card_name=self.settings['MA5 Card']
                #self.yaml = os.path.join(self.settings['cwd'],self.yamlname)
                #yamldir=self.settings['cwd']
                ma5_card_dir=os.getcwd()
                self.ma5_card = os.path.join(ma5_card_dir,self.ma5_card_name)
            
            
            else:
                self.ma5_card = self.settings['MA5 Card']
                self.ma5_card_name=os.path.split(self.ma5_card)[-1]
                #yamldir=os.path.split(os.path.dirname(self.yaml))[0]  ### I no longer know why we are doing this
                #yamldir = tpath
        elif 'Analyses' in self.settings: ## asseble from list of analyses, like atlas_susy_2016_07
            self.ma5_card=os.path.join(os.getcwd(),self.sample_name+'.card')
            with open(self.ma5_card,'w') as f:
                for analysis in self.settings['Analyses']:
                    if analysis not in self.analysismap:
                        continue
                    f.write(analysis+'    '+self.analysistypemap[analysis]+'    on    '+self.analysismap[analysis]+'   #\n')
        else:
            raise NameError("Need to name an MA5 Card or a list of analyses to use!")

        if not os.path.exists(self.ma5_card):
            raise NameError('Could not find the specified MA5 Card. Either specify the full path to it, or put it in the run directory.')


        

    def run(self, spc_file, temp_dir, log,data_point):
        
        #Run MadGraph

        MGinfo=self.MGrunner.run(spc_file,temp_dir,log,data_point)

        hepmc_events=MGinfo['hepmc']
        if hepmc_events is None:
            log.error('No hepmc events produced')
            raise

        ## set up the MA5 run directory
        tempMA5dir=os.path.join(temp_dir,'MA5')
        templaunchdir=os.path.join(tempMA5dir,'launch')
        if os.path.exists(tempMA5dir):
            shutil.rmtree(tempMA5dir)

        os.makedirs(templaunchdir)

        ## create the sample file
        
        ma5proc=os.path.join(templaunchdir,self.sample_name)

        """
            rm -f ma5.proc
    echo "set main.recast = on" > ma5.proc
    echo "import $pythf as ewset" >> ma5.proc
    echo "set main.recast.store_root = False" >> ma5.proc
    echo "set main.recast.card_path = MA5card.dat" >> ma5.proc
    echo "set ewset.xsection = $xsection" >> ma5.proc
    echo "submit MA5i/$cardstub" >> ma5.proc
    echo "N" >> ma5.proc
        """


        if os.path.exists(ma5proc):
            os.remove(ma5proc)
        with open(ma5proc, 'w') as f:
            f.write('set main.recast = on\n')
            f.write('set main.recast.store_root = False\n')
            f.write('import '+hepmc_events+' as '+self.sample_name+'\n')
            f.write('set main.recast.card_path = '+self.ma5_card+'\n')
            f.write('set '+self.sample_name+'.xsection = %.4e\n' %float(MGinfo['xs']))
            f.write('submit MA5/'+self.sample_name+'\n')
            f.write('N\n')
        
        ## Now to run
        debug.command_line_log(self.ma5_cmd+' -R < '+ma5proc, log)
        
        ## Do something with the results ... take only the 'best' exclusions for now

        clsresultfile=os.path.join(tempMA5dir,self.sample_name,'Output','SAF','CLs_output_summary.dat')

        if os.path.exists(clsresultfile):
            try:
                MA5results=ReadMASummary(clsresultfile)
                ### { 'Analysis':bestana, 'Region':bestregion, 'CLs': bestCLs, 'sig95': bestsig95}
            except Exception as e:
                log.debug('Failed to read MA5 results: '+str(e))

            
            """ now append to the spectrum file """
            OF = open(spc_file,"a")

            OF.write('BLOCK MADANALYSISPAD # Best analysis:'+MA5results['Analysis']+' \n')
            OF.write(' 1 %10.6E # Best CLs \n' % float(MA5results['CLs']))
            sig95=float(MA5results['sig95'])
            OF.write(' 2 %10.6E # Best sigma 95 \n' % sig95)

            if sig95 <= 0:
                xsMinusSigma95=-MGinfo['xs']
                rMinus1=-1
            else:
                xsMinusSigma95=MGinfo['xs']-sig95
                rMinus1=xsMinusSigma95/sig95
            OF.write(' 3 %10.6E # Cross-section - sigma95: > 0 if excluded, < 0 if allowed\n' % (xsMinusSigma95))
            OF.write(' 4 %10.6E # r - 1: > 0 if excluded, < 0 if allowed\n' % (rMinus1))
            OF.write('#10 %s # Best region \n' % MA5results['Region'])
            

        
            OF.close()

            if data_point.spc is None:
                data_point.spc=zslha.read(spc_file)
            else:
            #if data_point.spc is not None:
                data_point.spc.blocks['MADANALYSISPAD']={'1': float(MA5results['CLs']), '2':float(MA5results['sig95']), '3':xsMinusSigma95, '4':rMinus1}
                data_point.spc.blockcomments['MADANALYSISPAD']={'1': ' Best CLs', '2':'Best sigma 95', '3':'Cross-section - sigma95: > 0 if excluded, < 0 if allowed','4':'r - 1: > 0 if excluded, < 0 if allowed'}
               
                    

                



        # Cutflows can be picked up by 'Store Everything'

       
            

