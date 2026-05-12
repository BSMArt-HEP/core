"""
Tool to run the pipeline madgraph to create lhe file then hackanalysis to drive pythia and perform analysis. 


Need to specify in tools:

.. code-block:: text

    "MadGraphHackAnalysis" : {
        "HackAnalysis" : "/path/to/executable", (e.g. analysePYTHIA_LHE.exe)
        "YAML file" "<NAME>.yaml"    (file for configuring HackAnalysis run  )
    }

plus, in addition, settings for running madgraph which can also be found in ../collider/MGaux.py
These include as obligatory:

* "MadGraph" : "/path/to/modeldir/bin/generate_events",   (if you have already generated the output directory for your process)

OR

* "Proc Card": "proc_card.dat",    ( for generating the madgraph directory based on a series of commands to provide to madgraph)
* "MadGraph Path": "/path/to/MG5_aMC_v3.5.XX/bin/mg5_aMC",     (so that we can run MadGraph )
* "Output Directory": "/path/to/where/you/want/to/put/the/output", (this is OPTIONAL: if you have the command output <outputname> in your proc_card.dat then BSMArt will find it)

other optional MadGraph settings include:

.. code-block:: text

    {
        "python2": "True",  (if you need to run madgraph via python2)
        "MadGraph Directory": "/path/to/modeldir",  (if you want to specify the path to the modeldir, otherwise the code will infer it)
        "QNUMBERS file": "myQNUMBERs.slha", (for models that have particles unknown to pythia)
        "MadSpin": "True", (if you want to run MadSpin and have it handled by BSMArt)
        "k-factor": 1.2, (for multiplying the cross-section by a k-factor)
        "Cross-sections" : "/path/to/csv/file/containing/cross-sections.csv",      (if you have tabulated cross-sections as a function of *one* variable to use instead of the one computed in MadGraph; this is documented in the InterpolateXSection.py file.)
    }

---------------------------

25/03/13: Make more robust, in case for some reason HackAnalysis fails to run when doing batches, will give two more tries with the same point






"""

__meta__ = {
    "name": "MadGraphHackAnalysis",
    "requires": ["pyyaml"],
    "settings": {
        "HackAnalysis": "Path to HackAnalysis executable",
        "YAML file": "Path to YAML config file",
        "MadGraph": "Path to generate_events",
        "Proc Card": "Path to proc_card.dat",
        "MadGraph Path": "Path to mg5_aMC",
        "Output Directory": "Path to output directory",
        "Mode": "LHE or HEPMC",
        "GridPack": "Boolean",
        "Check Absolute Convergence": "Boolean",
        "Check Convergence": "Boolean",
        "Efficiency Threshold": "Float",
        "Convergence Margin": "Float",
        "Check Exclusion": "Boolean",
        "Batch Size": "Int",
        "Store Events": "Boolean",
        "Cores": "Int"
    }
}


import os
import shutil

from bsmart import debug
from bsmart.HEPRun import HepTool,DataPoint

from bsmart.collider.pythiacfg import PYTHIACFG
from bsmart.collider.splitlhe import splitlhe
from bsmart.collider.MGaux import MadGraphRunner

import yaml
import json
from collections import OrderedDict
from bsmart import zslha
#import fnmatch
import sys
#import math

try:
    import importlib
except Exception as e:
    print("Failed to import importlib, "+str(e))


###############################################################


references=r'''@article{Goodsell:2021iwc,
    author = "Goodsell, Mark D. and Priya, Lakshmi",
    title = "{Long dead winos}",
    eprint = "2106.08815",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    doi = "10.1140/epjc/s10052-022-10188-1",
    journal = "Eur. Phys. J. C",
    volume = "82",
    number = "3",
    pages = "235",
    year = "2022"
}

@article{Goodsell:2024aig,
    author = "Goodsell, Mark D.",
    title = "{HackAnalysis 2: A powerful and hackable recasting tool}",
    eprint = "2406.10042",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    month = "6",
    year = "2024"
}

@article{Corpe:2025sbw,
    author = "Corpe, Louie and Haddad, Abdelhamid and Goodsell, Mark",
    title = "{Recasting the ATLAS search for displaced hadronic jets in the ATLAS calorimeter with additional jets or leptons using surrogate models}",
     eprint = "2502.10231",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    month = "2",
    year = "2025"
}
'''


 
class NewTool(HepTool):
    """ MadGraphHackAnalysis runs madgraph and extracts the cross-section plus uncertainty 
    
    Various settings are possible:


    
    
    """
    def __init__(self, name, settings,global_settings=None):
        HepTool.__init__(self, name, settings,global_settings)
        self.citations=references
        ## Optional "Process" setting to label the process examined
        if 'Process' in self.settings:
            self.process = self.settings['Process'] ## to allow labelling the processes
        else:
            self.process='' 

        
        ## Optional "Mode" setting to choose "LHE" or "HEPMC", otherwise infer from the command name
        if 'Mode' in self.settings:
            
            tmode=self.settings['Mode'].upper()
            if tmode in ['LHE','HEPMC']:
                self.mode = tmode
            else:
                raise NameError("Operating mode not recognised; allowed values are LHE and HEPMC")
        else:
            if 'analysePYTHIA_LHE' in self.settings['HackAnalysis']:
                self.mode='LHE'
            elif 'analyseHEPMC' in self.settings['HackAnalysis']:
                self.mode='HEPMC'
            else:
                raise NameError("Executable not recognised and Mode not specified; allowed values are LHE and HEPMC")

        

        if global_settings is not None and global_settings.varnames is not None:
            self.varnames=global_settings.varnames
        else:
            self.varnames=None

        ## Require "YAML file" setting for running HackAnalysis
        if 'YAML file' not in self.settings:
            #self.log.error('No YAML file present!')
            raise NameError('No YAML file present!')
        tpath=os.path.dirname(self.settings['YAML file'])
        if tpath == '':
            
            self.yamlname=self.settings['YAML file']
            #self.yaml = os.path.join(self.settings['cwd'],self.yamlname)
            #yamldir=self.settings['cwd']
            yamldir=os.getcwd()
            self.yaml = os.path.join(yamldir,self.yamlname)
            
        else:
            self.yaml = self.settings['YAML file']
            self.yamlname=os.path.split(self.yaml)[-1]
            #yamldir=os.path.split(os.path.dirname(self.yaml))[0]  ### I no longer know why we are doing this
            yamldir = tpath
        # READ OTHER INFO FROM THE YAML!!!!
        try:
            with open(self.yaml, 'r') as stream:        
                self.yaml_data = yaml.safe_load(stream)
        except:
            raise NameError('Failed to open HackAnalysis yaml file!')
        
        self.qnumbers_name=None
        if 'QNUMBERS file' in self.yaml_data:
            ## need to know whether we have an absolute path for the QNUMBERS file, or no path:
            ## if there's no path then we need to copy it to the temp directory for running HackAnalysis
            ## We can extract it from the MG runner 
            tpath=os.path.dirname(self.yaml_data['QNUMBERS file'])
            if tpath == '':  # need to copy the QNUMBERS file to the working directory
                self.qnumbers_name=self.yaml_data['QNUMBERS file']
                if 'QNUMBERS file' not in self.settings:
                    self.settings['QNUMBERS file'] = self.yaml_data['QNUMBERS file']
            #else:
            #    self.qnumbers_name=None
                
        ## "GridPack": "True" enables gridpack mode, perfect for running large numbers of events
        ## GridPack not yet supported for HEPMC run because we only run single core HEPMC in HA2
        if 'GridPack' in self.settings and self.mode=='LHE':
            self.useGridPack=eval(self.settings['GridPack'])
        else:
            self.useGridPack=False


        ## "Check Absolute Convergence" setting when running batches of events, to check that cutflows converge to fixed precision
        if "Check Absolute Convergence" in self.settings:
            self.check_convergence = eval(self.settings['Check Absolute Convergence'])
        else:
            self.check_convergence = False

        ## "Check Convergence" setting when running batches of events, to check that cutflows converge to suitable precision depending on experimental sensitivity
        if "Check Convergence" in self.settings:
            self.check_exp_convergence = eval(self.settings['Check Convergence'])
            if self.check_exp_convergence:
                self.check_convergence = True
                
        else:
            self.check_exp_convergence=False

        if "Efficiency Threshold" in self.settings:
            self.convergence_threshold=float(self.settings['Efficiency Threshold'])
        else:
            self.convergence_threshold=1e-5

        if "Convergence Margin" in self.settings:
            self.convergence_error=float(self.settings['Convergence Margin'])
            self.safety_ratio=self.convergence_error
        else:
            self.convergence_error=0.1
            self.safety_ratio=0.25

        HAdir=os.path.dirname(settings['HackAnalysis'])
        #HApythonpath=os.path.join(os.path.dirname(settings['HackAnalysis']),'Statistics')
        HApythonpath=os.path.join(HAdir,'Statistics')
        if not os.path.exists(HApythonpath): # presumably using new cmake, so the stats directory is different
            if 'bin' in HAdir:
                HAdir=os.path.dirname(HAdir)
                HApythonpath=os.path.join(HAdir,'Statistics')

        # Modification 7/2/2025
        #HAstatspath=os.path.join(HApythonpath,'HAStats')
        HAstatspath=HApythonpath
        
        if "Check Exclusion" in self.settings:
            self.check_exclusion = eval(self.settings['Check Exclusion'])   
            if self.check_exclusion:
                #print('Adding stats path '+HAstatspath)

                if os.path.exists(HAstatspath):
                    if HApythonpath not in sys.path:
                        sys.path.insert(0,HApythonpath)
                #self.exclusions_module=importlib.import_module('HAstats') 
                self.exclusions_module=importlib.import_module('HAStats.HAstats')
                self.analysishandlers={}
                infopath=os.path.join(HApythonpath,"HAStats","Analyses")
                for analysis in self.yaml_data['analyses']:
                    analysisUPPER=analysis.upper()
                    try:
                        
                        with open(os.path.join(infopath,analysisUPPER+'_info.json'),'r') as thisinfo:
                            analysis_infodict=json.load(thisinfo)
                        
                        self.analysishandlers[analysis] = self.exclusions_module.AnalysisStats(analysis_infodict)

                    except Exception as e:
                        print("Failed to initialise analysis %s: %s" %(analysis,str(e)))
                        if analysis in self.analysishandlers:
                            print("Still somehow initialised the analysis, deleting")
                            self.analysishandlers.pop(analysis)
                        continue
        else:
            self.check_exclusion=False



        self.nevents=None
        
        ## We get nevents from the yaml file if it is there, otherwise use default of 10000
        if 'nevents' in self.yaml_data['settings']:
            self.nevents=int(self.yaml_data['settings']['nevents'])

        if self.nevents is None and self.useGridPack:
            
            self.nevents=10000
        
        self.combiner = None
        ## "Batch Size" setting if you want batch mode
        if 'Batch Size' in self.settings:
            self.batch_size=int(self.settings['Batch Size'])
            if self.batch_size < self.nevents:
                # need to modify the yaml file to only run the batch-size number of events
                self.yaml_data['settings']['nevents'] = self.batch_size
                # Import routines for merging json files etc
                
                """
                HApythonpath=os.path.join(os.path.dirname(settings['HackAnalysis']),'Statistics')
                if os.path.exists(HApythonpath):
                    sys.path.insert(0,HApythonpath)
                """
                
                if os.path.exists(HApythonpath):
                    if HApythonpath not in sys.path:
                        sys.path.insert(0,HApythonpath)
                    self.combiner=importlib.import_module('combine_ha')
                    if self.check_convergence:
                        self.convergence_module=importlib.import_module('ha_convergence')
                        self.convergence=self.convergence_module.HAconverger()
                else:
                    raise NameError("Cannot find Statistics path "+HApythonpath)
        else:
            self.batch_size=self.nevents

        if 'JSON file' in self.yaml_data['settings']:
            self.HAjson=self.yaml_data['settings']['JSON file']
        else:
            self.HAjson='hackanalysis.json'

        
        
        #Set up MadGraph

        


        #if self.mode=='LHE' and 'Extra Commands' not in self.settings:
        if 'Extra Commands' not in self.settings:    
            # need to make sure showering is turned off in madgraph, if using default settings
            # Otherwise we trust the user to know what they are doing!!
            if self.mode=='LHE':
                commandlist=[['shower=OFF'],['detector=OFF'],['analysis=OFF']]
            else:
                commandlist=[['shower=Pythia8'],['detector=OFF'],['analysis=OFF']]
            
            
            if 'MadSpin' in self.settings and eval(self.settings['MadSpin']):
                commandlist.append(['madspin=ON'])
            else:
                commandlist.append(['madspin=OFF'])
            commandlist.append(['done'])
            if 'run_card Commands' in self.settings:
                ## 'run_card Commands' in settings are commands passed to madgraph to update the run_card or parameters
                # These are in addition to the automatically created ones (e.g. setting the number of events) 
                for rccommand in self.settings['run_card Commands']:
                    commandlist.append(rccommand)
            else:
                commandlist.append(['update dependent'])

            #if 'MadSpin' in self.settings and eval(self.settings['MadSpin']):  
                #self.settings['Extra Commands']=[['shower=OFF'],['detector=OFF'],['analysis=OFF'],['madspin=ON'],['done'],[eventsline],['done']]
            #else:
                #self.settings['Extra Commands']=[['shower=OFF'],['detector=OFF'],['analysis=OFF'],['madspin=OFF'],['done'],[eventsline],['done']]
            
            # if we want to use the GridPack mode but in one big batch, we only need a minimum number of initial events
            # Otherwise we set the initial events equal to the batch size (which could be equal to the total number of events in one-shot runs)
            if self.useGridPack and self.batch_size >= self.nevents: 
                commandlist.append(['set nevents 10000']) ## basic number of events for gridpack run
            else:
                commandlist.append(['set nevents '+str(self.batch_size)])
            """

            if self.useGridPack:
                commandlist.append(['set nevents 10000']) ## basic number of events for gridpack run
            else:
                if 'nevents' in self.yaml_data['settings']:
                    commandlist.append(['set nevents '+str(self.yaml_data['settings']['nevents'])])
            """

            commandlist.append(['done'])
            self.settings['Extra Commands']=commandlist

        
        
        if 'Store Events' in self.settings:
            self.store_events=eval(self.settings['Store Events'])
        else:
            self.store_events = False


        self.hackanalysis_cmd = self.settings['HackAnalysis']+' '+self.yamlname

        
        if 'Cores' in self.settings:
            self.n_cores=int(self.settings['Cores'])
            self.yaml_data['settings']['cores'] = self.n_cores
        else:
            if 'cores' in self.yaml_data['settings']:
                self.n_cores = int(self.yaml_data['settings']['cores'])
                self.settings['Cores']=self.n_cores
            else:
                self.n_cores = 1
                self.yaml_data['settings']['cores'] = 1

        ## We use the temporary directory for the madgraph model directory if not told otherwise
        if global_settings is not None:
            #print("have global settings!")
            if hasattr(global_settings,'work_dir') and global_settings.work_dir is not None:
                #print("Setting TempDir!")
                self.settings['TempDir'] = global_settings.work_dir
                #print(self.settings['TempDir'])

        try:
            self.MGrunner=MadGraphRunner(self.settings)
        except Exception as e:
            raise NameError('Failed to initialise MadGraph: '+str(e))


        if 'LHE file' in self.yaml_data['settings']:
            self.LHEfile=self.yaml_data['settings']['LHE file']
        else:
            self.LHEfile='Split/split'

        if self.mode == 'LHE':
            if 'Config file' not in self.yaml_data['settings']:
                raise NameError('No cfg file for pythia present!')

            self.cfgname=self.yaml_data['settings']['Config file']
            #self.cfgfile = os.path.join(self.settings['cwd'],self.cfgname)
            self.cfgfile =os.path.join(yamldir,self.cfgname)
            self.pythia_cfg = PYTHIACFG()
            try:
                self.pythia_cfg.read(self.cfgfile)
            except Exception as e:
                raise NameError("Failed to initialise pythia configuration file, " + str(e))
        
        if 'Efficiency Filename' in self.yaml_data['settings']:
            self.effname = self.yaml_data['settings']['Efficiency Filename']
        else:
            self.effname = 'test.eff'
            self.yaml_data['settings']['Efficiency Filename'] = self.effname

        if 'Cutflow Filename' in self.yaml_data['settings']:
            self.cutflowname = self.yaml_data['settings']['Cutflow Filename']
        else:
            self.cutflowname = 'test_cf.eff'
            self.yaml_data['settings']['Cutflow Filename'] = self.cutflowname
        
        if 'Histogram Filename' in self.yaml_data['settings']:
            self.histoname = self.yaml_data['settings']['Histogram Filename']
        else:
            self.histoname = 'test.yoda'
            self.yaml_data['settings']['Histogram Filename'] = self.histoname
        


    def run(self, spc_file, temp_dir, log,data_point):
        ## set up madgraph configuration
        #var_dict = data_point.get_var_dict(self.varnames)
        ## HackAnalysis yaml file:
        local_yaml = os.path.join(temp_dir, self.yamlname)
        
        if self.useGridPack:
            # First remove any existing gridpack if it already exists
            gridpacktarball=os.path.join(self.MGrunner.MG_Model_Dir,'run_01_gridpack.tar')
            if os.path.exists(gridpacktarball):
                os.remove(gridpacktarball)
            gridpacktarball=os.path.join(self.MGrunner.MG_Model_Dir,'run_01_gridpack.tar.gz')
            if os.path.exists(gridpacktarball):
                os.remove(gridpacktarball)

        MGinfo=self.MGrunner.run(spc_file,temp_dir,log,data_point)
        
        # if the reported cross-section is too small, is there any point in running?
        if self.check_exp_convergence:
            if MGinfo['xs']*1e3*139 < 1:
                if self.check_exclusion: # point clearly not excluded
                    ### ADD CODE HERE!
                    log.info("Not trivially Excluded!") 

                #return
        ## now split lhe file, if needed
        # otherwise we write the full path to the lhe file in the yaml file
        
        if self.mode == 'LHE':
                lhe_events=None
                for tfile in os.listdir(self.MGrunner.MG_Run_Dir):
                    if tfile[-7:] == '.lhe.gz' or tfile[-4:] == '.lhe':
                        lhe_events=os.path.join(self.MGrunner.MG_Run_Dir,tfile)                    
                        break
                if not os.path.exists(lhe_events):
                    raise NameError("Failed to generate lhe events!")
                
                if self.n_cores > 1:
                    splitlhe(self.n_cores,lhe_events,os.path.join(temp_dir,self.LHEfile))
                else:
                    self.yaml_data['settings']['LHE file'] = lhe_events
                    #if os.path.exists(os.path.join(self.MG_Run_Dir,))
        else: ## HEPMC 
                hepmc_events=None
                for tfile in os.listdir(self.MGrunner.MG_Run_Dir):
                    if tfile[-9:] == '.hepmc.gz' or tfile[-6:] == '.hepmc':
                        hepmc_events=os.path.join(self.MGrunner.MG_Run_Dir,tfile)                    
                        break
                if not os.path.exists(hepmc_events):
                    raise NameError("Failed to generate hepmc events!") 
                self.yaml_data['settings']['HEPMC file'] = hepmc_events
            ## Now for procesing the files to set up for HA

        with open(local_yaml,'w') as stream:
                    whatever = yaml.dump(self.yaml_data,stream)


        var_dict = data_point.get_all_dict()
        # Generate pythia config file each time (since we might need to, e.g. for changing qcut )

        ## Now for HackAnalysis stuff

        if self.mode == 'LHE':
            local_cfg = os.path.join(temp_dir, self.cfgname)
            #self.pythia_cfg.write(local_cfg,data_point.inputs,self.varnames)
        
            self.pythia_cfg.write_dict(local_cfg,var_dict)

        ## Now to run
        if os.path.exists(self.HAjson):
            os.remove(self.HAjson)
        if os.path.exists(self.effname):
            os.remove(self.effname)
        if os.path.exists(self.cutflowname):
            os.remove(self.cutflowname)
        if self.histoname is not None:
            if os.path.exists(self.histoname):
                os.remove(self.histoname)

        os.chdir(temp_dir)
        if self.qnumbers_name is not None:
            if not os.path.exists(self.qnumbers_name):
                with open(self.qnumbers_name,'w') as QNF:
                    QNF.write(self.MGrunner.QNUMBERSfile)
        
        failed_attempts=0
        try:
            debug.command_line_log(self.hackanalysis_cmd, log)
        except Exception as e:
            failed_attempts=1
            log.error('Error running HackAnalysis: '+str(e))
        

        if not os.path.exists(self.effname) or not os.path.exists(self.cutflowname):
                ## point has failed
                failed_attempts=1
                log.error('Failed to run HackAnalysis, efficiency or cutflow missing')
                raise NameError('Failed to run HackAnalysis, efficiency or cutflow missing')
        output_json_data=OrderedDict()
        if os.path.exists(self.HAjson):
                    with open(self.HAjson,'r') as jsonfile:
                        output_json_data=json.load(jsonfile,object_pairs_hook=OrderedDict)
        else:
            failed_attempts=1
            log.error('Json file missing')
            #raise NameError('Failed to run HackAnalysis, json missing')

        ### If gridpack mode, now we need the gridpack & batch run
        needed_extra_events=self.nevents
        improve_regions={}
        converged=False
        doneevents=0
        ## here we use the number of  events we *should* have generated. Maybe better to read the number from the output

        if failed_attempts ==0:
            if self.useGridPack and self.batch_size >= self.nevents: 
                doneevents=10000
            else:
                doneevents=self.batch_size

            if self.check_convergence: # check convergence or exclusion
                        if self.check_exp_convergence:
                            #check_convergence_exclusion(self,hackanalysis_output,nloxs=None,input_max_batch_size=None,safety_ratio=0.25,verbose=False)
                            needed_extra_events,improve_regions,_=self.convergence.check_convergence_exclusion(output_json_data,None,self.batch_size,self.safety_ratio,True)
                        else:
                            #check_convergence_fixed_goal(self,efficiencydata,rel_uncertainty=0.1,eff_threshold=1e-5,input_max_batch_size=None,verbose=True)
                            #needed_extra_events,improve_regions,_=self.convergence.check_convergence_fixed_goal(output_json_data,1.e-1,1.e-5,self.batch_size)
                            needed_extra_events,improve_regions,_=self.convergence.check_convergence_fixed_goal(output_json_data,self.convergence_error,self.convergence_threshold,self.batch_size,True)

                        if needed_extra_events > 0:
                            outstr=''
                            log.debug('Estimate need to generate %d more events' %needed_extra_events)
                            for analysis in improve_regions:
                                linestring=','.join(improve_regions[analysis])
                                outstr=outstr+analysis+': '+linestring+'\n'
                            log.debug('Improve regions: \n'+outstr)

                        else:
                            converged=True

        if self.useGridPack and not converged:
            # First need to generate the GridPack
            try:
                if self.MGrunner.python2:
                    debug.command_line_log('python2 '+os.path.join(self.MGrunner.MG_Model_Dir,'bin','madevent')+' create_gridpack',log)
                else:
                    debug.command_line_log(os.path.join(self.MGrunner.MG_Model_Dir,'bin','madevent')+' create_gridpack',log)
            except Exception as e:
                log.error('Could not run madevent to create gridpack: '+str(e))
            
            # then check gridpack was created in the model directory (that's where it gets put ...)
            # since we remove the run_01 should always be called run_01
            gridpacktarball=os.path.join(self.MGrunner.MG_Model_Dir,'run_01_gridpack.tar.gz')
            if not os.path.exists(gridpacktarball):
                gridpacktarball=os.path.join(self.MGrunner.MG_Model_Dir,'run_01_gridpack.tar')
                if not os.path.exists(gridpacktarball):
                    raise NameError('No gridpack found at '+gridpacktarball)

            # now extract it
            log.debug('Extracting gridpack '+gridpacktarball)
            log.debug('Changing directory to '+temp_dir)
            os.chdir(temp_dir) # just to make sure ...

            madeventdir=os.path.join(temp_dir,'madevent')
            if os.path.exists(madeventdir): # clean up if it is already there
                debug.command_line_log('chmod a+w -R madevent',log)
                shutil.rmtree(madeventdir)

            #now unpack it    
            debug.command_line_log('tar -xf '+gridpacktarball,log)
            
            os.chdir(temp_dir)

            # make read only
            debug.command_line_log('chmod a-w -R madevent',log)

            # if python2 -> correct the run.sh script
            if self.MGrunner.python2:
                OFrun=open('run.sh','w')
                OFrun.write('#!/bin/bash\n')
                OFrun.write('basedir=%s\n'%temp_dir)
                OFrun.write('export  LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${basedir}/madevent/lib:${basedir}/HELAS/lib\n')
                OFrun.write('python2 ${basedir}/madevent/bin/gridrun $1 $2 1\n')
                OFrun.write('mv ./Events/GridRun_${2}/unweighted_events.lhe.gz events.lhe.gz\n')
                OFrun.write('rm -rf Events Cards P* *.dat randinit &> /dev/null\n')
                OFrun.close()

            ## remove unnecessary files? 
            
            # Now create the directories for the gridpacks to run
            # New: run directories are created each time we run the gridpack

            import random
            import subprocess
            #seedstart=random.randint(1,sys.maxsize-1)
            seedstart=random.randint(1,1000000)
            log.debug('Setting random seed start to %d' % seedstart)
            if self.batch_size < self.nevents:
                coreevents=int(self.batch_size/self.n_cores)+1
            else:
                coreevents=int(self.nevents/self.n_cores)+1


            var_dict = data_point.get_all_dict()

            if self.mode == 'LHE':
                local_cfg = os.path.join(temp_dir, self.cfgname)
                #self.pythia_cfg.write(local_cfg,data_point.inputs,self.varnames)
                self.pythia_cfg.write_dict(local_cfg,var_dict)

            def rungp(np,seed):
                log.debug('starting run for core %d, requesting %d events with random seed %d' %(np,coreevents,seed+np))
                ## added the /dev/null to clear out the mess
                command='../run.sh %d %d' %(coreevents,seed+np)
                #command='../run.sh %d %d' %(coreevents,seedstart+np)
                wdir=os.path.join(temp_dir,'gpr_%d' %np)
                if os.path.exists(wdir):
                    shutil.rmtree(wdir)
                os.mkdir(wdir)
                """
                command_line = subprocess.Popen(command, shell=True, cwd=wdir,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
                out, err = command_line.communicate()
                """

                # new 22/07/24, don't want the output! Direct to null in a platform indep way
                command_line = subprocess.Popen(command, shell=True, cwd=wdir,
                                    stdout=subprocess.DEVNULL,
                                    stderr=subprocess.PIPE)
                _, err = command_line.communicate()
                

                if err != b'':
                    log.error(err.decode("utf-8"))
                # check events exist
                eventfile=os.path.join(wdir,'events.lhe.gz')
                ending='.lhe.gz'
                if not os.path.exists(eventfile):
                    eventfile=os.path.join(wdir,'events.lhe')
                    ending='.lhe'
                    if not os.path.exists(eventfile):
                        raise NameError('Failed to produce events for core %d' %np)
                    
                # move file to desired location
                
                shutil.move(eventfile,self.LHEfile+'_'+str(np)+ending)
                ## clean up
                shutil.rmtree(wdir)
                # return the name of the file produced
                return self.LHEfile+'_'+str(np)+ending

            # now have already set doneevents and converged check
            #doneevents=0
            #converged = False

            while doneevents < self.nevents and not converged:
                seedstart+=self.n_cores
                doneevents += self.batch_size

                if self.n_cores > 1:
                    import multiprocessing as mp
                    import time
                
                    lhe_events=None
                    for tfile in os.listdir(self.MGrunner.MG_Run_Dir):
                        if tfile[-7:] == '.lhe.gz' or tfile[-4:] == '.lhe':
                            lhe_events=os.path.join(self.MGrunner.MG_Run_Dir,tfile)                    
                            break
                    processes=[]
                    splitdir=os.path.dirname(self.LHEfile)
    
                    if splitdir != '':
                        if not os.path.exists(splitdir):
                            os.makedirs(splitdir)
                        
                        # move events to desired location
                        


                    for np in range(self.n_cores):
                        #processes.append(mp.Process(target=debug.command_line_log_PID,args=('../run.sh %d %d' %(coreevents,seedstart+np),os.path.join(temp_dir,'gpr_%d' %np))))
                        processes.append(mp.Process(target=rungp,args=(np,seedstart)))
                    for p in processes:
                        p.start()
                        time.sleep(0.1)
                    
                    while True:
                        time.sleep(10)
                        allExit = True
                        for t in processes:
                            if t.exitcode is None:
                                allExit = False
                                break

                        if allExit:
                            break
                    for p in processes:
                        p.join()               
                else:
                    lhe_events=rungp(0,seedstart)
                    self.yaml_data['settings']['LHE file'] = lhe_events
                os.chdir(temp_dir)
                ## Now to run
                if os.path.exists(self.HAjson):
                    os.remove(self.HAjson)
                if os.path.exists(self.effname):
                    os.remove(self.effname)
                if os.path.exists(self.cutflowname):
                    os.remove(self.cutflowname)
                if self.histoname is not None:
                    if os.path.exists(self.histoname):
                        os.remove(self.histoname)
                #self.yaml_data['settings']['LHE file'] = lhe_events
                with open(local_yaml,'w') as stream:
                    whatever = yaml.dump(self.yaml_data,stream)

                debug.command_line_log(self.hackanalysis_cmd, log)
        
                if not os.path.exists(self.effname) or not os.path.exists(self.cutflowname):
                    ## point has failed
                    log.error('Failed to run HackAnalysis')
                    raise NameError('Failed to run HackAnalysis')
                
                if os.path.exists(self.HAjson):
                    with open(self.HAjson,'r') as jsonfile:
                        new_json_data=json.load(jsonfile,object_pairs_hook=OrderedDict)
                    if self.combiner is None:
                        output_json_data=new_json_data
                    else:
                        output_json_data=self.combiner.addresults(output_json_data,new_json_data)

                # write merged files so the user can keep track (or in case of crashes!)
            
                if self.batch_size < self.nevents:
                    with open('temp_HA_output.json','w') as fp:
                        json.dump(output_json_data,fp,indent=4)
                    if self.combiner is not None:
                        self.combiner.write_efficiency_file(output_json_data,'temp_HA_efficiency_file.eff')

                if self.check_convergence: # check convergence or exclusion
                    if self.check_exp_convergence:
                        #check_convergence_exclusion(self,hackanalysis_output,nloxs=None,input_max_batch_size=None,safety_ratio=0.25,verbose=False)
                        needed_extra_events,improve_regions,_=self.convergence.check_convergence_exclusion(output_json_data,None,self.batch_size,self.safety_ratio)
                    else:
                        #check_convergence_fixed_goal(self,efficiencydata,rel_uncertainty=0.1,eff_threshold=1e-5,input_max_batch_size=None,verbose=True)
                        #needed_extra_events,improve_regions,_=self.convergence.check_convergence_fixed_goal(output_json_data,1.e-1,1.e-5,self.batch_size)
                        needed_extra_events,improve_regions,_=self.convergence.check_convergence_fixed_goal(output_json_data,self.convergence_error,self.convergence_threshold,self.batch_size)

                    if needed_extra_events > 0:
                        outstr=''
                        log.debug('Estimate need to generate %d more events' %needed_extra_events)
                        for analysis in improve_regions:
                            linestring=','.join(improve_regions[analysis])
                            outstr=outstr+analysis+': '+linestring+'\n'
                        log.debug('Improve regions: \n'+outstr)

                    else:
                        converged=True



            os.chdir(temp_dir)

            # If we have multiple batches, need to merge the outputs into one file,
            # but otherwise we will just keep the original
            if self.batch_size < self.nevents:
                with open(self.HAjson,'w') as jsonfile:
                    json.dump(output_json_data,jsonfile,indent=4)
                if self.combiner is not None:
                    self.combiner.write_efficiency_file(output_json_data,self.effname)

                
            if os.path.exists('temp_HA_output.json'):
                os.remove('temp_HA_output.json')
            if os.path.exists('temp_HA_efficiency_file.eff'):
                os.remove('temp_HA_efficiency_file.eff')

            debug.command_line_log('chmod a+w -R madevent',log)
            shutil.rmtree('madevent')
                
        ####################################################  
        elif self.batch_size < self.nevents: # batch mode but not gridpack
            ## first write the temporary files before first loop, since we already ran HA
            # Let's now harden the tool since on my computer MadGraph has a tendency to randomly crash

            if failed_attempts ==0:
                with open('temp_HA_output.json','w') as fp:
                            json.dump(output_json_data,fp,indent=4)
                if self.combiner is not None:
                    self.combiner.write_efficiency_file(output_json_data,'temp_HA_efficiency_file.eff')
            
            while doneevents < self.nevents and not converged  and failed_attempts < 3: # this is for 3 *consecutive* failures we give up
                # first generate more events
                
                try:
                    MGinfo=self.MGrunner.run(spc_file,temp_dir,log,data_point)
                except Exception as e:
                    log.debug('Problem running MadGraph: '+str(e))
                    failed_attempts+=1
                    continue
                if self.mode == 'LHE':
                    lhe_events=None
                    for tfile in os.listdir(self.MGrunner.MG_Run_Dir):
                        if tfile[-7:] == '.lhe.gz' or tfile[-4:] == '.lhe':
                            lhe_events=os.path.join(self.MGrunner.MG_Run_Dir,tfile)                    
                            break
                    if not os.path.exists(lhe_events):
                        # Allow for failure
                        failed_attempts+=1
                        continue
                        
                        #raise NameError("Failed to generate lhe events!")
                    if self.n_cores > 1:
                        splitlhe(self.n_cores,lhe_events,os.path.join(temp_dir,self.LHEfile))
                    else:
                        self.yaml_data['settings']['LHE file'] = lhe_events
                else: ## HEPMC 
                    hepmc_events=None
                    for tfile in os.listdir(self.MGrunner.MG_Run_Dir):
                        if tfile[-9:] == '.hepmc.gz' or tfile[-6:] == '.hepmc':
                            hepmc_events=os.path.join(self.MGrunner.MG_Run_Dir,tfile)                    
                            break
                    if not os.path.exists(hepmc_events):
                        failed_attempts+=1
                        continue
                        #raise NameError("Failed to generate hepmc events!") 
                    self.yaml_data['settings']['HEPMC file'] = hepmc_events

                if os.path.exists(self.HAjson):
                    os.remove(self.HAjson)
                if os.path.exists(self.effname):
                    os.remove(self.effname)
                if os.path.exists(self.cutflowname):
                    os.remove(self.cutflowname)
                if self.histoname is not None:
                    if os.path.exists(self.histoname):
                        os.remove(self.histoname)
                #self.yaml_data['settings']['LHE file'] = lhe_events
                with open(local_yaml,'w') as stream:
                    whatever = yaml.dump(self.yaml_data,stream)

                attempt_ok = True
                try:
                    debug.command_line_log(self.hackanalysis_cmd, log)
                    if not os.path.exists(self.effname) or not os.path.exists(self.cutflowname):
                        ## point has failed
                        log.error('Failed to run HackAnalysis')
                        attempt_ok=False
                        #raise NameError('Failed to run HackAnalysis')
                
                    if os.path.exists(self.HAjson):
                        with open(self.HAjson,'r') as jsonfile:
                                new_json_data=json.load(jsonfile,object_pairs_hook=OrderedDict)
                        output_json_data=self.combiner.addresults(output_json_data,new_json_data)
                    else:
                        attempt_ok=False
                except Exception as e:
                    log.error('Error running HackAnalysis: '+str(e))
                    log.error('That makes %d consecutive failed attempts' %failed_attempts)
                    attempt_ok=False
        
                
                if attempt_ok: 
                    failed_attempts=0   
                    doneevents += self.batch_size
                    # write merged files so the user can keep track (or in case of crashes!)
            
                    if self.batch_size < self.nevents:
                        with open('temp_HA_output.json','w') as fp:
                            json.dump(output_json_data,fp,indent=4)
                        self.combiner.write_efficiency_file(output_json_data,'temp_HA_efficiency_file.eff')
                else:
                    failed_attempts+=1

                

                if self.check_convergence and attempt_ok: # check convergence or exclusion
                    if self.check_exp_convergence:
                        #check_convergence_exclusion(self,hackanalysis_output,nloxs=None,input_max_batch_size=None,safety_ratio=0.25,verbose=False)
                        needed_extra_events,improve_regions,_=self.convergence.check_convergence_exclusion(output_json_data,None,self.batch_size,self.safety_ratio)
                    else:
                        #check_convergence_fixed_goal(self,efficiencydata,rel_uncertainty=0.1,eff_threshold=1e-5,input_max_batch_size=None,verbose=True)
                        #needed_extra_events,improve_regions,_=self.convergence.check_convergence_fixed_goal(output_json_data,1.e-1,1.e-5,self.batch_size)
                        needed_extra_events,improve_regions,_=self.convergence.check_convergence_fixed_goal(output_json_data,self.convergence_error,self.convergence_threshold,self.batch_size)

                    if needed_extra_events > 0:
                        outstr=''
                        log.debug('Estimate need to generate %d more events' %needed_extra_events)
                        for analysis in improve_regions:
                            linestring=','.join(improve_regions[analysis])
                            outstr=outstr+analysis+': '+linestring+'\n'
                        log.debug('Improve regions: \n'+outstr)

                    else:
                        log.debug('Point has converged!')
                        converged=True

            # finished event loop for non-gridpack mode, dump the final file
            with open(self.HAjson,'w') as jsonfile:
                    json.dump(output_json_data,jsonfile,indent=4)
            self.combiner.write_efficiency_file(output_json_data,self.effname)    
            if os.path.exists('temp_HA_output.json'):
                os.remove('temp_HA_output.json')
            if os.path.exists('temp_HA_efficiency_file.eff'):
                os.remove('temp_HA_efficiency_file.eff')

        ####################################################  
        ## Here we have finished running 

        # Now to append results to the spc file

        HA_spc=zslha.read(self.effname)
        if self.check_exclusion:
            for analysis in output_json_data:
                if analysis not in self.analysishandlers:  # have to check that it was initialised correctly...
                    continue
                
                #getstatsfromdict(self,effdict,inputxs=None,with_significance=True,with_upper_limits=False)
                try:
                    exclusion_results=self.analysishandlers[analysis].getstatsfromdict(output_json_data[analysis],None,False,True)
                except Exception as e:
                    print("Exception %s running exclusion results for analyis %s, giving mu = 10" %(str(e),analysis))
                    exclusion_results={}
                    exclusion_results['mu95exp']=10.0
                    exclusion_results['mu95obs']=10.0
                    exclusion_results['CLs exp']=-1.0
                    exclusion_results['CLs obs']=-1.0
                    #continue
                print("============ %s ============== " %analysis)
                print(exclusion_results)
                HA_spc.analyses[analysis.upper()]['PROCESS']['MU95EXP']=float(exclusion_results['mu95exp'])
                HA_spc.analyses[analysis.upper()]['PROCESS']['MU95OBS']=float(exclusion_results['mu95obs'])
                #HA_spc.blockcomments['PROCESS'][analysis.upper()]['MU95EXP']=analysis+' (expected)'
                #HA_spc.blockcomments['PROCESS'][analysis.upper()]['MU95OBS']=analysis+' (observed)'
                if float(exclusion_results['mu95obs']) < 1:
                    HA_spc.analyses[analysis.upper()]['PROCESS']['EXCLUDED']=0
                    #HA_spc.blockcomments['PROCESS'][analysis.upper()]['EXCLUDED']='Point Excluded (0)'
                else:
                    HA_spc.analyses[analysis.upper()]['PROCESS']['EXCLUDED']=1

                HA_spc.analyses[analysis.upper()]['PROCESS']['CLSEXP']=float(exclusion_results['CLs exp'])
                HA_spc.analyses[analysis.upper()]['PROCESS']['CLSOBS']=float(exclusion_results['CLs obs'])
                
                if 'bayes' in exclusion_results:
                    HA_spc.analyses[analysis.upper()]['PROCESS']['BAYES']=float(exclusion_results['bayes'])
                    #HA_spc.blockcomments['PROCESS'][analysis.upper()]['EXCLUDED']='Point Allowed (1)'
            ## overwrite the efficiency file so we can append it to spc
            #HA_spc.write(self.effname)
            """
            mubest=10.0
            muexpbest=10.0
            #analysisresults=[]
            newblock={}
            newblockcomments={}

            newindex=10
            for analysis in output_json_data:
                #getstatsfromdict(self,effdict,inputxs=None,with_significance=True,with_upper_limits=False)
                exclusion_results=self.analysishandlers[analysis].getstatsfromdict(output_json_data[analysis],None,False,True)
                print("============ %s ============== " %analysis)
                print(exclusion_results)
                newblock[str(newindex)]=exclusion_results['mu95exp']
                newblock[str(newindex+1)]=exclusion_results['mu95obs']
                newblockcomments[str(newindex)]=analysis+' (expected)'
                newblockcomments[str(newindex+1)]=analysis+' (observed)'
                #analysisresults.append([analysis.upper(),exclusion_results['mu95exp'],exclusion_results['mu95obs']])
                mubest=min(exclusion_results['mu95obs'],mubest)
                muexpbest=min(exclusion_results['mu95exp'],muexpbest)
                newindex=newindex+2
            if mubest < 1:
                newblock['1']=0
                newblockcomments['1']='Point Excluded (0)'
            else:
                newblock['1']=1
                newblockcomments['1']='Point Allowed (0)'
            newblock['2']=mubest
            newblockcomments['2']='Observed limit on signal strength'
            newblock['3']=muexpbest
            newblockcomments['3']='Expected limit on signal strength'
            data_point.spc.blocks['HACKANALYSIS '+self.process]=
            """
        if data_point.spc is None:
            data_point.spc=HA_spc
        else:
            data_point.spc.add(HA_spc)
        
       
        with open(spc_file,'a') as FO:
            for analysis in output_json_data:
                upanalysis=analysis.upper()
                FO.write('BLOCK PROCESS %s\n' % upanalysis)
                FO.write('XS    %.4e   # (fb)\n' % HA_spc.analyses[upanalysis]['PROCESS']['XS'])
                FO.write('Events   %d   \n' % HA_spc.analyses[upanalysis]['PROCESS']['Events'])
                if self.check_exclusion and 'MU95EXP' in HA_spc.analyses[upanalysis]['PROCESS']:
                    FO.write('MU95EXP    %.4e   # \n' % HA_spc.analyses[upanalysis]['PROCESS']['MU95EXP'])
                    FO.write('MU95OBS    %.4e   # \n' % HA_spc.analyses[upanalysis]['PROCESS']['MU95OBS'])
                    FO.write('CLSEXP    %.5f   # \n' % HA_spc.analyses[upanalysis]['PROCESS']['CLSEXP'])
                    FO.write('CLSOBS    %.5f   # \n' % HA_spc.analyses[upanalysis]['PROCESS']['CLSOBS'])

                    excl=int(HA_spc.analyses[upanalysis]['PROCESS']['EXCLUDED'])

                    if excl == 0:
                        FO.write('EXCLUDED    %d   # Point excluded (0) \n' %excl)
                    else:
                        FO.write('EXCLUDED    %d   # Point allowed (1) \n' %excl)
                        
                if 'BAYES' in HA_spc.analyses[upanalysis]['PROCESS']:
                    FO.write('BAYES    %.4e   # \n' % HA_spc.analyses[upanalysis]['PROCESS']['BAYES'])

                FO.write('BLOCK EFFICIENCIES %s # name, eff, rel. uncert. [systup   systdown]\n' % analysis)
                for region in output_json_data[analysis]["Regions"]:
                    if "sysup" in output_json_data[analysis]["Regions"][region]:
                        FO.write('%s   %.8e   %.8e   %.8e   %.8e\n' %(region,output_json_data[analysis]["Regions"][region]["efficiency"],output_json_data[analysis]["Regions"][region]["stat"]*0.01,output_json_data[analysis]["Regions"][region]["sysup"]*0.01,output_json_data[analysis]["Regions"][region]["sysdown"]*0.01))
                    else:
                        FO.write('%s   %.8e   %.8e\n' %(region,output_json_data[analysis]["Regions"][region]["efficiency"],output_json_data[analysis]["Regions"][region]["stat"]*0.01))
        #else:
            #debug.command_line_log("cat " + self.effname + " >> " + os.path.join(spc_file), log)   

            #for b in HA_spc.blocks:
            #    zslha.write_block_head(b,SF)

        #debug.command_line_log("cat " + self.effname + " >> " + os.path.join(spc_file), log)
        # Now do we store the cutflows and yoda files too? These can be picked up by copying everything

        ###Clean up the split directory -- otherwise it will be copied too!!!
        if not self.store_events and self.mode == 'LHE' and (self.n_cores > 1 or self.useGridPack):
            if lhe_events[-3:]=='.gz':
                ending = '.lhe.gz'
            else:
                ending = '.lhe'
            for x in range(self.n_cores):
                try:
                    os.remove(self.LHEfile+'_'+str(x)+ending)
                except Exception as e:
                    log.warning('Failed to remove split file, '+str(e))

