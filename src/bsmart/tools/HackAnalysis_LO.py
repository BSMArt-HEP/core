"""
Tool to run HackAnalysis with event generation handled internally by pythia.


Need to specify in tools:

.. code-block:: text

    "HackAnalysis_LO" : {
        "Command" : "/path/to/executable", (i.e. analysePYTHIA.exe)
        "YAML file" "<NAME>.yaml"    (file for configuring HackAnalysis run  )
    }

Note that the pythia configuration file must also be present in the launch directory, and you must
specify it in the YAML file.


"""


__meta__ = {
    "name": "HackAnalysis_LO",
    "requires": ["pyyaml"],
    "settings": {
        "Command": "Path to executable",
        "YAML file": "Configuration YAML file"
    }
}

import os
import shutil

from bsmart import debug
from bsmart.HEPRun import HepTool, DataPoint

from bsmart.collider.pythiacfg import PYTHIACFG

import fnmatch
import sys
import math
from bsmart import zslha 

import yaml

class NewTool(HepTool):
    """ Runs madgraph and extracts the cross-section plus uncertainty """
    def __init__(self, name, settings,global_settings=None):
        HepTool.__init__(self, name, settings,global_settings)

        

        
        if 'YAML file' not in self.settings:
            #self.log.error('No YAML file present!')
            raise NameError('No YAML file present!')

        #if 'Store Cutflows' in self.settings:
        #    self.store_cutflows=eval(self.settings['Store Cutflows'])

        #if self.store_cutflows:
            

        if global_settings is not None and global_settings.varnames is not None:
            self.varnames=global_settings.varnames
        else:
            self.varnames=None
        
        #tpath=os.path.split(self.settings['YAML file'])
        tpath=os.path.dirname(self.settings['YAML file'])
        ## determine if absolute path is given, or if the YAML file is in the launch directory
        
        #if len(tpath) == 1:
        if tpath == '':
            self.yamlname=self.settings['YAML file']
            #self.yaml = os.path.join(self.settings['cwd'],self.yamlname)
            #yamldir=self.settings['cwd']
            yamldir=os.getcwd()
            self.yaml = os.path.join(yamldir,self.yamlname)
            
        else:
            self.yaml = self.settings['YAML file']
            self.yamlname=os.path.split(self.yaml)[-1]
            yamldir=os.path.split(os.path.dirname(self.yaml))[0]
        #print(self.yaml)
        #### READ OTHER INFO FROM THE YAML!!!!
        try:
            with open(self.yaml, 'r') as stream:        
                self.yaml_data = yaml.safe_load(stream)
        except:
            raise NameError('Failed to open HackAnalysis yaml file!')

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
        #if 'Config file' not in self.yaml_data['settings']:
        #    raise NameError('No cfg file for pythia present!')

        #self.cfgname=self.yaml_data['settings']['Config file']
        #self.cfgfile = os.path.join(self.settings['cwd'],self.cfgname)
        #self.cfgfile =os.path.join(yamldir,self.cfgname)
        #print(self.cfgfile)
        
        """
        if 'cfg file' not in self.settings:
            #log.error('No cfg file for pythia present!')
            raise NameError('No cfg file for pythia present!')
        tpath=os.path.split(self.settings['cfg file'])
        ## determine if absolute path is given, or if the YAML file is in the launch directory
        if len(tpath) == 0:
            self.cfgname=self.settings['YAML file']
            self.cfgfile = os.path.join(self.settings['cwd'],self.cfgname)
            
        else:
            self.cfgfile = self.settings['cfg file']
            self.cfgname=os.path.split(self.cfgfile)[-1]
        """
        self.command = self.settings['Command']+' '+self.yamlname

        if 'Cores' in self.settings:
            self.n_cores=int(self.settings['Cores'])
            self.yaml_data['settings']['cores'] = self.n_cores
        else:
            if 'cores' in self.yaml_data['settings']:
                self.n_cores = int(self.yaml_data['settings']['cores'])
            else:
                self.n_cores = 1
                self.yaml_data['settings']['cores'] = 1

        if 'Efficiency Filename' in self.yaml_data['settings']:
            self.effname = self.yaml_data['settings']['Efficiency Filename']
        else:
            self.effname = 'test.slha'

        if 'Cutflow Filename' in self.yaml_data['settings']:
            self.cutflowname = self.yaml_data['settings']['Cutflow Filename']
        else:
            self.cutflowname = 'test.eff'
        
        if 'Histogram Filename' in self.yaml_data['settings']:
            self.histoname = self.yaml_data['settings']['Histogram Filename']
        else:
            self.histoname = 'test.yoda'
            
    def run(self, spc_file, temp_dir, log,data_point):
        local_yaml = os.path.join(temp_dir, self.yamlname)

        var_dict = data_point.get_var_dict(self.varnames)

        with open(local_yaml,'w') as stream:
            whatever = yaml.dump(self.yaml_data,stream)

        local_cfg = os.path.join(temp_dir, self.cfgname)
        #self.pythia_cfg.write(local_cfg,data_point.inputs,self.varnames)
        self.pythia_cfg.write_dict(local_cfg,var_dict)
        try:
            ## also append the SLHA to the cfg file
            OF = open(local_cfg,"a")
            OF.write('SLHA:file ='+spc_file+'\n')
            OF.close()
        except Exception as e:
            raise NameError('Failed to edit  pythia cfg file!' + str(e))

        """
        if not os.path.exists(local_yaml):
            try:
                shutil.copyfile(self.yaml,local_yaml)
            except:
                log.error('Failed to copy YAML file! '+self.yaml)
                raise
        local_cfg = os.path.join(temp_dir, self.cfgname)
        if not os.path.exists(local_cfg):
            try:
                shutil.copyfile(self.cfgfile,local_cfg)
                ## also append the SLHA to the cfg file?
                OF = open(local_cfg,"a")
                OF.write('SLHA:file ='+spc_file+'\n')
                OF.close()
            except:
                log.error('Failed to copy cfg file!')
                raise
                #debug.command_line_log('echo "SLHA:file =' + spc_file + '" >> '+local_cfg,log)
        """
        ## Now to run
        if os.path.exists(self.effname):
            os.remove(self.effname)
        if os.path.exists(self.cutflowname):
            os.remove(self.cutflowname)
        if self.histoname is not None:
            if os.path.exists(self.histoname):
                os.remove(self.histoname)
        debug.command_line_log(self.command, log)
        
        if not os.path.exists(self.effname) or not os.path.exists(self.cutflowname):
            ## point has failed
            log.error('Failed to run HackAnalysis')
            raise NameError('Failed to run HackAnalysis')

        HA_spc=zslha.read(self.effname)
        if data_point.spc is None:
            data_point.spc=HA_spc
        else:
            data_point.spc.add(HA_spc)

        # Now to append results to the spc file
        debug.command_line_log("cat " + self.effname + " >> " + os.path.join(spc_file), log)

        # Now do we store the cutflows and yoda files too? These can be picked up by copying everything
