"""
Tool to run SModelS version 2 (kept as legacy)



"""

__meta__ = {
    "name": "SModelS2",
    "requires": ["smodels", "unum", "pyslha"],
    "settings": {
        "Path": "Path to SModelS",
        "Cache": "Path to cache",
        "Parameter file": "Path to parameters",
        "QNUMBERS file": "Path to QNUMBERS"
    }
}


import os
import sys
import shutil
#import subprocess
from bsmart import debug
import importlib
from bsmart.HEPRun import HepTool, DataPoint
try:
    from imp import reload
except:
    try:
        from importlib import reload
    except Exception as e:
        raise NameError("Failed to import SModelS reload from importlib: "+str(e))

from bsmart import zslha

# def noshell_command_log_withpid(command, wdir, log):
#     log.debug('Terminal command: %s' % command)
#     command_line = subprocess.Popen(command, shell=False, cwd=wdir,
#                                     stdout=subprocess.PIPE,
#                                     stderr=subprocess.PIPE)
#     pid=command_line.pid
#     out, err = command_line.communicate()
#     # log.debug(out)
#     if err != b'':
#         log.error(err.decode("utf-8"))
#     return str(pid)


class NewTool(HepTool):
    """ overload the init to initialise the DM limit calculator, but name and settings are already given in HepTool """
    def __init__(self, name, settings,global_settings=None):
        HepTool.__init__(self, name, settings,global_settings)

        if 'Path' in self.settings:
            #self.SModelS_path=self.settings['Path']
            #sys.path.insert(0,self.SModelS_path)
            sys.path.insert(0,self.settings['Path'])

        if 'Cache' in self.settings:
            os.environ['SMODELS_CACHEDIR']=self.settings['Cache']

        #import logging
        #logging.basicConfig(level=logging.CRITICAL)
        try:
            
            #from smodels.installation import installDirectory, version
            
            #from smodels.tools import modelTester
            self.modelTester = importlib.import_module('smodels.tools.modelTester')
            #from smodels.tools import crashReport
            #self.crashReport=importlib.import_module('smodels.tools.crashReport')
            #from smodels.tools import smodelsLogging
            self.smodelsLogging=importlib.import_module('smodels.tools.smodelsLogging')
            
            
            #from smodels.tools import runtime
            self.runtime=importlib.import_module('smodels.tools.runtime')
            #from smodels import particlesLoader
            self.runtime.modelFile=self.settings['QNUMBERS file']
            self.particlesLoader= importlib.import_module('smodels.particlesLoader')
            
            
        except Exception as e:
            raise NameError('Cannot import necessary SModelS packages, error '+str(e))
        
        self.smodelsLogging.setLogLevel('error')
        """ Read and check parameter file, exit parameterFile does not exist """
        #self.parser = modelTester.getParameters(self.settings['Parameter file'])
        self.parameterFile=self.settings['Parameter file']
        self.parser = self.modelTester.getParameters(self.parameterFile)

        """ Determine particles module from ini file, if necessary """
        #if parser.has_option("particles","model"):
        #    runtime.modelFile = parser.get( "particles", "model" )
        #    reload(particlesLoader)
        self.runtime.modelFile=self.settings['QNUMBERS file']
        reload(self.particlesLoader)
        db=None
        """ Check database location and load database, exit if not found """
        self.database, self.databaseVersion = self.modelTester.loadDatabase(self.parser, db)

        """ Get list of input files to be tested """
        #fileList, inDir = modelTester.getAllInputFiles(inFile)

        """ Create output directory if missing """
        #if not os.path.isdir(outputDir): os.mkdir(outputDir)

        """ Load analysis database, print list """
        self.listOfExpRes = self.modelTester.loadDatabaseResults(self.parser, self.database)

        self.timeout=120 ## allow one minute to run ...
        self.development=False
        #self.smodelsLogging.setLogLevel('error')
        
    ## Add path to smodels here

    ## import database



    
        
    def run(self, spc_file, temp_dir, log,data_point):

        """ Test all input points """
        #modelTester.testPoints(fileList, inDir, outputDir, parser, databaseVersion,
        #         listOfExpRes, timeout, development, parameterFile )
        #fileList, inDir = modelTester.getAllInputFiles(spc_file)
        #modelTester.testPoints(fileList, inDir, temp_dir, self.parser, self.databaseVersion,
        #         self.listOfExpRes, self.timeout, self.development, sefl.parameterFile )
        #runSingleFile(inputFile, outputDir, parser, databaseVersion, listOfExpRes,
        #          timeout, development, parameterFile):
        try:
            #self.smodelsLogging.setLogLevel('error')
            smodels_results=self.modelTester.runSingleFile(spc_file, temp_dir, self.parser, self.databaseVersion, self.listOfExpRes,
                  self.timeout, self.development, self.parameterFile)
        except Exception as e:
            raise NameError('SModelS exception '+str(e))
        #import logging
        #logging.basicConfig(level=logging.CRITICAL)
        rmax=0
        best_analysis='None'
        
        
        #for rrr in smodels_results[spc_file]['python']:
        #    print('Smodels result: '+str(rrr))
        #    #print(smodels_results[rrr])
        if spc_file in smodels_results and 'python' in smodels_results[spc_file] and 'ExptRes' in smodels_results[spc_file]['python']:
            for rrr in smodels_results[spc_file]['python']['ExptRes']:
                #print(rrr['AnalysisID']+str(': ')+str(rrr['r']))
                if 'r' in rrr:
                    if rrr['r'] > rmax:
                        rmax=rrr['r']
                        best_analysis=rrr['AnalysisID']
        
        """ now append to the spectrum file """
        with open(spc_file,"a") as OF:

            OF.write('BLOCK SMODELS # Best analysis:'+best_analysis+' \n')
            OF.write(' 1 %10.6E # Best R ratio (observed)\n' % float(rmax))
            OF.write(' 2 %10.6E # r - 1: > 0 if excluded, < 0 if allowed\n' % (rmax-1.0))      
            

        if data_point.spc is None:
            data_point.spc=zslha.read(spc_file)
        else:
            #if data_point.spc is not None:
            data_point.spc.blocks['SMODELS']={'1': float(rmax), '2':float(rmax-1.0)}
            data_point.spc.blockcomments['SModelS']={'1': 'Best R ratio (observed)','2': 'r - 1: > 0 if excluded, < 0 if allowed'}