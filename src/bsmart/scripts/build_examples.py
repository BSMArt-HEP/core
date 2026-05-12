#!/usr/bin/env python3

"""

Script to build all runnable examples!

"""


import os,sys,json
from bsmart.scripts import prepare_model
from bsmart.scripts import prepare_bsmart as prepareBSMArt
import argparse
import subprocess
import shutil
import collections
import importlib.resources

def havemodel(ModelName,extra_commands_list,run_if_missing=True):
    # extra_commands is now a list of args
    modelname=str(ModelName).replace('/','-')
    cwd = os.getcwd()
    scriptdir=os.path.join(cwd,'ScriptsAndLogs')
    bsmartdatafile=os.path.join(scriptdir,'BSMArt_data_'+modelname+'.json')

    if not os.path.exists(bsmartdatafile) and run_if_missing:
        # command_to_execute='python3 prepareModel.py '+ModelName+' '+extra_commands
        # call prepare_model.main
        
        args = ['prepare_model.py', ModelName] + extra_commands_list
        sys.argv = args
        try:
            prepare_model.main()
        except SystemExit as e:
             if e.code is not None and e.code != 0:
                  print("Error running prepare_model for "+ModelName)
                  return None
        except Exception as e:
             print("Error running prepare_model for "+ModelName + ": "+str(e))
             return None

    if os.path.exists(bsmartdatafile):
        try:
            with open(bsmartdatafile,'r') as bsmdf:
                bsmart_data=json.load(bsmdf)
            return bsmart_data
        except:
            return None
    else:
        return None

def main():

    try:
        parser = argparse.ArgumentParser(
            description='Build model examples. Either allow all examples to be built (default), or select for a specific model.')
        parser.add_argument("--Model", help="Build for only this model/these models",nargs="+",action="store",
                            default=None)
        # Are --short and --csv really necessary any more? 
        parser.add_argument("--NoSARAH", help="Don't run SARAH: only build examples for the cases where the data already exists or is not needed.",
                            action="store_true")
        args = parser.parse_args()
    except:
        # Don't need a default!
        raise SystemExit

    ## need this to locate the BSMArt install path!
    ## need this to locate the BSMArt install path!
    jsonFILE='HEPtoolpaths.json'
    if not os.path.exists(jsonFILE):
         try:
             from bsmart import get_heptools_path_file
             jsonFILE = get_heptools_path_file()
         except ImportError:
             pass

    if os.path.exists(jsonFILE):
        try:
            with open(jsonFILE) as json_data:
                paths_dict = json.load(json_data)
        except Exception as e:
            #log.error('Failed to load json file '+file)
            print('Failed to load json file '+jsonFILE)
            print('Json exception given as ' +str(e))
            raise SystemExit
    else:
        print('Failed to load json file '+jsonFILE)
        print("Please run BSMArt-InstallHEPTools or ensure HEPtoolpaths.json is present in the current directory or configuration folder.")
        raise SystemExit

    # if 'BSMArt' not in paths_dict:
    #     print('Could not locate BSMArt installation, exiting')
    #     raise SystemExit

    # examples_dir=os.path.join(paths_dict['BSMArt'],'InputExamples')
    
    # Use resources for examples directory
    # We need to copy examples out? 
    # prepare_bsmart.updateBSMArtTemplate takes a template FILE path.
    # So we need to access the file path of the json example.
    
    # If using importlib.resources.files, specific files can be accessed as path.
    # Directory access is also possible.
    
    examples_resource = None
    try:
         from importlib.resources import files
         examples_resource = files('bsmart.data').joinpath('InputExamples')
    except:
         print("Could not find InputExamples in package!")
         sys.exit(1)
         
    
    
    cwd = os.getcwd()
    scriptdir=os.path.join(cwd,'ScriptsAndLogs')



    allcodes={ 'SPheno':True, 'FlexibleSUSY': False, 'MicrOMEGAs':True, 'Vevacious++':True, 'HiggsSignals':True, 'HiggsBounds':True, 'HiggsTools':True, 'SModelS': True}
    ## Do MSSM ones first


    forcebuild = not args.NoSARAH
    #models=['MSSM','SMSSM','DiracGauginos','E6AFH','SM-SQQ','NoSARAH']
    # E6AFH is too slow ...
    models=['MSSM','SMSSM','DiracGauginos','SM-SQQ','NoSARAH']


    if args.Model is not None:
        models = args.Model

    def get_example_path(subpath_list):
         # subpath_list i.e. ['Minimise', 'Fitting_MSSM_mh.json']
         # We need to extract it to a temp file or get its path if it is on disk.
         # For simplicity, if we are installed as a directory, we can use the path.
         # If zip, we have to extract.
         # importlib.resources.as_file context manager gives a path.
         p = examples_resource
         for s in subpath_list:
              p = p.joinpath(s)
         
         # We have to keep the file alive if it is temp?
         # as_file returns a path we can use inside the block.
         # But we iterate over multiple files.
         # We will assume installed as standard package directory (since user is on linux and standard pip install usually unzips or we can assume it works for now, or copy to temp).
         # But as_file is safer.
         return p


    def update_template(model, data, paths, codes, subpath_list):
         p = get_example_path(subpath_list)
         with importlib.resources.as_file(p) as fpath:
              prepareBSMArt.updateBSMArtTemplate(model, data, paths, codes, str(fpath))
    
    def copy_file(subpath_list, dest):
         p = get_example_path(subpath_list)
         with importlib.resources.as_file(p) as fpath:
              shutil.copyfile(fpath, dest)


    if 'MSSM' in models:
        #mssmdata=havemodel('MSSM','--All',forcebuild)
        mssmdata=havemodel('MSSM',['--debug', '--SPheno', '--BSMArt', '--MicrOMEGAs', '--HiggsTools', '--HiggsBounds', '--HiggsSignals', '--VevaciousPlusPlus'],forcebuild)
        if mssmdata is not None:
            # Fitting
        
            """
            with open(os.path.join(examples_dir,'Minimise',template_name),'r') as IF:
                    template_dict = json.load(IF,object_pairs_hook=collections.OrderedDict)
            prepareBSMArt.updateBSMArtTemplate('MSSM',mssmdata,paths_dict,includecodes,template_dict,template_name)
            """
            update_template('MSSM',mssmdata,paths_dict,allcodes,['Minimise','Fitting_MSSM_mh.json'])
            update_template('MSSM',mssmdata,paths_dict,allcodes,['ValidityCheck','validity_MSSM.json'])
            update_template('MSSM',mssmdata,paths_dict,allcodes,['MLS','MLS_MSSM.json'])

    if 'SMSSM' in models:
        smssmdata=havemodel('SMSSM',['--SPheno', '--BSMArt'],forcebuild) 
        
        if smssmdata is not None:
            update_template('SMSSM',smssmdata,paths_dict,allcodes,['SMSSM','SMSSM_Grid.json'])

    if 'DiracGauginos' in models:
        #mdgssmdata=havemodel('DiracGauginos','--All',forcebuild)
        mdgssmdata=havemodel('DiracGauginos',['--debug', '--SPheno', '--BSMArt', '--MicrOMEGAs', '--HiggsTools', '--HiggsBounds', '--HiggsSignals', '--VevaciousPlusPlus'],forcebuild)

        if mdgssmdata is not None:
            update_template('DiracGauginos',mdgssmdata,paths_dict,allcodes,['DiracGauginos','DiracGauginos_MCMC.json'])
            update_template('DiracGauginos',mdgssmdata,paths_dict,allcodes,['DiracGauginos','DiracGauginos_Random.json'])
            update_template('DiracGauginos',mdgssmdata,paths_dict,allcodes,['DiracGauginos','DiracGauginos_AL.json'])





    if 'E6AFH' in models:
        E6AFH_data=havemodel('E6AFH',['--debug', '--SPheno', '--BSMArt', '--MicrOMEGAs', '--HiggsTools', '--HiggsBounds', '--HiggsSignals', '--VevaciousPlusPlus'],forcebuild) 
        if E6AFH_data is not None:
            update_template('E6AFH',E6AFH_data,paths_dict,allcodes,['E6AFH','CSV_High_E6AFH.json'])
            
            dest = os.path.join(os.getcwd(),'BSMArt_E6AFH','E6AFH_points.csv')
            copy_file(['E6AFH','E6AFH_points.csv'], dest)
            # shutil.copyfile(os.path.join(examples_dir,'E6AFH','E6AFH_points.csv'),os.path.join(os.getcwd(),'BSMArt_E6AFH','E6AFH_points.csv'))

    if 'SMSQQ' in models:
        smsqqdata=havemodel('SM-SQQ',['--debug', '--SPheno', '--BSMArt', '--MicrOMEGAs', '--HiggsTools', '--HiggsBounds', '--HiggsSignals', '--VevaciousPlusPlus'],forcebuild) 

        if smsqqdata is not None:
            update_template('SM-SQQ',smsqqdata,paths_dict,allcodes,['SMSQQ','SMSQQ_MCMC.json'])
            update_template('SM-SQQ',smsqqdata,paths_dict,allcodes,['SMSQQ','SMSQQ_Random.json'])
            update_template('SM-SQQ',smsqqdata,paths_dict,allcodes,['SMSQQ','SMSQQ_Grid.json'])
            update_template('SM-SQQ',smsqqdata,paths_dict,allcodes,['SMSQQ','SMSQQ_AL.json'])

    if 'NoSARAH' in models: 
        """
            Models that don't require SARAH. 
            First let's do the fake MLS scan
        """
        fakeMLSdir=os.path.join(os.getcwd(),'BSMArt_FakeMLS')
        if not os.path.exists(fakeMLSdir):
            os.mkdir(fakeMLSdir)
        copy_file(['MLS','FakeMLS.json'],os.path.join(fakeMLSdir,'FakeMLS.json'))
        copy_file(['MLS','MLS_LesHouches.in.MSSM'],os.path.join(fakeMLSdir,'MLS_LesHouches.in.MSSM'))
        
        # shutil.copyfile(os.path.join(examples_dir,'MLS','FakeMLS.json'),os.path.join(fakeMLSdir,'FakeMLS.json'))
        # shutil.copyfile(os.path.join(examples_dir,'MLS','MLS_LesHouches.in.MSSM'),os.path.join(fakeMLSdir,'MLS_LesHouches.in.MSSM'))
        fakeTooldir=os.path.join(fakeMLSdir,'Tools')
        if not os.path.exists(fakeTooldir):
            os.mkdir(fakeTooldir)
        copy_file(['MLS','Tools','fakefunc.py'],os.path.join(fakeTooldir,'fakefunc.py'))
        # shutil.copyfile(os.path.join(examples_dir,'MLS','Tools','fakefunc.py'),os.path.join(fakeTooldir,'fakefunc.py'))
        
        # BSMArtexec=os.path.join(paths_dict['BSMArt'],'bin','BSMArt')
        # Use simple command
        BSMArtexec = 'BSMArt'
        
        shellscript=os.path.join(fakeMLSdir,'runBSMArt_FakeMLS.sh')
        OF3=open(shellscript,'w')
        OF3.write('#!/usr/bin/env bash\n\n')
        OF3.write(BSMArtexec + ' --debug FakeMLS.json \n')
        OF3.close()
        os.system('chmod u+x '+shellscript)


if __name__ == "__main__":
    main()
