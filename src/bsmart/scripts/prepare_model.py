#!/usr/bin/env python3

"""
Script to prepare a model using SARAH for use with selected tools, including SPheno, MicrOmegas, HiggsBounds, HiggsSignals, HiggsTools, flavio, MultiNest and BSMArt

2025/01/29
"""


import os
from datetime import datetime
from socketserver import ThreadingUnixStreamServer
import sys

import shutil
import json
import collections

import argparse

import subprocess

import logging

try:
    from bsmart.scripts import prepare_bsmart as prepareBSMArt
except ImportError:
    try:
        import prepareBSMArt
    except:
        print('Could not import prepareBSMArt, necessary for creating the BSMArt template!')
        sys.exit(1)


log = None

def get_lib_path(basedir):
    dirs_to_try=[os.path.join(basedir,'lib64'),os.path.join(basedir,'lib')]

    pylibdir=''
    for trydir in dirs_to_try:
        try:    
            if os.path.exists(trydir):
                return trydir
        except: 
            pass
    
    raise NameError('No library directory?')    



def shell_command(command_to_execute,working_dir,errors_are_fatal=False,log_everything=False):
    global log

    try:

        command_line = subprocess.Popen(command_to_execute,shell=True, cwd=working_dir,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        out, err = command_line.communicate()
    except Exception as e:
        raise NameError(e)
    res = out.decode().strip()
    errmess=err.decode().strip()
    returnval=command_line.returncode
    #log.debug(out.decode("utf-8").strip())

    if returnval > 0: ## failed to run at all
        if err !=b'':
            log.error('Error message(s) from '+command_to_execute+': '+errmess)
        if errors_are_fatal:
            raise NameError(command_to_execute + ' failed: '+errmess)
    else:
        if err != b'':
            log.debug('Warning messages from '+command_to_execute+': '+errmess)
                
    if log_everything:
        log.info('Command '+str(command_to_execute)+', produced: \n'+str(res))
    else:
        log.debug('Command '+str(command_to_execute)+', produced: \n'+str(res))
    return res

def main():
    global log
    info_message='Please specify the tool(s) you want to use. SPheno will always be included.\n Valid options are: All, MicrOmegas, HiggsBounds, HiggsSignals, BSMArt.'


    parser = argparse.ArgumentParser(
            description=info_message)
    parser.add_argument('ModelName',
                        metavar='<Model Name>', type=str,
                        help='Input Model name')
    parser.add_argument("--All", help="Setup all tools",
                        action="store_true")
    parser.add_argument("--HiggsBounds", help="Configure for HiggsBounds",
                        action="store_true")
    parser.add_argument("--HiggsSignals", help="Configure for HiggsSignals",
                        action="store_true")
    parser.add_argument("--HiggsTools", help="Configure for HiggsTools",
                        action="store_true")
    parser.add_argument("--MicrOMEGAs", help="Configure for MicrOMEGAs",
                        action="store_true")

    parser.add_argument("--BSMArt", help="Configure for BSMArt",
                        action="store_true")
    parser.add_argument("--flavio", help="Configure for flavio",
                        action="store_true")
    parser.add_argument("--VevaciousPlusPlus", help="Configure Vevacious++",
                        action="store_true")
    parser.add_argument("--UFO", help="Create UFO file",
                        action="store_true")                    
    parser.add_argument("--anyBSM",help="Configure for anyBSM", action="store_true")
    parser.add_argument("--FlexibleSUSY",help="Model name for FlexibleSUSY",action="store",default='',type=str)
    parser.add_argument("--SPheno",help="Skip SPheno",
                        action="store_true",default=False)
    parser.add_argument("--Options", help="Options to pass to MakeSPheno",action="store",
                        default='',type=str)
    parser.add_argument("--MOOptions", help="Options to pass to MakeCHep",action="store",
                        default='',type=str)
    parser.add_argument("--Init", help="Extra commands for the Mathematica script, to be executed before MakeSPheno",action="store",
                        default='',type=str)
    parser.add_argument("--SkipFlavorKit",help="Skip FlavorKit",
                        action="store_true",default=False)
    parser.add_argument("--debug",help="Store extra debug info",
                        action="store_true",default=False)
    parser.add_argument("--math", help="Mathematica executable name (to run on command line)",nargs="?",action="store",
                        default='math')
    parser.add_argument('--local',help='Store data in the current directory/subdirectories only. Otherwise, will store data in the venv/user directory',action='store_true',default=False)


    try:
        args = parser.parse_args()
    except Exception as e:
        print(str(e))
        parser.print_help()
        #print(info_message +str(e))
        raise SystemExit

    model=args.ModelName
    modelname=str(args.ModelName).replace('/','-')

    mathexec = args.math

    runFlexibleSUSY=False
    if args.FlexibleSUSY != '':  # FlexibleSUSY may have a different name to the SARAH one
            runFlexibleSUSY=True

    if args.All:
        args.MicrOMEGAs = True
        args.HiggsBounds = True
        args.HiggsSignals = True
        args.HiggsTools = True
        args.VevaciousPlusPlus = True
        args.BSMArt = True
        args.UFO = True
        args.flavio = True
        args.anyBSM = True
        #args.FlexibleSUSY = True
        runFlexibleSUSY=True
        args.SPheno = True
        if args.FlexibleSUSY == '':  # FlexibleSUSY may have a different name to the SARAH one
            args.FlexibleSUSY = modelname
        

    ## Set up directories
    cwd = os.getcwd()
    scriptdir=os.path.join(cwd,'ScriptsAndLogs')

    if not os.path.exists(scriptdir):
        os.makedirs(scriptdir)

    
        
    ## Create the mathscript
    scriptname=os.path.join(scriptdir,'init'+modelname+'.m')

        
    """ Set up logging """
    log = logging.getLogger(model)
    logfile=os.path.join(scriptdir,'prep_'+modelname+'.log')
    if os.path.exists(logfile):
        os.remove(logfile)
    if args.debug:
        log.setLevel(logging.DEBUG)
    else:
        log.setLevel(logging.INFO)
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    fh = logging.FileHandler(logfile)
    if args.debug:
        fh.setLevel(logging.DEBUG)
    else:
        fh.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    log.addHandler(fh)
    log.addHandler(ch)

    allcmdarguments=' '.join(sys.argv[1:])
    log.info('Invoked with arguments: '+allcmdarguments)

    ## Load the paths
    ## Load the paths
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
                paths_dict = json.load(json_data, object_pairs_hook=collections.OrderedDict)
        except Exception as e:
            #log.error('Failed to load json file '+file)
            log.error('Failed to load json file '+jsonFILE)
            log.error('Json exception given as ' +str(e))
            raise SystemExit
    else:
        log.error('Failed to load json file '+jsonFILE)
        log.error('Please run BSMArt-InstallHEPTools first or ensure HEPtoolpaths.json is in the current directory or configuration folder.')
        raise SystemExit

    ## now check for paths
    if 'SARAH' not in paths_dict:
        log.error('Could not find SARAH path!')
        raise SystemExit

    if args.SPheno or args.All:
        if 'SPheno' not in paths_dict:
            log.warning('Could not find SPheno path!')
            args.SPheno = False    

    if args.FlexibleSUSY != '':
        if 'FlexibleSUSY' not in paths_dict:
            log.warning('Could not find FlexibleSUSY path!')
            args.FlexibleSUSY= '' 
            runFlexibleSUSY=False

    if args.MicrOMEGAs or args.All:
        if 'MicrOMEGAs' not in paths_dict:
            log.warning('Could not find MicrOMEGAs path!')
            args.MicrOMEGAs = False

    if args.HiggsBounds or args.All:
        if 'HiggsBounds' not in paths_dict:
            log.warning('Could not find HiggsBounds path!')
            args.HiggsBounds = False

    if args.HiggsSignals or args.All:
        if 'HiggsSignals' not in paths_dict:
            log.warning('Could not find HiggsSignals path!')
            args.HiggsSignals = False  

    if args.HiggsTools or args.All:
        if 'HiggsBoundsRepo' not in paths_dict or 'HiggsSignalsRepo' not in paths_dict:
            log.warning('Could not find paths for HiggsTools repositories!')
            args.HiggsTools = False  

        else:
            try: 
                import Higgs
            except:
                args.HiggsTools = False
        
    if args.anyBSM:
        try:
            import anyBSM
        except:
            log.warning('Could not import anyBSM: disabling configuration')
            args.anyBSM = False    

    if args.VevaciousPlusPlus or args.All:
        if 'Vevacious++' not in paths_dict:
            log.warning('Could not find Vevacious++ path!')
            args.VevaciousPlusPlus = False  

    try:
        OF=open(scriptname,'w')
    except:
        log.error("Could not open script " + scriptname + " for writing")
        raise SystemExit

    OF.write('<<\"'+str(os.path.join(paths_dict["SARAH"],"SARAH.m"))+'\";\n')
    OF.write('Start[\"'+model+'\"];\n')
    #OF.write('DirectoryNamesFile=OpenWrite[ToFileName["'+scriptdir+'",modelname+"_DirectoryNames.txt"]];\n')
    #OF.write('WriteModelDirectories=True;')
    if len(args.Init) > 0:
        OF.write(args.Init+'\n')
    if args.SkipFlavorKit:
        OF.write('SkipFlavorKit=True;\n')
    #print(args.Options)

    if args.SPheno:
        OF.write('MakeSPheno['+args.Options+'];\n')

    # BSMArt path logic for AuxFiles
    # Updated to use importlib.resources for package compatibility
    
    aux_script_name = 'SARAH_BSMArt_extra.m'
    aux_script_path = None
    
    try:
        import importlib.resources
        from importlib.resources import files
        
        # Locate the resource
        # We need to extract it to a file accessible by Mathematica
        # We copy it to the script directory (where we are generating scripts)
        
        # Assuming bsmart.data.AuxFiles is a package/resource location
        # Check if we can find it
        ref = files('bsmart.data.AuxFiles').joinpath(aux_script_name)
        if ref.is_file():
             dest_path = os.path.join(scriptdir, aux_script_name)
             # Read content and write to dest, or use as_file to copy?
             # Simple read/write is safest for text files
             with importlib.resources.as_file(ref) as f:
                 shutil.copy(str(f), dest_path)
             aux_script_path = dest_path.replace('\\','/') # Ensure forward slashes for Mathematica if needed
    except Exception as e:
        log.warning('Could not extract ' + aux_script_name + ' from resources: ' + str(e))
        # Fallback to local check or paths_dict if legacy
        if 'BSMArt' in paths_dict:
             p = os.path.join(paths_dict["BSMArt"],'AuxFiles',aux_script_name)
             if os.path.exists(p):
                  aux_script_path = p
        
    if not aux_script_path:
        # Final fallback, assume it's in CWD or let SARAH fail if not found
        aux_script_path = aux_script_name
        
    # Check existence
    if not os.path.exists(aux_script_path) and not os.path.isfile(aux_script_path):
         log.error('Could not find auxiliary SARAH commands ('+aux_script_name+')!')
         raise SystemExit

    OF.write('<<\"'+aux_script_path+'\";\n')


    bsmartdatafile=os.path.join(scriptdir,'BSMArt_data_'+modelname+'.json')

    OF.write('fname=OpenWrite[\"'+bsmartdatafile+'\"];\n')
    OF.write('WriteString[fname,"{"];\n')
    OF.write('WriteBSMArtData[fname];\n')

    if (args.HiggsSignals or args.HiggsBounds):
            doHBHS=True
            OF.write('WriteHBData[fname];\n')
    elif args.HiggsTools:
            OF.write('WriteHBData[fname];\n')

        

    OF.write('WriteString[fname,"}"];\n')

    if args.SPheno:
        qnumbers_slha=os.path.join(scriptdir,'QNUMBERS_'+modelname+'.slha')
        OF.write('MakeQNUMBERSBSMArt["'+qnumbers_slha+'",If[Head[ExtraOddParticles] == List, ExtraOddParticles, {}]];\n')
        
    OF.write('Close[fname];\n')


    #if args.MicrOMEGAs or args.All:
    if args.MicrOMEGAs:
        OF.write('MakeCHep['+args.MOOptions+'];\n')

    #if args.VevaciousPlusPlus or args.All:
    if args.VevaciousPlusPlus:
        OF.write('MakeVevacious[Version->"++"];')

    if args.UFO:
        if args.anyBSM:
            log.info('Configure UFO to make all vertices')
            OF.write('MakeUFO[Exclude -> {}];')
        else:
            OF.write('MakeUFO[];')
        
    #OF.write('Close[DirectoryNamesFile];\n')
    OF.close()

        
    ### Now run the script


    log.info('Running SARAH')
    try:
        if args.debug:
            extra_debug_info=shell_command(mathexec + ' < '+scriptname,cwd,True)
            log.debug('SARAH output: '+extra_debug_info)
        else:
            _ = shell_command(mathexec+' < '+scriptname,cwd,True)
    except:
        log.critical('Error running Mathematica script')
        raise SystemExit

    try:
        with open(bsmartdatafile) as bsmart_json_data:
            bsmart_dict = json.load(bsmart_json_data, object_pairs_hook=collections.OrderedDict)
    except Exception as e:
        #log.error('Failed to load json file '+file)
        log.error('Failed to load json file '+bsmartdatafile)
        log.error('Json exception given as ' +str(e))
        raise SystemExit    

    ### Now start moving stuff around and compiling

    #sarahsphdir=os.path.join(paths_dict["SARAH"],'Output',modelname,'EWSB','SPheno')

    if runFlexibleSUSY:
        log.info('Configuring FlexibleSUSY for model')
        try:
            #with open(bsmartdatafile,'w') as bsmart_json_data:
            #    json.dump(bsmart_dict, bsmart_json_data, indent=4)
            if 'MathematicaUserBaseDirectory' in paths_dict:
                os.environ["MATHEMATICA_USERBASE"] = paths_dict['MathematicaUserBaseDirectory']
                os.environ["WOLFRAM_USERBASE"] = paths_dict['MathematicaUserBaseDirectory']
            #os.environ["MATHEMATICA_USERBASE"] = paths_dict['FlexibleSUSY']
            os.chdir(paths_dict['FlexibleSUSY'])
            #_ = shell_command('math < '+scriptname,cwd,True)
            
            _ = shell_command('./createmodel --name=%s -f' % args.FlexibleSUSY,paths_dict['FlexibleSUSY'],True,True) # errors fatal, log output
            #os.system('./createmodel --name=%s -f' % args.FlexibleSUSY)
            
            try:
                looptoollib=get_lib_path(paths_dict['LoopTools'])
            except Exception as e:
                log.critical('Failed to get python library dir for LoopTools. '+str(e))
            try:
                collierlib=get_lib_path(paths_dict['COLLIER'])
            except Exception as e:
                log.critical('Failed to get python library dir for Collier. ' +str(e))


            inclist=' --with-loop-libraries=collier,looptools --with-looptools-libdir='+looptoollib+' --with-looptools-incdir='+os.path.join(paths_dict['LoopTools'],'include')
            inclist+=' --with-collier-libdir='+collierlib+' --with-collier-incdir='+os.path.join(paths_dict['COLLIER'],'include')
            cmd='./configure --with-math-cmd="%s" --with-models=%s' %(mathexec,args.FlexibleSUSY)
            cmd+= inclist
            
            #_ = shell_command(cmd,paths_dict['FlexibleSUSY'],True)
            #os.system(cmd)
            _ = shell_command(cmd,paths_dict['FlexibleSUSY'],True,True)
            try:
                import multiprocessing as mp
                cpus=mp.cpu_count()

                makeline='make -j%d' %cpus
            except:
                makeline='make'
            #os.system(makeline)
            _ = shell_command(makeline,paths_dict['FlexibleSUSY'],True,True)
            bsmart_dict['FSmodelname']=args.FlexibleSUSY

            os.chdir(cwd)
        except Exception as e:
            log.error('Failed to build FlexibleSUSY model: '+str(e))


    internalmodelname=str(bsmart_dict['BSMArtmodelname'])

    if args.SPheno:
        sarahsphdir=bsmart_dict["SPhenoOutputDir"]
        sphenotarget=os.path.join(paths_dict["SPheno"],modelname)
        if not os.path.exists(sarahsphdir):
            log.error("SPheno output not generated!")
            raise SystemExit

        if os.path.exists(sphenotarget):
            ## remove
            shutil.rmtree(sphenotarget)
        
        shutil.copytree(sarahsphdir,sphenotarget)

        log.info('Building SPheno')
        try:
            _ = shell_command("make Model="+modelname,paths_dict["SPheno"],True)
        except:
            log.critical('Error compiling SPheno script')
            raise SystemExit
    else: # remove the SPhenoOutputDir from the json data file, that way we know SPheno wasn't built
        if 'SPhenoOutputDir' in bsmart_dict:
            bsmart_dict.pop('SPhenoOutputDir')



        
    #if args.MicrOMEGAs or args.All:
    if args.MicrOMEGAs:
        log.info('Building MicrOMEGAs')
        sarahmodir=os.path.join(paths_dict["SARAH"],'Output',modelname,'EWSB','CHep')
        targetmodir=os.path.join(paths_dict["MicrOMEGAs"],'SARAH_'+modelname)
        
        if os.path.exists(targetmodir):
            shutil.rmtree(targetmodir)
        
        
        
        try:
            _ = shell_command('./newProject SARAH_'+modelname,paths_dict["MicrOMEGAs"],True)
        except:
            log.error('Failed to create new MicrOMEGAs project')
            raise SystemExit
                    
        
        #for ff in ['prtcls1.mdl','func1.mdl','extlib1.mdl','vars1.mdl','lgrng1.mdl']:
        for ff in ['prtcls1.mdl','func1.mdl','vars1.mdl','lgrng1.mdl']:
            shutil.copy(os.path.join(sarahmodir,ff),os.path.join(targetmodir,'work','models'))

        
        MOver=os.path.basename(paths_dict["MicrOMEGAs"])
        if len(MOver) < 11 or MOver[10] !='_': ## wtf? Assume latest version
            MOver='6.0.0'
        else:
            MOver=MOver[11:]
        
        if MOver[1] !='.':
            MOver='6.0.0'
        MOver=float(MOver[0]+MOver[2:])
        
        if args.BSMArt or args.All:
            # Determined executale/cpp name
            if MOver<52:
                moexec='MicrOmegas_v5.0_BSMArt'
            elif MOver<53.41:
                moexec='MicrOmegas_v5.2_BSMArt'
            elif MOver<60:
                moexec='MicrOmegas_v5.3_BSMArt_xsections'
            else:
                moexec='MicrOmegas_v6.1_BSMArt_xsections'
            
            mocpp=moexec+'.cpp'

            # Logic to find mocpp file using importlib.resources
            src_mocpp = None
            try:
                import importlib.resources
                from importlib.resources import files
                ref = files('bsmart.data.AuxFiles').joinpath(mocpp)
                if ref.is_file():
                    with importlib.resources.as_file(ref) as f:
                        shutil.copy(str(f),targetmodir)
                    # It's now in targetmodir as mocpp (shutil.copy copies filename)
                    src_mocpp = os.path.join(targetmodir, mocpp)
            except Exception as e:
                log.warning('Could not find ' + mocpp + ' in resources: ' + str(e))
                # Fallback
                if 'BSMArt' in paths_dict:
                    p = os.path.join(paths_dict["BSMArt"],'AuxFiles',mocpp)
                    if os.path.exists(p):
                        shutil.copy(p,targetmodir)
                        src_mocpp = os.path.join(targetmodir, mocpp)

            if not src_mocpp or not os.path.exists(src_mocpp):
                 log.error('Could not find BSMArt micrOMEGAs driver: '+mocpp)
                 raise SystemExit

            command='make main='+mocpp
            MOexec=os.path.join(targetmodir,moexec)
            
        else:
            MOexec=os.path.join(targetmodir,'CalcOmega_with_DDetection_MOv5')
            command='make main=CalcOmega_with_DDetection_MOv5.cpp'
            shutil.copy(os.path.join(sarahmodir,'CalcOmega_with_DDetection_MOv5.cpp'),targetmodir)
            shutil.copy(os.path.join(sarahmodir,'CalcOmega_MOv5.cpp'),targetmodir)

        try:    
            _ = shell_command(command,targetmodir,False)
        except:
            log.error('Failed to compile MicrOMEGAs project')
            raise SystemExit
        
        if not os.path.exists(MOexec):
            log.error('Error compiling MicrOMEGAs code (no executable found)')
            raise SystemExit
        
        #if MOver > 53.4 and args.BSMArt:
        #    MOexec=MOexec+' --xsections'

        bsmart_dict['MicrOMEGAs executable'] = MOexec
        bsmart_dict['MicrOMEGAs version'] = MOver

    # This was at the end of micromegas; now I will just do it automatically because of updates e.g. in FlexibleSUSY 
    with open(bsmartdatafile, 'w') as fp:
            json.dump(bsmart_dict, fp, indent=4)
    ### Now set up BSMArt scan directory, if required

    if not args.local:
        from bsmart import get_bsmart_config_dir
        global_config_dir = get_bsmart_config_dir()
        if os.path.isdir(global_config_dir):
            global_bsmartdatafile=os.path.join(global_config_dir,'BSMArt_data_'+modelname+'.json')
            with open(global_bsmartdatafile, 'w') as fp:
                json.dump(bsmart_dict, fp, indent=4)
    #if not args.BSMArt and not args.All:
    if not args.BSMArt:
        sys.exit()

    log.info("Setting up BSMArt")



    includecodes={ 'SPheno':args.SPheno, 'FlexibleSUSY': runFlexibleSUSY, 'MicrOMEGAs':args.MicrOMEGAs, 'Vevacious++':args.VevaciousPlusPlus, 'HiggsSignals':args.HiggsSignals, 'HiggsBounds':args.HiggsBounds, 'HiggsTools':args.HiggsTools, 'SModelS': True, 'Flavio': True}

    prepareBSMArt.makeBSMArtTemplates(modelname,bsmart_dict,paths_dict,includecodes)    
        

if __name__ == "__main__":
    main()
