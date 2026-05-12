#!/usr/bin/env python3

import argparse
import os
from datetime import datetime
import sys
import shutil
import json
import collections
import importlib
import importlib.util

# Package imports
import bsmart
import bsmart.tools
from bsmart import debug
from bsmart.scans import get_available_scans

# from bsmart import __version__
# from bsmart import __build_date__

# Global variables for compatibility
MAINpath = os.getcwd() # Or appropriate default

# Read the input file and extract all necessary information
def parse_input(file):
    
    try:
        with open(file) as json_data:
            d = json.load(json_data, object_pairs_hook=collections.OrderedDict)
    except Exception as e:
        #log.error('Failed to load json file '+file)
        print('Failed to load json file '+file)
        print('Json exception given as ' +str(e))
        raise SystemExit
    #log.debug('Content of input file: %s' % str(d))

    if 'Observables' not in d:
        d['Observables'] = collections.OrderedDict()
    if 'StoreAllPoints' not in d['Setup']:
        d['Setup']['StoreAllPoints'] = "False"
    if 'Cores' not in d['Setup']:
        d['Setup']['Cores'] = 1
    return d

def print_scan_and_tool_help(help_name:str=None):
    """
        New feature to help users: list scans and tools if no argument is given, otherwise print help about that scan if it is 
    """
    
    local_scanpath=os.path.join(os.getcwd(),'Scans')

    # This function walks the subdirectories recursively and gives information about the scans
    allscans=get_available_scans(None,with_meta = True)
    locals=get_available_scans(local_scanpath,with_meta = True)
    #allscans.extend(locals)
    allscans.update(locals) # dictionary now

    # main_tooldir=os.path.join(bsmartpath,'tools')
    # Using relative import logic or entry_points for tools might be better, 
    # but for now let's assume valid package structure finding tools.
    # We might need a better way to find the installed 'tools' directory.
    # existing logic: bsmartpath was '.../bsmart' (the package dir)
    # new bsmartpath: '.../site-packages/bsmart'
    main_tooldir = bsmart.tools

    local_toolpath=os.path.join(os.getcwd(),'Tools')

    

    if help_name is not None and help_name in allscans: # print all info
        scan = allscans[help_name]
        print("Scan %s"% help_name)
        print('Path: %s' % scan['path'])
        if 'meta' in scan and len(scan['meta']) >0 : # metadata provided
            print('---------')
            print('Metadata:')
            print('---------')
            for key,val in scan['meta'].items():
                # if val is a dict, print it in a better way
                if isinstance(val, dict):
                    print("   * %s:" % key)
                    for k,v in val.items():
                        print("     * %s: %s" % (k,v))
                else:
                    print("   * %s: %s" % (key,val))
        if 'documentation' in scan and scan['documentation'] is not None:
            print('--------------')
            print("Documentation:")
            print('--------------')
            print(scan['documentation'])
        sys.exit(0)
    
    alltools = get_available_scans(main_tooldir,with_meta = True)
    alltools.update(get_available_scans(local_toolpath,with_meta = True))

    if help_name is not None and help_name in alltools: # print all info
        tool = alltools[help_name]
        print("Tool %s"% help_name)
        print('Path: %s' % tool['path'])
        if 'meta' in tool and len(tool['meta']) >0 : # metadata provided
            print('---------')
            print('Metadata:')
            print('---------')
            for key,val in tool['meta'].items():
                # if val is a dict, print it in a better way
                if isinstance(val, dict):
                    print("   * %s:" % key)
                    for k,v in val.items():
                        print("     * %s: %s" % (k,v))
                else:
                    print("   * %s: %s" % (key,val))
        if 'documentation' in tool and tool['documentation'] is not None:
            print('--------------')
            print("Documentation:")
            print('--------------')
            print(tool['documentation'])
        sys.exit(0)

    # Print everything

    if help_name is not None:
        print('Scan or Tool %s not found' %help_name)

    print('----------------')
    print("Available scans:")
    print('----------------')
    for id,scan in allscans.items():
        print(" - %s" % id)
        
        # if 'meta' in scan and len(scan['meta']) >0 : # metadata provided
        #     for key,val in scan['meta'].items():
        #         print("   * %s : %s" % (key,val))
        

    print('----------------')
    print("Available tools:")
    print('----------------')
    for id,tool in alltools.items():
        print(" - %s" % id)
        
        # if 'meta' in tool and len(tool['meta']) >0 : # metadata provided
        #     for key,val in tool['meta'].items():
        #         print("   * %s : %s" % (key,val))
        
    print('----------------')
    print("Type BSMArt --Settings <scan name> or BSMArt --Settings <tool name> for more information")
    print('----------------')

    

    
    # print("\nAvailable tools:")
    # main_tooldir=os.path.dirname(bsmart.tools.__file__)
    # local_toolpath=os.path.join(os.getcwd(),'Tools')
    # alltooldirs=[main_tooldir,local_toolpath]
    # for tooldir in alltooldirs:
    #     if os.path.isdir(tooldir):
    #         for tool in sorted(os.listdir(tooldir)):
    #             if tool.endswith('.py') and not tool.startswith('_'):
    #                 print(" - %s" % tool[:-3])
    # 

def main():
    try:
        parser = argparse.ArgumentParser(
            description='Please give the name of the input file.')
        parser.add_argument('inputfile',
                            metavar='File', type=str,
                            help='Input file name',default=None, nargs='?')
        # Are --short and --csv really necessary any more? 
        parser.add_argument("--short", help="Store output in short tabbed form",
                            action="store_true")
        parser.add_argument("--csv", help="Store output in csv form",
                            action="store_true")
        parser.add_argument("--debug", help="write debug information",
                            action="store_true")
        parser.add_argument("--NoMPI", help="Do not check for MPI",
                            action="store_true")
        """
            New feature to help users (not yet implemented) 
        """
        parser.add_argument("--Settings", help="List all settings for your selected scan and tools!",
                            action="store_true")
        args = parser.parse_args()
    except Exception as e:
        print("Error parsing arguments", e)
        # print("No input file provided\n")
        # print_scan_and_tool_help()
        # raise SystemExit
        return

    if not args.inputfile and not args.Settings:
        print("No input file provided\n")
        print_scan_and_tool_help()
        return

    # set up paths
    cwd = os.getcwd()
    
    #BSMArtVersion='PRIVATEVERSION'
    #BSMArtDate='RELEASEDATE'

    # BSMArtVersion=__version__
    # BSMArtDate=__build_date__
    
    
    if args.NoMPI == True:
        MPIsize=0
        MPIrank=0
    else:
        try:
            from mpi4py.MPI import COMM_WORLD
            MPIsize=COMM_WORLD.Get_size() 
            if MPIsize > 1: # need parallel capabilities
                MPIrank=COMM_WORLD.Get_rank()
            else:
                MPIrank=0
        except:			
            MPIsize=0
            MPIrank=0
    
    
    if MPIrank ==0:
        now = datetime.now()
        timestamp = now.strftime('%Y%m%d_%H_%M_%S.%f')
        debug.print_BSMArt_splash()
        # print('                                                                ')
        # print('888888b.    .d8888b.  888b     d888        d8888         888    ')
        # print('888  "88b  d88P  Y88b 8888b   d8888       d88888         888    ')
        # print('888  .88P  Y88b.      88888b.d88888      d88P888         888    ')
        # print('8888888K.   "Y888b.   888Y88888P888     d88P 888 888d888 888888 ')
        # print('888  "Y88b     "Y88b. 888 Y888P 888    d88P  888 888P"   888    ')
        # print('888    888       "888 888  Y8P  888   d88P   888 888     888    ')
        # print('888   d88P Y88b  d88P 888   "   888  d8888888888 888     Y88b.  ')
        # print('8888888P"   "Y8888P"  888       888 d88P     888 888      "Y888 ')
        # print('                                                                ')
        # print('Version %s [%s]' %(BSMArtVersion,BSMArtDate))
        # print('By M. Goodsell                                                  ')
        # print('With contributions from:                                        ')
        # print(' * A. Joury (R. Moutafis)                                       ')
        # print(' * Luc Darmé                                                    ')
        # print(' * Martin Gabelmann                                             ')
        # print(' * Johannes Braathen                                            ')
        # print(' * Farid Ibrahimov                                              ')
        # print(' * Wojciech Kotlarski                                           ')
        # print('                                                                ')
    else:
        now = None
        timestamp=None
    
    if MPIsize > 1:
        now = COMM_WORLD.bcast(now,root=0)
        timestamp = COMM_WORLD.bcast(timestamp,root=0)
        COMM_WORLD.Barrier()
    
    
    
    if args.Settings:
        print_scan_and_tool_help(args.inputfile)
        return
    
    
    
    if not os.path.isfile(args.inputfile):
        print_scan_and_tool_help(args.inputfile)
        return
    
    ### Now load in the input file
    try:
        inputs=parse_input(args.inputfile)
    except:
        print('Failed to load input file %s!' %args.inputfile)
        #print_scan_and_tool_help(args.inputfile)
        raise SystemExit
    
    
    # Now set up configuration and the scan
    
    
    if args.short:
        inputs['Setup']['Short'] = "True"
    
    if args.csv:
        inputs['Setup']['csv'] = "True"
    
    if args.debug:
        inputs['Setup']['debug'] = "True"
    
    
    if 'Temporary Directory' in inputs['Setup']:
        tempstub=inputs['Setup']['Temporary Directory']
    else:
        tempstub=os.path.join(cwd,"Temp",timestamp)
    
    if MPIsize > 1:
        #tempstub=os.path.join(cwd,"Temp",timestamp,"MPI_"+str(MPIrank))
        tempstub=os.path.join(tempstub,"MPI_"+str(MPIrank))
    #else:
    #    tempstub=os.path.join(cwd,"Temp",timestamp)
    
    
    if 'debug' in inputs['Setup'] and (eval(inputs['Setup']['debug'])):
        temporary_dir = tempstub
    else:
        inputs['Setup']['debug'] = "False"
        if ('Temporary Directory' not in inputs['Setup']) and os.path.isdir('/dev/shm'):
            # use ramdisk to exchange files
            if MPIsize > 1:
                temporary_dir = os.path.join('/dev/shm/BSMART_Temp', timestamp, 'MPI_'+str(MPIrank))
            else:   
                temporary_dir = os.path.join('/dev/shm/BSMART_Temp', timestamp + '_PID' + str(os.getpid()))
        else:
            temporary_dir = tempstub
    
    os.makedirs(temporary_dir)    
    
    if MPIsize > 1:
        if MPIrank ==0:
            MAINlog = debug.new_logger(args.debug,'BSMArt main 0',os.path.join(temporary_dir,'BSMArt.log'))
        else:
            MAINlog = debug.new_logger(False,'BSMArt main '+str(MPIrank),os.path.join(temporary_dir,'BSMArt.log'))
    else:
        MAINlog = debug.new_logger(args.debug,'BSMArt main',os.path.join(temporary_dir,'BSMArt.log'))
        
    if not inputs['Setup']['RunName']:
        inputs['Setup']['RunName'] = 'BSMArt_run'
            
    inputs['Setup']['cwd'] = cwd
    inputs['Setup']['TempDir'] = temporary_dir
    inputs['Setup']['MAINpath']= MAINpath
    inputs['Setup']['PackagePath']= os.path.dirname(bsmart.__file__)
    
    """
    Allow setting of environment variables
    """
    if 'Environment Variables' in inputs['Setup']:
        for var,val in inputs['Setup']['Environment Variables'].items():
            os.environ[var]=val
    
    
    
    
    
    if 'MPI' not in inputs:
        inputs['MPI'] = {}
    inputs['MPI']['Size'] = MPIsize
    inputs['MPI']['Rank'] = MPIrank
     
    if MPIsize > 1:
        inputs['Setup']['Cores'] = 1
    
    
    # checking for the correct scan
    scan_initialised = False
    
    scan_type=inputs['Setup']['Type']
    
    # Do a quick check for the default scan specification. If it is not there, then we will get the scan by more sophisticated means
    
    scan_filename=scan_type+'.py'
    if importlib.util.find_spec('bsmart.scans.'+scan_type):
        import_string='bsmart.scans.'+scan_type
    else:
        allscans=get_available_scans(None,with_meta = False)
        if scan_type in allscans:
            import_string = allscans[scan_type]['module']
        else:
            
            local_scanpath=os.path.join(os.getcwd(),'Scans')
            
            locals=get_available_scans(local_scanpath,with_meta = False)
            try:
                import_string = locals[scan_type]['module']
                sys.path.append(local_scanpath)
            except:
            
                MAINlog.error('Fatal error: could not find file %s for scan type %s '%(scan_filename,scan_type))
                raise SystemExit
    
    try:
        
        # importlib.import_module(import_string)
        scan_class=importlib.import_module(import_string)
        scan = scan_class.NewScan(inputs, MAINlog)
        MAINlog.debug('Found scan ' + scan_type)
        scan_initialised = True
    except Exception as e:
        MAINlog.error('Fatal error, %s: could not initialise scan type %s '%(str(e),scan_type))
        raise SystemExit
                    
    
    ## Scan should have been initialised and is now class scan
    
    MAINlog.info('Starting scan')
    scan.launch()

if __name__ == "__main__":
    main()




                        
    
