
import os
from datetime import datetime
import sys

import shutil
import json
import collections

import importlib



"""

MAINpath = os.path.split(os.path.dirname(os.path.realpath( __file__ )))[0]



if not os.path.isdir(MAINpath):
    sys.exit('Problem locating BSMArt root directory:\n' + MAINpath)

## Insert root directory to path so that we can import our modules ...
sys.path.insert(0, MAINpath)

## Maybe this is dangerous
bsmartpath=os.path.join(MAINpath,'bsmart')
sys.path.insert(0, bsmartpath)
"""
import bsmart.zslha as zslha
#### Import package files

import bsmart.debug as debug

import bsmart.HEPRun as HEPRun
from bsmart.BSMlikelihood import MakeLikelihoods, safe_float, MakeGlobalLikelihood

import bsmart.BSMutil as BSMutil



def get_blackbox(file,args_dict=None,post_process=None,likelihood_return_type='Likelihood'):
    if args_dict is None:
        args={'NoMPI':True,'short':True,'csv':True,'debug':True}
    else:
        args=args_dict

    
    try:
        with open(file) as json_data:
            inputs = json.load(json_data, object_pairs_hook=collections.OrderedDict)
    except Exception as e:
        #log.error('Failed to load json file '+file)
        print('Failed to load json file '+file)
        print('Json exception given as ' +str(e))
        raise SystemExit
    #log.debug('Content of input file: %s' % str(d))
    
    if 'Observables' not in inputs:
        inputs['Observables'] = collections.OrderedDict()
    if 'StoreAllPoints' not in inputs['Setup']:
        inputs['Setup']['StoreAllPoints'] = "False"
    if 'Cores' not in inputs['Setup']:
        inputs['Setup']['Cores'] = 1


    if 'NoMPI' not in args: # must have been passed to the function by user
         args['NoMPI'] = False
    

    # set up paths
    cwd = os.getcwd()
    
    if args['NoMPI']:
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
    else:
        now = None
        timestamp=None

    if MPIsize > 1:
        now = COMM_WORLD.bcast(now,root=0)
        timestamp = COMM_WORLD.bcast(timestamp,root=0)
        COMM_WORLD.Barrier()

    # Now set up configuration and the scan 
    
    if 'short' in args and args["short"]:
        inputs['Setup']['Short'] = "True"
    
    if 'short' not in inputs['Setup']:
        inputs['Setup']['Short'] = "False"


    if 'csv' in args and args["csv"]:
        inputs['Setup']['csv'] = "True"
    
    if 'csv' not in inputs['Setup']:
        inputs['Setup']['csv'] = "False"

    if 'debug' in args and args["debug"]:
        inputs['Setup']['debug'] = "True"
    if 'debug' not in inputs['Setup']:
        inputs['Setup']["debug"] = "False"

    do_debug=False
    if 'debug' in inputs['Setup'] and (eval(inputs['Setup']['debug'])):
         do_debug=True
    
    if 'Temporary Directory' in inputs['Setup']:
        tempstub=inputs['Setup']['Temporary Directory']
    else:
        
        tempstub=os.path.join(cwd,"Temp",timestamp)
    
    if MPIsize > 1:
        #tempstub=os.path.join(cwd,"Temp",timestamp,"MPI_"+str(MPIrank))
        tempstub=os.path.join(tempstub,"MPI_"+str(MPIrank))
    #else:
    #    tempstub=os.path.join(cwd,"Temp",timestamp)
   
    temporary_dir = ""
    if do_debug:
        
        temporary_dir = tempstub
    else:
        
        inputs['Setup']['debug'] = "False"
        
        if ('Temporary Directory' not in inputs['Setup']) and os.path.isdir('/dev/shm'):
            # use ramdisk to exchange files
            if MPIsize > 1:
                temporary_dir = os.path.join('/dev/shm/BSMART_Temp', timestamp, 'MPI_'+str(MPIrank))
            else:   
                temporary_dir = os.path.join('/dev/shm/BSMART_Temp', timestamp)
        else:
            temporary_dir = tempstub
    
    os.makedirs(temporary_dir)    

    

    if MPIsize > 1:
        if MPIrank ==0:
            MAINlog = debug.new_logger(do_debug,'BSMArt main 0',os.path.join(temporary_dir,'BSMArt.log'))
        else:
            MAINlog = debug.new_logger(False,'BSMArt main '+str(MPIrank),os.path.join(temporary_dir,'BSMArt.log'))
    else:
        MAINlog = debug.new_logger(do_debug,'BSMArt main',os.path.join(temporary_dir,'BSMArt.log'))
        
    if not inputs['Setup']['RunName']:
        inputs['Setup']['RunName'] = 'BSMArt_run'
            
    inputs['Setup']['cwd'] = cwd
    inputs['Setup']['TempDir'] = temporary_dir
    #inputs['Setup']['MAINpath']= MAINpath
    #inputs['Setup']['PackagePath']= bsmartpath

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


    ## Now to import HEPRun etc

    
    runsettings=HEPRun.RunSettings(inputs)
    varnames = runsettings.varnames
    
    doplots=False
            
    if 'Plotting' in inputs:
        plotting=inputs['Plotting']
        plots=inputs['Plotting']['Plots']
        if plotting['Strategy']=='All' or plotting['Strategy']=='Valid':
            runsettings.store_points=True

        
        if len(plots) > 0:
            #import BSMplots
            doplots=True
            #BSMplots.say_hello(self)
            plotdir=os.path.join(inputs['Setup']['cwd'], inputs['Setup']['RunName'], 'Plots')


    # Create output and inputs directories 
            
    if os.path.exists(runsettings.output_dir):
        shutil.rmtree(runsettings.output_dir)
                
    os.makedirs(runsettings.output_dir)

    if runsettings.store_everything:
                if os.path.exists(runsettings.all_dir):
                    shutil.rmtree(runsettings.all_dir)
                
                os.makedirs(runsettings.all_dir)

            
    if os.path.exists(runsettings.inputs_dir):
                shutil.rmtree(runsettings.inputs_dir)    
                
    if runsettings.store_inputs:
                os.makedirs(runsettings.inputs_dir)

    if doplots:
                if inputs['MPI']['Rank'] == 0: ### Only one set of plots ... 
                    if os.path.exists(plotdir):
                        shutil.rmtree(plotdir)
                    os.makedirs(plotdir)		

    MAINlog.info('Parse Settings')
    MAINlog.info('Reading bare spc: '+str(runsettings.input_filename))

    
    try:
        inputspc= zslha.read(runsettings.input_filename,verbose=False)
    except Exception as e:
        MAINlog.error('Problem reading slha file '+runsettings.input_filename+', ' + str(e))
        inputspc = zslha.SLHA()

    
        

    
    
    

    

    def write_lh_file(point, dir, name,values_dict=None,Overloads=None):
        """ write Les Houches input file for given parameter point """
        #global inputspc,varnames  ## I think this binds the values of varnames and inputspc
        lhname=os.path.join(dir,name)
        
        #runspc=inputspc.deepcopy()
        runspc = inputspc
        if Overloads is not None:  
            """
                Allow overloading of block values, e.g. during marginalisation
                Format: a list of form
                [["SPHENOINPUTS",[1],1.0],...]
            """
            runspc = zslha.SLHA(inputspc) #
            for overload in Overloads:
                if len(overload) < 3:
                    continue
                block=overload[0].upper()
                blockkey=str(overload[1])[1:-1].replace(" ", "")
                value=overload[2]
                if block in runspc.blocks:
                    runspc.blocks[block][blockkey] = value

        if values_dict is None:
            zslha.write_lh(runspc, point, lhname,varnames)
        else:
            zslha.write_lh_dict(runspc, values_dict, lhname)


    if post_process is None:
        import math
        #use_loglike=eval(inputs['Setup'].get('LogLike',"True"))
        if runsettings.store_points_in_memory and runsettings.store_invalid_points:
            naive = False 
            runsettings.invalid_return_value = 0
        else:
            naive = True  # we treat invalid points as bad
            runsettings.invalid_return_value = []


            
        # if likelihood_return_type.upper() == 'LIKELIHOOD':
        #     use_loglike=False
        # else:
        #     use_loglike=True
          
        # Create the likelihood functions based on 'scaling' property of variables, if they have it; otherwise use gaussian likelihood    


        
        get_likelihood,_=MakeGlobalLikelihood(inputs['Observables'],return_type=likelihood_return_type)
        def postprocess(Point, observables, data_point,temp_dir,log, lock=None):
            """ return the likelihood; we won't get this far if the point failed to be generated """
            # For this algorithm we are minimising the objective function
            return get_likelihood(observables)
        # else:
        #     likelihood_fns, observable_masks = MakeLikelihoods(inputs["Observables"], loglike=True)

        #     loglike = True
        #     multiplier = 1.0
        #     if likelihood_return_type is not None:
        #         if likelihood_return_type.upper() == "LIKELIHOOD":
        #             loglike = False
                     
        #         elif likelihood_return_type.upper() == "LL":
        #             loglike = True
        #         elif likelihood_return_type.upper() == "NLL":
        #             loglike = True
        #             multiplier = -1.0

            


        #     if loglike:
        #         import sys
        #         maxloss = sys.float_info.max
        #         def smooth_cap_loss(x):
        #             """
        #             Caps the loss by applying a sigmoid.
        #             This is useful for losses that are unbounded.                   

        #             """
        #             return maxloss*math.expm1(x/maxloss) # 


        #         def postprocess(Point, observables, data_point,temp_dir,log, lock=None):
        #             """ return the likelihood; we won't get this far if the point failed to be generated """
        #             # For this algorithm we are minimising the objective function
        #             likeit=iter(likelihood_fns)

        #             return multiplier*math.sum([self.smooth_cap_loss((next(likeit))(val)) if mask and not math.isnan(val := safe_float(v)) else float((next(likeit) and False) or maxloss) 
        #                                     for v, mask in zip(observables, observable_masks) if mask])


        #     else:
        #         def postprocess(Point, observables, data_point,temp_dir,log, lock=None):
        #             """ return the likelihood; we won't get this far if the point failed to be generated """
        #             # For this algorithm we are minimising the objective function
        #             likeit=iter(likelihood_fns)

        #             lh = 1.0

        #             for v, mask in zip(observables, observable_masks):
        #                 if mask:
        #                     val =  next(likeit)(float(v))
        #                     if math.isnan(val): # fail -> zero likelihood
        #                         return 0.0
        #                     else:
        #                         lh *= val

                    
        #             return lh

        

        RunManager=HEPRun.HEPRun(runsettings, write_lh_file, postprocess,MAINlog)

    else:

        RunManager=HEPRun.HEPRun(runsettings, write_lh_file, post_process,MAINlog)



    return RunManager,inputs

