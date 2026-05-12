"""
Simple MCMC scan


Simple MCMC scan implementing the Metropolis-Hastings algorithm

Launches a seperate scan for each core, with the option to merge the
results at the end via "Merge Results":"True" in "Setup". Of course, will
also run in single-core mode.

For variables we require a range and a variance for the jumps:

.. code-block:: text

    "Variables": {
         "lambda": { "RANGE": [0.01,1.0], 
                     "VARIANCE": 0.2}
         },


Several possibilities for likelihood functions are provided -- see the BSMlikelihood file.

E.g. gaussian, Log likelihood, "USER", or "OFF" for observables that should not be taken into account for the likelihood calculation (but the user wants to retain them in the summary file for other purposes). These are set by providing e.g.:

.. code-block:: text

    "Observables": {
    	"mh" : { "SLHA": ["MASS", [25]], "SCALING": "GAUSS",
    		"MEAN": 125.0,
    		"VARIANCE":3.0
    	      }


or "SCALING":"LOG", or "SCALING":"OFF"

It is not necessary to specify the scaling explicitly: gaussian is taken by default

The results are summarised in a csv file (or several) in the Results subdirectory of the scan, which contains only the points retained by the scan (along with an ID for each point); the raw points will as usual be retained in the Spectrum_Files directory if so desired. 

"""

__meta__ = {
    "name": "MCMC",
    "requires": ["numpy"],
    "settings": {
        "Points": "Int",
        "Valid Points": "Int",
        "Merge Results": "Boolean",
        "Initial Sample": "Path to initial sample file",
        "Sample LogLike Name": "Column name for LogLike (if initial sample and relevant)",
        "Sample NLL Name": "Column name for NLL (if initial sample and relevant)",
        "Sample Likelihood Name": "Column name for Likelihood (if initial sample and relevant)"
    }
}



from bsmart.core import Scan as Scan

import numpy as np
import math
from bsmart import debug
import multiprocessing as mp

import queue
import time
import os
import shutil

from bsmart.HEPRun import UpdateProgress, DataPoint

from bsmart.BSMlikelihood import MakeLikelihoods, MakeGlobalLikelihood






class NewScan(Scan):
    """Scanner class for multicore MCMC Scans"""
    

    def __init__(self, inputs, log):
        Scan.__init__(self, inputs, log)
        self.header=[]
        self.header.append('ID')
        self.variances = []
        self.lock = None
        for v in self.inputs['Variables']:
            self.variances.append(self.inputs['Variables'][v]['VARIANCE'])
            self.header.append(v) ## append the variable name ... 

        # Create the likelihood functions based on 'scaling' property of variables, if they have it; otherwise use gaussian likelihood    
        
        for o in self.inputs['Observables']:
            self.header.append(o)
        
       
        #self.use_loglike=eval(self.inputs['Setup'].get('LogLike',"False"))

        #self.likelihoodfns=MakeLikelihoods(self.inputs['Observables'],self.use_loglike)
        #self.get_likelihood,_=MakeGlobalLikelihood(self.inputs['Observables'],self.use_loglike)
        self.get_likelihood,_=MakeGlobalLikelihood(self.inputs['Observables'],return_type="Likelihood")
        self.header.append('Likelihood')
        
        #self.headerlinestring='# '+', '.join(self.header)
        self.headerlinestring=','.join(self.header)
        if 'Points' in self.inputs['Setup']:
            self.maxpoints=self.inputs['Setup']['Points']
        else:
            self.maxpoints=None

        if 'Valid Points' in self.inputs['Setup']:
            self.maxvalidpoints=self.inputs['Setup']['Valid Points']
        else:
            self.maxvalidpoints=None

        if not self.maxvalidpoints and not self.maxpoints:
            log.error('No Points target or Valid Points target')
            raise

        if self.maxvalidpoints:
            self.target=self.maxvalidpoints
        else:
            self.target=self.maxpoints
        
        self.resultsdir=os.path.join(self.inputs['Setup']['cwd'], self.inputs['Setup']['RunName'],'Results')

        self.merge_results  = self.inputs['Setup'].get('Merge Results',False)

        """ 
            For resuming, or using an initial sample, we need just to have a starting point per core. 

            The initial sample should also be copied to the results
        
        """

        self.starting_points=[] # list of input values and their likelihoods
        self.merged_filename=os.path.join(self.resultsdir,'MCMC_results')
        self.number_of_prior_points=0
        self.number_of_prior_valid=0
        if self.resume: # this handled by core.py
            if os.path.exists(os.path.join(self.resultsdir,'MCMC_results')): # Previous scan actually completed!
                self.inputs['Setup']['Initial Sample'] = self.merged_filename
                self.number_of_prior_valid = max(0,len(list(open(self.merged_filename, 'r'))) - 1)
            elif all([os.path.exists(os.path.join(self.resultsdir,'MCMC_run_'+str(core))) for core in range(len(self.RunManager.CoreRuns))]): 
                # results files from all cores exist; read from each sample to get starting point, and merge the files into a results file if required
                totaldf = None
                for core in range(len(self.RunManager.CoreRuns)):
                    resultfile=os.path.join(self.resultsdir,'MCMC_run_'+str(core))
                    core_inputs,core_likes=self.get_previous_sample(resultfile)
                    if len(core_inputs) > 0:
                        self.starting_points.append([list(core_inputs[-1]),float(core_likes[-1])])
                        self.number_of_prior_valid += len(core_inputs)
                    else: # empty file
                        self.starting_points.append([[],0])

                              
            elif os.path.exists(self.RunManager.csv_filename):
                self.inputs['Setup']["Initial Sample"] = self.RunManager.csv_filename
                # we just treat this as an initial sample

            if os.path.exists(self.RunManager.csv_filename):
                # count number of lines in self.RunManager.csv_filename
                self.number_of_prior_points = max(0,len(list(open(self.RunManager.csv_filename, 'r'))) - 1)


                


        if  len(self.starting_points) == 0 and 'Initial Sample' in self.inputs['Setup'] and os.path.exists(self.inputs['Setup']['Initial Sample']):
            self.prior_sample, self.prior_likes = self.get_previous_sample(self.inputs['Setup']['Initial Sample'])
            # now must make some choices of where to start ...

            from bsmart.ml.seed_points import select_seeds

            # select seeds tries to minimis the nll
            selected_point_ids=list(select_seeds(self.prior_sample,-np.log(self.prior_likes),nseeds=self.runsettings.cores))

            self.starting_points=[ [self.prior_sample[id],self.prior_likes[id]] for id in selected_point_ids]

            if len(self.starting_points) < self.runsettings.cores: 
                self.starting_points = self.starting_points + [[],0] * (self.runsettings.cores - len(self.starting_points) )

        


        # set up results directory and files

        if os.path.exists(self.resultsdir) and not self.resume:
            shutil.rmtree(self.resultsdir)

        if not os.path.exists(self.resultsdir):
            os.makedirs(self.resultsdir)

        
        # initialise the scans
        self.ScanThreads=[]    
        for x in range(len(self.RunManager.CoreRuns)):
            self.ScanThreads.append(MCMC_onecore(self,x))
    
        if len(self.starting_points) > 0:
            for numthread,pt in enumerate(self.starting_points):
                log.info("Adding initial point " +str(pt))
                self.ScanThreads[numthread].last_point = list(pt[0])
                self.ScanThreads[numthread].likelihood_last_step = float(pt[1])

    def get_previous_sample(self,sample_file): 
        """
            Routine to load data from a previous sample and return a couple of numpy arrays of the values
        """
        self.log.info('Attempting to load prior sample '+sample_file)
        try:

            import pandas as pd
            df = pd.read_csv(sample_file)
            columns=df.columns.tolist()
            sample_inputs= df[self.inputs['Variables'].keys()].to_numpy()
            if  'Sample LogLike Name' in self.inputs["Setup"]:
                    likes = np.exp(df[self.inputs["Setup"]["Sample LogLike Name"]].to_numpy())
            elif  'Sample NLL Name' in self.inputs["Setup"]:
                    likes = np.exp(-df[self.inputs["Setup"]["Sample NLL Name"]].to_numpy())
            elif  'Sample Likelihood Name' in self.inputs["Setup"]:
                    likes = df[self.inputs["Setup"]["Sample Likelihood Name"]].to_numpy()
            elif  'Likelihood' in columns:
                    likes = df['Likelihood'].to_numpy()
            elif  'Result' in columns:
                    #if eval(self.inputs["Setup"].get("LogLike","True")):
                    #    likes = df['Result'].to_numpy()
                    #else:
                    #    likes = df['Likelihood'].to_numpy()
                    # Assume this comes from a prior MCMC scan
                    likes = df['Result'].to_numpy()
            else:
                    raise ValueError("Missing required columns for likelihood calculation")
            self.log.info('Loaded Initial sample '+sample_file)
            return sample_inputs,likes
        except Exception as e:
            self.log.error('Failed to load initial sample from file, will generate new one;  '+str(e))
            return np.empty(),np.empty()
     
    def get_next_pointnumber(self):
        if self.lock:
            self.lock.acquire()
            nextid=self.next_id.value
            if not self.maxpoints or (nextid < self.maxpoints):
                
                self.next_id.value=self.next_id.value + 1
            else:
                self.lock.release()
                raise # raise an exception to show we're done
            self.lock.release()
            
        else:
            nextid=self.next_id
            if not self.maxpoints or (self.next_id < self.maxpoints):
                self.next_id=self.next_id + 1
            else:
                raise # raise an exception to show we're done

        if not self.maxvalidpoints: # target is number of *valid* points in chain
            UpdateProgress(nextid, self.target)
        return nextid

        
            
    def propose_point(self,last_point):
            okay = False
            while okay == False:
                next_step = [np.random.normal(x, y)
                                  for x, y in zip(last_point, self.variances)]
                okay =  True 
                # check if point is in the given range
                for v, r in zip(next_step, self.inputs['Variables'].values()):
                    if v < min(r['RANGE'][0],r['RANGE'][1]):
                        okay = False 
                    if v > max(r['RANGE'][0],r['RANGE'][1]):
                        okay = False 
            return next_step

    def get_random_point(self):
            return [np.random.uniform(x['RANGE'][0], x['RANGE'][1])
                    for x in self.inputs['Variables'].values()]
          
    def run(self):
        start_time=time.time()
        self.log.debug('Started MCMC scan')
        if len(self.ScanThreads) > 1:
            
            manager = mp.Manager()
            List_all = manager.list()
            List_valid = manager.list() 
            self.lock=manager.Lock()
            #self.next_id=mp.Value('i',0) ## next_id is an integer (hence 'i') 
            self.next_id=mp.Value('i',self.number_of_prior_points) ## next_id is an integer (hence 'i') 
            #self.valid_points=mp.Value('i',0)
            self.valid_points=mp.Value('i',self.number_of_prior_valid)
                
            processes = [mp.Process(target=thread.run,args=(List_all, List_valid)) for thread in self.ScanThreads]
            for p in processes:
                p.start()
            for p in processes:
                p.join()
            self.all_points_run=list(List_all)
            self.valid_points_run=list(List_valid)
            
        else:
            self.all_points_run=[]
            self.valid_points_run=[]
            self.next_id=self.number_of_prior_points
            self.valid_points=self.number_of_prior_valid
            self.ScanThreads[0].run(self.all_points_run,self.valid_points_run)
        self.RunManager.all_points=self.all_points_run
        self.RunManager.valid_points=self.valid_points_run
        end_time=time.time()
        self.log.debug('Finished MCMC scan after %.3f s' %(end_time-start_time))
        UpdateProgress(self.target, self.target, message='',finish=True)

        ## merge output files if that is asked for
        if self.inputs['Setup']['Cores'] > 1:
            if 'Merge Results' in self.inputs['Setup']:
                if eval(self.inputs['Setup']['Merge Results']):
                    merged_file=os.path.join(self.resultsdir,'MCMC_results')
                    if not os.path.exists(merged_file): # don't overwrite previous results!
                        with open(merged_file,'w') as FO:
                            FO.write(self.headerlinestring+'\n')

                    for thread in self.ScanThreads:
                        debug.command_line_log('tail -n +2 '+thread.resultsfile+' >> '+merged_file,self.log)
                        os.remove(thread.resultsfile)
        
        
        
        
        

    #def postprocess(self,Point, observables, outslha,temp_dir,log, lock=None):
    def postprocess(self,Point, observables, data_point,temp_dir,log, lock=None):
        """ return the likelihood; we won't get this far if the point failed to be generated """
        """ calculate likelihood """

        res=0.0 # for the MCMC we need likelihoods, not logL



        #if self.use_loglike:
        #    try:
        #        res = math.exp(self.get_likelihood(observables))
        #    except:
        #        pass
        #else:
        try:
                res = self.get_likelihood(observables)

        except:
                res = 0.0

        #print(f"Got likelihood: {res}")
        return res
        '''
        lh = 1.
        
        for  v, f in zip(observables,self.likelihoodfns):
            lh = lh * f(float(v))
        
        return lh
        '''
        
        


class MCMC_onecore():
    def __init__(self,ParentScan,CoreNumber):
        
        self.ParentScan=ParentScan
        self.runner=ParentScan.RunManager.CoreRuns[CoreNumber]

        self.last_point = []
        self.likelihood_last_step = 0
        #self.steps = 1
        #self.points_total = 1
        self.converged = False

        # set up results file

        self.resultsfile=os.path.join(ParentScan.resultsdir,'MCMC_run_'+str(CoreNumber))

        
        with open(self.resultsfile,'w') as FO: # even if we resume, we restart here 
            FO.write(ParentScan.headerlinestring+'\n')


    def add_valid_point(self):
        if self.ParentScan.lock:
            self.ParentScan.lock.acquire()
            
            if not self.ParentScan.maxvalidpoints or (self.ParentScan.valid_points.value < self.ParentScan.maxvalidpoints):
                self.ParentScan.valid_points.value=self.ParentScan.valid_points.value+1

            else:
                self.ParentScan.lock.release()
                raise ## raises an error indicating we're done
            self.ParentScan.lock.release()
            
        else:
            # single core mode
            if not self.ParentScan.maxvalidpoints or (self.ParentScan.valid_points < self.ParentScan.maxvalidpoints):
                self.ParentScan.valid_points=self.ParentScan.valid_points+1

        if self.ParentScan.maxvalidpoints is not None:
            UpdateProgress(self.ParentScan.valid_points.value,self.ParentScan.maxvalidpoints)
                
    def run(self,all_list,valid_list):

        self.runner.log.debug('Starting MCMC scan')
       
        # Initialise the random seed in each thread, otherwise we get identical points!!
        np.random.seed()
        # start with random point
        found_start=False


        # For resume!!
        if self.last_point != [] and self.likelihood_last_step > 0:
            found_start = True
            self.runner.log.info("Starting scan at point "+str(self.last_point)+" with likelihood "+str(self.likelihood_last_step))
        
        while not found_start:
            
            try:
                tempID=self.ParentScan.get_next_pointnumber()
            except:
                ## have already exhausted the number of points to try
                return
            
            point = self.ParentScan.get_random_point()
            
            try:
                #outlh=self.runner.run_point(point,tempID,lock=self.ParentScan.lock,make_point_container=True)
                outlh=self.runner.run_point(point,tempID,lock=self.ParentScan.lock)
                """
                if self.ParentScan.RunManager.settings.store_points: 
                    tpoint=DataPoint(self.runner.last_point.ID,self.runner.last_point.inputs)
                    tpoint.outputs=self.runner.last_point.outputs
                    tpoint.obs_dict=self.runner.last_point.obs_dict
                    all_list.append(tpoint) 
                """
                #if self.ParentScan.RunManager.settings.store_points:  
                #    valid_list.append(self.runner.last_point) 
                self.likelihood_last_step = outlh
                self.add_valid_point()
                
                found_start=True
                self.last_point=point
            except:
                pass
            #if self.ParentScan.RunManager.settings.store_points: 
            #    all_list.append(self.runner.last_point) 

        # have a valid starting point, now we can engage the general algorithm
        self.reached_end=False
        while not self.reached_end:
            
            # First check in case one of the other threads just completed the list
            if self.ParentScan.maxvalidpoints is not None and (self.ParentScan.valid_points.value >= self.ParentScan.maxvalidpoints):
                self.reached_end = True
                break

            try:
                tempID=self.ParentScan.get_next_pointnumber()
            except:
                self.reached_end=True
                return
            point=self.ParentScan.propose_point(self.last_point)
           
            try:
                #outlh=self.runner.run_point(point,tempID,self.ParentScan.lock,make_point_container=True) 
                 outlh=self.runner.run_point(point,tempID,self.ParentScan.lock)   
            except:
                continue

            """
            if self.ParentScan.RunManager.settings.store_points: 
                tpoint=DataPoint(self.runner.last_point.ID,self.runner.last_point.inputs)
                tpoint.outputs=self.runner.last_point.outputs
                tpoint.obs_dict=self.runner.last_point.obs_dict
                all_list.append(tpoint) 
                #all_list.append(self.runner.last_point)   
            """
            # check whether to keep the point
            #print(f"checking {outlh} for {point} against {self.likelihood_last_step} for {self.last_point}")
            if (outlh > self.likelihood_last_step) or (np.random.uniform(0,1)*self.likelihood_last_step < outlh):
                self.likelihood_last_step = outlh
                #if self.ParentScan.RunManager.settings.store_points: 
                #    valid_list.append(self.runner.last_point)      
                """
                if self.ParentScan.RunManager.settings.store_points: 
                    valid_list.append(tpoint)
                """
                self.last_point=point
                # write the output to our file
                outline=[str(tempID)]
                outline.extend(map(str,point))
                outline.extend(map(str,self.runner.last_observables))
                outline.append(str(outlh))
                with open(self.resultsfile,'a') as FO:
                    FO.write(','.join(outline)+'\n')
                
                try:
                    self.add_valid_point()
                    
                except:
                    # have reached the end
                    print("Reached end of scan!")
                    self.reached_end=True
                    return

            """ # this is redundant, right?!!     
            else:
                continue
            """

            
        return

            
         

       
