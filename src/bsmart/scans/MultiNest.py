"""
MultiNest scan


Written by M. Goodsell

Shamelessly adapted from parts of PyMultiNest by Johannes Buchner
https://johannesbuchner.github.io/PyMultiNest

Based on MultiNest (and using the library) from:
https://github.com/farhanferoz/MultiNest

Requires:
	MultiNest (built with MPI if required)
	[Optional] mpi4py (for MPI operation)


For variables we require a range, e.g.:

"Variables": {
	"lambda": { "RANGE": [0.01,1.0] }
	},

Various possibilities for likelihood functions are provided: GAUSS, BIAS, SIGMA, USER, LOGUSER
or "OFF" for observables that should not be taken into account for the likelihood calculation 
(but the user wants to retain them in the summary file for other purposes). These are set by providing e.g.:

"Observables": {
	"mh" : { "SLHA": ["MASS", [25]], "SCALING": "GAUSS",
		"MEAN": 125.0,
		"VARIANCE":3.0
		  }

or "SCALING":"BIAS", or "SCALING":"OFF"

It is not necessary to specify the scaling explicitly: gaussian is taken by default

We take a uniform prior on all variables. This could be easily modified if there are convenient choices.

The MultiNest output files are stored in 
RunName/MultiNest

It is up to you to read them!!!

The raw points will as usual be retained in the Spectrum_Files directory if so desired, and 
merged if "Merge Results" is True:

"Setup" : {
	"Type" : "MultiNest",
	"Merge Results" : "True",
	"MultiNestLib": '/path/to/MultiNest/lib'
	...
}

Note that 


Can run either with or without MPI:

For MPI runs, you should run 

mpiexec -n <ncores> BSMArt <json file> 

or, slightly safer:

mpiexec -n <ncores> python3 -m mpi4py BSMArt <json file>


"""

__meta__ = {
    "name": "MultiNest",
    "requires": ["numpy","mpi4py"],
    "settings": {
        "MultiNestLib": "Path to MultiNest library",
        "RunName": "String",
        "Live Points": "Int",
        "Points": "Int",
        "Sampling Efficiency": "Float",
        "Evidence Tolerance": "Float"
    }
}



from bsmart.core import Scan as Scan

import numpy as np
import math
from bsmart import debug


import os
import shutil

import bsmart.BSMutil as BSMutil

from bsmart.BSMlikelihood import MakeLikelihoods

from ctypes import cdll
from ctypes.util import find_library 
import sys, os, threading

"""
try: # detect if run through mpiexec/mpirun
        from mpi4py import MPI
except:
        pass
"""


def _load_library(libpath,libname,log):
	libname = {
		#'darwin' : libname[len('lib'):],
		'win32'  : libname + '.dll',
		'cygwin' : libname + '.dll',
	}.get(sys.platform, libname + '.so')

	try:
		#if sys.platform == 'darwin':
		#	libname = find_library(libname)
		if libpath == '':
			return cdll.LoadLibrary(libname)
		else:
			return cdll.LoadLibrary(os.path.join(libpath,libname))
	except OSError as e:
		message = str(e)
		if ('%s: cannot open shared object file: No such file or directory' % libname) in message:

			log.critical('Could not load MultiNest library "%s"' % libname)
			log.critical('You have to build it first,')
			log.critical('and tell BSMArt where it is via Setup -> MultiNestLib')
			log.critical('Or add it to your LD_LIBRARY_PATH environment variable')
			
		elif message.endswith('cannot open shared object file: No such file or directory'):

			log.critical('Could not load MultiNest library: %s' % message.split(':')[0])
			log.critical('You have to build MultiNest,')
			log.critical('and tell BSMArt where it is via Setup -> MultiNestLib')
			log.critical('Or add it to your LD_LIBRARY_PATH environment variable')

		elif 'undefined symbol: mpi_' in message:

			log.critical('You tried to compile MultiNest linked with MPI,')
			log.critical('but now when running, MultiNest can not find the MPI linked libraries.')

		# the next if is useless because we can not catch symbol lookup errors (the executable crashes)
		# but it is still there as documentation.
		elif 'symbol lookup error' in message and 'mpi' in message:

			log.critical('You are trying to get MPI to run, but MPI failed to load.')
			log.critical('Specifically, mpi symbols are missing in the executable.')
			log.critical('Let me know if this is a problem of running python or a compilation problem.')
		# what if built with MPI, but don't have MPI
		else:
			log.critical('problem:', e)
		raise SystemExit

from ctypes import *
from numpy.ctypeslib import as_array
import signal, sys
#import inspect

def interrupt_handler(recvsignal, frame):
	sys.stderr.write('ERROR: Interrupt received: Terminating\n')
	os.kill(os.getpid(), signal.SIGTERM)



"""
#### We're supposed to use log likelihoods


# bias function to enhance likelihoods  
def bias(val, m, s):
	return s*math.log(val/m)

# gauss likelihood
def gauss(val, m, s):
	#print('%f %f %f' %(val,m,s))
	return -(val-m)**2/(2.*(s**2))

# Log likelihood
#def logL(val, m, s):
#	return 1./(1.+(val-m)**2/(2.*s)**2)

## Sigmoid for upper or lower bounds. The variance 's' controls how sharply it should turn off
def Sigmoid(val, s): 
	return - math.log(1.+math.exp(-val/s))




def FuncMe(scaling,mean,var):
	tmean = mean
	tvar=var
	if scaling == 'OFF':
		return lambda x : 0.0 # dummy function that returns 1
	#elif scaling == 'LOG':
#	return lambda x : logL(x, tmean, tvar)
	elif scaling == 'UPPER':
		return lambda x : Sigmoid(tmean-x, tvar)
	elif scaling == 'LOWER':
		return lambda x : Sigmoid(x-tmean, tvar)
	elif scaling == 'BIAS':
		return lambda x : bias(x,tmean,tvar)
	elif scaling == 'USER': #### User supplied likelihood/pval, e.g. from MicrOmegas
		return lambda x : x
	else:
		return lambda x : gauss(x, tmean, tvar)

"""


class NewScan(Scan):
	"""Scanner class for MultiNest Scans"""
	
	def __init__(self, inputs, log):
		Scan.__init__(self, inputs, log)
		

		self.downscalers,self.upscalers = BSMutil.create_scalers(self.inputs,log)

		# Create the likelihood functions based on 'scaling' property of variables, if they have it; otherwise use gaussian likelihood    
		
		self.likelihoodfns=MakeLikelihoods(self.inputs['Observables'],True) ## use log likelihood
		"""
		self.likelihoodfns=[]
		for o in self.inputs['Observables']:
			#self.header.append(o)
			obs=self.inputs['Observables'][o]
			scaling='GAUSS'
			if 'SCALING' in obs:
				scaling=obs['SCALING']
			if 'MEAN' in obs:
				tmean = obs['MEAN']
			else:
				tmean=0.0
			if 'VARIANCE' in obs:
				tvar = obs['VARIANCE']
			else:
				tvar = 0.0
				
			self.likelihoodfns.append(FuncMe(scaling,tmean,tvar))   
		"""
		# set up results directory and files

		if 'MultiNestLib' not in self.inputs['Setup']:
			self.log.warning('"MultiNestLib" should be in "Setup" options, giving path to the MultiNest Library\n Let us hope you have put the path to it in the LD_LIBRARY_PATH')
			self.multinestlib=''
		else:
			self.multinestlib= self.inputs['Setup']['MultiNestLib']


		self.multinestdir=os.path.join(self.inputs['Setup']['cwd'], self.inputs['Setup']['RunName'],'MultiNest')

		if self.inputs['MPI']['Rank'] ==0:
			if os.path.exists(self.multinestdir):
				shutil.rmtree(self.multinestdir)

			os.makedirs(self.multinestdir)
		
		

	def postprocess(self,Point, observables, data_point,temp_dir,log, lock=None):
		""" return the likelihood; we won't get this far if the point failed to be generated """
		""" calculate likelihood """
		lh = 0.
		
		for  v, f in zip(observables,self.likelihoodfns):
			lh = lh + f(float(v))
		

		return lh
		
	


	def run(self):
		n_params = None
		n_clustering_params = None
		wrapped_params = None
		importance_nested_sampling = True
		multimodal = True
		const_efficiency_mode = False
		#n_live_points = 400
		if 'Live Points' in self.inputs['Setup']:
			n_live_points = int(self.inputs['Setup']['Live Points'])
		elif 'Points' in self.inputs['Setup']:
			n_live_points = int(self.inputs['Setup']['Points'])
		else:
			n_live_points = 400

		
		evidence_tolerance = 0.5
		sampling_efficiency = 0.8
		n_iter_before_update = 100
		null_log_evidence = -1e90
		max_modes = 100
		mode_tolerance = -1e90
		#outputfiles_basename = "chains/1-"
		seed = -1
		verbose = False
		resume = True
		context = 0
		write_output = True
		log_zero = -1e100 
		max_iter = 0
		init_MPI = False
		dump_callback = None
		use_MPI = False
		self.log.debug('Started MultiNest scan')
		lib = _load_library(self.multinestlib,'libmultinest',self.log)
		lib_mpi = None

		if self.inputs['MPI']['Size'] > 1:
			try: # detect if run through mpiexec/mpirun
				from mpi4py import MPI
				if MPI.COMM_WORLD.Get_size() > 1: # need parallel capabilities
					self.log.info('MPI running! Core' + str(MPI.COMM_WORLD.Get_rank()))
					lib_mpi = _load_library(self.multinestlib,'libmultinest_mpi',self.log)
					use_MPI = True
				else:			
					self.log.info('Not run through mpiexec/mpirun!')
					use_MPI=False
			except ImportError as e:
				self.log.debug(e)
				use_MPI=False
	


		n_params = len(self.upscalers)
		n_dims=n_params
		n_clustering_params = n_dims
		wrapped_params = [0] * n_dims

		WrappedType = c_int * len(wrapped_params)
		wraps = WrappedType(*wrapped_params)

		if 'Sampling Efficiency' in self.inputs['Setup']: ### have already set it to 0.8 by default
		#sampling_efficiency = 'parameter'
			se=self.inputs['Setup']['Sampling Efficiency']
			try:
				val = float(se)
				sampling_efficiency = val
			except:
				if se == 'parameter':
					sampling_efficiency = 0.8
				elif se == 'model':
					sampling_efficiency = 0.3
				else:
					sampling_efficiency = 0.8
		
		if 'Evidence Tolerance' in self.inputs['Setup']: ### have already set it to 0.5 by default
			try:
				evidence_tolerance = float(self.inputs['Setup']['Evidence Tolerance'])
			except:
				evidence_tolerance = 0.5
			
		loglike_type = CFUNCTYPE(c_double, POINTER(c_double),c_int, c_int, c_double, c_void_p)
		
		dumper_type  = CFUNCTYPE(c_void_p, c_int, c_int, c_int,POINTER(c_double),POINTER(c_double),POINTER(c_double),
						c_double,c_double,c_double,c_void_p)

		# check if threads are involved
		# in that case we can't use the signal handler
		
		is_thread = threading.active_count() > 1
	

		#### Use flat prior for now: include possibilities of others later on; usually we convolve
		#### the cube -> Prior -> Point
		#### Here we need loglike to run cube -> upscaler -> run_point
		def loglike(cube, ndim, nparams, lnew, nullcontext):
			point = [f(x) for x,f in zip(cube,self.upscalers)]
			try:
				outlh=self.RunManager.CoreRuns[0].run_point(point)
			except: ## On exception we return a low value rather than failing. 
				outlh=-1e100 
			return outlh


		def dumper(nSamples,nlive,nPar,
			   physLive,posterior,paramConstr,
			   maxLogLike,logZ,logZerr,nullcontext):
			if dump_callback:
			# It's not clear to me what the desired PyMultiNest dumper callback
			# syntax is... but this should pass back the right numpy arrays,
			# without copies. Untested!
				pc =  as_array(paramConstr,shape=(nPar,4))
			
				dump_callback(nSamples,nlive,nPar,
					as_array(physLive,shape=(nPar+1,nlive)).T,
					as_array(posterior,shape=(nPar+2,nSamples)).T, 
					(pc[:,0],pc[:,1],pc[:,2],pc[:,3]), # (mean,std,bestfit,map)
					maxLogLike,logZ,logZerr, 0)
		if not is_thread:
			prev_handler = signal.signal(signal.SIGINT, interrupt_handler)
	
		# to avoid garbage collection of these ctypes, which leads to NULLs
		# we need to make local copies here that are not thrown away

		### Using temp_dir here proves that the writing is only in rank 0 process, so
		### therefore it is properly doing MPI:
		#outputfiles_basename = self.temp_dir +'/'
		## But we actually want the results in a proper results directory:
		outputfiles_basename = self.multinestdir +'/'
		s = outputfiles_basename.encode()

		if len(outputfiles_basename + "post_equal_weights.txt") < 100:
			# 100 characters work for all MultiNest versions
			sb = create_string_buffer(s, 100)
		elif len(outputfiles_basename + "post_equal_weights.txt") > 1000:
			# file name is too long
			raise ValueError("Filenames must be less than 1000 characters long.")
		else:
			# try 1000 character length file name (this may cause MultiNest failure
			# or filename truncation if not using the latest MultiNest version)
			sb = create_string_buffer(s, 1000)

		argtypes = [c_bool, c_bool, c_bool, 
			c_int, c_double, c_double, 
			c_int, c_int, c_int, c_int, 
			c_int, c_double, 
			lambda x: x, c_int, lambda x: x, 
			c_bool, c_bool, c_bool, c_bool, 
			c_double, c_int, loglike_type, dumper_type, c_int
			]
		args = [importance_nested_sampling, multimodal, const_efficiency_mode,
			n_live_points, evidence_tolerance, sampling_efficiency, 
			n_dims, n_params, n_clustering_params, max_modes, 
			n_iter_before_update, mode_tolerance, 
			sb, seed, wraps,
			verbose, resume, write_output, init_MPI,
			log_zero, max_iter, loglike, dumper, context]
		args_converted = [converter(v) for v, converter in zip(args, argtypes)]
		if use_MPI and lib_mpi is not None:
			from mpi4py import MPI
			lib_mpi.run(*args_converted)
			# wait for all processes to return
			# Otherwise rank=0 is still writing output files
			# while the others already try to read them.
			MPI.COMM_WORLD.Barrier()
		else:
			lib.run(*args_converted)
		if not is_thread:
			signal.signal(signal.SIGINT, prev_handler)
		#assert len(args) == len(argtypes) # to make sure stuff is still here
			#   outlh=self.runner.run_point(point,tempID,lock=self.ParentScan.lock,make_point_container=True)

		
