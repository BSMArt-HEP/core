"""
Diver scan

:author: M. Goodsell




Some code adapted from PyMultiNest

Uses the Diver library from:
https://diver.hepforge.org

You will have to install that library yourself, there is no cmake and, at least on my machine, at least one of the
parameters in the makefile needed modifying

Requires:
	Diver (built with MPI if required)
	[Optional] mpi4py (for MPI operation)


For variables we require a range, e.g.:

.. code-block:: json

    "Variables": {
    	"lambda": { "RANGE": [0.01,1.0] }
    	},


Various possibilities for likelihood functions are provided: GAUSS, BIAS, SIGMA, USER, 
or "OFF" for observables that should not be taken into account for the likelihood calculation 
(but the user wants to retain them in the summary file for other purposes). These are set by providing e.g.:

.. code-block:: json

    "Observables": {
    	"mh" : { "SLHA": ["MASS", [25]], "SCALING": "GAUSS",
    		"MEAN": 125.0,
    		"VARIANCE":3.0
    		  }

or "SCALING":"BIAS", or "SCALING":"OFF"

It is not necessary to specify the scaling explicitly: gaussian is taken by default

If you want to use your own likelihood function then add it as a tool, set the scaling to OFF for all observables except for the output of 
that tool, which you can set to "USER"

NB we compute LogLikelihoods here and subtract them all from zero to get a MinusLogLikelihood.



We take a uniform prior on all variables. This could be easily modified if there are convenient choices.

The Diver output files are stored in 
RunName/Diver

It is up to you to read them!!!

The raw points will as usual be retained in the Spectrum_Files directory if so desired, and 
merged if "Merge Results" is True:

.. code-block:: text

    "Setup" : {
    	"Type" : "Diver",
    	"Merge Results" : "True",
    	"DiverLib": '/path/to/Diver/lib'
    	...
    }


NB DiverLib just points to the lib directory. On linux machines we will use libdiver.so


Can run either with or without MPI:

For MPI runs, you should run 

mpiexec -n <ncores> BSMArt <json file> 

or, slightly safer:

mpiexec -n <ncores> python3 -m mpi4py BSMArt <json file>



We use the routine cdiver from cwrapper.f90:

.. code-block:: fortran

    subroutine cdiver(minusloglike_in, nPar, lowerbounds, upperbounds, path, nDerived, nDiscrete, discrete,  &
    					  partitionDiscrete, maxciv, maxgen, NP, nF, F, Cr, lambda, current, expon,              &
    					  bndry, jDE, lambdajDE,  convthresh, convsteps, removeDuplicates, doBayesian, prior_in, &
    					  maxNodePop, Ztolerance, savecount, resume, outputSamples, init_population_strategy,    &
    					  discard_unfit_points, max_initialisation_attempts, max_acceptable_value, seed,         &
    					  context, verbose) bind(c)

    	use iso_c_binding, only: c_int, c_bool, c_double, c_char, c_funptr, c_ptr, C_NULL_CHAR
    	use de, only: diver

    	integer, parameter :: maxpathlen = 300

    	type(c_funptr),  intent(in), value :: minusloglike_in, prior_in
    	type(c_ptr),     intent(inout)     :: context
    	integer(c_int),  intent(in), value :: nPar, nDerived, nDiscrete, maxciv, maxgen, NP, nF, bndry, convsteps, savecount, verbose
    	integer(c_int),  intent(in), value :: init_population_strategy, max_initialisation_attempts, seed
    	integer(c_int),  intent(in), target:: discrete(nDiscrete)
    	logical(c_bool), intent(in), value :: partitionDiscrete, current, expon, jDE, lambdajDE, removeDuplicates, doBayesian, resume
    	logical(c_bool), intent(in), value :: outputSamples, discard_unfit_points
    	real(c_double),  intent(in), value :: Cr, lambda, convthresh, maxNodePop, Ztolerance, max_acceptable_value
    	real(c_double),  intent(in)        :: lowerbounds(nPar), upperbounds(nPar), F(nF)
    	character(kind=c_char,len=1), dimension(maxpathlen), intent(in) :: path


There are therefore loads of obscure parameters to deal with. Gah.


"""

__meta__ = {
    "name": "Diver",
    "requires": ["numpy","mpi4py"],
    "settings": {
        "DiverLib": "Path to Diver library",
        "RunName": "Name of the run",
        "Convergence Threshold": "Float",
        "Max Civilisations": "Int",
        "Verbose": "Int",
        "Number of parameters": "Int",
        "NP": "Int"
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

def _load_library(libpath,libname,log):
	libname = {
		'darwin' : libname[len('lib'):],
		'win32'  : libname + '.dll',
		'cygwin' : libname + '.dll',
	}.get(sys.platform, libname + '.so')

	try:
		if sys.platform == 'darwin':
			libname = find_library(libname)
		if libpath == '':
			return cdll.LoadLibrary(libname)
		else:
			return cdll.LoadLibrary(os.path.join(libpath,libname))
	except OSError as e:
		message = str(e)
		if ('%s: cannot open shared object file: No such file or directory' % libname) in message:

			log.critical('Could not load library "%s"' % libname)
			log.critical('You have to build it first,')
			log.critical('and tell BSMArt where it is via Setup -> DiverLib')
			log.critical('Or add it to your LD_LIBRARY_PATH environment variable')
			
		elif message.endswith('cannot open shared object file: No such file or directory'):

			log.critical('Could not load library: %s' % message.split(':')[0])
			log.critical('You have to build Diver,')
			log.critical('and tell BSMArt where it is via Setup -> DiverLib')
			log.critical('Or add it to your LD_LIBRARY_PATH environment variable')

		elif 'undefined symbol: mpi_' in message:

			log.critical('You tried to compile Diver linked with MPI,')
			log.critical('but now when running, Diver can not find the MPI linked libraries.')

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
	return -math.log(1.+math.exp(-val/s))




def FuncMe(scaling,mean,var):
	tmean = mean
	tvar=var
	if scaling == 'OFF':
		return lambda x : 0.0 # dummy function that returns 0 (-loglikelihood(1) ...)
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
	"""Scanner class for Diver Scans"""
	
	def __init__(self, inputs, log):
		Scan.__init__(self, inputs, log)
		
		## Don't need scalers for Diver -> it takes ranges
		#self.downscalers,self.upscalers = BSMutil.create_scalers(self.inputs,log)

		# Create the likelihood functions based on 'scaling' property of variables, if they have it; otherwise use gaussian likelihood
		self.lowerbounds=[]
		self.upperbounds=[]
		for varname in inputs['Variables']:
			if 'RANGE' in inputs['Variables'][varname]:
				varmax=max(inputs['Variables'][varname]['RANGE'][0],inputs['Variables'][varname]['RANGE'][1])
				varmin=min(inputs['Variables'][varname]['RANGE'][0],inputs['Variables'][varname]['RANGE'][1])
				if varmin > varmax:
					tt=varmax
					varmax=varmin
					varmin=tt

				self.lowerbounds.append(varmin)
				self.upperbounds.append(varmax)
			else:
				self.log.critical('Diver requires RANGE be defined for each variable, none given for '+str(varname))
				raise NameError('No RANGE defined for variable '+str(varname))

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

		if 'DiverLib' not in self.inputs['Setup']:
			self.log.warning('"DiverLib" should be in "Setup" options, giving path to the Diver Library\n Let us hope you have put the path to it in the LD_LIBRARY_PATH')
			self.diverlib=''
		else:
			self.diverlib= self.inputs['Setup']['DiverLib']


		self.diverdir=os.path.join(self.inputs['Setup']['cwd'], self.inputs['Setup']['RunName'],'Diver')

		if self.inputs['MPI']['Rank'] ==0:
			if os.path.exists(self.diverdir):
				shutil.rmtree(self.diverdir)

			os.makedirs(self.diverdir)
		

		

	def postprocess(self,Point, observables, data_point,temp_dir,log, lock=None):
		""" return the likelihood; we won't get this far if the point failed to be generated """
		""" calculate likelihood """
		""" We need - LogLike for Diver """
		lh = 0.
		
		for  v, f in zip(observables,self.likelihoodfns):
			lh = lh - f(float(v))
		return lh
		
	


	def run(self):
		
		self.log.debug('Started Diver scan')
		use_MPI = False
		#lib = _load_library(self.multinestlib,'libmultinest',self.log)
		lib = None

		try:
			lib = _load_library(self.diverlib,'libdiver',self.log)
			#use_MPI = True
		except:
			self.log.debug(e)
			raise SystemExit 
		if self.inputs['MPI']['Size'] > 1:
			try: # detect if run through mpiexec/mpirun
				from mpi4py import MPI
				if MPI.COMM_WORLD.Get_size() > 1: # need parallel capabilities
					self.log.info('MPI running! Core' + str(MPI.COMM_WORLD.Get_rank()))
					use_MPI=True
				else:			
					self.log.info('Not run through mpiexec/mpirun!')
					use_MPI=False
			except ImportError as e:
				self.log.debug(e)
				use_MPI=False
	
		

		n_params = len(self.lowerbounds)
		assert n_params == len(self.upperbounds)

		# check if threads are involved
		# in that case we can't use the signal handler
		is_thread = threading.active_count() > 1

				
		## CFUNCTYPE(return type, arguments)
		minusloglike_type = CFUNCTYPE(c_double,POINTER(c_double),c_int,POINTER(c_int),  c_bool, c_bool, c_void_p)

		#### Use flat prior for now: include possibilities of others later on; usually we convolve
		#### the cube -> Prior -> Point
		#### Here we need loglike to run cube -> upscaler -> run_point
		validvector=True

		fquit = False

		self.call=0
		def minusloglike(cube,tempnps, fcall, fquit, fvalidvector, nullcontext):
			## NB fcall is the number of the call. It needs to be updated here
			## Passed as a pointer so we access its value through fcall[0]
			tcall = int(fcall[0])
			tcall=tcall+1
			fcall[0] = c_int(tcall)
			pp = [float(cube[n]) for n in range(tempnps)]
			fquit=c_bool(False)
			try:
				outlh=self.RunManager.CoreRuns[0].run_point(pp)
				fvalidvector=c_bool(True)
			except: ## On exception we return a high value rather than failing. 
				outlh=1e10 
				fvalidvector=c_bool(False)
			return outlh

		
		
		if not is_thread:
			prev_handler = signal.signal(signal.SIGINT, interrupt_handler)
	
		# to avoid garbage collection of these ctypes, which leads to NULLs
		# we need to make local copies here that are not thrown away

		outputfiles_basename = os.path.join(self.diverdir,'diverout')
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


		
		WrappedType = c_double * n_params
		zeroes=[0]*n_params
		
		nPar = n_params
		#NP=10*n_params
		#maxciv = 10
		
		if 'Convergence Threshold' in self.inputs['Setup']:
			self.convthresh = float(self.inputs['Setup']['Convergence Threshold'])
		else:
			self.convthresh = 0.001   ### This seems very low; it's the default but seems to never terminate for me.

		if 'Max Civilisations' in self.inputs['Setup']:
			self.maxcivs = int(self.inputs['Setup']['Max Civilisations'])
		else:
			self.maxcivs = 10   ### The default in the code was 10000. This seems to never terminate for me.

		if 'Verbose' in self.inputs['Setup']:
			self.verbose = int(self.inputs['Setup']['Verbose'])
		else:
			self.verbose = 1   ### Could be 1,2,3

		if 'Number of parameters' in self.inputs['Setup']:
			self.NP = int(self.inputs['Setup']['Number of parameters'])
		elif 'NP' in self.inputs['Setup']:
			self.NP = int(self.inputs['Setup']['Number of parameters'])
		else:
			self.NP = 10*n_params   ### Default choice in the code. 

		

		args_and_types = [[ minusloglike, minusloglike_type], #minusloglike
						[nPar, c_int], # nPar
						[WrappedType(*self.lowerbounds),lambda x: x], # lowerbounds
						[WrappedType(*self.upperbounds),lambda x: x], # upperbounds
						[sb, lambda x: x], # path
						[0, c_int], # nDerived
						[0, c_int], # nDiscrete
						[None, lambda x: x], # discrete
						[ False, c_bool], # partitionDiscrete
						[self.maxcivs, c_int], #  maxciv
						[300, c_int], # maxgen
						[self.NP, c_int], # NP
						[0, c_int], # nF
						[None, lambda x: x], # F(nF)
 						[0.9, c_double], # Cr
						[0.0, c_double], # lambda
						[False, c_bool], # current
						[False, c_bool], # expon
						[1, c_int], #bdnry, can be 0 (none),1 (brick wall),2 (random) or 3 (reflection)
						[True, c_bool], # jDE
						[True, c_bool], #lambdajDE
						[self.convthresh, c_double], #convthresh
						[10, c_int], #convsteps
						[True, c_bool], # removeDuplicates
						[False, c_bool], # doBayesian
						[None, lambda x: x], # prior_in
						[1.0, c_double], # maxNodePop
						[0.01, c_double], # Ztolerance
						[1, c_int], # savecount
						[False, c_bool], # resume
						[True, c_bool], # outputSamples
						[0, c_int], # init_population_strategy
						[False, c_bool], # discard_unfit_points
						[self.maxcivs, c_int], # max_initialisation_attempts
						[1e6, c_double], # max_acceptable_value
						[-1, c_int], # seed
						[None, lambda x: x], # context
						[self.verbose, c_int] # verbose
					  	]

		args_converted = [arg[1](arg[0]) for arg in args_and_types]
		
				
		if use_MPI and lib is not None:
			lib.cdiver(*args_converted)
			# wait for all processes to return
			# Otherwise rank=0 is still writing output files
			# while the others already try to read them.
			from mpi4py import MPI
			MPI.COMM_WORLD.Barrier()
		else:
			lib.cdiver(*args_converted)
		if not is_thread:
			signal.signal(signal.SIGINT, prev_handler)
		
		
