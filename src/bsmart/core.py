""" 
----------------------------------------------------------
 BSMArt Core module to define the basic Scan class
----------------------------------------------------------
 All scans inherit the basic features here.


Useful settings:
	'StoreEverything': will store all output from all codes in separate directories in All_Outputs
	'StoreAllPoints': Keep all points. Default is to concatenate them into one file stored in Results
	'StoreSeparateFiles': Separate the resulting spectrum files, stored in Results
	'Store In Memory': Store *all* the points in memory during running; useful for plotting and for user-defined scans (defualt False)
	'Output File': name for output file collecting spectra/csv/tabbed output
	'Spectrum File': name for the spectrum file. Only really used internally
	'Merge Results' : Merge the spectrum/csv files if MPI operation is used, also for MCMC Results.
	'csv' : Store results in a csv file
"""



import os
import shutil
import time
import datetime

import math
import itertools

import importlib
from six import string_types

from bsmart.HEPRun import HEPRun, RunSettings, CoreRunner, HepTool
from collections import OrderedDict

try:
	from bsmart import zslha
except:
	raise SystemExit('Failed to import zslha')


from bsmart import debug   


# raw string is allowed in old versions of python
BSMArt_references=r'''
@article{Goodsell:2023iac,
    author = "Goodsell, Mark D. and Joury, Ari",
    title = "{BSMArt: Simple and fast parameter space scans}",
    eprint = "2301.01154",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    doi = "10.1016/j.cpc.2023.109057",
    journal = "Comput. Phys. Commun.",
    volume = "297",
    pages = "109057",
    year = "2024"
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

@article{Faraggi:2023jzm,
    author = "Faraggi, Alon E. and Goodsell, Mark D.",
    title = "{$M_W$ in string derived Z$^\prime $ models}",
    eprint = "2312.13411",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    doi = "10.1140/epjc/s10052-024-12900-9",
    journal = "Eur. Phys. J. C",
    volume = "84",
    number = "6",
    pages = "589",
    year = "2024"
}

'''



class Scan:
	""" Parent Scan class from which all child NewScan classes are born """
	def __init__(self, inputs, log): 
		log.info('Initialise scan: %s' % inputs['Setup']['RunName'])
		#self.show_credits() # optional message to display credits for the scan!       
		self.inputs = inputs # store all inputs. NB These may be modified a little during the creation of self.runsettings
		# MAINpath is the working directory for the scan
		self.MAINpath=inputs['Setup']['MAINpath']
		self.temp_dir=inputs['Setup']['TempDir']
		
		self.log=log

		self.citations : str = None   # allow these to be overloaded in the "initialise" phase of the scan
		
		self.resume = eval(inputs['Setup'].get('Resume',"False"))



		# This may throw an exception, which should be caught by whatever initialises the scan.
		self.runsettings=RunSettings(inputs)
		self.varnames = self.runsettings.varnames
		

		'''
		## now setup the run manager
		self.runsettings=RunSettings(inputs['Setup']['RunName'])
		self.runsettings.cores=inputs['Setup']['Cores']
		self.runsettings.work_dir=inputs['Setup']['TempDir']
		self.runsettings.debug=eval(inputs['Setup']['debug'])
		self.runsettings.packagepath=inputs['Setup']['PackagePath']

		self.runsettings.no_slha_output=False

		self.resume = eval(inputs['Setup'].get('Resume',"False"))

		if inputs['MPI']['Size'] > 1:
			self.runsettings.use_MPI=True

		if 'SLHA_Output' in inputs['Setup']:
			self.runsettings.no_slha_output=not eval(str(inputs['Setup']['SLHA_Output']))
		
		if 'Codes' not in inputs:
			self.log.error('No codes found in input file')
			raise SystemExit
		else:
			self.runsettings.codes=inputs['Codes']
		
			### add observables defined within the code to the global observables list if present
			for code in inputs['Codes']:
				if eval(inputs['Codes'][code]['Run']) and 'Observables' in inputs['Codes'][code]:
					for obs in inputs['Codes'][code]['Observables']:
						inputs['Observables'][obs] = inputs['Codes'][code]['Observables'][obs]

		
		if 'StoreInputs' in inputs['Setup']:
			self.runsettings.store_inputs=eval(str(inputs['Setup']['StoreInputs']))
		else:
			self.runsettings.store_inputs=False

		if 'Short' in inputs['Setup']:
			self.runsettings.tabbed_output = eval(str(inputs['Setup']['Short']))
		
		if 'csv' in inputs['Setup']:
			self.runsettings.csv_output = eval(str(inputs['Setup']['csv']))
		if 'csv name' in inputs['Setup']:
			self.runsettings.csv_name=str(inputs['Setup']['csv name'])
		
		if 'Precision' in inputs['Setup']:
			self.runsettings.precision = int(inputs['Setup']['Precision'])
			
		if 'StoreAllPoints' in self.inputs['Setup']:
			self.runsettings.store_outputs=eval(self.inputs['Setup']['StoreAllPoints'])
		else:
			self.runsettings.store_outputs=True
			
		if 'Store In Memory' in self.inputs['Setup']:
			self.runsettings.store_points=eval(self.inputs['Setup']['Store In Memory'])

		
		if 'StoreSeparateFiles' in self.inputs['Setup']:
			self.runsettings.store_separate_files=eval(str(self.inputs['Setup']['StoreSeparateFiles']))
			
		self.runsettings.store_everything=False
		if 'StoreEverything' in self.inputs['Setup']:
			self.runsettings.store_everything=eval(str(self.inputs['Setup']['StoreEverything']))
			if self.runsettings.store_everything:
				self.runsettings.store_outputs=True;
				if inputs['MPI']['Size'] > 1:
					self.runsettings.all_dir=os.path.join(self.inputs['Setup']['cwd'], self.inputs['Setup']['RunName'],'All_Outputs','MPI_'+str(inputs['MPI']['Rank']))
				else:
					self.runsettings.all_dir=os.path.join(self.inputs['Setup']['cwd'], self.inputs['Setup']['RunName'], 'All_Outputs')
		if inputs['MPI']['Size'] > 1:    
			self.runsettings.inputs_dir = os.path.join(self.inputs['Setup']['cwd'], self.inputs['Setup']['RunName'], 'Input_Files','MPI_'+str(inputs['MPI']['Rank']))
			self.runsettings.output_dir = os.path.join(self.inputs['Setup']['cwd'], self.inputs['Setup']['RunName'],'Spectrum_Files','MPI_'+str(inputs['MPI']['Rank']))
		else:           
			self.runsettings.inputs_dir = os.path.join(self.inputs['Setup']['cwd'], self.inputs['Setup']['RunName'], 'Input_Files')
			self.runsettings.output_dir = os.path.join(self.inputs['Setup']['cwd'], self.inputs['Setup']['RunName'], 'Spectrum_Files')
		self.runsettings.output_filename='SLHA_output'
		if 'Output File' in self.inputs['Setup']:
			self.runsettings.output_filename=self.inputs['Setup']['Output File']
		
		if 'Spectrum File' in self.inputs['Setup']:
			self.runsettings.spc_filename=self.inputs['Setup']['Spectrum File']
		else:
			self.runsettings.spc_filename='SLHA_output'
		
		## the input file should be the input for the first code that we run, presumably SPheno. Note we imported the json as an OrderedDict
		firstcode=next(iter(self.inputs['Codes']))
		self.runsettings.input_filename=self.inputs['Codes'][firstcode]['InputFile']
		if 'Variables' in self.inputs:
			self.runsettings.variables=self.inputs['Variables']
			self.varnames = list(self.inputs['Variables'].keys())
			self.runsettings.varnames=self.varnames
		else:
			self.runsettings.variables=None
			self.varnames=None
			self.runsettings.varnames=None

		if 'Derived Variables' in self.inputs:
			"""
				Allow for functions of the variables, which should be provided as a dictionary like
				{"dvar1": "math.exp(mh)", ...}
				These will be treated like variables when writing the input files etc, but won't show up in the csv file
			"""
			self.runsettings.derived_variables = self.inputs['Derived Variables']

		if 'Observables' in self.inputs:
			self.runsettings.observables=self.inputs['Observables']

		if 'Fitting' in self.inputs:
			self.runsettings.fitting_inputs=self.inputs['Fitting']
		'''
			
		# Allow the scan to modify initial runsettings if needed
		self.initialise()


		#if self.runsetting.fitting_inputs is not None and 'Options' in self.runsetting.fitting_inputs and  'Overloads' in self.runsetting.fitting_inputs['Options']:
		#	self.overload_spc = 

		self.doplots=False
		
		if 'Plotting' in self.inputs:
			self.plotting=self.inputs['Plotting']
			self.plots=self.inputs['Plotting']['Plots']
			if self.plotting['Strategy']=='All' or self.plotting['Strategy']=='Valid':
				self.runsettings.store_points_in_memory=True
			if len(self.plots) > 0:
				#import BSMplots
				self.doplots=True
				#BSMplots.say_hello(self)
				self.plotdir=os.path.join(self.inputs['Setup']['cwd'], self.inputs['Setup']['RunName'], 'Plots')

		
		self.make_dirs()
		
		

		self.configure_input_spc()

		# store all valid and invalid points, although the user should decide what to do with them; if self.runsettings.store_points_in_memory is set then this is stored in the run manager object anyway
		# I am not sure that these are ever used -> remove?

		"""
		self.all_data = []
		self.all_valid=[]
		self.all_invalid=[]

		# I am not sure that these are ever used -> remove?
		self.pointID=0  ### so that each point gets a unique ID; we currently let the run manager take care of this
		"""


		# now initialise the run manager, where the work of running the codes is done
		
		if self.runsettings.no_slha_output:
			self.RunManager=HEPRun(self.runsettings, self.write_lh_file, self.postprocess_no_slha,self.log)
		else:
			self.RunManager=HEPRun(self.runsettings, self.write_lh_file, self.postprocess,self.log)


		self.write_citations()
		
		### END __init__
		


	def write_citations(self):
		if self.inputs['MPI']['Size'] > 1: ## only do this once!
			return 
		
		with open(os.path.join(self.inputs['Setup']['cwd'], self.inputs['Setup']['RunName'],'references.bib'),'w') as BO:
			BO.write(BSMArt_references)
			if self.citations is not None:
				BO.write(self.citations+'\n')

			for tool in self.RunManager.tools:
				if tool.citations is not None:
					BO.write(tool.citations+'\n')
			

	def initialise(self):
		""" method to allow the user scan to overload run settings etc during the initialisation process """
		pass
	
	def configure_input_spc(self):
		""" 
		We store a template input spectrum file and use zslha to write our input files 
		for different codes by default. However, the user can overload write_lh_file 
		or setlhvals as necessary if a different format is needed 
		"""
		self.log.info('Parse Settings')
		self.log.info('Reading bare spc: '+str(self.runsettings.input_filename))

		## in case there is no template and we want to use json input
		try:
			self.inputspc= zslha.read(self.runsettings.input_filename,verbose=False)
		except Exception as e:
			self.log.error('Problem reading slha file '+self.runsettings.input_filename+', ' + str(e))
			self.inputspc = zslha.SLHA()
		
		# for compatibility with json input
		if 'Blocks' in self.inputs:
			try:
				for current_block in self.inputs['Blocks']:
					
					if current_block not in self.inputspc.blocks:
						self.log.debug(current_block + ' not present ' )
						self.inputspc.blocks[current_block]={}
						
					for key in self.inputs['Blocks'][current_block].keys():
						self.log.debug('putting value of ' + key + ' to ' +str(self.inputs['Blocks'][current_block][key]))
						
						self.inputspc.blocks[current_block][key] = self.inputs['Blocks'][current_block][key]
						self.log.debug('Set ' + str(self.inputspc.blocks[current_block][key]))
			except Exception as e:
				self.log.debug('Exception reading json block input: '+str(e))
		

	def make_dirs(self):
		""" Create output and inputs directories"""
		
		
		if not self.resume and os.path.exists(self.runsettings.output_dir):
				shutil.rmtree(self.runsettings.output_dir)

		if not os.path.exists(self.runsettings.output_dir):	
			os.makedirs(self.runsettings.output_dir)

		if self.runsettings.store_everything:
			if not self.resume and os.path.exists(self.runsettings.all_dir):
				shutil.rmtree(self.runsettings.all_dir)
			
			if not os.path.exists(self.runsettings.all_dir):
				os.makedirs(self.runsettings.all_dir)

		
		if not self.resume and os.path.exists(self.runsettings.inputs_dir):
			shutil.rmtree(self.runsettings.inputs_dir)    
			
		if self.runsettings.store_inputs and not os.path.exists(self.runsettings.inputs_dir):
			os.makedirs(self.runsettings.inputs_dir)

		if self.doplots:
			if self.inputs['MPI']['Rank'] == 0: ### Only one set of plots ... 
				if os.path.exists(self.plotdir):
					shutil.rmtree(self.plotdir)
				os.makedirs(self.plotdir)
			
	
	def setlhvals(self,point,values_dict=None):
		###### USER MAY OVERLOAD THIS FUNCTION AS THEY WANT
		
		# If we want to modify this should create a deepcopy:
		# resspc = zslha.SLHA(self.inputspc)
		resspc=self.inputspc
		return resspc
	
	def write_lh_file(self,point, dir, name,values_dict=None,Overloads=None):
		""" write Les Houches input file for given parameter point """

		lhname=os.path.join(dir,name)
		runspc=self.setlhvals(point,values_dict)
		if Overloads is not None:  
			"""
				Allow overloading of block values, e.g. during marginalisation

				Format: a list of form
				[["SPHENOINPUTS",[1],1.0],...]
			"""
			runspc = zslha.SLHA(runspc) #
			for overload in Overloads:
				if len(overload) < 3:
					continue
				block=overload[0].upper()
				blockkey=str(overload[1])[1:-1].replace(" ", "")
				value=overload[2]
				if block in runspc.blocks:
					runspc.blocks[block][blockkey] = value

                            
		if values_dict is None:
			zslha.write_lh(runspc, point, lhname,self.varnames)
		else:
			zslha.write_lh_dict(runspc, values_dict, lhname)




		
	def start_scan(self):
		self.start_time = time.time()
		self.log.info('Running scan %s' % str(self.inputs['Setup']['RunName']))

	def finish_scan(self):
		self.end_time = time.time()
		elapsed=self.end_time-self.start_time
		self.log.info('Scan %s finished after %s ' % (str(self.inputs['Setup']['RunName']),str(datetime.timedelta(seconds=elapsed))))

		if self.inputs['MPI']['Size'] > 1:
			try: # detect if run through mpiexec/mpirun
				from mpi4py import MPI
			except ImportError as e:
				self.log.debug(e)
				
                        
			if 'Merge Results' in self.inputs['Setup']:
				if eval(self.inputs['Setup']['Merge Results']):
					self.log.info('Merging Results')

					#if self.runsettings.tabbed_output:

					if self.runsettings.store_outputs and not self.runsettings.store_separate_files:
						specfiles=MPI.COMM_WORLD.gather(self.RunManager.output_filename)
						if self.inputs['MPI']['Rank'] == 0:
							if 'Output File' in self.inputs['Setup']:
								merge_spectrum_file=os.path.join(self.inputs['Setup']['cwd'], self.inputs['Setup']['RunName'], 'Spectrum_Files',self.inputs['Setup']['Output File'])
							else:
								merge_spectrum_file=os.path.join(self.inputs['Setup']['cwd'], self.inputs['Setup']['RunName'], 'Spectrum_Files','SLHA_Output')
							if os.path.isfile(merge_spectrum_file):
								os.remove(merge_spectrum_file)

							for spec in specfiles:
								if os.path.isfile(spec):
									debug.command_line_log('cat '+spec+' >> '+merge_spectrum_file,self.log)
									os.remove(spec)
									
								

							## Change this so that the plot functions find it:
							self.RunManager.output_filename=merge_spectrum_file

					if self.runsettings.csv_output:
						csvfiles=MPI.COMM_WORLD.gather(self.RunManager.csv_filename)
						#print(csvfiles)
						if self.inputs['MPI']['Rank'] == 0:
							if 'csv name' in self.inputs['Setup']:
								merge_csv = os.path.join(self.inputs['Setup']['cwd'], self.inputs['Setup']['RunName'], 'Spectrum_Files',self.inputs['Setup']['csv name'])
							else:
								if 'Output File' in self.inputs['Setup']:
									merge_csv = os.path.join(self.inputs['Setup']['cwd'], self.inputs['Setup']['RunName'], 'Spectrum_Files',self.inputs['Setup']['Output File']+'.csv')
								else:
									merge_csv = os.path.join(self.inputs['Setup']['cwd'], self.inputs['Setup']['RunName'], 'Spectrum_Files','SLHA_Output.csv')

							## First do the first line of the first file, then cat all of them
							if os.path.isfile(merge_csv):
								os.remove(merge_csv)

							have_first_csv=False	
							#debug.command_line_log('head -n 1 '+self.RunManager.csv_filename+' >> '+merge_csv,self.log)
							for spec in csvfiles:
								if os.path.isfile(spec):
									if not have_first_csv:
										debug.command_line_log('head -n 1 '+spec+' >> '+merge_csv,self.log)
										have_first_csv=True
									debug.command_line_log('tail -n +2 '+spec+' >> '+merge_csv,self.log)
									os.remove(spec)
							## Change this so that the plot functions find it:
							self.RunManager.csv_filename=merge_csv

		### remove temporary directories
		if not self.runsettings.debug:
			shutil.rmtree(self.temp_dir)
		

	def launch(self):
		self.start_scan()
		self.run()
		self.finish_scan()
		if self.doplots and self.inputs['MPI']['Rank'] == 0:
			os.chdir(self.plotdir)
			from bsmart import BSMplots
			BSMplots.make_plots(self)
			
	'''
	# This prototype is not needed
	def generate_parameter_points(self):
		raise NotImplementedError("You need to implement the function 'generate_parameter_points' in the Scan class" )
	'''
	

	def run(self):
		raise NotImplementedError("You need to implement the function 'run' in the Scan class" ) 

	'''
	def postprocess_no_slha(self,temp_dir,log, lock=None):
		""" Default postprocessing method for use without SLHA output, you should overload this for more interesting/sophisticated scans """
		if os.path.exists(os.path.join(temp_dir,self.runsettings.spc_filename)):
			return True
		else:
			return False
	'''
	

	def postprocess(self,Point, observables, data_point,temp_dir,log, lock=None):
		""" Default postprocessing method for use with SLHA output, you should overload this for more interesting/sophisticated scans """
		if os.path.exists(os.path.join(temp_dir,self.runsettings.spc_filename)):
			return True
		else:
			return False
		
	def postprocess_no_slha(self,Point, observables, data_point,temp_dir,log, lock=None):
		""" 
			
			Default postprocessing method for use without SLHA output, you should overload this for more interesting/sophisticated scans

			By default I will take this to be the 
		 
	    """
		return self.postprocess(Point, observables, data_point,temp_dir,log, lock)
	##### HELPER FUNCTIONS 


	def lhsetcommonspheno(self):
		""" Sets useful sphenooptions """
		self.inputspc.blocks['SPHENOINPUT']['11']=0
		self.inputspc.blocks['SPHENOINPUT']['16']=0
		self.inputspc.blocks['SPHENOINPUT']['57']=0
		self.inputspc.blocks['SPHENOINPUT']['75']=0
		self.inputspc.blocks['SPHENOINPUT']['76']=0
		self.inputspc.blocks['SPHENOINPUT']['79']=0
		self.inputspc.blocks['SPHENOINPUT']['440']=0
		self.inputspc.blocks['SPHENOINPUT']['441']=0
		self.inputspc.blocks['SPHENOINPUT']['530']=0
		self.inputspc.blocks['SPHENOINPUT']['550']=0
	 
	def lhsetdiagonal(self,block,val):
		""" sets the given block to a diagonal value val """
		if block in self.inputspc.blocks:
			for key in self.inputspc.blocks[block].keys():
				parsed=list(map(int,key.split(',')))
				x0=parsed[0]
				for x in parsed[1:]:
					if x != x0:
						self.inputspc.blocks[block][key]=0.0
						break
				if x == x0:
					self.inputspc.blocks[block][key]=val

		
	

	def lhsetSUSYScalars(self,m0):
		SUSYblocks=['MSD2IN','MSE2IN','MSl2IN','MSQ2IN','MSU2IN']
		if(isinstance(m0,string_types)):
			ms2=m0+'*'+m0
		else:
			ms2=m0**2
		for block in SUSYblocks:
			if block in self.inputspc.blocks:
				for x in range(1,4):
					keystr=str(x)+','+str(x)
					self.inputspc.blocks[block][keystr]=ms2

	def lhcommonsusyscale(self,mSUSY):
		""" sets a common SUSY scale for use with low scale input """
		self.lhsetSUSYScalars(mSUSY)
		if(isinstance(mSUSY,string_types)):
			ms2=mSUSY+'*'+mSUSY
		else:
			ms2=mSUSY**2
		if 'EXTPAR' in self.inputspc.blocks:
			self.inputspc.blocks['EXTPAR']['1']=mSUSY
			self.inputspc.blocks['EXTPAR']['2']=mSUSY
			self.inputspc.blocks['EXTPAR']['3']=mSUSY
			self.inputspc.blocks['EXTPAR']['24']=mSUSY   ## mu term
			self.inputspc.blocks['EXTPAR']['25']=ms2 ## MA2
		
	def lhsetSUSY(self,m0,m12,T0):
		""" Sets common SUSY soft scalar masses, gaugino masses and stop trilinear """

		self.lhsetSUSYScalars(m0)

		if(isinstance(m0,string_types)):
			ms2=m0+'*'+m0
		else:
			ms2=m0**2
		
		if 'EXTPAR' in self.inputspc.blocks:
			self.inputspc.blocks['EXTPAR']['1']=m12
			self.inputspc.blocks['EXTPAR']['2']=m12
			self.inputspc.blocks['EXTPAR']['3']=m12
			self.inputspc.blocks['EXTPAR']['24']=m0   ## mu term
			self.inputspc.blocks['EXTPAR']['25']=ms2 ## MA2

			
		if 'TUIN' in self.inputspc.blocks:
			self.inputspc.blocks['TUIN']['3,3']=T0
			
			


