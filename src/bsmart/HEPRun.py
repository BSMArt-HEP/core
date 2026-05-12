"""
--------------------------------------------------------------------------------------------------
HEPRun class for handling running and IO with HEP tools.

Written by Mark Goodsell

Part of BSMArt, but it is self-contained (nb requires zslha and the classes for the different tools)
and you can import it into your own programs

In that case the useful classes are:

HEPRun (which manages everything)
CoreRunner (which manages running on one core)
DataPoint

You can use HEPRun to run batches of points (using multiprocessing!)
or CoreRunner.run_point(point) to run individual points

This is useful for MPI running mode, where the multiprocessing is handled outside the program.




Allowed settings in the input dictionary ("Setup" block unless otherwise specified):
    - RunName: Name of the run.
    - Cores: Number of cores to use.
    - TempDir: Temporary working directory.
    - debug: Enable debug mode (boolean).
    - SLHA_Output: Boolean, control SLHA output generation.
    - Skip SLHA Files: Boolean, for pure python scans (skips I/O).
    - Output File: Name for collected outputs.
    - Spectrum File: Name for intermediate spectrum file.
    - StoreInputs: Boolean, whether to store input Les Houches files.
    - StoreAllPoints: Boolean, whether to store all points (default True).
    - Store Output Files: Boolean, alias for StoreAllPoints.
    - Short: Boolean, enable tabbed output.
    - csv: Boolean, enable CSV output.
    - csv name: Name for the CSV file.
    - Precision: Integer, precision for observables.
    - Store Invalid: Boolean, store invalid points.
    - Store In Memory: Boolean, store DataPoint objects in memory.
    - StoreSeparateFiles: Boolean, store separate spectrum files.
    - StoreEverything: Boolean, retain all inputs and outputs for each point.
    - Use TQDM: Boolean, use tqdm progress bar.

Other blocks:
    - Codes: Dictionary of codes to run.
    - Observables: Dictionary of observables.
    - Variables: Dictionary of variables.
    - Derived Variables: Dictionary of derived variables.
    - Fitting: Dictionary for fitting inputs.
    - MPI: Dictionary containing MPI settings (Size, Rank).

--------------------------------------------------------------------------------------------------
"""
# For use with traceback.print_exc() in debugging:
#import traceback


# standard packages
import os
import sys
from bsmart import debug
import logging

import shutil
import multiprocessing as mp

try:
    if mp.get_start_method() != 'fork':
        mp.set_start_method('fork', force=True)
except (ValueError, RuntimeError): 
    pass


from collections import OrderedDict
from bsmart import zslha
import importlib
import queue
import time

import math

import copy

from typing import Optional, Union, Callable

try:
    from tqdm import tqdm

    have_tqdm = True
except:
    have_tqdm = False


##############################################################################################################################
#
# Define here the RunSettings class
#
##############################################################################################################################

class RunSettings:
    """
    Container class to store settings to be passed to the running module
    Ordinarily these are filled in in core.py
    """

    def __init__(self, inputs: Union[str, dict, OrderedDict] = None):

        self.name: str = ""
        self.cores: int = 1  # Number of cores
        self.work_dir: str = ""
        self.output_dir: str = ""
        self.inputs_dir: str = ""
        self.store_outputs: bool = (
            True  # Whether to store output spc files, either concatenated or as separate files
        )
        self.store_inputs: bool = (
            False  # Whether to store input Les Houches files, false by default
        )
        self.debug: bool = False  # Whether to enable additional debugging
        self.MAINpath: str = ""
        #self.packagepath: str = ""  
        self.codes: dict = {}
        self.observables = OrderedDict()
        self.input_filename: str = ""
        self.output_filename: str = "SLHA_output"
        self.spc_filename: str = "SLHA_output"
        self.output_file_suffix: str = ".slha"
        self.input_file_suffix: str = ".in"
        self.no_slha_output: bool = (
            True  # Tells the scan not to expect an slha file output. This might be the case for a pythonic tool.
        )
        self.no_slha_input: bool = (
            False  # Tells the scan not to write an input file. Intended either for fake scans, or ones where all codes are handled by python
        )
        self.csv_output: bool = False
        self.csv_name: str = ""
        self.tabbed_output: bool = False
        self.tabbed_name: str = ""
        self.store_invalid_points: bool = False  # store invalid points
        self.store_separate_files = False
        self.store_points_in_memory = False  ## store DataPoint objects of the points in memory. NB not including the whole spectrum file ...
        self.deep_store_points_in_memory = False  ## Even store the spectrum files in memory. Ignored if store_points is false
        self.invalid_return_value = []  # This should be specified by the scan!
        self.store_everything = False  # Retain all inputs and outputs for each point
        self.variables = None
        self.varnames = None
        self.marginal_variables = None
        self.fitting_inputs = None
        self.use_MPI = False
        self.derived_variables = None
        self.precision = None  # numerical precision to store observables
        self.use_tqdm = have_tqdm

        if inputs is not None:
            if isinstance(inputs, str):  # have just provided a name and nothing else
                self.name = inputs
            else:
                self.name = inputs["Setup"]["RunName"]
                self.cores = inputs["Setup"]["Cores"]
                self.work_dir = inputs["Setup"]["TempDir"]
                self.debug = eval(inputs["Setup"]["debug"])
                #self.packagepath = inputs["Setup"]["PackagePath"]
                self.no_slha_output = (
                    False  # above it is set to True, but in usual running we want false
                )

                if inputs["MPI"]["Size"] > 1:
                    self.use_MPI = True

                if "SLHA_Output" in inputs["Setup"]:
                    self.no_slha_output = not eval(str(inputs["Setup"]["SLHA_Output"]))

                if (
                    "Skip SLHA Files" in inputs["Setup"]
                ):  # For scans with purely python tools, where we don't want to keep i/o and just use datapoints
                    self.no_slha_output = eval(str(inputs["Setup"]["Skip SLHA Files"]))
                    self.no_slha_input = self.no_slha_output

                if "Codes" not in inputs:
                    raise NameError("No codes found in input file")
                else:
                    self.codes = inputs["Codes"]

                    ### add observables defined within the code to the global observables list if present
                    for code in inputs["Codes"]:
                        if (
                            eval(inputs["Codes"][code]["Run"])
                            and "Observables" in inputs["Codes"][code]
                        ):
                            for obs in inputs["Codes"][code]["Observables"]:
                                inputs["Observables"][obs] = inputs["Codes"][code][
                                    "Observables"
                                ][obs]

                
                firstcode = next(iter(inputs["Codes"]))
                if (
                    "Output File" in inputs["Setup"]
                ):  # This sets the name for collected outputs: spectrum files, csv or tabbed outputs
                    self.output_filename = inputs["Setup"]["Output File"]

                if "Spectrum File" in inputs["Setup"] and inputs["Setup"]["Spectrum File"] != "":  # This sets the name for the intermediate spectrum file
                    self.spc_filename = inputs["Setup"]["Spectrum File"]
                elif "OutputFile" in inputs["Codes"][firstcode] and inputs["Codes"][firstcode]["OutputFile"] != "":  
                    # take default from the output of the first code, if it has one. This is fine if it's SPheno, for example
                    self.spc_filename = inputs["Codes"][firstcode]["OutputFile"]

                ## the input file should be the input for the first code that we run, presumably SPheno. Note we imported the json as an OrderedDict
                
                if not self.no_slha_input:
                    self.input_filename = inputs["Codes"][firstcode]["InputFile"]
                else:
                    self.input_filename = self.spc_filename
                # If the

                if "StoreInputs" in inputs["Setup"]:
                    self.store_inputs = eval(str(inputs["Setup"]["StoreInputs"]))
                else:
                    self.store_inputs = False

                self.store_outputs = True
                if "StoreAllPoints" in inputs["Setup"]:

                    self.store_outputs = eval(inputs["Setup"]["StoreAllPoints"])

                if "Store Output Files" in inputs["Setup"]:
                    self.store_outputs = eval(inputs["Setup"]["Store Output Files"])

                if "Short" in inputs["Setup"]:
                    self.tabbed_output = eval(str(inputs["Setup"]["Short"]))

                if "csv" in inputs["Setup"]:
                    self.csv_output = eval(str(inputs["Setup"]["csv"]))
                if "csv name" in inputs["Setup"]:
                    self.csv_name = str(inputs["Setup"]["csv name"])

                if "Precision" in inputs["Setup"]:
                    self.precision = int(inputs["Setup"]["Precision"])

                if "Store Invalid" in inputs["Setup"]:
                    self.store_invalid_points = eval(inputs["Setup"]["Store Invalid"])

                if "Store In Memory" in inputs["Setup"]:
                    self.store_points_in_memory = eval(
                        inputs["Setup"]["Store In Memory"]
                    )

                if "StoreSeparateFiles" in inputs["Setup"]:
                    self.store_separate_files = eval(
                        str(inputs["Setup"]["StoreSeparateFiles"])
                    )

                # self.store_everything=False  # this is already set above by default
                if "StoreEverything" in inputs["Setup"]:
                    self.store_everything = eval(
                        str(inputs["Setup"]["StoreEverything"])
                    )
                    if self.store_everything:
                        self.store_outputs = True
                        if inputs["MPI"]["Size"] > 1:
                            self.all_dir = os.path.join(
                                inputs["Setup"]["cwd"],
                                inputs["Setup"]["RunName"],
                                "All_Outputs",
                                "MPI_" + str(inputs["MPI"]["Rank"]),
                            )
                        else:
                            self.all_dir = os.path.join(
                                inputs["Setup"]["cwd"],
                                inputs["Setup"]["RunName"],
                                "All_Outputs",
                            )
                if inputs["MPI"]["Size"] > 1:
                    self.inputs_dir = os.path.join(
                        inputs["Setup"]["cwd"],
                        inputs["Setup"]["RunName"],
                        "Input_Files",
                        "MPI_" + str(inputs["MPI"]["Rank"]),
                    )
                    self.output_dir = os.path.join(
                        inputs["Setup"]["cwd"],
                        inputs["Setup"]["RunName"],
                        "Spectrum_Files",
                        "MPI_" + str(inputs["MPI"]["Rank"]),
                    )
                else:
                    self.inputs_dir = os.path.join(
                        inputs["Setup"]["cwd"],
                        inputs["Setup"]["RunName"],
                        "Input_Files",
                    )
                    self.output_dir = os.path.join(
                        inputs["Setup"]["cwd"],
                        inputs["Setup"]["RunName"],
                        "Spectrum_Files",
                    )
                # self.output_filename='SLHA_output'  # already set above

                # else:  # already set above
                #    self.spc_filename='SLHA_output'

                if "Variables" in inputs:
                    self.variables = inputs["Variables"]
                    if hasattr(inputs["Variables"], "keys"):
                        self.varnames = list(inputs["Variables"].keys())
                    else:
                        self.varnames = list(inputs["Variables"])

                if "Derived Variables" in inputs:
                    """
                    Allow for functions of the variables, which should be provided as a dictionary like
                    {"dvar1": "math.exp(mh)", ...}
                    These will be treated like variables when writing the input files etc, but won't show up in the csv file
                    """
                    self.derived_variables = inputs["Derived Variables"]

                if "Observables" in inputs:
                    self.observables = inputs["Observables"]

                if "Fitting" in inputs:
                    self.fitting_inputs = inputs["Fitting"]

                if "Use TQDM" in inputs["Setup"]:
                    self.use_tqdm = have_tqdm and eval(
                        inputs["Setup"]["Use TQDM"]
                    )  # use it if we have it!



##############################################################################################################################
#
# Define here the DataPoint class
#
##############################################################################################################################

class DataPoint:
    """Class to hold info about each point.
    We need at least an ID.
     A list of inputs is the first step to filling.
     observables are a list of floats
     outputs are a list of *results* from the computation, intended to be returned to the scan. E.g. a likelihood
     valid is whether the point passes validity checks etc


    """

    def __init__(self, ID: int, inputs: list = None):
        self.ID = ID
        if inputs is None:
            self.inputs = []
        else:
            self.inputs = inputs

        # self.varnames: list=None # list of names of the variables
        self.outputs: list = []
        self.observables: list = []
        self.valid = False
        self.var_dict = None
        self.varnames = None
        self.marginal_variables = None
        self.obs_dict = OrderedDict()  # future use
        self.spc = None

    def clear_obs(self):
        self.outputs = []
        self.observables = []
        self.valid = False
        # self.var_dict=None
        self.obs_dict = OrderedDict()  # future use
        self.spc = None

    def make_var_dict(self, var_names: list):
        self.var_dict = {
            var_name: var_val for var_name, var_val in zip(var_names, self.inputs)
        }
        self.var_dict["VARIABLE"] = self.inputs
        self.var_dict["math"] = math

    def get_var_dict(self, var_names: list = None):
        if self.var_dict is not None:
            return self.var_dict

        if var_names is not None:
            var_dict = {
                var_name: var_val for var_name, var_val in zip(var_names, self.inputs)
            }
            var_dict["VARIABLE"] = (
                self.inputs
            )  ## Added for backward compatibility, maybe remove?
        else:
            if (
                self.var_dict is not None
            ):  # don't want to overwrite it if we already evaluated it!
                return self.var_dict
            var_dict = {
                "VARIABLE": self.inputs
            }  ## Added for backward compatibility, maybe remove?
        var_dict["math"] = (
            math  ## to allow maths functions such as sqrt or cos/sin/tan!
        )
        return var_dict

    # def calculate_derived_vars(self,derived_vars: dict):
    #    for var in

    def get_all_dict(self) -> dict:
        if self.var_dict is not None:
            res = self.var_dict
        else:
            res = self.get_var_dict()
        if self.obs_dict is not None:
            res.update(self.obs_dict)
        return res

    def to_dict(self, all_data: bool = False):
        """
        Creates a dictionary representation.
        If all_data is True, we export the spc object too
        """

        # return dict(self.__dict__)
        res = {"ID": self.ID, "inputs": self.inputs}

        # for var,val in self.var_dict.items():
        #    if var != 'math' and var != '__builtins__':
        #        print(f"{var}")
        res["variables"] = {
            var_name: var_val
            for var_name, var_val in self.var_dict.items()
            if var_name != "math" and var_name != "__builtins__"
        }
        if self.obs_dict is not None:
            res["variables"].update(self.obs_dict)
        res["observables"] = list(self.observables)
        res["results"] = self.outputs
        res["valid"] = self.valid
        # res['observables'] = copy.deepcopy(self.observables)
        # res['results'] = copy.deepcopy(self.outputs)
        # res['valid'] = self.valid
        if all_data and self.spc is not None:
            res["spc"] = copy.deepcopy(self.spc.to_dict())
        return res


def UpdateProgress(number, length, message="", finish=False):
    """
    Ridiculously basic progress bar
    Implemented here so that we don't need to install extra packages

    If tqdm is installed, you can use it in this function. See the example below:
    from tqdm import tqdm

    with tqdm(total=length) as pbar:
        for i in range(number):
            # Do stuff...
            pbar.update()

    """
    ## below is not compatible with python < 3.6
    # string=f'['
    if length > 0:
        string = "["
        done = number * 50 // length
        todo = 50 - done
        percentage = number * 100 // length
        if done != 0:
            done = done - 1
        string = (
            string
            + ("=" * done)
            + ">"
            + (" " * todo)
            + "] "
            + str(percentage)
            + "% "
            + message
        )
        #    if number == length:
        if finish:
            print(string)
        else:
            print(string, end="\r")


def _normalize_valid_bound(val):
    """
    Normalize bounds used in VALID intervals to floats or None.

    Accepts:
      - None                     -> None
      - float/int (incl. +/-inf) -> float(val)
      - strings like 'inf', '+inf', 'infty', 'infinity', '-inf' -> +/-inf
      - other strings convertible to float -> float(value)
    """
    if val is None:
        return None

    if isinstance(val, (int, float)):
        return float(val)

    if isinstance(val, str):
        s = val.strip().lower()
        if s in ("inf", "+inf", "infty", "infinity", "+infty", "+infinity"):
            return float("inf")
        if s in ("-inf", "-infty", "-infinity"):
            return float("-inf")
        # fall back to numeric parse
        return float(s)

    return float(val)


def create_validator(obs, range):
    # Normalise bounds so that VALID can safely use numeric or string +/-inf
    norm_range = [_normalize_valid_bound(v) for v in range]
    varmax = max(norm_range)
    varmin = min(norm_range)

    def tfunc(obs_dict):
        val = float(eval(obs, obs_dict))
        if val < varmin:
            return False
        if val > varmax:
            return False
        return True

    return tfunc


def create_condition(cond_string):
    def tfunc(obs_dict):
        # print(f"Checking condition {cond_string}:")
        # print(f"{bool(eval(cond_string,obs_dict))}")
        return bool(eval(cond_string, obs_dict))

    return tfunc


"""
def smart_format(value,precision):
    if isinstance(value, (float)):
        return "{0:.{1}g}".format(value, precision)
    return str(value)
"""

"""
def format_variable(value):
    if isinstance(value, (float)):
        return str(round(value,8))
    return str(value)
"""


def format_variable(value):
    if isinstance(value, (int)):
        return str(value)
    elif isinstance(value, (float)):
        return "{0:.{1}g}".format(value, 9)
    return str(value)


def create_smart_format(precision):
    def tfunc(value):
        if isinstance(value, (int)):
            return str(value)
        elif isinstance(value, (float)):
            return "{0:.{1}g}".format(value, precision)
        return str(value)

    return tfunc




##############################################################################################################################
#
# Define here the HepTool class
#
##############################################################################################################################

class HepTool:
    """
    Template class to run the different HEP tools
    Codes are supposed to place their output in spc_file which will get read at the end

    At init we may also provide global settings so that the code can take what they need (e.g. number of cores)
    These should not be stored by default

    But they may want access at intermediate stages to things like the input variables
    or details of the spectrum. For that we can use the DataPoint class

    """
    
    
    def __init__(self, name, settings, global_settings=None):
        self.name = name
        self.citations = None  # allow tools to specify references
        self.settings = settings
        self.check_validity = False
        self.validity_checks = []
        if "Observables" in self.settings:
            for obs in self.settings["Observables"]:
                if "VALID" in self.settings["Observables"][obs]:
                    self.validity_checks.append(
                        create_validator(
                            obs, self.settings["Observables"][obs]["VALID"]
                        )
                    )
                    self.check_validity = True

        if "Conditions" in self.settings:
            for tcond in self.settings["Conditions"]:
                self.validity_checks.append(create_condition(tcond))
                self.check_validity = True

        # whether to look for files
        self.no_slha_output = False
        if global_settings is not None and global_settings.no_slha_output:
            self.no_slha_output = True

    def validity_check(self, obs_dict):
        return all([x(obs_dict) for x in self.validity_checks])

    def launch(self, spc_file, temp_dir, log, data_point):
        log.debug("Running %s " % self.name)
        os.chdir(temp_dir)

        try:
            self.run(spc_file, temp_dir, log, data_point)
        except:
            log.debug("Failed to run " + self.name)
            raise

        if not self.no_slha_output:
            # Should get the tool to check this if needed. We no longer force the generation of an output file
            if not os.path.exists(os.path.join(temp_dir, spc_file)):
                log.warning("No output produced for code " + self.name)
                raise

        """ Should we put the onus on the tool to do this instead? """
        if (
            data_point.spc is not None and "Observables" in self.settings
        ):  ## collect observables into the datapoint
            # my_observables={ obs: data_point.spc.Value(val['SLHA'][0], val['SLHA'][1]) for obs,val in zip(self.settings['Observables'].keys(),self.settings['Observables'].values()}
            # data_point.obs_dict.update(my_observables)
            for obs, val in zip(
                self.settings["Observables"].keys(),
                self.settings["Observables"].values(),
            ):
                if (
                    "SLHA" in val
                ):  # New: allow for observables that are not SLHA values but may be functions.
                    data_point.obs_dict[obs] = float(
                        data_point.spc.safeValue(val["SLHA"][0], val["SLHA"][1])
                    )
                elif "FUNCTION" in val:
                    data_point.obs_dict[obs] = eval(
                        val["FUNCTION"], data_point.get_all_dict()
                    )  # allow for derived variables

        if self.check_validity:
            if not self.validity_check(data_point.get_all_dict()):
                raise NameError("Tool " + self.name + " failed validity check")

    def run(self, spc_file, log):
        ### This routine should clean any output files first, then run the code.
        raise NotImplementedError(
            "You need to implement the function 'run' in the HepTool class for "
            + self.name
        )



##############################################################################################################################
#
# Define here the write_inputfile function
#
##############################################################################################################################


def _create_default_write_inputfile(
    input_filename : str,
    log: logging.Logger,
    block_info: Union[dict, OrderedDict] = None,
):
    """
    Routine to create a function that will write an slha input file from a template and/or list of blocks. 
    """
    log.info('Reading bare spc: '+str(input_filename))

    
    try:
        inputspc= zslha.read(input_filename,verbose=False)
    except Exception as e:
        log.error('Problem reading slha file '+input_filename+', ' + str(e))
        inputspc = zslha.SLHA()
        

    # in case there is no template and we want to use json input
    if 'Blocks' in block_info:
        try:
            for current_block in block_info['Blocks']:
                
                if current_block not in inputspc.blocks:
                    inputspc.blocks[current_block]={}
                        
                for key in block_info['Blocks'][current_block].keys():
                    inputspc.blocks[current_block][key] = block_info['Blocks'][current_block][key]
        except Exception as e:
            log.debug('Exception reading json block input: '+str(e))

    # now create the function to write the input file
    def write_lh_file(
        point: list, 
        dir: str, 
        name: str,
        values_dict: dict=None,
        Overloads: list = None
    ):
        """ write Les Houches input file for given parameter point """

        lhname=os.path.join(dir,name)
        if Overloads is not None:  
            """
                Allow overloading of block values, e.g. during marginalisation

                Format: a list of form
                [["SPHENOINPUTS",[1],1.0],...]
            """
            # obtain a new copy of the input spc
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
            zslha.write_lh(runspc, point, lhname,self.varnames)
        else:
            zslha.write_lh_dict(runspc, values_dict, lhname)

    return write_lh_file





##############################################################################################################################
#
# Define here the HEPRun class
#
##############################################################################################################################


class HEPRun:
    """ 
    This is the main class to handle running. 
    
    """

    def __init__(self, 
    runsettings: RunSettings = None,
    write_inputfile: Callable = None,
    postprocess: Callable = None,
    log: logging.Logger = None,
    inputs: Union[dict,OrderedDict] = {},
    ):
        if runsettings is None: 
            # use the inputs instead!
            self.settings = RunSettings(inputs)
        else:
            self.settings = runsettings

        if write_inputfile is None:
            self.write_inputfile = _create_default_write_inputfile(
                self.settings.input_filename,
                self.log,
                inputs,
            )
        else:
            self.write_inputfile = write_inputfile

        self.postprocess = postprocess
        self.work_dir = runsettings.work_dir
        self.use_tqdm = runsettings.use_tqdm
        if log == None:
            self.log = debug.new_logger(
                self.settings.debug,
                self.settings.name,
                os.path.join(self.settings.work_dir, self.settings.name),
            )
            self.log.debug("Initialise HEPRun")
        else:
            log.debug("Initialise HEPRun")
            self.log = log

        if self.settings.store_points_in_memory:
            self.all_points = []
            self.valid_points = []

        self.all_batch_points = []
        self.valid_batch_points = []

        self.nextpointID = 1

        if self.settings.cores > 1:
            self.input_points = mp.Queue()
            self.output_fileheaders = mp.Queue()
            self.invalid_fileheaders = mp.Queue()
            self.progress_messages = mp.Queue()
        else:
            self.input_points = queue.Queue()
            self.output_fileheaders = queue.Queue()
            self.invalid_fileheaders = queue.Queue()
            self.progress_messages = queue.Queue()

        if self.settings.fitting_inputs is not None:

            self.do_marginalisation = True
            if "Method" in self.settings.fitting_inputs:
                self.marginalisation_method = self.settings.fitting_inputs["Method"]
            else:
                self.marginalisation_method = None
            self.set_up_marginalisation()

            if "Options" in self.settings.fitting_inputs:
                ## should be dict, but need to evaluate if there is a "True"
                self.marginalisation_options = {
                    x: eval(str(self.settings.fitting_inputs["Options"][x]))
                    for x in self.settings.fitting_inputs["Options"]
                }

            else:
                self.marginalisation_options = {"ftol": 1e-1, "disp": False, "eps": 0.1}
                # print('Options:')
                # print(dict(eval(self.marginalisation_options)))
                # raise

        else:
            self.do_marginalisation = False

        self.tools = []
        self.set_up_codes()

        ## Set up a run handler for each core
        self.CoreRuns = []

        for x in range(self.settings.cores):
            self.CoreRuns.append(CoreRunner(self, "Core_" + str(x)))

        self.all_tools_run_state()

        ## set up an output file and how we format output if we want to store csv or tab values
        if self.settings.csv_output:

            if self.settings.csv_name == "":
                # self.csv_filename=os.path.join(self.settings.output_dir,self.settings.output_filename)
                # if self.settings.store_outputs:
                log.info(
                    "CSV output without a csv file name -> appending .csv to output filename"
                )
                # self.csv_filename=self.csv_filename+'.csv'
                csv_filename = self.settings.output_filename + ".csv"
            else:
                csv_filename = self.settings.csv_name
            self.csv_filename = os.path.join(self.settings.output_dir, csv_filename)

            # else:
            #    self.csv_filename=os.path.join(self.settings.output_dir,self.settings.csv_name)

            if os.path.isfile(self.csv_filename):
                os.remove(self.csv_filename)
            self.output_fileheaders.put([self.csv_filename, ","])

            if self.settings.store_invalid_points:  # keep a csv file for invalid points
                self.invalid_csv_filename = os.path.join(
                    self.settings.output_dir, "invalid_" + csv_filename
                )
                self.invalid_fileheaders.put([self.invalid_csv_filename, ","])

        if self.settings.tabbed_output:
            if self.settings.tabbed_name == "":
                tabbed_name = self.settings.output_filename + ".dat"
                self.tabbed_filename = os.path.join(
                    self.settings.output_dir, self.settings.output_filename
                )
                if self.settings.store_outputs:
                    log.info(
                        "Tabbed output without a tab file name -> appending .dat to output filename"
                    )
                    self.tabbed_filename = self.tabbed_filename + ".dat"
            else:
                tabbed_name = self.settings.tabbed_name
                self.tabbed_filename = os.path.join(
                    self.settings.output_dir, self.settings.tabbed_name
                )
            if os.path.isfile(self.tabbed_filename):
                os.path.remove(self.tabbed_filename)
            self.output_fileheaders.put([self.tabbed_filename, "\t"])
            if self.settings.store_invalid_points:  # keep a csv file for invalid points
                self.invalid_tabbed_filename = os.path.join(
                    self.settings.output_dir, "invalid_" + tabbed_name
                )
                self.invalid_fileheaders.put([self.invalid_tabbed_filename, "\t"])

        # Define functions for formatting variables and observables
        # In zslha the variables are printed out only to 8 significant figures anyway.
        """
        def format_variable(value):
            if isinstance(value, (float)):
                return "{0:.{1}g}".format(value, 8)
            return str(value)
        if self.settings.precision is not None:
            def format_observable(self,value):
                print(value)
                if isinstance(value, (float)):
                    return "{0:.{1}g}".format(value, self.settings.precision)
                return str(value)    
        else:
            def format_observable(self,value):
                print(value)
                return str(value)    
        """

        self.format_variable = format_variable  # defined above, rounds to 8 decimal places, since that is what is done in zslha
        if self.settings.precision is not None:
            self.format_observable = create_smart_format(self.settings.precision)
        else:
            self.format_observable = lambda x: str(x)
        if self.settings.store_inputs:
            self.input_filename_prefix = os.path.join(
                self.settings.inputs_dir, self.settings.name + "_"
            )

        if self.settings.store_outputs:
            if self.settings.store_separate_files:
                # self.output_filename_prefix=os.path.join(self.settings.output_dir,self.settings.name+'_')
                self.output_filename_prefix = os.path.join(
                    self.settings.output_dir, self.settings.output_filename + "_"
                )
            else:
                self.output_filename = os.path.join(
                    self.settings.output_dir, self.settings.output_filename
                )

            if self.settings.store_everything:
                # self.output_directory_prefix=os.path.join(self.settings.all_dir,self.settings.name+'_')
                # self.output_filename_prefix=os.path.join(self.settings.all_dir,self.settings.output_filename+'_')
                self.output_directory_prefix = os.path.join(
                    self.settings.all_dir, self.settings.output_filename + "_"
                )

    def set_up_marginalisation(self):
        import bsmart.marginalisation as marginalisation

        if "Variables" not in self.settings.fitting_inputs:
            raise NameError("Fitting must include Variables!")

        self.marginal_variables = list(self.settings.fitting_inputs["Variables"].keys())
        # print('Marginal variables: ' +str(self.marginal_variables))
        self.marginal_downscalers, self.marginal_upscalers, self.marginal_ranges = (
            marginalisation.create_scalers_marginal(
                self.settings.fitting_inputs, self.log
            )
        )
        self.number_tools_run_marginalisation = 0
        # print('here')
        if "Observables" not in self.settings.fitting_inputs:
            raise NameError("Fitting must include Observables!")

        self.marginal_observables = list(
            self.settings.fitting_inputs["Observables"].keys()
        )
        if (
            (
                self.marginalisation_method == None
                or self.marginalisation_method == "Brentq"
                or self.marginalisation_method != "Minimise"
            )
            and len(self.marginal_variables) == 1
            and len(self.marginal_observables) == 1
        ):
            self.oned_marginal = True
            self.marginal_chisquared = marginalisation.create_diff_marginal(
                self.settings.fitting_inputs["Observables"]
            )
        else:
            if self.marginalisation_method == "Minimise":
                self.marginalisation_method = None
            self.oned_marginal = False
            self.marginal_chisquared = marginalisation.create_chisquare_marginal(
                self.settings.fitting_inputs["Observables"]
            )

    def set_up_codes(self):
        """Set up the HEP Tools needed in the scan"""
        # classpath = os.path.join(self.settings.packagepath, "tools")
        found_observables = []
        for toolname in self.settings.codes:
            self.log.debug("Setting up " + toolname)
            if eval(self.settings.codes[toolname]["Run"]):
                """
                Allow definition of observables inside the tool
                """
                if self.do_marginalisation and len(found_observables) < len(
                    self.marginal_observables
                ):
                    self.number_tools_run_marginalisation += 1
                if "Observables" in self.settings.codes[toolname]:
                    for obs in self.settings.codes[toolname]["Observables"]:
                        # print(obs)
                        self.settings.observables[obs] = self.settings.codes[toolname][
                            "Observables"
                        ][obs]

                        if self.do_marginalisation and obs in self.marginal_observables:
                            found_observables.append(obs)

                # Try to import from bsmart.tools
                tool_module = None
                try:
                    tool_module = importlib.import_module("bsmart.tools." + toolname)
                except ImportError:
                    pass

                if tool_module:
                    try:
                        self.tools.append(
                            tool_module.NewTool(
                                toolname, self.settings.codes[toolname], self.settings
                            )
                        )
                    except Exception as e:
                        self.log.critical(
                            "Failed to initialise " + toolname + " from bsmart.tools, error " + str(e)
                        )
                        raise SystemExit
                else:
                    self.log.warning(
                        "Required tool "
                        + toolname
                        + " not found in bsmart.tools, trying user directory ./Tools"
                    )
                    ## try a user directory
                    classfile = os.path.join(os.getcwd(), "Tools", toolname + ".py")
                    if not os.path.isfile(classfile):
                        self.log.critical(
                            "Required tool not present "
                            + toolname
                            + ", should be in file "
                            + classfile
                        )
                        raise SystemExit
                    try:
                        if os.path.join(os.getcwd(), "Tools") not in sys.path:
                            sys.path.append(os.path.join(os.getcwd(), "Tools"))
                        new_class = importlib.import_module(toolname)
                        self.tools.append(
                            new_class.NewTool(
                                toolname, self.settings.codes[toolname], self.settings
                            )
                        )
                    except Exception as e:
                        self.log.critical(
                            "Failed to initialise " + toolname + " ,error " + str(e)
                        )
                        raise SystemExit
        ## If we didn't find all the observables in the tools
        if self.do_marginalisation and len(found_observables) < len(
            self.marginal_observables
        ):
            self.extra_observable_step = True
        else:
            self.extra_observable_step = False

    def all_tools_run_state(self, state=True):
        for tool in self.tools:
            for core in self.CoreRuns:
                core.run_tools[tool.name] = state

    def set_tool_run_state(self, tool, state=True):
        for core in self.CoreRuns:
            core.run_tools[tool] = state

    def newpoint(self, inputs):
        nextID = self.nextpointID
        self.nextpointID = self.nextpointID + 1
        return DataPoint(nextID, inputs)

    def new_point_to_stack(self, inputs):
        nextID = self.nextpointID
        self.nextpointID = self.nextpointID + 1
        # self.input_points.put(DataPoint(nextID,inputs))
        self.input_points.put(
            {"ID": nextID, "inputs": inputs}
        )  # since we're not storing anything here

    def run_batch(self, sample):
        sample_size = len(sample)
        sample_start = self.nextpointID
        for point in sample:
            self.new_point_to_stack(point)
        for _ in range(len(self.CoreRuns)):
            self.input_points.put(None)  ## Add a 'sentinel' to indicate termination
        self.log.info(
            "Starting batch run of "
            + str(sample_size)
            + " points over "
            + str(self.settings.cores)
            + " cores"
        )
        start_time = time.time()

        if not self.use_tqdm:
            UpdateProgress(0, sample_size, (" 0 of " + str(sample_size) + " done"))

        self.all_batch_points = []
        self.valid_batch_points = []
        if self.settings.cores > 1:
            with mp.Manager() as manager:
                # List_all = manager.list()   # better to replace this with a queue!
                # Using a manager queue is more robust, especially with 'spawn' or 'forkserver' start methods
                # which are becoming default in newer Python versions.
                result_queue = manager.Queue()
                lock = manager.Lock()

                # processes = [mp.Process(target=corerun.run_all,args=(List_all, List_valid, sample_start, sample_size,lock)) for corerun in self.CoreRuns]
                processes = [
                    mp.Process(
                        target=corerun.run_all,
                        args=(result_queue, sample_start, sample_size, lock),
                    )
                    for corerun in self.CoreRuns
                ]
                for p in processes:
                    p.start()

                if self.use_tqdm:
                    # with tqdm(total=sample_size,ncols=50) as pbar:
                    pbar = tqdm(total=sample_size, ncols=50)
                else:
                    pbar = None
                while True:
                    try:
                        item = result_queue.get(timeout=0.1)

                        # print(f"Got item {item}")
                        if item["valid"]:
                            self.valid_batch_points.append(item)
                        self.all_batch_points.append(item)
                        if self.use_tqdm:
                            pbar.update(1)
                    except (
                        queue.Empty
                    ):  # Queue.Empty, not necessarily finished! Check! Could just have timed out, weird!
                        if all(not w.is_alive() for w in processes):
                            break
                        continue
                    except KeyboardInterrupt:  # something went wrong or interrupted
                        self.log.warning("Running interrupted")
                        for p in processes:
                            p.terminate()
                        for p in processes:
                            p.join()
                        raise
                    except (
                        Exception
                    ) as e:  # Something else when wrong. Should not ignore
                        self.log.error("Running hit exception " + str(e))
                        for p in processes:
                            p.terminate()
                        for p in processes:
                            p.join()
                        raise
                    # except KeyboardInterrupt:
                    # except: # something went wrong
                    #    continue

                for p in processes:
                    p.join()

                nfinal = len(self.all_batch_points)
                nvalid = len(self.valid_batch_points)

        else:
            result_queue = queue.Queue()
            if self.use_tqdm:
                pbar = tqdm(total=sample_size, ncols=50)
            else:
                pbar = None
            self.CoreRuns[0].run_all(result_queue, sample_start, sample_size, Lock=None, pbar=pbar)

            while True:
                try:
                    item = result_queue.get_nowait()
                except queue.Empty:  # Queue.Empty, not necessarily finished! Check!
                    break

                except Exception as e:
                    raise NameError("Error collecting results, " + str(e))

                if item["valid"]:
                    self.valid_batch_points.append(item)
                self.all_batch_points.append(item)
                if self.use_tqdm:
                    pbar.update(1)
            nfinal = len(self.all_batch_points)
            nvalid = len(self.valid_batch_points)

        if self.settings.store_points_in_memory:
            self.all_points = self.all_points + self.all_batch_points
            self.valid_points = self.valid_points + self.valid_batch_points

        message = "Batch finished, %d of %d points evaluated are valid" % (
            nvalid,
            nfinal,
        )
        if not self.use_tqdm:
            UpdateProgress(sample_size, sample_size, message, True)
        self.log.info("Batch finished ...")
        self.log.info(
            "... with %d points in total, of which %d are valid" % (nfinal, nvalid)
        )
        end_time = time.time()
        self.log.info("Batch run took %f seconds " % (end_time - start_time))

        ### Now return the results, just the result values in the correct order:
        if self.settings.cores > 1:
            return_vals = [[] for x in range(len(self.all_batch_points))]
            reordered_all_batch = [None] * len(return_vals)
            for point in self.all_batch_points:
                return_vals[point["ID"] - sample_start] = point["results"]
                reordered_all_batch[point["ID"] - sample_start] = point
            self.all_batch_points = reordered_all_batch
            # no need to reorder the valid list, since we know it will have holes?

            # otherwise could do
            # self.valid_batch_points.sort(key=lambda x: x["ID"])
            #

            # return_vals[point.ID-sample_start] = point.outputs  ## use the ID to order the points
            return return_vals

        else:  ## points already correctly ordered
            return [point["results"] for point in self.all_batch_points]
            # return [point.outputs for point in self.all_batch_points] ### note we set the output to self.settings.invalid_return_value



##############################################################################################################################
#
# Define here the CoreRunner class
#
##############################################################################################################################

class CoreRunner:
    """Class to handle running on one core, takes HEPRun as parent, which becomes the self.RunHandler"""

    def __init__(self, ParentRunner, name):
        ParentRunner.log.debug("Setting up runner name " + name)
        ## set up temp directory and output directory
        self.RunHandler = ParentRunner
        # self.output_file=ParentRunner.settings.output_filename
        self.spc_file = ParentRunner.settings.spc_filename

        self.name = name
        self.temp_dir = os.path.join(self.RunHandler.work_dir, name)
        if os.path.isdir(self.temp_dir):
            shutil.rmtree(self.temp_dir)
        os.makedirs(self.temp_dir)
        self.log = debug.new_core_logger(
            self.RunHandler.settings.debug, name, os.path.join(self.temp_dir, name)
        )

        self.last_observables = None
        # self.last_spc=None   # this won't pickle

        self.last_point_dict = None  # We need this to pass all point information back and forth, in case we are storing it.

        self.run_tools = (
            {}
        )  ### defined to allow tools to not necessarily be run every time
        # for code in self.ParentRunner.tools

    def __del__(self):
        """destructor method: delete working directory to clear up?"""

        if not self.RunHandler.settings.debug:
            if os.path.isdir(
                self.temp_dir
            ):  ## needed because we may have already deleted it
                shutil.rmtree(self.temp_dir)

    # def run_all(self, all_list, valid_list, sample_start, sample_size, Lock=None):
    def run_all(self, run_queue, sample_start, sample_size, Lock=None, pbar=None):
        """
        Runs points from a list; usually invoked by run_batch in the RunHandler.
        """
        while True:  # not self.RunHandler.input_points.empty():
            try:
                # point=self.RunHandler.input_points.get()
                self.last_point_dict = self.RunHandler.input_points.get(timeout=1.0)
                # print(f"got {self.last_point_dict}")
            except:
                # self.RunHandler.progress_messages.put('')
                return 0

            if self.last_point_dict is None:
                return 0

            message = ""

            point_number_in_batch = int(self.last_point_dict["ID"]) - sample_start + 1
            try:



                if self.RunHandler.settings.store_points_in_memory:
                    # Do not need to pass the point data, it is already specified above
                    # output=self.run_point(None,None,lock=Lock,point_data=point_dict) # should be able to access this object
                    output = self.run_point(None, None, lock=Lock)
                else:
                    output = self.run_point(
                        self.last_point_dict["inputs"],
                        self.last_point_dict["ID"],
                        lock=Lock,
                    )

                # point_dict = copy.deepcopy(self.last_point_dict)
                # point_dict = self.last_point_dict
                self.last_point_dict["valid"] = True
                if "results" not in self.last_point_dict:
                    self.last_point_dict["results"] = output


                message = (
                    "Batch point "
                    + str(point_number_in_batch)
                    + "/"
                    + str(sample_size)
                    + " is valid! "
                )
                if self.RunHandler.settings.debug:
                    message = message + str(self.last_point_dict["inputs"])
                    # message=message+str(point.inputs)

            except Exception as e:
                #traceback.print_exc() 
                self.log.debug("Running hit exception " + str(e))
                # point_dict = self.last_point_dict
                self.last_point_dict["valid"] = False
                self.last_point_dict["results"] = (
                    self.RunHandler.settings.invalid_return_value
                )

                message = (
                    "Batch point "
                    + str(point_number_in_batch)
                    + "/"
                    + str(sample_size)
                    + " is not valid! "
                )
                if self.RunHandler.settings.debug:
                    message = message + str(self.last_point_dict["inputs"])


            run_queue.put(self.last_point_dict)


            # self.RunHandler.progress_messages.put(message)
            if not self.RunHandler.use_tqdm:
                UpdateProgress(point_number_in_batch, sample_size, message)
            elif pbar: # only provided in single core mode.
                pbar.update(1)
            
    def run_point(self, Point=None, ID=None, lock=None):
        """
        Run a single point: create LH file, set codes to work on it; uses a single core for this
        If ID is None, we assume we are running in single core mode

        Point is just a list of floats corresponding to the input variables.
        """

        if Point is None:  # get data indirectly via the dictionary.
            Point = self.last_point_dict[
                "inputs"
            ]  # if this is missing, get an exception -> ok. Still need this further down
            tID = self.last_point_dict["ID"]
            point_container = DataPoint(self.last_point_dict["ID"], Point)
        else:
            if not ID:
                """running in single core mode, no batch, but we may want to store the point/files for which we need a unique ID"""
                point_container = self.RunHandler.newpoint(Point)
                tID = point_container.ID
            else:
                tID = ID
                """ 
                        If we are handling the cores within the scan, we might want to still store the points (e.g. in the MCMC scan). 
                        Then we need to create a point container for the point

                        Change of plan: now we always create the point so that we can pass info between the tools
                    """
                # if self.RunHandler.settings.store_points and make_point_container:
                #   point_container=DataPoint(tID,Point)
                point_container = DataPoint(tID, Point)

        """ 
            If we have no point ID, then we are running a single point in single core mode, 
            so we should stick it on the stack if we are storing points. BUT this is no good in multicore mode: 
            the points will be lost at the end. So instead we just keep the *last* point and leave it up to the 
            scan to handle the storage in that case.
        """
        # if self.RunHandler.settings.store_points and ((not ID) or (make_point_container)):

        if self.RunHandler.settings.varnames is not None:
            point_container.make_var_dict(self.RunHandler.settings.varnames)

        if self.RunHandler.settings.derived_variables is not None:
            """
            Derived variables should be a function of the variables.
            We expect a dictionary of the form
            {"name": "expression", ... }

            If you need something that depends on the observables too, then I will treat that as an observable
            and you should define it as an observable with the 'FUNCTION' property
            """
            vardict = point_container.get_var_dict()
            derived = {
                var: eval(self.RunHandler.settings.derived_variables[var], vardict)
                for var in self.RunHandler.settings.derived_variables
            }
            point_container.var_dict.update(vardict)
        # if self.RunHandler.settings.store_points and (not ID):

        # This doesn't work (pickling)
        # self.last_point=point_container

        self.log.debug("Running point with input parameters: %s" % str(Point))

        point_container.valid = True
        os.chdir(self.temp_dir)

        if self.RunHandler.do_marginalisation:
            # print('Let us do fit')
            try:
                self.do_marginalisation(point_container)
                point_container.clear_obs()  # Entire toolchain will be run again with marginalised values
                point_container.valid = True
            except:  # failure to marginalise ...
                self.log.warning("Marginalisation failed")
                point_container.valid = False
                # if self.RunHandler.settings.store_points_in_memory:
                if self.RunHandler.settings.store_invalid_points:

                    point_container.marginal_variables = [0.0] * len(
                        self.RunHandler.marginal_variables
                    )
                    if self.RunHandler.settings.store_points_in_memory:
                        self.last_point_dict = point_container.to_dict()
                    """
                    self.last_point_dict = point_container.to_dict()
                    point_container.clear_obs()
                    point_container.valid = False
                    self.last_point_dict = point_container.to_dict()
                    self.last_observables=[]
                    point_all_dict=point_container.get_all_dict()
                    for obs,val in zip(self.RunHandler.settings.observables.keys(),self.RunHandler.settings.observables.values()):
                            #self.last_observables.append('NaN') # for invalid points we will be missing observables, declare them as NaN
                            self.last_observables.append(math.nan)
                    """
                else:
                    raise
            ## Now the results should be stored in the point_container.var_dict

        if (
            not self.RunHandler.settings.no_slha_output
        ):  # clear output files if we need them
            if os.path.exists(self.spc_file):
                try:
                    os.remove(self.spc_file)
                except Exception:
                    self.log.debug("Failed to remove file", exc_info=True)

        if (
            point_container.valid
        ):  # haven't tried to do marginalisation or it didn't fail
            if (
                not self.RunHandler.settings.no_slha_input
            ):  # check if we need an input file, and create it

                if os.path.exists(self.RunHandler.settings.input_filename):
                    try:
                        os.remove(self.RunHandler.settings.input_filename)
                    except Exception:
                        self.log.debug("Failed to remove file", exc_info=True)

                self.RunHandler.write_inputfile(
                    Point,
                    self.temp_dir,
                    self.RunHandler.settings.input_filename,
                    point_container.var_dict,
                )

                if self.RunHandler.settings.store_inputs:
                    shutil.copyfile(
                        self.RunHandler.settings.input_filename,
                        self.RunHandler.input_filename_prefix
                        + str(tID)
                        + self.RunHandler.settings.input_file_suffix,
                    )

            # moved this up, here we check whether to go on because if the marginalisation failed
            # point_container.valid=True
            for code in self.RunHandler.tools:
                # print(self.run_tools)
                if self.run_tools[code.name]:
                    self.log.debug("Launching " + str(code))
                    # point_container.valid=False
                    try:
                        code.launch(
                            self.spc_file, self.temp_dir, self.log, point_container
                        )
                        # point_container.valid=True
                    except Exception:
                        self.log.debug(
                            self.name + " failed to run " + code.name, exc_info=True
                        )
                        point_container.valid = False
                        # if self.last_point_dict is not None and self.RunHandler.settings.store_points_in_memory: # try to keep whatever observables we already gleaned
                        if self.RunHandler.settings.store_points_in_memory:
                            self.last_point_dict = point_container.to_dict()
                            # self.last_point_dict = point_data
                        break
                        # raise #### return with error

        if (
            not point_container.valid
            and not self.RunHandler.settings.store_invalid_points
        ):
            # no need to continue
            raise

        """

        if self.RunHandler.settings.no_slha_output:
            try:
                result=self.RunHandler.postprocess(Point,self.temp_dir, self.log, lock)
                self.last_observables = None
            except:
                self.log.debug(self.name+' failed to produce point')
                raise

        else:
        """

        if not self.RunHandler.settings.no_slha_output:  
            # read in the spc if we expect it and haven't already done so

            ## slha output -> read slha file first and extract observables
            if point_container.spc is None:
                # self.last_slha=zslha.read(self.output_file)
                if os.path.isfile(self.spc_file):
                    point_container.spc = zslha.read(self.spc_file)
            else:
                self.last_slha = point_container.spc

                # self.RunHandler.settings.observables contains all observables.
                # point_container should already have some values loaded into it. We just need to update them if not present

                # self.last_observables=[self.last_slha.Value(obs['SLHA'][0], obs['SLHA'][1]) for obs in self.RunHandler.settings.observables.values()]

        if len(self.RunHandler.settings.observables) > 0:
            try:
                self.last_observables = []
                point_all_dict = point_container.get_all_dict()
                for obs, val in zip(
                    self.RunHandler.settings.observables.keys(),
                    self.RunHandler.settings.observables.values(),
                ):
                    if (
                        obs in point_container.obs_dict
                    ):  # Trust that the correct value has been entered by the tools into the point container
                        self.last_observables.append(point_container.obs_dict[obs])
                    elif (
                        point_container.valid
                    ):  # observables that are not specified in the tool, but directly in the settings: extract these after all tools have run.
                        obsval = 0.0
                        if (
                            "SLHA" in val
                        ):  # New: allow for observables that are not SLHA values but may be functions.
                            obsval = float(
                                point_container.spc.Value(
                                    val["SLHA"][0], val["SLHA"][1]
                                )
                            )
                        elif "FUNCTION" in val:
                            obsval = eval(
                                val["FUNCTION"], point_all_dict
                            )  # allow for derived variables
                        else:
                            raise NameError(
                                "Unknown observable information: "
                                + str(obs)
                                + ", "
                                + str(val)
                            )
                        self.last_observables.append(obsval)
                        point_all_dict[obs] = obsval
                        point_container.obs_dict[obs] = obsval
                    else:
                        # self.last_observables.append('NaN') # for invalid points we will be missing observables, declare them as NaN
                        self.last_observables.append(math.nan)
            except Exception as e:
                self.log.debug(
                    self.name
                    + " point "
                    + str(tID)
                    + " failed to read observables, "
                    + str(e)
                )
                point_container.valid = False
                if not self.RunHandler.settings.store_invalid_points:
                    raise

        if point_container.valid:  # point ran successfully
            try:
                # result=self.RunHandler.postprocess(Point,self.last_observables,self.last_slha,self.temp_dir, self.log, lock)
                result = self.RunHandler.postprocess(
                    Point,
                    self.last_observables,
                    point_container,
                    self.temp_dir,
                    self.log,
                    lock,
                )
                point_container.outputs = result
            except:
                self.log.debug(
                    self.name + " point " + str(tID) + " failed postprocessing"
                )
                point_container.valid = False
                if not self.RunHandler.settings.store_invalid_points:
                    raise
        else:
            result = self.RunHandler.settings.invalid_return_value

        if (
            self.RunHandler.settings.csv_output
            or self.RunHandler.settings.tabbed_output
        ):
            # output_fileheaders will be a list with [0] being the filename and [1] being the separator.

            if (
                point_container.valid and not self.RunHandler.output_fileheaders.empty()
            ) or (
                not point_container.valid
                and not self.RunHandler.invalid_fileheaders.empty()
            ):
                headerline = []
                if tID:
                    headerline.append("ID")

                if self.RunHandler.settings.variables is not None:
                    for v in self.RunHandler.settings.variables:
                        headerline.append(v)
                else:
                    for x in range(len(Point)):
                        headerline.append("VARIABLE[" + str(x) + "]")

                if self.RunHandler.do_marginalisation:
                    for v in self.RunHandler.marginal_variables:
                        headerline.append(v)

                for obs in self.RunHandler.settings.observables:
                    headerline.append(obs)

                if point_container.valid:
                    if hasattr(result, "__len__") and len(result) > 1:
                        for x in range(len(result)):
                            headerline.append("RESULT[" + str(x) + "]")
                    else:
                        if str(result) != "":
                            headerline.append("Result")

                if point_container.valid:
                    getlist = self.RunHandler.output_fileheaders
                else:
                    getlist = self.RunHandler.invalid_fileheaders
                # while not self.RunHandler.output_fileheaders.empty():
                #    header_to_write=self.RunHandler.output_fileheaders.get()
                while not getlist.empty():
                    header_to_write = getlist.get()
                    # header_to_write is [filename,separator]
                    headerlinestring = header_to_write[1].join(headerline)

                    ## check to see if file already exists and prepend the header
                    if lock:
                        lock.acquire()
                    if not os.path.isfile(header_to_write[0]):
                        debug.command_line_log(
                            "echo "
                            + headerlinestring
                            + " >> "
                            + str(header_to_write[0]),
                            self.log,
                        )
                    else:  # either we have some results already, or we are resuming an old scan
                        with open(header_to_write[0], "r+") as f:
                            s = f.read()
                            if not s.startswith(headerlinestring):
                                f.seek(0)
                                f.write(headerlinestring + "\n" + s)
                    if lock:
                        lock.release()
            # end if not self.RunHandler.output_fileheaders.empty():
            outline = []

            if tID:
                outline.append(str(tID))

            # outline.extend(map(str,Point))
            outline.extend(map(self.RunHandler.format_variable, Point))

            if self.RunHandler.do_marginalisation:
                # outline.extend(map(str,point_container.marginal_variables))
                outline.extend(
                    map(
                        self.RunHandler.format_variable,
                        point_container.marginal_variables,
                    )
                )

            # outline.extend(map(str,self.last_observables))
            if self.last_observables is not None:
                outline.extend(
                    map(self.RunHandler.format_observable, self.last_observables)
                )

            if hasattr(result, "__len__") and len(result) > 1:
                # outline.extend(map(str,result))
                outline.extend(map(self.RunHandler.format_observable, result))
            else:
                if str(result) != "":
                    # outline.append(str(result))
                    outline.append(self.RunHandler.format_observable(result))

            if self.RunHandler.settings.csv_output:
                outlinestr = ",".join(outline)
                if point_container.valid:
                    outf = self.RunHandler.csv_filename
                else:
                    outf = self.RunHandler.invalid_csv_filename
                if lock:
                    lock.acquire()
                debug.command_line_log(
                        "echo " + outlinestr + " >> " + outf, self.log
                    )
                    # debug.command_line_log("echo " + outlinestr + " >> " + self.RunHandler.csv_filename, self.log)

                if lock:
                    lock.release()

            if self.RunHandler.settings.tabbed_output:
                outlinestr = "\t".join(outline)
                if point_container.valid:
                    outf = self.RunHandler.tabbed_filename
                else:  # we already quit if we are not storing data
                    outf = self.RunHandler.invalid_tabbed_filename

                if lock:
                    lock.acquire()
                debug.command_line_log(
                        "echo " + outlinestr + " >> " + outf, self.log
                    )
                    # debug.command_line_log("echo " + outlinestr + " >> " + self.RunHandler.tabbed_filename, self.log)
                if lock:
                    lock.release()

        if self.RunHandler.settings.store_everything:
            ## copy everything except the log which is called <name>, and anything else that might be onerous

            shutil.copytree(
                self.temp_dir,
                self.RunHandler.output_directory_prefix + str(tID),
                ignore=shutil.ignore_patterns("*.out", self.name),
            )

        if self.RunHandler.settings.store_outputs:

            # check first whether the point is valid. If so, then there must be an spc file.
            # If the point is not valid, only store the file if we are told to -- and if it exists!
            if point_container.valid or (
                self.RunHandler.settings.store_invalid_points
                and os.path.isfile(self.spc_file)
            ):
                if self.RunHandler.settings.store_separate_files:
                    shutil.copyfile(
                        self.spc_file,
                        self.RunHandler.output_filename_prefix
                        + str(tID)
                        + self.RunHandler.settings.output_file_suffix,
                    )

                else:
                    if lock:
                        lock.acquire()
                    debug.command_line_log(
                        "cat "
                        + self.spc_file
                        + " >> "
                        + self.RunHandler.output_filename,
                        self.log,
                    )
                    debug.command_line_log(
                        'echo "ENDOFPARAMETERPOINT" >> '
                        + self.RunHandler.output_filename,
                        self.log,
                    )
                    if lock:
                        lock.release()

        # if self.RunHandler.settings.store_points and ((not ID) or (make_point_container)):
        #    point_container.outputs=result
        #    point_container.valid=True

        if (
            self.last_observables is not None
        ):  ## if no slha then no observables to store ...
            point_container.observables = self.last_observables

        # if not point_container.valid and not self.RunHandler.settings.store_invalid_points:
        #    raise

        # point_container.valid=True
        # if self.last_point_dict is not None and self.RunHandler.settings.store_points_in_memory:
        if self.RunHandler.settings.store_points_in_memory:
            self.last_point_dict = point_container.to_dict()
            # self.last_point_dict = point_data
        if not point_container.valid:
            raise

        # self.last_point=point_container
        # self.last_point_dict = point_container.to_dict() # by default not including the spc file

        return result  # this is for scans that just want a function that returns a number or some quantity

    def do_marginalisation(self, point_container):
        import scipy

        # print('hello!')
        # var_dict=point_container.make_var_dict(self.RunHandler.settings.varnames)
        # var_dict=point_container.get_var_dict(self.RunHandler.settings.varnames)
        var_dict = point_container.var_dict

        # print(var_dict)
        def fun(x0):  ## function to pass to optimisation routine
            try:
                point_container.clear_obs()
                # print("h1"+str(x0))
                if os.path.exists(self.RunHandler.settings.input_filename):
                    try:
                        os.remove(self.RunHandler.settings.input_filename)
                    except Exception:
                        self.log.debug("Failed to remove file", exc_info=True)

                if os.path.exists(self.spc_file):
                    try:
                        os.remove(self.spc_file)
                    except Exception:
                        self.log.debug("Failed to remove file", exc_info=True)
                # if len(self.RunHandler.marginal_variables) == 1: # special handling for 1d case
                if self.RunHandler.oned_marginal:
                    for var, x in zip(self.RunHandler.marginal_variables, x0):
                        var_dict[var] = x
                else:
                    for var, x, f in zip(
                        self.RunHandler.marginal_variables,
                        x0,
                        self.RunHandler.marginal_downscalers,
                    ):
                        var_dict[var] = f(x)

                if (
                    "Overloads" in self.RunHandler.marginalisation_options
                ):  # allow overriding values in the input file, e.g. turning off decays in SPheno
                    self.RunHandler.write_inputfile(
                        point_container.inputs,
                        self.temp_dir,
                        self.RunHandler.settings.input_filename,
                        var_dict,
                        self.RunHandler.marginalisation_options["Overloads"],
                    )
                else:
                    self.RunHandler.write_inputfile(
                        point_container.inputs,
                        self.temp_dir,
                        self.RunHandler.settings.input_filename,
                        var_dict,
                    )

                for ncode in range(self.RunHandler.number_tools_run_marginalisation):
                    code = self.RunHandler.tools[ncode]
                    self.log.debug("Launching " + str(code))
                    try:
                        code.launch(
                            self.spc_file, self.temp_dir, self.log, point_container
                        )
                    except Exception:
                        self.log.debug(
                            self.name + " failed to run " + code.name, exc_info=True
                        )
                        raise  #### return with error

                if self.RunHandler.extra_observable_step:
                    if point_container.obs_dict is None:
                        point_container.obs_dict = OrderedDict()

                    for obs, val in zip(
                        self.RunHandler.settings.observables.keys(),
                        self.RunHandler.settings.observables.values(),
                    ):
                        # point_container.obs_dict[obs] = float(point_container.spc.Value(val['SLHA'][0], val['SLHA'][1]))
                        if (
                            "SLHA" in val
                        ):  # New: allow for observables that are not SLHA values but may be functions.
                            point_container.obs_dict[obs] = float(
                                point_container.spc.Value(
                                    val["SLHA"][0], val["SLHA"][1]
                                )
                            )
                        elif "FUNCTION" in val:
                            point_container.obs_dict[obs] = eval(
                                val["FUNCTION"], point_container.get_all_dict()
                            )  # allow for derived variables

                out_observables = [
                    float(eval(x, point_container.obs_dict))
                    for x in self.RunHandler.marginal_observables
                ]

                res = self.RunHandler.marginal_chisquared(out_observables)
                # print('res: '+str(res))
                return res
            except:
                return 1e32

        # result=scipy.optimize.minimize(fun, x0, args=(), method=None, jac=None, hess=None, hessp=None, bounds=None, constraints=(), tol=None)
        # print('Try optimize')
        try:
            # result=scipy.optimize.minimize(fun, [0.0]*len(self.RunHandler.marginal_variables), tol=1.0)
            # result=scipy.optimize.minimize(fun, [0.0]*len(self.RunHandler.marginal_variables), method=self.RunHandler.marginalisation_options, options={'xatol': 1,'gtol': 1e-22, 'disp': True})
            # result=scipy.optimize.minimize(fun, [0.0]*len(self.RunHandler.marginal_variables), method=self.RunHandler.marginalisation_method, options=dict(eval(self.RunHandler.marginalisation_options)))

            if self.RunHandler.oned_marginal:
                import bsmart.marginalisation as marginalisation

                result = marginalisation.oned_usebrentq(
                    fun,
                    self.RunHandler.marginal_ranges,
                    self.log,
                    options=self.RunHandler.marginalisation_options,
                )
            else:

                result = scipy.optimize.minimize(
                    fun,
                    [0.0] * len(self.RunHandler.marginal_variables),
                    method=self.RunHandler.marginalisation_method,
                    options=self.RunHandler.marginalisation_options,
                )

            # result=scipy.optimize.basinhopping(fun, [0.0]*len(self.RunHandler.marginal_variables))
        except Exception as e:
            self.log.debug("Failed to optimize:" + str(e))
            raise
        # print('ran'+str(result.x))
        results = []

        if self.RunHandler.oned_marginal:
            results.append(result)
            for var in self.RunHandler.marginal_variables:
                point_container.var_dict[var] = result
        else:

            for var, x, f in zip(
                self.RunHandler.marginal_variables,
                result.x,
                self.RunHandler.marginal_downscalers,
            ):
                res = f(x)
                # print(var+':'+str(res))
                try:
                    # print(point_container.var_dict)
                    point_container.var_dict[var] = res
                except Exception as e:
                    self.log.debug("Failed to set point container " + str(e))
                results.append(res)
        point_container.marginal_variables = results
        # print('Results: '+str(results))

        return
