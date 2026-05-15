Json Reference
==============

The BSMArt JSON input file is the main input file for a scan. It is a JSON file that contains all the necessary information for a scan. The JSON file is divided into several sections, each of which contains information about a specific aspect of the scan. 

The JSON file is divided into the following sections:

* `Codes`: Contains information about the codes to be run.
* `Setup`: Contains options to be passed to the scan. This includes information about the running of the tools, as well as scan specific settings. 
* `Variables`: Contains information about the variables to be scanned. The keys of this section are the names of the variables, and the values are dictionaries containing information about the variable. These can be as simple as defining a range, as in the schematic below.
* `Derived Variables`: Contains information about variables that are derived from the scanned variables via python functions, including mathematical operations such as `math.exp(variable)`.
* `Observables`: Observables that are collected at the end of running all tools. This is for backward compatibility; it is generally more useful to define observables within each tool's entry in the `Codes` section.
* `Fitting`: Contains information about fitting a selection of variables (not listed in the `Variables` section) based on a selection of observables. For example, fitting a quartic coupling to the Higgs mass. 
* `Plotting`: Contains information about the plots to be generated.

A typical JSON file is schematically like:
```json
{
    "Codes": {
        "SPheno": {
            "Command": "/path/to/SPheno/bin/SPhenoMSSM",
            "InputFile": "fast.in",
            "OutputFile": "SPheno.spc.MSSM",
            "Observables": {
                "mh": {
                    "SLHA": ["MASS",[25]],
                    "SCALING": "GAUSS",
                    "MEAN": 125.09,
                    "VARIANCE": 3.0
                }
            },
            "Run": "True"
        }
    },
    "Setup": {
        "RunName": "MyRun",
        "Type": "Grid",
        "Cores": 1,
        "StoreEverything": "False",
        "StoreAllPoints": "True",
        "Store In Memory": "False",
        "Store Invalid": "False",
        "StoreSeparateFiles": "False",
        "Store Inputs": "False",
        "csv": "True",
        "Precision": 8,
        "Spectrum File": "SPheno.spc.MSSM",
        "OutputFile": "SLHA_output"
  },
    "Variables": {
        "m0" : { "RANGE": [200,4000]},
        "m12" : { "RANGE": [200,4000]}		
    },
    "Derived Variables": {
        "m122": "m12**2"
    },
    "Observables": {},
    "Fitting": {
        "Variables": { "mA" : {"RANGE": [-6000,6000]} },
        "Observables": {"mh" :{"TARGET": 125,"VARIANCE":0.5}},
        "Options" : {"eps": 0.1, "ftol" : 1.0, "disp" : "False", "Overloads": [["SPHENOINPUT",[11],0]]}
    },
    "Plotting": {
        "Strategy": "csv",
        "Plots": {
            "m0_vs_mh": {
                "Labels": ["$m_0$ (GeV)","$m_{h}$ (GeV)"],
                "Values": ["m0","mh"],
                "Ranges": [[200,3500],[118,128]]
            }
        }
    }
}
```


Running options
---------------

The following options are available for running tools in the `Setup` section of the JSON file:

* RunName: Name of the run. All output files will be placed in a folder named of this name.
* Type: The type of scan to run, e.g. Grid, Random, MCMC, ...
* Cores: Number of cores to use, for parallelisation via python multiprocessing. Behaviour will differ (or be ignored) for running via `mpi`.
* TempDir: Temporary working directory.
* debug: Enable debug mode (boolean).
* SLHA_Output: Boolean, whether to generate SLHA files (True by default).
* Skip SLHA Files: Boolean, for pure python scans (skips I/O).
* Output File: Name for the final output file. This will be used to name the final outputs: the csv/tabbed output will be named this plus `.csv` or `.tab`. If the SLHA files are stored separately, they will be named `<Output File>_1`, `<Output File>_2`, etc. If the SLHA files are stored together, they will just be named `<Output File>`. 
* Spectrum File: Name for *intermediate* spectrum file, passed between the different codes in the temporary run directory. If this option is absent it will be taken by default to be the OutputFile of the first code; for example, if SPheno is the first tool, then it would be `SPheno.spc.<Modelname>`. Can be set to anything except a blank string!
* StoreInputs: Boolean, whether to store input Les Houches files.
* StoreAllPoints: Boolean, whether to store all points (default True).
* Store Output Files: Boolean, alias for StoreAllPoints.
* Short: Boolean, enable tabbed output.
* csv: Boolean. If True, variables and observables will be stored in a csv file.
* csv name: Name for the CSV file.
* Precision: Integer, precision for csv output.
* Store Invalid: Boolean. If True, invalid points from the scan will be stored, with data kept in a separate csv file.
* Store In Memory: Boolean, store DataPoint objects in memory. Necessary for extracting information from failed points (when using hierarchical losses).
* StoreSeparateFiles: Boolean. If True, the spectrum file for each point will be stored in a separate file.
* StoreEverything: Boolean, If True, all data from each point (input files and outputs) will be stored in separate directories within the `Spectrum_Files` subdirectory.
* Use TQDM: Boolean, use tqdm progress bar.



Fitting
-------

The `Fitting` section is used to fit a selection of variables to a selection of observables. 

Variables to be marginalised over should *not* appear in the `Variables` section of the JSON file.

The `Observables` section in the `Fitting` section *must* already appear as an observable from a tool or in the global `Observables` section. Here we require a `TARGET` value and a `VARIANCE` (to search for values within that tolerance) which is information for the fitting routine only (not likelihood computations).

Options:

* eps: Epsilon for the minimiser (default 0.1)
* ftol: Tolerance for the minimiser (default 1.0)
* disp: Display the minimiser output (default False)
* Overloads: allows overloading of options in the Les Houches Input files during the minimisation process,  e.g. [["SPHENOINPUT",[11],0]] tells the minimisation routine to change the SLHA input value 11 to 0 during this step. After the minimisation is complete, the original value will be restored. This is particularly useful for turning off slow features of the code, for example disabling calculation of decay spectra from SPheno when we are fitting values to the Higgs mass (as in the example given here).