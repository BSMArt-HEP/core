Installation
============

From version 2.0 `BSMArt` is installable as a package, requiring python 3.8 and higher.
A legacy version for older machines compatible with python 3.5 will remain available but is not covered by this documentation. 

Since `BSMArt` is now installable as a package, it can be installed using pip. We recommend using a virtual environment for this:

```bash
python3 -m venv BSMArt_venv
source BSMArt_venv/bin/activate
pip install bsmart
```
After installation, some shells (e.g. `tsch`) require the `rehash` command to update paths and find the new `BSMArt` executable. 

If installing as a user without a virtual environment, you may need to use the `--user` flag with `pip`, e.g.: `pip install --user bsmart`.

Many scans/tools require additional packages to be installed. To install as complete a list as possible, you can use:

```bash
pip install bsmart[all]
```
Attention: this will include some large packages such as `pytorch`. One hot tip, if you already have `pytorch` installed globally, is to create the virtual environment with `python3 -m venv --system-site-packages BSMArt_venv`. The disadvantage of this is that there may be conflicts with the global installation (not from `BSMArt` itself, but perhaps among dependencies). So this can be changed later by editing `BSMArt_venv/pyvenv.cfg` and setting `include-system-site-packages = false`.


The cleaner option is to use `BSMArt-CheckPackages` to check which packages are needed for your scan/tools and install them:
```bash
BSMArt-CheckPackages [--install] <jsonfile>
```
or specify which packages you want to check:
```bash
BSMArt-CheckPackages [--install] --scan <scan1> [scan2 ...] --tool <tool1> [tool2 ...]
```
The flag `--install` will install the packages if they are not installed.

You can then proceed to install necessary HEP tools using `BSMArt-InstallHEPTools` (see next section) or skip ahead to building examples



Installing HEP tools
---------------------

To install the necessary tools, use `BSMArt-InstallHEPTools`:

```bash
BSMArt-InstallHEPTools --All
```
Alternatively individual tools can be installed using e.g.:

```bash
BSMArt-InstallHEPTools --SPheno --MicrOMEGAs --HiggsTools --VevaciousPlusPlus --SModelS
```

This will download and install the necessary tools, creating a file `HEPtoolpaths.json`. By default, if you are using a virtual environment, this file will be created in the virtual environment, in 



```bash
BSMArt_venv/etc/bsmart/HEPtoolpaths.json
```

Otherwise, it will be created in the current directory; this can be forced by the option `--local`. It is a simple dictionary of paths:

```json
{
    "MicrOMEGAs": "/path/to/MicrOMEGAs",
    "SPheno": "/path/to/SPheno",
    "HiggsTools": "/path/to/HiggsTools",
    "VevaciousPlusPlus": "/path/to/VevaciousPlusPlus",
    "SModelS": "/path/to/SModelS"
}
```

If you have already installed the necessary tools for your scan, you can just create a `HEPtoolpaths.json` file with the paths to the tools and place it in the current directory or configuration folder and skip this step.

Quickstart
----------

A quickstart script is provided to help you get started. It will prepare the MSSM using SARAH and BSMArt, including all installed tools, and create example scans (one lightning fast, the other longer):

```bash
BSMArt-QuickStart
```

This will create a directory `BSMArt_QuickStart` with two example scans, `QuickStart` and `Lightning`. You can then launch them by typing:

```bash
cd BSMArt_QuickStart
BSMArt QuickStart_MSSM.json
```

or

```bash
cd BSMArt_QuickStart
BSMArt Lightning_MSSM.json
```

In case of problems, it is recommended to use the `--debug` flag to see warning messages; otherwise, the logs will be stored on either `/dev/shm/BSMArt_Temp` (if `/dev/shm` exists) or in subdirectory `Temp`.


Example scans
-------------

There are several example scans available included in the package. A subset of these are buildable so that you can run them directly on your system: you can do this via

```bash
BSMArt-BuildExamples
```

This will launch `SARAH` to build code for several different models, and create several scan directories containing scans over them. It will also create some examples with toy tools (not requiring `SPheno`). 

Community examples on github
----------------------------

To share scans, tools and example `json` files, we created a new github repository:

[https://github.com/BSMArt-HEP/examples](github.com/BSMArt-HEP/examples)

Instructions can be found there about downloading and using them for your own scans. 

We strongly encourage you to contribute by adding your own scans, tools and example `json` files to this repository via a pull request!


Setting up your own scan
------------------------

If you need to run `SARAH` for your model to generate the different codes, the `BSMArt-PrepareModel` script can be used. It will prepare the model using SARAH and configure it for all the included codes:

```bash
BSMArt-PrepareModel --All <modelname>
```

Alternatively, only a subset of codes can be configured by using different flags. 

This script will generate the code, compile the necessary outputs, and create a `BSMArt_modelname` directory with the necessary files to run scans, including an empty template `json` file and a template input file for SPheno, if required. Once you have configured the scan (according to the type of scan, tools desired, parameters to scan over, plots to make, etc) you can launch it with:

```bash
cd BSMArt_<modelname>
BSMArt <template>.json
```

and that's it! If you have problems, add the --debug tag to see warning messages and keep more log info in a Temp directory; otherwise, the logs will be stored on either /dev/shm/BSMArt_Temp (if /dev/shm exists) or in subdirectory Temp. This behaviour can be changed by setting the environment variable `BSMArt_Temp` to a different directory.



