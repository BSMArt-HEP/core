#!/usr/bin/env python3

"""
Quickstart Script to prepare the MSSM using SARAH and BSMArt, including all installed tools,
and create example scans (one lightning fast, the other longer)
"""


import os
from datetime import datetime
import sys

import shutil
import json
import collections

import argparse

import subprocess
import importlib.resources
from bsmart.scripts import prepare_model

def main():

    parser = argparse.ArgumentParser(
        description="Quickstart Script to prepare the MSSM using SARAH and BSMArt."
    )
    # Add any arguments if needed, for now just help is sufficient
    args_parsed = parser.parse_args()


    ## Load the paths

    cwd = os.getcwd()
        

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
            print('Failed to load json file '+jsonFILE)
            print('Json exception given as ' +str(e))
            raise SystemExit
    else:
         print('Failed to load json file '+jsonFILE)
         print('Please run BSMArt-InstallHEPTools first or ensure HEPtoolpaths.json is in the current directory or configuration folder.')
         raise SystemExit

    ### Now run the script

    # command='python3 prepareModel.py --debug --SPheno --BSMArt --MicrOMEGAs --HiggsTools --HiggsBounds --HiggsSignals --VevaciousPlusPlus'
    ## For debugging use ReadLists->True
    # command=command + ' --SkipFlavorKit --Options="IncludeLoopDecays->False" MSSM' 

    # We call prepare_model.main() via simulating args
    
    
    # But we need to handle quoting of Options properly if passing via sys.argv manually
    # argparse handles quoting from shell. Here I am constructing the list directly.
    # '--Options="IncludeLoopDecays->False"' might pass "IncludeLoopDecays->False" (with quotes) or not depending on parsing.
    # Usually we pass '--Options=IncludeLoopDecays->False'.
    
    sys.argv = ['prepare_model.py', '--debug', '--SPheno', '--BSMArt', '--MicrOMEGAs', '--HiggsTools', '--HiggsBounds', '--HiggsSignals', '--VevaciousPlusPlus', '--SkipFlavorKit', '--Options=IncludeLoopDecays->False', 'MSSM']
    
    try:
        prepare_model.main()
    except SystemExit as e:
        if e.code is not None and e.code != 0:
             print('Error preparing model!')
             raise e
    except Exception as e:
        print('Error preparing model! '+str(e))
        raise SystemExit


    ### Now start moving stuff around and compiling



    print("Setting up the QuickStart files")

    BSMArtScanDir=os.path.join(cwd,'BSMArt_'+'MSSM')
    QuickStartScanDir=os.path.join(cwd,'BSMArt_QuickStart')



    if os.path.exists(QuickStartScanDir):
        shutil.rmtree(QuickStartScanDir)

    # QuickStartTemplate=os.path.join(cwd,'QuickStart.in.MSSM')
    # LightningTemplate=os.path.join(cwd,'lightning.in.MSSM')
    
    # Use resources
    try:
        from importlib.resources import files
        QuickStartTemplate = files('bsmart.data').joinpath('QuickStart.in.MSSM')
        LightningTemplate = files('bsmart.data').joinpath('lightning.in.MSSM')
    except:
         print("Could not find template files in package!")
         raise SystemExit

    if not os.path.exists(QuickStartScanDir):    
        os.makedirs(QuickStartScanDir)

    # Copy files
    if QuickStartTemplate.is_file():
         with importlib.resources.as_file(QuickStartTemplate) as f:
             shutil.copy(f,QuickStartScanDir)
    if LightningTemplate.is_file():
         with importlib.resources.as_file(LightningTemplate) as f:
             shutil.copy(f,QuickStartScanDir)

    shutil.copy(os.path.join(BSMArtScanDir,'QNUMBERS_MSSM.slha'),QuickStartScanDir)

    modelname='MSSM'
    internalmodelname='MSSM'
    #internalmodelname=str(bsmart_dict['BSMArtmodelname'])

    Qtemplate='QuickStart.in.'+internalmodelname
    Ltemplate='lightning.in.'+internalmodelname


    """
    Now set up the json files, starting from the one created by prepareModel!
    """


    try:
        with open(os.path.join('BSMArt_MSSM','BSMArt_MSSM.json')) as template_json_data:
            template_dict = json.load(template_json_data, object_pairs_hook=collections.OrderedDict)
    except Exception as e:
        #log.error('Failed to load json file '+file)
        print('Failed to load json file '+'BSMArt_MSSM/BSMArt_MSSM.json')
        print('Json exception given as ' +str(e))
        raise SystemExit


    """
    Quickstart -> full scan
    lightning -> very quick scan

    """
    quickstart_dict=template_dict # could do .copy()


    ### first the full scan

    if 'SPheno' in quickstart_dict['Codes']:
        quickstart_dict['Codes']['SPheno']['Run']="True"
        quickstart_dict['Codes']['SPheno']['InputFile']=Qtemplate
        quickstart_dict['Codes']['SPheno']['Observables']={"mh" : { "SLHA": ["MASS", [25]], "SCALING": "GAUSS", "MEAN":125.0, "VARIANCE": 3.0}}


    if 'HiggsBounds' in quickstart_dict['Codes'] and eval(quickstart_dict['Codes']['HiggsBounds']["Run"]):
        quickstart_dict['Codes']['HiggsBounds']['Observables']={"HBRatio" : {"SLHA": ["HIGGSBOUNDS",[101]], "SCALING":"UPPER","MEAN":1.0,"VARIANCE":0.1},"HBOK" : {"SLHA": ["HIGGSBOUNDS",[5]],"SCALING":"OFF"}}

    if 'HiggsSignals' in quickstart_dict['Codes'] and eval(quickstart_dict['Codes']['HiggsSignals']["Run"]):
        quickstart_dict['Codes']['HiggsSignals']['Observables']={"HSpval" : {"SLHA": ["HIGGSSIGNALS",[101]], "SCALING":"LOWER","MEAN":0.001,"VARIANCE":0.1}}

    if 'HiggsTools' in quickstart_dict['Codes'] and eval(quickstart_dict['Codes']['HiggsTools']["Run"]):
        quickstart_dict['Codes']['HiggsTools']['Observables']={
            "HBresult" : {"SLHA": ["HIGGSTOOLS",[1]],"SCALING":"OFF"},
            "HBobsr" : {"SLHA": ["HIGGSTOOLS",[2]], "SCALING":"UPPER","MEAN":1.0,"VARIANCE":0.1},
            "HSchi2" : {"SLHA": ["HIGGSTOOLS",[11]],"SCALING":"OFF"},
            "HSnobs" : {"SLHA": ["HIGGSTOOLS",[12]],"SCALING":"OFF"},
            "HSpval" : {"SLHA": ["HIGGSTOOLS",[13]], "SCALING":"LOWER","MEAN":0.001,"VARIANCE":0.1}
            }

    if 'MicrOmegas' in quickstart_dict['Codes'] and eval(quickstart_dict['Codes']['MicrOmegas']["Run"]):
        quickstart_dict['Codes']['MicrOmegas']['Observables']={
            "Oh2" : { "SLHA": ["DARKMATTER", [1]],"SCALING":"LOG", "MEAN":0.112,"VARIANCE":0.05},
            "CDM1" : { "SLHA": ["DARKMATTER", [2]],"SCALING":"OFF"},
            "DMpval" : { "SLHA": ["DARKMATTER", [401]],"SCALING":"OFF"}
            }


    if 'SModelS' in quickstart_dict['Codes'] and eval(quickstart_dict['Codes']['MicrOmegas']["Run"]) and eval(quickstart_dict['Codes']['MicrOmegas']["Cross-sections"]):
        quickstart_dict['Codes']['SModelS']['Observables']={
            "SModelSr": {"SLHA":["SMODELS",[1]], "SCALING": "UPPER", "MEAN": 1.0, "VARIANCE": 0.2}, "SModelSrm1": {"SLHA":["SMODELS",[2]], "SCALING":"OFF"}
            }

    ## Create the json files
        

    if 'Vevacious++' in quickstart_dict['Codes'] and eval(quickstart_dict['Codes']['Vevacious++']["Run"]):
        quickstart_dict['Codes']['Vevacious++']['Observables']={
                "Vstab": {"SLHA": ["VEVACIOUSSTABILITY",[1]], "SCALING":"OFF"}
            }


    quicksetup={
            "RunName": "QuickStart",
            "Type": "MCMC",
            "csv": "True",
            "Cores": 1,
            "Merge Results": "True",
            "StoreEverything": "False",
            "Store In Memory": "True",
            "StoreAllPoints": "True",
            "StoreAllPoints": "True",
            "StoreSeparateFiles": "False",
            "Spectrum File": "SPheno.spc.MSSM",
            "Output File": "MSSM_out",
            "Points":100
    }

    for setting,value in quicksetup.items():
        quickstart_dict['Setup'][setting]=value

    quickstart_dict['Variables']={
        "m0": { "RANGE": [100,5000.0],"VARIANCE": 500.0},
        "m12": { "RANGE": [100,5000.0],"VARIANCE": 500.0}
        }

    quickstart_dict['Plotting']['Strategy']='SLHA'
    quickstart_dict['Plotting']['Plots']={
        "m0_vs_mh": { "Labels": ["$m_0$ (GeV)","$m_h$ (GeV)"], "Values": [["MINPAR",[1]],["MASS",[25]]], "Ranges":[[100,5000],[118,127]]},
        "m12_vs_mh": { "Labels": ["$m_{1/2}$ (GeV)","$m_h$ (GeV)"], "Values": [["MINPAR",[2]],["MASS",[25]]], "Ranges":[[100,5000],[118,127]]}
        }

    with open(os.path.join(QuickStartScanDir,'QuickStart_MSSM.json'), 'w') as fp:
        json.dump(quickstart_dict, fp, indent=4)



    """

    Now for the lightning run

    """
    lightning_dict=quickstart_dict


    ### Set all codes except SPheno to false
    for code in lightning_dict['Codes']:
        lightning_dict['Codes'][code]['Run']='False'

    lightning_dict['Codes']['SPheno']['Run']="True"
    lightning_dict['Codes']['SPheno']['InputFile']="lightning.in.MSSM"

    light_setup={
        "RunName": "Lightning",
        "Type": "Random",
        "csv": "True",
        "Cores": 1,
        "Merge Results": "True",
        "StoreEverything": "False",
        "StoreAllPoints": "True",
        "StoreSeparateFiles": "False",
    #    "Spectrum File": "SPheno.spc.MSSM",  ## don't need to change these, they are inherited
    #    "Output File": "MSSM_out",
        "Points":100
        }

    for setting,value in light_setup.items():
        lightning_dict['Setup'][setting]=value



    lightning_dict['Plotting']['Strategy']='csv'
    lightning_dict['Plotting']['Plots']={
        "m0_vs_mh": { "Labels": ["$m_0$ (GeV)","$m_h$ (GeV)"], "Values": ["m0","mh"], "Ranges":[[100,3000],[118,122]]},
        "m12_vs_mh": { "Labels": ["$m_{1/2}$ (GeV)","$m_h$ (GeV)"], "Values": ["m12","mh"], "Ranges":[[100,3000],[118,122]]}
        }

    with open(os.path.join(QuickStartScanDir,'lightning_MSSM.json'), 'w') as fp:
        json.dump(lightning_dict, fp, indent=4)







    print("All done!")

    print('To run the lightning Random MSSM scan, type:\n    cd BSMArt_QuickStart\n   BSMArt lightning_MSSM.json\n')
    print('To run the example MCMC scan, type:\n     cd BSMArt_QuickStart\n   BSMArt QuickStart_MSSM.json\n')

if __name__ == "__main__":
    main()
