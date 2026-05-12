#!/usr/bin/env python3
from __future__ import annotations
"""

Script to create or update template json files for BSMArt.

Use after prepareModel
This is used by the example scripts too, so you can have multiple scans per model without rerunning all the codes


The only relevant argument is the model name, which should correspond to a file in the ScriptsAndLogs directory (BSMArt_data_<modelname>.json)

2025/01/29

"""



import os
import sys

import shutil
import json
import collections

import argparse

#import compactJSON





class CompactJSONEncoder(json.JSONEncoder):
    """A JSON Encoder that puts small containers on single lines."""

    CONTAINER_TYPES = (list, tuple, dict)
    """Container datatypes include primitives or other containers."""

    MAX_WIDTH = 70
    """Maximum width of a container that might be put on a single line."""

    MAX_ITEMS = 10
    """Maximum number of items in container that might be put on single line."""

    def __init__(self, *args, **kwargs):
        # using this class without indentation is pointless
        if kwargs.get("indent") is None:
            kwargs["indent"] = 4
        super().__init__(*args, **kwargs)
        self.indentation_level = 0

    def encode(self, o):
        """Encode JSON object *o* with respect to single line lists."""
        if isinstance(o, (list, tuple)):
            return self._encode_list(o)
        if isinstance(o, dict):
            return self._encode_object(o)
        if isinstance(o, float):  # Use scientific notation for floats
            return format(o, "g")
        return json.dumps(
            o,
            skipkeys=self.skipkeys,
            ensure_ascii=self.ensure_ascii,
            check_circular=self.check_circular,
            allow_nan=self.allow_nan,
            sort_keys=self.sort_keys,
            indent=self.indent,
            separators=(self.item_separator, self.key_separator),
            default=self.default if hasattr(self, "default") else None,
        )

    def _encode_list(self, o):
        if self._put_on_single_line(o):
            return "[" + ", ".join(self.encode(el) for el in o) + "]"
        self.indentation_level += 1
        output = [self.indent_str + self.encode(el) for el in o]
        self.indentation_level -= 1
        return "[\n" + ",\n".join(output) + "\n" + self.indent_str + "]"

    def _encode_object(self, o):
        if not o:
            return "{}"

        # ensure keys are converted to strings
        o = {str(k) if k is not None else "null": v for k, v in o.items()}

        if self.sort_keys:
            o = dict(sorted(o.items(), key=lambda x: x[0]))

        if self._put_on_single_line(o):
            return (
                "{ "
                + ", ".join(
                    f"{json.dumps(k)}: {self.encode(el)}" for k, el in o.items()
                )
                + " }"
            )

        self.indentation_level += 1
        output = [
            f"{self.indent_str}{json.dumps(k)}: {self.encode(v)}" for k, v in o.items()
        ]
        self.indentation_level -= 1

        return "{\n" + ",\n".join(output) + "\n" + self.indent_str + "}"

    def iterencode(self, o, **kwargs):
        """Required to also work with `json.dump`."""
        return self.encode(o)

    def _put_on_single_line(self, o):
        return (
            self._primitives_only(o)
            and len(o) <= self.MAX_ITEMS
            and len(str(o)) - 2 <= self.MAX_WIDTH
        )

    def _primitives_only(self, o: list | tuple | dict):
        if isinstance(o, (list, tuple)):
            return not any(isinstance(el, self.CONTAINER_TYPES) for el in o)
        elif isinstance(o, dict):
            return not any(isinstance(el, self.CONTAINER_TYPES) for el in o.values())

    @property
    def indent_str(self) -> str:
        if isinstance(self.indent, int):
            return " " * (self.indentation_level * self.indent)
        elif isinstance(self.indent, str):
            return self.indentation_level * self.indent
        else:
            raise ValueError(
                f"indent must either be of type int or str (is: {type(self.indent)})"
            )

########################################################################################################


#def updateBSMArtTemplate(modelname,bsmart_dict,paths_dict,includecodes,template_dict,template_name):
def updateBSMArtTemplate(modelname,bsmart_dict,paths_dict,in_includecodes,template,template_name=None):
   
    with open(template,'r') as IF:
        template_dict = json.load(IF,object_pairs_hook=collections.OrderedDict)

        # need the name of the template file -> we'll keep it
    if template_name is None:
        template_path,template_name = os.path.split(template)
    else:
        template_path,_ = os.path.split(template)
    cwd=os.getcwd()
    includecodes={}
    for code in in_includecodes:
        if in_includecodes[code] and code in paths_dict:
            includecodes[code] = True
        else:
            includecodes[code] = False
    
    if 'VevaciousPlusPlus' in template_dict:
        includecodes['VevaciousPlusPlus'] = includecodes['Vevacious++']
        if 'SARAH' not in paths_dict:
            includecodes['VevaciousPlusPlus'] = False
    else:
        includecodes['VevaciousPlusPlus'] = False

    outjson=template_dict
    if 'SPhenoOutputDir' in bsmart_dict and 'BSMArtmodelname' in bsmart_dict and includecodes['SPheno']:
        sarahsphdir=bsmart_dict["SPhenoOutputDir"]
        internalmodelname=str(bsmart_dict['BSMArtmodelname'])
        #sphenotarget=os.path.join(paths_dict["SPheno"],modelname)
    else:
        sarahsphdir=None
        internalmodelname=modelname
        #sphenotarget=None
        includecodes['SPheno'] = False

    BSMArtScanDir=os.path.join(cwd,'BSMArt_'+modelname)

    if not os.path.exists(BSMArtScanDir):    
        os.makedirs(BSMArtScanDir)    

    if 'SPheno' in template_dict['Codes']:
        if includecodes['SPheno']:
            outjson['Codes']['SPheno']['Command'] = os.path.join(paths_dict['SPheno'],'bin','SPheno'+internalmodelname)


            #template='LesHouches.in.'+internalmodelname
            #template_low='LesHouches.in.'+internalmodelname+'_low'
            inputfile=os.path.basename(template_dict['Codes']['SPheno']['InputFile'])
            if inputfile== '':
                inputfile='LesHouches.in.'+internalmodelname
            paramfile=os.path.join(template_path,inputfile)

            if not os.path.exists(paramfile):  # take default :)
                paramfile =os.path.join(sarahsphdir,'Input_Files','LesHouches.in.'+internalmodelname)


            if os.path.exists(paramfile):
                shutil.copy(paramfile,BSMArtScanDir)
            outputfile=''
            #if 'FlexibleSUSY' in template_dict['Codes']:  # why???
            #    outputfile=os.path.basename(template_dict['Codes']['FlexibleSUSY'].get('OutputFile',''))

            outputfile=os.path.basename(template_dict['Codes']['SPheno'].get('OutputFile',''))
                
                
            if outputfile == '':
                outputfile='SPheno.spc.'+internalmodelname

            paramfile2=os.path.join(sarahsphdir,'Input_Files','LesHouches.in.'+internalmodelname+'_low')

            if os.path.exists(paramfile2):
                shutil.copy(paramfile2,BSMArtScanDir)

            QNUMBERS_slha='QNUMBERS_'+internalmodelname+'.slha'

            if os.path.exists(os.path.join(sarahsphdir,'Input_Files',QNUMBERS_slha)):
                shutil.copy(os.path.join(sarahsphdir,'Input_Files',QNUMBERS_slha),BSMArtScanDir) 
            else:
                QNUMBERS_slha='QNUMBERS_'+modelname+'.slha'
                if os.path.exists(os.path.join(cwd,'ScriptsAndLogs',QNUMBERS_slha)): ## in this case, older version of SARAH so no duplicate
                    shutil.copy(os.path.join(cwd,'ScriptsAndLogs',QNUMBERS_slha),BSMArtScanDir) 
            outjson['Codes']['SPheno']['InputFile']=inputfile
            outjson['Codes']['SPheno']['OutputFile'] = outputfile
            if 'Spectrum File' not in outjson['Setup'] or outjson['Setup']['Spectrum File'] == '':
                outjson['Setup']['Spectrum File'] = outputfile
            #QNUMBERS_slha=os.path.join(BSMArtScanDir,QNUMBERS_slha)
        else:
            #outjson['Codes']['SPheno']["Run"] = "False"
            outjson['Codes'].pop('SPheno') # delete!!

    if 'FlexibleSUSY' in template_dict['Codes']:
        if includecodes['FlexibleSUSY'] and 'FSmodelname' in bsmart_dict:
            FSmodelname=bsmart_dict['FSmodelname']
            FSmodeldir=os.path.join(paths_dict['FlexibleSUSY'],'models',FSmodelname)
            outjson['Codes']['FlexibleSUSY']['Command'] = os.path.join(FSmodeldir,'run_'+FSmodelname+'.x')


            #template='LesHouches.in.'+internalmodelname
            #template_low='LesHouches.in.'+internalmodelname+'_low'
            inputfile=os.path.basename(template_dict['Codes']['FlexibleSUSY']['InputFile'])
            if inputfile== '':
                inputfile='LesHouches.in.'+FSmodelname
            
            outputfile=os.path.basename(template_dict['Codes']['FlexibleSUSY'].get('OutputFile',''))
            if outputfile == '':
                outputfile='FlexibleSUSY.spc.'+FSmodelname
            

            paramfile=os.path.join(template_path,inputfile)
            #if os.path.exists(inputfile):
            #    paramfile=inputfile
            #else:
            if not os.path.exists(paramfile):
                paramfile =os.path.join(FSmodeldir,inputfile)


            if os.path.exists(paramfile):
                destparamfile=os.path.join(BSMArtScanDir,inputfile)
                if os.path.exists(destparamfile):
                    # if initial template actually contained the correct input then we are using both SPheno & FS, so ok
                    # Otherwise we need to rename
                    if template_dict['Codes']['FlexibleSUSY']['InputFile'] == '':
                        inputfile+='.FS'
                        destparamfile=os.path.join(BSMArtScanDir,inputfile)
                        shutil.copy(paramfile,destparamfile)
                else:
                    shutil.copy(paramfile,destparamfile)
            outjson['Codes']['FlexibleSUSY']['InputFile']=inputfile
            outjson['Codes']['FlexibleSUSY']['OutputFile'] = outputfile

            ## Cannot run SPheno and FlexibleSUSY at the same time by default (needs special care)
            if 'SPheno' in outjson['Codes'] and eval(outjson['Codes']['SPheno']['Run']):
                outjson['Codes']['FlexibleSUSY']['Run'] = "False"
        else:
            #outjson['Codes']['FlexibleSUSY']["Run"] = "False"
            outjson['Codes'].pop('FlexibleSUSY') # delete!

    #if args.HiggsBounds or args.HiggsSignals or args.All:
    if includecodes['HiggsBounds'] or includecodes['HiggsSignals']:
        if 'Hneutral' not in bsmart_dict or 'Hcharged' not in bsmart_dict:
            includecodes['HiggsBounds']=False
            includecodes['HiggsSignals']=False
        else:
            HBdata=str(bsmart_dict['Hneutral'])+' '+str(bsmart_dict['Hcharged'])
    
        


    if includecodes['MicrOMEGAs']:
        targetmodir=os.path.join(paths_dict["MicrOMEGAs"],'SARAH_'+modelname)
        #MOexec=os.path.join(targetmodir,'MicrOmegas_v5.3_BSMArt')

        if 'MicrOMEGAs executable' in bsmart_dict:
            
            MOexec=bsmart_dict['MicrOMEGAs executable']
        else:
            #print('Missing MO exec!')
        #if not os.path.exists(MOexec):
            includecodes['MicrOMEGAs'] = False

    ## Move vevacious files and create the init files
    
    #if args.VevaciousPlusPlus or args.All:

    if includecodes['VevaciousPlusPlus']:
        sarahVVdir=os.path.join(paths_dict["SARAH"],'Output',modelname,'Vevacious')

        if not os.path.exists(sarahVVdir):
            includecodes['VevaciousPlusPlus'] = False


    if includecodes['VevaciousPlusPlus']:
        sarahVVdir=os.path.join(paths_dict["SARAH"],'Output',modelname,'Vevacious')

        BSMVevaciousdir=os.path.join(BSMArtScanDir,'Vevacious')
        if os.path.exists(BSMVevaciousdir):
            shutil.rmtree(BSMVevaciousdir)
        
        os.makedirs(BSMVevaciousdir)
        shutil.copy(os.path.join(sarahVVdir,internalmodelname+'.vin'),BSMVevaciousdir)
        shutil.copy(os.path.join(sarahVVdir,'ScaleAndBlock.xml'),os.path.join(BSMVevaciousdir,internalmodelname+'.xml'))
        #shutil.copy(os.path.join(BSMVevaciousdir,'ScaleAndBlock.xml'),BSMVevaciousdir)
        #for ff in [internalmodelname+'.vin','ScaleAndBlock.xml']:
        #    shutil.copy(os.path.join(BSMVevaciousdir,ff),BSMArtScanDir)
        VVPotentialFile=os.path.join(BSMVevaciousdir,'PotentialFunctionInitialization.xml')

        with open(VVPotentialFile,'w') as VVP:
            VVP.write('<VevaciousPlusPlusPotentialFunctionInitialization>\n')
            VVP.write('<LagrangianParameterManagerClass>\n')
            VVP.write('<ClassType>\n SlhaCompatibleWithSarahManager\n </ClassType>\n')
            VVP.write('<ConstructorArguments>\n <ScaleAndBlockFile>\n')
            VVP.write(os.path.join(BSMVevaciousdir,internalmodelname+'.xml'))
            #VVP.write(os.path.join(BSMVevaciousdir,'ScaleAndBlock.xml'))
            VVP.write(' </ScaleAndBlockFile></ConstructorArguments>\n') 
            VVP.write('</LagrangianParameterManagerClass>\n')
            VVP.write('<PotentialFunctionClass>\n<ClassType> \n FixedScaleOneLoopPotential \n </ClassType>\n <ConstructorArguments>\n')
            VVP.write('<ModelFile>\n')
            VVP.write(os.path.join(BSMVevaciousdir,internalmodelname+'.vin'))
            VVP.write('</ModelFile>\n')
            VVP.write('<AssumedPositiveOrNegativeTolerance>\n 0.5 \n </AssumedPositiveOrNegativeTolerance> \n </ConstructorArguments>\n </PotentialFunctionClass>\n')
        
            VVP.write('</VevaciousPlusPlusPotentialFunctionInitialization>\n')

        ### NEW: Copy over the initialisation files instead of using the THDM ones
    
        VVPotentialInitFile=os.path.join(BSMVevaciousdir,'PotentialMinimizerInitialization.xml')
        shutil.copy(os.path.join(paths_dict['Vevacious++'],'InitializationFiles','THDMInitializationFiles','PotentialMinimizerInitialization.xml'),VVPotentialInitFile)

        ## edit the line to use 

        VVTunnelingFile=os.path.join(BSMVevaciousdir,'TunnelingCalculatorInitialization.xml')
        shutil.copy(os.path.join(paths_dict['Vevacious++'],'InitializationFiles','THDMInitializationFiles','TunnelingCalculatorInitialization.xml'),VVTunnelingFile)
        
        VVinitfile=os.path.join(BSMVevaciousdir,internalmodelname+'Initialization.xml')
        with open(VVinitfile,'w') as VVF:
            VVF.write('<VevaciousPlusPlusObjectInitialization>\n')
            VVF.write('  <PotentialFunctionInitializationFile>\n')
            VVF.write(VVPotentialFile+'\n') 
            VVF.write('  </PotentialFunctionInitializationFile>\n')
    
            VVF.write('  <PotentialMinimizerInitializationFile>\n')
            ## Use THDM
            VVF.write(VVPotentialInitFile+'\n')
            #VVF.write(os.path.join(paths_dict['Vevacious++'],'InitializationFiles','THDMInitializationFiles','PotentialMinimizerInitialization.xml')+'\n')
            #/ada1/lpthe/goodsell/Scratch/SPheno/Vevacious++/VevaciousPlusPlus/InitializationFiles/MSSMInitialization/PotentialMinimizerInitialization.xml
            VVF.write('  </PotentialMinimizerInitializationFile>\n')
    
            VVF.write('  <TunnelingCalculatorInitializationFile>\n')
            VVF.write(VVTunnelingFile+'\n')

            ## Use THDM
            #VVF.write(os.path.join(paths_dict['Vevacious++'],'InitializationFiles','THDMInitializationFiles','TunnelingCalculatorInitialization.xml')+'\n')
            #/ada1/lpthe/goodsell/Scratch/SPheno/Vevacious++/VevaciousPlusPlus/InitializationFiles/MSSMInitialization/TunnelingCalculatorInitialization.xml
            VVF.write('  </TunnelingCalculatorInitializationFile>\n')
            VVF.write('  </VevaciousPlusPlusObjectInitialization>\n')
            

        
    #OF2=open(os.path.join(BSMArtScanDir,'Template_'+modelname+'.json'),'w')
        
    #if args.HiggsBounds or args.All:
    if 'HiggsBounds' in template_dict['Codes']:
        if includecodes['HiggsBounds'] and 'HiggsBounds' in paths_dict:

            outjson['Codes']['HiggsBounds']['Command'] = os.path.join(paths_dict['HiggsBounds'],'HiggsBounds')

            if 'Neutral Higgs' not in template_dict['Codes']['HiggsBounds'] or template_dict['Codes']['HiggsBounds']['Neutral Higgs'] is None:
                outjson['Codes']['HiggsBounds']["Neutral Higgs"] = bsmart_dict['Hneutral']
                outjson['Codes']['HiggsBounds']["Charged Higgs"]= bsmart_dict['Hcharged']
        else:
            #outjson['Codes']['HiggsBounds']["Run"] = "False"
            outjson['Codes'].pop('HiggsBounds')

    #if args.HiggsSignals or args.All:
    if 'HiggsSignals' in template_dict['Codes']:
        if includecodes['HiggsSignals'] and 'HiggsSignals' in paths_dict:
            outjson['Codes']['HiggsSignals']['Command'] = os.path.join(paths_dict['HiggsSignals'],'HiggsSignals')

            
            if  'Neutral Higgs' not in template_dict['Codes']['HiggsSignals'] or template_dict['Codes']['HiggsSignals']['Neutral Higgs'] is None or template_dict['Codes']['HiggsTools']['Neutral Higgs'] == []:
                    outjson['Codes']['HiggsSignals']["Neutral Higgs"] = bsmart_dict['Hneutral']
                    outjson['Codes']['HiggsSignals']["Charged Higgs"]= bsmart_dict['Hcharged']
                    outjson['Codes']['HiggsSignals']["Options"]  = 'latestresults 2 effC '+HBdata
        else:
            #outjson['Codes']['HiggsSignals']["Run"] = "False"
            outjson['Codes'].pop('HiggsSignals')

    if 'Neutral Higgs' not in bsmart_dict or "Charged Higgs" not in bsmart_dict or 'HiggsBoundsRepo' not in paths_dict or 'HiggsSignalsRepo' not in paths_dict:
        includecodes['HiggsTools'] = False

    if 'HiggsTools' in template_dict['Codes']:
        if includecodes['HiggsTools']:
            outjson['Codes']['HiggsTools']["HiggsBounds Dataset"] =  paths_dict['HiggsBoundsRepo']
            outjson['Codes']['HiggsTools']["HiggsSignals Dataset"] = paths_dict['HiggsSignalsRepo']

            if 'Neutral Higgs' not in template_dict['Codes']['HiggsTools'] or template_dict['Codes']['HiggsTools']['Neutral Higgs'] is None or template_dict['Codes']['HiggsTools']['Neutral Higgs'] == []:
                outjson['Codes']['HiggsTools']["Neutral Higgs"] = bsmart_dict['Neutral Higgs']
                outjson['Codes']['HiggsTools']["Charged Higgs"] = bsmart_dict['Charged Higgs']
        else:
            #outjson['Codes']['HiggsTools']["Run"] = "False"
            outjson['Codes'].pop('HiggsTools')


    if 'MicrOmegas' in template_dict['Codes']:
        if includecodes['MicrOMEGAs']:
            outjson['Codes']['MicrOmegas']['Command'] = MOexec
            
            if "InputFile" not in outjson['Codes']['MicrOmegas'] or outjson['Codes']['MicrOmegas']['InputFile'] is None or outjson['Codes']['MicrOmegas']['InputFile'] == '':
                outjson['Codes']['MicrOmegas']["InputFile"] = "SPheno.spc."+internalmodelname
        else:
            #outjson['Codes']['MicrOmegas']['Run'] = "False"
            outjson['Codes'].pop('MicrOmegas')

    if 'flavio' in template_dict['Codes']:
        outjson['Codes']['flavio']['InputFile']= "WCXF.SMEFT."+internalmodelname

    #if args.VevaciousPlusPlus or args.All:
    if 'VevaciousPlusPlus' in template_dict['Codes']:   
        if includecodes['VevaciousPlusPlus']:
            outjson['Codes']['VevaciousPlusPlus']['Command'] = os.path.join(paths_dict['Vevacious++'],'bin','VevaciousPlusPlus')
            outjson['Codes']['VevaciousPlusPlus']["Initialization File"] = VVinitfile
        else:
            #outjson['Codes']['VevaciousPlusPlus']["Run"] = "False"
            outjson['Codes'].pop('VevaciousPlusPlus')
            if 'Vevacious++' in outjson['Codes']:
                outjson['Codes'].pop('Vevacious++')

    if 'HackAnalysis_LO' in template_dict['Codes']:
        if 'hackanalysis' in paths_dict:
            outjson['Codes']['HackAnalysis_LO']["Command"] = os.path.join(paths_dict['hackanalysis'],'analysePYTHIA.exe')
            if "YAML file" not in outjson['Codes']['HackAnalysis_LO'] or outjson['Codes']['HackAnalysis_LO']['YAML file'] is None or outjson['Codes']['HackAnalysis_LO']['YAML file'] == '':
                outjson['Codes']['HackAnalysis_LO']["YAML file"] = "LOpythia.yaml"
                outjson['Codes']['HackAnalysis_LO']["Config file"]= "LO.cfg"

            
            if "Efficiency Filename" not in outjson['Codes']['HackAnalysis_LO'] or outjson['Codes']['HackAnalysis_LO']['Efficiency Filename'] is None or outjson['Codes']['HackAnalysis_LO']['Efficiency Filename'] == '':    
                outjson['Codes']['HackAnalysis_LO']["Efficiency Filename"] = internalmodelname+".eff"
                outjson['Codes']['HackAnalysis_LO']["Cutflow Filename"] = internalmodelname+".cutflow"
                outjson['Codes']['HackAnalysis_LO']["Histogram Filename"] = internalmodelname+".yoda"
        else:
            #outjson['Codes']['HackAnalysis_LO']["Run"] = "False"
            outjson['Codes'].pop('HackAnalysis_LO')

    if 'MadGraphHackAnalysis' in template_dict['Codes']:
        if 'hackanalysis' in paths_dict and 'madgraph' in paths_dict: 
            outjson['Codes']['MadGraphHackAnalysis']["HackAnalysis"] = os.path.join(paths_dict['hackanalysis'],'analysePYTHIA_LHE.exe')
            outjson['Codes']['MadGraphHackAnalysis']["MadGraph"] =os.path.join(paths_dict['madgraph'],internalmodelname,'bin','generate_events')
            if "YAML file" not in outjson['Codes']['MadGraphHackAnalysis'] or outjson['Codes']['MadGraphHackAnalysis']['YAML file'] is None or outjson['Codes']['MadGraphHackAnalysis']['YAML file'] == '':
                outjson['Codes']['MadGraphHackAnalysis']["YAML file"] = "LOpythia.yaml"
                outjson['Codes']['MadGraphHackAnalysis']["Config file"]= "LO.cfg"

            
            if "Efficiency Filename" not in outjson['Codes']['MadGraphHackAnalysis'] or outjson['Codes']['MadGraphHackAnalysis']['Efficiency Filename'] is None or outjson['Codes']['MadGraphHackAnalysis']['Efficiency Filename'] == '':    
                outjson['Codes']['MadGraphHackAnalysis']["Efficiency Filename"] = internalmodelname+".eff"
                outjson['Codes']['MadGraphHackAnalysis']["Cutflow Filename"] = internalmodelname+".cutflow"
                outjson['Codes']['MadGraphHackAnalysis']["Histogram Filename"] = internalmodelname+".yoda"
        else:
            #outjson['Codes']['MadGraphHackAnalysis']["Run"] = "False"
            outjson['Codes'].pop('MadGraphHackAnalysis')

    if 'ZPEED' in template_dict['Codes']:
        if 'ZPEED' in paths_dict and includecodes['MicrOMEGAs']:
            outjson['Codes']['ZPEED']["Path"] = paths_dict['ZPEED']
        else:
            #outjson['Codes']['ZPEED']["Run"] = "False"
            outjson['Codes'].pop('ZPEED')
            includecodes['ZPEED'] = False

    if 'Zprime' in template_dict['Codes']:
        if not includecodes['MicrOMEGAs']:
            outjson['Codes']['Zprime']['Run'] = "False"
            outjson['Codes'].pop('Zprime')

    if 'SModelS' in template_dict['Codes']:
        if includecodes['SModelS'] and 'SModelS' in paths_dict:
            
            Defaultsmodelsinifile=os.path.join(paths_dict['SModelS'],'parameters.ini')
            templatesmodelsinifile=""
            if 'Parameter file' in template_dict['Codes']['SModelS']:
                templatesmodelsinifile=template_dict['Codes']['SModelS']['Parameter file']
            #else:
            #    templatesmodelsinifile=Defaultsmodelsinifile
            #else:
            #    templatesmodelsinifile=""  # base template contains it but defines it to be ""

            if templatesmodelsinifile != "" and os.path.exists(templatesmodelsinifile):   #  For the case where the parameters.ini is provided as a template! Just need to copy it
                #smodelsinifile=os.path.join(BSMArtScanDir,templatesmodelsinifile)#
                smodelsinifile=os.path.join(BSMArtScanDir,os.path.basename(templatesmodelsinifile))

                #smodelsinifile=os.path.join(BSMArtScanDir,'smodels_parameters.ini')
                outjson['Codes']['SModelS']['Parameter file'] = smodelsinifile   # will always rename it to smodels_parameters.ini
                try:
                    shutil.copy(templatesmodelsinifile,smodelsinifile)
                except: # probably target already exists!
                    pass
            else:
                if templatesmodelsinifile != "": # have named the file, but it's not present
                    smodelsinifile=os.path.join(BSMArtScanDir,os.path.basename(templatesmodelsinifile))
                else:
                    smodelsinifile=os.path.join(BSMArtScanDir,'smodels_parameters.ini')
                outjson['Codes']['SModelS']['Parameter file'] = smodelsinifile
                if os.path.exists(Defaultsmodelsinifile):
                    IF=open(Defaultsmodelsinifile,'r')
                    
                    OFS=open(smodelsinifile,'w')
                    for line in IF:
                        if line.startswith('model=share.models.mssm'):
                            OFS.write('; '+line)
                            OFS.write('model='+os.path.join(BSMArtScanDir,QNUMBERS_slha)+'\n')
                        else:
                            OFS.write(line)
                    OFS.close()
                    IF.close()

            outjson['Codes']['SModelS']['Path'] = paths_dict['SModelS'] 
            outjson['Codes']['SModelS']['Cache'] = paths_dict['SModelS Cache']
            #outjson['Codes']['SModelS']["QNUMBERS file"] = QNUMBERS_slha
            outjson['Codes']['SModelS']["QNUMBERS file"] = os.path.join(BSMArtScanDir,QNUMBERS_slha)
        else:
            #outjson["Codes"]["SModelS"]["Run"] = "False"
            outjson['Codes'].pop('SModelS')

    if 'MultiNest' in paths_dict and 'MultiNestLib' in template_dict['Setup']:
        outjson['Setup']['MultiNestLib']=paths_dict['MultiNest']
    if 'Diver' in paths_dict and 'DiverLib' in template_dict['Setup']:
        outjson['Setup']['DiverLib']=paths_dict['Diver']
    
    
    with open(os.path.join(BSMArtScanDir,template_name),'w') as OF2:
        json.dump(outjson, OF2, cls=CompactJSONEncoder, indent=4)
        #json.dump(outjson, OF2, indent=4)




    if 'BSMArt' in paths_dict:
        bsmartpath=paths_dict['BSMArt']
    else: ## find BSMArt run path from where the script is located
        # bsmartpath = os.path.split(os.path.dirname(os.path.realpath( __file__ )))[0]
        # In package mode, we don't rely on this path for the executable.
        bsmartpath = None

    if template_name[-5:] == '.json':
        runname='_'+template_name[:-5]
    else:
        runname=''
    
    # BSMArtexec=os.path.join(bsmartpath,'bin','BSMArt')
    BSMArtexec = 'BSMArt' # Assume it is in PATH
    if bsmartpath and os.path.exists(os.path.join(bsmartpath, 'bin', 'BSMArt')):
         BSMArtexec=os.path.join(bsmartpath,'bin','BSMArt')

    shellscript=os.path.join(BSMArtScanDir,'runBSMArt'+runname+'.sh')
    OF3=open(shellscript,'w')
    OF3.write('#!/usr/bin/env bash\n\n')
    # OF3.write('python3 '+BSMArtexec + ' --debug '+template_name+' \n')
    if BSMArtexec == 'BSMArt':
         OF3.write(BSMArtexec + ' --debug '+template_name+' \n')
    else:
         OF3.write('python3 '+BSMArtexec + ' --debug '+template_name+' \n')

    OF3.close()
    os.system('chmod u+x '+shellscript)


##########################################################################################################
def makeBSMArtTemplates(modelname,bsmart_dict,paths_dict,includecodes):

    basetemplate=collections.OrderedDict()
    jsonbasefile='BSMArt_basetemplate.json'
    
    # Check current directory
    if os.path.exists(jsonbasefile):
        template_path = jsonbasefile
    else:
        # Check explicit BSMArt path
        if 'BSMArt' in paths_dict and os.path.exists(os.path.join(paths_dict['BSMArt'],'AuxFiles',jsonbasefile)):
            template_path = os.path.join(paths_dict['BSMArt'],'AuxFiles',jsonbasefile)
        else:
            # Check package resources
            try:
                import importlib.resources
                # Assume AuxFiles were moved to bsmart.data.AuxFiles or bsmart.data?
                # I moved AuxFiles to src/bsmart/data. So the file is at bsmart.data.BSMArt_basetemplate.json (if I flat moved it?)
                # Wait, I did `git mv AuxFiles src/bsmart/data/`. So it is `src/bsmart/data/AuxFiles/BSMArt_basetemplate.json`.
                # So resource is `bsmart.data.AuxFiles`
                # Python < 3.9 resources might be tricky with subdirectories.
                # Let's try files() interface.
                from importlib.resources import files
                template_path = files('bsmart.data').joinpath('AuxFiles', jsonbasefile)
                if not template_path.is_file():
                     # Maybe I moved it differently?
                     # I did `git mv AuxFiles src/bsmart/data/` -> `src/bsmart/data/AuxFiles`
                     pass
            except ImportError:
                 # Fallback for older python?
                 pass 
            except Exception as e:
                 pass

    # Actually resolving the file path for updateBSMArtTemplate which expects a path string
    if not os.path.exists(str(template_path)):
         # Last ditch: try finding it relative to this script if user is running from source structure
         # But we are installed.
         raise NameError('Could not find BSMArt_basetemplate.json, essential for constructing the templates')

    # updateBSMArtTemplate(modelname,bsmart_dict,paths_dict,includecodes,jsonbasefile,'BSMArt_'+modelname+'.json')
    # We pass the full path
    updateBSMArtTemplate(modelname,bsmart_dict,paths_dict,includecodes,str(template_path),'BSMArt_'+modelname+'.json')


def main():
    info_message='Script to just set up a BSMArt directory and templates for an existing model. Will not generate any code. Just supply the Model Name.'

    
    parser = argparse.ArgumentParser(
            description=info_message)
    parser.add_argument('ModelName',
                        metavar='<Model Name>', type=str,
                        help='Input Model name')    
    parser.add_argument('--template', type=str,help='Read and modify an example template instead of creating a new one',default=None)
    parser.add_argument('--scriptdir', type=str,help='Directory where the BSMArt data file for the model has been stored. By default, the script assumes the subdirectory ScriptsAndLogs of the current directory.',default=None)
    try:
        args = parser.parse_args()
    except Exception as e:
        print(str(e))
        parser.print_help()
        #print(info_message +str(e))
        raise SystemExit


    ## Set up directories
    cwd = os.getcwd()

    

    if args.scriptdir is not None and os.path.isdir(args.scriptdir):
        scriptdir = args.scriptdir
    else:
        scriptdir=os.path.join(cwd,'ScriptsAndLogs')
        
    bsmartdatafile=os.path.join(scriptdir,'BSMArt_data_'+modelname+'.json') 

    if not os.path.exists(bsmartdatafile):
        from bsmart import get_bsmart_config_dir
        global_config_dir = get_bsmart_config_dir()
        if os.path.isdir(global_config_dir):
            bsmartdatafile=os.path.join(global_config_dir,'BSMArt_data_'+modelname+'.json')
    
    #if not os.path.isdir(scriptdir):
    #    scriptdir=cwd


    model=args.ModelName
    modelname=str(args.ModelName).replace('/','-')
    #scriptname=os.path.join(scriptdir,'init'+modelname+'.m')


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



    includecodes={ 'SPheno':True, 'FlexibleSUSY': True, 'MicrOMEGAs':True, 'Vevacious++':True, 'HiggsSignals':True, 'HiggsBounds':True, 'HiggsTools':True, 'SModelS': True}

    for code in includecodes:
        if code not in paths_dict:
            includecodes[code] = False
   
    includecodes['flavio'] = True
    if 'SARAH' not in paths_dict:
        includecodes['Vevacious++'] = False

    #bsmartdatafile=os.path.join(scriptdir,'BSMArt_data_'+modelname+'.json')     
    try:
        with open(bsmartdatafile) as bsmart_json_data:
            bsmart_dict = json.load(bsmart_json_data, object_pairs_hook=collections.OrderedDict)

        if 'Hneutral' not in bsmart_dict or 'Hcharged' not in bsmart_dict:
            includecodes['HiggsBounds'] = False
            includecodes['HiggsSignals'] = False
        if 'Neutral Higgs' not in bsmart_dict or 'Charged Higgs' not in bsmart_dict:
            includecodes['HiggsTools'] = False
    

    except Exception as e:
        #log.error('Failed to load json file '+file)
        print('Failed to load json file '+bsmartdatafile)
        print('Json exception given as ' +str(e))
        print('Data for SPheno/HiggsTools will not be updated')
        bsmart_dict={} 
        includecodes['SPheno'] = False
        includecodes['HiggsBounds'] = False
        includecodes['HiggsSignals'] = False
        includecodes['HiggsTools'] = False

    if args.template is None:
        makeBSMArtTemplates(modelname,bsmart_dict,paths_dict,includecodes)
    else:
        updateBSMArtTemplate(modelname,bsmart_dict,paths_dict,includecodes,args.template)

if __name__ == "__main__":
    main()
