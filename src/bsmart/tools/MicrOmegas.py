"""
Tool to run MicrOMEGAs



"""

__meta__ = {
    "name": "MicrOmegas",
    "requires": [],
    "settings": {
        "Command": "Path to MicrOmegas executable",
        "LSP": "List of allowed LSP PIDs",
        "DD_Limits": "Boolean",
        "Cross-sections": "Boolean",
        "OutputFile": "Path to output file",
        "InputFile": "Path to input SLHA file"
    }
}


import os
import shutil
#import subprocess
from bsmart import debug
from bsmart.HEPRun import HepTool, DataPoint
from bsmart.dmlimits import DMlimits
from bsmart import zslha

# def noshell_command_log_withpid(command, wdir, log):
#     log.debug('Terminal command: %s' % command)
#     command_line = subprocess.Popen(command, shell=False, cwd=wdir,
#                                     stdout=subprocess.PIPE,
#                                     stderr=subprocess.PIPE)
#     pid=command_line.pid
#     out, err = command_line.communicate()
#     # log.debug(out)
#     if err != b'':
#         log.error(err.decode("utf-8"))
#     return str(pid)

references=r'''@article{Alguero:2023zol,
    author = "Alguero, G. and Belanger, G. and Boudjema, F. and Chakraborti, S. and Goudelis, A. and Kraml, S. and Mjallal, A. and Pukhov, A.",
    title = "{micrOMEGAs 6.0: N-component dark matter}",
    eprint = "2312.14894",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    doi = "10.1016/j.cpc.2024.109133",
    journal = "Comput. Phys. Commun.",
    volume = "299",
    pages = "109133",
    year = "2024"
}

@article{Alguero:2022inz,
    author = "Alguero, Gael and Belanger, Genevieve and Kraml, Sabine and Pukhov, Alexander",
    title = "{Co-scattering in micrOMEGAs: A case study for the singlet-triplet dark matter model}",
    eprint = "2207.10536",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    doi = "10.21468/SciPostPhys.13.6.124",
    journal = "SciPost Phys.",
    volume = "13",
    pages = "124",
    year = "2022"
}

@article{Belanger:2021smw,
    author = "Belanger, Genevi\`eve and others",
    title = "{Leptoquark manoeuvres in the dark: a simultaneous solution of the dark matter problem and the $ {R}_{D^{\left(\ast \right)}} $ anomalies}",
    eprint = "2111.08027",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    reportNumber = "KIAS-Q21013, MSUHEP-21-031, IFIC/21-46, TTK-21-44, P3H-21-091",
    doi = "10.1007/JHEP02(2022)042",
    journal = "JHEP",
    volume = "02",
    pages = "042",
    year = "2022"
}

@article{Belanger:2020gnr,
    author = "Belanger, Genevieve and Mjallal, Ali and Pukhov, Alexander",
    title = "{Recasting direct detection limits within micrOMEGAs and implication for non-standard Dark Matter scenarios}",
    eprint = "2003.08621",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    doi = "10.1140/epjc/s10052-021-09012-z",
    journal = "Eur. Phys. J. C",
    volume = "81",
    number = "3",
    pages = "239",
    year = "2021"
}

@article{Belanger:2018ccd,
    author = "B\'elanger, Genevi\`eve and Boudjema, Fawzi and Goudelis, Andreas and Pukhov, Alexander and Zaldivar, Bryan",
    title = "{micrOMEGAs5.0 : Freeze-in}",
    eprint = "1801.03509",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    doi = "10.1016/j.cpc.2018.04.027",
    journal = "Comput. Phys. Commun.",
    volume = "231",
    pages = "173--186",
    year = "2018"
}

@article{Barducci:2016pcb,
    author = "Barducci, D. and Belanger, G. and Bernon, J. and Boudjema, F. and Da Silva, J. and Kraml, S. and Laa, U. and Pukhov, A.",
    title = "{Collider limits on new physics within micrOMEGAs$\_$4.3}",
    eprint = "1606.03834",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    doi = "10.1016/j.cpc.2017.08.028",
    journal = "Comput. Phys. Commun.",
    volume = "222",
    pages = "327--338",
    year = "2018"
}

@article{Belanger:2014vza,
    author = "B\'elanger, G. and Boudjema, F. and Pukhov, A. and Semenov, A.",
    title = "{micrOMEGAs4.1: two dark matter candidates}",
    eprint = "1407.6129",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    doi = "10.1016/j.cpc.2015.03.003",
    journal = "Comput. Phys. Commun.",
    volume = "192",
    pages = "322--329",
    year = "2015"
}

@article{Belanger:2013oya,
    author = "Belanger, G. and Boudjema, F. and Pukhov, A. and Semenov, A.",
    title = "{micrOMEGAs$\_$3: A program for calculating dark matter observables}",
    eprint = "1305.0237",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    reportNumber = "LAPTH-023-13",
    doi = "10.1016/j.cpc.2013.10.016",
    journal = "Comput. Phys. Commun.",
    volume = "185",
    pages = "960--985",
    year = "2014"
}

@article{Belanger:2010gh,
    author = "Belanger, G. and Boudjema, F. and Brun, P. and Pukhov, A. and Rosier-Lees, S. and Salati, P. and Semenov, A.",
    title = "{Indirect search for dark matter with micrOMEGAs2.4}",
    eprint = "1004.1092",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    reportNumber = "IRFU-10-24, LAPTH-012-10.",
    doi = "10.1016/j.cpc.2010.11.033",
    journal = "Comput. Phys. Commun.",
    volume = "182",
    pages = "842--856",
    year = "2011"
}
'''



class NewTool(HepTool):
    """ overload the init to initialise the DM limit calculator, but name and settings are already given in HepTool """
    def __init__(self, name, settings,global_settings=None):
        HepTool.__init__(self, name, settings,global_settings)

        self.citations=references
        """ 
            Need to set the hidden envirnoment variable 'nParProc' to 1 or less 
            so that MicrOMEGAs doesn't run many processes at once.
            But we leave it unset if the user specifies single core running, so that
            they can generate the libraries as quickly as possible.
        """
        if (hasattr(global_settings,'cores') and global_settings.cores > 1) or (hasattr(global_settings,'use_MPI') and global_settings.use_MPI):
            os.environ['nParProc'] = '0'

        self.MOworkdir=os.path.join(os.path.dirname(self.settings['Command']),'work')
        self.MOstub='_'+os.uname()[1]+'_'

        self.check_LSP=False
        if 'LSP' in settings:
            self.list_LSP=list(settings['LSP'])
            if len(self.list_LSP) > 0:
                self.check_LSP=True
              
        
        self.ddlimits=False
        if 'DD_Limits' in self.settings:
            if eval(self.settings['DD_Limits']):
                self.ddlimits=True
                self.limitclass=DMlimits()
            else:
                self.ddlimits=False
        else:
            self.ddlimits=False
        self.command=os.path.expanduser(self.settings['Command'])
        if not os.path.isfile(self.command):
            #fdir=os.path.dirname
            raise NameError("Cannot identify MicrOmegas executable. Please give it (without additional command line parameters) under MicrOmegas->settings->Command")
        
        if 'Cross-sections' in self.settings:
            if eval(self.settings['Cross-sections']):
                self.command=[self.command,'--xsections']

        if 'OutputFile' in self.settings:
            self.output_file=self.settings['OutputFile']
        else:
            self.output_file='omg.out'
        
    def run(self, spc_file, temp_dir, log,data_point):
        if self.check_LSP and data_point.spc is not None: ### check that the LSP is in the allowed list, e.g. if it is not charged ... 
            if 'LSP' in data_point.spc.blocks:
                LSP=int(data_point.spc.Value('LSP',[1]))
                if LSP not in self.list_LSP:
                    log.debug('Forbidden LSP found %d! Not running MircOMEGAs' % LSP)
                    raise
                else:
                    log.debug('LSP %d ok, running MicrOMEGAs' % LSP)
                    
        if os.path.exists(self.output_file):
            os.remove(self.output_file) 
        #print(spc_file+', '+self.settings['InputFile'])
        if spc_file != self.settings['InputFile']:
            try:
                shutil.copyfile(spc_file,self.settings['InputFile'])
            except:
                log.debug('File '+spc_file+' not copied')
        # Run MO

        ppid=debug.noshell_command_log_PID(self.command, temp_dir, log)

        #ppid=debug.command_line_log_PID(self.command, temp_dir, log)
        ## do clearup on temporary MO dir 
        
        tdir=os.path.join(self.MOworkdir,self.MOstub+ppid+'_')
        log.debug("Guessing temporary dir: "+tdir)

        if os.path.isdir(tdir):
            log.debug("MO temp directory still exists! Removing it ...")
            shutil.rmtree(tdir)
        
        
        # Glue output to the SPheno file
        if os.path.exists(self.output_file):

            debug.command_line_log('sed -i "1 i\\BLOCK DARKMATTER #" '
                                   + self.output_file, log)

            micrOMEGAS_spc=zslha.read(self.output_file)
            
            if self.ddlimits:
                try:
                    Mdm=float(micrOMEGAS_spc.Value("DARKMATTER",[3]))
                    Oh2=float(micrOMEGAS_spc.Value("DARKMATTER",[1]))
                    SDp=float(micrOMEGAS_spc.Value("DARKMATTER",[202]))
                    SDn=float(micrOMEGAS_spc.Value("DARKMATTER",[204]))
                    SIp=float(micrOMEGAS_spc.Value("DARKMATTER",[201]))
                    SIn=float(micrOMEGAS_spc.Value("DARKMATTER",[203]))
                    SI=(54.0*SIp+(131.0-54.0)*SIn)*1e-36/131.0
                    
                    #print(f"{Mdm}, {Oh2}, {SI}, {SDp}, {SDn}")
                    exclusion=self.limitclass.DMexclusion(Mdm,Oh2,SI,SDp,SDn)
                    micrOMEGAS_spc.blocks['DARKMATTER']['501'] = exclusion
                    micrOMEGAS_spc.blockcomments['DARKMATTER']['501'] = 'BSMArt Direct detection exclusion'
                    strex='%10.6E' % float(exclusion)
                    #debug.command_line_log('echo \"501  '+strex+' # BSMArt Direct detection exclusion \" >> '+ self.output_file, log)
                except Exception as e:
                    log.debug('Exception in runnng ddlimits: '+str(e)) 
                    pass
                #ata_point.spc=micrOMEGAS_spc

            """
            if self.ddlimits:
                try:
                    exclusion=self.limitclass.excludefile(self.output_file)
                    strex='%10.6E' % float(exclusion)
                    debug.command_line_log('echo \"501  '+strex+' # BSMArt Direct detection exclusion \" >> '+ self.output_file, log)
                except Exception as e:
                    log.debug('Exception in runnng ddlimits: '+str(e)) 
                    pass

            ## lazy. We should do this further up
            
            micrOMEGAS_spc=zslha.read(self.output_file)
            if data_point.spc is None:
                data_point.spc=zslha.read(spc_file)
                #ata_point.spc=micrOMEGAS_spc
            """
            #else:
            if data_point.spc is None:
                data_point.spc=zslha.read(spc_file)
            data_point.spc.add(micrOMEGAS_spc)
            


            if self.check_LSP: ### check that the LSP is in the allowed list, e.g. if it is not charged ... 
                LSP=int(micrOMEGAS_spc.Value('DARKMATTER',[2]))
                """
                LSP=0
                with open(self.output_file,'r') as IF:
                    for line in IF:
                        if line.startswith('#'):
                            continue
                        data=line.split()
                        if len(data) < 2:
                            continue
                        if data_point.spc is not None:
                            data_point.blocks["DARKMATTER"][data[0]] = float(data[1])
                        if '# CDM1' in line:
                            #data=line.split()
                            if int(data[0]) == 2:
                                LSP=int(data[1])
                                break
                """
                if LSP not in self.list_LSP:
                    log.debug('Forbidden LSP found %d!' % LSP)
                    raise
                else:
                    log.debug('LSP %d ok' % LSP)
            #debug.command_line_log("echo \"Block DARKMATTER # \" >> "
            #                       + os.path.join(temp_dir, spc_file), log)
            #if os.path.exists(os.path.join(temp_dir, spc_file)):
            #    # append result to it

            # append results to spc file

            
            debug.command_line_log("cat " + self.output_file + " >> "
                                    + os.path.join(temp_dir, spc_file), log)

            
        else:
            log.error("MicrOmegas output not written!")
            raise

        # if self.ddlimits:
        #     try:
        #         exclusion=self.limitclass.excludefile(output_file)
        #         debug.command_line_log('echo \"501  '+str(exclusion)+' # Direct detection exclusion \" >> '+ os.path.join(temp_dir, spc_file), log)
        #     except:
        #         pass
