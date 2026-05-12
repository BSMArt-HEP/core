"""
Tool to run resummino


It will be necessary to provide an input file for resummino

EITHER you define an input file in 'Input File', containing all the process info. That file must contain the input filename and output json file that can then be read.
OR the code will write it itself, potentially for several processes, for which you need to define 'Processes'
 


"""

__meta__ = {
    "name": "resummino",
    "requires": [],
    "settings": {
        "Command": "Path to resummino executable",
        "Input Files": "List of input files",
        "Processes": "List of processes",
        "pdf sets": "List of PDF sets",
        "CMS": "Center of mass energy"
    }
}


import os
import shutil
from bsmart import debug
from bsmart import zslha
from bsmart.HEPRun import HEPRun, RunSettings, CoreRunner, HepTool, DataPoint

import json


def writeresumminofile(fname,process,slhain,jsonout="BSMArt_resummino.json",pdfsets=["NNPDF31_lo_as_0118","NNPDF31_nlo_as_0118"],cms=13000.0):

    OF=open(fname,'w')
    OF.write('collider_type         = proton-proton  # proton-proton or proton-antiproton\n')
    OF.write('center_of_mass_energy = %.0f\n' %cms)
    OF.write('particle1 = %d\n' % int(process[0]))
    OF.write('particle2 = %d\n' % int(process[1]))
    OF.write('result = total  # total, pt, ptj or m.\n')
    OF.write('M      = auto  # auto = sqrt((p1 + p2)^2)\n')
    OF.write('pt     = auto\n')
    OF.write('slha = %s\n' %slhain)
    OF.write('output = %s\n' %jsonout)
    
    OF.write('pdf_format = lhgrid\n')
    OF.write('pdf_lo = %s\n' % pdfsets[0])
    OF.write('pdfset_lo = 0\n')
    OF.write('pdf_nlo = %s\n' % pdfsets[1])
    OF.write('pdfset_nlo = 0\n')

    OF.write('mu_f = 1.0\n')
    OF.write('mu_r = 1.0\n')

    OF.write('precision = 0.01  # desired precision\n')
    OF.write('max_iters = 5  # maximum iterations\n')

    OF.write('Minv_min = auto\n')
    OF.write('Minv_max = auto\n')

    OF.write('xmin = auto\n')

    OF.close()


class NewTool(HepTool):
    def __init__(self, name, settings,global_settings=None):
        
        HepTool.__init__(self, name, settings,global_settings)
        
        self.command=settings['Command']

        self.input_filenames=[]
        self.processes=[]
        self.pdfsets=["NNPDF31_lo_as_0118","NNPDF31_nlo_as_0118"]
        self.jsonout="BSMArt_resummino.json"
        self.write_inputs=True
        self.cms_energy=13000
        self.vary_cms=False
        if 'json output' in settings:
            self.jsonout=settings['json output']

        if 'Processes' in settings:
            for nproc,proc in enumerate(settings['Processes']):
                self.processes.append(list(proc))
                self.input_filenames.append('BSMArt_resummino_%d.in' %nproc)
            if 'pdf sets' in settings:
                self.pdfsets=list(settings['pdf sets'])

            if 'CMS' in settings: ## choose to vary cms energy
                
                try:
                    self.cms_energy=float(settings['CMS'])
                except: # assum it's a variable!
                    self.vary_cms = True
                    self.cms_energy=settings['CMS']    
                           
        elif 'Input Files' in settings:
            self.write_inputs=False
            cwd=os.getcwd()
            for ff in settings('Input Files'):
                
                if os.path.dirname(ff) == "":
                    ff=os.path.join(cwd,ff)                
                if os.path.exists(ff):
                    self.input_filenames.append(ff)
                else:
                    raise NameError("Cannot find input file "+ff)
                
            self.input_filenames=list(settings['Input Files'])
    
        #self.nprocesses=max(len(self.processes),len(self.input_filenames))

    def run(self, spc_file, temp_dir, log,data_point):
        entries={}
        comments={}
        ## relevant if this is the first tool
        if data_point.spc is None:
            data_point.spc= zslha.read(spc_file)
        for nproc,resumfile in enumerate(self.input_filenames):
            if self.write_inputs:
                if self.vary_cms:
                    cms=self.cms_energy
                    if data_point is not None:
                        var_dict = data_point.get_all_dict()
                    else:
                        var_dict={}
                    try:
                        cms=float(eval(self.cms_energy,var_dict))
                    except Exception as e:
                        raise NameError("Failed to evalue CMS energy, " +str(e))
                    writeresumminofile(resumfile,self.processes[nproc],os.path.join(temp_dir,spc_file),self.jsonout, self.pdfsets,cms)
                else:
                    if not os.path.exists(resumfile):
                        writeresumminofile(resumfile,self.processes[nproc],os.path.join(temp_dir,spc_file),self.jsonout, self.pdfsets,self.cms_energy)
            else:
                if not os.path.exists(resumfile):
                    raise NameError("Could not find resummino input file!")

                
            if os.path.exists(self.jsonout):
                os.remove(self.jsonout)

            ## now run resummino!

            debug.command_line_log(self.command+' '+resumfile,log)

            if not os.path.exists(self.jsonout):
                raise NameError("Failed to find output json after running resummino")
            
            ## now read the results:
            with open(self.jsonout,'r') as json_data:
                resummino_info = json.load(json_data)
            
            #res_lo=
            #res_nlo=float(resummino_info['nlo'])
            #res_nll=float(resummino_info['nll'])
            entries[str(nproc*10+1)]=float(resummino_info['lo'])
            entries[str(nproc*10+2)]=float(resummino_info['nlo'])
            entries[str(nproc*10+3)]=float(resummino_info['nll'])
            if self.write_inputs:
                comments[str(nproc*10+1)]='lo for process %d,%d' % (self.processes[nproc][0],self.processes[nproc][1])
                comments[str(nproc*10+2)]='nlo for process %d,%d' %(self.processes[nproc][0],self.processes[nproc][1])
                comments[str(nproc*10+3)]='nll for process %d,%d' % (self.processes[nproc][0],self.processes[nproc][1])
            else:
                comments[str(nproc*10+1)]='lo for process %d' % nproc
                comments[str(nproc*10+2)]='nlo for process %d' %nproc
                comments[str(nproc*10+3)]='nll for process %d' %nproc 
            log.debug('Obtained nll cross-section of %.4e for process %s' %(entries[str(nproc*10+3)],comments[str(nproc*10+3)]))
        data_point.spc.blocks['RESUMMINO']=entries
        data_point.spc.blockcomments['RESUMMINO']=comments

        with open(spc_file,'a') as OF:
            OF.write('BLOCK RESUMMINO \n')
            for entry in entries:
                OF.write(' %s  %10.8e # %s\n' %(entry,entries[entry],comments[entry]))