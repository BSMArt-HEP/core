###############################################################
## Aux functions for madgraph reading
## Part of BSMArt

import os
import fnmatch
import math


import sys
import shutil
from bsmart import debug
from bsmart.HEPRun import DataPoint
from bsmart import zslha
from bsmart.collider import InterpolateXSection

def quadrature(input):
    res=0.0
    for x in input:
        res=res+float(x)**2
    return math.sqrt(res)

def readbanner(rundir):
    infile = None
    for myfile in os.listdir(rundir):
        if fnmatch.fnmatch(myfile, 'run_*_banner.txt'):
            infile=myfile
            break

    if infile is None:
        raise
    IF=open(rundir+'/'+infile,'r')
    xs=None
    for line in IF:
        if 'Integrated weight' in line:
            sline=line.strip()
            data=sline.replace(' ','').split(':')
            xs=float(data[1])
            break
    IF.close()
    if xs:
        return([xs,0.0,0.0])
    else:
        raise

def readmerged(rundir):
    infile=None
    for myfile in os.listdir(rundir):
        if fnmatch.fnmatch(myfile, 'tag_*_merged_xsecs.txt'):
            infile=myfile
            break
    
    if infile is None:
        raise  ## not a genuine error but it found nothing

    IF=open(os.path.join(rundir,infile),'r')
    xs=None
    for line in IF:
        ## take first line that looks like 3 numbers: scale, xsection, uncertainty
        
        sline=line.strip()
        data=(sline.replace('\t',' ')).replace('  ',' ').split()
        if len(data) > 2:
            try:
                xs=float(data[1])
                break
            except:
                continue
    IF.close()
    if xs:
        return([xs,0.0,0.0])
    else:
        raise


def readparton(rundir):
    infile=os.path.join(rundir,'parton_systematics.log')
    if not os.path.isfile(infile):
        res=readbanner(rundir)
        return res

    IF=open(infile,'r')

    for line in IF:
        sline=line.strip()
        rhs=sline.split(':')
        if 'original cross-section' in sline:
            xs=float(rhs[1].replace(' ',''))
        if 'scale variation' in sline:
            scalev=rhs[1].replace(' ','').split('%')
        if 'central scheme variation' in sline:
            central=rhs[1].replace(' ','').split('%')
        if 'PDF variation' in sline:
            pdf=rhs[1].replace(' ','').split('%')
    if not xs or not scalev or not central or not pdf:
        raise
    tplus=quadrature([scalev[0],central[0], pdf[0]])
    tminus=quadrature([scalev[1],central[1], pdf[1]])
    return([xs,tplus*0.01,tminus*0.01])

def readxsdata(rundir):
    """  will try the above routines to see which works """

    try:
        xsdata=readmerged(rundir)
        return xsdata
    except:
        pass

    try:
        xsdata=readparton(rundir)
        return xsdata
    except:
        pass

    try:
        xsdata=readbanner(rundir)
        return xsdata
    except:
        raise NameError('Could not read either merged, parton systematics or banner for cross-section info')
    ## shouldn't be here
    raise 

#class MG_EXEC_LIST():
    """ 
    Class to deal with writing a list of commands to madgraph at execution time
    e.g. bin/generate_events < filename.mg

    Format is as a list of lists. Each entry commands[0], commands[1] is a list of strings
    that will be joined by spaces if they are not null

    e.g. commands[0]=['set decay 1000022','1.96e-15/declen'] or commands[1] = ['set xqcut','mg/6']

    For now, only the second one will be evaluated against the list of variables.
    
    """
#    def __init__(self,commands):
#        self.entries=[]


def write_command_file(fname,commands,var_dict):
    """ function to write a set of commands for madgraph with the possibility of replacing """
    with open(fname, 'w') as f:
        for command in commands:
            outstr=command[0] 
            if len(command) > 1:
                try:
                    outstr=outstr+(' %12.8e' % eval(command[1],var_dict))
                except:
                    outstr=outstr+' ' +str(command[1])
            if len(command) > 2:
                outstr=outstr+' '+str(command[2])
            f.write(outstr+'\n')



"""

Now we define a class to run MadGraph as a subset of another tool!

"""


class MadGraphRunner():
    def __init__(self, settings):

        #if hasattr(global_settings,'cores') and global_settings.cores > 1:
        #    raise NameError("MadGraph tool cannot be run in multicore mode!")

        """
        Initialisation.

        Will give the option to generate the working directory if no command is given and a proc card is given instead!
        
        """
        if 'Process' in settings:
            self.process = settings['Process'] ## to allow labelling the processes
        else:
            self.process='HAPROCESS' 

        if 'python2' in settings:
            self.python2 = eval(settings['python2'])
        else:
            self.python2 = False

        if 'MadGraph' not in settings:
            for setting in ['Proc Card','MadGraph Path']:
                if setting not in settings:
                    raise NameError('Neither Command nor '+setting+' specified for MadGraph; you need to supply the Proc Card and Madgraph Path, or just a command for the working directory!')
            proc_card=settings['Proc Card']
            ## try to read proc card to get the output directory path
            if 'Output Directory' in settings:
                #self.madgraph_cmd=os.path.join(settings['Output Directory'],'bin','generate_events')
                mgdir=settings['Output Directory']

            else:
                mgdir=None
                with open(proc_card,'r') as f:
                    for line in f:
                        if line.startswith('output'):
                            mgdir=line.strip()[7:]
                            break

                if mgdir is None:
                    if 'TempDir' in settings:
                        mgdir=os.path.join(settings['TempDir'],self.process)
                        print("setting madgraph dir: "+mgdir)
                    else:
                        mgdir=self.process
                    with open(proc_card,'a') as f:
                        f.write('output '+mgdir+'\n')
                          
            
            tpath=os.path.dirname(mgdir)
            if tpath == '':
                outdir=os.path.join(os.getcwd(),mgdir)
            else:
                outdir=mgdir
            self.madgraph_cmd=os.path.join(outdir,'bin','generate_events')
            if not os.path.exists(outdir):  
                ## no output directory -> now run the proc card
                try:
                    if self.python2:
                        print('python2 '+settings['MadGraph Path']+' '+proc_card)
                        os.system('python2 '+settings['MadGraph Path']+' '+proc_card)
                    else:
                        os.system(settings['MadGraph Path']+' '+proc_card)

                except:
                    raise NameError('Could not create MadGraph working directory')
                if not os.path.exists(outdir):
                    raise NameError('Could not find or generate Madgraph output directory')
                    
                ## Stop the automatic html opening
                with open(os.path.join(outdir,'Cards','me5_configuration.txt'),'a') as f:
                    f.write('automatic_html_opening=False\n')
                    f.write('auto_update = 0\n')
                    if 'Cores' in settings:
                        f.write('nb_core=%d\n' % int(settings['Cores']))



        else:
            self.madgraph_cmd = settings['MadGraph']



        # Settings
        
        if 'MadGraph Directory' in settings:
            self.MG_Model_Dir = settings['MadGraph Directory']
        else:
            self.MG_Model_Dir = os.path.split(os.path.dirname(self.madgraph_cmd))[0]
        if not os.path.exists(self.MG_Model_Dir):
            raise NameError('Could not find the necessary Madgraph directory')
        self.MG_Event_Dir = os.path.join(self.MG_Model_Dir,'Events')
        self.MG_Run_Dir = os.path.join(self.MG_Event_Dir,'run_01')
        self.MG_hepmc_Dir = self.MG_Run_Dir
        self.MG_HTML_Dir = os.path.join(self.MG_Model_Dir,'HTML','run_01')
        self.MG_Card_Dir = os.path.join(self.MG_Model_Dir,'Cards')
        self.MG_param_card = os.path.join(self.MG_Card_Dir,'param_card.dat')
        self.MG_run_card = os.path.join(self.MG_Card_Dir,'run_card.dat')
        self.MG_pythia_card = os.path.join(self.MG_Card_Dir,'pythia_card.dat')

        
        # self.run_pythia=bool(eval(settings.get('Pythia','False')))
        # if self.run_pythia:
        #     self.cfgfile=settings['Pythia Card']
        #     #self.cfgfile = os.path.join(self.settings['cwd'],self.cfgname)
        #     self.pythia_cfg = PYTHIACFG()
        #     try:
        #         self.pythia_cfg.read(self.cfgfile)
        #     except Exception as e:
        #         raise NameError("Failed to initialise pythia configuration file, " + str(e)) 
            

#        if 'Model_Dir' not in self.settings:
#            log.error('No model directory defined for Madgraph')
#            raise
#        self.MG_Model_Dir = self.settings['Model_Dir']

        if 'QNUMBERS file' in settings and os.path.exists(settings['QNUMBERS file']):
            with open(settings['QNUMBERS file'],'r') as qnfile:
                self.QNUMBERSfile=qnfile.read()
        else:
            self.QNUMBERSfile=None

        """ MG model dir is one step up from the bin directory, use standard trick to get there """
        

        self.run_madspin=False

        if 'MadSpin' in settings:
            self.run_madspin=eval(settings['MadSpin'])


        ## always want the extra commands to make sure we run pythia etc
        ## But we also might want to change settings in the run card such as the xqcut etc, which can be done here

        self.write_mg_extra = True
        self.extra_mg_file = 'extra.mg'
        
       
        mgfname=''
        for cha in self.process:
            if cha.isalnum():
                mgfname += cha
        if len(mgfname) > 0:
            self.extra_mg_file = mgfname+'.mg'
        self.madgraph_cmd = self.madgraph_cmd + ' < '+self.extra_mg_file
        if 'Extra Commands' in settings:
            
            self.simple_mg = True
            
            
            self.extra_commands=list(settings['Extra Commands'])
            for command in self.extra_commands:
                if 'madspin' in command[0]:
                    tempcmd=command[0].replace(' ','')
                    if tempcmd == 'madspin=OFF':
                        self.run_madspin=False
                    else:
                        self.run_madspin=True
                if len(command) > 1:
                    self.simple_mg = False
                    #break
        else:
            self.simple_mg = False
            if self.run_madspin:
                self.extra_commands=[['shower=Pythia8'],['detector=OFF'],['analysis=OFF'],['madspin=ON'],['done'],['done']]
            else:
                self.extra_commands=[['shower=Pythia8'],['detector=OFF'],['analysis=OFF'],['madspin=OFF'],['done'],['done']]


        ## madgraph puts the events in a different directory if it runs madspin
        if self.run_madspin:
            self.MG_hepmc_Dir = self.MG_Run_Dir+'_decayed_1'
            
        """ Allow for a k-factor """

        if 'k-factor' in settings:
            self.kfactor=float(settings['k-factor'])
        else:
            self.kfactor=None


        """ Allow for interpolated cross-sections instead of extracting from MadGraph """
        if 'Cross-sections' in settings:
            try:
                self.XSdata=InterpolateXSection.XS(settings['Cross-sections'])
            except:
                self.XSdata=None

        else:
            self.XSdata=None

    def run(self, spc_file, temp_dir, log,data_point):
        local_file = os.path.join(temp_dir, spc_file)
        #print("local file: "+local_file)
        
        if os.path.exists(os.path.join(self.MG_Model_Dir,'RunWeb')):
            os.remove(os.path.join(self.MG_Model_Dir,'RunWeb'))

        ## new: remove *all* events in the event directory
            
        for filename in os.listdir(self.MG_Event_Dir):
            filepath = os.path.join(self.MG_Event_Dir, filename)
            try:
                shutil.rmtree(filepath)
            except: ## probably a file instead
                try:
                    os.remove(filepath)
                except:
                    raise NameError('Failed to remove run dir '+str(filepath)+'; '+str(e))
                    
        """
        # redundant
        if os.path.exists(self.MG_Run_Dir):
            try:
                shutil.rmtree(self.MG_Run_Dir)
            except Exception as e:
                raise NameError('Failed to remove run dir! '+str(e))
        """
        if os.path.exists(self.MG_HTML_Dir):
            try:
                shutil.rmtree(self.MG_HTML_Dir)
            except Exception as e:
                log.warning('Failed to remove HTML dir! '+str(e))

        if self.run_madspin: ## mg creates an additional directory with _decayed_1 added:
            if os.path.exists(self.MG_hepmc_Dir):
                try:
                    shutil.rmtree(self.MG_hepmc_Dir)
                except Exception as e:
                    raise NameError('Failed to remove run dir! '+str(e))
            if os.path.exists(self.MG_HTML_Dir+'_decayed_1'):
                try:
                    shutil.rmtree(self.MG_HTML_Dir+'_decayed_1')
                except Exception as e:
                    log.warning('Failed to remove HTML dir! '+str(e))

        
        log.debug('Copying '+local_file+' to '+self.MG_param_card)    
        shutil.copyfile(local_file,self.MG_param_card)

        ## append QNUMBERS if defined

        if self.QNUMBERSfile is not None:
            with open(self.MG_param_card,'a') as PC:
                PC.write(self.QNUMBERSfile)

        """  For modifying the pythia card ... not implemented, probably not required because use the command settings
        local_cfg = os.path.join(temp_dir, self.cfgname)
        #self.pythia_cfg.write(local_cfg,data_point.inputs,self.varnames)
        self.pythia_cfg.write_dict(local_cfg,var_dict)
        """

        ## set up madgraph configuration
        
        #var_dict = data_point.get_var_dict(self.varnames)
        if data_point is not None:
            var_dict = data_point.get_all_dict()
        else:
            var_dict={}
        

        if self.write_mg_extra:
            ## if it doesn't exist -> write it 
            ## if it is not simple -> replace it!
            
            #
            #alternative: write it anyway each time! We might have several different MadGraph processes
            # instead we will give different names for different processes
            if not self.simple_mg or not os.path.exists(self.extra_mg_file):
                write_command_file(self.extra_mg_file,self.extra_commands,var_dict)

        """  MG produces too many warnings that fail ...
        debug.command_line_log(self.madgraph_cmd, log)
        """
        
        try:
            if self.python2:
                debug.command_line_log('python2 '+self.madgraph_cmd, log)
            else:
                debug.command_line_log(self.madgraph_cmd, log)
        except: 
            pass # MG produces too many warnings that fail ...
        
        xsdata=None
        if self.run_madspin:
            try:
                xsdata=readxsdata(self.MG_Run_Dir+'_decayed_1')
            except:
                pass
        if xsdata is None:
            try:
                xsdata=readxsdata(self.MG_Run_Dir)
            except:
                log.debug('Failed to read xs data')
                raise

        ### Check whether we want k-factors or imported XSection data

        if self.kfactor is not None:
            xsdata[0]=xsdata[0]*self.kfactor

        if self.XSdata is not None:
            nloxs,nloerr=self.XSdata.get_xs_info(var_dict)
            
            if nloxs > 0:
                newxsdata=[nloxs,nloerr/nloxs,nloerr/nloxs]
            else:
                newxsdata=[nloxs,0.0,0.0]
            xsdata=newxsdata+xsdata
                                 
        


        """ now append to the local file """
        OF = open(local_file,"a")

        OF.write('BLOCK MADGRAPH '+self.process+'\n')
        OF.write(' 1 %10.6E # Cross-section (pb) \n' % xsdata[0])
        OF.write(' 2 %10.4E # + Fractional uncertainty \n' % xsdata[1])
        OF.write(' 3 %10.4E # - Fractional uncertainty \n' % xsdata[2])

        if len(xsdata) >= 6:  ### Have done an interpolation, store originals
            OF.write(' 4 %10.6E # MadGraph computed Cross-section (pb) \n' % xsdata[3])
            OF.write(' 5 %10.4E # + MadGraph computed Fractional uncertainty \n' % xsdata[4])
            OF.write(' 6 %10.4E # - MadGraph computed Fractional uncertainty \n' % xsdata[5])

        OF.close()

        if data_point.spc is None:
            data_point.spc=zslha.read(local_file)
        else:
        #if data_point.spc is not None:
            data_point.spc.blocks['MADGRAPH']={'1': xsdata[0], '2':xsdata[1], '3':xsdata[2]}
            data_point.spc.blockcomments['MADGRAPH']={'1': ' Cross-section (pb)', '2':'+ Fractional uncertainty', '3':'- Fractional uncertainty'}
            if len(xsdata) >= 6:  ### Have done an interpolation, store originals
                data_point.spc.blocks['MADGRAPH'].update({'4': xsdata[3], '5':xsdata[4], '6':xsdata[5]})
                data_point.spc.blockcomments['MADGRAPH'].update({'4': 'MadGraph computed Cross-section (pb)', '5': '+ MadGraph computed Fractional uncertainty', '6':'- MadGraph computed Fractional uncertainty' })
       
        
        hepmc_events=None
        for tfile in os.listdir(self.MG_hepmc_Dir):
            if tfile[-9:] == '.hepmc.gz' or tfile[-6:] == '.hepmc':
                hepmc_events=os.path.join(self.MG_hepmc_Dir,tfile)                    
                break

        
        results={'xs':xsdata[0], '+ Fractional uncertainty':xsdata[1], '- Fractional uncertainty':xsdata[2],'hepmc':hepmc_events}

        return results