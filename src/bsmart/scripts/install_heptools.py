#!/usr/bin/env python3

"""
Simple install script for:

SARAH
SPheno
MicrOMEGAs
HiggsBounds
HiggsSignals
HiggsTools
FlexibleSUSY
VevaciousPlusPlus
SModelS
ZPEED
MultiNest
Diver
BSMArt

Everything will be installed as a subdirectory from where the script is called. 
Location info will be stalled in a json file which can later be used by the 

prepareModel.py 

script to set up a specific model

version of 29/01/2025
"""

# Default urls:

default_urls={"SARAH": "https://sarah.hepforge.org/downloads/?f=SARAH-4.15.4.tar.gz", 
              "SPheno": "https://spheno.hepforge.org/downloads/?f=SPheno-4.0.5.tar.gz",
              "HiggsBounds":"https://gitlab.com/higgsbounds/higgsbounds/-/archive/5.10.2/higgsbounds-5.10.2.tar.gz", 
              "HiggsSignals":"https://gitlab.com/higgsbounds/higgssignals/-/archive/2.6.2/higgssignals-2.6.2.tar.gz",
              "MicrOMEGAs":"https://zenodo.org/records/13376690/files/micromegas_6.1.15.tgz?download=1",
              "Vevacious++":"https://github.com/JoseEliel/VevaciousPlusPlus/archive/refs/heads/master.zip",
              "Diver":"https://diver.hepforge.org/downloads/?f=diver-1.0.5.tar.gz", 
              "MultiNest":"https://github.com/farhanferoz/MultiNest/archive/refs/heads/master.zip",
              "HiggsTools": "https://gitlab.com/higgsbounds/higgstools/-/archive/main/higgstools-main.tar.gz",
              "HiggsBoundsRepo":"https://gitlab.com/higgsbounds/hbdataset/-/archive/master/hbdataset-master.tar.gz", 
              "HiggsSignalsRepo":"https://gitlab.com/higgsbounds/hsdataset/-/archive/main/hsdataset-main.tar.gz",
              "SModelS":"https://github.com/SModelS/smodels/archive/refs/heads/main.zip",
              "ZPEED":"https://github.com/kahlhoefer/ZPEED/archive/refs/heads/master.zip",
              "FlexibleSUSY":"https://github.com/FlexibleSUSY/FlexibleSUSY/archive/refs/tags/v2.8.0.tar.gz",
              "COLLIER":"https://collier.hepforge.org/downloads/collier-1.2.8.tar.gz",
              "LoopTools":"https://feynarts.de/looptools/LoopTools-2.16.tar.gz"}


import os
from datetime import datetime
import sys

import shutil
import json
import collections
from collections import OrderedDict
import argparse
import requests

try:
    from bsmart import get_heptools_path_file
except ImportError:
    # Fallback if bsmart is not installed/importable in this context (e.g. running standalone script)
    def get_heptools_path_file(create_dir=False):
        # Default to CWD if we can't determine global config
        return os.path.join(os.getcwd(), 'HEPtoolpaths.json')

import wget
import tarfile
import subprocess
import tempfile

import ssl
ssl._create_default_https_context = ssl._create_unverified_context  

import logging

# Globals
reinstall=False
clean_tarballs=False

cwd = os.getcwd()
tarballdir=os.path.join(cwd,'Downloads')

# Globals populated in main or functions
args = None
log = logging.getLogger('HEPToolsInstall')
latest_urls = default_urls
paths_dict = {}
mathexec = 'math'
installed_something=False
installed_tools=[]

if not os.path.exists(tarballdir):
    os.makedirs(tarballdir)

def shell_command(command_to_execute,working_dir):
    global log


    try:

        command_line = subprocess.Popen(command_to_execute,shell=True, cwd=working_dir,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        out, err = command_line.communicate()
    except Exception as e:
        raise NameError(e)
    res = out.decode().strip()
    #log.debug(out.decode("utf-8").strip())
    
    if err != b'':
        log.debug('Warning messages from '+command_to_execute+': '+err.decode("utf-8").strip())
    return res


def use_requests(temp_url,temp_filename):
    """ This for fastjet, which seems to have blocked the python wget """
    
    try:
        import requests
        
    except Exception as e:
        raise NameError("Failed with requests import " + str(e))

    notrack=False
    try:
        from rich.progress import track
    except:
        notrack=True
    
        
    log.debug('Trying to download ' +temp_url+' by a different method')
    try:
        
        r = requests.get(temp_url,stream=True, headers={'User-agent': 'Mozilla/5.0'})


        r.raise_for_status()
        with open(temp_filename, 'wb') as f:
            total_length = int(r.headers.get('content-length'))
            #print("Length: " +str(total_length))
            #with Bar('Downloading %s' % temp_filename,max=total_length) as bar:
                    #for chunk in progress.bar(r.iter_content(chunk_size=1024), expected_size=(total_length/1024) + 1):
            if notrack:
                chunktotal=(total_length/1024) + 1
                ichunktotal=int(chunktotal)
                #nn = 0
                for nn,chunk in enumerate(r.iter_content(chunk_size=1024)):
                    if chunk:
                        f.write(chunk)
                        f.flush()
                        #nn=nn+1
                        done=nn*50//ichunktotal
                        todo=50-done
                        percentage=nn*100//ichunktotal
                    
                        string='['+('='*done)+(' '*todo)+'] Downloading '+str(percentage)+'% '
                        print(string,end='\r')
                print('['+('='*50)+'] Downloading 100%')
            else:
                for chunk in track(r.iter_content(chunk_size=1024),total=(total_length/1024) + 1,description='Downloading '+temp_filename):
                    if chunk:
                        f.write(chunk)
                        f.flush()
                        
            
            return temp_filename
    except Exception as e:
        raise NameError("Failed to download using requests, "+str(e))


def get_file(temp_url,user_filename=None):
    global tarballdir,cwd
    #temp_filename=os.path.basename(temp_url).strip('?f=')
    if user_filename is None:
        temp_filename=os.path.basename(temp_url).replace('?f=','')
    else:
        temp_filename=user_filename
    log.info("Getting %s " % temp_filename)
    
    os.chdir(tarballdir)
            
    if os.path.exists(temp_filename) and not clean_tarballs:
        log.debug('Found '+temp_filename)
        os.chdir(cwd)
        return temp_filename
    else:
        if clean_tarballs and os.path.exists(temp_filename):
            os.remove(temp_filename)
            
        try:
            filename = wget.download(temp_url,temp_filename)
            os.chdir(cwd)
            return filename
        except Exception as e:
            os.chdir(cwd)
            log.info("Could not download " +str(temp_url) +", " +str(e))
            
    try:
        os.chdir(tarballdir)
        filename = use_requests(temp_url,temp_filename)
        os.chdir(cwd)
        return filename
    except Exception as ee:
        os.chdir(cwd)
        raise NameError(str(ee))


def extract_tar(tarball):
    """ extract the tarball and return the base directory """
    global reinstall,tarballdir
    tar = tarfile.open(os.path.join(tarballdir,tarball))
    #firstentry=tar.next()
    ## This to ignore stupid MacOS files found in MicrOMEGAs release 
    allentries=[ x for x in tar.getnames() if not x.startswith('._')]
    
    #allentries=tar.list()
    baseprefix=os.path.commonprefix(allentries)
    if os.path.exists(baseprefix):
        if reinstall:
            shutil.rmtree(baseprefix)
        else:
            return baseprefix,True
            #raise NameError(baseprefix)
            
#    print(baseprefix)
#    print ("my name = %s "%tarball.strip('.tar.gz'))
    if baseprefix == '':
        baseprefix = tarball.strip('.tar.gz')
        tar.extractall(path=baseprefix)
        tar.close()
    else:
        tar.extractall()
        tar.close()
    return baseprefix,False
    

def install_SARAH():
    global latest_urls
    sarahurl=latest_urls["SARAH"]
    tfile=get_file(sarahurl)
    try:
        basedir,pre_existing=extract_tar(tfile)
    except Exception as e:
        #basedir=str(e)
        raise
    return basedir

"""
def install_BSMArt():
    global latest_urls
    tempurl=latest_urls["BSMArt"]
    tfile=get_file(tempurl)
    try:
        basedir,pre_existing=extract_tar(tfile)
    except Exception as e:
        #basedir=str(e)
        raise
    return basedir

"""

def install_ZPEED():
    global latest_urls,reinstall,tarballdir,cwd
    baseprefix='ZPEED-master'
    if os.path.exists(baseprefix):
        if reinstall:
            shutil.rmtree(baseprefix)
        else:
            return baseprefix
            #raise NameError(baseprefix)
    tempurl=latest_urls["ZPEED"]
    #os.chdir(tarballdir)
    zipname='zpeed-main.zip'
    try:
        zipfile=get_file(tempurl,zipname)
    except Exception as e:
        raise NameError("Failed to find or download SModelS, "+str(e))
    
    ## should be a zip file
    os.chdir(cwd)
    # now try to extract
    log.debug('Exctracting ZPEED zip')
    #_=shell_command("unzip "+os.path.join(tarballdir,zipname),cwd)
    _=shell_command("unzip "+os.path.join(tarballdir,zipfile),cwd)

    return baseprefix 


def install_SPheno(F90):
    global latest_urls
    sphenourl=latest_urls["SPheno"]
    tfile=get_file(sphenourl)
    try:
        basedir,pre_existing=extract_tar(tfile)
        if pre_existing:
            return basedir
    except Exception as e:
        raise NameError("Error extracting SPheno "+str(e))
    
        
    makepath=os.path.join(cwd,basedir,"Makefile")
    ## modify the Makefile to select the compiler
    IF=open(makepath,'r')
    OF=open("temp",'w')
    OF.write("F90 = %s\n" % F90)
    for line in IF:
        if "F90" in line:
            #print(line.strip()[:2])
            if line.strip().startswith('F90'):
                OF.write('# '+ line)
            else:
                OF.write(line)
        else:
            OF.write(line)
    IF.close()
    OF.close()

    shutil.move("temp",makepath)
    ## Compile SPheno
    #os.chdir(os.path.join(cwd,basedir))
    _ = shell_command("make",os.path.join(cwd,basedir))

    os.chdir(cwd)
    
    return basedir

def install_MicrOMEGAs():
    global latest_urls
    tempurl=latest_urls["MicrOMEGAs"]
    tfile=get_file(tempurl)
    try:
        basedir,pre_existing=extract_tar(tfile)
        if pre_existing:
            return basedir
    except Exception as e:
        
        raise NameError("Error extracting MicrOMEGAs "+str(e))

    ## Compile

    log.info("Building MicrOMEGAs")
    _ = shell_command("make",os.path.join(cwd,basedir))

    os.chdir(cwd)
    
    return basedir


def install_LoopTools():
    global latest_urls
    lturl=default_urls["LoopTools"]
    tfile=get_file(lturl)
    try:
        basedir,pre_existing=extract_tar(tfile)
        if pre_existing and not reinstall:
            return os.path.join(cwd, basedir) ## installation directory
        
        basedir_src = basedir + '_src'
        shutil.move(basedir, basedir_src)
    except Exception as e:
        raise NameError("Error extracting tarball for " +str(e))
    installdir=os.path.join(cwd, basedir)

    
    
    builddir=os.path.join(cwd,basedir_src)

    log.debug("Building LoopTools")
    #installdir=os.path.join(cwd, basedir)
    _=shell_command('FFLAGS="-O3 -fPIC" ./configure --prefix=%s' %installdir, builddir)
    _=shell_command("make",builddir)
    _=shell_command("make install", builddir)

    os.chdir(cwd)
    shutil.rmtree(basedir_src)

    return installdir

def install_COLLIER():
    global latest_urls
    clurl=default_urls["COLLIER"]
    tfile=get_file(clurl)
    try:
        basedir,pre_existing=extract_tar(tfile)
        if pre_existing and not reinstall:
            return os.path.join(cwd, basedir) ## installation directory
        basedir_src = basedir + '_src'
        shutil.move(basedir, basedir_src)
    except Exception as e:
        raise NameError("Error extracting tarball for " +str(e))

    installdir=os.path.join(cwd,basedir)
    
    builddir=os.path.join(cwd,basedir_src,'build')
    if os.path.exists(builddir):
        shutil.rmtree(builddir)
    os.makedirs(builddir)

    log.debug("Building COLLIER")
    #installdir=os.path.join(cwd,basedir)
    _=shell_command('cmake .. -Dstatic=ON -DCMAKE_POSITION_INDEPENDENT_CODE=ON -DCMAKE_INSTALL_PREFIX=%s' %installdir, builddir)
    _=shell_command("make",builddir)
    _=shell_command("make install",builddir)

    os.chdir(cwd)
    shutil.rmtree(basedir_src)

    return installdir

def install_FlexibleSUSY():
    global latest_urls
    fsurl=default_urls["FlexibleSUSY"]
    tfile=get_file(fsurl)
    try:
        basedir,pre_existing=extract_tar(tfile)
    except Exception as e:
        raise NameError("Error extracting tarball for" +str(e))

    return basedir

def install_HiggsBoundsSignals(NAME):
    global latest_urls
    tempurl=latest_urls[NAME]
    tfile=get_file(tempurl)
    try:
        basedir,pre_existing=extract_tar(tfile)
        if pre_existing:
            builddir=os.path.join(basedir,'build')
            if os.path.exists(builddir):
                return builddir
            
    except Exception as e:
        raise NameError("Error extracting tarball for "+NAME+", "+str(e))

    ## Make build directory
        
    ## Compile HiggsBounds
    #os.chdir(os.path.join(cwd,basedir))
    
    builddir=os.path.join(cwd,basedir,'build')
    #print("now to build " +builddir)
    if os.path.exists(builddir):
        shutil.rmtree(builddir)
    os.makedirs(builddir)

    log.debug("Building "+NAME)
    _=shell_command("cmake ..",builddir)
    _=shell_command("make",builddir)
    
    os.chdir(cwd)
    
    return builddir


def install_Vevaciouspp():
    global latest_urls,reinstall,tarballdir,cwd
    baseprefix='VevaciousPlusPlus-master'
    if os.path.exists(baseprefix):
        if reinstall:
            shutil.rmtree(baseprefix)
        else:
            return baseprefix
            #raise NameError(baseprefix)

    tempurl=latest_urls["Vevacious++"]
    
    
    #os.chdir(tarballdir)
    zipname='Vevacious.zip'
    try:
        zipfile=get_file(tempurl,zipname)
    except Exception as e:
        raise NameError("Failed to find or download Vevacious, "+str(e))
    """
    if os.path.exists(zipname):
        if reinstall:
            os.remove(zipname)
    else:

        
        try:
            wget.download(tempurl,zipname)
        except Exception as e:
            raise NameError("Could not download Vevacious++, "+str(e))
    """
    ## should be a zip file
    os.chdir(cwd)
    # now try to extract
    log.debug('Exctracting Vevacious++ zip')
    #_=shell_command("unzip "+os.path.join(tarballdir,zipname),cwd)
    _=shell_command("unzip "+os.path.join(tarballdir,zipfile),cwd)
    
    
    #print("Now to create build dir")
    builddir=os.path.join(cwd,baseprefix,'build')
    #print("now to build " +builddir)
    if os.path.exists(builddir):
        shutil.rmtree(builddir)
    os.makedirs(builddir)
    log.debug('Building Vevacious++')
    _ = shell_command("cmake ..",builddir)
    _ = shell_command("make", builddir)
    
    os.chdir(cwd)
    
    return baseprefix 


def install_Diver(FF):
    global latest_urls
    diverurl=latest_urls["Diver"]
    tfile=get_file(diverurl)
    try:
        basedir,pre_existing=extract_tar(tfile)
        if pre_existing:
            return basedir
    except Exception as e:
        raise NameError("Error extracting SPheno "+str(e))
    
        
    makepath=os.path.join(cwd,basedir,"makefile")

    ## Need to get gcc verion:
    fortver=-1
    
    if FF=='gfortran' or FF=='mpif90':
        tver=shell_command(FF+' -dumpversion',os.path.join(cwd,basedir))
        try:
            ver=int(tver)
            fortver=ver
        except:
            pass
    
    ## modify the Makefile to select the compiler
    print('WARNING: I am using the default mpif90 compiler and the options that allow diver to compile on my machine.\nYou may need to modify them to get Diver to compile')
    IF=open(makepath,'r')
    OF=open("temp",'w')
    #OF.write("F90 = %s\n" % F90)
    for line in IF:
        sline=line.strip()
        if sline.startswith('#'):
            OF.write(line)
        elif sline.startswith('FOPT='):
            if fortver > 9:
                OF.write('# '+ line)
                OF.write('FOPT=-g -ffixed-line-length-none -DMPI -cpp -fallow-argument-mismatch -Wno-conversion #-Wall -fcheck=all\n')
            else:
                OF.write(line)
        elif sline.startswith('FF='):
            OF.write('# '+ line)
            OF.write('FF='+str(FF)+'\n')
        else:
            OF.write(line)
    IF.close()
    OF.close()

    shutil.move("temp",makepath)
    log.debug("Making Diver")
    _ = shell_command("make libdiver.so",os.path.join(cwd,basedir))

    os.chdir(cwd)
    
    return basedir

    
def install_MultiNest():
    global latest_urls,reinstall,cwd
    baseprefix='MultiNest-master'
    os.chdir(cwd)
    if os.path.exists(baseprefix):
        if reinstall:
            shutil.rmtree(baseprefix)
        else:
            dirs=os.listdir(os.path.join(cwd,baseprefix))
            cmakedir=''
            for tdir in dirs:
                if '_CMake' in tdir:
                    cmakedir=str(tdir)
                    break
            if cmakedir == '':
                cmakedir= 'MultiNest_v3.12_CMake'

            libpath=os.path.join(cwd,baseprefix,cmakedir,'multinest','lib') 
            if os.path.exists(libpath):
                return libpath
                ## otherwise we need to carry on
            #raise NameError(baseprefix)

    tempurl=latest_urls["MultiNest"]
    #print(tempurl)
    zipname='MultiNest.zip'
    try:
        zipfile=get_file(tempurl,zipname)
    except Exception as e:
        raise NameError("Failed to find or download MultiNest, "+str(e))
    """
    if os.path.exists(zipname):
        if reinstall:
            os.remove(zipname)
    else:
        try:
            wget.download(tempurl,zipname)
        except Exception as e:
            raise NameError("Could not download MultiNest, "+str(e))
    
    ## should be a zip file
    
    # now try to extract
    log.debug("Unzipping multinest")
    _ = shell_command("unzip "+zipname,cwd)
    """
    # now try to extract
    log.debug("Unzipping multinest")
    _ = shell_command("unzip "+os.path.join(tarballdir,zipfile),cwd)

    #### When I write this, the version is MultiNest_v3.12_CMake
    #### Hopefully if this gets updated, then the format will stay the same
    dirs=os.listdir(os.path.join(cwd,baseprefix))
    cmakedir=''
    for tdir in dirs:
        if '_CMake' in tdir:
            cmakedir=str(tdir)
            break
    if cmakedir == '':
        cmakedir= 'MultiNest_v3.12_CMake'


    builddir=os.path.join(cwd,baseprefix,cmakedir,'multinest','build')
    #print("now to build " +builddir)
    if os.path.exists(builddir):
        shutil.rmtree(builddir)
    os.makedirs(builddir)

    log.info("Building MultiNest")
    _ = shell_command("cmake ..",builddir)
    _ = shell_command("make",builddir)
    
    os.chdir(cwd)
    
    return os.path.join(cwd,baseprefix,cmakedir,'multinest','lib') 


def install_SModelS():
    global latest_urls,reinstall,tarballdir,cwd
    baseprefix='smodels-main'
    if os.path.exists(baseprefix):
        if reinstall:
            shutil.rmtree(baseprefix)
        else:
            return baseprefix
            #raise NameError(baseprefix)

    tempurl=latest_urls["SModelS"]
    
    
    #os.chdir(tarballdir)
    zipname='smodels-main.zip'
    try:
        zipfile=get_file(tempurl,zipname)
    except Exception as e:
        raise NameError("Failed to find or download SModelS, "+str(e))
    
    ## should be a zip file
    os.chdir(cwd)
    # now try to extract
    log.debug('Exctracting SModelS zip')
    #_=shell_command("unzip "+os.path.join(tarballdir,zipname),cwd)
    _=shell_command("unzip "+os.path.join(tarballdir,zipfile),cwd)
    
    ## make sure requirements are satisfied

    _=shell_command("pip3 install -r "+os.path.join(cwd,'smodels-main','smodels','share','requirements.txt'),cwd)
    
    return baseprefix 


def add_SARAH_paths(fakedir=None):
    """
    Add the path for the SARAH installation to the Mathematica directory.
    If fakedir is given, we create a local environment 
    """

    global paths_dict,mathexec
    
    if 'SARAH' not in paths_dict:
        raise NameError('SARAH not already installed!')

    """
    tempfile defaults to w+b so that it can be read/write without being closed. But then we can only write
    bytestrings to it. Instead we can use mode='w' so that I can just write to the file
    """
    
    
    #with tempfile.NamedTemporaryFile(suffix='.m',delete_on_close=False,mode='w') as FO:
    with tempfile.NamedTemporaryFile(suffix='.m',mode='w',delete=False) as FO:  # delete_on_close is new! Damnit
        if fakedir is not None: # check/create a fake directory to kid mathematica into
            if not os.path.exists(fakedir):
                raise NameError('Provided mathematica launch directory does not exist: ' +str(fakedir))
            userbase=None
            userbaseline=r'Print[$UserBaseDirectory];'+'\n'
            #userbaseline=b'Print[$UserBaseDirectory];\n'
            FO.write(userbaseline)
            FO.close()
            command=mathexec +' -script '+FO.name
            print("Running command: "+command)
            try:
                #os.system('cat '+FO.name)
                #os.system(command)
                userbase = shell_command(command,cwd)
                print("Found %s (userbase) " %userbase)
                
                #mathlicensedir=os.path.join(userbase,".Mathematica","Licensing")
                mathlicensedir=os.path.join(userbase,"Licensing")
                if os.path.exists(mathlicensedir):
                    os.system('ln -s %s %s' %(mathlicensedir,fakedir))
                if not os.path.exists(os.path.join(fakedir,'Kernel')):
                    os.mkdir(os.path.join(fakedir,'Kernel'))
                with open(os.path.join(fakedir,'Kernel','init.m'),'w') as OF:
                    OF.write(r'PrependTo[$Path,"'+paths_dict['SARAH']+'"];\n')
                paths_dict['MathematicaUserBaseDirectory']=fakedir
                log.info("Added path to Mathematica!")
            except Exception as e:
                log.warn("Failed to get Mathematica to update its own path!" +str(e))         
                
        else:   # Just add path to default mathematica kernel
            #FO, tempscript=tempfile.mkstemp(suffix='.m')
            cmdline='With[{newPath = {"'+paths_dict['SARAH']+'"}},PutAppend[Unevaluated[$Path = Join[newPath,$Path]],FileNameJoin[{$UserBaseDirectory, "Kernel", "init.m"}]]]\n'
            FO.write(cmdline)
            
            command=mathexec+' -script '+FO.name
            try:
                _ = shell_command(command,cwd)
                log.info("Added path to Mathematica!")
            except Exception as e:
                log.warn("Failed to get Mathematica to update its own path!" +str(e))
        # clean up file
        os.remove(FO.name)  # so frustrating
    
        


def go_SARAH():
    global installed_something,mathexec
    try:
        sarah_path=install_SARAH()
    except Exception as e:
        log.critical('Failed to install SARAH, '+str(e))
        raise NameError ('Failed to install SARAH, '+str(e))
    #print("SARAH path: %s " %sarah_path)
    paths_dict['SARAH']=os.path.join(cwd,sarah_path)
    installed_something=True
    if args.AddPaths:
        ## First need to get the Mathematica $UserBaseDirectory
        try:
            add_SARAH_paths()
        except Exception as e:
            log.critical('Failed to update paths as requested, '+str(e))
            raise
    else:
        log.info("Installed SARAH!")
        print("Consider adding:")
        print('$Path = Join[{"'+paths_dict['SARAH']+'"}, $Path]\nto your Mathematica Kernel/init.m')
        log.debug("Installed SARAH!\n Consider adding:")
        log.debug('$Path = Join[{"'+paths_dict['SARAH']+'"}, $Path]\nto your Mathematica Kernel/init.m')
        installed_tools.append('SARAH')
    

def go_SPheno():
    global installed_something
    log.info('Installing SPheno')
    try:
        SPheno_path=install_SPheno(args.Fortran)
    except Exception as e:
        raise NameError('Failed to install SPheno, '+str(e))
    paths_dict['SPheno']=os.path.join(cwd,SPheno_path)
    installed_something=True
    log.info("Installed SPheno")
    
    installed_tools.append('SPheno')
    

def go_MicrOMEGAs():
    global installed_something
    try:
        MO_path=install_MicrOMEGAs()
    except Exception as e:
        log.error('Problem installing MicrOMEGAs: ' +str(e))
        raise 
    
    paths_dict['MicrOMEGAs']=os.path.join(cwd,MO_path)
    installed_something=True
    log.info("Installed MicrOMEGAs")
    installed_tools.append('MicrOMEGAs')


def go_HiggsBounds():
    global installed_something
    try:
        HBbuild_path=install_HiggsBoundsSignals('HiggsBounds')
    except Exception as e:
        print('Problem installing HiggsBounds: ' +str(e))
        raise 

    paths_dict['HiggsBounds']=os.path.join(cwd,HBbuild_path)
    installed_something=True
    print("Installed HiggsBounds")
    installed_tools.append('HiggsBounds')


def go_HiggsSignals():
    global installed_something
    log.info('Installing HiggsSignals')
    try:
        HSbuild_path=install_HiggsBoundsSignals('HiggsSignals')
    except Exception as e:
        log.error('Problem installing HiggsSignals: ' +str(e))
        raise 
    paths_dict['HiggsSignals']=os.path.join(cwd,HSbuild_path)
    installed_something=True
    log.info("Installed HiggsSignals")
    installed_tools.append('HiggsSignals')


def go_HiggsTools():
    global installed_something
    try:
        HTbuild_path=install_HiggsBoundsSignals('HiggsTools')
        HTmainpath=os.path.split(os.path.join(cwd,HTbuild_path))[0]        
    except Exception as e:
        log.error('Problem installing HiggsTools: ' +str(e))
        raise SystemExit
    ### Add the path to the tool:

    paths_dict['HiggsTools']=os.path.join(cwd,HTbuild_path)

    
    
    
    
    """ now install repos """
    tempurl=latest_urls['HiggsBoundsRepo']
    tfile=get_file(tempurl)
    try:
        basedir,pre_existing=extract_tar(tfile)
        paths_dict['HiggsBoundsRepo'] = os.path.join(cwd,basedir)
    except:
        log.error('Problem downloading HiggsBounds repository')
        raise 
    tempurl=latest_urls['HiggsSignalsRepo']
    tfile=get_file(tempurl)
    try:
        basedir,pre_existing=extract_tar(tfile)
        paths_dict['HiggsSignalsRepo'] = os.path.join(cwd,basedir)
    except:
        log.error('Problem downloading HiggsSignals repository')
        raise 
    installed_something=True
    log.info('Installing python interface to HiggsTools')
    if reinstall:
        try:
            import Higgs
            shell_command('pip3 uninstall Higgs',HTmainpath)

        except: ## higgs not installed anyway
            pass
    try: 
        import Higgs
        #import pylha
    except:

        try:
            shell_command('pip3 install .',HTmainpath)
            #shell_command('pip3 install pylha',HTmainpath)
        except Exception as e:
            log.error('Failed to install HiggsTools python interface '+str(e))
            raise 

    
    log.info("Installed HiggsTools!")
    installed_tools.append('HiggsTools')


def go_VevaciousPlusPlus():
    log.info('Installing Vevacious++')
    try:
        VVbuild_path=install_Vevaciouspp()
        paths_dict['Vevacious++']=os.path.join(cwd,VVbuild_path)
        installed_something=True
        log.info("Installed Vevacious++")
        installed_tools.append('Vevacious++')
    except Exception as e:
        log.error('Problem installing Vevacious++: ' +str(e))
        raise


def go_SModelS():
    log.info('Installing SModelS')
    try:
        smodels_path=install_SModelS()
        paths_dict['SModelS']=os.path.join(cwd,smodels_path)
        paths_dict['SModelS Cache']=os.path.join(cwd,smodels_path,'CACHE')
        #if os.path.exists(paths_dict['SModelS Cache']):
        #    shutil.rmtree(paths_dict['SModelS Cache'])
        if not os.path.exists(paths_dict['SModelS Cache']):
            os.makedirs(paths_dict['SModelS Cache'])
        installed_something=True
        log.info("Downloaded SModelS")
        installed_tools.append('SModelS')
    except Exception as e:
        log.error('Problem Downloading SModelS: ' +str(e))
        raise


def go_ZPEED():
    log.info('Installing ZPEED')
    try:
        zpeed_path=install_ZPEED()
        paths_dict['ZPEED']=os.path.join(cwd,zpeed_path)

        installed_something=True
        log.info("Downloaded ZPEED")
        installed_tools.append('ZPEED')
    except Exception as e:
        log.error('Problem Downloading ZPEED: ' +str(e))
        raise


def go_MultiNest():
    log.info('Installing MultiNest')
    try:
        MultiNest_build_path=install_MultiNest()
        paths_dict['MultiNest']=MultiNest_build_path
        installed_something=True
        log.info("Installed MultiNest")
        installed_tools.append('MultiNest')
    except Exception as e:
        log.error('Problem installing MultiNest: ' +str(e))
        raise

def go_Diver():
    try:
        Diver_path=install_Diver(args.FF)
        paths_dict['Diver']=os.path.join(cwd,Diver_path,'lib')
        installed_something=True
        log.info("Installed Diver") 
        installed_tools.append('Diver')
    except Exception as e:
        log.error('Failed to install Diver! You may need to check the settings. '+str(e))
        raise

def go_COLLIER():
    try:
        fs_path=install_COLLIER()
        installed_something=True
        log.info("Installed COLLIER")
        installed_tools.append('COLLIER')
        paths_dict['COLLIER'] = fs_path
    except Exception as e:
        log.error('Failed to install COLLIER! You may need to check the settings. '+str(e))
        raise

def go_LoopTools():
    try:
        fs_path=install_LoopTools()
        installed_something=True
        paths_dict['LoopTools'] = fs_path
        log.info("Installed LoopTools")
        installed_tools.append('LoopTools')
        return fs_path
    except Exception as e:
        log.error('Failed to install LoopTools! You may need to check the settings. '+str(e))
        raise

def go_FlexibleSUSY():
    try:
        fs_path=install_FlexibleSUSY()
        fs_path=os.path.join(cwd,fs_path)
        paths_dict['FlexibleSUSY'] = fs_path
        installed_something=True
        log.info("Installed FlexibleSUSY") 
        installed_tools.append('FlexibleSUSY')
        
        # Now to set up paths so that FlexibleSUSY can find SARAH via Needs`
        if not args.AddPaths:
            add_SARAH_paths(fs_path)
        log.info('Created special directories')
        

    except Exception as e:
        log.error('Failed to install FlexibleSUSY! You may need to check the settings. '+str(e))
        raise

"""
def go_BSMArt():
    log.info('Installing BSMArt')
    try:
        build_path=install_BSMArt()
        paths_dict['BSMArt']=os.path.join(cwd,build_path)
        installed_something=True
        log.info("Installed BSMArt!")
        log.debug("Consider adding:\n"+os.path.join(paths_dict['BSMArt'],'bin')+"\n to your environment path, or creating a link to BSMArt.")
        print("Consider adding:\n"+os.path.join(paths_dict['BSMArt'],'bin')+"\n to your environment path, or creating a link to BSMArt.")
        installed_tools.append('BSMArt')
    except Exception as e:
        log.error('Failed to install BSMArt! '+str(e))
        raise 
"""

def main():
    global args, reinstall, clean_tarballs, mathexec, log, latest_urls, paths_dict, installed_something, installed_tools

    info_message='Please specify the tool(s) you want to install.\n Valid options are: All, SARAH, SPheno, MicrOmegas, HiggsBounds, HiggsSignals, VevaciousPlusPlus, SModelS, HiggsTools, ZPEED, FlexibleSUSY'

    try:
        parser = argparse.ArgumentParser(
            description=info_message)
        parser.add_argument("--All", help="Install all tools",
                            action="store_true")
        parser.add_argument("--SARAH",  help="Install SARAH",
                            action="store_true")
        parser.add_argument("--SPheno", help="Install SPheno",
                            action="store_true")
        parser.add_argument("--HiggsBounds", help="Install HiggsBounds",
                            action="store_true")
        parser.add_argument("--HiggsSignals", help="Install HiggsSignals",
                            action="store_true")
        parser.add_argument("--HiggsTools", help="Install HiggsTools",
                            action="store_true")
    #    parser.add_argument("--git", help="Download using git rather than zip", action="store_true")
        parser.add_argument("--MicrOMEGAs", help="Install MicrOMEGAs",
                            action="store_true")
        parser.add_argument("--VevaciousPlusPlus", help="Install Vevacious++",
                                action="store_true")
        parser.add_argument("--ZPEED", help="Install ZPEED",
                                action="store_true")
        parser.add_argument("--MultiNest", help="Install MultiNest",
                            action="store_true")
        parser.add_argument("--Diver", help="Install Diver",
                            action="store_true")         
        parser.add_argument("--SModelS", help="Install SModelS",
                            action="store_true")
        parser.add_argument("--FlexibleSUSY", help="Install FlexibleSUSY",
                            action="store_true")                    
        parser.add_argument("--COLLIER", help="Install COLLIER",
                            action="store_true")           
        #parser.add_argument("--BSMArt", help="Install BSMArt",
        #                    action="store_true")
        parser.add_argument("--math", help="Mathematica executable name (to run on command line)",nargs="?",action="store",
                            default='math')
        parser.add_argument("--Fortran", help="Fortran compiler name",nargs="+",action="store",
                            default='gfortran')
        parser.add_argument("--FF", help="mpi Fortran compiler name for Diver, e.g. mpif90",nargs="+",action="store",
                            default='mpif90')                    
        parser.add_argument("--Reinstall", help="Reinstall packages",
                            action="store_true")
        parser.add_argument("--Redownload",help="Delete local copies of tarballs before reinstallation",action="store_true")
        parser.add_argument("--AddPaths", help="Configure paths for SARAH and/or BSMArt",
                            action="store_true")
        parser.add_argument("--LocalUrls", help="Don't try to update the urls from the internet",
                            action="store_true")
        parser.add_argument("--local", help="Install HEPtoolpaths.json locally in CWD",
                        action="store_true")

        args = parser.parse_args()
    except Exception as e:
        print(info_message +str(e))
        raise SystemExit

    if args.Reinstall:
        reinstall=True

    if args.Redownload:
        reinstall=True
        clean_tarballs=True

    mathexec=args.math

    #print(cwd)


    scriptdir=os.path.join(cwd,'ScriptsAndLogs')

    if not os.path.exists(scriptdir):
        os.makedirs(scriptdir)

    """ Set up logging """
    # log is already defined as global, we configure it here
    logfile=os.path.join(scriptdir,'HEPToolsInstall.log')
    if os.path.exists(logfile):
        os.remove(logfile)
    log.setLevel(logging.DEBUG)
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    fh = logging.FileHandler(logfile)
    fh.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    log.addHandler(fh)
    log.addHandler(ch)

    ### Now try to update urls from the website

    online_urls='http://goodsell.pages.in2p3.fr/bsmart/json/HEPurls.json'
    urls_json='HEPurls.json'

    if not args.LocalUrls:
        log.info("Updating paths ...")
        try:
            os.chdir(scriptdir)
            filename=wget.download(online_urls)
            #print(filename)
            if filename !=urls_json:
                shutil.move(filename,urls_json)
        except: ## doesn't matter at this point
            os.chdir(cwd)
            print('Could not update paths')

    # Global latest_urls already refers to default_urls

    if not args.LocalUrls and os.path.exists(urls_json):
        try:
            with open(urls_json) as json_data:
                http_paths = json.load(json_data)
            latest_urls=http_paths
            
        except Exception as e:
            #log.error('Failed to load json file '+file)
            print('Failed to load json file '+urls_json)
            print('Json exception given as ' +str(e))
            latest_urls=default_urls
    else:
        latest_urls=default_urls

    os.chdir(cwd)
    ## Load the paths
    if args.local:
        jsonFILE = 'HEPtoolpaths.json'
        print('Using local configuration file: HEPtoolpaths.json')
        # We might want to warn if a global one exists?
    else:
        jsonFILE = get_heptools_path_file(create_dir=True)
        print('Using configuration file: ' + jsonFILE)

    if os.path.exists(jsonFILE):

        try:
            with open(jsonFILE) as json_data:
                # paths_dict global
                paths_dict.update(json.load(json_data, object_pairs_hook=collections.OrderedDict))
        except Exception as e:
            #log.error('Failed to load json file '+file)
            print('Failed to load json file '+jsonFILE)
            print('Json exception given as ' +str(e))
            raise SystemExit
    
    # func_dict population
    func_dict=OrderedDict()

    if args.All or args.SARAH:
        func_dict['SARAH']=go_SARAH

    if args.All or args.SPheno:
        func_dict['SPheno']=go_SPheno

    if args.All or args.MicrOMEGAs:
        func_dict['MicrOMEGAs']=go_MicrOMEGAs

    if args.HiggsBounds:
        func_dict['HiggsBounds']=go_HiggsBounds
    
    if args.HiggsSignals:
        func_dict['HiggsSignals']=go_HiggsSignals

    if args.All or args.HiggsTools:
        func_dict['HiggsTools']=go_HiggsTools

    if args.All or args.VevaciousPlusPlus:
        func_dict['VevaciousPlusPlus']=go_VevaciousPlusPlus   

    if args.All or args.SModelS:
        func_dict['SModelS']=go_SModelS
    
    if args.All or args.ZPEED:
        func_dict['ZPEED']=go_ZPEED

    if args.All or args.MultiNest:
        func_dict['MultiNest']=go_MultiNest          

    if args.All or args.Diver:
        func_dict['Diver']=go_Diver           
    
    if args.All or args.FlexibleSUSY:
        func_dict['COLLIER']=go_COLLIER
        func_dict['LoopTools']=go_LoopTools
        func_dict['FlexibleSUSY']=go_FlexibleSUSY

    #if args.All or args.BSMArt:
    #    func_dict['BSMArt']=go_BSMArt     

    for tool in func_dict:
        #print('Calling: '+tool)
        try:
            func_dict[tool]()
        except Exception as e:
            log.error('Failed to install '+tool+', exception: '+str(e))
        #print(installed_tools)
        #print(installed_something)

        
    if len(installed_tools) > 0:
        os.chdir(cwd)
        with open(jsonFILE, 'w') as fp:
            json.dump(paths_dict, fp, indent=4)
        print('Tools installed:')
        for tool in installed_tools:
            print(' * '+tool)
    print("All done!")

if __name__ == "__main__":
    main()
