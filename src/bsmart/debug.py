""" Debugging routines modified from xBit by F. Staub """


import logging
import subprocess
import sys


def print_BSMArt_splash():
    from bsmart import __version__
    from bsmart import __build_date__
    BSMArtVersion=__version__
    BSMArtDate=__build_date__
    print('                                                                ')
    print('888888b.    .d8888b.  888b     d888        d8888         888    ')
    print('888  "88b  d88P  Y88b 8888b   d8888       d88888         888    ')
    print('888  .88P  Y88b.      88888b.d88888      d88P888         888    ')
    print('8888888K.   "Y888b.   888Y88888P888     d88P 888 888d888 888888 ')
    print('888  "Y88b     "Y88b. 888 Y888P 888    d88P  888 888P"   888    ')
    print('888    888       "888 888  Y8P  888   d88P   888 888     888    ')
    print('888   d88P Y88b  d88P 888   "   888  d8888888888 888     Y88b.  ')
    print('8888888P"   "Y8888P"  888       888 d88P     888 888      "Y888 ')
    print('                                                                ')
    print('Version %s [%s]' %(BSMArtVersion,BSMArtDate))
    print('By M. Goodsell                                                  ')
    print('With contributions from:                                        ')
    print(' * A. Joury (R. Moutafis)                                       ')
    print(' * Luc Darmé                                                    ')
    print(' * Martin Gabelmann                                             ')
    print(' * Johannes Braathen                                            ')
    print(' * Farid Ibrahimov                                              ')
    print(' * Wojciech Kotlarski                                           ')
    print(' * Miguel Crispim Romao                                         ')
    print(' * Fernando Abreu de Souza                                      ')
    print('                                                                ')



""" Test existence of psutil """

psutil_available=False
try:
    import psutil
    psutil_available=True
except:
    pass


# def new_logger(debug, name, file):
#     '''Create Logger'''
    
#     #logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s',level=logging.CRITICAL)
#     #logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s',level='critical')
#     BSMArtlog = logging.getLogger(name)
#     BSMArtlog.propagate=False ### don't send messages up to the root logger
#     BSMArtlog.setLevel(logging.DEBUG)
#     # file output
#     BSMfh = logging.FileHandler(file)
#     BSMch = logging.StreamHandler()
#     if debug:
#         BSMch.setLevel(logging.INFO)
#         BSMArtlog.setLevel(logging.DEBUG)
#     else:
#         BSMch.setLevel(logging.CRITICAL)
#         BSMfh.setLevel(logging.INFO)
#         #ch.setLevel(logging.ERROR)
        

#     ## could add threadname, but it's redundant
#     BSMArtformatter = \
#             logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
#     BSMfh.setFormatter(BSMArtformatter)
#     BSMch.setFormatter(BSMArtformatter)
    
#     BSMArtlog.addHandler(BSMfh)
#     BSMArtlog.addHandler(BSMch)
#     BSMArtlog.debug('Logger %s initialised' % name)
#     return BSMArtlog


def new_logger(debug, name, file):
    """ 
        Create the logger for scans
    """

    
    if debug:
        console_level=logging.INFO
        logfile_level=logging.DEBUG
    else:
        console_level=logging.CRITICAL
        logfile_level=logging.INFO
    
    return _new_general_logger(console_level,logfile_level,name,file)


def new_core_logger(debug, name, file):
    """ 
        Create the logger for the core running: we don't need to write much from the console in this case, keep it in the files. 
        Otherwise we clutter the screen with lots of error messages.
    """

    
    if debug:
        console_level=logging.CRITICAL
        logfile_level=logging.DEBUG
    else:
        console_level=logging.CRITICAL
        logfile_level=logging.INFO   # This could be changed
    
    return _new_general_logger(console_level,logfile_level,name,file)



def _new_general_logger(console_level,logfile_level,name,file):
    """ Private method for creating a logger with specified logging levels """
    BSMArtlog = logging.getLogger(name)
    BSMArtlog.propagate=False ### don't send messages up to the root logger
    BSMArtlog.setLevel(logging.DEBUG)
    # file output
    BSMfh = logging.FileHandler(file)
    BSMch = logging.StreamHandler()
    
    BSMch.setLevel(console_level)
    BSMArtlog.setLevel(logfile_level)
    
    ## could add threadname, but it's redundant
    BSMArtformatter = \
            logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    BSMfh.setFormatter(BSMArtformatter)
    BSMch.setFormatter(BSMArtformatter)
    
    BSMArtlog.addHandler(BSMfh)
    BSMArtlog.addHandler(BSMch)
    BSMArtlog.debug('Logger %s initialised' % name)
    return BSMArtlog






# def command_line_log(command, log, ctimeout=None):
#     log.debug('Terminal command: %s' % command)
#     command_line = subprocess.Popen(command, shell=True,
#                                     stdout=subprocess.PIPE,
#                                     stderr=subprocess.PIPE)
#     try:
#         out, err = command_line.communicate(timeout=ctimeout)
#         # log.debug(out)
#         if err != b'':  ## not necessarily an actual error?
#             log.error(err.decode("utf-8"))
#     except Exception as e:
#         command_line.kill()  ## necessary, otherwise child process can continue forever!
#         raise NameError(command+' timed out: '+str(e))


def command_line_log(command, log, ctimeout=None):
    """
    Execute commands in a shell, logging both the command and any errors.
    22/07/2024: discard the output of stdout since we never use it. 
    """
    log.debug('Terminal command: %s' % command)
    command_line = subprocess.Popen(command, shell=True,
                                    stdout=subprocess.DEVNULL,
                                    stderr=subprocess.PIPE)
    try:
        _, err = command_line.communicate(timeout=ctimeout)
        # log.debug(out)
        if err != b'':  ## not necessarily an actual error?
            log.error(err.decode("utf-8"))
    except Exception as e:
        command_line.kill()  ## necessary, otherwise child process can continue forever!
        raise NameError(command+' timed out: '+str(e))


def command_line_log_PID(command, wdir, log):
    log.debug('Terminal command: %s' % command)
    command_line = subprocess.Popen(command, shell=True, cwd=wdir,
                                    stdout=subprocess.DEVNULL,
                                    stderr=subprocess.PIPE)
    pid=command_line.pid
    _, err = command_line.communicate()
    # log.debug(out)
    if err != b'':
        log.error(err.decode("utf-8"))
    return str(pid)

def noshell_command_log_PID(command, wdir, log):
    log.debug('Terminal command: %s' % command)
    command_line = subprocess.Popen(command, shell=False, cwd=wdir,
                                    stdout=subprocess.DEVNULL,
                                    stderr=subprocess.PIPE)
    pid=command_line.pid
    _, err = command_line.communicate()
    # log.debug(out)
    if err != b'':
        log.error(err.decode("utf-8"))
    return str(pid)

def onecore_command_log(command, log, ctimeout=None):
    log.debug('Terminal command: %s' % command)
    command_line = psutil.Popen(command, shell=True,
                                    stdout=subprocess.DEVNULL,
                                    stderr=subprocess.PIPE)
    pid=command_line.pid
    cpuN=command_line.cpu_num()
    log.debug("Runnning on CPU number: "+str(cpuN))
    command_line.cpu_affinity([cpuN])
    # log.debug(out)
    try:
        _, err = command_line.communicate(timeout=ctimeout)
        # log.debug(out)
        if err != b'':  ## not necessarily an actual error
            log.error(err.decode("utf-8"))
    except Exception as e:
        command_line.kill()  ## necessary, otherwise child process can continue forever!
        raise NameError(command+' timed out: '+str(e))

