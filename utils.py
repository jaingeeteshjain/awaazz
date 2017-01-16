#!/usr/bin/env python2
'''
# ================================================================================#
#-- PROJECT NAME : Diarization package using i-vector clustering. 
#-- BACKEND : LIUM Speaker Diarization
#-- TASK : This utility file contains modules for subprocess, threading and file checking.

#-- Author : Sruthi.S
#-- Date : September 27th, 2016
# ================================================================================#
'''

from . import VConf
import os
import subprocess

CONFIGURATION = VConf()

def alive_threads(t_dict):
    """Check how much threads are running and alive in a thread dictionary

    :type t: dictionary
    :param t: thread dictionary like  { key : thread_obj, ... }"""
    num = 0
    for thr in t_dict:
        if t_dict[thr].is_alive():
            num += 1
    return num

def start_subprocess(commandline):
    """Start a subprocess using the given commandline and check for correct
    termination.

    :type commandline: string
    :param commandline: the command to run in a subprocess"""
    print commandline
    proc = subprocess.Popen(commandline, shell = True)
    retval = proc.wait()
    if retval != 0:
        err = OSError("Subprocess %s closed unexpectedly [%s]" % (str(proc),
                                                                  commandline))
        err.errno = retval
        raise err

def check_cmd_output(command):
    "Run a shell command and return the result as string"
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT,
                               universal_newlines=True)
    output = process.communicate()
    retcode = process.poll()
    if retcode:
        raise subprocess.CalledProcessError(retcode,
                                            command)
                                            # output=output[0]???
    return output

def ensure_file_exists(filename):
    """Ensure file exists and is not empty, otherwise raise an IOError.

    :type filename: string
    :param filename: file to check"""
    if not os.path.exists(filename):
        raise IOError("File %s doesn't exist or not correctly created" 
                      % filename)
    if not (os.path.getsize(filename) > 0):
        raise IOError("File %s empty" % filename)
    
    (shortname, extension) = os.path.splitext(filename)

def is_good_wave(filename):
    """Check if the wave is in correct format for LIUM.

    :type filename: string
    :param filename: file to check"""
    import wave
    par = None
    try:
        w_file = wave.open(filename)
        par = w_file.getparams()
        w_file.close()
    except wave.Error, exc:
        print exc
        return False
    if par[:3] == (1, 2, 16000) and par[-1:] == ('not compressed',):
        return True
    else:
        return False

def check_deps():
    """Check for dependencies."""
    ensure_file_exists(CONFIGURATION.LIUM_JAR)

    dir_m = os.path.join(CONFIGURATION.DB_DIR, "M")
    dir_f = os.path.join(CONFIGURATION.DB_DIR, "F")
    dir_u = os.path.join(CONFIGURATION.DB_DIR, "U")
    ensure_file_exists(CONFIGURATION.UBM_PATH)
    if not os.path.exists(CONFIGURATION.DB_DIR):
        os.makedirs(CONFIGURATION.DB_DIR)
    if os.listdir(CONFIGURATION.DB_DIR) == []:
        print("WARNING: Gmm db directory found in %s is empty" 
                % CONFIGURATION.DB_DIR)
    if not os.path.exists(dir_m):
        os.makedirs(dir_m)
    if not os.path.exists(dir_f):
        os.makedirs(dir_f)
    if not os.path.exists(dir_u):
        os.makedirs(dir_u)

def humanize_time(secs):
    """Convert seconds into time format.

    :type secs: integer
    :param secs: the time in seconds to represent in human readable format
           (hh:mm:ss)"""
    mins, secs = divmod(secs, 60)
    hours, mins = divmod(mins, 60)
    return '%02d:%02d:%02d,%s' % (hours, mins, int(secs), 
                                  str(("%0.3f" % secs))[-3:])
