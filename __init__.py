#!/usr/bin/env python2
'''
# ================================================================================#
#-- PROJECT NAME : Diarization package using i-vector clustering. 
#-- BACKEND : LIUM Speaker Diarization
#-- TASK : This is the configuration file, implemented as singleton. 
		   This file has the paths for all dependencies

#-- Author : Sruthi.S
#-- Date : September 27th, 2016
# ================================================================================#
'''

class VConf(object):
    """"Configuration for Voiceid, implemented as singleton"""
    __instance = None

    def __new__(cls, *args, **kwargs):
        if not cls.__instance:
            cls.__instance = super(VConf, cls).__new__(
                                cls, *args, **kwargs)
        return cls.__instance

    def __init__(self, *args, **kwargs):
        object.__init__(self, *args, **kwargs)
        self.QUIET_MODE = False
        self.VERBOSE = False
        self.KEEP_INTERMEDIATE_FILES = True
        self.LIUM_JAR = './LIUM_SpkDiarization-8.4.1.jar'
        self.UBM_PATH = './models/ubm.gmm'
        self.DB_DIR = './gmm_db'
        self.GENDER_GMMS = './models/gender.gmms'
        self.SMS_GMMS = './models/sms.gmms'
        self.S_GMMS = './models/s.gmms'
        self.OUTPUT_FORMAT = 'srt'
