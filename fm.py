#!/usr/bin/env python2
'''
# ================================================================================#
#-- PROJECT NAME : Diarization package using i-vector clustering. 
#-- BACKEND : LIUM Speaker Diarization
#-- TASK : This module contains the the low level file manipulation functions

#-- Author : Sruthi.S
#-- Date : September 27th, 2016
# ================================================================================#
'''

import os
import re
import struct
from __init__ import VConf
import utils
import sys

CONFIGURATION = VConf()

JAVA_MEM1 = '1G'
JAVA_MEM2 = '2024'
JAVA_EXE = 'java'

def wave_duration(wavfile):
    """Extract the duration of a wave file in sec.

    :type wavfile: string
    :param wavfile: the wave input file"""
    import wave
    w_file = wave.open(wavfile)
    par = w_file.getparams()
    w_file.close()
    return par[3] / par[2]


def merge_waves(input_waves, wavename):
    """Take a list of waves and append them to a brend new destination wave.

    :type input_waves: list
    :param input_waves: the wave files list

    :type wavename: string
    :param wavename: the output wave file to be generated"""
    #if os.path.exists(wavename):
            #raise Exception("File gmm %s already exist!" % wavename)
    waves = [w_names.replace(" ", "\ ") for w_names in input_waves]
    w_names = " ".join(waves)
    commandline = "sox " + str(w_names) + " " + str(wavename)
    
    utils.start_subprocess(commandline)


def file2wav(filename):
    """Take any kind of video or audio and convert it to a
    "RIFF (little-endian) data, WAVE audio, Microsoft PCM, 16 bit,
    mono 16000 Hz" wave file using gstreamer. If you call it passing a wave it
    checks if in good format, else it converts the wave in the good format.

    :type filename: string
    :param filename: the input audio/video file to convert"""
    name, ext = os.path.splitext(filename)
    if ext == '.wav' and utils.is_good_wave(filename):
        pass
    else:
        if ext == '.wav':
            name += '_'
        utils.start_subprocess("gst-launch filesrc location='" + filename
           + "' ! decodebin ! audioresample ! 'audio/x-raw-int,rate=16000' !"
           + " audioconvert !"
           + " 'audio/x-raw-int,rate=16000,depth=16,signed=true,channels=1' !"
           + " wavenc ! filesink location=" + name + ".wav ")
    utils.ensure_file_exists(name + '.wav')
    return name + ext


def merge_gmms(input_files, output_file):
    """Merge two or more gmm files to a single gmm file with more voice models.

    :type input_files: list
    :param input_files: the gmm file list to merge

    :type output_file: string
    :param output_file: the merged gmm output file"""
    num_gmm = 0
    gmms = ''
    for ifile in input_files:
        try:
            current_f = open(ifile, 'rb')
        except (IOError, OSError):
            continue
        kind = current_f.read(8)
        if kind != 'GMMVECT_':
            raise Exception('different kinds of models!')
        
        line = current_f.read(4)
        num = struct.unpack('>i', line)
        num_gmm += int(num[0])
        byte = current_f.read(1)
        all_other = ""
        
        while byte:
            # Do stuff with byte.
            all_other += byte
            byte = current_f.read(1)
            
        gmms += all_other
        current_f.close()

    num_gmm_string = struct.pack('>i', num_gmm)
    os.system('echo "hello" | sudo chmod 777 '+output_file)
    new_gmm = open(output_file, 'w')
    new_gmm.write("GMMVECT_")
    new_gmm.write(num_gmm_string)
    new_gmm.write(gmms)

    new_gmm.close()


def get_gender(input_file):
    """Return gender of a given gmm file.

    :type input_file: string
    :param input_file: the gmm file
    """
    gmm = open(input_file, 'rb')
    gmm.read(8)  # kind
    num_gmm_string = gmm.read(4)
    num_gmm = struct.unpack('>i', num_gmm_string)

    if num_gmm != (1,):
        raise Exception('Loop needed for gmms')

    gmm.read(8)  # gmm_1
    gmm.read(4)  # nothing

    str_len = struct.unpack('>i', gmm.read(4))
    gmm.read(str_len[0]) # name

    gender = gmm.read(1)
    return gender


def _read_gaussian(g_file):
    """Read a gaussian in binary format"""
    full = 0
    g_key = g_file.read(8)     # read string of 8bytes kind
    if g_key != 'GAUSS___':
        raise Exception("Error: the gaussian is not"
                        + " of GAUSS___ key  (%s)" % g_key)
    g_id = g_file.read(4)
    g_length = g_file.read(4)  # readint 4bytes representing the name len
    g_name = g_file.read(int(struct.unpack('>i', g_length)[0]))
    g_gender = g_file.read(1)
    g_kind = g_file.read(4)
    g_dim = g_file.read(4)
    g_count = g_file.read(4)
    g_weight = g_file.read(8)
    dimension = int(struct.unpack('>i', g_dim)[0])
    g_header = g_key + g_id + g_length + g_name + g_gender + g_kind + g_dim
    g_header = g_header + g_count + g_weight
    datasize = 0
    if g_kind == full:
        for jay in range(dimension) :
            datasize += 1
            tee = jay
            while tee < dimension :
                datasize += 1
                tee += 1
    else:
        for jay in range(dimension) :
            datasize += 1
            tee = jay
            while tee < jay + 1 :
                datasize += 1
                tee += 1
    return g_header + g_file.read(datasize * 8)


def _read_gaussian_container(g_file):
    """Read a gaussian container in a binary format"""
    # gaussian container
    chk = g_file.read(8)   # read string of 8bytes
    if chk != "GAUSSVEC":
        raise Exception("Error: the gaussian container" +
                        " is not of GAUSSVEC kind %s" % chk)
    gcs = g_file.read(4)   # readint 4bytes representing the size of
    # the gaussian container
    stuff = chk + gcs
    for index in range(int(struct.unpack('>i', gcs)[0])):
        stuff += _read_gaussian(g_file)
    return stuff


def _read_gmm(g_file):
    """Read a gmm in a binary format"""
    myfile = {}
    # single gmm
    kind = g_file.read(8)     # read string of 8bytes kind
    if kind != "GMM_____":
        raise Exception("Error: Gmm section doesn't match GMM_____ kind")
    hash_ = g_file.read(4)  # readint 4bytes representing the hash (compatib)
    length = g_file.read(4)  # readint 4bytes representing the name length
    # read string of length bytes
    name = g_file.read(int(struct.unpack('>i', length)[0]))
    myfile['name'] = name
    gen = g_file.read(1)     # read a char representing the gender
    gaussk = g_file.read(4) # readint 4bytes representing the gaussian kind
    dim = g_file.read(4)   # readint 4bytes representing the dimension
    comp = g_file.read(4)  # readint 4bytes representing the number of comp
    gvect_header = kind + hash_ + length + name + gen + gaussk + dim + comp
    myfile['header'] = gvect_header
    myfile['content'] = _read_gaussian_container(g_file)
    return myfile


def split_gmm(input_file, output_dir=None):
    """Split a gmm file into gmm files with a single voice model.

    :type input_file: string
    :param input_file: the gmm file containing various voice models

    :type output_dir: string
    :param output_dir: the directory where is splitted the gmm input file"""


    g_file = open(input_file, 'rb')
    key = g_file.read(8)
    if key != 'GMMVECT_':  # gmm container
        raise Exception('Error: Not a GMMVECT_ file!')

    size = g_file.read(4)
    num = int(struct.unpack('>i', size)[0])  # number of gmms
    main_header = key + struct.pack('>i', 1)
    files = []
    for num in range(num):
        files.append(_read_gmm(g_file))
    g_file.close()
    file_basename = input_file[:-4]
    index = 0
    basedir, filename = os.path.split(file_basename)
    if output_dir != None:
        basedir = output_dir
        for g_file in files:
            newname = os.path.join(basedir, "%s%04d.gmm" % (filename, index))
            gmm_f = open(newname, 'w')
            gmm_f.write(main_header + g_file['header'] + g_file['content'])
            gmm_f.close()
            index += 1


def rename_gmm(input_file, identifier):
    """Rename a gmm with a new speaker identifier(name) associated.

    :type input_file: string
    :param input_file: the gmm file to rename

    :type identifier: string
    :param identifier: the new name or identifier of the gmm model"""
    gmm = open(input_file, 'rb')
    new_gmm = open(input_file + '.new', 'w')
    kind = gmm.read(8)
    new_gmm.write(kind)
    num_gmm_string = gmm.read(4)
    num_gmm = struct.unpack('>i', num_gmm_string)
    if num_gmm != (1,):
        raise Exception('Loop needed for gmms')
    new_gmm.write(num_gmm_string)
    gmm_1 = gmm.read(8)
    new_gmm.write(gmm_1)
    nothing = gmm.read(4)
    new_gmm.write(nothing)
    str_len = struct.unpack('>i', gmm.read(4))
    name = gmm.read(str_len[0])
    new_len = struct.pack('>i', len(identifier))
    new_gmm.write(new_len)
    new_gmm.write(identifier)
    all_other = gmm.read()
    new_gmm.write(all_other)
    gmm.close()
    new_gmm.close()

def build_gmm(filebasename, identifier):
    """Build a gmm (Gaussian Mixture Model) file from a given wave with a
    speaker identifier  associated.

    :type filebasename: string
    :param filebasename: the input file basename

    :type identifier: string
    :param identifier: the name or identifier of the speaker"""

    diarization_standard(filebasename)
    ident_seg(filebasename, identifier)
    _train_init(filebasename)
    _train_map(filebasename)

def seg2trim(filebasename):
    """Take a wave and splits it in small waves in this directory structure
    <file base name>/<cluster>/<cluster>_<start time>.wav

    :type filebasename: string
    :param filebasename: filebasename of the seg and wav input files"""
    segfile = filebasename + '.ev_is.seg'
    seg = open(segfile, 'r')
    for line in seg.readlines():
        if not line.startswith(";;"):
            arr = line.split()
            clust = arr[7]
            start = float(arr[2]) / 100
            end = float(arr[3]) / 100
            try:
                mydir = os.path.join(filebasename, clust)
                os.makedirs(mydir)
            except os.error, err:
                if err.errno == 17:
                    pass
                else:
                    raise os.error
            wave_path = os.path.join(filebasename, clust,
                                     "%s_%07d.%07d.wav" % (clust, int(start),
                                                           int(end)))
            
            commandline = "sox %s.wav %s trim  %s %s" % (filebasename,
                                                         wave_path,
                                                         start, end)    
                        
            utils.start_subprocess(commandline)
            utils.ensure_file_exists(wave_path)
    seg.close()

def seg2srt(segfile):
    """Take a seg file and convert it in a subtitle file (srt).

    :type segfile: string
    :param segfile: the segmentation file to convert"""
    def readtime(aline):
        "Help function for sort, to extract time from line"
        return int(aline[2])

    basename = os.path.splitext(segfile)[0]
    seg = open(segfile, 'r')
    lines = []
    for line in seg.readlines():
        if not line.startswith(";;"):
            arr = line.split()
            lines.append(arr)
    seg.close()
    lines.sort(key=readtime, reverse=False)
    fileoutput = basename + ".srt"
    srtfile = open(fileoutput, "w")
    row = 0
    for line in lines:
        row = row + 1
        start = float(line[2]) / 100
        end = start + float(line[3]) / 100
        srtfile.write(str(row) + "\n")
        srtfile.write(utils.humanize_time(start) + " --> "
                      + utils.humanize_time(end) + "\n")
        srtfile.write(line[7] + "\n")
        srtfile.write("\n")
    srtfile.close()
    utils.ensure_file_exists(basename + '.srt')

def ident_seg(filebasename, identifier):
    """Substitute cluster names with speaker names ang generate a
    "<filebasename>.ident.seg" file."""
    ident_seg_rename(filebasename, identifier, filebasename + '.ident')


def ident_seg_rename(filebasename, identifier, outputname):
    """Take a seg file and substitute the clusters with a given name or
    identifier."""
    seg_f = open(filebasename + '.seg', 'r')
    clusters = []
    lines = seg_f.readlines()
    for line in lines:
        for key in line.split():
            if key.startswith('cluster:'):
                clu = key.split(':')[1]
                clusters.append(clu)
    seg_f.close()
    output = open(outputname + '.seg', 'w')
    clusters.reverse()
    for line in lines:
        for clu in clusters:
            line = line.replace(clu, identifier)
        output.write(line)
    output.close()
    utils.ensure_file_exists(outputname + '.seg')

def srt2subnames(filebasename, key_value):
    """Substitute cluster names with real names in subtitles."""

    def replace_words(text, word_dic):
        """Take a text and replace words that match a key in a dictionary with
        the associated value, return the changed text"""
        rec = re.compile('|'.join(map(re.escape, word_dic)))

        def translate(match):
            "not documented"
            return word_dic[match.group(0)] + '\n'

        return rec.sub(translate, text)

    file_original_subtitle = open(filebasename + ".srt")
    original_subtitle = file_original_subtitle.read()
    file_original_subtitle.close()
    key_value = dict(map(lambda (key, value): (str(key) + "\n", value),
                         key_value.items()))
    text = replace_words(original_subtitle, key_value)
    out_file = filebasename + ".ident.srt"
    # create a output file
    fout = open(out_file, "w")
    fout.write(text)
    fout.close()
    utils.ensure_file_exists(out_file)

def file2trim(filename):
    """Take a video or audio file and converts it into smaller waves according
    to the diarization process.

    :type filename: string
    :param filename: the input file audio/video"""
    if not CONFIGURATION.QUIET_MODE:
        print "*** converting video to wav ***"
    file2wav(filename)
    file_basename = os.path.splitext(filename)[0]
    if not CONFIGURATION.QUIET_MODE:
        print "*** diarization ***"
    diarization(file_basename)
    if not CONFIGURATION.QUIET_MODE:
        print "*** trim ***"
    seg2trim(file_basename)


def diarization_standard(filebasename):
    """Take a wave file in the correct format and build a segmentation file.
    The seg file shows how much speakers are in the audio and when they talk.

    :type filebasename: string
    :param filebasename: the basename of the wav file to process"""

    utils.start_subprocess(JAVA_EXE +' -Xmx' + JAVA_MEM1 + ' -jar '
                           + CONFIGURATION.LIUM_JAR
                           + ' fr.lium.spkDiarization.system.Diarization '
                           +'--fInputMask=./%s.wav --sOutputMask=./%s.seg --doCEClustering '
                           +filebasename)

    utils.ensure_file_exists(filebasename + '.seg')


def diarization(filebasename, h='3', c='75'):

    fDescStart="audio16kHz2sphinx,1:1:0:0:0:0,13,0:0:0"
    fDesc="sphinx,1:1:0:0:0:0,13,0:0:0"
    fDescD="sphinx,1:3:2:0:0:0,13,0:0:0:0"
    fDescLast="sphinx,1:3:2:0:0:0,13,1:1:0:0"
    fDescCLR="sphinx,1:3:2:0:0:0,13,1:1:300:4"

    if not os.path.isdir(filebasename):
        os.mkdir(filebasename)

    print 'compute the MFCC'

    utils.start_subprocess(JAVA_EXE +' -Xmx' + JAVA_MEM1 + ' -classpath '
                           + CONFIGURATION.LIUM_JAR
                           + ' fr.lium.spkDiarization.tools.Wave2FeatureSet --help '
                           + '--fInputMask=%s.wav --fInputDesc='+fDescStart 
                           + ' --fOutputMask=%s.mfcc --fOutputDesc='+fDesc
                           + ' '+filebasename)

    utils.ensure_file_exists(filebasename + '.mfcc')

    print 'check the MFCC'

    utils.start_subprocess(JAVA_EXE +' -Xmx' + JAVA_MEM1 + ' -classpath '
                           + CONFIGURATION.LIUM_JAR
                           + ' fr.lium.spkDiarization.programs.MSegInit   --help '
                           + '--fInputMask=%s.mfcc --fInputDesc='+fDesc 
                           + ' --sInputMask= --sOutputMask=%s.i.seg '
                           + filebasename)

    utils.ensure_file_exists(filebasename + '.i.seg')

    print 'GLR based segmentation, make small segments'

    utils.start_subprocess(JAVA_EXE +' -Xmx' + JAVA_MEM1 + ' -classpath '
                           + CONFIGURATION.LIUM_JAR
                           + ' fr.lium.spkDiarization.programs.MSeg --kind=FULL '
                           + '--sMethod=GLR  --help --fInputMask=%s.mfcc '
                           + '--fInputDesc='+fDesc+' --sInputMask=%s.i.seg '
                           + '--sOutputMask=%s.s.seg ' + filebasename)

    utils.ensure_file_exists(filebasename + '.s.seg')

    print 'Segmentation: linear clustering'
    l = '2'
    utils.start_subprocess(JAVA_EXE +' -Xmx' + JAVA_MEM1 + ' -classpath '
                           + CONFIGURATION.LIUM_JAR
                           + ' fr.lium.spkDiarization.programs.MClust --help '
                           + '--fInputMask=%s.mfcc --fInputDesc='+fDesc 
                           + ' --sInputMask=%s.s.seg --sOutputMask=%s.l.seg '
                           + '--cMethod=l --cThr='+l+ ' '+filebasename)

    utils.ensure_file_exists(filebasename + '.l.seg')

    print 'hierarchical clustering'

    utils.start_subprocess(JAVA_EXE +' -Xmx' + JAVA_MEM1 + ' -classpath '
                           + CONFIGURATION.LIUM_JAR
                           + ' fr.lium.spkDiarization.programs.MClust --help '
                           + '--fInputMask=%s.mfcc --fInputDesc='+fDesc
                           + ' --sInputMask=%s.l.seg --sOutputMask=%s.h.seg '
                           + ' --cMethod=h --cThr='+h+ ' '+ filebasename)

    utils.ensure_file_exists(filebasename + '.h.seg')

    print 'initialize GMM'

    utils.start_subprocess(JAVA_EXE +' -Xmx' + JAVA_MEM1 + ' -classpath '
                           + CONFIGURATION.LIUM_JAR
                           + ' fr.lium.spkDiarization.programs.MTrainInit --help '
                           + '--nbComp=8 --kind=DIAG --fInputMask=%s.mfcc '
                           + '--fInputDesc='+fDesc+ ' --sInputMask=%s.h.seg '
                           + '--tOutputMask=%s.init.gmms '+filebasename)

    utils.ensure_file_exists(filebasename + '.init.gmms')

    print 'EM computation'

    utils.start_subprocess(JAVA_EXE +' -Xmx' + JAVA_MEM1 + ' -classpath '
                           + CONFIGURATION.LIUM_JAR
                           + ' fr.lium.spkDiarization.programs.MTrainEM --help '
                           + '--nbComp=8 --kind=DIAG --fInputMask=%s.mfcc '
                           + '--fInputDesc='+fDesc+' --sInputMask=%s.h.seg '
                           + '--tOutputMask=%s.gmms --tInputMask=%s.init.gmms '
                           + filebasename)

    utils.ensure_file_exists(filebasename + '.gmms')

    print 'Viterbi decoding'

    utils.start_subprocess(JAVA_EXE +' -Xmx' + JAVA_MEM1 + ' -classpath '
                           + CONFIGURATION.LIUM_JAR
                           + ' fr.lium.spkDiarization.programs.MDecode --help '
                           + '--fInputMask=%s.mfcc --fInputDesc='+fDesc
                           + ' --sInputMask=%s.h.seg --sOutputMask=%s.d.seg '
                           + '--dPenality=250 --tInputMask=%s.gmms '
                           + filebasename)

    utils.ensure_file_exists(filebasename + '.d.seg')

    print 'Speech/Music/Silence segmentation'

    utils.start_subprocess(JAVA_EXE +' -Xmx' + JAVA_MEM1 + ' -classpath '
                           + CONFIGURATION.LIUM_JAR
                           + ' fr.lium.spkDiarization.programs.MDecode --help '
                           + '--fInputDesc='+fDescD+' --fInputMask=%s.mfcc '
                           + '--sInputMask=%s.i.seg --sOutputMask=%s.pms.seg '
                           + '--dPenality=10,10,50 --tInputMask='+CONFIGURATION.SMS_GMMS
                           + ' '+filebasename)

    utils.ensure_file_exists(filebasename + '.pms.seg')

    print 'filter spk segmentation according pms segmentation'

    utils.start_subprocess(JAVA_EXE +' -Xmx' + JAVA_MEM1 + ' -classpath '
                           + CONFIGURATION.LIUM_JAR
                           + ' fr.lium.spkDiarization.tools.SFilter --help '
                           + '--fInputDesc='+fDescD+' --fInputMask=%s.mfcc '
                           + '--fltSegMinLenSpeech=150 --fltSegMinLenSil=25 --sFilterClusterName=j '
                           + '--fltSegPadding=25 --sFilterMask=%s.pms.seg --sInputMask=%s.d.seg '
                           + '--sOutputMask=%s.flt.seg '+filebasename)

    utils.ensure_file_exists(filebasename + '.flt.seg')

    print 'Set gender and bandwith'

    utils.start_subprocess(JAVA_EXE +' -Xmx' + JAVA_MEM1 + ' -classpath '
                           + CONFIGURATION.LIUM_JAR
                           + ' fr.lium.spkDiarization.programs.MScore --help '
                           + '--sGender --sByCluster --fInputDesc='+fDescLast
                           + ' --fInputMask=%s.mfcc --sInputMask=%s.flt.seg '
                           + '--sOutputMask=%s.g.seg --tInputMask='+CONFIGURATION.GENDER_GMMS
                           + ' '+filebasename)

    utils.ensure_file_exists(filebasename + '.g.seg')

    print 'ILP Clustering'

    utils.start_subprocess(JAVA_EXE +' -Xmx' + JAVA_MEM1 + ' -classpath '
                           + CONFIGURATION.LIUM_JAR
                           + ' fr.lium.spkDiarization.programs.ivector.ILPClustering --cMethod=es_iv '
                           + '--ilpThr='+c+' --help --sInputMask=%s.g.seg --sOutputMask=%s.ev_is.seg '
                           + '--fInputMask=%s.mfcc --fInputDesc='+fDescLast+' --tInputMask=./ubm/wld.gmm '
                           + '--nEFRMask=mat/wld.efn.xml --ilpGLPSolProgram=glpsol ' 
                           + '--nMahanalobisCovarianceMask=./mat/wld.mahanalobis.mat '
                           + '--tvTotalVariabilityMatrixMask=./mat/wld.tv.mat '
                           + '--ilpOutputProblemMask=%s.ilp.problem.txt '
                           + '--ilpOutputSolutionMask=%s.ilp.solution.txt '+filebasename)

    utils.ensure_file_exists(filebasename + '.ev_is.seg')

    if not CONFIGURATION.KEEP_INTERMEDIATE_FILES:

        f_list = os.listdir(filebasename)

        for ext in f_list:
            os.remove(ext)

def _train_init(filebasename):
    """Train the initial speaker gmm model."""
    utils.start_subprocess(JAVA_EXE +' -Xmx' + JAVA_MEM2 + 'm -cp '
                           + CONFIGURATION.LIUM_JAR
                           + ' fr.lium.spkDiarization.programs.MTrainInit --help '
                           + '--sInputMask=%s.ident.seg --fInputMask=%s.wav '
                           + '--fInputDesc="audio16kHz2sphinx,1:3:2:0:0:0,13,1:1:300:4" '
                           + '--emInitMethod=copy --tInputMask='+CONFIGURATION.UBM_PATH
                           + ' --tOutputMask=%s.init.gmm '+filebasename)

    utils.ensure_file_exists(filebasename + '.init.gmm')

def _train_map(filebasename):
    """Train the speaker model using a MAP adaptation method."""
    utils.start_subprocess(JAVA_EXE +' -Xmx' + JAVA_MEM2 + 'm -cp '
                           + CONFIGURATION.LIUM_JAR
                           + ' fr.lium.spkDiarization.programs.MTrainMAP --help '
                           + '--sInputMask=%s.ident.seg --fInputMask=%s.wav '
                           + '--fInputDesc="audio16kHz2sphinx,1:3:2:0:0:0,13,1:1:300:4" '
                           + '--tInputMask=%s.init.gmm --emCtrl=1,5,0.01 --varCtrl=0.01,10.0 '
                           + '--tOutputMask=%s.gmm '+filebasename)

def wav_vs_gmm(filebasename, gmm_file, gender, custom_db_dir=None):
    """Match a wav file and a given gmm model file and produce a segmentation
    file containing the score obtained.

    :type filebasename: string
    :param filebasename: the basename of the wav file to process

    :type gmm_file: string
    :param gmm_file: the path of the gmm file containing the voice model

    :type gender: char
    :param gender: F, M or U, the gender of the voice model

    :type custom_db_dir: None or string
    :param custom_db_dir: the voice models database to use"""

    database = CONFIGURATION.DB_DIR
    
    if custom_db_dir != None:
        database = custom_db_dir
    gmm_name = os.path.split(gmm_file)[1]

    print 'database', database
    print 'gmm_name', gmm_name

    # utils.start_subprocess(JAVA_EXE +' -Xmx256M -cp '+CONFIGURATION.LIUM_JAR
    #                        + ' fr.lium.spkDiarization.programs.MScore --sInputMask=%s.ev_is.seg '
    #                        + '--fInputMask=%s.wav --sOutputMask=%s.ident.' + gender + '.'
    #                        + gmm_name + '.ev_is.seg --sOutputFormat=seg,UTF8 '
    #                        + '--fInputDesc=audio16Khs2sphinx,1:3:2:0:0:0,13,1:0:300:4 '
    #                        + '--tInputMask=' + database + '/' + gender + '/' + gmm_file
    #                        + ' --sTop=8,' + CONFIGURATION.UBM_PATH
    #                        + '  --sSetLabel=add --sByCluster ' + filebasename)

    utils.start_subprocess('echo "hello" | sudo -S '+JAVA_EXE +' -Xmx256M -cp '
                           + CONFIGURATION.LIUM_JAR
                           +' fr.lium.spkDiarization.programs.MScore --sInputMask=%s.ev_is.seg '
                           + '--fInputMask=%s.wav --sOutputMask=%s.ident.' + gender + '.'
                           + gmm_name + '.ev_is.seg --sOutputFormat=seg,UTF8 '
                           + '--fInputDesc=audio2sphinx,1:3:2:0:0:0,13,1:0:300:4 --tInputMask='
                           + database + '/' + gender + '/' + gmm_file+' --sTop=8,'
                           +CONFIGURATION.UBM_PATH+ ' --sSetLabel=add --sByCluster '+filebasename)

    utils.ensure_file_exists(filebasename + '.ident.'
                             + gender + '.' + gmm_name + '.ev_is.seg')
