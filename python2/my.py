#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import argparse
import collections
from ConfigParser import SafeConfigParser
from contextlib import closing
import csv
from glob import glob
import gzip
import job
import math
import re
import subprocess
import tempfile
from time import localtime, strftime
from urllib import urlretrieve
from warnings import warn
from zipfile import ZipFile
try:
    import pysam
except ImportError:
    raise ImportError,"The Python pysam module is required to run this program. Try 'pip install pysam'."
try:
    from bs4 import BeautifulSoup
except ImportError:
    raise ImportError,"The Python BeautifulSoup module is required to run this program. Try 'pip install beautifulsoup4'."
try:
    from ftputil import FTPHost
except ImportError:
    raise ImportError,"The Python ftputil module is required to run this program. Try 'pip install ftputil'."
try:
    import MySQLdb
    import MySQLdb.cursors
except ImportError:
    raise ImportError,"The Python MySQLdb module is required to run this program. Try 'pip install MySQL-python'."
try:
    import requests
except ImportError:
    raise ImportError,"The Python requests module is required to run this program. Try 'pip install requests'."
try:
    import sh
    from sh import git
except ImportError:
    raise ImportError,"The Python sh module is required to run this program. Try 'pip install sh'."

def file_len(fname):
    p = subprocess.Popen(['wc', '-l', fname], stdout=subprocess.PIPE, 
                         stderr=subprocess.PIPE)
    result, err = p.communicate()
    if p.returncode != 0:
        raise IOError(err)
    return int(result.strip().split()[0])

class QseqError(Exception): pass
class MyError(Exception): pass
class GetMetadataError(Exception): pass
class RunJobError(Exception): pass
class CheckFilesError(Exception): pass
class FastqError(Exception): pass
class SqlColumnsError(Exception): pass

class LaneRead(object):
    def __init__(self, read=None, lane_qseq=None, files_to_merge=None,
                 tile_files_compressed=None, merge_cmd=None, ):
        self.read = read
        self.lane_qseq = lane_qseq
        self.files_to_merge = files_to_merge if files_to_merge else []
        self.tile_files_compressed = tile_files_compressed
        self.merge_cmd = merge_cmd

class LaneReads(dict):
    def __init__(self, **kwds):
        super(LaneReads,self).__init__(**kwds)

class Lane(object):
    def __init__(self, machine = None, run = None, lane = None,
                 barcodes = None, readgroups = None, lane_qseqs_dir = None, reads = None,
                 tile_qseqs = None, lane_qseqs_demux_dir=None, demultiplex_command = None,
                 lane_qseqs_demux_move_commands = None, lane_qseqs_remove_intermediates_command = None):

        self.machine = machine
        self.run = run
        self.lane = lane
        self.barcodes = barcodes if barcodes else []
        self.readgroups = readgroups if readgroups else []
        self.lane_qseqs_dir = lane_qseqs_dir
        self.reads = reads if reads else LaneReads()
        self.tile_qseqs = tile_qseqs if tile_qseqs else []
        self.lane_qseqs_demux_dir = lane_qseqs_demux_dir
        self.demultiplex_command = demultiplex_command
        self.lane_qseqs_demux_move_commands = lane_qseqs_demux_move_commands if lane_qseqs_demux_move_commands else []
        self.lane_qseqs_remove_intermediates_command = lane_qseqs_remove_intermediates_command

class Lanes(dict):
    def __init__(self, **kwds):
        super(Lanes,self).__init__(**kwds)

ReadgroupPk = collections.namedtuple('ReadgroupPk', 'machine, run, lane, barcode')

class ReadGroup(object):
    def __init__(self, do_not_align=False, sample=None, original_sample=None, library=None, readgroup_id=None,
                 machine=None, run=None, lane=None, sample_sheet_barcode=None, barcode=None,
                 owner=None, email=None, collaborator=None, collaborator_email=None,
                 study=None, cohort=None, family=None, individual=None, father=None, mother=None, sex=None,
                 tumor_status=None, species=None, sample_source=None,
                 multiplexed_library=None, library_protocol=None, bioanalyzer_mean_bp=None, predicted_median_insert_size=None, median_insert_size=None,
                 capture_bait_interval_list=None, genomic_target_interval_list=None,
                 barcode_tag_list=None, sequencing_library=None,
                 adapter_1=None, adapter_2=None, experiment_protocol=None, experiment=None,
                 sequencing_center=None, platform=None, instrument_model=None,
                 run_date=None, run_folder=None, flowcell=None, flowcell_type=None,
                 index_format=None, index1=None, index2=None, platform_unit=None,
                 read_1_length=None, read_2_length=None,
                 readgroup_description=None, additional_readgroup_description=None,
                 run_top_dir=None, run_dir=None, fastq_top_dir=None, fastq_glob=None,
                 #readgroup_qseqs=None, readgroup_qseqs_read_lengths=None, tile_qseq_created_date=None,
                 novoalign_cmds=None, readgroup_bam_parts=None, readgroup_bam=None,  novoalign_bam_validate=None,
                 fixmate_bam_parts=None, fixmate_bam=None, fixmate_bams=None, fixmate_baits=None, fixmate_bai_parts=None, fixmate_bai=None, fixmate_bam_validate_parts=None, fixmate_bam_validate=None,
                 library_bam=None, library_bai=None, library_bam_validate=None,
                 library_markdup_bam=None, library_markdup_bai=None, library_markdup_bam_validate=None,
                 sample_bam=None, sample_bai=None, sample_bam_validate=None,
                 study_bam=None, study_bai=None, study_bam_validate=None,
                 fastq_reads=None, fastq_peek=None,
                 RG=None):

        self.do_not_align = bool(do_not_align)
        self.sample = sample
        self.original_sample = original_sample
        self.library = library
        self.readgroup_id = readgroup_id
        self.machine = machine
        self.run = run
        self.lane = lane
        self.sample_sheet_barcode = sample_sheet_barcode
        self.barcode = barcode
        self.owner = owner
        self.email = email
        self.collaborator = collaborator
        self.collaborator_email = collaborator_email
        self.study = study
        self.cohort = cohort
        self.family = family
        self.individual = individual
        self.father = father
        self.mother = mother
        self.sex = sex
        self.tumor_status = tumor_status
        self.species = species
        self.sample_source = sample_source
        self.multiplexed_library = multiplexed_library
        self.library_protocol = library_protocol
        self.bioanalyzer_mean_bp = bioanalyzer_mean_bp
        self.predicted_median_insert_size = predicted_median_insert_size
        self.median_insert_size = median_insert_size
        self.capture_bait_interval_list = capture_bait_interval_list
        self.genomic_target_interval_list = genomic_target_interval_list
        self.barcode_tag_list = barcode_tag_list
        self.sequencing_library = sequencing_library
        self.adapter_1 = adapter_1
        self.adapter_2 = adapter_2
        self.experiment_protocol = experiment_protocol
        self.experiment = experiment
        self.sequencing_center = sequencing_center
        self.platform = platform
        self.instrument_model = instrument_model
        self.run_date = run_date
        self.run_folder = run_folder
        self.flowcell = flowcell
        self.flowcell_type = flowcell_type
        self.index_format = index_format
        self.index1 = index1
        self.index2 = index2
        self.platform_unit = platform_unit
        self.read_1_length = read_1_length
        self.read_2_length = read_2_length
        self.readgroup_description = readgroup_description
        self.additional_readgroup_description = additional_readgroup_description
        self.run_top_dir = run_top_dir
        self.run_dir = run_dir
        self.fastq_top_dir = fastq_top_dir
        self.fastq_glob = fastq_glob
        #self.qseq_glob = qseq_glob
        #self.readgroup_qseqs = readgroup_qseqs if readgroup_qseqs else []
        #self.readgroup_qseqs_read_lengths = readgroup_qseqs_read_lengths if readgroup_qseqs_read_lengths else []
        #self.tile_qseq_created_date = tile_qseq_created_date
        self.novoalign_cmds = novoalign_cmds if novoalign_cmds else []
        self.novoalign_bam_validate = novoalign_bam_validate
        self.readgroup_bam_parts = readgroup_bam_parts if readgroup_bam_parts else {}
        self.readgroup_bam = readgroup_bam
        self.fixmate_bam_parts = fixmate_bam_parts if fixmate_bam_parts else {}
        self.fixmate_bam = fixmate_bam
        self.fixmate_bams = fixmate_bams if fixmate_bams else []
        self.fixmate_baits = fixmate_baits if fixmate_baits else set()
        self.fixmate_bai_parts = fixmate_bai_parts if fixmate_bai_parts else {}
        self.fixmate_bai = fixmate_bai
        self.fixmate_baits = fixmate_baits if fixmate_baits else []
        self.fixmate_bam_validate_parts = fixmate_bam_validate_parts if fixmate_bam_validate_parts else {}
        self.fixmate_bam_validate = fixmate_bam_validate
        self.library_bam = library_bam
        self.library_bai = library_bai
        self.library_bam_validate = library_bam_validate
        self.library_markdup_bam = library_markdup_bam
        self.library_markdup_bai = library_markdup_bai
        self.library_markdup_bam_validate = library_markdup_bam_validate
        self.sample_bam = sample_bam
        self.sample_bai = sample_bai
        self.sample_bam_validate = sample_bam_validate
        self.study_bam = study_bam
        self.study_bai = study_bai
        self.study_bam_validate = study_bam_validate
        self.fastq_reads = fastq_reads if fastq_reads else {}
        self.fastq_peek = fastq_peek
        self.RG = RG

class ReadGroups(dict):
    def __init__(self, **kwds):
        super(ReadGroups,self).__init__(**kwds)

def makedir(dir):
    try:
        os.makedirs(dir)
    except OSError:
        if os.path.isdir(dir):
            pass
        else:
            raise
    finally:
        return dir

def split_file(file, lines, suffix_length=8, prefix=None, repeat_header=False, compress=False):
    #returns list of output files
    if not file_exists(file):
        raise MyError("split_file can't find {}".format(file))
    if not is_number(lines) or lines < 1:
        raise MyError("split_file lines {} invalid".format(lines))
    if not is_number(suffix_length) or suffix_length < 1:
        raise MyError("split_file suffix_length {} invalid".format(suffix_length))

    #return value
    parts = []

    if not prefix:
        prefix = file
    file_suffix = prefix+'{:0='+str(suffix_length)+'}'+('.gz' if compress else '')
    lines_out = 0
    suffix = 0
    file_out = None
    header = ''
    if repeat_header:
        with open_gz_or_text(file) as f:
            for line in f:
                if line.startswith('#'):
                    header += line
                else:
                    break
    with open_gz_or_text(file) as f:
        for line in f:
            if not file_out or lines_out >= lines:
                lines_out = 0
                if file_out:
                    file_out.close()
                filename_out = file_suffix.format(suffix)
                suffix += 1
                parts += [filename_out]
                if compress:
                    file_out = gzip.open(filename_out, 'w')
                else:
                    file_out = open(filename_out, 'w')
                if repeat_header:
                    file_out.write(header)
            if repeat_header and line.startswith('#'):
                continue
            file_out.write(line)
            lines_out += 1

    if file_out:
        file_out.close()

    return parts

def get_immediate_subdirectory_paths(dir):
    return [os.path.join(dir, name) for name in os.listdir(dir) if os.path.isdir(os.path.join(dir, name))]

def get_immediate_subdirectory_names(dir):
    return [name for name in os.listdir(dir) if os.path.isdir(os.path.join(dir, name))]

def is_gz(filename):
    GZIP_MAGIC = b"\x1F\x8B"
    with open(filename, 'rb') as fh:
        magic = fh.read(len(GZIP_MAGIC))
    if magic == GZIP_MAGIC:
        return True
    else:
        return False

def open_gz_or_text(filename):
    if is_gz(filename):
        return gzip.open(filename, 'r')
    else:
        return open(filename, 'r')

def open_bam(filename):
    if is_gz(filename):
        return pysam.Samfile(filename, "rb")
    else:
        return pysam.Samfile(filename, "r")

def unique_dir(prefix='', dir=os.path.join(os.path.expanduser('~'),'tmp')):
    makedir(dir)
    return tempfile.mkdtemp(prefix=prefix+"_", dir=dir)

#http://stackoverflow.com/questions/183480/is-this-the-best-way-to-get-unique-version-of-filename-w-python
def unique_file(file_name):
    dirname, filename = os.path.split(file_name)
    prefix, suffix = os.path.splitext(filename)
    fd, filename = tempfile.mkstemp(suffix, prefix+"_", dirname)
    return os.fdopen(fd,'w'), filename

#deprecated
def list_files(pathnames):
    return sorted(set(flatten([glob(os.path.expanduser(os.path.expandvars(p))) for p in pathnames])))

def unglob(globs,sort_unique=True):
    if isinstance(globs, (basestring)):
        globs = [globs]
    if sort_unique:
        return sorted(set(flatten([glob(os.path.expanduser(os.path.expandvars(p))) for p in globs]))) if globs and isinstance(globs, (list, tuple, set)) else []
    else:
        return flatten([glob(os.path.expanduser(os.path.expandvars(p))) for p in globs]) if globs and isinstance(globs, (list, tuple, set)) else []

def get_files(names):
    return (file for file in names if os.path.isfile(file))

def is_int(s):
    try:
        n = int(s)
        if (math.isinf(n) or math.isnan(n)): return False
        return True
    except ValueError:
        return False

def is_number(s):
    try:
        n = float(s)
        if (math.isinf(n) or math.isnan(n)): return False
        return True
    except ValueError:
        return False

#http://code.activestate.com/recipes/363051/
#http://rightfootin.blogspot.com/2006/09/more-on-python-flatten.html
def flatten(l, ltypes=(list, tuple)):
    ltype = type(l)
    l = list(l)
    i = 0
    while i < len(l):
        while isinstance(l[i], ltypes):
            if not l[i]:
                l.pop(i)
                i -= 1
                break
            else:
                l[i:i + 1] = l[i]
        i += 1
    return ltype(l)

#http://kogs-www.informatik.uni-hamburg.de/~meine/python_tricks
def flattenX(x):
    """flatten(sequence) -> list

    Returns a single, flat list which contains all elements retrieved
    from the sequence and all recursively contained sub-sequences
    (iterables)."""

    result = []
    for el in x:
        #if isinstance(el, (list, tuple)):
        if hasattr(el, "__iter__") and not isinstance(el, basestring):
            result.extend(flatten(el))
        else:
            result.append(el)
    return result

def file_exists(file_path):
    if isinstance(file_path, basestring) and os.path.isfile(file_path):
        try:
            open(file_path)
            return True
        except IOError as e: return False
    else: return False

def is_exe(file_path):
    return file_exists(file_path) and os.access(file_path, os.X_OK)

def list_existing_files(files):
    """list_existing_files(sequence of paths) -> list of paths

    Returns a single, flat list of all paths that openable files
    or None if all files are missing"""

    result = []
    for f in flatten(files):
        if file_exists(f): result.append(f)
    return None if len(result) == 0 else result

def list_missing_files(files):
    """list_missing_files(sequence of paths) -> list of paths

    Returns a single, flat list of all paths that are not files or cannot be opened
    or None if all files are present and accessible"""

    result = []
    for f in flatten(files):
        if not file_exists(f): result.append(f)
    return None if len(result) == 0 else result

def print_log(s, log_file=None, print_this=True):
    t = strftime("%Y-%m-%d %H:%M:%S", localtime())
    if print_this:
        log_str = '{}\t{}'.format(t, (s if isinstance(s, basestring) else '\n'.join(s)).rstrip('\n'))
        print log_str
    if log_file:
        if isinstance(s, basestring):
            s = s.split('\n')
        for line in s:
            log_str = '{}\t{}'.format(t, line.rstrip('\n'))
            with open(log_file, 'a') as l:
                l.write(log_str+'\n')

def localtime_stamp():
    return '{}'.format(strftime("%Y-%m-%d %H:%M:%S", localtime()))

def localtime_squish():
    return '{}'.format(strftime("%Y%m%d%H%M%S", localtime()))

def print_err(s):
    sys.stderr.write('[{}] {}\n'.format(strftime("%Y-%m-%d %H:%M:%S", localtime()),s.rstrip('\n')))

def query_yes_no(question, default="yes"):
    """Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is one of "yes" or "no".
    """
    valid = {"yes": True, "y": True, "ye": True,
             "no": False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = raw_input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "
                             "(or 'y' or 'n').\n")
def rcut(string, cut):
    return re.sub()

def bam_strip(string):
    s = string[:-3] if string.endswith('.gz') else string
    return s[:-4] if s.endswith('.bam') or s.endswith('.sam') else s

def qseq_strip(string):
    s = string[:-3] if string.endswith('.gz') else string
    return s[:-9] if s.endswith('_qseq.txt') else s

def vcf_strip(string):
    s = string[:-3] if string.endswith('.gz') else string
    return s[:-4] if s.endswith('.vcf') else s

def r_strip(string, strip):
    return string[:-len(strip)] if string.endswith(strip) else string

def l_strip(string, strip):
    return string[len(strip):] if string.startswith(strip) else string

def swap_ext(string, old_ext, new_ext):
    return r_strip(string, old_ext)+new_ext

QseqTrio = collections.namedtuple('QseqTrio', 'read1 read2 read3')
FilegroupInfo = collections.namedtuple('FilegroupInfo', 'machine, run, lane')

def get_qseqs(dir, lane, read, gz=None):
    #return a list of qseqs
    qseqs = glob(os.path.join(dir,'*qseq*'))
    qseqs_tmp = set([q for q in qseqs if (lane == '*' or qseq_peek(q)['lane'] == int(lane)) and (read == '*' or qseq_peek(q)['read'] == int(read))])
    if gz == None:
        return sorted(qseqs_tmp)
    elif gz:
        return sorted(q for q in qseqs_tmp if is_gz(q))
    else:
        return sorted(q for q in qseqs_tmp if not is_gz(q))

def qseq_read_dict(qseq_line, path=None):
    (MACHINE,RUN,LANE,TILE,X,Y,INDEX,READ,SEQUENCE,QUALITY,FILTER) = range(11)
    info = {'machine': None, 'run': None, 'lane': None, 'tile': None, 'x': None, 'y': None, 'index': None, 'read': None, 'sequence': None, 'quality': None, 'filter': None}
    fields = qseq_line.rstrip('\n').split('\t')
    if ((len(fields) == FILTER+1) and (is_number(fields[LANE])) and (int(fields[LANE]) in range(1,9))
        and (is_number(fields[TILE])) and (is_number(fields[X])) and (is_number(fields[Y]))
        and (is_number(fields[INDEX])) and (is_number(fields[READ])) and (int(fields[READ]) in range(1,4))
        and (re.search(r"[^ACGTN.]", fields[SEQUENCE], re.I) == None) and (len(fields[SEQUENCE]) == len(fields[QUALITY]))
        and (is_number(fields[FILTER])) and (int(fields[FILTER]) in range(0,2))):
        info['machine'] = fields[MACHINE]
        info['run'] = fields[RUN]
        info['lane'] = int(fields[LANE])
        info['tile'] = int(fields[TILE])
        info['x'] = int(fields[X]) if '.' not in fields[X] else int(10*round(float(fields[X]))+1000)
        info['y'] = int(fields[Y]) if '.' not in fields[Y] else int(10*round(float(fields[Y]))+1000)
        info['index'] = int(fields[INDEX])
        info['read'] = int(fields[READ])
        info['sequence'] = fields[SEQUENCE]
        info['quality'] = fields[QUALITY]
        info['filter'] = int(fields[FILTER])
        info['read_length'] = len(fields[SEQUENCE])
        info['path'] = path
        return info
    else:
        return None

def qseq_read_valid(dict_or_qseq_line):
    info = qseq_read_dict(dict_or_qseq_line) if isinstance(dict_or_qseq_line, basestring) else dict_or_qseq_line #if not str, assume dict
    if info == None or info['machine'] == None:
        return False
    else:
        return True

def qseq_all_reads_valid(path_to_qseq_file):
    try:
        with open_gz_or_text(path_to_qseq_file) as f:
            for read in f:
                (MACHINE,RUN,LANE,TILE,X,Y,INDEX,READ,SEQUENCE,QUALITY,FILTER) = range(11)
                fields = read.rstrip('\n').split('\t')
                if not(((len(fields) == FILTER+1) and (is_number(fields[LANE])) and (int(fields[LANE]) in range(1,9))
                        and (is_number(fields[TILE])) and (is_number(fields[X])) and (is_number(fields[Y]))
                        and (is_number(fields[INDEX])) and (is_number(fields[READ])) and (int(fields[READ]) in range(1,4))
                        and (re.search(r"[^ACGTN.]", fields[SEQUENCE], re.I) == None) and (len(fields[SEQUENCE]) == len(fields[QUALITY]))
                        and (is_number(fields[FILTER])) and (int(fields[FILTER]) in range(0,2)))):
                    return False
    except:
        return False
    else:
        return True

def qseq_reads_match(read1, read2):
    #input may be string or qseq_line_dict
    p1 = qseq_read_dict(read1) if isinstance(read1, basestring) else read1
    p2 = qseq_read_dict(read2) if isinstance(read2, basestring) else read2
    if p1 == None or p2 == None:
        return False
    else:
        return p1['machine'] == p2['machine'] and p1['run'] == p2['run'] and p1['lane'] == p2['lane'] and p1['tile'] == p2['tile'] and p1['x'] == p2['x'] and p1['y'] == p2['y']

def extimate_read_length(qseq):
    pass

#@HWI-ST0860:79:C03YLACXX:5:1101:1395:2165 1:N:0:ACTTGA
#@HWI-ST0860:79:C03YLACXX:5:1101:1395:2165 1:N:0:
__sequence_identifier_rx = re.compile(r'^@(?P<instrument>[a-zA-Z0-9_-]+):(?P<run>[0-9.eE-]+):(?P<flowcell>[a-zA-Z0-9]+):(?P<lane>[0-9.eE-]+):(?P<tile>[0-9.eE-]+):(?P<x_pos>[0-9.eE-]+):(?P<y_pos>[0-9.eE-]+) (?P<read>[0-9.eE-]+):(?P<is_filtered>[nyNY]+):(?P<control_number>[0-9]+):(?P<index_sequence>[ACGT]+)?$', re.I)
__sequence_rx = re.compile(r'^[ACGTN.-]+$', re.I)
def fastq_peek(fastq):
    #returns dict with metadata (machine, run, lane, barcode, read) from first four lines of a fastq
    with open_gz_or_text(fastq) as file:
        sequence_identifier = file.readline().rstrip('\n')
        match = __sequence_identifier_rx.match(sequence_identifier)
        if not match:
            raise FastqError ('"{}" in file {} line 1 is not a proper fastq sequence identifier'.format(sequence_identifier, fastq))
        peek = match.groupdict('')
        peek['machine'] = peek['instrument']
        peek['barcode'] = peek['index_sequence']
        sequence = file.readline().rstrip('\n')
        if not __sequence_rx.match(sequence):
            raise FastqError ('"{}" in file {} line 2 is not a proper fastq sequence'.format(fastq, sequence))
        peek['sequence'] = sequence
        quality_score_identifier = file.readline().rstrip('\n')
        if not quality_score_identifier == '+':
            raise FastqError ('"{}" in file {} line 3 is not a proper fastq quality score identifier'.format(fastq, quality_score_identifier))
        quality_score = file.readline().rstrip('\n')
        if not len(quality_score) == len(sequence):
            raise FastqError ('"{}" in file {} line 4 is not a proper fastq quality score because it is not the same length as sequence "{}"'.format(fastq, quality_score, sequence))
        peek['quality_score'] = quality_score
    peek['path'] = fastq
    peek['compressed'] = is_gz(fastq)
    return peek

def get_fastq_readgroups(fastqs, metadata):
    #returns metadata with fastq_reads:paths and fastq_peek
    #from list of fastq globs and metadata file
    fastqs = unglob(fastqs)
    readgroups = get_readgroups_from_metadata(metadata)
    if len(fastqs) == 0 or len(readgroups) == 0:
        return None
    for p in [fastq_peek(f) for f in fastqs]:
        meta = readgroups.get((p['machine'],p['run'],p['lane'],p['barcode']))
        if not meta:
            raise FastqError ('fastq file "{}" not found in metadata'.format(p['path']))
        #fastq_reads is a dict with key (read,path) and value fastq_peek data from first fastq record in theache file
        meta.fastq_reads[(p['read'],p['path'])] = [p] if meta.fastq_reads.get(p['read']) == None else meta.fastq_reads[p['read']].append(p)
        #meta.fastq_reads[p['read']] = [p] if meta.fastq_reads.get(p['read']) == None else meta.fastq_reads[p['read']].append(p)
    return {r: readgroups[r] for r in readgroups if readgroups[r].fastq_reads}

def fastqs_hierarchy(fastqs):
    #returns {machines:runs:lanes:barcodes:reads:peek}
    #from list of fastq globs
    fastqs = unglob(fastqs)
    if len(fastqs) == 0:
        return None
    peeks = [fastq_peek(f) for f in fastqs]
    machines = {x['machine']: None for x in peeks}
    for machine in machines:
        machines[machine] = {x['run']: None for x in peeks if x['machine'] == machine}
        for run in machines[machine]:
            machines[machine][run] = {x['lane']: None for x in peeks if x['machine'] == machine and x['run'] == run}
            for lane in machines[machine][run]:
                machines[machine][run][lane] = {x['barcode']: None for x in peeks if x['machine'] == machine and x['run'] == run and x['lane'] == lane}
                for barcode in machines[machine][run][lane]:
                    machines[machine][run][lane][barcode] = {x['read']: x for x in peeks if x['machine'] == machine and x['run'] == run and x['lane'] == lane and x['barcode'] == barcode}
    return machines

def qseq_peek(qseq, barcode_re=None):
    #returns dict with metadata (machine, run, lane, read) from first line of a qseq and barcode from file path
    with open_gz_or_text(qseq) as file:
        read = file.readline()
    peek = qseq_read_dict(read, qseq)
    peek['compressed'] = is_gz(qseq)
    match = re.search(barcode_re, peek['path']) if barcode_re else None
    peek['barcode'] = match.group(1) if match else ''
    return peek

def qseqs_hierarchy(qseqs, barcode_re=None):
    #returns {machines:runs:lanes:barcodes:reads:peek}
    #from list of qseq globs and barcode regular expression
    qseqs = unglob(qseqs)
    if len(qseqs) == 0:
        return None
    peeks = [qseq_peek(q, barcode_re) for q in qseqs]
    machines = {x['machine']: None for x in peeks}
    for machine in machines:
        machines[machine] = {x['run']: None for x in peeks if x['machine'] == machine}
        for run in machines[machine]:
            machines[machine][run] = {x['lane']: None for x in peeks if x['machine'] == machine and x['run'] == run}
            for lane in machines[machine][run]:
                machines[machine][run][lane] = {x['barcode']: None for x in peeks if x['machine'] == machine and x['run'] == run and x['lane'] == lane}
                for barcode in machines[machine][run][lane]:
                    machines[machine][run][lane][barcode] = {x['read']: x for x in peeks if x['machine'] == machine and x['run'] == run and x['lane'] == lane and x['barcode'] == barcode}
    return machines

def qseq_tiles_hierarchy(qseqs):
    #returns {machines:runs:lanes:reads:tiles:peek}
    #from list of qseq globs and barcode regular expression
    qseqs = unglob(qseqs)
    if len(qseqs) == 0:
        return None
    peeks = [qseq_peek(q) for q in qseqs]
    machines = {x['machine']: None for x in peeks}
    for machine in machines:
        machines[machine] = {x['run']: None for x in peeks if x['machine'] == machine}
        for run in machines[machine]:
            machines[machine][run] = {x['lane']: None for x in peeks if x['machine'] == machine and x['run'] == run}
            for lane in machines[machine][run]:
                machines[machine][run][lane] = {x['read']: None for x in peeks if x['machine'] == machine and x['run'] == run and x['lane'] == lane}
                for read in machines[machine][run][lane]:
                    machines[machine][run][lane][read] = {x['tile']: x for x in peeks if x['machine'] == machine and x['run'] == run and x['lane'] == lane and x['read'] == read}
    return machines

def qseq_tiles_reads_hierarchy(qseqs):
    #returns {machines:runs:lanes:tiles:reads:peek}
    #from list of qseq globs and barcode regular expression
    qseqs = unglob(qseqs)
    if len(qseqs) == 0:
        return None
    peeks = [qseq_peek(q) for q in qseqs]
    machines = {x['machine']: None for x in peeks}
    for machine in machines:
        machines[machine] = {x['run']: None for x in peeks if x['machine'] == machine}
        for run in machines[machine]:
            machines[machine][run] = {x['lane']: None for x in peeks if x['machine'] == machine and x['run'] == run}
            for lane in machines[machine][run]:
                machines[machine][run][lane] = {x['tile']: None for x in peeks if x['machine'] == machine and x['run'] == run and x['lane'] == lane}
                for tile in machines[machine][run][lane]:
                    machines[machine][run][lane][tile] = {x['read']: x for x in peeks if x['machine'] == machine and x['run'] == run and x['lane'] == lane and x['tile'] == tile}
    return machines

def qseq_read1_matches(q1,q2):
    #checks for match of first record
    return qseq_reads_match(qseq_peek(q1), qseq_peek(q2))

def qseq_all_reads_match(q1,q2,q3=None):
    #checks for matches of all records
    f1 = open_gz_or_text(q1)
    f2 = open_gz_or_text(q2)
    f3 = open_gz_or_text(q2) if q3 else None
    l1 = f1.readline()
    l2 = f2.readline()
    if q3:
        l3 = f3.readline()
    else:
        l3 = False
    if not (bool(l1) or bool(l2) or bool(l3)):
        f1.close()
        f2.close()
        if q3:
            f3.close
        return True
    if (bool(l1) != bool(l2)) or (q3 and bool(l1) != bool(l3)):
        #files with different number of records
        return False
    if not qseq_reads_match(l1,l2):
        return False
    elif q3 and not qseq_reads_match(l1,l3):
        return False

def qseqs2qseqtrios(qseqs):
    #check that the first reads in each file match,
    #non-duplicate paired end and barcode read paths, typically (pe1, barcode, pe2)
    #returns three lists of sorted, matchedfiles, some possibly empty
    qseqs1 = sorted(set([q for q in qseqs if qseq_peek(q)['read'] == 1]))
    qseqs2 = sorted(set([q for q in qseqs if qseq_peek(q)['read'] == 2]))
    qseqs3 = sorted(set([q for q in qseqs if qseq_peek(q)['read'] == 3]))
    if not qseqs1:
        raise QseqError('at line {} no files for read 1 in {}'.format(line_count, qseqs))
    if len(qseqs1) == len(qseqs2) and len(qseqs1) == len(qseqs3):
        #paried end reads with barcode
        for i in range(len(qseqs1)):
            if not qseq_read1_matches(qseqs1[i], qseqs3[i]) or not qseq_read1_matches(qseqs1[i], qseqs2[i]):
                raise QseqError('at line {} reads do not match in [{}]  [{}] and [{}]'.format(line_count, qseqs1[i], qseqs2[i], qseqs3[i]))
    elif len(qseqs1) == len(qseqs2) and not qseqs3:
        #paired end reads without barcode
        for i in range(len(qseqs1)):
            if not qseq_read1_matches(qseqs1[i], qseqs2[i]):
                raise QseqError('at line {} reads do not match in [{}] and [{}]'.format(line_count, qseqs1[i], qseqs2[i]))
    else:
        raise QseqError('at line {} number of files differs for paired end reads in {}'.format(line_count, qseqs))
    return (qseqs1, qseqs2, qseqs3)

def bam_peek(path_to_sam_or_bam_file):
    #returns dict of header; h['RG'] will be one or more dicts of RG keys and values
    with open_bam(path_to_sam_or_bam_file) as samfile:
        h = samfile.header
    return h

def get_metadata(path_to_metadata, primary_key=('machine','run','lane','barcode')):
    if not file_exists(path_to_metadata):
        raise GetMetadataError('metadata file {} does not exist'.format(path_to_metadata))
    null_key = tuple(['' for m in range(len(primary_key))])
    #get rows from metadata file as a dictionary with keys from first row column headers
    #ignore blank primary keys
    md = [m for m in csv.DictReader(open(path_to_metadata,'rb'), dialect=csv.excel_tab)
          if tuple([m[k] for k in primary_key]) != null_key]
    #create dictionary of metadata with primary key
    metadata = {tuple([m[k] for k in primary_key]): m for m in md}
    #check for duplicate primary keys
    if len(metadata) != len(md):
        raise GetMetadataError('primary key (column(s) {}) duplicated in metadata'.format(','.join(primary_key)))
    else:
        return metadata

def get_readgroups_from_metadata(path_to_metadata):
    #get metadata file to map machine, run, lane, barcode to readgroup info
    metadata = get_metadata(path_to_metadata)
    if not metadata:
        raise Qseqs2BamsError('no records selected from metadata file {}'.format(args.metadata))
    readgroups = ReadGroups()
    for k in metadata.keys():
        m = metadata[k]
        readgroups[k] = ReadGroup(
            do_not_align=m.get('do_not_align', None),
            sample=m.get('sample', None),
            original_sample=m.get('original_', None),
            library=m.get('library', None),
            readgroup_id=m.get('readgroup_id', None),
            machine=m.get('machine', None),
            run=m.get('run', None),
            lane=m.get('lane', None),
            sample_sheet_barcode=m.get('sample_sheet_barcode', None),
            barcode=m.get('barcode', None),
            owner=m.get('owner', None),
            email=m.get('email', None),
            collaborator=m.get('collaborator', None),
            collaborator_email=m.get('collaborator_email', None),
            study=m.get('study', None),
            cohort=m.get('cohort', None),
            family=m.get('family', None),
            individual=m.get('individual', None),
            father=m.get('father', None),
            mother=m.get('mother', None),
            sex=m.get('sex', None),
            tumor_status=m.get('tumor_status', None),
            species=m.get('species', None),
            multiplexed_library=m.get('multiplexed_library', None),
            library_protocol=m.get('library_protocol', None),
            bioanalyzer_mean_bp=m.get('bioanalyzer_mean_bp', None),
            predicted_median_insert_size=m.get('predicted_median_insert_size', None),
            median_insert_size=m.get('median_insert_size', None),
            capture_bait_interval_list=m.get('capture_bait_interval_list', None),
            genomic_target_interval_list=m.get('genomic_target_interval_list', None),
            barcode_tag_list=m.get('barcode_tag_list', None),
            sequencing_library=m.get('sequencing_library', None),
            adapter_1=m.get('adapter_1', None),
            adapter_2=m.get('adapter_2', None),
            experiment_protocol=m.get('experiment_protocol', None),
            experiment=m.get('experiment', None),
            sequencing_center=m.get('sequencing_center', None),
            platform=m.get('platform', None),
            instrument_model=m.get('instrument_model', None),
            run_date=m.get('run_date', None),
            run_folder=m.get('run_folder', None),
            flowcell=m.get('flowcell', None),
            index1=m.get('index1', None),
            index2=m.get('index2', None),
            platform_unit=m.get('platform_unit', None),
            read_1_length=m.get('read_1_length', None),
            read_2_length=m.get('read_2_length', None),
            readgroup_description=m.get('readgroup_description', None),
            additional_readgroup_description=m.get('additional_readgroup_description', None),
            run_top_dir=m.get('run_top_dir', None),
            run_dir=m.get('run_dir', None),
            fastq_top_dir=m.get('fastq_top_dir', None),
            fastq_glob=m.get('fastq_glob', None),
        )
    return readgroups

def default_parser():
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--top_dir', '-d', #required=True,
                        help='top-level directory for all analysis')
    parser.add_argument('--parent_job_dir',
                        help='top-level directory for all analysis')
    parser.add_argument('--email', nargs='*',
                        help='email address for job completion notifications')
    parser.add_argument('--synchronous', action='store_true', default=False,
                        help='run job synchronously')
    parser.add_argument('--hold_jid', nargs='*',
                        help='space-separated list of job ids that must complete before this job runs')
    parser.add_argument('--config',
                        help='additional (overriding?) configuration file (default=[pypeline.cfg, ~/pypeline.cfg])')
    return parser

def create_default_directory_structure(top_dir):
    dirs = {}
    dirs['top'] = top_dir
    makedir(dirs['top'])
    dirs['jobs'] = os.path.join(dirs['top'],'jobs')
    makedir(dirs['jobs'])
    dirs['gatk_jobs'] = os.path.join(dirs['jobs'],'gatk_jobs')
    makedir(dirs['gatk_jobs'])
    dirs['tmp'] = os.path.join(dirs['top'],'tmp')
    makedir(dirs['tmp'])
    #dirs['qseqs'] = os.path.join(dirs['top'],'qseqs')
    #makedir(dirs['qseqs'])
    dirs['fastqs'] = os.path.join(dirs['top'],'fastqs')
    makedir(dirs['fastqs'])
    dirs['linkdatagen'] = os.path.join(dirs['top'],'linkdatagen')
    makedir(dirs['linkdatagen'])
    dirs['plink'] = os.path.join(dirs['top'],'plink')
    makedir(dirs['plink'])
    dirs['bams'] = os.path.join(dirs['top'],'bams')
    makedir(dirs['bams'])
    dirs['readgroup_bam_parts'] = os.path.join(dirs['bams'],'readgroup_bam_parts')
    makedir(dirs['readgroup_bam_parts'])
    dirs['readgroup_bams'] = os.path.join(dirs['bams'],'readgroup_bams')
    makedir(dirs['readgroup_bams'])
    dirs['library_bams'] = os.path.join(dirs['bams'],'library_bams')
    makedir(dirs['library_bams'])
    dirs['sample_bams'] = os.path.join(dirs['bams'],'sample_bams')
    makedir(dirs['sample_bams'])
    dirs['study_bams'] = os.path.join(dirs['bams'],'study_bams')
    makedir(dirs['study_bams'])
    dirs['vcfs'] = os.path.join(dirs['top'],'vcfs')
    makedir(dirs['vcfs'])
    dirs['mpileups'] = os.path.join(dirs['vcfs'],'mpileups')
    makedir(dirs['mpileups'])
    dirs['gatk'] = os.path.join(dirs['top'],'gatk')
    makedir(dirs['gatk'])
    dirs['metrics'] = os.path.join(dirs['top'],'metrics')
    makedir(dirs['metrics'])
    dirs['fastq_metrics'] = os.path.join(dirs['metrics'],'fastq_metrics')
    makedir(dirs['fastq_metrics'])
    dirs['novoalign_metrics'] = os.path.join(dirs['metrics'],'novoalign_metrics')
    makedir(dirs['novoalign_metrics'])
    dirs['picard_bam_metrics'] = os.path.join(dirs['metrics'],'picard_bam_metrics')
    makedir(dirs['picard_bam_metrics'])
    dirs['readgroup_metrics'] = os.path.join(dirs['picard_bam_metrics'],'readgroup_metrics')
    makedir(dirs['readgroup_metrics'])
    dirs['library_metrics'] = os.path.join(dirs['picard_bam_metrics'],'library_metrics')
    makedir(dirs['library_metrics'])
    dirs['library_pre-markdup_metrics'] = os.path.join(dirs['library_metrics'],'library_pre-markdup_metrics')
    makedir(dirs['library_pre-markdup_metrics'])
    dirs['library_duplication_metrics'] = os.path.join(dirs['library_metrics'],'library_duplication_metrics')
    makedir(dirs['library_duplication_metrics'])
    dirs['library_post-markdup_metrics'] = os.path.join(dirs['library_metrics'],'library_post-markdup_metrics')
    makedir(dirs['library_post-markdup_metrics'])
    dirs['sample_metrics'] = os.path.join(dirs['picard_bam_metrics'],'sample_metrics')
    makedir(dirs['sample_metrics'])
    dirs['depth_of_coverage'] = os.path.join(dirs['metrics'],'depth_of_coverage')
    makedir(dirs['depth_of_coverage'])
    dirs['quality_score_recalibration_metrics'] = os.path.join(dirs['metrics'],'quality_score_recalibration_metrics')
    makedir(dirs['quality_score_recalibration_metrics'])
    dirs['unified_genotyper_metrics'] = os.path.join(dirs['metrics'],'unified_genotyper_metrics')
    makedir(dirs['unified_genotyper_metrics'])
    dirs['haplotype_caller_metrics'] = os.path.join(dirs['metrics'],'haplotype_caller_metrics')
    makedir(dirs['haplotype_caller_metrics'])
    dirs['vcf_metrics'] = os.path.join(dirs['metrics'],'vcf_metrics')
    makedir(dirs['vcf_metrics'])
    return dirs

def get_config(args=None, config_file='pypeline.cfg'):
    config = SafeConfigParser()
    config.readfp(open(os.path.join(os.path.dirname(__file__), config_file)))
    cfg_files = [os.path.expanduser('~/.'+config_file)]
    if args and 'config' in dir(args) and args.config:
        cfg_files.append(args.config)
    config.read(cfg_files)
    return config

def run_job(commands, jobName, job_dir, qout_dir=None, email=None,
            synchronous=False, processors=None, memory='6G', queue=None, hold_jid=None):
    j = job.Job(jobName=jobName, outputPath=qout_dir if qout_dir else job_dir,
                errorPath=qout_dir if qout_dir else job_dir, workingDirectory=job_dir,
                email=email, blockEmail=not bool(email), processors=processors,
                memory=memory, queue=queue, hold_jid=hold_jid, exclusive=True)
    try:
        j.executeCommands(commands, synchronous=synchronous)
    except Exception as e:
        raise RunJobError('There was an error running job {}\n{}\nintermediate files are in {}\n'.format(jobName, e, job_dir))
    else:
        with open(os.path.join(job_dir,'{}_job_status_report.txt'.format(jobName)), 'w') as s:
            s.write(format('\n'.join(j.completionMsg)))
        return j

def check_files(files_exist=None,files_contain=None,files_not_before=None):
    """files_exist = arrray of files
       files_contain = array of tuples (file,expected)
       files_not_before = array of tuples (file1,file2)
       returns True if all files in files_exist exist
       and all files in files_contain[0] exist and have contents equal to files_contain[1]
       and all files in files_not_before exist and m_time of files_not_before[1] >= m_time of files_not_before[0]"""
    if isinstance(files_exist, basestring):
        files_exist = [files_exist]
    if files_exist and not isinstance(files_exist, (list, tuple, set)):
        raise CheckFilesError("files_exist must be list, tuple, or set")
    if files_contain:
        if isinstance(files_contain, (list, tuple, set)):
            for fc in files_contain:
                if not isinstance(fc,(list, tuple)) or len(fc)<2:
                    raise CheckFilesError("files_contain must contain list or tuple")
        else:
            raise CheckFilesError("files_contain must be list, tuple, or set")
    if files_not_before:
        if isinstance(files_not_before, (list, tuple, set)):
            for fnb in files_not_before:
                if not isinstance(fnb,(list, tuple)) or len(fnb)<2:
                    raise CheckFilesError("files_not_before must contain list or tuple")
        else:
            raise CheckFilesError("files_not_before must be list, tuple, or set")
    try:
        if files_exist:
            for file in files_exist:
                if not file_exists(file):
                    return False
        if files_contain:
            for file,expected in files_contain:
                if not file_exists(file):
                    return False
                with open(file) as f:
                    f_contents = f.read()
                    if f_contents != expected:
                        return False
        if files_not_before:
            for file1,file2 in files_not_before:
                if not file_exists(file1) or not file_exists(file2):
                    return False
                if os.path.getmtime(file2) < os.path.getmtime(file1):
                    return False
    except:
        return False
    return True

#sql columns
def get_sql_column_spec(columns):
    return {c: {'type': None, 'size_min': None, 'size_max': None, 'min': None, 'max': None} for c in columns}

def get_sql_data_dict(columns, values):
    return {columns[i]: values[i] for i in range(len(values))}

def update_sql_column_spec(sql_column_spec, sql_data_dict):
    for c in sql_data_dict.keys():
        this_value = sql_data_dict[c]
        spec = sql_column_spec.get(c)
        if spec and this_value != '':
            spec['size_min'] = len(this_value) if spec['size_min'] == None else min(len(this_value),spec['size_min'])
            spec['size_max'] = len(this_value) if spec['size_max'] == None else max(len(this_value),spec['size_max'])
            if is_number(this_value):
                spec['min'] = float(this_value) if spec['min'] == None else min(float(this_value),spec['min'])
                spec['max'] = float(this_value) if spec['max'] == None else max(float(this_value),spec['max'])
            if not spec['type'] == 'varchar':
                if not is_number(this_value):
                    spec['type'] = 'varchar'
                elif not spec['type'] == 'float':
                    if not is_int(this_value):
                        spec['type'] = 'float'
                    else:
                        spec['type'] = 'int'
#microsoft sql server
def get_microsoft_sql_server_scripts(
    data_file=None,
    database='database',
    schema='schema',
    table_name='table',
    index_base_name=None,
    primary_key=[],
    clustered_index=[],
    unique_indexes=[],
    indexes=[],
    columns_out=[],
    columns_out_spec={},
    rows_to_delete=0
    ):
    #returns [create_table_script, import_data_script, index_scripts]

    if not data_file:
        data_file = r'C:\path\to\{}.txt'.format(r_strip(table_name, '.txt'))

    if data_file[1] == ':':
        nix_data_file = '~'+data_file.replace('\\','/')[2:]
        windows_data_file = data_file
    else:
        nix_data_file = data_file
        windows_data_file = 'C:'+('' if data_file[0] in ('\\','/') else '\\')+data_file.replace('/','\\')

    table_name = r_strip(table_name, '.txt').replace('.','_')[:128]

    if not index_base_name:
        index_base_name = table_name

    if primary_key:
        if not isinstance(primary_key, (list, tuple, set)):
            primary_key = primary_key.split(',')
        elif isinstance(primary_key, (list, tuple, set)) and len(primary_key) == 1 and ',' in primary_key[0]:
            primary_key = primary_key[0].split(',')
    else:
        primary_key = []

    if clustered_index:
        if not isinstance(clustered_index, (list, tuple, set)):
            clustered_indexi = clustered_index.split(',')
        elif isinstance(clustered_index, (list, tuple, set)) and len(clustered_index) == 1 and ',' in clustered_index[0]:
            clustered_index = clustered_index[0].split(',')
    else:
        clustered_index = []

    if unique_indexes:
        for i in range(len(unique_indexes)):
            if not isinstance(unique_indexes[i], (list, tuple, set)):
                unique_indexes[i] = unique_indexes[i].split(',')    
    else:
        unique_indexes = []

    if indexes:
        for i in range(len(indexes)):
            if not isinstance(indexes[i], (list, tuple, set)):
                indexes[i] = indexes[i].split(',')    
    else:
        indexes = []

    #CREATE TABLE expression
    create_table_script = """SET ANSI_NULLS ON
GO
SET QUOTED_IDENTIFIER ON
GO
SET ANSI_PADDING ON
GO

--if necessary, substitute your [database]  for [{0}]
USE [{0}];

--substitute your [schema].[table_name] for [{1}].[{2}]
IF OBJECT_ID(N'[{0}].[{1}].[{2}]', N'U') IS NOT NULL 
DROP TABLE [{0}].[{1}].[{2}];
GO
CREATE TABLE [{0}].[{1}].[{2}] (
""".format(database, schema, table_name)

    #specify columns
    for c in columns_out:
        spec = columns_out_spec.get(c)
        if spec:
            if c == columns_out[0] and c.startswith('#'):
                c = c[1:]
            if spec['type'] == 'varchar':
                if spec['size_min'] == spec['size_max'] and spec['size_max']<=1000:
                    ms_type = '[CHAR] ({})'.format(spec['size_max'])
                else:
                    ms_type = '[VARCHAR] ({})'.format(spec['size_max'] if spec['size_max']<=1000 else 'MAX')
            elif spec['type'] == 'float':
                if spec['min'] >= -3.4e+38 and spec['max'] <= 3.4e+38:
                    ms_type = '[REAL]'
                elif spec['min'] >= -1.79e+308 and spec['max'] <= 1.79e+308:
                    ms_type = '[FLOAT]'
                else:
                    ms_type = '[VARCHAR] ({})'.format(spec['size_max'] if spec['size_max']<=1000 else 'MAX')
            elif spec['type'] == 'int':
                if spec['min'] >= 0 and spec['max'] <= 1:
                    ms_type = '[BIT]'
                elif spec['min'] >= 0 and spec['max'] <= 255:
                    ms_type = '[TINYINT]'
                elif spec['min'] >= -32768 and spec['max'] <= 32767:
                    ms_type = '[SMALLINT]'
                elif spec['min'] >= -2147483648 and spec['max'] <= 2147483647:
                    ms_type = '[INT]'
                elif spec['min'] >= -9223372036854775808 and spec['max'] <= 9223372036854775807:
                    ms_type = '[BIGINT]'
                else:
                    ms_type = '[VARCHAR] ({})'.format(spec['size_max'] if spec['size_max']<=1000 else 'MAX')
            else:
                if c in ('compound_het'):
                    ms_type = '[BIT]'
                else:
                    ms_type = '[CHAR] (1)'
            ms_spec = '\t[{0}] {1} {2}NULL{3}\n'.format(c, ms_type, 'NOT ' if c in primary_key or c in clustered_index else '', ',' if c != columns_out[-1] else '')

            #add column spec to CREATE TABLE expression 
            create_table_script += ms_spec

    #used to assure unique index names
    index_names = []

    #add primary key, if any
    if primary_key:
        col_part = '_'.join(primary_key)[:128]
        prefix_part = 'PK_'+index_base_name+'_'[:128-len(col_part)]
        ix_name = prefix_part + col_part
        suffix_part = 2
        while ix_name in index_names:
            ix_name = ix_name[:128-len(str(suffix_part))] + str(suffix_part)
            suffix_part += 1
        index_names.append(ix_name)
        ix_cols = ',\n\t'.join(['['+c+'] ASC' for c in primary_key])
        clustered = '' if clustered_index else 'CLUSTERED'

        #add primary key to CREATE TABLE expression
        create_table_script += """, CONSTRAINT [{0}] PRIMARY KEY {1}
(
	{2}
)WITH (PAD_INDEX  = OFF, STATISTICS_NORECOMPUTE  = OFF, IGNORE_DUP_KEY = OFF, ALLOW_ROW_LOCKS  = ON, ALLOW_PAGE_LOCKS  = ON) ON [PRIMARY]
""".format(ix_name, clustered, ix_cols)

    #end of CREATE TABLE expression
    create_table_script += ') ON [PRIMARY]\nGO\n\nSET ANSI_PADDING OFF\nGO\n'

    #windows import data script
    import_data_script = """
--if necessary, substitute your [database].[schema].[table_name] for [{0}].[{1}].[{2}]
--if necessary, substitute your path to unzipped file for {3}
--delete first {4} rows of comments before importing

TRUNCATE TABLE [{0}].[{1}].[{2}];

DECLARE @bulk_cmd varchar(MAX)
SET @bulk_cmd = 'BULK INSERT [{0}].[{1}].[{2}]
FROM ''{3}'' 
WITH (ROWTERMINATOR = '''+CHAR(10)+''',
FIRSTROW=2)'
EXEC(@bulk_cmd)
GO
""".format(database, schema, table_name, windows_data_file, rows_to_delete)

    #no *nix import script until we get the ability for *nix to communicate with ms sql

    #create index scripts, if any
    index_scripts = ''

    #clustered_index is a list of columns
    if clustered_index:
        col_part = '_'.join(clustered_index)[:128]
        prefix_part = 'IX_'+index_base_name+'_'[:128-len(col_part)]
        ix_name = prefix_part + col_part
        suffix_part = 2
        while ix_name in index_names:
            ix_name = ix_name[:128-len(str(suffix_part))] + str(suffix_part)
            suffix_part += 1
        index_names.append(ix_name)
        ix_cols = ',\n\t'.join(['['+c+'] ASC' for c in clustered_index])
        index_scripts += """
IF  EXISTS (SELECT * FROM sys.indexes WHERE object_id = OBJECT_ID(N'[{2}].[{3}]') AND name = N'{0}')
DROP INDEX [{0}] ON [{2}].[{3}] WITH ( ONLINE = OFF )
GO
CREATE CLUSTERED INDEX [{0}] ON [{1}].[{2}].[{3}]
(
    {4}
)WITH (PAD_INDEX  = OFF, STATISTICS_NORECOMPUTE  = OFF, SORT_IN_TEMPDB = OFF, IGNORE_DUP_KEY = OFF, ONLINE = OFF, ALLOW_ROW_LOCKS  = ON, ALLOW_PAGE_LOCKS  = ON) ON [PRIMARY]
GO
""".format(ix_name, database, schema, table_name, ix_cols)

    #unique_indexes is a list of unique indexes containing lists of columns
    for i in unique_indexes:
        col_part = '_'.join(i)[:128]
        prefix_part = 'IX_'+index_base_name+'_'[:128-len(col_part)]
        ix_name = prefix_part + col_part
        suffix_part = 2
        while ix_name in index_names:
            ix_name = ix_name[:128-len(str(suffix_part))] + str(suffix_part)
            suffix_part += 1
        index_names.append(ix_name)
        ix_cols = ',\n\t'.join(['['+c+'] ASC' for c in i])
        index_scripts += """
IF  EXISTS (SELECT * FROM sys.indexes WHERE object_id = OBJECT_ID(N'[{2}].[{3}]') AND name = N'{0}')
DROP INDEX [{0}] ON [{2}].[{3}] WITH ( ONLINE = OFF )
GO
CREATE UNIQUE NONCLUSTERED INDEX [{0}] ON [{1}].[{2}].[{3}]
(
    {4}
)WITH (PAD_INDEX  = OFF, STATISTICS_NORECOMPUTE  = OFF, SORT_IN_TEMPDB = OFF, IGNORE_DUP_KEY = OFF, ONLINE = OFF, ALLOW_ROW_LOCKS  = ON, ALLOW_PAGE_LOCKS  = ON) ON [PRIMARY]
GO
""".format(ix_name, database, schema, table_name, ix_cols)

    #indexes is a list of nonclustered indexes containing lists of columns
    for i in indexes:
        col_part = '_'.join(i)[:128]
        prefix_part = 'IX_'+index_base_name+'_'[:128-len(col_part)]
        ix_name = prefix_part + col_part
        suffix_part = 2
        while ix_name in index_names:
            ix_name = ix_name[:128-len(str(suffix_part))] + str(suffix_part)
            suffix_part += 1
        index_names.append(ix_name)
        ix_cols = ',\n\t'.join(['['+c+'] ASC' for c in i])
        index_scripts += """
IF  EXISTS (SELECT * FROM sys.indexes WHERE object_id = OBJECT_ID(N'[{2}].[{3}]') AND name = N'{0}')
DROP INDEX [{0}] ON [{2}].[{3}] WITH ( ONLINE = OFF )
GO
CREATE NONCLUSTERED INDEX [{0}] ON [{1}].[{2}].[{3}]
(
    {4}
)WITH (PAD_INDEX  = OFF, STATISTICS_NORECOMPUTE  = OFF, SORT_IN_TEMPDB = OFF, IGNORE_DUP_KEY = OFF, ONLINE = OFF, ALLOW_ROW_LOCKS  = ON, ALLOW_PAGE_LOCKS  = ON) ON [PRIMARY]
GO
""".format(ix_name, database, schema, table_name, ix_cols)

    return [create_table_script, import_data_script, index_scripts]

def get_mysql_scripts(
    data_file=None,
    database='database',
    schema='schema',
    table_name='table',
    index_base_name=None,
    primary_key=[],
    clustered_index=[],
    unique_indexes=[],
    indexes=[],
    columns_out=[],
    columns_out_spec={},
    rows_to_delete=0,
    include_mysql_schema=False,
    engine='MyISAM',
    ):
    #returns [create_table_script, import_data_script_unix, import_data_script_windows]

    if not data_file:
        data_file = '\path\to\{}.txt'.format(r_strip(table_name, '.txt'))

    if data_file[1] == ':':
        nix_data_file = '~'+data_file.replace('\\','/')[2:]
        windows_data_file = data_file
    else:
        nix_data_file = data_file
        windows_data_file = 'C:'+('' if data_file[0] in ('\\','/') else '\\')+data_file.replace('/','\\')

    forbidden = re.compile(r'[^0-9,a-z,A-Z$_]')
    table_name = forbidden.sub('_',r_strip(table_name, '.txt'))[:64]

    if not index_base_name:
        index_base_name = table_name

    if primary_key:
        if not isinstance(primary_key, (list, tuple, set)):
            primary_key = primary_key.split(',')
        elif isinstance(primary_key, (list, tuple, set)) and len(primary_key) == 1 and ',' in primary_key[0]:
            primary_key = primary_key[0].split(',')
        primary_key = [forbidden.sub('_', n)[:64] for n in  primary_key]
    else:
        primary_key = []

    #MyISAM doesn't support clustered indexes
    #TODO: support for non-MyISAM clustered indexes; for now treat as regular indexes
    if clustered_index:
        if indexes == None:
            indexes = [clustered_index]
        else:
            indexes.append(clustered_index)

    if unique_indexes:
        for i in range(len(unique_indexes)):
            if not isinstance(unique_indexes[i], (list, tuple, set)):
                unique_indexes[i] = unique_indexes[i].split(',')
                for j in range(len(unique_indexes[i])): 
                    unique_indexes[i][j] = forbidden.sub('_', unique_indexes[i][j])[:64]
    else:
        unique_indexes = []

    if indexes:
        for i in range(len(indexes)):
            if not isinstance(indexes[i], (list, tuple, set)):
                indexes[i] = indexes[i].split(',')
                for j in range(len(indexes[i])): 
                    indexes[i][j] = forbidden.sub('_', indexes[i][j])[:64]
    else:
        indexes = []

    #CREATE TABLE expression
    create_table_script = """# *nix CREATE TABLE usage:
#mysql -h <host> -P <port> -u <user> -p "{0}" < {2}.mysql;

# Windows CREATE TABLE usage"
#C:\\Program Files\\MySQL\\MySQL Server 5.5\\bin\\mysql -h host -P port --u user -p "{0}" < {3}.mysql;

#if executing the CREATE TABLE statement within MySQL, uncomment the next line:
#USE `{0}`;

delimiter $$
DROP TABLE IF EXISTS `{1}`$$
CREATE TABLE `{1}` (
""".format(schema if schema else database, table_name, nix_data_file, windows_data_file)

    #specify columns
    table_specs = []
    for c in columns_out:
        spec = columns_out_spec.get(c)
        if spec:
            if c == columns_out[0] and c.startswith('#'):
                c = c[1:]
            if spec['type'] == 'varchar':
                if spec['size_max']<=500:  #row size limit of 65535 bytes, therefore text preferred over varchar, choosing arbitrary limit
                    my_type = 'varchar({})'.format(spec['size_max'])
                else:
                    my_type = 'text'
            elif spec['type'] == 'float':
                my_type = 'float'
            elif spec['type'] == 'int':
                if spec['min'] >= 0 and spec['max'] <= 1:
                    my_type = 'tinyint(1)'
                elif spec['min'] >= 0 and spec['max'] <= 255:
                    my_type = 'tinyint unsigned'
                elif spec['min'] >= -32768 and spec['max'] <= 32767:
                    my_type = 'smallint'
                elif spec['min'] >= -8388608 and spec['max'] <= 8388607:
                    my_type = 'mediumint'
                elif spec['min'] >= -2147483648 and spec['max'] <= 2147483647:
                    my_type = 'int'
                elif spec['min'] >= -9223372036854775808 and spec['max'] <= 9223372036854775807:
                    my_type = 'bigint'
                elif spec['size_max']<=500:
                    my_type = 'varchar({})'.format(spec['size_max'])
                else:
                    my_type = 'text'
            else:
                if c in ('compound_het'):
                    my_type = 'tinyint(1)'
                else:
                    my_type = 'char(1)'
            my_spec = '\t`{}` {} NULL'.format(forbidden.sub('_',c)[:64], my_type)
            table_specs.append(my_spec)

    #add column spec to CREATE TABLE expression 
    create_table_script += ',\n'.join(table_specs)

    keys = []
    if primary_key or indexes:
        create_table_script += ','

    #used to assure unique index names
    index_names = []

    #add primary key, if any
    if primary_key:
        #ix_name = forbidden.sub('_','PK_'+index_base_name+'_'+'_'.join(primary_key))[:64]
        ix_cols = '`,`'.join(primary_key)
        keys.append("""
PRIMARY KEY (`{0}`)""".format(ix_cols))

    #unique_indexes is a list of indexes containing lists of columns
    if unique_indexes:
        for i in unique_indexes:
            col_part = forbidden.sub('_','_'.join(i))[:64]
            prefix_part = forbidden.sub('_','IX_'+index_base_name+'_')[:64-len(col_part)]
            ix_name = prefix_part + col_part
            suffix_part = 2
            while ix_name in index_names:
                ix_name = ix_name[:64-len(str(suffix_part))] + str(suffix_part)
                suffix_part += 1
            index_names.append(ix_name)
            ix_cols = '`,`'.join(i)
            keys.append("""
UNIQUE KEY `{0}` (`{1}`)""".format(ix_name, ix_cols))

    #indexes is a list of indexes containing lists of columns
    if indexes:
        for i in indexes:
            col_part = forbidden.sub('_','_'.join(i))[:64]
            prefix_part = forbidden.sub('_','IX_'+index_base_name+'_')[:64-len(col_part)]
            ix_name = prefix_part + col_part
            suffix_part = 2
            while ix_name in index_names:
                ix_name = ix_name[:64-len(str(suffix_part))] + str(suffix_part)
                suffix_part += 1
            index_names.append(ix_name)
            ix_cols = '`,`'.join(i)
            keys.append("""
KEY `{0}` (`{1}`)""".format(ix_name, ix_cols))

    create_table_script += ','.join(keys)

    #end of CREATE TABLE expression
    create_table_script += """
) ENGINE={0} DEFAULT CHARSET=latin1$$

DELIMITER ;""".format(engine)

    #mysqlimport
    import_data_script_unix = """# *nix import:
#mysqlimport --host <host> --port <port> --user <user> -p --delete --local --verbose --ignore-lines={0} {1} "{2}"
""".format(rows_to_delete, schema if schema else database, nix_data_file)

    import_data_script_windows = """# Windows import:
#"C:\\Program Files\\MySQL\\MySQL Server 5.5\\bin\\mysqlimport" --host <host> --port <port> --user <user> --password=<password> --delete --local --verbose --ignore-lines={0} {1} "{2}"
""".format(rows_to_delete, schema if schema else database, windows_data_file)

    return [create_table_script, import_data_script_unix, import_data_script_windows]

def done_file(path):
    return os.path.join(os.path.dirname(path), '.'+os.path.basename(path)+'.done')

class cd():
    """Context manager for changing the current working directory

    Usage:
        with cd(new_dir):
            do something in new_dir, return to current directory when done
    """
    def __init__(self, newPath):
        self.newPath = newPath

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)

    def __exit__(self, etype, value, traceback):
        self.connection(close)

def mysql(cursor, query, print_query=False):
    """Execute each MySQL statement in query"""
    if isinstance(query, basestring):
        query = [query]
    for q in query:
        cursor.execute(q)
        q.replace('\n', ' ')
        if print_query:
            print('Executed: {}'.format(q))

def set_key_value(cursor, table, key_col, key, value_col, value):
    """Insert or update a key-value pair"""
    query = """INSERT INTO `{table}` (`{key_col}`, `{value_col}`) VALUES {key, value} ON DUPLICATE KEY UPDATE {`key_col`} = '{value}'""".format(table,key,value)
    mysql(cursor, query)

def get_key_value(cursor, table, key_col, key, value_col, value):
    """Get the value from a key-value pair"""
    query = """SELECT `{table}`.`{value_col}`) FROM `{table}` WHERE {`key_col`} = '{key}'""".format(table,key,value)
    mysql(cursor, query)
    return cursor.fetchone()

def get_key_value_dict(host, user, password, port, database, table, key_col, value_col,):
    """Get a key-value table from MySQL and return as a dict"""
    query = """SELECT `{key}`, `{value}` FROM `{database}`.`{table}`""".format(key=key_col, value=value_col, database=database, table=table)
    with MySQLdb.connect(host=host, user=user, passwd=password, port=port, db=database, cursorclass=MySQLdb.cursors.DictCursor) as cursor:
        try:
            cursor.execute(query)
            data = cursor.fetchall()
        except MySQLdb.Error, e:
            msg = ''
            try:
                msg = 'MySQL Error [{}]: {}'.format(e.args[0], e.args[1])
            except IndexError:
                msg = 'MySQL Error {}'.format(e)
            raise MyError('Database connection failed. {}'.format(msg))
        except:
            data = []
    if data == None:
        return {}
    else:
        return {d[key_col]: d[value_col] for d in data}

def set_key_value_dict(host, user, password, port, database, table, key_col, value_col, data):
    """Insert/update a MySQL key-value table from a dict"""
    with MySQLdb.connect(host=host, user=user, passwd=password, port=port, db=database,) as cursor:
        for key, value in data.items():
            query = """INSERT INTO `{database}`.`{table}` (`{key_col}`, `{value_col}`)
            VALUES ('{key}', '{value}')
            ON DUPLICATE KEY UPDATE `{value_col}` = '{value}'""".format(database=database, table=table, key_col=key_col, value_col=value_col, key=key, value=value)
            cursor.execute(query)

def list_hrefs(url, pattern=None, case_insensitive=True):
    """Scrape a web page for hrefs, filter by regex pattern"""
    if pattern:
        rx = re.compile(pattern, re.I) if case_insensitive else re.compile(pattern)
    with closing(requests.get(url)) as r:
        html = r.text
    soup = BeautifulSoup(html)
    hrefs = []
    for link in soup.find_all('a'):
        href =  link.get('href')
        if href and (pattern == None or rx.search(href)):
            hrefs.append(href)
    return hrefs

def list_ftp_files(host, directory, pattern=None, case_insensitive=True):
    """Get list of the filenames in an FTP directory"""
    with FTPHost(host, 'anonymous', '') as h:
        h.chdir(directory)
        names = h.listdir(h.curdir)
    if pattern:
        rx = re.compile(pattern, re.I) if case_insensitive else re.compile(pattern)
        names = [n for n in names if rx.match(n)]
    return names

def download_ftp_dir(host, remote_directory, local_directory, pattern=None, case_insensitive=True):
    """Download all files from an FTP directory; return list of downloaded files"""
    makedir(local_directory)
    done_file = os.path.join(local_directory, '.'+os.path.basename(remote_directory)+'.done')
    if file_exists(done_file):
        print('{} {} already downloaded; skipping. To reinstall "rm {}"'.format(host, os.path.basename(remote_directory), done_file))
    else:
        print('Downloading {} to {}.'.format(os.path.join(host, remote_directory), local_directory))
        with FTPHost(host, 'anonymous', '') as h:
            h.chdir(remote_directory)
            names = h.listdir(h.curdir)
            if pattern:
                rx = re.compile(pattern, re.I) if case_insensitive else re.compile(pattern)
                names = [n for n in names if rx.match(n)]
            if not names:
                print('WARNING: nothing to download from {}.'.format(os.path.join(host, remote_directory)))
            else:
                for name in names:
                    if h.path.isfile(name):
                        name_done_file = os.path.join(local_directory, '.'+name+'.done')
                        if file_exists(name_done_file):
                            print('{} {} already downloaded; skipping. To reinstall "rm {}"'.format(host, name, name_done_file))
                        else:
                            print('Downloading {} to {}.'.format(os.path.join(host, remote_directory, name), os.path.join(local_directory, name)))
                            h.download(name, os.path.join(local_directory, name))
                            with open(name_done_file, 'w'):
                                pass
                            print('Downloaded {} to {}.'.format(os.path.join(host, remote_directory, name), os.path.join(local_directory, name)))
                with open(done_file, 'w'):
                    pass
                print('Downloaded {} to {}.'.format(os.path.join(host, remote_directory), local_directory))
                return [os.path.join(local_directory, n) for n in names]

def getzip(url, zipfile, unzipdir): 
    """Download zipfile from url, extract contents to unzipdir, and remove zipfile."""
    done_file = os.path.join(unzipdir, '.'+os.path.basename(zipfile)+'.done')
    if file_exists(done_file):
        print('{} already downloaded and extracted; skipping. To reinstall "rm {}"'.format(os.path.basename(zipfile), done_file))
    else:
        print('Downloading  {} as {}.'.format(url, zipfile))
        urlretrieve(url, zipfile)
        print('Extracting {} into {}.'.format(zipfile, unzipdir))
        with ZipFile(zipfile, 'r') as zip:
            zip.extractall(unzipdir)
        os.remove(zipfile)
        with open(done_file, 'w'):
            pass

def getzip_requests(url, zipfile, unzipdir): 
    """(Deprecated; use getzip) Download zipfile from url, extract contents to unzipdir, and remove zipfile. Uses requests.get. Slow."""
    with closing(requests.get(url, stream=True)) as r:
        if r.headers.get('content-type') == None or r.headers.get('content-type') != 'application/zip':
            warning = "{} doesn't seem to be a zip file. Unzipping may fail.".format(url)
            warn(warning)
        with open(zipfile, 'wb') as fd:
            for chunk in r.iter_content():
                fd.write(chunk)
    with ZipFile(zipfile, 'r') as zip:
        zip.extractall(unzipdir)
    os.remove(zipfile)

def number_list(l):
    """Return a copy of a list with 1-based ordinal numbers before each item."""
    return ['{i:>{s}}. {v}'.format(s=len(str(len(l))), i=i+1, v=l[i]) for i in range(len(l))]

def select_item(items, default=None, title='Items', prompt='item'):
    """Print a numbered list and return user's choice"""
    selected_item = None
    print """
*** {} ***
{}
""".format(title, '\n'.join(number_list(items)))
    if default and default in items:
        default_index = items.index(default)+1
    else:
        default_index = None
    while selected_item == None:
        try:
            selected_index = raw_input('{}. Enter number (1 to {}) {}. 0 for none. :'
                                       .format(prompt, len(items), '' if default_index == None else '[{}]'.format(default_index)))
            selected_index = selected_index.strip().rstrip('.')
            if default_index != None and selected_index == '':
                selected_item = items[default_index-1]
            elif is_int(selected_index):
                selected_index = int(selected_index)
                if selected_index == 0:
                    return None
                elif (selected_index > 0 and selected_index <= len(items)):
                    selected_item = items[selected_index-1]
        except:
            pass
    return selected_item

def parse_vcf_meta_info(vcf_file):
    CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,GT = range(10)
    info_re = re.compile(r'INFO=<ID=(?P<id>[^,]+),Number=(?P<number>[^,]+),Type=(?P<type>[^,]+),Description="(?P<description>[^"]+)">', re.I)
    filter_re = re.compile(r'FILTER=<ID=(?P<id>[^,]+),Description="(?P<description>[^"]+)">', re.I)
    format_re = re.compile(r'FORMAT=ID=(?P<id>[^,]+),Number=(?P<number>[^,]+),Type=(?P<type>[^,]+),Description="(?P<description>[^"]+)">', re.I)
    alt_re = re.compile(r'ALT=ID=(?P<id>[^,]+),Description=(?P<description>[^>]+)>', re.I)
    contig_re = re.compile(r'contig=ID=(?P<id>[^,]+),URL=(?P<url>[^>]+)>', re.I)
    sample_re = re.compile(r'SAMPLE=ID=(?P<id>[^,]+),Genomes=(?P<genomes>[^,]+),Mixture=(?P<mixture>[^,]+),Description=(?P<description>[^>]+)>', re.I)
    result = {'column_names': [], 'GTS': [], 'INFO': {}, 'FILTER': {}, 'FORMAT': {}}
    with open_gz_or_text(vcf_file) as v:
        for line in v:
            line = line.rstrip('\n')
            if line == '':
                continue
            if not line.startswith('#'):
                break
            if line.upper().startswith('#CHROM'):
                result['column_names'] = line[1:].split('\t')
                if len(result['column_names']) > GT:
                    result['GTS'] = result['column_names'][GT:]
            if line.startswith('##'):
                line = line.lstrip('#')
                if '=' not in line:
                    result[line] = None
                else:
                    if '=<' in line:
                        k,v = line.split('=<')
                    else:
                        k,v = line.split('=')
                    if k.upper() == 'INFO':
                        m = info_re.match(line)
                        if m:
                            result['INFO'][m.group('id')] = {'Number': m.group('number'), 'Type': m.group('type'), 'Description': m.group('description')}
                    elif k.upper() == 'FILTER':
                        m = filter_re.match(line)
                        if m:
                            result['FORMAT'][m.group('id')] = {'Description': m.group('description')}
                    elif k.upper() == 'FORMAT':
                        m = format_re.match(line)
                        if m:
                            result['FILTER'][m.group('id')] = {'Number': m.group('number'), 'Type': m.group('type'), 'Description': m.group('description')}
                    elif k.upper() == 'ALT':
                        m = alt_re.match(line)
                        if m:
                            result['ALT'][m.group('id')] = {'Number': m.group('number'), 'Description': m.group('description')}
                    elif k.upper() == 'contig':
                        m = contig_re.match(line)
                        if m:
                            result['contig'][m.group('id')] = {'URL': m.group('url').split(',')}
                    elif k.upper() == 'SAMPLE':
                        m = sample_re.match(line)
                        if m:
                            result['SAMPLE'][m.group('id')] = {'Genomes': m.group('genomes').split(';'), 'Mixture': m.group('mixture').split(';'), 'Description': m.group('description').split(';')}
                    elif k.upper() == 'PEDIGREE':
                        names_ids = line.lstrip('<').rstrip('>').split(',')
                        names_ids = [ni.split('=') for ni in names_ids]
                        result['PEDIGREE'] = {ni[0]: ni[1] for ni in names_ids}
                    else:
                        result[k] = v
    return result

