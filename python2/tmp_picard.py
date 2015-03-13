#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import argparse
from ConfigParser import SafeConfigParser #configparser in python 3
import shutil
import my
import job

class PicardError(Exception): pass

#-t BuildBamIndex INPUT=/scratch0/tmp/myourshaw/gmd/bams/test/GMD1A_XT_Illumina.2011-05-11.HWI-ST430.243.6.ATCACGA.novoalign.bam
#-t BamIndexStats --job_dir /scratch0/tmp/myourshaw/gmd1/jobs/bam_metrics_BHIqqb --picard_args INPUT=/scratch0/tmp/myourshaw/gmd1/bams/readgroup_bams/HWI-ST430.243.6.ATCACGA.novoalign.fixmate.bam --redirect_stdout /scratch0/tmp/myourshaw/gmd1/metrics/bam_metrics/testHWI-ST430.243.6.ATCACGA.novoalign.fixmate.bam.metrics.indexstats

def __picard_info(reference=None):
    if not reference:
        config = SafeConfigParser()
        config.readfp(open(os.path.join(os.path.dirname(__file__),'pypeline.cfg')))
        cfg_files = [os.path.expanduser('~/.pypeline.cfg')]
        config.read(cfg_files)
        reference = config.get('reference','default_reference')
    return {
    'ValidateSamFile': {'defaults': {'INPUT': None, 'OUTPUT': None, 'MODE': 'VERBOSE', 'REFERENCE_SEQUENCE': reference, 'VALIDATE_INDEX': 'true',},
        'run_parameters': {'processors': 8, 'memory': '6G'}},
    'FixMateInformation': {'defaults': {'INPUT': None, 'CREATE_INDEX': 'true', 'CREATE_MD5_FILE': 'true', 'SORT_ORDER': 'coordinate',},
        'run_parameters': {'processors': 8, 'memory': '6G'}},
    'MergeSamFiles': {'defaults': {'OUTPUT': None, 'CREATE_INDEX': 'true', 'CREATE_MD5_FILE': 'true', 'SORT_ORDER': 'coordinate', 'ASSUME_SORTED': 'false', 'MERGE_SEQUENCE_DICTIONARIES': 'false', 'USE_THREADING': 'true',},
        'run_parameters': {'processors': 8, 'memory': '6G'}},
    'MarkDuplicates': {'defaults': {'INPUT': None, 'OUTPUT': None, 'CREATE_INDEX': 'true', 'CREATE_MD5_FILE': 'true', 'REMOVE_DUPLICATES': 'false', 'MAX_FILE_HANDLES_FOR_READ_ENDS_MAP': 8000, 'READ_NAME_REGEX': '"[^_]+_[^_]+_[0-9]+_([0-9]+)_([0-9.]+)_([0-9.]+).*"', 'OPTICAL_DUPLICATE_PIXEL_DISTANCE': 100, 'METRICS_FILE': None,},
        'run_parameters': {'processors': 8, 'memory': '6G'}},
    'BuildBamIndex': {'defaults': {'INPUT': None},
        'run_parameters': {'processors': 8, 'memory': '6G'}},
    'BamIndexStats': {'defaults': {'INPUT': None},
        'run_parameters': {'processors': 8, 'memory': '6G'}},
    'CalculateHsMetrics': {'defaults': {'INPUT': None, 'OUTPUT':None, 'BAIT_INTERVALS':None, 'TARGET_INTERVALS':None,},
        'run_parameters': {'processors': 8, 'memory': '6G'}},
    'CollectMultipleMetrics': {'defaults': {'INPUT': None, 'OUTPUT':None, 'REFERENCE_SEQUENCE': reference,},
        'run_parameters': {'processors': 8, 'memory': '6G'}},
    'CollectGcBiasMetrics': {'defaults': {'INPUT': None, 'OUTPUT':None, 'REFERENCE_SEQUENCE': reference, 'CHART_OUTPUT':None, 'SUMMARY_OUTPUT':None,},
        'run_parameters': {'processors': 8, 'memory': '6G'}},
    'EstimateLibraryComplexity': {'defaults': {'INPUT': None, 'OUTPUT':None, 'READ_NAME_REGEX':r'"[^_]+_[^_]+_[0-9]+_([0-9]+)_([0-9.]+)_([0-9.]+).*"', 'OPTICAL_DUPLICATE_PIXEL_DISTANCE':100,},
        'run_parameters': {'processors': 8, 'memory': '6G'}},
}

def run(tool, picard_args, multi_valued_arg_strings=None, redirect_stdout=None, job_dir=None,
        reference=None, picard_current=None, email=None, synchronous=False, hold_jid=None):
    config = my.get_config()
    picard_info = __picard_info(reference)
    tool_info = picard_info.get(tool)
    if not tool_info:
        raise PicardError('unable to run {}'.format(tool))
    out_dir = os.path.dirname(picard_args['OUTPUT']) if picard_args.get('OUTPUT') else  None
    if out_dir:
        my.makedir(out_dir)
    if not job_dir:
        job_dir = my.unique_dir('picard_'+tool,out_dir) if out_dir else my.unique_dir('picard_'+tool)
    my.makedir(job_dir)
    tmp_dir = picard_args['TMP_DIR'] if picard_args.get('TMP_DIR') else my.unique_dir('picard_tmp_'+tool, job_dir)
    my.makedir(tmp_dir)
    if not picard_current:
        picard_current = config.get('picard','picard_current')
    jar = os.path.join(picard_current, tool+'.jar')
    if not my.file_exists(jar):
        raise PicardError("Can't find {}".format(jar))
    picard_options = tool_info['defaults']
    if picard_args:
        for k in picard_args.keys():
            if isinstance(picard_args[k], (list, tuple)):
                picard_options[k] = ' {}='.format(k).join(picard_args[k])
            else:
                picard_options[k] = picard_args[k]
    for o in picard_options.keys():
        if picard_options[o] == None:
            raise PicardError('missing value for picard {} option {}'.format(tool, o))
    option_str = picard_options = ' '.join(['{}={}'.format(k,picard_options[k]) for k in picard_options.keys()])
    if isinstance(multi_valued_arg_strings, (list, tuple)):
        multi_valued_arg_strings = ' '.join(multi_valued_arg_strings)
    if multi_valued_arg_strings:
        option_str += ' '+multi_valued_arg_strings
    cmd = 'java -Xmx5g -Djava.io.tmpdir={tmp_dir} -jar {jar} TMP_DIR={tmp_dir} VERBOSITY=INFO QUIET=false VALIDATION_STRINGENCY=STRICT COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=1000000 {option_str};'.format(tmp_dir=tmp_dir, jar=jar, option_str=option_str)
    if redirect_stdout:
        cmd = cmd.rstrip(';')+' > '+redirect_stdout+';'
    jobName = 'picard_{}{}'.format(tool, '_'+os.path.basename(picard_args['OUTPUT']) if picard_args.get('OUTPUT') else '')
    qout_dir = job_dir
    processors = tool_info['run_parameters']['processors']
    memory = tool_info['run_parameters']['memory']
    job = my.run_job(cmd, jobName, job_dir, qout_dir, email, synchronous=synchronous, processors=processors, memory=memory, hold_jid=hold_jid)
    if synchronous:
        shutil.rmtree(tmp_dir, ignore_errors=True)
    return job

def main():
    
    #command line arguments
    parser = argparse.ArgumentParser(#parents=[my.default_parser()],
        description = 'Wrapper to run picard commands',
        epilog = 'pypeline.picard version 1.0β1 ©2011-2012 Michael Yourshaw all rights reserved')
    parser.add_argument('--tool', '-t', required=True, choices=__picard_info().keys(),
        help='name of picard tool',)
    parser.add_argument('--picard_current',
        help='path/to/picard/directory',)
    parser.add_argument('--job_dir',
        help='path/to/job/directory',)
    parser.add_argument('--email', nargs='*', default=None,
        help='email address(es) to receive job completion notification',)
    parser.add_argument('--synchronous', action='store_true', default=False,
        help='run job synchronously')
    parser.add_argument('--hold_jid', nargs='*', default=None,
        help='space-separated list of job ids that must complete before this job runs')
    parser.add_argument('--picard_args', nargs='*',
        help='keyword=value arguments specific to the picard tool, e.g., INPUT=foo/bar.bam',)
    parser.add_argument('--multi_valued_arg_strings', nargs='*',
        help='keyword=value arguments where a keyword is used multiple times',)
    parser.add_argument('--redirect_stdout', default=None,
        help="file to receive stdout for tools such as BamIndexStats that don't specify OUTPUT",)
    parser.add_argument('--reference', '-r',
        help='path to reference sequence fasta file (default: from config->reference.default_reference)')
    args = parser.parse_args()
    
    picard_args = {k.split('=')[0]: k.split('=')[1] for k in args.picard_args}
    job = run(tool=args.tool, picard_args=picard_args, multi_valued_arg_strings=args.multi_valued_arg_strings,
        redirect_stdout=args.redirect_stdout, job_dir=args.job_dir, reference=args.reference,
        picard_current=args.picard_current, email=args.email, synchronous= args.synchronous,
        hold_jid=args.hold_jid)
    

if __name__ == "__main__": sys.exit(main())
