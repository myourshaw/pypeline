#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import argparse
from ConfigParser import SafeConfigParser
import re
from glob import glob, iglob
import my
import job
import qseq2rg_novobarcode_and_cat
import qseq_metrics

#
#110715_SN973_0041_AC041EACXX
#-d /scratch1/tmp/myourshaw/gmd -q /scratch1/tmp/myourshaw/gmd/qseqs/HWI-ST973/41/7/s_7_[123]_qseq.txt.gz
#-d /scratch1/tmp/myourshaw/gmd -q /scratch1/tmp/myourshaw/gmd/qseqs/HWI-ST973/41/8/s_8_[123]_qseq.txt.gz

#110812_SN860_0079_2011-143_AC03YLACXX
#-d /scratch1/tmp/myourshaw/gmd -q /scratch1/tmp/myourshaw/gmd/qseqs/HWI-ST0860/79/1/s_1_[123]_qseq.txt.gz
#-d /scratch1/tmp/myourshaw/gmd -q /scratch1/tmp/myourshaw/gmd/qseqs/HWI-ST0860/79/2/s_2_[123]_qseq.txt.gz
#-d /scratch1/tmp/myourshaw/gmd -q /scratch1/tmp/myourshaw/gmd/qseqs/HWI-ST0860/79/3/s_3_[123]_qseq.txt.gz
#-d /scratch1/tmp/myourshaw/gmd -q /scratch1/tmp/myourshaw/gmd/qseqs/HWI-ST0860/79/4/s_4_[123]_qseq.txt.gz
#-d /scratch1/tmp/myourshaw/gmd -q /scratch1/tmp/myourshaw/gmd/qseqs/HWI-ST0860/79/5/s_5_[123]_qseq.txt.gz
#-d /scratch1/tmp/myourshaw/gmd -q /scratch1/tmp/myourshaw/gmd/qseqs/HWI-ST0860/79/6/s_6_[123]_qseq.txt.gz
#-d /scratch1/tmp/myourshaw/gmd -q /scratch1/tmp/myourshaw/gmd/qseqs/HWI-ST0860/79/7/s_7_[123]_qseq.txt.gz
#-d /scratch1/tmp/myourshaw/gmd -q /scratch1/tmp/myourshaw/gmd/qseqs/HWI-ST0860/79/8/s_8_[123]_qseq.txt.gz

#-d /scratch1/tmp/myourshaw/gmd -q /scratch1/tmp/myourshaw/gmd/qseqs/*/*/*/s_*_[123]_qseq.txt.gz

class DemultiplexLaneQseqsError(Exception): pass

def run(config, dirs, metadata, qseqs, barcodefile, novobarcode):

    novobarcode = novobarcode if novobarcode and my.is_exe(novobarcode) else config.get('novocraft','novobarcode')
    if not my.is_exe(novobarcode):
        raise DemultiplexLaneQseqsError ("Can't find novoalign executable")
    
    barcodefile = barcodefile if barcodefile and my.file_exists(barcodefile) else config.get('novocraft','truseq_xt_barcodes')
    if not my.file_exists(barcodefile):
        raise DemultiplexLaneQseqsError ("Can't find barcode tag file")
    barcodes = ['NC']
    with open(barcodefile, 'r') as b:
        for line in b:
            line = line.strip()
            if not line:
                continue
            fields = line.split()
            if fields[0].lower() in ('distance','format'):
                continue
            if len(fields) > 1 and not re.search(r"[^ACGTN.]",fields[1]):
                barcodes += [fields[1]]

    my_name = '{}_{}'.format('demultiplex_lane_qseqs', my.localtime_squish())
    job_dir = my.unique_dir(my_name, dirs['jobs'])

    #get list of qseq files for each lane and examine first read in each file to get machine, run, lane, read
    qseqs = set(my.flatten([glob(q) for q in qseqs]))
    peeks = [my.qseq_peek(q) for q in qseqs]

    #for each machine, run, lane get matching qseq read files; read 2 is the read of the barcode
    print 'Demultiplexing lane qseqs'
    print '\t'.join(['machine', 'run', 'lane', 'demultiplexed_qseq_files_folder'])
    cmds = []
    dirs = []
    machines = sorted({x['machine'] for x in peeks})
    for machine in machines:
        runs = sorted({x['run'] for x in peeks if x['machine'] == machine})
        for run in runs:
            lanes = sorted({x['lane'] for x in peeks if x['machine'] == machine and x['run'] == run})
            for lane in lanes:
                reads = sorted({x['read'] for x in peeks if x['machine'] == machine and x['run'] == run and x['lane'] == lane})
                if 1 in reads:
                    file1 = [x['path'] for x in peeks if x['machine'] == machine and x['run'] == run and x['lane'] == lane and x['read'] == 1][0]
                elif 3 in reads:
                    file1 = [x['path'] for x in peeks if x['machine'] == machine and x['run'] == run and x['lane'] == lane and x['read'] == 3][0]
                else:
                    continue
                if 2 in reads:
                    qseqtagfile = [x['path'] for x in peeks if x['machine'] == machine and x['run'] == run and x['lane'] == lane and x['read'] == 2][0]
                else:
                    continue
                if 3 in reads and 1 in reads:
                    file2 = [x['path'] for x in peeks if x['machine'] == machine and x['run'] == run and x['lane'] == lane and x['read'] == 3][0]
                else:
                    file2 = ''
                output_folder = os.path.dirname(file1)
                dirs += [output_folder]
                print '\t'.join([str(machine), str(run), str(lane), str(output_folder)])
                cmds += ['{0} -b {1} -d {2} -F QSEQ -f {3} {4} -i {5} --ILQ_SKIP --QSEQ_OUT --NC_OFF;'\
                    .format(novobarcode, barcodefile, output_folder, file1, file2, qseqtagfile)]

    job_name = my_name
    print 'Running job {}. Job_info, STDOUT, and STDERR are in {}. ...'.format(job_name, job_dir)    
    job = my.run_job(cmds, job_name, job_dir, synchronous=True)
    print 'Job {} done. Job_info, STDOUT, and STDERR are in {}.'.format(job_name, job_dir)
    print 'Compressing and renaming demultiplexed qseq files. ...'
    print '\t'.join(['machine', 'run', 'lane', 'read', 'barcode', 'demultiplexed_qseq_file'])
    cmds = []
    for dir in dirs:
        for barcode in barcodes:
            barcode_dir = os.path.join(dir, barcode)
            if os.path.isdir(barcode_dir):
                for file in os.listdir(barcode_dir):
                    qseq = os.path.join(barcode_dir, file)
                    if os.path.isfile(qseq):
                        peek = my.qseq_peek(qseq)
                        if peek:
                            compressed_file = os.path.join(dir, '{}.{}.{}.{}.{}_qseq.txt.gz'
                                .format(peek['machine'], peek['run'], peek['lane'], peek['read'], barcode))
                            print '\t'.join([str(peek['machine']), str(peek['run']), str(peek['lane']), str(peek['read']), str(barcode), str(compressed_file)])
                            cmds += ['gzip {}; mv {} {};'.format(qseq, qseq+'.gz', compressed_file)]
    job_name = my_name+'_compress'
    print 'Running job {}. Job_info, STDOUT, and STDERR are in {}. ...'.format(job_name, job_dir)    
    job = my.run_job(cmds, job_name, job_dir, synchronous=True)
    print 'Job {} done. Job_info, STDOUT, and STDERR are in {}.'.format(job_name, job_dir)

def main():

    #command line arguments
    parser = argparse.ArgumentParser(parents=[my.default_parser()],
        description = 'list of qseq lane files (all multiplexed with the same barcode set) > qseq_dir/machine.run.lane.read.barcode_qseq.txt.gz',
        epilog = 'pypeline.demultiplex_lane_qseqs version 1.0β1 ©2011-2012 Michael Yourshaw all rights reserved')
    parser.add_argument('--qseqs', '-q', nargs='+', required=True,
        help='list of qseq files')
    parser.add_argument('--barcodefile', '-b',
        help='novobarcode-style tag file (default: TruSeq/SureSelectXT set of 12 barcodes)')
    parser.add_argument('--novobarcode','-n',
        help='path to novobarcode executable (default: from config->novobarcode)')
    args = parser.parse_args()
    dirs = my.create_default_directory_structure(args.top_dir)
    config = my.get_config(args)
    
    job = run(config, dirs, args.qseqs, args.barcodefile, args.novobarcode)

if __name__ == "__main__": sys.exit(main())
