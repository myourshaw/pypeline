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
#-d /scratch1/tmp/myourshaw/gmd -q /data/storage-1-02/solexa/110715_SN973_0041_AC041EACXX/Data/Intensities/BaseCalls/QSEQ/s_[78]_[123]_*_qseq.txt.gz

#110812_SN860_0079_2011-143_AC03YLACXX
#-d /scratch1/tmp/myourshaw/gmd -q /data/storage-1-00/HiSeq2k.Freimerlab/110812_SN860_0079_2011-143_AC03YLACXX/Data/Intensities/BaseCalls/QSEQ/s_[1-8]_[123]_*_qseq.txt.gz


class Qseqs2RgsError(Exception): pass

def main():

    #command line arguments
    parser = argparse.ArgumentParser(parents=[my.default_parser()],
        description = 'list of qseq files > topdir/qseqs/machine.run.lane.read[.barcode]_qseq.txt.gz',
        epilog = 'pypeline.qseq2rg version 1.0β1 ©2011-2012 Michael Yourshaw all rights reserved')
    parser.add_argument('--qseqs', '-q', nargs='+', required=True,
                                            help='list of qseq files')
    parser.add_argument('--barcode_file', '-b',
                                            help='novobarcode-style tag file (default: TruSeq/SureSelectXT set of 12 barcodes)')
    parser.add_argument('--novobarcode','-n',
                                            help='path to novobarcode executable (default: from config->novobarcode)')
    args = parser.parse_args()
    dirs = my.create_default_directory_structure(args.top_dir)
    config = my.get_config(args)
    
    my_name = '{}_{}'.format('qseqs2rgs', my.localtime_squish())
    job_dir = my.unique_dir(my_name, dirs['jobs'])
    
    #get list of qseq files for each tile and examine first read in each file to get machine, run, lane, read
    qseqs = []
    for q in args.qseqs:
        qseqs += glob(q)
    qseqs = {q for q in qseqs}
    peeks = [my.qseq_peek(q) for q in qseqs]
    
    #merge qseq tile files to get one file per machine, run, lane, read
    cmds= []
    machines = sorted({x['machine'] for x in peeks})
    for machine in machines:
        runs = sorted({x['run'] for x in peeks if x['machine'] == machine})
        for run in runs:
            lanes = sorted({x['lane'] for x in peeks if x['machine'] == machine and x['run'] == run})
            for lane in lanes:
                reads = sorted({x['read'] for x in peeks if x['machine'] == machine and x['run'] == run and x['lane'] == lane})
                #check for matching tile files for each read
                read_dict = {read: {x['tile'] for x in peeks if x['machine'] == machine and x['run'] == run and x['lane'] == lane and x['read'] == read} for read in reads}
                u = set()
                for r in read_dict.keys():
                    u |= read_dict[r]
                for r in read_dict.keys():
                    if read_dict[r] != u:
                        raise Qseqs2RgsError("Read sets do not have identical tiles")
                merged_files_dir = os.path.join(dirs['qseqs'],str(machine),str(run),str(lane))
                my.makedir(merged_files_dir)
                for read in reads:
                    merged_file = os.path.join(merged_files_dir,'s_{}_{}_qseq.txt.gz'.format(lane,read))
                    files_to_merge = []
                    compression = set()
                    tiles = sorted({x['tile'] for x in peeks if x['machine'] == machine and x['run'] == run and x['lane'] == lane and x['read'] == read})
                    for tile in tiles:
                        files_to_merge += [x['path'] for x in peeks if x['machine'] == machine and x['run'] == run and x['lane'] == lane and x['read'] == read and x['tile'] == tile]
                        compression |= {x['compressed'] for x in peeks if x['machine'] == machine and x['run'] == run and x['lane'] == lane and x['read'] == read and x['tile'] == tile}
                    if True in compression and False in compression:
                        raise Qseqs2RgsError("Compressed and uncompressed qseqs cannot be mixed")
                    if True in compression:
                        cmds += ['gzip -cd {} | cat - | gzip > {}'.format(' '.join(files_to_merge), merged_file)]
                    else:
                        cmds += ['cat {} | gzip > {}'.format(' '.join(files_to_merge), merged_file)]
    job_name = my_name+'_merge'
    job = my.run_job(cmds, job_name, job_dir, synchronous=True)
    print job
    
    novobarcode = args.novobarcode if args.novobarcode and my.is_exe(args.novobarcode) else config.get('novocraft','novobarcode')
    
    qseq_dir = os.getcwd() if args.qseq_dir == None else os.path.abspath(args.qseq_dir)
    exptrun_match = re.search(r"(?P<date>\d+)_(?P<instrument>[^_]+)_(?P<run>\d+)_([^_]+_)?(?P<flowcell>[^_/]+)",qseq_dir)
    exptrun = args.exptrun if args.exptrun else exptrun_match.group(0) if exptrun_match else re.sub(r"[^a-zA-Z0-9_.]",'_', qseq_dir).strip('_')
    readgroup_tmpdir = my.unique_dir('readgroup_{}_tmp'.format(exptrun), args.tmp_dir)
    qout_dir = os.path.join(readgroup_tmpdir, 'qout')
    my.makedir(qout_dir)
    
    lane_barcodes = {}
    NUM_LANES = 8
    if not args.lane_barcodes:
        for l in range(1, NUM_LANES+1):
            lane_barcodes[str(l)] = set()
    else:
        lane = None
        for lb in args.lane_barcodes:
            if lb.isdigit():
                lane = lb
                lane_barcodes.setdefault(lane,set())
            elif lane:
                if lb == '*' or my.file_exists(lb):
                    lane_barcodes[lane] = lb
                elif re.search(r"[^ACGT]",lb,re.I):
                    raise Qseq2RgError('{} is not a valid barcode or barcode file'.format(lb))
                else:
                    lane_barcodes[lane].add(lb)
            else:
                raise Qseq2RgError('{} is not a valid lane number'.format(lb))
    cmds = []
    novobarcode_output_dir_list = []
    for lane, barcodes  in lane_barcodes.items():
        qseq_reads = []
        for read in ('1','2','3'):
            compressed = False
            qseq_reads.append(my.get_qseqs(qseq_dir, lane, read, gz=False))
            if not qseq_reads:
                compressed = True
                qseq_reads.append(my.get_qseqs(qseq_dir, lane, read, gz=True))
        if not qseq_reads:
            raise Qseq2RgError('no files for lane {} in {}'.format(lane, qseq_dir))
        elif not barcodes:
            for qseqs in qseq_reads:
                if compressed:
                    cmds.append('gzip -cd {} | cat - | gzip > {}'.format(' '.join(qseqs), os.path.join(args.readgroup_qseqs_dir,'{}.{}.{}_qseq.txt.gz'.format(exptrun, lane, read))))
                else:
                    cmds.append('cat {} | gzip > {}'.format(' '.join(qseqs), os.path.join(args.readgroup_qseqs_dir,'{}.{}.{}_qseq.txt.gz'.format(exptrun, lane, read))))
        else:
            #read 2 [1] is the barcode
            if not qseq_reads[1]:
                raise Qseq2RgError('no barcode reads for lane {} in {}'.format(lane, qseq_dir))
            elif len(qseq_reads[0]) != len(qseq_reads[1]) or (qseq_reads[2] and len(qseq_reads[2]) != len(qseq_reads[1])):
                raise Qseq2RgError('number of qseq files for lane {} in {} does not match number of barcode qseq files'.format(lane, qseq_dir))
            else:
                #possible lane-specific barcode file
                default_barcode_file = lane_barcodes[lane] if len(lane_barcodes[lane]) == 1 and my.file_exists(lane_barcodes[lane]) else args.default_barcode_file
                cmds.append('python {} --file1s {} --file2s {} --qseqtagfiles {} --exptrun {} --lane {} --barcodes {} --default_barcode_file {} --novobarcode {} --readgroup_qseqs_dir {} {} {} {} {} {} {}'\
                 .format(qseq2rg_novobarcode_and_cat.__file__, ' '.join(qseq_reads[0]), ' '.join(qseq_reads[2]), ' '.join(qseq_reads[1]), exptrun, lane, ' '.join(barcodes), args.default_barcode_file, novobarcode, args.readgroup_qseqs_dir,  '--email '+args.email if args.email else '',	'--sendemailforintermediatejobs' if args.sendemailforintermediatejobs else '',	'--deletetempfiles' if args.deletetempfiles else '',	'--restart' if args.restart else '',	'--config '+args.config if args.config else '', '--tmpdir '+args.tmpdir if args.tmpdir else ''))
    
    jobName = 'qseq2rg_{}_{}'.format(exptrun, my.localtime_squish())
    my.print_log('job {} started; commands:\n{}'.format(jobName, '\n'.join(cmds)))
    j = job.Job(jobName = jobName,
    outputPath = args.qout_dir, errorPath = args.qout_dir, workingDirectory = readgroup_tmpdir,
    email = args.email, blockEmail = not args.sendemailforintermediatejobs)
    
    try:
        j.executeCommands(cmds, synchronous=True)
    except Exception as e:
        my.print_err('There was an error running the job {}:\n{}\nintermediate files are in {}'.format(jobName, e, readgroup_tmpdir))
    else:
        my.print_log('job {} finished; status:\n{}'.format(jobName, '\n'.join(j.completionMsg)))
        with open(os.path.join(readgroup_tmpdir,'{}_job_status_report.txt'.format(jobName)), 'w') as s:
            s.write(format('\n'.join(j.completionMsg)))
    
            #qseq statistics
            my.print_log('calculating qseq metrics started')
            qseq_metrics.gather_qseq_stats(my.qseqs2qseqtrios(my.get_qseqs(args.readgroup_qseqs_dir,'*','*')), os.path.join(args.qseq_metrics_dir, '{}.qseq_metrics'.format(jobName)))
            my.print_log('calculating qseq metrics finished; results are in [{}]'.format(args.qseq_metrics_dir))
    
        if args.deletetempfiles:
            shutil.rmtree(readgroup_tmpdir, ignore_errors=True)

if __name__ == "__main__": sys.exit(main())
