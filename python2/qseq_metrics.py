#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import shutil
from glob import glob
import re
import my
import job
import qseq_metrics_calculator
import qseq_metrics_merger

#--tiles -d /scratch0/tmp/myourshaw/gmd1 -j /scratch0/tmp/myourshaw/gmd1/jobs/qseq_metrics_tiles_test -o /scratch0/tmp/myourshaw/gmd1/jobs/qseq_metrics_tile_test/qseq_metrics_tiles_test.txt -q /data/storage-1-04/archive/illumina/110616_SN860_0065_2011-101_B817EAABXX/Data/Intensities/BaseCalls/QSEQ/s_[1-7]_[123]_*_qseq.txt /data/storage-1-02/archive/myourshaw/illumina/110623_SN860_0067_2011-100R_A81MVKABXX/Data/Intensities/BaseCalls/QSEQ/s_[1-8]_[123]_*_qseq.txt /data/storage-1-00/HiSeq2k.Freimerlab/110812_SN860_0079_2011-143_AC03YLACXX/Data/Intensities/BaseCalls/QSEQ/s_[1-8]_[123]_*_qseq.txt.gz /data/storage-1-02/solexa/110511_SN430_0243_B817FLABXX/Data/Intensities/BaseCalls/s_[1-7]_[123]_*_qseq.txt /data/storage-1-02/solexa/110715_SN973_0041_AC041EACXX/Data/Intensities/BaseCalls/QSEQ/s_[78]_[123]_*_qseq.txt.gz

#MMS
#python /home/myourshaw/lab/pypeline/python/qseq_metrics.pyc -d /scratch0/tmp/myourshaw/mms1 -q /scratch0/tmp/myourshaw/mms1/qseqs/HWUSI-EAS335/61DRC/4/HWUSI-EAS335.61DRC.4.1_qseq.txt.gz /scratch0/tmp/myourshaw/mms1/qseqs/HWUSI-EAS335/61DRC/3/HWUSI-EAS335.61DRC.3.1_qseq.txt.gz /scratch0/tmp/myourshaw/mms1/qseqs/HWUSI-EAS172/61DRF/8/HWUSI-EAS172.61DRF.8.1_qseq.txt.gz /scratch0/tmp/myourshaw/mms1/qseqs/HWUSI-EAS172/61DRE/6/HWUSI-EAS172.61DRE.6.1_qseq.txt.gz -o /scratch0/tmp/myourshaw/mms1/metrics/qseq_metrics/readgroup_qseq_metrics_20110922093617.txt -b '\.([ACGT]+)_qseq.txt'

class QseqMetricsError(Exception): pass


def main():
    import argparse
    parser = argparse.ArgumentParser(parents=[my.default_parser()],
        description = 'Calculate pass/fail stats and validate matching read ids of qseq files',
        epilog = 'pypeline.qseq_metrics version 1.0β1 ©2011-2012 Michael Yourshaw all rights reserved')
    parser.add_argument('--qseqs', '-q', required=True, nargs='+',
        help='list of qseq files')
    parser.add_argument('--output', '-o',
        help='output file')
    parser.add_argument('--barcode_re', '-b', default=r"\.([ACGT]+)_qseq.txt",
        help='regular expression for barcode in qseq path/filename (default("\.([ACGT]+)_qseq.txt")')
    parser.add_argument('--tiles', action='store_true', default=False,
        help='qseqs are raw tile qseqs (default: qseqs are readgroup qseqs)')
    parser.add_argument('--job_dir', '-j',
        help='parent job directory')
    args = parser.parse_args()

    dirs = my.create_default_directory_structure(args.top_dir)
    config = my.get_config(args)
    
    my_name = '{}_{}'.format('qseq_metrics', my.localtime_squish())
    my_job_dir = my.unique_dir(my_name, args.job_dir) if args.job_dir else my.unique_dir(my_name, dirs['jobs'])

    args.output = args.output if args.output else os.path.join(dirs['qseq_metrics'], 'qseq_metrics.{}.txt'.format(my.localtime_squish()))
    my.makedir(os.path.dirname(args.output))

    qseq_info = my.qseqs_hierarchy(args.qseqs, args.barcode_re) if not args.tiles else my.qseq_tiles_reads_hierarchy(args.qseqs)
    if not qseq_info or len(qseq_info) == 0:
        raise QseqMetricsError('No qseq files')
    
    job_dir = my_job_dir
    cmds= []
    metrics_tmps = []
    #calculate metrics for the read files of a group
    for machine in sorted(qseq_info):
        for run in sorted(qseq_info[machine]):
            for lane in sorted(qseq_info[machine][run]):
                for barcode_tile in sorted(qseq_info[machine][run][lane]):
                    qseq_readgroup = []
                    for read in sorted(qseq_info[machine][run][lane][barcode_tile]):
                        qseq_readgroup += [qseq_info[machine][run][lane][barcode_tile][read]['path']]
                    metrics_tmp_fh, metrics_tmp = my.unique_file(os.path.join(my_job_dir, 'metrics.{}.{}.{}.{}'.format(machine,run,lane,barcode_tile)))
                    metrics_tmps += [metrics_tmp]
                    cmds += ['python {} -q {} -o {}'.format(
                        qseq_metrics_calculator.__file__, ' '.join(my.flatten(qseq_readgroup)), metrics_tmp)]
    job_name = my_name+'_calculator'
    job = my.run_job(cmds, job_name, job_dir, synchronous=False, processors=8)
    qseq_metrics_calculator_jobid = job.jobId

    #merge all the results into a single file
    cmd = 'python {} -i {} -o {}'.format(
        qseq_metrics_merger.__file__, ' '.join(metrics_tmps), args.output)
    job_name = my_name+'_merger'
    job = my.run_job(cmd, job_name, job_dir, synchronous=False, processors=8, hold_jid=qseq_metrics_calculator_jobid)
    qseq_metrics_merger_jobid = job.jobId


if __name__ == "__main__": sys.exit(main())
