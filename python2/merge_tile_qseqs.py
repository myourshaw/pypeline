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

#python /home/myourshaw/lab/pypeline/python/merge_tile_qseqs.py -d /scratch0/tmp/myourshaw/gmd \
#-q /data/storage-1-04/archive/illumina/110616_SN860_0065_2011-101_B817EAABXX/Data/Intensities/BaseCalls/QSEQ/s_[1-7]_[123]_*_qseq.txt \
#/data/storage-1-02/archive/myourshaw/illumina/110623_SN860_0067_2011-100R_A81MVKABXX/Data/Intensities/BaseCalls/QSEQ/s_[1-7]_[123]_*_qseq.txt;

#MMS
#-d /scratch0/tmp/myourshaw/mms -q /home/solexa/solexa_datasets/100528_HWUSI-EAS335_61DRC/Data/Intensities/BaseCalls/s_3_1_*_qseq.txt* /home/solexa/solexa_datasets/100528_HWUSI-EAS335_61DRC/Data/Intensities/BaseCalls/s_4_1_*_qseq.txt* /home/solexa/solexa_datasets/100608_HWUSI-EAS172_61DRF/Data/Intensities/BaseCalls/s_8_1_*_qseq.txt* /home/solexa/solexa_datasets/100616_HWUSI-EAS172_61DRE/Data/Intensities/BaseCalls/s_6_1_*_qseq.txt*

class MergeTileQseqsError(Exception): pass

def run(config, dirs, metadata, qseqs, hold_jid):
    my_name = '{}_{}'.format('merge_tile_qseqs', my.localtime_squish())
    job_dir = my.unique_dir(my_name, dirs['jobs'])
    
    #get list of qseq files for each tile and examine first read in each file to get machine, run, lane, read
    qseqs = my.qseq_tiles_hierarchy(qseqs)

    #merge qseq tile files to get one file per machine, run, lane, read
    print 'Merging tile qseqs into lane qseqs'
    print '\t'.join(['machine', 'run', 'lane', 'read', 'merged_qseq_file'])
    merged_files = []
    cmds= []
    for machine in sorted(qseqs):
        for run in sorted(qseqs[machine]):
            for lane in sorted(qseqs[machine][run]):
                tiles = [set(qseqs[machine][run][lane][read]) for read in qseqs[machine][run][lane]]
                u = set()
                for t in tiles:
                    u |= t
                for t in tiles:
                    if t != u:
                        raise MergeTileQseqsError("Read sets do not have identical tiles")
                merged_files_dir = os.path.join(dirs['qseqs'], str(machine), str(run), str(lane))
                my.makedir(merged_files_dir)
                for read in sorted(qseqs[machine][run][lane]):
                    merged_file = os.path.join(merged_files_dir,'{}.{}.{}.{}_qseq.txt.gz'.format(machine, run, lane, read))
                    merged_files += [merged_file]
                    print '\t'.join([str(machine), str(run), str(lane), str(read), str(merged_file)])
                    files_to_merge = []
                    compression = set()
                    for tile in sorted(qseqs[machine][run][lane][read]):
                        peek = qseqs[machine][run][lane][read][tile]
                        files_to_merge += [peek['path']]
                        compression.add(peek['compressed'])
                    if True in compression and False in compression:
                        raise MergeTileQseqsError("Compressed and uncompressed qseqs cannot be mixed")
                    if True in compression:
                        cmds += ['gzip -cd {} | cat - | gzip > {}'.format(' '.join(files_to_merge), merged_file)]
                    else:
                        cmds += ['cat {} | gzip > {}'.format(' '.join(files_to_merge), merged_file)]
    job_name = my_name
    print 'Running job {}. Job_info, STDOUT, and STDERR will be in {}'.format(job_name, job_dir)    
    job = my.run_job(cmds, job_name, job_dir, synchronous=False)
    return {'job': job, 'merged_files': merged_files}
    

def main():

    #command line arguments
    parser = argparse.ArgumentParser(parents=[my.default_parser()],
        description = 'list of qseq files for tiles > topdir/qseqs/machine/run/lane/machine.run.lane.read_qseq.txt.gz',
        epilog = 'pypeline.merge_tile_qseqs version 1.0β1 ©2011-2012 Michael Yourshaw all rights reserved')
    parser.add_argument('--metadata', '-m', required=True,
        help='metadata file to map machine, run, lane, barcode to adapters, library, sample and other readgroup info',)
    parser.add_argument('--qseqs', '-q', nargs='+', required=True,
        help='list of qseq files')
    args = parser.parse_args()
    dirs = my.create_default_directory_structure(args.top_dir)
    config = my.get_config(args)
    
    results = run(config, dirs, args.metadata, args.qseqs, hold_jid)
    

if __name__ == "__main__": sys.exit(main())
