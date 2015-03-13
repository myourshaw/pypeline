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

#-q /scratch0/tmp/myourshaw/gmd/qseqs/*/*/*/*.*.*.[123].*_qseq.txt.gz
#-q /home/solexa/solexa_datasets/100528_HWUSI-EAS335_61DRC/Data/Intensities/BaseCalls/s_3_1_0001_qseq.txt* /home/solexa/solexa_datasets/100528_HWUSI-EAS335_61DRC/Data/Intensities/BaseCalls/s_4_1_0001_qseq.txt* /home/solexa/solexa_datasets/100608_HWUSI-EAS172_61DRF/Data/Intensities/BaseCalls/s_8_1_0001_qseq.txt* /home/solexa/solexa_datasets/100616_HWUSI-EAS172_61DRE/Data/Intensities/BaseCalls/s_6_1_0001_qseq.txt*
#-q /home/solexa/solexa_datasets/100616_HWUSI-EAS172_61DRE/Data/Intensities/BaseCalls/s_6_1_0001_qseq.txt*
#-q /data/storage-1-00/HiSeq2k.Freimerlab/110929_SN860_0095_B_2011-173_D0A1YACXX/Data/Intensities/BaseCalls/QSEQ/s_[456]_[123]_*_qseq.txt.gz

class QseqInfoError(Exception): pass

def main():
    import argparse
    parser = argparse.ArgumentParser(
        description = 'List machine, run, lane, read, barcode for qseq files',
        epilog = 'pypeline.qseq_info version 1.0β1 ©2011-2012 Michael Yourshaw all rights reserved')
    parser.add_argument('--qseqs', '-q', required=True, nargs='+',
        help='list of qseq files')
    args = parser.parse_args()

    print 'Readgroup info from first read in each qseq file ...'
    
    qseq_files = my.unglob(args.qseqs)
    files_count = len(qseq_files)
    qseqs_count = 0

    qseqs = my.qseq_tiles_reads_hierarchy(qseq_files)
    if not qseqs or len(qseqs) == 0:
        raise QseqInfoError('No qseq files')

    print '\t'.join(['machine', 'run', 'lane', 'tile', 'read', 'read_length', 'qseq'])
    for machine in sorted(qseqs):
        for run in sorted(qseqs[machine]):
            for lane in sorted(qseqs[machine][run]):
                for tile in sorted(qseqs[machine][run][lane]):
                    for read in sorted(qseqs[machine][run][lane][tile]):
                        qseqs_count +=1
                        peek = qseqs[machine][run][lane][tile][read]
                        path = peek['path']
                        read_length = peek['read_length']
                        print '{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(machine, run, lane, tile, read, read_length, path)
    print 'This dataset contains {} files, {} distinct readgroups'.format(files_count, qseqs_count)

if __name__ == "__main__": sys.exit(main())
