#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import shutil
from glob import glob
import re
import my
import job

class QseqMetricsMergerError(Exception): pass

#-o /scratch0/tmp/myourshaw/gmd1/metrics/qseq_metrics/tile_qseq_metrics_20110924144100.txt -i /scratch0/tmp/myourshaw/gmd1/jobs/pypeline_20110923180426_EzgO2l/tile_qseq_metrics_l7LiBd/qseq_metrics_20110923180749_N6bHf_/metrics_*
#-o /scratch0/tmp/myourshaw/gmd1/metrics/qseq_metrics/readgroup_qseq_metrics_20110924144100.txt -i /scratch0/tmp/myourshaw/gmd1/jobs/pypeline_20110923180426_EzgO2l/tile_qseq_metrics_l7LiBd/qseq_metrics_20110923180749_N6bHf_/metrics_*

def main():
    import argparse
    parser = argparse.ArgumentParser(
        description = 'Merge qseq metrics for multiple groups of read files into a single file',
        epilog = 'pypeline.qseq_metrics_merger version 1.0β1 ©2011-2012 Michael Yourshaw all rights reserved')
    parser.add_argument('--input', '-i', required=True, nargs='+',
        help='list of temporary metrics files to merge')
    parser.add_argument('--output', '-o', required=True,
        help='output file')
    args = parser.parse_args()

    header_written = False
    with open (args.output, 'w') as o:
        for file in args.input:
            with open(file, 'r') as f:
                for line in f:
                    if line.startswith('#') and not header_written:
                        header_written = True
                        o.write(line)
                    elif not line.startswith('#') and not line.strip() == '':
                        o.write(line)
    for file in args.input:
        os.remove(file)


if __name__ == "__main__": sys.exit(main())
