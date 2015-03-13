#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import argparse
import gzip
import re
from collections import Counter
import my

class FastqRemoveBarcodeMismatchesError(Exception): pass

#-i /home/vchang/HiSeq2k.Xinmin/Sample_GAR30/GAR30_CGATGT_L001_R1_009.fastq.gz,CGATGT -d /scratch1/tmp/myourshaw/sarc/fastqs

def run(input, output_dir=None, delimiter=','):
    if output_dir:
        my.makedir(output_dir)
    #@HWI-ST0860:79:C03YLACXX:5:1101:1395:2165 1:N:0:ACTTGA
    #@HWI-ST0860:79:C03YLACXX:5:1101:1395:2165 1:N:0:
    sequence_identifier_rx = re.compile(r'^@(?P<instrument>[a-zA-Z0-9_-]+):(?P<run>[0-9.eE-]+):(?P<flowcell>[a-zA-Z0-9]+):(?P<lane>[0-9.eE-]+):(?P<tile>[0-9.eE-]+):(?P<x_pos>[0-9.eE-]+):(?P<y_pos>[0-9.eE-]+) (?P<read>[0-9.eE-]+):(?P<is_filtered>[nyNY]+):(?P<control_number>[0-9]+)(:(?P<index_sequence>[ACGTN]+))?$', re.I)
    sequence_rx = re.compile(r'^[ACGTN.-]+$', re.I)
    for i in input:
        fastq_barcode = i.split(delimiter)
        fastq = fastq_barcode[0]
        expected_barcode = fastq_barcode[1] if len(fastq_barcode) > 1 else None
        if not my.file_exists(fastq):
            raise FastqRemoveBarcodeMismatchesError('File {} does not exist.'.format(fastq))
        counts = Counter()
        with my.open_gz_or_text(fastq) as f:
            line_number = 0
            copy = False
            if expected_barcode != None:
                if output_dir and output_dir != os.path.dirname(fastq):
                    output_fastq = os.path.join(output_dir, os.path.basename(fastq))
                else:
                    output_fastq = fastq+'.single_barcode.fastq.gz'
                out = gzip.open(output_fastq, 'w')
            for line in f:
                line_number += 1
                if line .startswith('@'):
                    sequence_identifier = line.rstrip('\n')
                    match = sequence_identifier_rx.match(sequence_identifier)
                    if not match:
                        raise FastqRemoveBarcodeMismatchesError ('"{}" in file {} line {} is not a proper fastq sequence identifier'.format(sequence_identifier, fastq, line_number))
                    peek = match.groupdict('')
                    peek['machine'] = peek['instrument']
                    peek['barcode'] = peek['index_sequence']
                    counts[peek['index_sequence']] += 1
                    copy = expected_barcode and peek['index_sequence'] == expected_barcode
                if copy:
                    out.write(line)
            print('{}:'.format(fastq))
            for c in sorted(counts.items()):
                print('\t{}: {}'.format(c[0], c[1]))
        out.close()

def main():

    #command line arguments
    parser = argparse.ArgumentParser(
        description = 'count barcodes and remove mismatches from fastq files',
        epilog = 'fastq_remove_barcode_mismatches 1.0β1 ©2011-2014 Michael Yourshaw all rights reserved')
    parser.add_argument('--input', '-i', required=True, nargs='+',
        help='list of fastq files[,optional expected barcode] (required; example: /dir/foo.fastq.gz,ACTGTG /dir/bar.fastq.gz,CCAGGG)')
    parser.add_argument('--output_dir', '-d',
        help='optional output directory for cleaded up fastqs')
    parser.add_argument('--delimiter', default = ',',
        help='delimiter between fastq and barcode (default: comma)')
    args = parser.parse_args()

    run(input=args.input, output_dir=args.output_dir, delimiter=args.delimiter)
    

if __name__ == "__main__": sys.exit(main())
