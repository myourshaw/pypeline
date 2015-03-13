#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import argparse
import collections
import my

#-q /scratch0/tmp/myourshaw/gmd/qseqs/HWI-ST0860/79/1/HWI-ST0860.79.1.1.ACAGTGA_qseq.txt.gz /scratch0/tmp/myourshaw/gmd/qseqs/HWI-ST0860/79/1/HWI-ST0860.79.1.2.ACAGTGA_qseq.txt.gz /scratch0/tmp/myourshaw/gmd/qseqs/HWI-ST0860/79/1/HWI-ST0860.79.1.3.ACAGTGA_qseq.txt.gz -o /scratch0/tmp/myourshaw/gmd/jobs/qseq_metrics_20110913133544_grxSYL/metrics_WD2ZPO

def main():
    parser = argparse.ArgumentParser(
        description = 'Calculate pass/fail metrics for qseq files and validate that rows match in each file of a group of reads',
        epilog = 'pypeline.qseq_metrics_calculator version 1.0β1 ©2011-2012 Michael Yourshaw all rights reserved')
    parser.add_argument('--qseqs', '-q', nargs='+',
        help='list of 1 to 3 qseq files that are expected to have corresponding reads; 4+ ignored')
    parser.add_argument('--output', '-o', required=True,
        help='output file')
    args = parser.parse_args()
    
    #calculate_qseq_metrics
    a = len(args.qseqs) > 0
    b = len(args.qseqs) > 1
    c = len(args.qseqs) > 2
    if a and my.qseq_peek(args.qseqs[0]):
        aok = True
        a_qseq = my.open_gz_or_text(args.qseqs[0])
    else:
        aok = False
    if b and my.qseq_peek(args.qseqs[1]):
        bok = True
        b_qseq = my.open_gz_or_text(args.qseqs[1])
    else:
        bok = False
    if c and my.qseq_peek(args.qseqs[2]):
        cok = True
        c_qseq = my.open_gz_or_text(args.qseqs[2])
    else:
        cok = False
    
    results = {'a_file': args.qseqs[0] if a else '', 'b_file': args.qseqs[1] if b else '', 'c_file': args.qseqs[2] if c else '',
                'a_pass': 0, 'a_fail': 0, 'a_count': 0, 'a_pct_pass': 0, 'a_length': 0,
                'b_pass': 0, 'b_fail': 0, 'b_count': 0, 'b_pct_pass': 0, 'b_length': 0,
                'c_pass': 0, 'c_fail': 0, 'c_count': 0, 'c_pct_pass': 0, 'c_length': 0,
                'ab_first_mismatch_line': 0, 'ac_first_mismatch_line': 0, 'bc_first_mismatch_line': 0}
    line_count = 0

    while aok or bok or cok:
        if aok:
            linea = a_qseq.readline()
            if linea:
                a = linea.rstrip('\n').split('\t')
                results['a_count'] += 1
                if len(a) > 10:
                    results['a_length'] = max((results['a_length'],len(a[7])))
                    if a[10] == '1':
                        results['a_pass'] += 1
                    else:
                        results['a_fail'] += 1
            else:
                aok = False
        if bok:
            lineb = b_qseq.readline()
            if lineb:
                b = lineb.rstrip('\n').split('\t')
                results['b_count'] += 1
                if len(b) > 10:
                    results['b_length'] = max((results['b_length'],len(a[7])))
                    if b[10] == '1':
                        results['b_pass'] += 1
                    else:
                        results['b_fail'] += 1
            else:
                bok = False
        if cok:
            linec = c_qseq.readline()
            if linec:
                c = linec.rstrip('\n').split('\t')
                results['c_count'] += 1
                if len(c) > 10:
                    results['c_length'] = max((results['c_length'],len(a[7])))
                    if c[10] == '1':
                        results['c_pass'] += 1
                    else:
                        results['c_fail'] += 1
            else:
                cok = False
        line_count += 1
        if a and b and not results['ab_first_mismatch_line'] and (a[:6] != b[:6]):
            results['ab_first_mismatch_line'] = line_count
        if a and c and not results['ac_first_mismatch_line'] and (a[:6] != c[:6]):
            ac_first_mismatch_line = line_count
        if b and c and not results['bc_first_mismatch_line'] and (b[:6] != c[:6]):
            bc_first_mismatch_line = line_count
    if a and results['a_count'] > 0:
        results['a_pct_pass'] = float(results['a_pass'])/float(results['a_count'])
    if b and results['b_count'] > 0:
        results['b_pct_pass'] = float(results['b_pass'])/float(results['b_count'])
    if c and results['c_count'] > 0:
        results['c_pct_pass'] = float(results['c_pass'])/float(results['c_count'])

    #save results in output file, for upstream cat with other result sets
    result_fields =['a_file', 'b_file', 'c_file',
                    'a_pass', 'a_fail', 'a_count', 'a_pct_pass', 'a_length',
                    'b_pass', 'b_fail', 'b_count', 'b_pct_pass', 'b_length',
                    'c_pass', 'c_fail', 'c_count', 'c_pct_pass', 'c_length',
                    'ab_first_mismatch_line', 'ac_first_mismatch_line', 'bc_first_mismatch_line']
    with open(args.output, 'w') as o:
        #header
        o.write('#{}\n'.format('\t'.join(result_fields)))
        #results
        o.write('{{{}}}\n'.format('}\t{'.join(result_fields)).format(**results))


if __name__ == "__main__": sys.exit(main())
