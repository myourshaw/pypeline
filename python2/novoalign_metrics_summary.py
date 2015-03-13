#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import re
import my

#-i /scratch1/tmp/myourshaw/gmdjen_redo/metrics/novoalign_metrics/fastq2bam_novoalign_*.e* -o /scratch1/tmp/myourshaw/gmdjen_redo/metrics/novoalign_metrics/consolidated_novoalign_metrics.txt

class NovoalignMetricsSummaryError(Exception): pass

name = 'novoalign_metrics_summary'
version = 1.0
copyright = '©2011-2013 Michael Yourshaw all rights reserved'
run_time = my.localtime_stamp()


def run(input, output):
    Paired_reads_re = re.compile(r'#\s*Paired Reads:\s*(?P<value>\d+)')
    Pairs_Aligned_re = re.compile(r'#\s*Pairs Aligned:\s*(?P<value>\d+)')
    Read_Sequences_re = re.compile(r'#\s*Read Sequences:\s*(?P<value>\d+)')
    Aligned_re = re.compile(r'#\s*Aligned:\s*(?P<value>\d+)')
    Unique_Alignment_re = re.compile(r'#\s*Unique Alignment:\s*(?P<value>\d+)')
    Gapped_Alignment_re = re.compile(r'#\s*Gapped Alignment:\s*(?P<value>\d+)')
    Quality_Filter_re = re.compile(r'#\s*Quality Filter:\s*(?P<value>\d+)')
    Homopolymer_Filter_re = re.compile(r'#\s*Homopolymer Filter:\s*(?P<value>\d+)')
    Elapsed_Time_sec_re = re.compile(r'#\s*Elapsed Time:\s*(?P<value>[0-9.]+)')
    CPU_Time_min_re = re.compile(r'#\s*CPU Time:\s*(?P<value>[0-9.]+)')
    Fragment_header_re = re.compile(r'#\s+From\s+To\s+Count')
    Fragment_re = re.compile(r'#\s*(?P<from>\d+)\s*(?P<to>\d+)\s*(?P<count>\d+)')
    Fragment_Summary_re = re.compile(r'#\s*Mean\s*(?P<mean>[0-9.]+),\s*Std Dev\s*(?P<std_dev>[0-9.]+)')
    files = my.unglob(input)
    fragment_lengths = output.rstrip('.txt')+'.fragment_lengths.txt'
    my.makedir(os.path.dirname(output))
    with open(output, 'w') as o, open(fragment_lengths, 'w') as f:
        o.write('#Path\tFile\tPaired_reads\tPairs_Aligned\tRead_Sequences\tAligned\tUnique_Alignment\tGapped_Alignment\tQuality_Filter\tHomopolymer_Filter\tElapsed_Time_sec\tCPU_Time_min\tMean_Fragment_Size\tStdDev_Fragment_Size\n')
        f.write('#Path\tFile\tFrom\tTo\tCount\n')
        for file in files:
            data = dict()
            with open(file) as input:
                fragments = False
                for line in input:
                    if not fragments:
                        m = Paired_reads_re.match(line)
                        if m:
                            data['Paired_reads'] = m.group('value')
                            continue
                        m = Pairs_Aligned_re.match(line)
                        if m:
                            data['Pairs_Aligned'] = m.group('value')
                            continue
                        m = Read_Sequences_re.match(line)
                        if m:
                            data['Read_Sequences'] = m.group('value')
                            continue
                        m = Aligned_re.match(line)
                        if m:
                            data['Aligned'] = m.group('value')
                            continue
                        m = Unique_Alignment_re.match(line)
                        if m:
                            data['Unique_Alignment'] = m.group('value')
                            continue
                        m = Gapped_Alignment_re.match(line)
                        if m:
                            data['Gapped_Alignment'] = m.group('value')
                            continue
                        m = Quality_Filter_re.match(line)
                        if m:
                            data['Quality_Filter'] = m.group('value')
                            continue
                        m = Homopolymer_Filter_re.match(line)
                        if m:
                            data['Homopolymer_Filter'] = m.group('value')
                            continue
                        m = Elapsed_Time_sec_re.match(line)
                        if m:
                            data['Elapsed_Time_sec'] = m.group('value')
                            continue
                        m = CPU_Time_min_re.match(line)
                        if m:
                            data['CPU_Time_min'] = m.group('value')
                            continue
                        if Fragment_header_re.match(line):
                            fragments = True
                            continue
                    else: #fragmenmts
                        m = Fragment_re.match(line)
                        if m:
                            f.write(file+'\t'+os.path.basename(file)+'\t'+'\t'.join(m.groups())+'\n')
                            continue
                        else:
                            m = Fragment_Summary_re.match(line)
                            if m:
                                data['Mean_Fragment_Size'] = m.group('mean')
                                data['StdDev_Fragment_Size'] = m.group('std_dev')
                                continue
                if data:
                    o.write('\t'.join((file, os.path.basename(file),data['Paired_reads'], data['Pairs_Aligned'], data['Read_Sequences'], data['Aligned'], data['Unique_Alignment'], data['Gapped_Alignment'], data['Quality_Filter'], data['Homopolymer_Filter'], data['Elapsed_Time_sec'], data['CPU_Time_min'], data['Mean_Fragment_Size'], data['StdDev_Fragment_Size'])) + '\n')
                    
            
        
def main():
    
    #command line arguments
    parser = argparse.ArgumentParser(
        description = 'summarize metrics in Novoalign STDERR files',
        epilog = '{} version {} {}'.format(name, version, copyright))
    parser.add_argument('--input', '-i', nargs='+', required=True,
        help='input Novoalign STDERR file(s)')
    parser.add_argument('--output', '-o', required=True,
        help='output file tab delimeted (Fragment Length Distributions will be in output.fragment_lengths.txt)')
    args = parser.parse_args()

    run(input=args.input, output=args.output)

if __name__ == "__main__": sys.exit(main())
