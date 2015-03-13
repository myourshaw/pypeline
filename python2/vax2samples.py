#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import glob
import argparse
import re
import gzip
import my

class Vax2Samples(Exception): pass

#-i /scratch1/tmp/myourshaw/hlee/vcfs/hlee55gmd28bipolar20ocdtic18.analysis_ready.vcf.vax -s 121 -o /scratch1/tmp/myourshaw/hlee/vcfs/hlee55gmd28bipolar20ocdtic18.analysis_ready.vcf.vax.samples_variants_conseq7.txt.gz


def main():

    additional_columns = ['SAMPLE', 'ALLELE1', 'PHASE', 'ALLELE2', 'ZYGOSITY', 'compound_het', 'GT', 'GT_INFO']

    #command line arguments
    parser = argparse.ArgumentParser(
        description = 'expand VAX output file by samples, and filter by Consequence_rank',
        epilog = 'pypeline.vax2samples version 1.0β1 ©2011-2012 Michael Yourshaw all rights reserved')
    parser.add_argument('--input', '-i', required=True,
        help='input vep[.gz] or vax[.gz] file')
    parser.add_argument('--output', '-o',
        help='output file (default:<input>.samples.txt)')
    parser.add_argument('--sample_column_count', '-s', required=True, type=int,
        help='number of genotype columns')
    parser.add_argument('--consequence_threshhold', '-t', type=int, default=9,
        help='only output records where Consequence_rank <= consequence_threshhold (default: 9)')
    parser.add_argument('--output_uncalled_genotypes', action='store_true', default=False,
        help='output uncalled (./.) genotypes (default: False)')
    parser.add_argument('--output_homozygous_reference_genotypes', action='store_true', default=False,
        help='output homozygous_reference (0/0) genotypes (default: False)')
    parser.add_argument('--uncompressed', '-u', action='store_true', default=False,
        help='do not gzip output (default: False)')
    args = parser.parse_args()
    
    if not args.output:
        args.output = args.input+'.samples.txt'+'.gz' if not args.uncompressed else ''
    if args.uncompressed and args.output.endswith('.gz'):
        args.output = args.output[:-3]
    if not args.uncompressed and not args.output.endswith('.gz'):
        args.output = args.output+'.gz'

    columns_in = {}

    gt_re = re.compile(r'(?P<allele1>[0-9.]+)(?P<phase>[/|])(?P<allele2>[0-9.]+)', re.I)
    
    with my.open_gz_or_text(args.input) as i, gzip.open(args.output, 'w') if not args.uncompressed else open(args.output, 'w') as o:
        line_count = 0
        line_out_count = 0
        for line in i:
            #if line_out_count > 10000: break
            line_count += 1
            if line.startswith('##'):
                continue
            elif line.startswith('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'):
                col_list = line.rstrip('\n').split('\t')
                columns_in = {col_list[x]: x for x in range(len(col_list))}
                sample_columns = col_list[9:9+args.sample_column_count]
                columns_out = col_list + additional_columns
                o.write('\t'.join(columns_out)+'\n')
                continue
            elif not columns_in:
                raise Vax2Samples('data encountered before column names at line {}'.format(line_count))
            else:
                data_list = line.rstrip('\n').split('\t')
                data_in = {k: data_list[columns_in[k]] for k in columns_in.keys()}
                if my.is_int(data_in['Consequence_rank']) and int(data_in['Consequence_rank']) > args.consequence_threshhold:
                    continue
                data_in['compound_het'] = ''
                alt_list = data_in['ALT'].rstrip(',').split(',')
                ref = data_in['REF']
                these_alleles = [ref] + alt_list
                for alt in alt_list:
                    data_in['ALT'] = alt
                    for sample in sample_columns:
                        data_in['SAMPLE'] = sample
                        data_in['GT_INFO'] = data_in[sample]
                        m = gt_re.match(data_in[sample])
                        allele1x, phase, allele2x = m.groups() if m else ('.','/','.')
                        data_in['GT'] = '{}{}{}'.format(allele1x, phase, allele2x)
                        allele1 = these_alleles[int(allele1x)] if allele1x != '.' else '.'
                        data_in['ALLELE1'] = allele1
                        allele2 = these_alleles[int(allele2x)] if allele2x != '.' else '.'
                        data_in['ALLELE2'] = allele2
                        data_in['PHASE'] = phase
                        data_in['ZYGOSITY'] = '' if allele1 == '.' or allele2 == '.' else '0' if allele1 == ref and allele2 == ref else '2' if allele1 == alt and allele2 == alt else '1' if (allele1 == ref and allele2 == alt) or (allele1 == alt and allele2 == ref) else str(int(allele1x)+int(allele2x))
                        if not args.output_homozygous_reference_genotypes and data_in['ZYGOSITY'] == '0':
                            continue
                        if not args.output_uncalled_genotypes and data_in['ZYGOSITY'] == '':
                            continue
                        data_out = [data_in[c] for c in columns_out]
                        o.write('\t'.join(data_out)+'\n')
                        line_out_count += 1
    print 'done, wrote {} lines'.format(line_out_count)


if __name__ == "__main__": sys.exit(main())
