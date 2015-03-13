#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import argparse
from ConfigParser import SafeConfigParser #configparser in python 3
import gzip
import re
from collections import Counter, OrderedDict
from vcf import Reader
from string import maketrans
from warnings import warn
import my
import job
import sql_columns
import vax_merge_post_process


class Vcf2GtxError(Exception): pass

name = 'pypeline.vcf2gtx'
version = 75
copyright = '©2011-2014 Michael Yourshaw all rights reserved'
run_time = my.localtime_stamp()


# -i /share/apps/myourshaw/vax/test/t.vcf
# -i /share/apps/myourshaw/vax/test/multialleic_test.vcf


def plx(a1, a2):
    """complute index in PL field for a genotype"""
    return (a2*(a2+1)/2)+a1


def run(input, output=None, force=False):
    """extracts per-sample data and statistics from a VCF file and writes to a GTX file"""

    my_name = 'vcf2gtx'

    my.print_log('{} {} version {}. Copyright {}.'.format(run_time, name, version, copyright))
    my.print_log("""IMPORTANT:
This application will create a hidden file with a name like '.mydata.vcf.gtx.done'.
After a restart, files with an associated the .done file will not be recreated.
To force overwrite a an existing file, delete the associated .done file or use the --force parameter.
""")

    # output directories and files
    if output:
        gtx = output
    else:
        gtx = my.swap_ext(input, '.vcf', '.gtx')
    my.makedir(os.path.dirname(gtx))
    gtx_done = my.done_file(gtx)
    if my.file_exists(gtx) and my.file_exists(gtx_done) and not force:
        warn('{} already created. Skipping.'.format(gtx))
        return
    try:
        os.remove(gtx_done)
    except OSError:
        pass

    gtx_header_lines_top = [
        '## GTX file created at {} by {} version {}.'.format(run_time, name, version),
        '## Copyright {}.'.format(copyright),
        '## GTX file has one line per sample with allele and genotype statistics for each sample',
        '## GTX input VCF file {}'.format(input),
        '## GTX output GTX file {}'.format(gtx),
    ]
    gtx_header_lines_bottom = [
        '## AD1 ALLELE1 depth',
        '## AD2 ALLELE2 depth',
        '## ADREF REF allele depth',
        '## ALLELE1 The first allele of this genotype',
        '## ALLELE2 The second allele of this genotype',
        '## DP Approximate read depth (reads with MQ=255 or with bad mates are filtered)',
        '## GQ Genotype Quality (= second lowest PL, max 99; the PL of this sample\'s genotype is scaled to 0)',
        '## GT1 The index in the GT field of the first allele of this genotype',
        '## GT2 The index in the GT field of the second allele of this genotype',
        '## PHASED Whether the genotype is phased',
        '## PL The phred-scaled genotype likelihoods rounded to the closest integer',
        '## PLA1A1 Normalized, Phred-scaled likelihood for genotype ALLELE1/ALLELE1.',
        '## PLA1A2 Normalized, Phred-scaled likelihood for genotype ALLELE1/ALLELE2.',
        '## PLA2A2 Normalized, Phred-scaled likelihood for genotype ALLELE2/ALLELE2.',
        '## PLRA1 Normalized, Phred-scaled likelihood for genotype REF/ALLELE1.',
        '## PLRA2 Normalized, Phred-scaled likelihood for genotype REF/ALLELE2.',
        '## PLRR Normalized, Phred-scaled likelihood for genotype REF/REF.',
        '## SAMPLE Sample ID taken from genotype column heading in VCF file',
        '## VCFA1A1C number of samples with ALLELE1/ALLELE1 in called genotypes of the VCF samples',
        '## VCFA1A2C number of samples with ALLELE1/ALLELE2 in called genotypes of the VCF samples',
        '## VCFA2A2C number of samples with ALLELE2/ALLELE2 in called genotypes of the VCF samples',
        '## VCFAC total number of alleles in called genotypes of the VCF samples',
        '## VCFAC1 number of ALLELE1s in called genotypes of the VCF samples',
        '## VCFAC2 number of ALLELE2s in called genotypes of the VCF samples',
        '## VCFACREF number of REF alleles in called genotypes of the VCF samples',
        '## VCFRA1C number of samples with REF/ALLELE1 in called genotypes of the VCF samples',
        '## VCFRA2C number of samples with REF/ALLELE2 in called genotypes of the VCF samples',
        '## VCFRRC number of samples with REF/REF in called genotypes of the VCF samples',
        '## VCFSAMP number of called genotypes of the VCF samples',
    ]

    gtx_cols = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'SAMPLE', 'GT', 'GT1', 'GT2', 'ALLELE1',
                'ALLELE2', 'PHASED', 'ADREF', 'AD1', 'AD2', 'DP', 'GQ', 'PL', 'PLRR', 'PLRA1', 'PLA1A1', 'PLRA2',
                'PLA1A2', 'PLA2A2', 'VCFSAMP', 'VCFAC', 'VCFACREF', 'VCFAC1', 'VCFAC2', 'VCFRRC', 'VCFRA1C',
                'VCFA1A1C', 'VCFRA2C', 'VCFA1A2C', 'VCFA2A2C',]

    stats = Counter()

    my.print_log('Starting vcf2gtx')
    for l in gtx_header_lines_top:
        my.print_log(l.lstrip('# '))

    # open the VCF file with vcf.Reader
    v = Reader(open(input, 'r'))

    # the INFO field will be expanded to columns
    info_cols = ['i_' + c for c in v.infos.keys()]
    output_cols = gtx_cols + info_cols

    # metadata
    with open(gtx, 'w') as gtxfh:
        # write header lines
        gtxfh.write('\n'.join(gtx_header_lines_top + v._header_lines + gtx_header_lines_bottom)+'\n')
        gtxfh.write('#'+'\t'.join(output_cols)+'\n')

        # process each data record
        for r in v:
            stats['vcf_records'] += 1
            if stats['vcf_records'] % 10000 == 0:
                my.print_log('{}:{}\t{} records'.format(r.CHROM, r.POS, stats['vcf_records']))

            # DEBUG
            # if vcf_record_count > 100: break
            # continue

            # q is the set of samples with called genotypes
            q = filter(lambda x: x.called, r.samples)

            # expand data to one record per sample
            for s in r.samples:
                data = OrderedDict([(c, '') for c in output_cols])
                data['CHROM'] = r.CHROM
                data['POS'] = r.POS
                data['ID'] = ';'.join(r.ID) if isinstance(r.ID, (list, tuple)) else r.ID if r.ID else '.'
                data['REF'] = r.REF
                data['ALT'] = ','.join(map(str, r.ALT))
                data['QUAL'] = r.QUAL if r.QUAL else '.'
                data['FILTER'] = ';'.join(r.FILTER) if r.FILTER else '.'
                # data['INFO'] = r.INFO
                for k in r.INFO.keys():
                    if isinstance(r.INFO[k], bool) and r.INFO[k]:
                        data['i_'+k] = 1
                    elif isinstance(r.INFO[k], (list, tuple)):
                        data['i_'+k] = ','.join(map(str, r.INFO[k]))
                    elif data['i_'+k] is None:
                        data['i_'+k] = ''
                    else:
                        data['i_'+k] = r.INFO[k]
                # data['FORMAT'] = r.FORMAT
                data['SAMPLE'] = s.sample
                data['GT'] = s.gt_nums if s.called else './.'
                data['PHASED'] = 1 if s.phased else ''

                # statistics for all samples
                data['ADREF'] = s.data.AD[0]
                data['DP'] = s.data.DP
                data['GQ'] = s.data.GQ if s.data.GQ is not None else ''
                data['VCFSAMP'] = r.num_called
                data['VCFAC'] = 2*len(q)
                data['VCFACREF'] = (sum(s.gt_alleles[0] == '0' for s in q) + sum(s.gt_alleles[1] == '0' for s in q))

                # statistics for samples with called genotypes
                if s.called:
                    a1s, a2s = s.gt_alleles # e.g., ['0','1']
                    a1x, a2x = map(int, s.gt_alleles) # e.g., [0,1]
                    data['GT1'], data['GT2'] = map(int, s.gt_alleles)
                    data['ALLELE1'] = r.alleles[a1x]
                    data['ALLELE2'] = r.alleles[a2x]
                    data['AD1'] = s.data.AD[a1x]
                    data['AD2'] = s.data.AD[a2x]
                    data['PL'] = ','.join(map(str, s.data.PL))
                    data['PLRR'] = s.data.PL[plx(0, 0)]
                    data['PLRA1'] = s.data.PL[plx(0, a1x)]
                    data['PLA1A1'] = s.data.PL[plx(a1x, a1x)]
                    data['PLRA2'] = s.data.PL[plx(0, a2x)]
                    data['PLA1A2'] = s.data.PL[plx(a1x, a2x)]
                    data['PLA2A2'] = s.data.PL[plx(a2x, a2x)]
                    foo = (sum([s.gt_alleles[0] == a1s for s in q]) + sum([s.gt_alleles[1] == a1s for s in q]))
                    bar = (sum(s.gt_alleles[0] == a1s for s in q) + sum(s.gt_alleles[1] == a1s for s in q))
                    data['VCFAC1'] = (sum(s.gt_alleles[0] == a1s for s in q) + sum(s.gt_alleles[1] == a1s for s in q))
                    data['VCFAC2'] = (sum(s.gt_alleles[0] == a2s for s in q) + sum(s.gt_alleles[1] == a2s for s in q))
                    data['VCFRRC'] = sum(s.gt_alleles[0] == '0' and s.gt_alleles[1] == '0' for s in q)
                    data['VCFRA1C'] = sum(s.gt_alleles[0] == '0' and s.gt_alleles[1] == a1s for s in q)
                    data['VCFA1A1C'] = sum(s.gt_alleles[0] == a1s and s.gt_alleles[1] == a1s for s in q)
                    data['VCFRA2C'] = sum(s.gt_alleles[0] == '0' and s.gt_alleles[1] == a2s for s in q)
                    data['VCFA1A2C'] = sum((s.gt_alleles[0] == a1s and s.gt_alleles[1] == a2s)
                                           or (s.gt_alleles[0] == a2s and s.gt_alleles[1] == a1s) for s in q)
                    data['VCFA2A2C'] = sum(s.gt_alleles[0] == a2s and s.gt_alleles[1] == a2s for s in q)

                # write per-sample GTX output record
                output_record = '\t'.join(map(str, data.values()))+'\n'
                gtxfh.write(output_record)
                stats['gtx_records'] += 1

   #done
    with open(gtx_done, 'w'):
        pass
    my.print_log('Input vcf {} had {} records.'.format(input, stats['vcf_records']))
    my.print_log('Output gtx file {} had {} records.'.format(gtx, stats['gtx_records']))
    my.print_log('vcf2gtx done.')





def main():

    #command line arguments
    parser = argparse.ArgumentParser(parents=[my.default_parser()],
        description = 'create a version of a vcf file flattened by sample with statistics',
        epilog = '{} version {} {}'.format(name, version, copyright))
    parser.add_argument('--input', '-i', '-v', required=True,
        help='input vcf file')
    parser.add_argument('--output', '-o', '-g',
        help='output gtx file (default: <input>.gtx')
    parser.add_argument('--force', action='store_true',
        help='overwrite existing gtx file')

    args = parser.parse_args()

    run(input=args.input, output=args.output, force=args.force)


if __name__ == "__main__": sys.exit(main())
