#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from itertools import combinations_with_replacement
import argparse
from collections import Counter
try:
    import cyvcf as vcf
except ImportError:
    import vcf
import my


class NhlbiVcfs2DbError(Exception):
    pass

# -i /Volumes/scratch/broad/ExAC_release/release0.2/ExAC.r0.2.sites.vep.vcf.gz


def output_file_names(input):
    if input.endswith('.gz'):
        return (my.swap_ext(input, 'gz', 'alleles.gz'),
                my.swap_ext(input, 'gz', 'genotypes.gz'))
    elif input.endswith('.bz2'):
        return (my.swap_ext(input, 'bz2', 'alleles.bz2'),
                my.swap_ext(input, 'bz2', 'genotypes.bz2'))
    else:
        return (input + '.alleles',
                input + '.genotypes')


def writer(file):
    if file.endswith('.gz'):
        import gzip
        return gzip.open(file, 'w')
    elif file.endswith('.bz2'):
        import bz2
        return bz2.open(file, 'w')
    else:
        return open(file,'w')


def run(input, DEBUG=False):

    allele_file, genotype_file = output_file_names(input)
    allele_file_done, genotype_file_done = map(my.done_file, (allele_file, genotype_file))
    if not DEBUG and all(map(my.file_exists, (allele_file, genotype_file, allele_file_done, genotype_file_done))):
        print('{} already converted. To redo:\nrm {} {}'.format(input, allele_file_done, genotype_file_done))
        return

    counter = Counter()

    per_locus_keys = ['AN', 'AN_AFR', 'AN_AMR', 'AN_EAS', 'AN_FIN', 'AN_NFE', 'AN_OTH', 'AN_SAS', 'AN_Adj',]

    per_allele_keys = ['AC', 'AC_AFR', 'AC_AMR', 'AC_EAS', 'AC_FIN', 'AC_NFE', 'AC_OTH', 'AC_SAS', 'AC_Adj',
                       'AF', 'AF_AFR', 'AF_AMR', 'AF_EAS', 'AF_FIN', 'AF_NFE', 'AF_OTH', 'AF_SAS']

    per_genotype_keys = ['GC', 'GC_AFR', 'GC_AMR', 'GC_EAS', 'GC_FIN', 'GC_NFE', 'GC_OTH', 'GC_SAS']

    allele_file_cols = ['CHROM', 'POS', 'REF', 'allele'] + per_locus_keys + per_allele_keys

    genotype_file_cols = ['CHROM', 'POS', 'REF', 'allele1', 'allele2', 'GC_All'] + per_genotype_keys

    with writer(allele_file) as allele_fh, writer(genotype_file) as genotype_fh:
        allele_fh.write('\t'.join(allele_file_cols)+'\n')
        genotype_fh.write('\t'.join(genotype_file_cols)+'\n')

        print('Reading {}; writing {} and {}'.format(input, allele_file, genotype_file))
        vcf_records = vcf.Reader(open(input, 'r'))
        for r in vcf_records:
            counter['vcf_records_in'] += 1
            if counter['vcf_records_in'] % 10000 == 0:
                print(counter['vcf_records_in'])

            # if counter['vcf_records_in'] == 91:
            #     print 'foo'
            #genotypes and indexes thereunto
            gtx = [g for g in combinations_with_replacement(range(len(r.alleles)), 2)]
            genotypes = [(r.alleles[x[0]], r.alleles[x[1]]) for x in gtx]

            #insert allele count values for homozygous REF
            for p in ['', '_AFR', '_AMR', '_EAS', '_FIN', '_NFE', '_OTH', '_SAS', '_Adj']:
                r.INFO['AC'+p].insert(0, r.INFO['AN'+p] - sum(r.INFO['AC'+p]))

            r.INFO['AC_Hom'].insert(0, r.INFO['AN'] - sum(r.INFO['AC_Adj']))

            for p in ['_AFR', '_AMR', '_EAS', '_FIN', '_NFE', '_OTH', '_SAS',]:
                r.INFO['Hom'+p].insert(0, r.INFO['AN'+p] - sum(r.INFO['Hom'+p] + r.INFO['Het'+p]))

            #insert allele frequency values for homozygous REF
            r.INFO['AF'].insert(0, 1.0 - sum(r.INFO['AF']))

            #add within population allele frequencies
            for p in ['_AFR', '_AMR', '_EAS', '_FIN', '_NFE', '_OTH', '_SAS',]:
                r.INFO['AF'+p] = [0 if r.INFO['AN'+p] == 0
                                  else float(r.INFO['AC'+p][i])/float(r.INFO['AN'+p])
                                  for i in range(len(r.alleles))]

            #combine Hom and Het individuals' genotype counts into one list for population independent
            #  and one list per population
            r.INFO['GC'] = []
            het_start, het_end = 0, len(r.alleles) - 1
            for i in range(len(r.alleles)):
                r.INFO['GC'].append(r.INFO['AC_Hom'][i])
                r.INFO['GC'] += r.INFO['AC_Het'][het_start:het_end]
                het_start = het_end
                het_end += len(r.alleles) - (i+2)

            r.INFO['GC_All'] = sum(r.INFO['GC'])

            for p in ['_AFR', '_AMR', '_EAS', '_FIN', '_NFE', '_OTH', '_SAS']:
                r.INFO['GC'+p] = []
                het_start, het_end = 0, len(r.alleles) - 1
                for i in range(len(r.alleles)):
                    r.INFO['GC'+p].append(r.INFO['Hom'+p][i])
                    r.INFO['GC'+p] += r.INFO['Het'+p][het_start:het_end]
                    het_start = het_end
                    het_end += len(r.alleles) - (i+2)

            #write allele file records
            for i in range(len(r.alleles)):
                allele_data = [r.CHROM, r.POS, r.REF, r.alleles[i]]
                allele_data += [r.INFO[d] for d in per_locus_keys]
                allele_data += [r.INFO[d][i] for d in per_allele_keys]
                allele_fh.write('\t'.join(map(str, allele_data))+'\n')
                counter['allele_records_out'] += 1
                # if DEBUG:
                #     print allele_data

            #write genotype file records
            for i in range(len(gtx)):
                allele1x, allele2x = gtx[i]
                genotype_data = [r.CHROM, r.POS, r.REF, r.alleles[allele1x], r.alleles[allele2x], r.INFO['GC_All']]
                genotype_data += [r.INFO[d][i] for d in per_genotype_keys]
                genotype_fh.write('\t'.join(map(str, genotype_data))+'\n')
                counter['genotype_records_out'] += 1
                # if DEBUG:
                #     print genotype_data

            # if DEBUG and counter['vcf_records_in'] >= 1000:
            #     break

    with open(allele_file_done, 'w'):
        pass
    with open(genotype_file_done, 'w'):
        pass
    print('ExAC to database table files done.')
    print('vcf records input: {}'.format(counter['vcf_records_in']))
    print('allele records output: {}'.format(counter['allele_records_out']))
    print('genotype records output: {}'.format(counter['genotype_records_out']))


def main():

    parser = argparse.ArgumentParser(
        description = 'create a database/excel-friendly vcf file by expanding alleles and genotypes',
        epilog = 'pypeline.exac2db version 1.0β1 ©2011-2014 Michael Yourshaw all rights reserved')
    parser.add_argument('--input', '-i', required=True,
        help='input vcf[.gz] files')
    parser.add_argument('--DEBUG', default=False, action='store_true',
        help='input vcf[.gz] files')
    args = parser.parse_args()

    run(input=args.input, DEBUG=args.DEBUG)


if __name__ == "__main__": sys.exit(main())
