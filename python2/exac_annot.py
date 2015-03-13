#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = 'myourshaw'

import sys
import os
import argparse
import sqlite3
from collections import Counter, namedtuple

from my import is_int, done_file, localtime_stamp, open_gz_or_text, swap_ext


#-i /Volumes/scratch/vf_20140828/vcfs/gatk/vf_20140828.hc.analysis_ready.vcf.vax.head -d /Volumes/scratch/broad/ExAC_release/release0.2/exac
#-i /Volumes/scratch/vf_20140828/vcfs/gatk/vf_20140828.hc.analysis_ready.vcf.gtx.head --allele1_col --allele2_col -d /Volumes/scratch/broad/ExAC_release/release0.2/exac


class ExACAnnotError(Exception): pass


Key = namedtuple('Key', ['CHROM', 'POS', 'REF', 'ALT', 'ALLELE1', 'ALLELE2'])

def make_output_file_name(input, ext='exac'):
    if input.endswith('.gz'):
        return swap_ext(input, 'gz', ext+'.txt.gz')
    elif input.endswith('.txt'):
        return swap_ext(input, '.txt', ext+'.txt')
    else:
        return input + '.'+ext+'.txt'


def writer(name):
    if name.endswith('.gz'):
        import gzip
        return gzip.open(name, 'w')
    else:
        return open(name,'w')

metadata_text ='''##INFO=<ID=AC_Adj1,Number=1,Type=Integer,Description="Adjusted allele1 count">
##INFO=<ID=AC_Adj2,Number=1,Type=Integer,Description="Adjusted allele2 count">
##INFO=<ID=AC_AFR1,Number=1,Type=Integer,Description="American American allele1 count">
##INFO=<ID=AC_AFR2,Number=1,Type=Integer,Description="African/African American allele2 count">
##INFO=<ID=AC_AMR1,Number=1,Type=Integer,Description="American allele1 count">
##INFO=<ID=AC_AMR2,Number=1,Type=Integer,Description="American allele2 count">
##INFO=<ID=AC_EAS1,Number=1,Type=Integer,Description="East Asian allele1 count">
##INFO=<ID=AC_EAS2,Number=1,Type=Integer,Description="East Asian allele2 count">
##INFO=<ID=AC_FIN1,Number=1,Type=Integer,Description="Finnish allele1 count">
##INFO=<ID=AC_FIN2,Number=1,Type=Integer,Description="Finnish allele2 count">
##INFO=<ID=AC_NFE1,Number=1,Type=Integer,Description="Non-Finnish European allele1 count">
##INFO=<ID=AC_NFE2,Number=1,Type=Integer,Description="Non-Finnish European allele2 count">
##INFO=<ID=AC_OTH1,Number=1,Type=Integer,Description="Other allele1 count">
##INFO=<ID=AC_OTH2,Number=1,Type=Integer,Description="Other allele2 count">
##INFO=<ID=AC_SAS1,Number=1,Type=Integer,Description="South Asian allele1 count">
##INFO=<ID=AC_SAS2,Number=1,Type=Integer,Description="South Asian allele2 count">
##INFO=<ID=AC1,Number=1,Type=Integer,Description="Allele count in genotypes, for allele1">
##INFO=<ID=AC2,Number=1,Type=Integer,Description="Allele count in genotypes, for allele2">
##INFO=<ID=AF_AFR1,Number=1,Type=Float,Description="American American Allele Frequency, for allele1">
##INFO=<ID=AF_AFR2,Number=1,Type=Float,Description="American American Allele Frequency, for allele2">
##INFO=<ID=AF_AMR1,Number=1,Type=Float,Description="American Allele Frequency, for allele1">
##INFO=<ID=AF_AMR2,Number=1,Type=Float,Description="American Allele Frequency, for allele2">
##INFO=<ID=AF_EAS1,Number=1,Type=Float,Description="East Asian Allele Frequency, for allele1">
##INFO=<ID=AF_EAS2,Number=1,Type=Float,Description="East Asian Allele Frequency, for allele2">
##INFO=<ID=AF_FIN1,Number=1,Type=Float,Description="Finnish Allele Frequency, for allele1">
##INFO=<ID=AF_FIN2,Number=1,Type=Float,Description="Finnish Allele Frequency, for allele2">
##INFO=<ID=AF_NFE1,Number=1,Type=Float,Description="Finnish Allele Frequency, for allele1">
##INFO=<ID=AF_NFE2,Number=1,Type=Float,Description="Finnish Allele Frequency, for allele2">
##INFO=<ID=AF_OTH1,Number=1,Type=Float,Description="Other Allele Frequency, for allele1">
##INFO=<ID=AF_OTH2,Number=1,Type=Float,Description="Other Allele Frequency, for allele2">
##INFO=<ID=AF_SAS1,Number=1,Type=Float,Description="South Asian Allele Frequency, for allele1">
##INFO=<ID=AF_SAS2,Number=1,Type=Float,Description="South Asian Allele Frequency, for allele2">
##INFO=<ID=AF1,Number=1,Type=Float,Description="Allele Frequency, for allele1">
##INFO=<ID=AF2,Number=1,Type=Float,Description="Allele Frequency, for allele2">
##INFO=<ID=allele1,Number=1,Type=String,Description="The first alle of this genotype (will be REF if genotype includes REF, otherwise one of the ALTs)">
##INFO=<ID=allele2,Number=1,Type=String,Description="The second alle of this genotype (will be REF if genotype is homozygous REF, otherwise one of the ALTs)">
##INFO=<ID=AN_Adj,Number=1,Type=Integer,Description="Adjusted Chromosome Count">
##INFO=<ID=AN_AFR,Number=1,Type=Integer,Description="African/African American Chromosome Count">
##INFO=<ID=AN_AMR,Number=1,Type=Integer,Description="American Chromosome Count">
##INFO=<ID=AN_EAS,Number=1,Type=Integer,Description="East Asian Chromosome Count">
##INFO=<ID=AN_FIN,Number=1,Type=Integer,Description="Finnish Chromosome Count">
##INFO=<ID=AN_NFE,Number=1,Type=Integer,Description="Non-Finnish European Chromosome Count">
##INFO=<ID=AN_OTH,Number=1,Type=Integer,Description="Other Chromosome Count">
##INFO=<ID=AN_SAS,Number=1,Type=Integer,Description="South Asian Chromosome Count">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=GC_AFR,Number=1,Type=Integer,Description="African/African American Genotype Count">
##INFO=<ID=GC_AMR,Number=1,Type=Integer,Description="American  Genotype Count">
##INFO=<ID=GC_EAS,Number=1,Type=Integer,Description="East Asian  Genotype Count">
##INFO=<ID=GC_FIN,Number=1,Type=Integer,Description="Finnish  Genotype Count">
##INFO=<ID=GC_NFE,Number=1,Type=Integer,Description="Non-Finnish European  Genotype Count">
##INFO=<ID=GC_OTH,Number=1,Type=Integer,Description="Other  Genotype Count">
##INFO=<ID=GC_SAS,Number=1,Type=Integer,Description="South Asian  Genotype Count">
##INFO=<ID=GC,Number=1,Type=Integer,Description="Genotype Count for this genotype">
'''

def run(input, database_file, delim='\t', chrom_col='CHROM', pos_col='POS', ref_col='REF', alt_col='ALT',
        allele1_col='ALLELE1', allele2_col='ALLELE2', skip_string=None, ignore_lines=0, header_line='#',):

    header_line_number = int(header_line) if is_int(header_line) and int(header_line) > 0 else None
    header_chars = header_line if not is_int(header_line) else None

    output = make_output_file_name(input)
    output_done_file = done_file(output)

    if all(map(os.path.exists, (output, output_done_file))):
        print('Output {} already created. To redo:\nrm {}'.format(os.path.basename(output), output_done_file))
        return

    key_names = Key(chrom_col, pos_col, ref_col, alt_col, allele1_col, allele2_col)

    print('{} Creating {} with ExAC annotations'.format(localtime_stamp(), output))
    counter = Counter()
    header_found = False
    with open_gz_or_text(input) as ifh, writer(output) as ofh, sqlite3.connect(database_file) as db:
        db.row_factory = sqlite3.Row

        ofh.write('##ExAC ANNOT output produced {}\n'.format(localtime_stamp()))

        # get column names from SQLite tables
        # get per locus (counts of individuals) column names from allele1 (allele2's are identical)
        per_locus_keys = ['AN', 'AN_AFR', 'AN_AMR', 'AN_EAS', 'AN_FIN', 'AN_NFE', 'AN_OTH', 'AN_SAS', 'AN_Adj']
        # ignore column names for alleles, which will appear in a differend order in the output
        ignore_exac_columns = list(key_names) + ['allele', 'allele1', 'allele2'] + per_locus_keys
        alleles_column_names = [str(c['name']) for c in db.execute('PRAGMA table_info(alleles);')]
        genotypes_column_names = [str(c['name']) for c in db.execute('PRAGMA table_info(genotypes);')]
        per_locus_column_names = [c for c in alleles_column_names if c in per_locus_keys]
        allele1_column_names = [c+'1' for c in alleles_column_names if c not in ignore_exac_columns]
        allele2_column_names = [c+'2' for c in alleles_column_names if c not in ignore_exac_columns]
        genotype_column_names = [c for c in genotypes_column_names if c not in ignore_exac_columns]
        exac_column_names = (per_locus_column_names
                             + [j for i in zip(allele1_column_names, allele2_column_names) for j in i]
                             + genotype_column_names)

        # read input
        for line in ifh:
            counter['input_lines'] += 1
            if counter['input_lines'] % 10000 == 0:
                print ('{} {}'.format(localtime_stamp(), counter['input_lines']))

            # empty line
            if not bool(line.strip()):
                ofh.write(line)
                continue

            # header line?
            if not header_found:
                if ((header_line_number and header_line_number == counter['input_lines'])
                        or (header_chars and header_chars != '#' and line.startswith(header_chars))
                        or (header_chars and header_chars == '#' and line.startswith(header_chars) and line[1] != '#')):
                    header_found = True

                    # get input column names and keys for db queries
                    input_column_names = line.rstrip('\n').split(delim)
                    input_column_names[0] = input_column_names[0].lstrip('#')

                    # default allele column names
                    if not key_names.ALLELE1:
                        for a in ('allele1', 'allele_1'):
                            if a in map(str.lower, input_column_names):
                                key_names = key_names._replace(ALLELE1 = input_column_names[map(str.lower, input_column_names).index(a)])
                        if not key_names.ALLELE1:
                            key_names = key_names._replace(ALLELE1 = ref_col)
                    if not key_names.ALLELE2:
                        for a in ('allele2', 'allele_2'):
                            if a in map(str.lower, input_column_names):
                                key_names = key_names._replace(
                                    ALLELE2 = input_column_names[map(str.lower, input_column_names).index(a)])
                        if not key_names.ALLELE2:
                            key_names = key_names._replace(ALLELE2 = alt_col)
                    add_allele1_col = key_names.ALLELE1 not in input_column_names
                    add_allele2_col = key_names.ALLELE2 not in input_column_names
                    # are all key column names really in the input header
                    try:
                        key_x = Key(*[map(str.lower, input_column_names).index(x) for x in map(str.lower,key_names)])
                    except TypeError:
                        raise ExACAnnotError('Header {} at line {} does not contain all of {}.'
                                             .format(line, counter['input_lines'], key_names))

                    # add output column names
                    added_output_column_names = ([key_names.ALLELE1] if add_allele1_col else []
                                                + [key_names.ALLELE2] if add_allele2_col else []
                                                + [n if n.lower() not in map(str.lower, input_column_names)
                                                 else n+'_exac' for n in exac_column_names
                                                 if n.lower() not in map(str.lower, key_names)])
                    output_column_names = input_column_names + added_output_column_names

                    # write metadata and header row
                    ofh.write(metadata_text)
                    ofh.write(delim.join(output_column_names) + '\n')
                else:
                    ofh.write(line)
                continue

            # skippable line?
            elif ((skip_string and line.startswith(skip_string))
                or (ignore_lines and counter['input_lines'] <= ignore_lines)):
                ofh.write(line)
                continue

            # data line
            else:
                counter['input_data_records'] += 1
                input_data = line.rstrip('\n').split(delim)

                # get db query keys from input data
                chom, pos, ref, alt, allele1, allele2 = [input_data[x] for x in key_x]

                if ',' in allele1:
                    raise ExACAnnotError('Multiple alleles not permitted in {}. Error in line {} of {}.'
                                         .format(key_names.ALLELE1, counter['input_lines'], input))
                # output annotated record for each genotype
                allele2s = allele2.split(',')
                genotypes = [(allele1, a) for a in allele2s]
                for allele1, allele2 in genotypes:

                    # get locus- and allele-specific data from sqlite db for allele1
                    r = [d for d in db.execute(
                        '''SELECT * FROM alleles WHERE CHROM = ? and POS = ? and REF = ? AND allele = ?;''',
                        (chom, pos, ref, allele1))]
                    if r:
                        per_locus_data = ([allele1] if add_allele1_col else [] + [allele2] if add_allele2_col else [] +
                                          [r[0][str(n)] for n in per_locus_column_names])
                        allele1_data = [r[0][str(n)] for n in alleles_column_names if n not in ignore_exac_columns]
                    else:
                        per_locus_data = ['', ''] + ['' for n in per_locus_column_names]
                        allele1_data = ['' for n in alleles_column_names if n not in ignore_exac_columns]

                    # get allele-specific data from sqlite db for allele2
                    # skip double database query if genotype is homozygous
                    if allele1 == allele2:
                        allele2_data = allele1_data
                    else:
                        r = [d for d in db.execute(
                            '''SELECT * FROM alleles WHERE CHROM = ? and POS = ? and REF = ? AND allele = ?;''',
                            (chom, pos, ref, allele2))]
                        if r:
                            allele2_data = [r[0][str(n)] for n in alleles_column_names if n not in ignore_exac_columns]
                        else:
                            allele2_data = ['' for n in alleles_column_names if n not in ignore_exac_columns]

                    # get genotype data from sqlite db for allele1/2
                    r = [d for d in db.execute(
                        '''SELECT * FROM genotypes WHERE CHROM = ? and POS = ? and REF = ? AND allele1 = ? and allele2 = ?;''',
                        (chom, pos, ref, allele1, allele2))]
                    if r:
                        genotype_data = [r[0][str(n)] for n in genotypes_column_names if n not in ignore_exac_columns]
                    else:
                        genotype_data = ['' for n in genotypes_column_names if n not in ignore_exac_columns]

                    # compose output line and write
                    exac_data = (per_locus_data
                        + [j for i in zip(allele1_data, allele2_data) for j in i]
                        + [d for d in genotype_data])
                    output_data = map(str, input_data + exac_data)
                    ofh.write('\t'.join(output_data) + '\n')
                    counter['output_data_records'] += 1

    # finish
    with open(output_done_file, 'w'):
        pass
    print('{} Created {} with ExAC annotations'.format(localtime_stamp(), output))
    print('''
    input lines: {input_lines}
    input data records: {input_data_records}
    output data records: {output_data_records}'''.format(**counter))

def main():

    parser = argparse.ArgumentParser(
        description = 'Add ExAC allele and genotype data to a VCF or similar file',
        epilog = 'pypeline.exac2db version 1.0β1 ©2011-2014 Michael Yourshaw all rights reserved')
    parser.add_argument('--input', '-i', required=True,
        help='input file file (tab-separated, with columns for CHROM POS REF ALT')
    parser.add_argument('--database_file', '-d', required=True,
        help='SQLite database file')
    parser.add_argument('--delim', default='\t',
        help='input file field delimiter')
    parser.add_argument('--chrom_col', default='CHROM',
        help='header line CHROM column name (default: CHROM')
    parser.add_argument('--pos_col', default='POS',
        help='header line POS column name POS (default: POS')
    parser.add_argument('--ref_col', default='REF',
        help='header line REF column name (default: REF')
    parser.add_argument('--alt_col', default='ALT',
        help='header line ALT column name (default: ALT')
    parser.add_argument('--allele1_col',
        help='header line ALLELE1 column name (default: allele1 or allele_1, or ref_col')
    parser.add_argument('--allele2_col',
        help='header line ALLELE2 column name (default: allele2 or allele_2, or alt_col')
    parser.add_argument('--skip_string', type=str,
        help='skip lines that start with this string')
    parser.add_argument('--ignore_lines', type=int, default=0,
        help='ignore the first n lines')
    parser.add_argument('--header_line', default= '#',
        help='1-based row number of column names or # if the first row that starts with a single # character')
    args = parser.parse_args()

    run(input=args.input, database_file=args.database_file, delim=args.delim, chrom_col=args.chrom_col,
        pos_col=args.pos_col, ref_col=args.ref_col, alt_col=args.alt_col, allele1_col=args.allele1_col,
        allele2_col=args.allele2_col, skip_string=args.skip_string, ignore_lines=args.ignore_lines,
        header_line=args.header_line,)


if __name__ == "__main__": sys.exit(main())
