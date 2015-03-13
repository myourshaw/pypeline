#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import glob
import argparse
import re
import gzip
import my

class Vcfs2DbsnpPopulationIndividualError(Exception): pass

#-i /scratch1/tmp/myourshaw/resources/dbsnp135_20120118/ByPopulation/*.vcf.gz
#-i /scratch0/tmp/myourshaw/dbsnp/dbsnp137/human_9606_dbSNP137_downloaded_20120816/VCF/ByPopulation/*.vcf.gz

def main():

    #command line arguments
    parser = argparse.ArgumentParser(
        description = 'extract population and individual ids from dbSNP vcf files',
        epilog = 'pypeline.dbsnp_population_individual version 1.0β1 ©2011-2012 Michael Yourshaw all rights reserved')
    parser.add_argument('--input', '-i', nargs='+',
        help='input vcf[.gz] files')
    args = parser.parse_args()
    
    CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,GT = range(10)
    
    pop = {}
    vcfs = my.unglob(args.input)
    for vcf in vcfs:
        dbSNP_POP_ID = None
        dbSNP_LOC_POP_ID = None
        with my.open_gz_or_text(vcf) as vcf_in:
            for line in vcf_in:
                line = line.rstrip('\n')
                if not bool(line.strip()):
                    continue
                elif line.startswith('#'):
                    if line.startswith('##dbSNP_POP_ID='):
                        dbSNP_POP_ID = line.split('=')[1]
                    elif line.startswith('##dbSNP_LOC_POP_ID='):
                        dbSNP_LOC_POP_ID = line.split('=')[1].split('-')[1]
                    elif line.upper().startswith('#CHROM'):
                        if not (dbSNP_POP_ID and dbSNP_LOC_POP_ID):
                            raise Vcfs2DbsnpPopulationIndividualError('no dbSNP_POP_ID or dbSNP_LOC_POP_ID before [{}] in {}'.format(line, vcf))
                        if not line.upper().startswith('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO'):
                            raise Vcfs2DbsnpPopulationIndividualError('invalid header format [{}] in {}'.format(line, vcf))
                        elif line.upper().startswith('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'):
                            gt_cols = set(line.split('\t')[GT:])
                            if pop.get(dbSNP_POP_ID):
                                pop[dbSNP_POP_ID][1] |= gt_cols
                            else:
                                pop[dbSNP_POP_ID] = [dbSNP_LOC_POP_ID,gt_cols]
                        break
    print 'dbSNP_POP_ID\tdbSNP_LOC_POP_ID\tSAMPLE'
    for POP_ID in pop.keys():
        for indiv in sorted(pop[POP_ID][1]):
            print '{}\t{}\t{}'.format(POP_ID, pop[POP_ID][0], indiv)
    pass
        


if __name__ == "__main__": sys.exit(main())
