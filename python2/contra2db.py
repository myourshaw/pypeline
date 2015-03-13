#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import glob
import argparse
import re
import time
import warnings
import copy
import my
import sql_columns


#-i /scratch0/tmp/myourshaw/contra/GMD*/table/GMD*.* /scratch0/tmp/myourshaw/contra/JEN*/table/JEN*.* -o /scratch0/tmp/myourshaw/contra/merge/GMDJEN_Illumina
#-i /scratch0/tmp/myourshaw/contra/GMD154A/table/GMD*.* -o /scratch0/tmp/myourshaw/contra/merge/GMD154A_Agilent

class Contra2DbError(Exception): pass

def main():

    #command line arguments
    parser = argparse.ArgumentParser(
        description = 'merge contra output prior to database input',
        epilog = 'pypeline.contra2db version 1.0β1 ©2011-2013 Michael Yourshaw all rights reserved')
    parser.add_argument('--input', '-i', required = True, nargs='+',
        help='input files (files must begin with <sampleName>.)')
    parser.add_argument('--output', '-o', required = True,
        help='output prefix (dir/prefix)')
    args = parser.parse_args()
    
    output_dir = os.path.dirname(args.output)
    my.makedir(output_dir)
    prefix = os.path.basename(args.output)
    
    DetailsFILTERED = os.path.join(output_dir, prefix+'_DetailsFILTERED_merged.txt')
    TargetRegions = os.path.join(output_dir, prefix+'_TargetRegions_merged.txt')
    LargeDeletion = os.path.join(output_dir, prefix+'_LargeDeletion_merged.txt')
    vcf_file = os.path.join(output_dir, prefix+'_vcf_merged.txt')
    
    files = my.unglob(args.input)
    
    details_header = 'Targeted.Region.ID\tExon.Number\tGene.Sym\tChr\tOriStCoordinate\tOriEndCoordinate\tMean.of.LogRatio\tAdjusted.Mean.of.LogRatio\tSD.of.LogRatio\tMedian.of.LogRatio\tnumber.bases\tP.Value\tAdjusted.P.Value\tgain.loss\ttumour.rd\tnormal.rd\ttumour.rd.ori\tnormal.rd.ori\tMinLogRatio\tMaxLogRatio\tBinNumber'
    target_header = 'Targeted.Region.ID\tExon.Number\tGene.Sym\tChr\tOriStCoordinate\tOriEndCoordinate\tMean.of.LogRatio\tAdjusted.Mean.of.LogRatio\tSD.of.LogRatio\tMedian.of.LogRatio\tnumber.bases\tP.Value\tAdjusted.P.Value\tgain.loss\ttumour.rd\tnormal.rd\ttumour.rd.ori\tnormal.rd.ori\tMinLogRatio\tMaxLogRatio\tBinNumber'
    deletion_header = 'Chr\tTarget.Start\tTarget.End\tNumberOfTargets\tOriStCoordinate\tOriEndCoordinate\tCBS.Mean\tLogRatios\tAbove.PValues.Cutoff\tCalls'
    vcf_header = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO'
    
    with open (DetailsFILTERED, 'w') as details, open (TargetRegions, 'w') as target, open (LargeDeletion, 'w') as deletion, open (vcf_file, 'w') as vcf:
        details.write(details_header.replace('.','_')+'\tsample\n')
        target.write(target_header.replace('.','_')+'\tsample\n')
        deletion.write(deletion_header.replace('.','_')+'\tsample\n')
        vcf.write(vcf_header.replace('.','_')+'\tsample\n')
        
        for file in files:
            sample = os.path.basename(file).split('.')[0]
            output = False
            with open(file) as input:
                for line in input:
                    line = line.rstrip('\n')
                    if line.strip() == '':
                        continue
                    if line.startswith('##'): #VCF headers
                        continue
                    if not output:
                        #vcf files have trailing spaces in header
                        line = '\t'.join([f.strip() for f in line.split('\t')])
                        output = details if line == details_header and file.endswith('.DetailsFILTERED.txt') else target if line == target_header and not file.endswith('.DetailsFILTERED.txt') else deletion if line == deletion_header else vcf if line == vcf_header else None
                        if not output:
                            break
                        else:
                            print 'Merging {}'.format(file)
                    else:
                        output.write(line+'\t'+sample+'\n')
                        
    sql_columns.run(input=[DetailsFILTERED, TargetRegions, LargeDeletion, vcf_file], header_line=1)
    
if __name__ == "__main__": sys.exit(main())
