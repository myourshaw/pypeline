#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import glob
import re
import my

#-p 5:95751821 5:95748132 5:95735818 5:95733119 5:95730675 5:95768746 5:95751753 5:95746563 5:95746544 5:95735874 5:95735734 5:95746477 5:95746477 5:95734724 5:95757580 5:95748156 5:95751806 5:95746652 5:95751785 5:95728974 5:95728898 -i /scratch1/tmp/myourshaw/resources/1000genomes_release_v3/ALL.chr5.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz.vax.gz.flax.gz -o /home/myourshaw/resources_1000genomes_release_v3_ALL.chr5.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz.vax.gz.flax.gz.PCSK1.flax
#-p 5:95751821 5:95748132 5:95735818 5:95733119 5:95730675 5:95768746 5:95751753 5:95746563 5:95746544 5:95735874 5:95735734 5:95746477 5:95746477 5:95734724 5:95757580 5:95748156 5:95751806 5:95746652 5:95751785 5:95728974 5:95728898 -i /scratch0/tmp/myourshaw/1000genomes/phase1_integrated_release_version3_20120430/ALL.chr5.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz -o /home/myourshaw/1000genomes_phase1_integrated_release_version3_20120430_ALL.chr5.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz.PCSK1.vcf
#-p 5:95751821 5:95748132 5:95735818 5:95733119 5:95730675 5:95768746 5:95751753 5:95746563 5:95746544 5:95735874 5:95735734 5:95746477 5:95746477 5:95734724 5:95757580 5:95748156 5:95751806 5:95746652 5:95751785 5:95728974 5:95728898 -i /scratch0/tmp/myourshaw/dbsnp/dbsnp137/human_9606_dbSNP137_20120616/VCF/05-*-*.vcf.gz -o /home/myourshaw/dbsnp_dbsnp137_human_9606_dbSNP137_20120616_VCF_05-all-all.vcf.gz.PCSK1.vcf
#-p 5:95751821 5:95748132 5:95735818 5:95733119 5:95730675 5:95768746 5:95751753 5:95746563 5:95746544 5:95735874 5:95735734 5:95746477 5:95746477 5:95734724 5:95757580 5:95748156 5:95751806 5:95746652 5:95751785 5:95728974 5:95728898 -i /scratch0/tmp/myourshaw/dbsnp/dbsnp137/human_9606_dbSNP137_20120616/VCF/00-All.vcf.gz -o /home/myourshaw/dbsnp_dbsnp137_human_9606_dbSNP137_20120616_VCF_00-All.vcf.gz.PCSK1.vcf

def main():

    #command line arguments
    parser = argparse.ArgumentParser(
        description = 'grep vcf file(s) for certain positions',
        epilog = 'vcf_grep version 1.0β1 ©2011 Michael Yourshaw all rights reserved')
    parser.add_argument('--vcfs', '-i', nargs='+',
        help='input vcf file(s)')
    parser.add_argument('--output', '-o',
        help='input vcf file(s)')
    parser.add_argument('--positions', '-p', nargs='+',
        help='list of positions to grep (chrom:pos)')
    args = parser.parse_args()
    
    vcfs=my.unglob(args.vcfs)
    positions = [(p.split(':')[0],p.split(':')[1]) for p in args.positions]
    write_header = True
    with open(args.output,'w') as o:
        for vcf in vcfs:
            with my.open_gz_or_text(vcf) as v:
                o.write('##file='+vcf+'\n')
                for line in v:
                    if line.strip() == '':
                        continue
                    elif line.startswith('#'): #and write_header:
                        o.write(line)
                        write_header = False
                    else:
                        fields = line[:31].split('\t')
                        if (fields[0],fields[1]) in positions:
                            o.write(line)
                    
    
    
    
if __name__ == "__main__": sys.exit(main())
