#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
from glob import glob
import re

#-i /Volumes/storage/wsEVS/text/wsEVS_SNP_download_*.txt -o /Volumes/storage/wsEVS/text/wsEVS_SNP_download_ALL_db.txt

def flatten(l, ltypes=(list, tuple)):
    ltype = type(l)
    l = list(l)
    i = 0
    while i < len(l):
        while isinstance(l[i], ltypes):
            if not l[i]:
                l.pop(i)
                i -= 1
                break
            else:
                l[i:i + 1] = l[i]
        i += 1
    return ltype(l)

def main():
    import argparse
    parser = argparse.ArgumentParser(
        description = 'Merge text files downloaded from http://snp.gs.washington.edu/EVS/ into database-ready tab-delimited file',
        epilog = 'pypeline.wsevs_text2db version 1.0β1 ©2011-2012 Michael Yourshaw all rights reserved')
    parser.add_argument('--input', '-i', required=True, nargs='+',
        help='list of wsEVS_SNP_download_*.txt files to merge')
    parser.add_argument('--output', '-o', required=True,
        help='output file')
    args = parser.parse_args()

    out = {}
    HEADER = 'CHROM\tPOS\trsID\tALT\tREF\tEuropeanAmericanAlleleCount1\tEuropeanAmericanAlleleCount2\tAfricanAmericanAlleleCount1\tAfricanAmericanAlleleCount2\tAllAlleleCount1\tAllAlleleCount2\tMAF_EA\tMAF_AA\tMAF_All\tAvgSampleReadDepth\tGenes\tGeneAccession\tFunctionGVS\tAminoAcidChange\tProteinPos\tcDNAPos\tConservationScorePhastCons\tConservationScoreGERP\tGranthamScore\tPolyphen\tRefBaseNCBI37\tChimpAllele\tClinicalInfo\tFilterStatus\n'.format(**out)

    db = open(args.output,'w')
    db.write(HEADER)
    
    INPUT_COLS = '#base(NCBI.37) rsID Alleles EuropeanAmericanAlleleCount AfricanAmericanAlleleCount AllAlleleCount MAFinPercent(EA/AA/All) AvgSampleReadDepth Genes GeneAccession FunctionGVS AminoAcidChange ProteinPos cDNAPos ConservationScorePhastCons ConservationScoreGERP GranthamScore Polyphen RefBaseNCBI37 ChimpAllele ClinicalInfo FilterStatus'.split()
    base,rsID,Alleles,EuropeanAmericanAlleleCount,AfricanAmericanAlleleCount,AllAlleleCount,MAFinPercent,AvgSampleReadDepth,Genes,GeneAccession,FunctionGVS,AminoAcidChange,ProteinPos,cDNAPos,ConservationScorePhastCons,ConservationScoreGERP,GranthamScore,Polyphen,RefBaseNCBI37,ChimpAllele,ClinicalInfo,FilterStatus = range(len(INPUT_COLS))

    files = sorted(set(flatten([glob(f) for f in args.input])))
    for file in files:
        with open(file,'r') as f:
            line_count = 0
            for line in f:
                line_count += 1
                if line.startswith('#'):
                    continue
                fields = line.rstrip('\n').split()
                out['chrom'], out['pos'] = fields[base].split(':')
                out['rsID'] = fields[rsID]
                alleles = fields[Alleles].split('/')
                out['Allele1'] = alleles[0]
                out['Allele2'] = alleles[1]
                europeanAmericanAlleleCount = fields[EuropeanAmericanAlleleCount].split('/')
                out['EuropeanAmericanAlleleCount1'] = int(europeanAmericanAlleleCount[0])
                out['EuropeanAmericanAlleleCount2'] = int(europeanAmericanAlleleCount[1])
                africanAmericanAlleleCount = fields[AfricanAmericanAlleleCount].split('/')
                out['AfricanAmericanAlleleCount1'] = int(africanAmericanAlleleCount[0])
                out['AfricanAmericanAlleleCount2'] = int(africanAmericanAlleleCount[1])
                allAlleleCount = fields[AllAlleleCount].split('/')
                out['AllAlleleCount1'] = int(allAlleleCount[0])
                out['AllAlleleCount2'] = int(allAlleleCount[1])
                mAFinPercent = fields[MAFinPercent].split('/')
                out['MAF_EA'] = float(mAFinPercent[0])/100
                out['MAF_AA'] = float(mAFinPercent[1])/100
                out['MAF_All'] = float(mAFinPercent[2])/100
                out['AvgSampleReadDepth'] = int(fields[AvgSampleReadDepth])
                out['Genes'] = fields[Genes]
                out['GeneAccession'] = fields[GeneAccession]
                out['FunctionGVS'] = fields[FunctionGVS]
                out['AminoAcidChange'] = fields[AminoAcidChange]
                out['ProteinPos'] = fields[ProteinPos]
                out['cDNAPos'] = fields[cDNAPos]
                out['ConservationScorePhastCons'] = fields[ConservationScorePhastCons]
                out['ConservationScoreGERP'] = fields[ConservationScoreGERP]
                out['GranthamScore'] = fields[GranthamScore]
                out['Polyphen'] = fields[Polyphen]
                out['RefBaseNCBI37'] = fields[RefBaseNCBI37]
                out['ChimpAllele'] = fields[ChimpAllele]
                out['ClinicalInfo'] = fields[ClinicalInfo]
                out['FilterStatus'] = fields[FilterStatus]
                line_out = '{chrom}\t{pos}\t{rsID}\t{Allele1}\t{Allele2}\t{EuropeanAmericanAlleleCount1}\t{EuropeanAmericanAlleleCount2}\t{AfricanAmericanAlleleCount1}\t{AfricanAmericanAlleleCount2}\t{AllAlleleCount1}\t{AllAlleleCount2}\t{MAF_EA}\t{MAF_AA}\t{MAF_All}\t{AvgSampleReadDepth}\t{Genes}\t{GeneAccession}\t{FunctionGVS}\t{AminoAcidChange}\t{ProteinPos}\t{cDNAPos}\t{ConservationScorePhastCons}\t{ConservationScoreGERP}\t{GranthamScore}\t{Polyphen}\t{RefBaseNCBI37}\t{ChimpAllele}\t{ClinicalInfo}\t{FilterStatus}\n'.format(**out)
                db.write(line_out)

if __name__ == "__main__": sys.exit(main())
