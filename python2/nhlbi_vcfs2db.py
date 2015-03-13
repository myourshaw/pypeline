#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import glob
import argparse
import re
import time
import warnings
import gzip
import copy
import my
import xml.sax.saxutils

class NhlbiVcfs2DbError(Exception): pass

#-i /scratch0/tmp/myourshaw/nhlbi/ESP6500_20120903/ESP6500.snps.vcf/ESP6500.chrX.snps.vcf -o /scratch0/tmp/myourshaw/nhlbi/ESP6500_20120903/ESP6500.snps.vcf/test
#-i /scratch0/tmp/myourshaw/nhlbi/ESP6500SI_20121126/snps_indels.vcf/ESP6500SI.chr22.snps_indels.vcf -o /scratch0/tmp/myourshaw/nhlbi/ESP6500SI_20121126/snps_indels.vcf/ESP6500SI.chr22.snps_indels.vcf.db.test
#-i /data/storage-1-01/archive/myourshaw/scratch0/nhlbi/ESP6500SI_20121126/indel_test.vcf -o /data/storage-1-01/archive/myourshaw/scratch0/nhlbi/ESP6500SI_20121126/indel_test.vcf.db
#-i /data/storage-1-01/archive/myourshaw/scratch0/nhlbi/ESP6500SI_20121126/snps_indels.vcf/ESP6500SI.chrX.snps_indels.vcf -o /data/storage-1-01/archive/myourshaw/scratch0/nhlbi/ESP6500SI_20121126/snps_indels.vcf/ESP6500SI.chrX.snps_indels.vcf
#-i /data/storage-1-01/archive/myourshaw/scratch0/nhlbi/ESP6500SI_20121126/snps_indels.vcf/ESP6500SI.chrY.snps_indels.vcf -o /data/storage-1-01/archive/myourshaw/scratch0/nhlbi/ESP6500SI_20121126/ESP6500SI.chrY.snps_indels.vcf
#--uncompressed -i /scratch1/tmp/myourshaw/resources/nhlbi/ESP6500SI-V2_20131218/ESP6500SI-V2-SSA137.dbSNP138-rsIDs.snps_indels.vcf/ESP6500SI-V2-SSA137.updatedRsIds.chr22.snps_indels.vcf -o /scratch1/tmp/myourshaw/resources/nhlbi/ESP6500SI-V2_20131218/ESP6500SI-V2-SSA137.dbSNP138-rsIDs.snps_indels.vcf/ESP6500SI-V2-SSA137.updatedRsIds.chr22.snps_indels.vcf
#--uncompressed -i /scratch1/tmp/myourshaw/resources/nhlbi/ESP6500SI-V2_20131218/ESP6500SI-V2-SSA137.dbSNP138-rsIDs.snps_indels.vcf/ESP6500SI-V2-SSA137.updatedRsIds.chrY.snps_indels.vcf -o /scratch1/tmp/myourshaw/resources/nhlbi/ESP6500SI-V2_20131218/ESP6500SI-V2-SSA137.dbSNP138-rsIDs.snps_indels.vcf/ESP6500SI-V2-SSA137.updatedRsIds.chrY.snps_indels.vcf
#--uncompressed -i /scratch1/tmp/myourshaw/resources/nhlbi/ESP6500SI-V2_20131218/ESP6500SI-V2-SSA137.dbSNP138-rsIDs.snps_indels.vcf/ESP6500SI-V2-SSA137.updatedRsIds.chr*.snps_indels.vcf -o /scratch1/tmp/myourshaw/resources/nhlbi/ESP6500SI-V2_20131218/ESP6500SI-V2-SSA137.dbSNP138-rsIDs.snps_indels.vcf/ESP6500SI-V2-SSA137.updatedRsIds.chr_.snps_indels.vcf
def main():

    #command line arguments
    parser = argparse.ArgumentParser(
        description = 'create a database/excel-friendly vcf file by expanding alleles and genotypes',
        epilog = 'pypeline.nhlbi_vcfs2db version 1.0β1 ©2011-2012 Michael Yourshaw all rights reserved')
    parser.add_argument('--input', '-i', nargs='+',
        help='input vcf[.gz] files')
    parser.add_argument('--output', '-o', required=True,
        help='prefix for merged ouput files')
    parser.add_argument('--AApop', default='ASW',
        help='African American population ID (default: ASW)')
    parser.add_argument('--AAsuper', default='AFR',
        help='African American super-population ID (default: AFR)')
    parser.add_argument('--EApop', default='CEU',
        help='European American population ID (default: CEU)')
    parser.add_argument('--EAsuper', default='EUR',
        help='European American super-population ID (default: EUR)')
    parser.add_argument('--uncompressed', action='store_true', default=False,
        help='do not gzip output (default: False)')
    args = parser.parse_args()
    
    db_output_file = args.output + '_db.txt' + ('' if args.uncompressed else '.gz')
    allele_output_file = args.output + '_allele_counts.txt' + ('' if args.uncompressed else '.gz')
    genotype_output_file = args.output + '_genotype_counts.txt' + ('' if args.uncompressed else '.gz')

    CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,GT = range(10)
    
    #regular expressions for vcf header lines and genotypes
    fileformat_re = re.compile(r'##fileformat=[^0-9.]*(?P<version>[0-9.]+)', re.I)
    filter_re = re.compile(r'##FILTER=<ID=(?P<id>[^,]+),Description="(?P<description>[^"]+)">', re.I)
    format_re = re.compile(r'##FORMAT=<ID=(?P<id>[^,]+),Number=(?P<number>[^,]+),Type=(?P<type>[^,]+),Description="(?P<description>[^"]+)">', re.I)
    info_re = re.compile(r'##INFO=<ID=(?P<id>[^,]+),Number=(?P<number>[^,]+),Type=(?P<type>[^,]+),Description="(?P<description>[^"]+)">', re.I)
    metadata_re = re.compile(r'##(?P<key>[^=]+)="?(?P<value>.+)"?$', re.I)
    gt_re = re.compile(r'(?P<allele1>[0-9.]+)(?P<phase>[/|])(?P<allele2>[0-9.]+)', re.I)
    indel_gt_re = re.compile(r'(?P<AR1>[AR])(?P<IX1>\d*)(?P<AR2>[AR]*)(?P<IX2>\d*)',re.I)
    
    vcfs = my.unglob(args.input, sort_unique=False)

    #minimal format validation and get all info and gt columns
    has_genotypes = False
    format_ok = False
    info_cols = set()
    gt_cols = set()
    for vcf in vcfs:
        with my.open_gz_or_text(vcf) as vcf_in:
            for line in vcf_in:
                line = line.rstrip('\n')
                if not bool(line.strip()):
                    continue
                if not format_ok:
                    m = fileformat_re.match(line)
                    if not m:
                        warnings.warn('missing fileformat in {}'.format(vcf))
                    elif not float(m.group('version')) >= 4.0:
                        warnings.warn('obsolete fileformat [{}] in {}'.format(line, vcf))
                    else:
                        format_ok = True
                        continue
                elif line.startswith('#'):
                    if line.upper().startswith('#CHROM'):
                        if not line.upper().startswith('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO'):
                            raise NhlbiVcfs2DbError('invalid header format [{}] in {}'.format(line, vcf))
                        elif line.upper().startswith('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'):
                            has_genotypes = True
                            warnings.warn('input file {} has genotypes'.format(vcf))
                            gt_cols = gt_cols | set(line.split('\t')[GT:])
                        break
                    else:
                        info_m = info_re.match(line)
                        if info_m:
                            info_cols.add(info_m.group('id'))
    info_cols = [i+'_info' if i in ['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT'] else i for i in sorted(info_cols)]
    gt_cols = sorted(gt_cols)

    #read all vcfs and write flat file output
    counters = dict(
        vcf_files = 0,
        vcf_records_in =  0,
        total_vcf_records_in = 0,
        genotypes_out = 0,
        alleles_out = 0,
    )
    
    db_output_cols = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','allele1','allele2','oALT','EA_AC_allele1','EA_AC_allele2','AA_AC_allele1','AA_AC_allele2','TAC_allele1','TAC_allele2','MAF_EA','MAF_AA','MAF_All','EA_GTC_this','AA_GTC_this','GTC_this'] + info_cols
    genotype_output_cols = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','allele1','allele2','oALT','population','super_population','sample_count']
    allele_output_cols = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','allele','oALT','population','super_population','allele_count']
    with  open(db_output_file,'w') if args.uncompressed else gzip.open(db_output_file, 'w') as db_out, open(genotype_output_file,'w') if args.uncompressed else gzip.open(genotype_output_file, 'w') as genotypes_out, open(allele_output_file,'w') if args.uncompressed else gzip.open(allele_output_file, 'w') as alleles_out:
        db_out.write('\t'.join(db_output_cols)+'\n')
        genotypes_out.write('\t'.join(genotype_output_cols)+'\n')
        alleles_out.write('\t'.join(allele_output_cols)+'\n')
        for vcf in vcfs:
            counters['vcf_files'] += 1
            counters['vcf_records_in'] = 0
            with my.open_gz_or_text(vcf) as vcf_in:
                print 'reading {}'.format(vcf)
                vcf_time = time.strftime("%Y-%m-%d %H:%M:%S",time.localtime(os.path.getmtime(vcf)))
                vcf_line_count = 0
                for line in vcf_in:
                    vcf_line_count += 1
                    line = line.rstrip('\n')
                    if not bool(line.strip()):
                        continue
                    #meta-information records
                    if line.startswith('#'):
                        if line.upper().startswith('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO'):
                            vcf_cols = line.split('\t')
                            these_gt_cols = vcf_cols[GT:] if len(vcf_cols) > GT else []
                    #data records
                    else:
                        counters['vcf_records_in'] += 1
                        counters['total_vcf_records_in'] += 1
                        fields = [f.strip() for f in line.split('\t')]
                        chrom, pos, id, ref, oalt, qual, filter, info = fields[:FORMAT]
                        chrom = chrom.upper()
                        #unescape things like &amp;
                        info = xml.sax.saxutils.unescape(info)
                        format = fields[FORMAT] if len(fields) > FORMAT else '.'
                        this_gt_dict = {g:'' for g in these_gt_cols}
                        gts = fields[GT:] if len(fields) > GT else None
                        if gts:
                            for i in range(len(these_gt_cols)):
                                this_gt_dict[these_gt_cols[i]] = gts[i]
#TODO: VCF 4.1 allows additional info in ALT and ID
                        ids_list = id.rstrip(';').split(';')
                        filters_list = filter.rstrip(';').split(';')
                        alt_list = oalt.rstrip(',').split(',')
                        these_alleles = [ref] + alt_list
                        alts_ref = alt_list + [ref]
                        #oALT may have multiple ,-separated alleles
                        info = 'oALT='+oalt+';'+info if info != '.' else 'oALT='+oalt+';'
                        info_list = info.rstrip(';').split(';')
                        info_dict = {i:'' for i in info_cols}
                        format_list = format.split(':') if format else []
                        #assign values to info_dict
                        for this_info in info_list:
                            if this_info != '.':
                                info_kv = this_info.split('=',1)
                                if info_kv[0] in ['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']:
                                    info_kv[0]+='_info' 
                                info_key = info_kv[0]
                                #value or 1 for presence of flag
                                info_value = info_kv[1] if len(info_kv)>1 else '1'
                                info_dict[info_key] = info_value
                                #allele count in genotypes, for each ALT allele, in the same order as listed
                                if info_key == 'AC':
                                    info_dict[info_key] = info_value.split(',')
                                #allele frequency for each ALT allele in the same order as listed: use this when estimated from primary data, not called genotypes
                                elif info_key == 'AF':
                                    info_dict[info_key] = info_value.split(',')
                                if (info_key == 'CG' or info_key == 'CP') and info_value == 'NA':
                                    info_dict[info_key] = ''
                        info_data = [info_dict[i] for i in info_cols]
                        #flatten output by genotypes and alleles and alts
                        ac_list = info_dict.get('AC')
                        if ac_list and len(ac_list)>=a-1:
                            this_info_dict['AC'] = ac_list[a]
                        af_list = info_dict.get('AF')
                        if af_list and len(af_list)>=a-1:
                            this_info_dict['AF'] = af_list[a]

                        if len(alt_list)>1:
                            pass #DEBUG
                        #Observed Genotypes. For INDELs, A1, A2, or An refers to the N-th alternate allele while R refers to the reference allele.
                        #convert 'AA,AG,GG' to [['A','A'],['A','G'],['G','G']]
                        #chr X and Y can be hemizygous; convert 'A,G' to [['A','.'],['G','.']]
                        gts_list = info_dict.get('GTS').split(',')
                        has_indel = sum([len(s) for s in these_alleles])!=len(these_alleles)
                        if has_indel:
                            for i in range(len(gts_list)):
                                a1,a2 = '.','.'
                                indel_gt_m = indel_gt_re.match(gts_list[i])
                                if indel_gt_m:
                                    if indel_gt_m.group('AR1') == 'R': #ref
                                        a1 = these_alleles[0]
                                    if indel_gt_m.group('AR1') == 'A': #alt
                                        a1 = these_alleles[int(indel_gt_m.group('IX1'))]
                                    #chr X and Y can be hemizygous
                                    if indel_gt_m.group('AR2') == 'R': #ref
                                        a2 = these_alleles[0]
                                    if indel_gt_m.group('AR2') == 'A': #alt
                                        a2 = these_alleles[int(indel_gt_m.group('IX2'))]
                                    gts_list[i]=[a1,a2]
                                else:
                                    raise NhlbiVcfs2DbError('invalid indel GTS format [{}] in {} line {}'.format(info_dict.get('GTS'), vcf, vcf_line_count))
                        else:
                            gts_list = [list(s) if len(s)>1 else [s,'.'] for s in gts_list]
                        #African American Genotype Counts in the order of listed GTS
                        aa_gtc_list = info_dict['AA_GTC'].split(',')
                        #European American Genotype Counts in the order of listed GTS
                        ea_gtc_list = info_dict['EA_GTC'].split(',')
                        #Total Genotype Counts in the order of listed GTS
                        gtc_list = info_dict['GTC'].split(',')
                        #African American Allele Count in the order of AltAlleles,RefAllele.  For INDELs, A1, A2, or An refers to the N-th alternate allele while R refers to the reference allele.
                        aa_ac_list = info_dict.get('AA_AC').split(',')
                        #European American Allele Count in the order of AltAlleles,RefAllele.  For INDELs, A1, A2, or An refers to the N-th alternate allele while R refers to the reference allele.
                        ea_ac_list = info_dict.get('EA_AC').split(',')
                        #Total Allele Count in the order of AltAlleles,RefAllele.  For INDELs, A1, A2, or An refers to the N-th alternate allele while R refers to the reference allele.
                        tac_list = info_dict.get('TAC').split(',')
                        #Minor Allele Frequency in percent in the order of EA,AA,All
                        maf_ea,maf_aa,maf_All = info_dict['MAF'].split(',')
                        for this_alt in alt_list:
                            out_list = [chrom, pos, id, ref, this_alt, qual, filter, info]
                            for i in range(len(gts_list)):
                                allele1, allele2 = gts_list[i] if len(gts_list[i]) == 2 else [gts_list[i][0],'.']
                                if allele1 != ref and allele2 == ref:
                                    allele1, allele2 = (allele2, allele1)
                                ea_gtc = ea_gtc_list[i]
                                aa_gtc = aa_gtc_list[i]
                                gtc = gtc_list[i]
                                allele1x = alts_ref.index(allele1) if allele1 != '.' else None
                                allele2x = alts_ref.index(allele2) if allele2 != '.' else None
                                ea_ac_allele1 = ea_ac_list[allele1x] if allele1x != None else ''
                                ea_ac_allele2 = ea_ac_list[allele2x] if allele2x != None else ''
                                aa_ac_allele1 = aa_ac_list[allele1x] if allele1x != None else ''
                                aa_ac_allele2 = aa_ac_list[allele2x] if allele2x != None else ''
                                tac_allele1 = tac_list[allele1x] if allele1x != None else ''
                                tac_allele2 = tac_list[allele2x] if allele2x != None else ''
                                db_out_list = out_list + [allele1, allele2, oalt, ea_ac_allele1, ea_ac_allele2, aa_ac_allele1, aa_ac_allele2, tac_allele1, tac_allele2, maf_ea, maf_aa, maf_All, ea_gtc, aa_gtc, gtc] + info_data
                                db_out.write('\t'.join(db_out_list)+'\n')
                            for population in [(aa_gtc_list,args.AApop,args.AAsuper), (ea_gtc_list,args.EApop,args.EAsuper)]:
                                for i in range(len(gts_list)):
                                    allele1, allele2 = gts_list[i] if len(gts_list[i]) == 2 else [gts_list[i][0],'']
                                    pop,super_pop = population[1:3]
                                    sample_count = population[0][i]
                                    genotype_out_list = out_list + [allele1, allele2, oalt, pop, super_pop, sample_count]
                                    genotypes_out.write('\t'.join(genotype_out_list)+'\n')
                                    counters['genotypes_out']+=1
                            for population in [(aa_ac_list,args.AApop,args.AAsuper), (ea_ac_list,args.EApop,args.EAsuper)]:
                                for i in range(len(alts_ref)):
                                    allele = alts_ref[i]
                                    pop,super_pop = population[1:3]
                                    sample_count = population[0][i]
                                    allele_out_list = out_list + [allele, oalt, pop, super_pop, sample_count]
                                    alleles_out.write('\t'.join(allele_out_list)+'\n')
                                    counters['alleles_out']+=1
            print 'input vcf records: {}'.format(counters['vcf_records_in'])
    print 'done. vcf files processed: {}; total input vcf records: {}; total genotype records: {}; total allele records: {}'.format(counters['vcf_files'],counters['total_vcf_records_in'],counters['genotypes_out'],counters['alleles_out'])
    

if __name__ == "__main__": sys.exit(main())
