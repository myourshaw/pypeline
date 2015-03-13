#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import glob
import argparse
import re
import time
import warnings
import my

class Vcfs2DbError(Exception): pass

#--sites_only -i /scratch1/tmp/myourshaw/resources/00-All.vcf -o /scratch1/tmp/myourshaw/resources/00-All.txt
#-i /scratch0/tmp/myourshaw/niehs/niehs.88samples.polymorphic.filtered.vcf -o /scratch0/tmp/myourshaw/niehs/db/niehs.88samples
#--sites_only -i /scratch0/tmp/myourshaw/completegenomics/20110715_complete_genomics_public_variants/ALL.wgs.complete_genomic.20110715.snps.genotypes.vcf -o /scratch0/tmp/myourshaw/completegenomics/20110715_complete_genomics_public_variants/db/completegenomics_snps

#hane
#-o /scratch0/tmp/myourshaw/hlee_vcf -i /scratch0/tmp/hlee/Cohn/SMDS/output/genotype/smds.08232011.snps.recalibrated.vcf /scratch0/tmp/hlee/Cohn/SMDS/output/genotype/smds.08232011.indels.filtered.vcf /scratch0/tmp/hlee/Cohn/Spondylothoracic/combined/output/genotype/vqsr/Spondylothoracic.snps.recalibrated.vcf /scratch0/tmp/hlee/Cohn/Spondylothoracic/combined/output/genotype/vqsr/Spondylothoracic.indels.filtered.vcf /scratch0/tmp/hlee/Cohn/Opsismodysplasia-SED-PAP/combined/output/genotype/snps/Opsismodysplasia-SED-PAP.snps.recalibrated.vcf /scratch0/tmp/hlee/Cohn/Opsismodysplasia-SED-PAP/combined/output/genotype/indels/Opsismodysplasia-SED-PAP.indels.filtered.vcf /scratch0/tmp/hlee/Cohn/Acrodysostosis/variant/acro.08312011.snps.recalibrated.vcf /scratch0/tmp/hlee/Cohn/Acrodysostosis/variant/acro.08252011.indels.filtered.vcf /scratch0/tmp/hlee/variant/08082011/output/snps/08082011.snps.recalibrated.vcf /scratch0/tmp/hlee/variant/08082011/output/indels/08082011.indels.filtered.vcf /scratch0/tmp/hlee/variant/08122011/output/vqsr.filtered/08122011.snps.recalibrated.vcf /scratch0/tmp/hlee/variant/08122011/output/vqsr.filtered/08122011.indels.filtered.vcf

#NIEHS95
#-o /scratch1/tmp/myourshaw/resources/niehs/niehs_exome_data_09212011/niehs95_db -i /scratch1/tmp/myourshaw/resources/niehs/niehs_exome_data_09212011/niehs.95samples.polymorphic.filtered.indels.vcf /scratch1/tmp/myourshaw/resources/niehs/niehs_exome_data_09212011/niehs.95samples.polymorphic.filtered.snps.vcf

#GMD32
#-o /scratch0/tmp/myourshaw/gmd/vcfs/gmd32_db -i /scratch0/tmp/myourshaw/gmd/vcfs/gmd32.analysis_ready.vcf

#GMD32
#-o /scratch0/tmp/myourshaw/gmd/vcfs/gmd32_db -i /scratch0/tmp/myourshaw/gmd/vcfs/gmd32.analysis_ready.vcf

#1kg
#-i /scratch1/tmp/myourshaw/resources/ALL.wgs.merged_beagle_mach.20101123.snps_indels_svs.sites.vcf -o /scratch1/tmp/myourshaw/resources/ALL.wgs.merged_beagle_mach.20101123.snps_indels_svs.sites_db

#GMD32
#-o /scratch0/tmp/myourshaw/gmd/vcfs/gmd32_db -i /scratch0/tmp/myourshaw/gmd/vcfs/gmd32.analysis_ready.vcf

#gmd28mms
#qout=/scratch0/tmp/myourshaw/gmd28mms/vcfs/qout
#mkdir -p $qout;
#cmd="python /home/myourshaw/lab/pypeline/python/vcfs2db.py -i /scratch0/tmp/myourshaw/gmd28mms/vcfs/gmd28mms.analysis_ready.vcf -o /scratch0/tmp/myourshaw/gmd28mms/vcfs/gmd28mms_db";
#echo "$cmd" | qsub -q all.q@compute-4* -cwd -V -M myourshaw@ucla.edu -m eas -terse -hard -pe serial 8 -l mem_free=6G -N vcfs2db_gmd28mms -e $qout -o $qout;

#dbsnp135
#-o /scratch1/tmp/myourshaw/resources/dbsnp135/dbsnp135_db -i /scratch1/tmp/myourshaw/resources/dbsnp135/00-All.vcf
#qout=/scratch1/tmp/myourshaw/resources/dbsnp135/qout
#mkdir -p $qout;
#cmd="python /home/myourshaw/lab/pypeline/python/vcfs2db.py -i /scratch1/tmp/myourshaw/resources/dbsnp135/00-All.vcf /scratch1/tmp/myourshaw/resources/dbsnp135/dbsnp135_db";
#echo "$cmd" | qsub -q all.q@compute-4* -cwd -V -M myourshaw@ucla.edu -m eas -terse -hard -pe serial 8 -l mem_free=6G -N vcfs2db_dbsnp135 -e $qout -o $qout;

def main():

    #command line arguments
    parser = argparse.ArgumentParser(
      description = 'merge multiple vcfs into a set of normalized files for database entry',
      epilog = 'pypeline.vcfs2db version 1.0β1 ©2011-2012 Michael Yourshaw all rights reserved')
    parser.add_argument('--input', '-i', nargs='+',
        help='input vcf files')
    parser.add_argument('--out_dir', '-o', required=True,
        help='output directory')
    parser.add_argument('--prefix',
        help='optional output file name prefix (default: output directory name)')
    parser.add_argument('--sites_only', action='store_true', default=False,
        help=' emit site only, ignore genotypes')
    args = parser.parse_args()
    
    if not args.prefix:
        args.prefix = os.path.basename(args.out_dir)
    
    CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,GT = range(10)
    
    #regular expressions for vcf header lines and genotypes
    fileformat_re = re.compile(r'##fileformat=[^0-9.]*(?P<version>[0-9.]+)', re.I)
    filter_re = re.compile(r'##FILTER=<ID=(?P<id>[^,]+),Description="(?P<description>[^"]+)">', re.I)
    format_re = re.compile(r'##FORMAT=<ID=(?P<id>[^,]+),Number=(?P<number>[^,]+),Type=(?P<type>[^,]+),Description="(?P<description>[^"]+)">', re.I)
    info_re = re.compile(r'##INFO=<ID=(?P<id>[^,]+),Number=(?P<number>[^,]+),Type=(?P<type>[^,]+),Description="(?P<description>[^"]+)">', re.I)
    metadata_re = re.compile(r'##(?P<key>[^=]+)="?(?P<value>.+)"?$', re.I)
    gt_re = re.compile(r'(?P<allele1>[0-9.]+)(?P<phase>[/|])(?P<allele2>[0-9.]+)', re.I)
    
    my.makedir(args.out_dir)
        
    sites_out = open(os.path.join(args.out_dir, args.prefix+'.variant.txt'), 'w')
    sites_out.write('#VCF\tTIME\tCHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')

    acafan_out = open(os.path.join(args.out_dir, args.prefix+'.acafan.txt'), 'w')
    acafan_out.write('#VCF\tTIME\tCHROM\tPOS\tREF\tALT\tAC\tAF\tAN\tAFG\n')
    
    variant_filter_out = open(os.path.join(args.out_dir, args.prefix+'.variant_filter.txt'), 'w')
    variant_filter_out.write('#VCF\tTIME\tCHROM\tPOS\tID\tREF\tALT\tFILTER\n')

    variant_info_out = open(os.path.join(args.out_dir, args.prefix+'.variant_info.txt'), 'w')
    variant_info_out.write('#VCF\tTIME\tCHROM\tPOS\tID\tREF\tALT\tKEY\tVALUE\n')
 
    id_out = open(os.path.join(args.out_dir, args.prefix+'.id.txt'), 'w')
    id_out.write('#VCF\tTIME\tCHROM\tPOS\tID\tREF\tALT\n')
    
    #per-vcf output
    filter_out = open(os.path.join(args.out_dir, args.prefix+'.filter.txt'), 'w')
    filter_out.write('#VCF\tTIME\tID\tDESCRIPTION\n')
    info_out = open(os.path.join(args.out_dir, args.prefix+'.info.txt'), 'w')
    info_out.write('#VCF\tTIME\tID\tNUMBER\tTYPE\tDESCRIPTION\n')
    metadata_out = open(os.path.join(args.out_dir, args.prefix+'.metadata.txt'), 'w')
    metadata_out.write('#VCF\tTIME\tKEY\tVALUE\n')
    
    if not args.sites_only:
        ##vcf-like file with sample id column
        #vcf_out = open(args.output, 'w')
        #vcf_out.write('VCF\tTIME\t\tSAMPLE\tCHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tGT\n')
        variant_sample_out = open(os.path.join(args.out_dir, args.prefix+'.variant_sample.txt'), 'w')
        variant_sample_out.write('#VCF\tTIME\tSAMPLE\tCHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tGT\tALLELE1\tPHASE\tALLELE2\tZYGOSITY\n')
        variant_sample_format_out = open(os.path.join(args.out_dir, args.prefix+'.variant_sample_format.txt'), 'w')
        variant_sample_format_out.write('#VCF\tTIME\tSAMPLE\tCHROM\tPOS\tID\tREF\tALT\tALLELE1\tPHASE\tALLELE2\tKEY\tVALUE\n')
        format_out = open(os.path.join(args.out_dir, args.prefix+'.format.txt'), 'w')
        format_out.write('#VCF\tTIME\tID\tNUMBER\tTYPE\tDESCRIPTION\n')
    
    vcfs = sorted(my.flatten([glob.glob(v) for v in args.input]))
    all_lines_count = 0
    for vcf in vcfs:
        vcf_time = time.strftime("%Y-%m-%d %H:%M:%S",time.localtime(os.path.getmtime(vcf)))
        format_ok = False
        samples = []
        line_count = 0
        for line in open(vcf, 'r'):
            all_lines_count += 1
            line_count += 1
            line = line.rstrip('\n')
            if not bool(line.strip()):
                continue
            if not format_ok:
                m = fileformat_re.match(line)
                if not m:
                    raise Vcfs2DbError('missing fileformat in {} at line {}'.format(vcf, line_count))
                elif not float(m.group('version')) >= 4.0:
                    raise Vcfs2DbError('obsolete fileformat in {} at line {}, {}'.format(vcf, line_count, m.group('version')))
                else:
                    format_ok = True
                    continue
            elif line.startswith('#'):
                if samples:
                    continue
                filter_m, format_m, info_m, metadata_m = (filter_re.match(line), format_re.match(line), info_re.match(line), metadata_re.match(line))
                if filter_m:
                    filter_out.write('{}\t{}\t{}\t{}\n'.format(vcf, vcf_time, filter_m.group('id'), filter_m.group('description')))
                elif info_m:
                    info_out.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(vcf, vcf_time, info_m.group('id'), info_m.group('number'), info_m.group('type'), info_m.group('description')))
                elif not args.sites_only and format_m:
                    format_out.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(vcf, vcf_time, format_m.group('id'), format_m.group('number'), format_m.group('type'), format_m.group('description')))
                elif line.upper().startswith('#CHROM'):
                    s = line.split('\t')
                    samples = s[GT:] if len(s) > GT else None
                    if not args.sites_only and not samples:
                        warnings.warn('no samples; maybe this is a sites_only file, or possibly there is a missing column header for samples in {} at line {}'.format(vcf, line_count))
                elif metadata_m and not (filter_m or info_m or format_m):
                    metadata_out.write('{}\t{}\t{}\t{}\n'.format(vcf, vcf_time, metadata_m.group('key'), metadata_m.group('value')))
                continue
            #data records
            else:
                fields = line.split('\t')
                chrom, pos, id, ref, alt, qual, filter, info = fields[:FORMAT]
                format = fields[FORMAT] if len(fields) > FORMAT else None
                gts = fields[GT:] if len(fields) > GT else None
                ids_list = id.split(';')
                filters_list = filter.split(';')
                alt_list = alt.split(',')
                these_alleles = [ref] + alt_list
                info_list = info.split(';')
                #allele counts and allele frequencies
                ac, af, an = [], [], ''
                for i in info_list:
                    #allele count in genotypes, for each ALT allele, in the same order as listed
                    if i.startswith('AC='):
                        ac = i.lstrip('AC=').split(',')
                    #allele frequency for each ALT allele in the same order as listed: use this when estimated from primary data, not called genotypes
                    elif i.startswith('AF='):
                        af = i.lstrip('AF=').split(',')
                    #total number of alleles in called genotypes
                    elif i.startswith('AN='):
                        an = i.lstrip('AN=')
                if len(ac) > 0 or len(af) > 0 or an:
                    for i in range(len(alt_list)):
                        this_ac = ac[i] if len(ac)>0 and len(ac)>=i-1 else ''
                        this_af = af[i] if len(af)>0 and len(af)>=i-1 else ''
                        afg = str(float(this_ac)/float(an)) if my.is_number(this_ac) and my.is_number(an) else ''
                        acafan_out.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(vcf, vcf_time, chrom, pos, ref, alt_list[i], this_ac, this_af, an, afg))
                for this_alt in alt_list:
                    sites_out.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(vcf, vcf_time, chrom, pos, id, ref, this_alt, qual, filter, info))
                    for this_id in ids_list:
                        if this_id != '.':
                            id_out.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(vcf, vcf_time, chrom, pos, this_id, ref, this_alt))
                    for this_filter in filters_list:
                        if this_filter != '.':
                            variant_filter_out.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(vcf, vcf_time, chrom, pos, id, ref, this_alt, this_filter))
                    for this_info in info_list:
                        if this_info != '.':
                            info_kv = this_info.split('=')
                            info_key = info_kv[0]
                            info_value = info_kv[1] if len(info_kv)>1 else ''
                            variant_info_out.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(vcf, vcf_time, chrom, pos, id, ref, this_alt, info_key, info_value))
                if not args.sites_only and samples:
                    for i in range(len(samples)):
                        this_sample = samples[i]
                        this_gt = gts[i]
                        these_format_keys, these_gt_values = (format.split(':'), this_gt.split(':'))
                        this_GT = these_gt_values[0]
                        m = gt_re.match(this_GT)
                        allele1x, phase, allele2x = m.groups() if m else ('.','/','.')
                        allele1, allele2 = (these_alleles[int(allele1x)] if allele1x != '.' else '.', these_alleles[int(allele2x)] if allele2x != '.' else '.')
                        zygosity = '' if allele1 == '.' or allele2 == '.' else 0 if allele1 == ref and allele2 == ref else 2 if allele1 == alt and allele2 == alt else 1 if (allele1 == ref and allele2 == alt) or (allele1 == alt and allele2 == ref) else int(allele1x)+int(allele2x)
                        variant_sample_out.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(vcf, vcf_time, this_sample, chrom, pos, id, ref, alt, qual, filter, info, format, this_gt, allele1, phase, allele2, zygosity))
                        if len(these_format_keys) > 1:
                            gt_keys, gt_values = (these_format_keys[1:], these_gt_values[1:])
                            for j in range(len(gt_keys)):
                                if len(gt_values) > j:
                                    this_gt_key, this_gt_value = (gt_keys[j], gt_values[j])
                                    variant_sample_format_out.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(vcf, vcf_time, this_sample, chrom, pos, id, ref, alt, allele1, phase, allele2, this_gt_key, this_gt_value))

    sites_out.close()
    acafan_out.close()
    variant_info_out.close()
    variant_filter_out.close()
    filter_out.close()
    info_out.close()
    metadata_out.close()
    if not args.sites_only:
        #vcf_out.close()
        variant_sample_out.close()
        variant_sample_format_out.close()
        format_out.close()
    print '{} lines done'.format(all_lines_count)

	
if __name__ == "__main__": sys.exit(main())
