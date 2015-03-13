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

class DbsnpVcfs2DbError(Exception): pass

#--uncompressed -i /scratch1/tmp/myourshaw/resources/nhlbi/ESP6500SI-V2_20131218/ESP6500SI-V2-SSA137.dbSNP138-rsIDs.snps_indels.vcf/ESP6500SI-V2-SSA137.updatedRsIds.chr*.snps_indels.vcf -o /scratch1/tmp/myourshaw/resources/nhlbi/ESP6500SI-V2_20131218/ESP6500SI-V2-SSA137.dbSNP138-rsIDs.snps_indels.vcf/ESP6500SI-V2-SSA137.updatedRsIds.chr_.snps_indels.vcf
#--uncompressed -i /scratch1/tmp/myourshaw/resources/dbsnp/138/00-All.vcf.gz -o /scratch1/tmp/myourshaw/resources/dbsnp/138/dbsnp138_00_All.db
def run(input, output, uncompressed=False):
    
    vcfs = my.unglob(input, sort_unique=False)
    if not uncompressed and not output.endswith('.gz'):
        output = output+'.gz'
    CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT = range(9)
    
    #regular expressions for vcf header lines and genotypes
    fileformat_re = re.compile(r'##fileformat=[^0-9.]*(?P<version>[0-9.]+)', re.I)
    filter_re = re.compile(r'##FILTER=<ID=(?P<id>[^,]+),Description="(?P<description>[^"]+)">', re.I)
    format_re = re.compile(r'##FORMAT=<ID=(?P<id>[^,]+),Number=(?P<number>[^,]+),Type=(?P<type>[^,]+),Description="(?P<description>[^"]+)">', re.I)
    info_re = re.compile(r'##INFO=<ID=(?P<id>[^,]+),Number=(?P<number>[^,]+),Type=(?P<type>[^,]+),Description="(?P<description>[^"]+)">', re.I)
    metadata_re = re.compile(r'##(?P<key>[^=]+)="?(?P<value>.+)"?$', re.I)
    gt_re = re.compile(r'(?P<allele1>[0-9.]+)(?P<phase>[/|])(?P<allele2>[0-9.]+)', re.I)
    
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
                            raise DbsnpVcfs2DbError('invalid header format [{}] in {}'.format(line, vcf))
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
        db_records_out = 0,
    )
    
    db_output_cols = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','oALT'] + info_cols
    with  open(output,'w') if uncompressed else gzip.open(output, 'w') as db_out:
        db_out.write('\t'.join(db_output_cols)+'\n')
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
                        info_data = [info_dict[i] for i in info_cols]
                        #flatten output by alts

                        for this_alt in alt_list:
                            out_list = [chrom, pos, id, ref, this_alt, qual, filter, info]
                            db_out_list = out_list + [oalt] + info_data
                            db_out.write('\t'.join(db_out_list)+'\n')
                            counters['db_records_out'] += 1
            print 'input vcf records: {}'.format(counters['vcf_records_in'])
    print 'done. vcf files processed: {}; total input vcf records: {}; output db records: {}'.format(counters['vcf_files'],counters['total_vcf_records_in'],counters['db_records_out'])
    

def main():

    #command line arguments
    parser = argparse.ArgumentParser(
        description = 'create a database/excel-friendly vcf file by expanding alleles',
        epilog = 'pypeline.vax_dbsnp_vcfs2db version 1.0β1 ©2011-2014 Michael Yourshaw all rights reserved')
    parser.add_argument('--input', '-i', nargs='+',
        help='input vcf[.gz] files')
    parser.add_argument('--output', '-o', required=True,
        help='prefix for merged ouput files')
    parser.add_argument('--uncompressed', action='store_true', default=False,
        help='do not gzip output (default: False)')
    args = parser.parse_args()
    
    run(input=args.input, output=args.output, uncompressed=args.uncompressed)


if __name__ == "__main__": sys.exit(main())
