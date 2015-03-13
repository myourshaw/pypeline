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

class Vcfs2AlleleFrequenciesError(Exception): pass

#dbsnp135
#-i /data/storage-1-03/archive/myourshaw/dbsnp135/ByPopulation/*.vcf.gz -o /data/storage-1-03/archive/myourshaw/dbsnp135/ByPopulation/ByPopulation_allele_frequencies.txt

def main():

    #command line arguments
    parser = argparse.ArgumentParser(
        description = 'create a database/excel-friendly vcf file by expanding ID, ALT, INFO, and GT',
        epilog = 'pypeline.vcfs2flat version 1.0β1 ©2011-2012 Michael Yourshaw all rights reserved')
    parser.add_argument('--input', '-i', nargs='+',
        help='input vcf[.gz] files')
    parser.add_argument('--output', '-o', required=True,
        help='merged ouput file')
    args = parser.parse_args()
    
    CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,GT = range(10)
    
    
    #regular expressions for vcf header lines and genotypes
    fileformat_re = re.compile(r'##fileformat=[^0-9.]*(?P<version>[0-9.]+)', re.I)
    filter_re = re.compile(r'##FILTER=<ID=(?P<id>[^,]+),Description="(?P<description>[^"]+)">', re.I)
    format_re = re.compile(r'##FORMAT=<ID=(?P<id>[^,]+),Number=(?P<number>[^,]+),Type=(?P<type>[^,]+),Description="(?P<description>[^"]+)">', re.I)
    info_re = re.compile(r'##INFO=<ID=(?P<id>[^,]+),Number=(?P<number>[^,]+),Type=(?P<type>[^,]+),Description="(?P<description>[^"]+)">', re.I)
    metadata_re = re.compile(r'##(?P<key>[^=]+)="?(?P<value>.+)"?$', re.I)
    gt_re = re.compile(r'(?P<allele1>[0-9.]+)(?P<phase>[/|])(?P<allele2>[0-9.]+)', re.I)
    
    vcfs = my.unglob(args.input)
    
    #minimal format validation and get all info columns
    has_genotypes = False
    format_ok = False
    info_cols = set()
    for vcf in vcfs:
        with my.open_gz_or_text(vcf) as vcf_in:
            for line in vcf_in:
                line = line.rstrip('\n')
                if not bool(line.strip()):
                    continue
                if not format_ok:
                    m = fileformat_re.match(line)
                    if not m:
                        raise Vcfs2AlleleFrequenciesError('missing fileformat in {}'.format(vcf))
                    elif not float(m.group('version')) >= 4.0:
                        raise Vcfs2AlleleFrequenciesError('obsolete fileformat [{}] in {}'.format(line, vcf))
                    else:
                        format_ok = True
                        continue
                elif line.startswith('#'):
                    if line.upper().startswith('#CHROM'):
                        if not line.upper().startswith('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO'):
                            raise Vcfs2AlleleFrequenciesError('invalid header format [{}] in {}'.format(line, vcf))
                        elif line.upper().startswith('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT'):
                            has_genotypes = True
                        break
                    else:
                        info_m = info_re.match(line)
                        if info_m:
                            info_cols.add(info_m.group('id'))
    info_cols = sorted(info_cols)

    #read all vcfs and write flat file output
    counters = dict(
        vcf_files = 0,
        vcf_records_in =  0,
        total_vcf_records_in = 0,
        records_out = 0,
    )
    output_cols = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','VCF','TIME'] + info_cols + ['sample_count', 'allele_count', 'ref_allele_count', 'alt_allele_count', 'het_sample_count', 'hom_sample_count', 'alt_allele_frequency']
    
    with open(args.output,'w') as variants_out:
        variants_out.write('\t'.join(output_cols)+'\n')
        for vcf in vcfs:
            counters['vcf_files'] += 1
            counters['vcf_records_in'] = 0
            with my.open_gz_or_text(vcf) as vcf_in:
                print 'reading {}'.format(vcf)
                vcf_time = time.strftime("%Y-%m-%d %H:%M:%S",time.localtime(os.path.getmtime(vcf)))
                for line in vcf_in:
                    line = line.rstrip('\n')
                    if not bool(line.strip()):
                        continue
                    #meta-information records
                    if line.startswith('#'):
                        if line.upper().startswith('#CHROM'):
                            vcf_cols = line.split('\t')
                            gt_cols = vcf_cols[GT:] if len(vcf_cols) > GT else None
                    #data records
                    else:
                        counters['vcf_records_in'] += 1
                        counters['total_vcf_records_in'] += 1
                        fields = [f.strip() for f in line.split('\t')]
                        chrom, pos, id, ref, alt, qual, filter, info = fields[:FORMAT]
                        format = fields[FORMAT] if len(fields) > FORMAT else '.'
                        gts = fields[GT:] if len(fields) > GT else None
                        ids_list = id.rstrip(';').split(';')
                        filters_list = filter.rstrip(';').split(';')
                        alt_list = alt.rstrip(',').split(',')
                        these_alleles = [ref] + alt_list
                        info_list = info.rstrip(';').split(';')
                        info_dict = {i:'' for i in info_cols}
                        format_list = format.split(':') if format else []
                        #assign values to info_dict
                        for this_info in info_list:
                            if this_info != '.':
                                info_kv = this_info.split('=')
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
                        #flatten output by IDs, ALTs, GTs
                        for this_id in ids_list:
                            for a in range(len(alt_list)):
                                this_info_dict = copy.deepcopy(info_dict)
                                this_alt = alt_list[a]
                                ac_list = info_dict.get('AC')
                                if ac_list and len(ac_list)>=a-1:
                                    this_info_dict['AC'] = ac_list[a]
                                af_list = info_dict.get('AF')
                                if af_list and len(af_list)>=a-1:
                                    this_info_dict['AF'] = af_list[a]
                                #flatten by sample
                                base_values_out = [chrom, pos, this_id, ref, this_alt, qual, filter, info]
                                info_values_out = [this_info_dict[k] for k in sorted(this_info_dict.keys())]
                                #genotypes present
                                if gt_cols:
                                    sample_count = 0
                                    allele_count = 0
                                    ref_allele_count = 0
                                    alt_allele_count = 0
                                    di_alt_allele_count = 0
                                    het_sample_count = 0
                                    hom_ref_sample_count = 0
                                    hom_alt_sample_count = 0
                                    di_het_sample_count = 0
                                    for s in range(len(gt_cols)):
                                        this_sample = gt_cols[s]
                                        this_gt = gts[s]
                                        m = gt_re.match(this_gt)
                                        allele1x, phase, allele2x = m.groups() if m else ('.','/','.')
                                        allele1, allele2 = (these_alleles[int(allele1x)] if allele1x != '.' else '.', these_alleles[int(allele2x)] if allele2x != '.' else '.')
                                        if allele1 != '.' and allele2 != '.':
                                            sample_count+=1
                                            allele_count+=2
                                            if allele1 == ref and allele2 == ref:
                                                hom_ref_sample_count+=1
                                                ref_allele_count+=2
                                            elif allele1 == this_alt and allele2 == this_alt:
                                                hom_alt_sample_count+=1
                                                alt_allele_count+=2
                                            elif (allele1 == this_alt and allele2 == ref) or (allele2 == this_alt and allele1 == ref):
                                                het_sample_count+=1
                                                alt_allele_count+=1
                                                ref_allele_count+=1
                                            elif (allele1 == this_alt and allele2 != this_alt and allele2 != ref) or (allele2 == this_alt and allele1 != this_alt and allele1 != ref):
                                                di_het_sample_count+=1
                                                alt_allele_count+=1
                                                di_alt_allele_count+=1
                                    ref_allele_frequency = float(ref_allele_count)/float(allele_count)
                                    alt_allele_frequency = float(alt_allele_count)/float(allele_count)
                                    values_out = base_values_out + [vcf, vcf_time] + info_values_out + [str(n) for n in [
                                        sample_count, allele_count, ref_allele_count, alt_allele_count, di_alt_allele_count,
                                        het_sample_count, hom_ref_sample_count, hom_alt_sample_count, di_het_sample_count,
                                        ref_allele_frequency, alt_allele_frequency]]
                                    line_out = '\t'.join(values_out)
                                    variants_out.write(line_out+'\n')
                                    counters['records_out'] += 1
                                #no genotypes
                                else:
                                    values_out = base_values_out + [vcf, vcf_time] + info_values_out + ['', '', '', '', '', '', '']
                                    line_out = '\t'.join(values_out)
                                    variants_out.write(line_out+'\n')
                                    counters['records_out'] += 1
            print 'input vcf records: {}'.format(counters['vcf_records_in'])
    print 'done. vcf files processed: {}; total input vcf records: {}; total output records: {}'.format(counters['vcf_files'],counters['total_vcf_records_in'],counters['records_out'])

	
if __name__ == "__main__": sys.exit(main())
