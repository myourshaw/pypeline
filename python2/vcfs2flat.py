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

class Vcfs2FlatError(Exception): pass

#dbsnp135
#-i /scratch1/tmp/myourshaw/resources/dbsnp135/00-All.vcf

#-o /scratch1/tmp/myourshaw/resources/dbsnp135/ByPopulationNoGeno/ByPopulationNoGeno.variants.txt -i /scratch1/tmp/myourshaw/resources/dbsnp135/ByPopulationNoGeno/*.vcf.gz

#-i /scratch1/tmp/myourshaw/resources/dbsnp135/ByPopulation/*.vcf.gz -o /scratch1/tmp/myourshaw/resources/dbsnp135/ByPopulation/ByPopulation.txt

#-i /scratch1/tmp/myourshaw/resources/dbsnp135/ByPopulation/MEX-12161-MT.vcf.gz /scratch1/tmp/myourshaw/resources/dbsnp135/ByPopulation/CHB-12157-MT.vcf.gz -o /scratch1/tmp/myourshaw/resources/dbsnp135/ByPopulation/test.txt

#--wsevs -i /Volumes/scratch/wsEVS/20111205/vcf/wsEVS_SNP_download_chr*.vcf -o /Volumes/scratch/wsEVS/20111205/vcf/wsEVS_SNP_download_db.txt

#1000genomes sites
#-i /Volumes/scratch/ncbi/1000genomes/ftp/release/20110521/ALL.wgs.phase1_integrated_calls.20101123.snps_indels_svs.sites.vcf.gz -o /Volumes/scratch/ncbi/1000genomes/ftp/release/20110521/ALL_wgs_phase1_integrated_calls_20101123_snps_indels_svs_sites_db.txt

#1000genomes genotypes
#-i /Volumes/scratch/ncbi/1000genomes/ftp/release/20110521/ALL.chr*.phase1_integrated_calls.20101123.snps_indels_svs.genotypes.vcf.gz -o /Volumes/scratch/ncbi/1000genomes/ftp/release/20110521/ALL_phase1_integrated_calls_20101123_snps_indels_svs_genotypes_db.txt

#1000 genomes with populations
#-i /scratch1/tmp/myourshaw/resources/1000genomes_release_v3/ALL.chr22.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz -o /scratch1/tmp/myourshaw/resources/1000genomes_release_v3/ALL.chr22.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz.flat --samples_populations /scratch1/tmp/myourshaw/resources/1000genomes_release_v3/samples_populations/20111108_1000genomes_samples_sex_population.txt
#-i /scratch1/tmp/myourshaw/resources/1000genomes_release_v3/sorted_vcf/ALL.chr*.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz -o /scratch1/tmp/myourshaw/resources/1000genomes_release_v3/ALL.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz.allele_frequencies.txt --samples_populations /scratch1/tmp/myourshaw/resources/1000genomes_release_v3/samples_populations/20111108_1000genomes_samples_sex_population.txt

#dbsnp135 with populations
#-i /scratch1/tmp/myourshaw/resources/dbsnp135_20120118/ByPopulation/*-*-*.vcf.gz -o /scratch1/tmp/myourshaw/resources/dbsnp135_20120118/dbsnp135_20120118_ByPopulation.allele_frequencies.txt --samples_populations /scratch1/tmp/myourshaw/resources/dbsnp135_20120118/samples_populations/hapmap_samples_sex_population.txt

def main():

    #command line arguments
    parser = argparse.ArgumentParser(
        description = 'create a database/excel-friendly vcf file by expanding ID, ALT, INFO, and GT',
        epilog = 'pypeline.vcfs2flat version 1.0β1 ©2011-2012 Michael Yourshaw all rights reserved')
    parser.add_argument('--input', '-i', nargs='+',
        help='input vcf[.gz] files')
    parser.add_argument('--output', '-o', required=True,
        help='merged ouput file')
    parser.add_argument('--uncompressed', action='store_true', default=False,
        help='do not gzip output (default: False)')
    parser.add_argument('--genotypes', action='store_true', default=False,
        help='include FORMAT and GT columns, if present in input')
    parser.add_argument('--qual_filter_info', action='store_true', default=False,
        help='include QUAL, FILTER, and INFO columns')
    parser.add_argument('--samples', action='store_true', default=False,
        help='create an output file flattened by sample, if present')
    parser.add_argument('--vcf_info', action='store_true', default=False,
        help='include VCF path and TIME columns')
    parser.add_argument('--wsevs', action='store_true', default=False,
        help='input is wsevs NHLBI exome vcf')
    parser.add_argument('--samples_populations',
        help='tab-delimeted file for stratification of samples by gender and population #sample\tsex\tpopulation\tsuper_population')
    args = parser.parse_args()
    
    if not args.uncompressed and not args.output.endswith('.gz'):
        args.output = args.output+'.gz'

    CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,GT = range(10)
    
    #regular expressions for vcf header lines and genotypes
    fileformat_re = re.compile(r'##fileformat=[^0-9.]*(?P<version>[0-9.]+)', re.I)
    filter_re = re.compile(r'##FILTER=<ID=(?P<id>[^,]+),Description="(?P<description>[^"]+)">', re.I)
    format_re = re.compile(r'##FORMAT=<ID=(?P<id>[^,]+),Number=(?P<number>[^,]+),Type=(?P<type>[^,]+),Description="(?P<description>[^"]+)">', re.I)
    info_re = re.compile(r'##INFO=<ID=(?P<id>[^,]+),Number=(?P<number>[^,]+),Type=(?P<type>[^,]+),Description="(?P<description>[^"]+)">', re.I)
    metadata_re = re.compile(r'##(?P<key>[^=]+)="?(?P<value>.+)"?$', re.I)
    gt_re = re.compile(r'(?P<allele1>[0-9.]+)(?P<phase>[/|])(?P<allele2>[0-9.]+)', re.I)
    
    vcfs = my.unglob(args.input)
    
    samples_populations = {}
    populations = set()
    super_populations = set()
    base_stats_cols = ['sample_count', 'allele_count', 'ref_allele_count', 'alt_allele_count', 'di_alt_allele_count',
            'het_sample_count', 'hom_ref_sample_count', 'hom_alt_sample_count', 'di_het_sample_count',]
    stats_cols = ['ref_allele_frequency','alt_allele_frequency'] + base_stats_cols
    if args.samples_populations:
        with open(args.samples_populations) as sp:
            for line in sp:
                if line.startswith("#"):
                    continue
                sample,sex,population,super_population = line.rstrip('\n').split('\t')
                samples_populations[sample] = ('male' if sex == '1' else 'female' if sex == '2' else None,population if population else None,super_population if super_population else None)
                populations.add(population)
                super_populations.add(super_population)
        stats_cols += ['{}_{}'.format(b,s) for b in base_stats_cols for s in ['male','female']]
        stats_cols += ['{}_{}{}'.format(b,p,s) for p in sorted(super_populations)+sorted(populations) for b in base_stats_cols for s in ['','_male','_female']]
    
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
                        raise Vcfs2FlatError('missing fileformat in {}'.format(vcf))
                    elif not float(m.group('version')) >= 4.0:
                        raise Vcfs2FlatError('obsolete fileformat [{}] in {}'.format(line, vcf))
                    else:
                        format_ok = True
                        continue
                elif line.startswith('#'):
                    if line.upper().startswith('#CHROM'):
                        if not line.upper().startswith('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO'):
                            raise Vcfs2FlatError('invalid header format [{}] in {}'.format(line, vcf))
                        elif line.upper().startswith('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'):
                            has_genotypes = True
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
        records_out = 0,
    )
    
    output_cols = ['#CHROM','POS','ID','REF','ALT']
    if args.qual_filter_info:
        output_cols += ['QUAL','FILTER','INFO']
    if has_genotypes and args.genotypes:
        output_cols += ['FORMAT'] + gt_cols
    if args.qual_filter_info:
        output_cols += info_cols
    if has_genotypes:
        output_cols += stats_cols
        if args.samples:
            output_cols += output_cols + ['SAMPLE','ALLELE1','ALLELE2']
    if args.wsevs:
        output_cols += [
                    'AA_AC_ALT', 'AA_AC_REF',
                    'EA_AC_ALT', 'EA_AC_REF',
                    'All_AC_ALT', 'All_AC_REF',
                    'AA_GTC_hom_alt', 'AA_GTC_het', 'AA_GTC_hom_ref', 'AA_GTC_di_alt',
                    'EA_GTC_hom_alt', 'EA_GTC_het', 'EA_GTC_hom_ref', 'EA_GTC_di_alt',
                    'GTC_hom_alt', 'GTC_het', 'GTC_hom_ref', 'GTC_di_alt',
                    'AA_ALT_AF', 'EA_ALT_AF', 'All_ALT_AF',
                    'AA_REF_AF', 'EA_REF_AF', 'All_REF_AF',
                    'AA_MAF', 'EA_MAF', 'All_MAF'
                    ]
    if args.vcf_info:
        output_cols += ['VCF','TIME']
    with gzip.open(args.output, 'w') if not args.uncompressed else open(args.output,'w') as variants_out:
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
                            these_gt_cols = vcf_cols[GT:] if len(vcf_cols) > GT else []
                    #data records
                    else:
                        counters['vcf_records_in'] += 1
                        counters['total_vcf_records_in'] += 1
                        fields = [f.strip() for f in line.split('\t')]
                        chrom, pos, id, ref, alt, qual, filter, info = fields[:FORMAT]
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
                        alt_list = alt.rstrip(',').split(',')
                        these_alleles = [ref] + alt_list
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
                                if args.wsevs and (info_key == 'CG' or info_key == 'CP') and info_value == 'NA':
                                    info_dict[info_key] = ''
                        #flatten output by IDs, ALTs, GTs
                        hom_ref= ref+ref
                        for this_id in ids_list:
                            for a in range(len(alt_list)):
                                this_info_dict = copy.deepcopy(info_dict)
                                this_alt = alt_list[a]
                                hom_alt = this_alt+this_alt
                                het = this_alt+ref
                                ac_list = info_dict.get('AC')
                                if ac_list and len(ac_list)>=a-1:
                                    this_info_dict['AC'] = ac_list[a]
                                af_list = info_dict.get('AF')
                                if af_list and len(af_list)>=a-1:
                                    this_info_dict['AF'] = af_list[a]
                                if args.wsevs:
                                    if len(alt_list)>1:
                                        pass
                                    gts_list = info_dict.get('GTS').split(',')
                                    aa_gtc_list = this_info_dict['AA_GTC'].split(',')
                                    ea_gtc_list = this_info_dict['EA_GTC'].split(',')
                                    gtc_list = this_info_dict['GTC'].split(',')
                                    #check for presence of het for 2 different alt alleles
                                    di_alt_list = [gts_list.index(x) for x in gts_list if x != hom_alt and x != het and this_alt in x]
                                    aa_ac_list = info_dict.get('AA_AC').split(',')
                                    if aa_ac_list:
                                        aa_ac_sum = sum([int(x) for x in aa_ac_list])
                                        aa_ac_min = min([int(x) for x in aa_ac_list])
                                        aa_ac_ref = aa_ac_list[-1]
                                        if len(aa_ac_list)>=a-1:
                                            aa_ac_alt = aa_ac_list[a]
                                        aa_gtc_di_alt = sum([int(aa_gtc_list[i]) for i in di_alt_list])
                                    ea_ac_list = info_dict.get('EA_AC').split(',')
                                    if ea_ac_list:
                                        ea_ac_sum = sum([int(x) for x in ea_ac_list])
                                        ea_ac_min = min([int(x) for x in ea_ac_list])
                                        ea_ac_ref = ea_ac_list[-1]
                                        if len(ea_ac_list)>=a-1:
                                            ea_ac_alt = ea_ac_list[a]
                                        ea_gtc_di_alt = sum([int(ea_gtc_list[i]) for i in di_alt_list])
                                    tac_list = info_dict.get('TAC').split(',')
                                    if tac_list:
                                        tac_sum = sum([int(x) for x in tac_list])
                                        tac_min = min([int(x) for x in tac_list])
                                        tac_ref = tac_list[-1]
                                        if len(tac_list)>=a-1:
                                            tac_alt = tac_list[a]
                                        gtc_di_alt = sum([int(gtc_list[i]) for i in di_alt_list])
                                    aa_gtc_hom_alt = aa_gtc_list[gts_list.index(hom_alt)]
                                    aa_gtc_het = aa_gtc_list[gts_list.index(het)]
                                    aa_gtc_hom_ref = aa_gtc_list[gts_list.index(hom_ref)]
                                    ea_gtc_hom_alt = ea_gtc_list[gts_list.index(hom_alt)]
                                    ea_gtc_het = ea_gtc_list[gts_list.index(het)]
                                    ea_gtc_hom_ref = ea_gtc_list[gts_list.index(hom_ref)]
                                    gtc_hom_alt = gtc_list[gts_list.index(hom_alt)]
                                    gtc_het = gtc_list[gts_list.index(het)]
                                    gtc_hom_ref = gtc_list[gts_list.index(hom_ref)]
                                    #in file maf is *100 and rounded to 4 digits
                                    #and really is the minor allele, not the always alt allele
                                    #ea_maf,aa_maf,all_maf = record_in[header_cols['MAF']].split(',')
                                    aa_alt_af = 0.0 if aa_ac_sum == 0 else float(aa_ac_alt)/aa_ac_sum
                                    ea_alt_af = 0.0 if ea_ac_sum == 0 else float(ea_ac_alt)/ea_ac_sum
                                    all_alt_af = 0.0 if tac_sum == 0 else float(tac_alt)/tac_sum
                                    aa_ref_af = 0.0 if aa_ac_sum == 0 else float(aa_ac_ref)/aa_ac_sum
                                    ea_ref_af = 0.0 if ea_ac_sum == 0 else float(ea_ac_ref)/ea_ac_sum
                                    all_ref_af = 0.0 if tac_sum == 0 else float(tac_ref)/tac_sum
                                    aa_maf = 0.0 if aa_ac_sum == 0 else float(aa_ac_min)/aa_ac_sum
                                    ea_maf = 0.0 if ea_ac_sum == 0 else float(ea_ac_min)/ea_ac_sum
                                    all_maf = 0.0 if tac_sum == 0 else float(tac_min)/tac_sum
                                    wsevs_fields = [str(f) for f in [
                                        aa_ac_alt,aa_ac_ref,
                                        ea_ac_alt,ea_ac_ref,
                                        tac_alt,tac_ref,
                                        aa_gtc_hom_alt,aa_gtc_het,aa_gtc_hom_ref,aa_gtc_di_alt,
                                        ea_gtc_hom_alt,ea_gtc_het,ea_gtc_hom_ref,ea_gtc_di_alt,
                                        gtc_hom_alt,gtc_het,gtc_hom_ref,gtc_di_alt,
                                        aa_alt_af,ea_alt_af,all_alt_af,
                                        aa_ref_af,ea_ref_af,all_ref_af,
                                        aa_maf,ea_maf,all_maf,
                                    ]]
                                    
                                values_out = []
                                base_values_out = [chrom, pos, this_id, ref, this_alt]
                                if args.qual_filter_info:
                                    base_values_out += [qual, filter, info]
                                info_values_out = [this_info_dict[k] for k in sorted(this_info_dict.keys())]
                                gt_values_out = [this_gt_dict[k] for k in sorted(this_gt_dict.keys())]
                                values_out += base_values_out
                                stats = {c:0 for c in stats_cols}
                                sample_count = 0
                                allele_count = 0
                                ref_allele_count = 0
                                alt_allele_count = 0
                                di_alt_allele_count = 0
                                het_sample_count = 0
                                hom_ref_sample_count = 0
                                hom_alt_sample_count = 0
                                di_het_sample_count = 0
                                #allele statics if genotypes present in this file
                                if has_genotypes and these_gt_cols:
                                    for s in range(len(these_gt_cols)):
                                        this_sample = these_gt_cols[s]
                                        sample_population = samples_populations.get(this_sample)
                                        this_gt = gts[s]
                                        m = gt_re.match(this_gt)
                                        allele1x, phase, allele2x = m.groups() if m else ('.','/','.')
                                        allele1, allele2 = (these_alleles[int(allele1x)] if allele1x != '.' else '.', these_alleles[int(allele2x)] if allele2x != '.' else '.')
                                        if allele1 != '.' and allele2 != '.':
                                            #sample_count+=1
                                            increment_stats(stats,'sample_count',1,sample_population)
                                            #allele_count+=2
                                            increment_stats(stats,'allele_count',2,sample_population)
                                            if allele1 == ref and allele2 == ref:
                                                #hom_ref_sample_count+=1
                                                increment_stats(stats,'hom_ref_sample_count',1,sample_population)
                                                #ref_allele_count+=2
                                                increment_stats(stats,'ref_allele_count',2,sample_population)
                                            elif allele1 == this_alt and allele2 == this_alt:
                                                #hom_alt_sample_count+=1
                                                increment_stats(stats,'hom_alt_sample_count',1,sample_population)
                                                #alt_allele_count+=2
                                                increment_stats(stats,'alt_allele_count',2,sample_population)
                                            elif (allele1 == this_alt and allele2 == ref) or (allele2 == this_alt and allele1 == ref):
                                                #het_sample_count+=1
                                                increment_stats(stats,'het_sample_count',1,sample_population)
                                                #alt_allele_count+=1
                                                increment_stats(stats,'alt_allele_count',1,sample_population)
                                                #ref_allele_count+=1
                                                increment_stats(stats,'ref_allele_count',1,sample_population)
                                            elif (allele1 == this_alt and allele2 != this_alt and allele2 != ref) or (allele2 == this_alt and allele1 != this_alt and allele1 != ref):
                                                #di_het_sample_count+=1
                                                increment_stats(stats,'di_het_sample_count',1,sample_population)
                                                #alt_allele_count+=1
                                                increment_stats(stats,'alt_allele_count',1,sample_population)
                                                #di_alt_allele_count+=1
                                                increment_stats(stats,'di_alt_allele_count',1,sample_population)
                                    stats['ref_allele_frequency'] = float(stats['ref_allele_count'])/float(stats['allele_count'])
                                    stats['alt_allele_frequency'] = float(stats['alt_allele_count'])/float(stats['allele_count'])

                                if has_genotypes and args.genotypes:
                                    values_out += [format] + gt_values_out
                                if args.qual_filter_info:
                                    values_out += info_values_out
                                if has_genotypes:
                                    values_out += [str(stats[c]) for c in stats_cols]
                                if args.wsevs:
                                    values_out += wsevs_fields
                                if args.vcf_info:
                                    output_cols += [vcf, vcf_time]
                                if not args.samples:
                                    line_out = '\t'.join(values_out)
                                    variants_out.write(line_out+'\n')
                                    counters['records_out'] += 1
                                else:
                                    if not has_genotypes or not these_gt_cols:
                                        values_out += output_cols + ['','','']
                                        line_out = '\t'.join(values_out)
                                        variants_out.write(line_out+'\n')
                                        counters['records_out'] += 1
                                    else:
                                        for s in range(len(these_gt_cols)):
                                            this_sample = these_gt_cols[s]
                                            this_gt = gts[s]
                                            m = gt_re.match(this_gt)
                                            allele1x, phase, allele2x = m.groups() if m else ('.','/','.')
                                            allele1, allele2 = (these_alleles[int(allele1x)] if allele1x != '.' else '.', these_alleles[int(allele2x)] if allele2x != '.' else '.')
                                            values_out += output_cols + [this_sample, allele1, allele2]
                                            line_out = '\t'.join(values_out)
                                            variants_out.write(line_out+'\n')
                                            counters['records_out'] += 1
            print 'input vcf records: {}'.format(counters['vcf_records_in'])
    print 'done. vcf files processed: {}; total input vcf records: {}; total output records: {}'.format(counters['vcf_files'],counters['total_vcf_records_in'],counters['records_out'])

def increment_stats(stats,stat,value,sample_population):   
    stats[stat]+=value
    if sample_population:
        sex,population,super_population = sample_population
        if sex:
            stats[stat+'_'+sex]+=value
        if population:
            stats[stat+'_'+population]+=value
            if sex:
                stats[stat+'_'+population+'_'+sex]+=value
        if super_population:
            stats[stat+'_'+super_population]+=value
            if sex:
                stats[stat+'_'+super_population+'_'+sex]+=value
    

if __name__ == "__main__": sys.exit(main())
