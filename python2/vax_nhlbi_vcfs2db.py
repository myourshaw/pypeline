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
import tarfile
import copy
import xml.sax.saxutils
import collections
import my
import sql_columns
try:
    import MySQLdb
    import MySQLdb.cursors
except ImportError:
    raise ImportError,"The Python MySQLdb module is required to run this program. Try 'pip install MySQL-python'."
try:
    import sh
    from sh import git
except ImportError:
    raise ImportError,"The Python sh module is required to run this program. Try 'pip install sh'."

class NhlbiVcfs2DbError(Exception): pass

#-i /scratch1/vax/75/nhlbi/ESP6500SI-V2-SSA137.protein-hgvs-update.snps_indels.vcf.tar.gz -o /scratch1/vax/75/nhlbi/ESP6500SI-V2-SSA137.protein-hgvs-update.snps_indels.vcf -d vax_test -p m.cha3ly

class NhlbiVcfs2Db:
    log_file = ''
    has_genotypes = False
    format_ok = False
    info_cols = set()
    gt_cols = set()
    metadata_done = False
    counters = collections.Counter()
    these_gt_cols = None

    db_out = None
    genotypes_out = None
    alleles_out = None
    
    #regular expressions for vcf header lines and genotypes
    fileformat_re = re.compile(r'##fileformat=[^0-9.]*(?P<version>[0-9.]+)', re.I)
    filter_re = re.compile(r'##FILTER=<ID=(?P<id>[^,]+),Description="(?P<description>[^"]+)">', re.I)
    format_re = re.compile(r'##FORMAT=<ID=(?P<id>[^,]+),Number=(?P<number>[^,]+),Type=(?P<type>[^,]+),Description="(?P<description>[^"]+)">', re.I)
    info_re = re.compile(r'##INFO=<ID=(?P<id>[^,]+),Number=(?P<number>[^,]+),Type=(?P<type>[^,]+),Description="(?P<description>[^"]+)">', re.I)
    metadata_re = re.compile(r'##(?P<key>[^=]+)="?(?P<value>.+)"?$', re.I)
    gt_re = re.compile(r'(?P<allele1>[0-9.]+)(?P<phase>[/|])(?P<allele2>[0-9.]+)', re.I)
    indel_gt_re = re.compile(r'(?P<AR1>[AR])(?P<IX1>\d*)(?P<AR2>[AR]*)(?P<IX2>\d*)',re.I)
    chr_re = re.compile(r'\.chr(?P<chrom>[0-9XY]{1,2})\.', re.I)
    #mySQL names
    forbidden = re.compile(r'[^0-9,a-z,A-Z$_]')
    
    CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,GT = range(10)
    
    def print_log(self, s):
        my.print_log(s, self.log_file)
        
    def run(self, input, output, compress=False, AApop='ASW', AAsuper='AFR', EApop='CEU', EAsuper='EUR',
            do_not_install_mysql_database=False, database=None, table=None, config=None, host=None, port=None, user=None, password=None,):
        self.AApop=AApop
        self.AAsuper=AAsuper
        self.EApop=EApop
        self.EAsuper=EAsuper
        
        self.log_file = output + '.log'
        try:
            os.remove(self.log_file)
        except OSError:
            pass
        
        if do_not_install_mysql_database:
            db_output_file = output + '_db.txt' + ('.gz' if compress else '')
            allele_output_file = output + '_allele_counts.txt' + ('.gz' if compress else '')
            genotype_output_file = output + '_genotype_counts.txt' + ('.gz' if compress else '')
        else:
            self.compress = False
            freq_table = table if table else my.r_strip(os.path.basename(output), '.vcf')
            freq_table = self.forbidden.sub('_',my.r_strip(freq_table, '.txt'))[:64]
            db_output_file = os.path.join(os.path.dirname(output), freq_table+'.txt')
            allele_count_table = freq_table+'_allele_counts'
            allele_count_table = self.forbidden.sub('_',my.r_strip(allele_count_table, '.txt'))[:64]
            allele_output_file = os.path.join(os.path.dirname(output), allele_count_table+'.txt')
            genotype_count_table = freq_table+'_genotype_counts'
            genotype_count_table = self.forbidden.sub('_',my.r_strip(genotype_count_table, '.txt'))[:64]
            genotype_output_file = os.path.join(os.path.dirname(output), genotype_count_table+'.txt')
            
        done_file = my.done_file(output)
        if my.file_exists(done_file):
            self.print_log('NHLBI-EVS file(s) already processed. To reprocess "rm {}"'.format(done_file))
        else:
            self.print_log('Creating NHLBI-EVS database text file(s). This will take ~3 hours.')
        
            self.print_log('DB output file: {}.'.format(db_output_file))
            self.print_log('Allele output file: {}.'.format(allele_output_file))
            self.print_log('Genotype output file: {}.'.format(genotype_output_file))
            self.db_out = gzip.open(db_output_file, 'w') if compress else open(db_output_file,'w')
            self.genotypes_out= gzip.open(genotype_output_file, 'w') if compress else open(genotype_output_file,'w')
            self.alleles_out = gzip.open(allele_output_file, 'w') if compress else open(allele_output_file,'w')
            
            inputs = my.unglob(input, sort_unique=False)
    
            #minimal format validation and get all info and gt columns
            #only metadata from the first vcf file is used
            for input_file in inputs:
                self.counters['input_files'] += 1
                if tarfile.is_tarfile(input_file):
                    self.print_log('Extracting data from tar file {}'.format(input_file))
                    with tarfile.open(input_file) as tar:
                        members = tar.getmembers()
                        names = [m.name for m in members]
                        #ESP6500SI-V2-SSA137.updatedProteinHgvs.chr1.snps_indels.vcf
                        file_order = {}
                        for member in members:
                            chr_m = self.chr_re.search(member.name)
                            if not chr_m:
                                warnings.warn('Missing chromosome in file name {}'.format(member.name))
                            else:
                                chrom = 23 if chr_m.group('chrom') == 'X' else 24 if chr_m.group('chrom') == 'Y' else int(chr_m.group('chrom'))
                                file_order[chrom] = member
                        for i in sorted(file_order.keys()):
                            this_member = file_order[i]
                            vcf = tar.extractfile(this_member)
                            self.these_gt_cols = None
                            self.counters['vcf_files'] += 1
                            self.counters['vcf_lines_in'] = 0
                            self.counters['vcf_records_in'] = 0
                            self.print_log('Processing data in {}'.format(this_member.name))
                            vcf_line_count = 0
                            for line in vcf:
                                self.process_line(line)
                            vcf.close()
                            self.print_log('VCF input file {} had {} lines, {} vcf records.'.format(this_member.name, self.counters['vcf_lines_in'], self.counters['vcf_records_in']))
        
                else:
                    with my.open_gz_or_text(input_file) as vcf:
                        self.counters['vcf_files'] += 1
                        self.counters['vcf_lines_in'] = 0
                        self.counters['vcf_records_in'] = 0
                        self.print_log('reading {}'.format(input_file))
                        self.these_gt_cols = None
                        for line in vcf:
                            self.process_line(line)
                    self.print_log('VCF input file {} had {} lines, {} vcf records.'.format(input_file, self.counters['vcf_lines_in'], self.counters['vcf_records_in']))
                    
            with open(done_file, 'w'):
                pass
            self.print_log('Created NHLBI-EVS database text file(s). vcf files processed: {}; total input vcf records: {}; total genotype records: {}; total allele records: {}'.format(self.counters['vcf_files'],self.counters['total_vcf_records_in'],self.counters['genotypes_out'],self.counters['alleles_out']))

        #install on mySQL database
        if do_not_install_mysql_database:
            self.print_log('Skipping mySQL database installation.')
        else:
            #create mySQL tables
            
            done_file = my.done_file(output+'.database')
            if my.file_exists(done_file):
                self.print_log('NHLBI-EVS database tables and procedures already installed. To reinstall "rm {}"'.format(done_file))
            else:
                self.print_log('Installing NHLBI-EVS database tables and procedures. This will take ~30 hours.')
                self.install_database(db_output_file, freq_table, allele_output_file, allele_count_table, genotype_output_file, genotype_count_table,
                                 database, host, port, user, password,)
                
            with open(done_file, 'w'):
                pass
            self.print_log('Installed NHLBI-EVS database tables and procedures.')
        return True
    
    def process_line(self, line):
        self.counters['vcf_lines_in'] += 1
        self.counters['total_vcf_lines_in'] += 1
        line = line.rstrip('\n')
        if not bool(line.strip()):
            return
        if not self.metadata_done:
            if self.process_metadata(line):
                return
            else:
                self.metadata_done = True
                self.info_cols = [i+'_info' if i in ['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT'] else i for i in sorted(self.info_cols)]
                self.gt_cols = sorted(self.gt_cols)
                self.write_headers()
        else:
            #meta-information records
            if line.startswith('#'):
                if line.upper().startswith('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO'):
                    vcf_cols = line.split('\t')
                    self.these_gt_cols = vcf_cols[self.GT:] if len(vcf_cols) > self.GT else []
            else:
                self.write_data(line)
        
        
    def process_metadata(self, line):
        if not self.format_ok:
            m = self.fileformat_re.match(line)
            if not m:
                warnings.warn('missing fileformat in {}'.format(vcf))
            elif not float(m.group('version')) >= 4.0:
                warnings.warn('obsolete fileformat [{}] in {}'.format(line, vcf))
            else:
                self.format_ok = True
                return True
        elif line.startswith('#'):
            if line.upper().startswith('#CHROM'):
                if not line.upper().startswith('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO'):
                    raise NhlbiVcfs2DbError('invalid header format [{}] in {}'.format(line, vcf))
                elif line.upper().startswith('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'):
                    self.has_genotypes = True
                    warnings.warn('input file {} has genotypes'.format(vcf))
                    self.gt_cols = self.gt_cols | set(line.split('\t')[self.GT:])
                return False
            else:
                info_m = self.info_re.match(line)
                if info_m:
                    self.info_cols.add(info_m.group('id'))
                return True
    
    def write_headers(self):
        db_output_cols = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','allele1','allele2','oALT_NHLBI','EA_AC_allele1','EA_AC_allele2','AA_AC_allele1','AA_AC_allele2','TAC_allele1','TAC_allele2','MAF_EA','MAF_AA','MAF_All','EA_GTC_this','AA_GTC_this','GTC_this'] + self.info_cols
        genotype_output_cols = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','allele1','allele2','oALT_NHLBI','population','super_population','sample_count']
        allele_output_cols = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','allele','oALT_NHLBI','population','super_population','allele_count']
        self.db_out.write('\t'.join(db_output_cols)+'\n')
        self.genotypes_out.write('\t'.join(genotype_output_cols)+'\n')
        self.alleles_out.write('\t'.join(allele_output_cols)+'\n')
    
    def write_data(self, line):
        self.counters['vcf_records_in'] += 1
        self.counters['total_vcf_records_in'] += 1
        fields = [f.strip() for f in line.split('\t')]
        chrom, pos, id, ref, oalt, qual, filter, info = fields[:self.FORMAT]
        chrom = chrom.upper()
        #unescape things like &amp;
        info = xml.sax.saxutils.unescape(info)
        format = fields[self.FORMAT] if len(fields) > self.FORMAT else '.'
        gts = fields[self.GT:] if len(fields) > self.GT else None
        if gts:
            this_gt_dict = {g:'' for g in self.these_gt_cols}
            for i in range(len(self.these_gt_cols)):
                this_gt_dict[self.these_gt_cols[i]] = gts[i]
#TODO: VCF 4.1 allows additional info in ALT and ID
        ids_list = id.rstrip(';').split(';')
        filters_list = filter.rstrip(';').split(';')
        alt_list = oalt.rstrip(',').split(',')
        these_alleles = [ref] + alt_list
        alts_ref = alt_list + [ref]
        #oALT may have multiple ,-separated alleles
        info = 'oALT='+oalt+';'+info if info != '.' else 'oALT='+oalt+';'
        info_list = info.rstrip(';').split(';')
        info_dict = {i:'' for i in self.info_cols}
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
        info_data = [info_dict[i] for i in self.info_cols]
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
                indel_gt_m = self.indel_gt_re.match(gts_list[i])
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
                self.db_out.write('\t'.join(db_out_list)+'\n')
            for population in [(aa_gtc_list,self.AApop,self.AAsuper), (ea_gtc_list,self.EApop,self.EAsuper)]:
                for i in range(len(gts_list)):
                    allele1, allele2 = gts_list[i] if len(gts_list[i]) == 2 else [gts_list[i][0],'']
                    pop,super_pop = population[1:3]
                    sample_count = population[0][i]
                    genotype_out_list = out_list + [allele1, allele2, oalt, pop, super_pop, sample_count]
                    self.genotypes_out.write('\t'.join(genotype_out_list)+'\n')
                    self.counters['genotypes_out']+=1
            for population in [(aa_ac_list,self.AApop,self.AAsuper), (ea_ac_list,self.EApop,self.EAsuper)]:
                for i in range(len(alts_ref)):
                    allele = alts_ref[i]
                    pop,super_pop = population[1:3]
                    sample_count = population[0][i]
                    allele_out_list = out_list + [allele, oalt, pop, super_pop, sample_count]
                    self.alleles_out.write('\t'.join(allele_out_list)+'\n')
                    self.counters['alleles_out']+=1

    def install_database(self, db_output_file, freq_table, allele_output_file, allele_count_table, genotype_output_file, genotype_count_table,
                         database, host, port, user, password,):
        #freq table
        text_file = db_output_file
        table = freq_table
        sql_file = text_file+'.mysql'
        sql_done_file = my.done_file(sql_file)
        if my.file_exists(sql_done_file):
            self.print_log('Script {} already created. To redo "rm {}"'.format(sql_file, sql_done_file))
        else:
            self.print_log('Making CREATE TABLE `{}` script. This will take ~6 hours..'.format(table))
            sql_columns.run(text_file, database=database, schema=database, table=table, clustered_index=['CHROM', 'POS', 'REF', 'ALT','allele1','allele2'], indexes=['ID'])
            with open(sql_done_file, 'w'):
                pass
            self.print_log('Made CREATE TABLE `{}` script.'.format(table))
            
        self.print_log('Creating `{}` table.'.format(table))
        with open(sql_file, 'r') as sql:
            sh.mysql('-h', host, '-P', port, '-u', user, '-p{}'.format(password), database, _in=sql)
        self.print_log('Created `{}` table.'.format(table))
        
        #import freq table data
        self.print_log('Importing data to `{}`.`{}` table. This will take ~3 minutes.'.format(database, table))
        sh.mysqlimport('-h', host, '-P', port, '-u', user, '-p{}'.format(password), '--delete', '--local', '--ignore-lines=1', database, text_file)
        self.print_log('Imported data to `{}`.`{}` table.'.format(database, table))
    
        #allele_counts table
        text_file = allele_output_file
        table = allele_count_table
        sql_file = text_file+'.mysql'
        sql_done_file = my.done_file(sql_file)
        if my.file_exists(sql_done_file):
            self.print_log('{} already created. To redo "rm {}"'.format(sql_file, sql_done_file))
        else:
            self.print_log('Making CREATE TABLE `{}` script. This will take ~2 hours.'.format(table))
            sql_columns.run(text_file, database=database, schema=database, table=table, clustered_index=['CHROM', 'POS', 'REF', 'ALT','allele','population'])
            with open(sql_done_file, 'w'):
                pass
            self.print_log('Made CREATE TABLE `{}` script.'.format(table))
            
        self.print_log('Creating `{}` table.'.format(table))
        with open(sql_file, 'r') as sql:
            sh.mysql('-h', host, '-P', port, '-u', user, '-p{}'.format(password), database, _in=sql)
        self.print_log('Created `{}` table.'.format(table))
        
        #import allele_counts table data
        self.print_log('Importing data to `{}`.`{}` table. This will take ~2 minutes.'.format(database, table))
        sh.mysqlimport('-h', host, '-P', port, '-u', user, '-p{}'.format(password), '--delete', '--local', '--ignore-lines=1', database, text_file)
        self.print_log('Imported data to `{}`.`{}` table.'.format(database, table))
    
        #allele_counts_all table
        #allele counts by population and also by AA+EA combined
        allele_counts_all_table = allele_count_table[:60]+'_all'
        table = allele_counts_all_table
        self.print_log('Creating `{}` table. This will take ~15 minutes.'.format(table))
        query = [
'DROP TABLE IF EXISTS `{0}`.`{1}`'.format(database, table),
"""CREATE TABLE `{0}`.`{1}`
AS
SELECT CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, allele, oALT_NHLBI, 'CEU+ASW' AS population, 'EUR+AFR' AS super_population, SUM(allele_count) AS allele_count
FROM `{0}`.`{2}`
GROUP BY CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, allele, oALT_NHLBI
UNION
SELECT *
FROM `{0}`.`{2}`
""".format(database, table, allele_count_table)
]
        with MySQLdb.connect(host,user,password,database,port=port) as cursor:
            my.mysql(cursor, query)
        self.print_log('Created `{}` table.'.format(table))
    
        #genotype_counts table
        text_file = genotype_output_file
        table = genotype_count_table
        sql_file = text_file+'.mysql'
        sql_done_file = my.done_file(sql_file)
        if my.file_exists(sql_done_file):
            self.print_log('Script {} already created. To redo "rm {}"'.format(sql_file, sql_done_file))
        else:
            self.print_log('Making CREATE TABLE `{}` script. This will take ~20 hours..'.format(table))
            sql_columns.run(text_file, database=database, schema=database, table=table, clustered_index=['CHROM', 'POS', 'REF', 'ALT','allele1','allele2', 'population'])
            with open(sql_done_file, 'w'):
                pass
            self.print_log('Made CREATE TABLE `{}` script.'.format(table))
            
        self.print_log('Creating `{}` table.'.format(table))
        with open(sql_file, 'r') as sql:
            sh.mysql('-h', host, '-P', port, '-u', user, '-p{}'.format(password), database, _in=sql)
        self.print_log('Created `{}` table.'.format(table))
        
        #import genotype_counts table data
        self.print_log('Importing data to `{}`.`{}` table. This will take ~3 minutes.'.format(database, table))
        sh.mysqlimport('-h', host, '-P', port, '-u', user, '-p{}'.format(password), '--delete', '--local', '--ignore-lines=1', database, text_file)
        self.print_log('Imported data to `{}`.`{}` table.'.format(database, table))
        
        #genotype_counts_all table
        #genotype counts by population and also by AA+EA combined
        genotype_counts_all_table = genotype_count_table[:60]+'_all'
        table = genotype_counts_all_table
        self.print_log('Creating `{}` table. This will take ~25 minutes.'.format(table))
        query = [
'DROP TABLE IF EXISTS `{0}`.`{1}`'.format(database, table),
"""CREATE TABLE `{0}`.`{1}`
AS
SELECT CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, allele1, allele2, oALT_NHLBI, 'CEU+ASW' AS population, 'EUR+AFR' AS super_population, sum(sample_count) AS sample_count
FROM `{0}`.`{2}`
GROUP BY CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, allele1, allele2, oALT_NHLBI
UNION
SELECT *
FROM `{0}`.`{2}`
""".format(database, table, genotype_count_table)
]
        with MySQLdb.connect(host,user,password,database,port=port) as cursor:
           my.mysql(cursor, query)
        self.print_log('Created `{}` table.'.format(table))
       
        #gtc_sum table
        #genotype count sums for freq table
        gtc_sum_table = freq_table[:56]+'_gtc_sum'
        table = gtc_sum_table
        self.print_log('Creating `{}` table. This will take ~3 minutes.'.format(table))
        query = [
'DROP TABLE IF EXISTS `{0}`.`{1}`'.format(database, table),
"""CREATE TABLE `{0}`.`{1}`
AS
SELECT CHROM,POS,REF,ALT
,SUM(EA_GTC_this) AS EA_GTC_total
,SUM(AA_GTC_this) AS AA_GTC_total
,SUM(GTC_this) AS GTC_total
FROM `{0}`.`{2}`
GROUP BY CHROM,POS,REF,ALT
""".format(database, table, freq_table)
]
        with MySQLdb.connect(host,user,password,database,port=port) as cursor:
           my.mysql(cursor, query)
        self.print_log('Created `{}` table.'.format(table))
       
        #vax table
        #add total genotype counts to freq table
        vax_table = freq_table[:60]+'_vax'
        table = vax_table
        self.print_log('Creating `{}` table. This will take ~8 minutes.'.format(table))
        query = [
'DROP TABLE IF EXISTS `{0}`.`{1}`'.format(database, table),
"""CREATE TABLE `{0}`.`{1}`
AS
SELECT e.*
,x.EA_GTC_total, x.AA_GTC_total,x.GTC_total
,CASE WHEN x.EA_GTC_total = 0 THEN NULL ELSE CAST(e.EA_AC_allele1/(2.0*x.EA_GTC_total) AS DECIMAL(10,6)) END AS EA_allele1_freq
,CASE WHEN x.EA_GTC_total = 0 THEN NULL ELSE CAST(e.EA_AC_allele2/(2.0*x.EA_GTC_total) AS DECIMAL(10,6)) END AS EA_allele2_freq
,CASE WHEN x.AA_GTC_total = 0 THEN NULL ELSE CAST(e.AA_AC_allele1/(2.0*x.AA_GTC_total) AS DECIMAL(10,6)) END AS AA_allele1_freq
,CASE WHEN x.AA_GTC_total = 0 THEN NULL ELSE CAST(e.AA_AC_allele2/(2.0*x.AA_GTC_total) AS DECIMAL(10,6)) END AS AA_allele2_freq
,CASE WHEN x.GTC_total = 0 THEN NULL ELSE CAST(e.TAC_allele1/(2.0*x.GTC_total) AS DECIMAL(10,6)) END AS EA_AA_allele1_freq
,CASE WHEN x.GTC_total = 0 THEN NULL ELSE CAST(e.TAC_allele2/(2.0*x.GTC_total) AS DECIMAL(10,6)) END AS EA_AA_allele2_freq
,CASE WHEN x.EA_GTC_total = 0 THEN NULL ELSE CAST(e.EA_GTC_this/x.EA_GTC_total AS DECIMAL(10,6)) END AS EA_GT_freq
,CASE WHEN x.AA_GTC_total = 0 THEN NULL ELSE CAST(e.AA_GTC_this/x.AA_GTC_total AS DECIMAL(10,6)) END AS AA_GT_freq
,CASE WHEN x.GTC_total = 0 THEN NULL ELSE CAST(e.GTC_this/x.GTC_total AS DECIMAL(10,6)) END AS EA_AA_GT_freq
FROM `vax_test`.`ESP6500SI_V2_SSA137_protein_hgvs_update_snps_indels` AS e
JOIN `vax_test`.`ESP6500SI_V2_SSA137_protein_hgvs_update_snps_indels_gtc_sum` AS x
ON e.chrom=x.chrom and e.pos=x.pos and e.ref=x.ref and e.alt=x.alt
""".format(database, table, freq_table, gtc_sum_table),
"""ALTER TABLE `{0}`.`{1}`
ADD PRIMARY KEY `IX_CHROM_POS_REF_ALT_allele1_allele2` (`CHROM` ASC, `POS` ASC, `REF` ASC, `ALT` ASC, `allele1` ASC, `allele2` ASC),
ADD INDEX `IX_ID` (`ID` ASC)
""".format(database, table),
]
        with MySQLdb.connect(host,user,password,database,port=port) as cursor:
           my.mysql(cursor, query)
        self.print_log('Created `{}` table.'.format(table))
       
        #create stored procedures
        
        #get_nhlbi_evs
        #CALL `get_nhlbi_evs`('9', 37783990, 'T', 'G');
        proc = 'get_nhlbi_evs'
        self.print_log('Creating `{}`.`{}` stored procedure'.format(database, proc))
        query = [
'DROP PROCEDURE IF EXISTS `{0}`.`{1}`'.format(database, proc),
"""CREATE DEFINER=CURRENT_USER PROCEDURE `{0}`.`{1}`(CHROM VARCHAR(2), POS INT(11), REF VARCHAR(255), ALT VARCHAR(255))
BEGIN
SELECT *
FROM `{0}`.`{2}` AS n
WHERE n.`CHROM` = CHROM
AND n.`POS` = POS
AND n.`REF` = REF
AND n.`ALT` = ALT
AND n.`allele1` = ref AND n.`allele2` = alt
;
END
""".format(database, proc, vax_table)
]
        with MySQLdb.connect(host,user,password,database,port=port) as cursor:
            my.mysql(cursor, query)
        self.print_log('Created `{}`.`{}` stored procedure'.format(database, proc))
            
        #get_nhlbi_evs_by_locus
        #CALL `vax_test`.`get_nhlbi_evs_by_locus`('9', 37783990);
        proc = 'get_nhlbi_evs_by_locus'
        self.print_log('Creating `{}`.`{}` stored procedure'.format(database, proc))
        query = [
'DROP PROCEDURE IF EXISTS `{0}`.`{1}`'.format(database, proc),
"""CREATE DEFINER=CURRENT_USER PROCEDURE `{0}`.`{1}`(CHROM VARCHAR(2), POS INT(11))
BEGIN
SELECT *
FROM `{0}`.`{2}` AS n
WHERE n.`CHROM` = CHROM
AND n.`POS` = POS
;
END
""".format(database, proc, vax_table)
]
        with MySQLdb.connect(host,user,password,database,port=port) as cursor:
            my.mysql(cursor, query)
        self.print_log('Created `{}`.`{}` stored procedure'.format(database, proc))
            
        #get_nhlbi_evs_by_ID
        #CALL `vax_test`.`get_nhlbi_evs_by_ID`('rs141138948');
        proc = 'get_nhlbi_evs_by_ID'
        self.print_log('Creating `{}`.`{}` stored procedure'.format(database, proc))
        query = [
'DROP PROCEDURE IF EXISTS `{0}`.`{1}`'.format(database, proc),
"""CREATE DEFINER=CURRENT_USER PROCEDURE `{0}`.`{1}`(ID VARCHAR(15))
BEGIN
SELECT *
FROM `{0}`.`{2}` AS n
WHERE n.`ID` = ID
;
END
""".format(database, proc, vax_table)
]
        with MySQLdb.connect(host,user,password,database,port=port) as cursor:
            my.mysql(cursor, query)
        self.print_log('Created `{}`.`{}` stored procedure'.format(database, proc))
    
        return True
    
def main():

    #command line arguments
    parser = argparse.ArgumentParser(
        description = 'create a database/excel-friendly vcf file by expanding alleles and genotypes',
        epilog = 'pypeline.nhlbi_vcfs2db version 1.0β1 ©2011-2014 Michael Yourshaw all rights reserved')
    parser.add_argument('--input', '-i', nargs='+',
        help='input vcf[.gz] files')
    parser.add_argument('--output', '-o', required=True,
        help='prefix for merged ouput files')
    parser.add_argument('--compress', action='store_true', default=False,
        help='gzip output (default: False)')
    parser.add_argument('--AApop', default='ASW',
        help='African American population ID (default: ASW)')
    parser.add_argument('--AAsuper', default='AFR',
        help='African American super-population ID (default: AFR)')
    parser.add_argument('--EApop', default='CEU',
        help='European American population ID (default: CEU)')
    parser.add_argument('--EAsuper', default='EUR',
        help='European American super-population ID (default: EUR)')
    #databases
    parser.add_argument('--do_not_install_mysql_database', action='store_true', default=False,
        help='Do not install mySQL database (default: False)')
    parser.add_argument('--database', '-d', type=str,
        help='MySQL database')
    parser.add_argument('--table', type=str,
        help='MySQL table')
    parser.add_argument('--host', '-H', type=str,
        help='MySQL database server hostname or ip address (default: config.db.HOST, example: cortex.local)')
    parser.add_argument('--port', '-P', type=int,
        help='MySQL database (default: config.db.PORT, example: 3306)')
    parser.add_argument('--user', '-u', type=str,
        help='MySQL database user with priveleges to install Ensembl (default: config.db.USER, example: sa)')
    parser.add_argument('--password', '-p', type=str,
        help='MySQL password for installing user (default: config.db.USER, example: None = enter silently at prompt)')
    args = parser.parse_args()
    config = my.get_config(args, 'vax.cfg')

    if not args.do_not_install_mysql_database:
        if not args.database:
            try:
                args.database = raw_input('MySQL vax database for NHLBI/EVS tables: ')
            except:
                raise NhlbiVcfs2DbError('A database is required, unless do_not_install_mysql_database is specified. Installer aborted.')
            else:
                #empty database
                if not args.database:
                    raise NhlbiVcfs2DbError('A database is required, unless do_not_install_mysql_database is specified. Installer aborted.')
        #MySQL install user and password
        if not args.user:
            try:
                args.user = config.get('db','USER')
            except:
                pass
            else:
                if not args.user:
                    try:
                        args.user = raw_input('MySQL user with install permissions: ')
                    except:
                        raise NhlbiVcfs2DbError('A user is required. Installer aborted.')
                    else:
                        #empty user
                        if not args.user:
                            raise NhlbiVcfs2DbError('A user is required. Installer aborted.')
        if not args.password:
            try:
                args.password = config.get('db','PASSWORD')
            except:
                pass
            else:
                if not args.password:
                    try:
                        args.password = getpass('MySQL password with install permissions: ')
                    except:
                        raise NhlbiVcfs2DbError('A password is required. Installer aborted.')
                    else:
                        #empty password
                        if not args.password:
                            raise NhlbiVcfs2DbError('A password is required. Installer aborted.')
        if not args.host:
            try:
                args.host = config.get('db','HOST')
            except:
                pass
            else:
                if not args.host:
                    try:
                        args.host = raw_input('MySQL host (e.g., 127.0.0.1, localhost, or cortex.local): ')
                    except:
                        raise NhlbiVcfs2DbError('A host is required. Installer aborted.')
                    else:
                        #empty host
                        if not args.host:
                            raise NhlbiVcfs2DbError('A host is required. Installer aborted.')
        if not args.port:
            try:
                args.port = int(config.get('db','PORT'))
            except:
                pass
            else:
                if not args.port:
                    try:
                        args.port = int(raw_input('MySQL port (e.g., 3306: '))
                    except:
                        raise NhlbiVcfs2DbError('A port is required. Installer aborted.')
                    else:
                        #empty port
                        if not args.port:
                            raise NhlbiVcfs2DbError('A port is required. Installer aborted.')    
    
    NHLBI = NhlbiVcfs2Db()
    NHLBI.run(input=args.input, output=args.output, compress=args.compress,
        AApop=args.AApop, AAsuper=args.AAsuper, EApop=args.EApop, EAsuper=args.EAsuper,
        do_not_install_mysql_database=args.do_not_install_mysql_database, database=args.database,
        table=args.table, password=args.password, config=config, host=args.host, port=args.port, user=args.user,)
    
if __name__ == "__main__": sys.exit(main())
