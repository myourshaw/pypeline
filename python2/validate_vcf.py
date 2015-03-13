#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import argparse
from warnings import warn
import re
import my

class ValidateVcfError(Exception): pass

#-v /home/myourshaw/lab/pypeline/vax23/gmd_test.vcf

def write_error(err_count, max_errors, err_msg, output_file):
    err_count+=1
    if err_count > max_errors:
        raise ValidateVcfError('VCF file has > {} errors; see {}'.format(max_errors, output_file.name))
    else:
        output_file.write(err_msg+'\n')
        return err_count

def run(input, output=None, max_errors=1):
    if not output:
        output = input+'.validate'

    CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT = range(9)
    valid_types = ('Integer', 'Float', 'Flag', 'Character', 'String')
    valid_alt_ids = ('DEL', 'INS', 'DUP', 'INV', 'CNV')

    gt_re = re.compile(r'(?P<allele1>[0-9.]+)(?P<phase>[/|])?(?P<allele2>[0-9.]+)?', re.I)
    ws_re = re.compile(r'[\s]', re.I)
    ws_colon_re = re.compile(r'[\s:]', re.I)
    ws_semicolon_re = re.compile(r'[\s;]', re.I)
    ws_comma_angle_brackets_re = re.compile(r'[\s,<>]', re.I)
    ws_semicolon_equal_sign_re = re.compile(r'[\s;=]', re.I)
    not_ACGTN_re = re.compile(r'[^ACGTN]', re.I)
    not_ACGTNX_re = re.compile(r'[^ACGTNX]', re.I)

    err_count = 0
    line_count = 0
    header_seen = False
    header_cols = None
    column_count = 0
    version = None
    that_chrom = None
    that_pos = 0
    data_seen = False
    
    output_file = open(output, 'w')
    try:
        with open(input) as v:
            for line in v:
                line_count += 1
                line = line.rstrip('\n')
                if line.startswith('#'):
                    if line_count == 1:
                        if not line.startswith('##fileformat='): ##fileformat=VCFv4.1
                            err_msg = 'line {} missing fileformat {}'.format(line_count, line)
                            err_count = write_error(err_count, max_errors, err_msg, output_file)
                        else:
                            version = line[17:].split('.')
                            if not my.is_int(version[0]) or int(version[0]) < 4:
                                warn('version < 4')
                    if header_seen:
                        err_msg = 'line {} comment after header {}'.format(line_count, line)
                        err_count = write_error(err_count, max_errors, err_msg, output_file)
                    elif line.startswith('##INFO'):
                        if not (line.startswith('##INFO=<') and line.endswith('>')):
                            err_msg = 'line {} invalid ##INFO header {}'.format(line_count, line)
                            err_count = write_error(err_count, max_errors, err_msg, output_file)
                        else:
                            info = line[8:-1]
                            info = info.split(',',3)
                            if not(len(info) == 4
                                   and info[0].startswith('ID=')
                                   and info[1].startswith('Number=') and len(info[1])>7 and (info[1][7:] in ('.','A','G') or my.is_int(info[1][7:]))
                                   and info[2].startswith('Type=') and len(info[2])>5 and info[2][5:] in valid_types
                                   and info[3].startswith('Description="') and info[3].endswith('"')):
                                err_msg = 'line {} invalid ##INFO header {}'.format(line_count, line)
                                err_count = write_error(err_count, max_errors, err_msg, output_file)
                    elif line.startswith('##FILTER'):
                        if not (line.startswith('##FILTER=<') and line.endswith('>')):
                            err_msg = 'line {} invalid ##FILTER header {}'.format(line_count, line)
                            err_count = write_error(err_count, max_errors, err_msg, output_file)
                        else:
                            filter = line[10:-1].split(',',1)
                            if not(len(filter) == 2
                                   and filter[0].startswith('ID=')
                                   and filter[1].startswith('Description="') and filter[1].endswith('"')):
                                err_msg = 'line {} invalid ##FILTER header {}'.format(line_count, line)
                                err_count = write_error(err_count, max_errors, err_msg, output_file)
                    elif line.startswith('##FORMAT'):
                        if not (line.startswith('##FORMAT=<') and line.endswith('>')):
                            err_msg = 'line {} invalid ##FORMAT header {}'.format(line_count, line)
                            err_count = write_error(err_count, max_errors, err_msg, output_file)
                        else:
                            format = line[10:-1].split(',',3)
                            if not(len(format) == 4
                                   and format[0].startswith('ID=')
                                   and format[1].startswith('Number=') and len(format[1])>7 and (format[1][7:] in ('.','A','G') or my.is_int(format[1][7:]))
                                   and format[2].startswith('Type=') and len(format[2])>5 and format[2][5:] in valid_types
                                   and format[3].startswith('Description="') and format[3].endswith('"')):
                                err_msg = 'line {} invalid ##FORMAT header {}'.format(line_count, line)
                                err_count = write_error(err_count, max_errors, err_msg, output_file)
                    elif line.startswith('##ALT'):
                        if not (line.startswith('##ALT=<') and line.endswith('>')):
                            err_msg = 'line {} invalid ##ALT header {}'.format(line_count, line)
                            err_count = write_error(err_count, max_errors, err_msg, output_file)
                        else:
                            alt = line[7:-1].split(',',1)
                            if not(len(alt) == 2
                                   and alt[0].startswith('ID=')
                                   and alt[0].split(':')[0] not in valid_alt_ids
                                   #TODO:VCF4.1 doesn't show Description requiring quotes but 1000 genomes uses quotes
                                   #and alt[1].startswith('Description=')):
                                   and alt[1].startswith('Description="') and alt[1].endswith('"')):
                                err_msg = 'line {} invalid ##ALT header {}'.format(line_count, line)
                                err_count = write_error(err_count, max_errors, err_msg, output_file)
                    elif line.startswith('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO'):
                        header_seen = True
                        header_cols = line.split('\t')
                        column_count = len(header_cols)
                        if column_count > 8:
                            if header_cols[8] != 'FORMAT':
                                err_msg = 'line {} no FORMAT column {}'.format(line_count, line)
                                err_count = write_error(err_count, max_errors, err_msg, output_file)
                                if column_count < 10:
                                    warn('no genotype columns')
                    else:
                        continue
                elif not header_seen:
                    err_msg = 'line {} missing header {}'.format(line_count, line)
                    err_count = write_error(err_count, max_errors, err_msg, output_file)
                else:
                    fields = line.split('\t')
                    field_count = len(fields)
                    if field_count != column_count:
                        err_msg = 'line {} has field count {} not equal to column count {} {}'.format(line_count, field_count, column_count, line)
                        err_count = write_error(err_count, max_errors, err_msg, output_file)
                    if ws_colon_re.search(fields[CHROM]):
                        err_msg = 'line {} CHROM has whitespace or colon {}'.format(line_count, line)
                        err_count = write_error(err_count, max_errors, err_msg, output_file)
                    if not my.is_int(fields[POS]):
                        err_msg = 'line {} POS {} is not numeric {}'.format(line_count, fields[POS], line)
                        err_count = write_error(err_count, max_errors, err_msg, output_file)
                    if fields[CHROM] == that_chrom:
                        if int(fields[POS]) < that_pos:
                            err_msg = 'line {} POS {} < {} {}'.format(line_count, fields[POS], that_pos, line)
                            err_count = write_error(err_count, max_errors, err_msg, output_file)
                        else:
                            that_pos = int(fields[POS])
                    else:
                        that_chrom = fields[CHROM]
                        that_pos = 0
                    if ws_re.search(fields[ID]):
                        err_msg = 'line {} ID {} has whitespace {}'.format(line_count, fields[ID], line)
                        err_count = write_error(err_count, max_errors, err_msg, output_file)
                    if not_ACGTN_re.search(fields[REF]):
                        err_msg = 'line {} REF {} has non-ACGTN {}'.format(line_count, fields[REF], line)
                        err_count = write_error(err_count, max_errors, err_msg, output_file)
                    alt = fields[ALT].split(',')
                    for a in alt:
                        #allow X in ALT to accomodate non-standard VCFs created by linkdatagen
                        if not (a == '.' or (a.startswith('<') and a.endswith('>')) or fields[INFO].find('SYTYPE=BND')!=-1 or fields[INFO].find('IMPRECISE')!=-1) and not_ACGTNX_re.search(a):
                            err_msg = 'line {} ALT {} has invalid character(s) {}'.format(line_count, a, line)
                            err_count = write_error(err_count, max_errors, err_msg, output_file)
                    if not(fields[QUAL] == '.' or my.is_number(fields[QUAL])):
                        err_msg = 'line {} QUAL {} is not numeric {}'.format(line_count, fields[QUAL], line)
                        err_count = write_error(err_count, max_errors, err_msg, output_file)
                    if fields[FILTER] == '0' or ws_semicolon_re.search(fields[ALT]):
                        err_msg = 'line {} FILTER {} is zero or has whitespace or semicolon {}'.format(line_count, fields[FILTER], line)
                        err_count = write_error(err_count, max_errors, err_msg, output_file)
                    info = fields[INFO].split(';')
                    for i in info:
                        kv = i.split('=',1)
                        if len(kv)>1 and ws_semicolon_equal_sign_re.search(kv[1]):
                            err_msg = 'line {} INFO {} has whitespace, semicolon, or equals-sign {}'.format(line_count, i, line)
                            err_count = write_error(err_count, max_errors, err_msg, output_file)
                    if field_count > 8:
                        format = fields[FORMAT].split(':')
                        if format[0] and format[0] != 'GT':
                            err_msg = 'line {} FORMAT {} does not start with GT {}'.format(line_count, fields[FORMAT], line)
                            err_count = write_error(err_count, max_errors, err_msg, output_file)
                    if field_count > 9 and format[0] and format[0] == 'GT':
                        for i in range(9,field_count-1):
                            gt = fields[i].split(':')[0]
                            m = gt_re.match(gt)
                            if not m:
                                err_msg = 'line {} invalid FORMAT {} {}'.format(line_count, fields[FORMAT], line)
                                err_count = write_error(err_count, max_errors, err_msg, output_file)
                    data_seen = True
        if header_seen and data_seen:
            output_file.write('OK')
        else:
            err_msg = 'line {} header and/or data missing {}'.format(line_count, line)
            err_count = write_error(err_count, max_errors, err_msg, output_file)
    except Exception as e:
        err_msg = 'error {} validating input file {}'.format(e, input)
        err_count = write_error(err_count, max_errors, err_msg, output_file)
    else:
        return 0 if err_count == 0 else 100
    finally:
        output_file.close()


    
def main():
   
    #command line arguments
    parser = argparse.ArgumentParser(
        description = 'validate vcf file',
        epilog = 'pypeline.validate_vcf version 1.0β1 ©2011-2012 Michael Yourshaw all rights reserved')
    parser.add_argument('--input', '-v', '-i', required=True,
        help='input vcf file to be validated')
    parser.add_argument('--output', '-o',
        help='output file(default: <input vcf>.validate; if input is valid file output contains "OK")')
    parser.add_argument('--max_errors', type=int, default=1,
        help='maximum number of errors to report (default: 1)')
    args = parser.parse_args()
    
    status = run(args.input, args.output, args.max_errors)
    
    return status


if __name__ == "__main__": sys.exit(main())
