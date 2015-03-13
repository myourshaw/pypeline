#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import stat
import argparse
from ConfigParser import SafeConfigParser #configparser in python 3
from glob import glob
from warnings import warn
import tempfile
from time import localtime, strftime
import re
from warnings import warn
import my
import job
import qseq2bam

class Qseqs2BamsError(Exception): pass

#1/12 lane
#-m /scratch0/tmp/myourshaw/gmd/metadata.txt -d /scratch0/tmp/myourshaw/gmd/ -q /scratch0/tmp/myourshaw/gmd/qseqs/HWI-ST430/243/6/HWI-ST430.243.6.[123].ATCACGA_qseq.txt.gz

#all
#-m /scratch0/tmp/myourshaw/gmd/metadata.txt -d /scratch0/tmp/myourshaw/gmd/ -q /scratch0/tmp/myourshaw/gmd/qseqs/*/*/*/*_qseq.txt.gz

def main():
    
    #command line arguments
    parser = argparse.ArgumentParser(parents=[my.default_parser()],
        description = 'From a metadata file, align multiple readgroups of qseq readgroup files with novoalign, fix mate information with picard',
        epilog = 'pypeline.qseqs2bams version 1.0β1 ©2011-2012 Michael Yourshaw all rights reserved')
    parser.add_argument('--metadata', '-m', required=True,
        help='metadata file to map machine, run, lane, barcode to adapters, library, sample and other readgroup info',)
    parser.add_argument('--qseqs', '-q', nargs='+',
        help='qseqs to align',)
    parser.add_argument('--barcode_re', '-b', default=r"\.([ACGT]+)_qseq.txt",
        help='regular expression for barcode in qseq path/filename (default("\.([ACGT]+)_qseq.txt")')
    parser.add_argument('--novoalign',
        help='path to novoalign executable (default: from config->novocraft.novoalign)')
    parser.add_argument('--index',
        help='path to novoindex-produced index (default: from config->novocraft.default_index)')
    parser.add_argument('--reference', '-r',
        help='path to reference sequence fasta file (default: from config->reference.default_reference)')
    parser.add_argument('--samtools',
        help='path to samtools executable (default: from config->samtools.samtools)')

    args = parser.parse_args()
    dirs = my.create_default_directory_structure(args.top_dir)
    config = my.get_config(args)
    
    args.novoalign = args.novoalign if args.novoalign else config.get('novocraft','novoalign')
    if not my.is_exe(args.novoalign):
        raise Qseqs2BamsError('{} is not a novoalign executable'.format(args.novoalign))
    args.index = args.index if args.index else config.get('novocraft','default_index')
    if not my.file_exists(args.index):
        raise Qseqs2BamsError('cannot find novolaign genome index file {}'.format(args.index))
    himem = bool(os.path.getsize(args.index) > 6500000000)
    args.reference = args.reference if args.reference else config.get('reference','default_reference')
    if not my.file_exists(args.reference):
        raise Qseqs2BamsError('cannot find reference sequence fasta file {}'.format(args.reference))
    args.samtools = args.samtools if args.samtools else config.get('samtools','samtools')
    if not my.is_exe(args.samtools):
        raise Qseqs2BamsError('{} is not a samtools executable'.format(args.samtools))

    my_name = '{}_{}'.format('qseqs2bams', my.localtime_squish())
    job_dir = my.unique_dir(my_name, dirs['jobs'])
    
    print 'Aligning reads in qseqs for ...'

    #get metadata file to map machine, run, lane, barcode to readgroup info
    metadata = my.get_metadata(args.metadata)
    if not metadata:
        raise Qseqs2BamsError('no records selected from metadata file {}'.format(args.metadata))

    qseqs = my.qseqs_hierarchy(args.qseqs, barcode_re=args.barcode_re)

    readgroup_bams = []
    cmds = []
    
    readgroups = {}
    
    print '\t'.join(['machine', 'run', 'lane', 'barcode', 'qseq1', 'qseq2', 'readgroup_bam'])
    
    #get reads 1 2 3 for each machine, run, lane, barcode
    #and add to readgroups
    for machine in qseqs:
        for run in qseqs[machine]:
            for lane in qseqs[machine][run]:
                for barcode in qseqs[machine][run][lane]:
                    pk = (str(machine), str(run), str(lane), str(barcode))
                    md = metadata[pk]
                    sample = md.get('sample', None)
                    if not sample:
                        raise Qseqs2BamsError('no sample specified in metadata file for ({})'.format(','.join(pk)))
                    adapter_1, adapter_2 = md.get('adapter_1', ''), md.get('adapter_2', '')
                    if not adapter_1 or not adapter_2:
                        raise Qseqs2BamsError('adapter_1 and/or adapter_2 missing in metadata file for ({})'.format(','.join(pk)))
                    #reads in qseqs[machine][run][lane][barcode] are represented as (path, read_length)
                    #first read
                    ql1 = qseqs[machine][run][lane][barcode].get(1, ())
                    #second read of unbarcoded paired end run, or barcode read of barcoded run
                    ql2 = qseqs[machine][run][lane][barcode].get(2, ())
                    #second read of barcoded paired end run
                    ql3 = qseqs[machine][run][lane][barcode].get(3, ())
                    q1 = ql1[0] if ql1 else ''
                    q2 = ql2[0] if ql2 else ''
                    q3 = ql3[0] if ql3 else ''
                    l1 = int(md['read_1_length']) if md['read_1_length'] else ql1[1] if q1l else 0
                    l2 = int(md['read_2_length']) if md['read_2_length'] else ql2[1] if ql2 else 0
                    l3 = int(md['read_3_length']) if md['read_3_length'] else ql3[1] if ql3 else 0
                    qseq1 = ''
                    qseq2 = ''
                    read_length_1 = 0
                    read_length_2 = 0
                    if len([r for r in [q1, q2, q3] if r]) == 0:
                        warn ("Can't find reads 1, 2 , or 3 for ({}).".format(','.join(pk)))
                    if len([r for r in [q1, q2, q3] if r]) == 1:
                        #align sole read regardless of read number
                        qseq1 = q1 if q1 else q2 if q2 else q3
                        read_length_1 = l1 if q1 else l2 if q2 else l3
                        if qseq1 != q1:
                            warn("Aligning {} as the only read file for ({}).".format(qseq1, ','.join(pk)))
                    if q1:
                        qseq1 = q1
                        read_length_1 = l1
                        if q3:
                            #barcoded paired end run
                            qseq2 = q3
                            read_length_2 = l3
                        elif q2:
                            #unbarcoded paired end run
                            #or barcoded single end run
                            if not barcode:
                                qseq2 = q2
                                read_length_2 = l2
                    elif q3:
                        #was barcoded paired end but first read not submitted
                        qseq1 = q3
                        read_length_1 = l3
                        warn("{} seems to be the second read of a barcoded paired end run but the first read was not submitted for ({}).".format(qseq1, ','.join(pk)))
                    if read_length_1 != read_length_2 or read_length_1 < 36 or read_length_2 < 36:
                        warn('Read lengths {} and/or {} are unequal or less than 36 bases for '.format(read_length_1, read_length_2, ','.join(pk)))
                    read_length = max((read_length_1, read_length_2))
                    library = md.get('library', sample)
                    library_protocol = md.get('library_protocol', '')
                    sequencing_library = md.get('sequencing_library', '')
                    sequencing_center = md.get('sequencing_center', 'UCLA')
                    run_date = md.get('run_date', None)
                    if not run_date:
                        run_date = strftime('%Y-%m-%d', localtime(os.path.getctime(qseq1)))
                    predicted_median_insert_size = md.get('predicted_median_insert_size', '')
                    platform = md.get('platform', 'ILLUMINA')
                    instrument_model = md.get('instrument_model', '')
                    flowcell = md.get('flowcell', '')
                    platform_unit = '{}.{}{}.{}{}'.format(machine, run, '.' + flowcell if flowcell else '', lane, '.' + barcode if barcode else '')
                    run_folder = md.get('run_folder', '')
                    readgroup_id = '{}.{}.{}{}'.format(machine, run, lane, '.' + barcode if barcode else '')
                    readgroup_description = 'sample={};library={};library_protocol={};sequencing_library={};sequencing_center={};platform={};instrument_model={};platform_unit={};run_date={};run_folder={};machine={};run={};flowcell={};lane={};barcode={};read_length={};predicted_median_insert_size={};adapter_1={};adapter_2={}'\
                        .format(sample, library, library_protocol, sequencing_library, sequencing_center, platform, instrument_model, platform_unit, run_date, run_folder, machine, run, flowcell, lane, barcode, read_length, predicted_median_insert_size, adapter_1, adapter_2)
                    addtional_readgroup_description = md.get('addtional_readgroup_description', None)
                    if addtional_readgroup_description:
                        readgroup_description = readgroup_description + ';' + addtional_readgroup_description
                    RG = '"@RG\\tID:{}\\tCN:{}\\tDT:{}\\tLB:{}\\tPI:{}\\tPL:{}\\tPU:{}\\tSM:{}\\tDS:{}"'.\
                        format(readgroup_id, sequencing_center, run_date, library, predicted_median_insert_size, platform, platform_unit, sample, readgroup_description)

                    output_name = '{}.{}.{}'.format(sample, library, readgroup_id)
                    readgroup_bam = os.path.join(dirs['readgroup_bams'], '{}.novoalign.bam'.format(output_name))
                    readgroup_bams += [readgroup_bam]

                    print '\t'.join([str(machine), str(run), str(lane), str(barcode), str(qseq1), str(qseq2), str(readgroup_bam)])

                    readgroups[pk] = {'qseqs': (qseq1, qseq2), 'readgroup_bam': readgroup_bam, 'RG': RG, 'read_length': read_length, 'adapters': (adapter_1, adapter_2)}

    print 'Queuing {} jobs. Job_info, STDOUT, and STDERR will be in {}. ...'.format(my_name, job_dir)    
    for pk in readgroups:
        rg = readgroups[pk]
        cmds += ['python {} --novoalign {} --index {} --reference {} --samtools {} --parent_job_dir {} --qseqs {} --readgroup_bam {} --RG {} --read_length {} --adapters {}'
            .format(qseq2bam.__file__, args.novoalign, args.index, args.reference, args.samtools, job_dir, ' '.join(rg['qseqs']), rg['readgroup_bam'], rg['RG'], rg['read_length'], ' '.join(rg['adapters']))]
        
    job_name = my_name
    job = my.run_job(cmds, job_name, job_dir, synchronous=False)
    print 'Job {}, jobId {}, queued. Job_info, STDOUT, and STDERR will be in {}.'.format(job_name, job.jobId, job_dir)


if __name__ == "__main__": sys.exit(main())
