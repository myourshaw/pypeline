#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import argparse
from ConfigParser import SafeConfigParser
import re
from glob import glob, iglob
from time import localtime, strftime
import shutil
import my
import job
import qseq_metrics
import picard
import verify_sam_file_validation
import bam_metrics
import merge_validate_qseqs

class Fastq2BamError(Exception): pass

#-d /scratch1/tmp/myourshaw/gmdjen -m /scratch1/tmp/myourshaw/gmdjen/gmdjen_metadata.txt --fastqs /data/storage-1-02/archive/myourshaw/HiSeq2k.Pathology/120417_SN973_0083_AC0P8TACXX/Data/Intensities/BaseCalls/Unaligned/Project_GMD/Sample_GMD73/GMD73_CGATGT_L004_R[12]_001.fastq.gz

def get_lane_info_from_fastqs(fastqs, readgroups):
    
    #lane_info (the return value) is a collection of Lane objects
    lane_info = my.Lanes()
    
    fastq_info = my.fastqs_hierarchy(fastqs)
    for machine in sorted(fastq_info):
        for run in sorted(fastq_info[machine]):
            for lane in sorted(fastq_info[machine][run]):
                        
                #l is a Lane object in to hold data about this lane
                lane_info[(machine,run,lane)] = my.Lane(machine = machine, run = run, lane = lane)
                l = lane_info[(machine,run,lane)]
                l.readgroups = [r for r in readgroups if r[0]==machine and r[1]==run and r[2]==str(lane)]
                
                #verify that this machine, run, lane is in metadata
                if not l.readgroups or len(l.readgroups) == 0:
                    raise Fastq2BamError('No metadata for machine {} run {} lane {}'.format(machine, run, lane))
                    
                #readgroup metadata tells us which (if any) barcodes are in this lane
                l.barcodes = sorted(set([r[3] for r in l.readgroups]))
                
                max_read_lengths = [0,0,0]
                for read in sorted(qseq_info[machine][run][lane]):
                    l.reads[read] = my.LaneRead(read=read,
                        lane_qseq = os.path.join(l.lane_qseqs_dir,'{}.{}.{}.{}_qseq.txt.gz'.format(machine, run, lane, read))
                    )
                    
                    #r is a Read object; there can be 1-3 read qseqs per tile
                    r = l.reads[read]
                    #many tile qseqs that will be merged into one lane qseq per read
                    compression = set()
                    for tile in sorted(qseq_info[machine][run][lane][read]):
                        r.files_to_merge += [qseq_info[machine][run][lane][read][tile]['path']]
                        max_read_lengths[read-1] = max(max_read_lengths[read-1], qseq_info[machine][run][lane][read][tile]['read_length'])
                        compression.add(qseq_info[machine][run][lane][read][tile]['compressed'])
                    tile_qseq_created_date = strftime('%Y-%m-%d', localtime(os.path.getctime(r.files_to_merge[0])))
                    
                    #verify that all tile qseqs are either compressed or uncompressed
                    if True in compression and False in compression:
                        raise MergeTileQseqsError("Compressed and uncompressed qseqs cannot be mixed")
                    r.tile_files_compressed = compression.pop()
                    
                    #command to merge tile qseqs for this read
                    #r.merge_cmd = 'gzip -cd {} | cat - | gzip > {}'.format(' '.join(r.files_to_merge), r.lane_qseq) \
                    #    if r.tile_files_compressed else 'cat {} | gzip > {}'.format(' '.join(r.files_to_merge), r.lane_qseq)
                    r.merge_cmd = '{} {} --qseqs {} --output {}'.format(python, merge_validate_qseqs.__file__, ' '.join(r.files_to_merge), r.lane_qseq, '--overwrite' if overwrite_merged_tile_qseqs else '')
                    
                #identify data reads and lengths thereof
                file1,file2 = '',''
                length1,length2 = 0,0
                
                #the readgroup qseq for this read is the entire lane qseq if the lane isn't barcoded
                if l.barcodes == ['']:
                    barcode = ''
                    this_rg = readgroups[(machine,run,str(lane),barcode)]
                    read1 = l.reads.get(1,None)
                    read2 = l.reads.get(2,None)
                    read3 = l.reads.get(3,None)
                    if read1:
                        file1 = read1.lane_qseq
                        length1 = max_read_lengths[0]
                        if read3:
                            file2 = read3.lane_qseq
                            length2 = max_read_lengths[2]
                        #is read2 a barcode?
                        elif read2 and max_read_lengths[1] > 32:
                            file2 = read2.lane_qseq
                            length2 = max_read_lengths[1]
                    elif read3:
                        file1 = read3.lane_qseq
                        length1 = max_read_lengths[2]
                    elif read2:
                        file1 = read2.lane_qseq
                        length1 = max_read_lengths[1]
                    this_rg.readgroup_qseqs += [file1]
                    this_rg.readgroup_qseqs_read_lengths += [length1]
                    if file2:
                        this_rg.readgroup_qseqs += [file2]
                        this_rg.readgroup_qseqs_read_lengths += [length2]
                    this_rg.tile_qseq_created_date = tile_qseq_created_date
                    
                #directories, files and commands to demultiplex a barcoded lane
                else: #if l.barcodes != ['']:
                    #check barcode_tag_list in metadata
                    barcode_tag_list = set()
                    for barcode in l.barcodes:
                        this_rg = readgroups[(machine,run,str(lane),barcode)]
                        
                    l.lane_qseqs_demux_dir = os.path.join(l.lane_qseqs_dir, 'demux')
                    my.makedir(l.lane_qseqs_demux_dir)
                    read1 = l.reads.get(1,None)
                    file1 = read1.lane_qseq if read1 else ''
                    if file1:
                        length1 = max_read_lengths[0]
                    read3 = l.reads.get(3,'')
                    file2 = read3.lane_qseq if read3 else ''
                    if file2:
                        length2 = max_read_lengths[2]
                    read2 = l.reads.get(2,'')
                    qseqtagfile = read2.lane_qseq if read2 else ''
                    if not qseqtagfile:
                        raise Fastq2BamError('qseq tag file (read 2) missing for machine {}, run {}, lane {}'
                            .format(machine, run, lane))
                    l.demultiplex_command = '{0} -b {1} -d {2} -F QSEQ -f {3} {4} -i {5} --ILQ_SKIP --QSEQ_OUT --GZIP;'\
                        .format(novobarcode, barcodefile, l.lane_qseqs_demux_dir, file1, file2, qseqtagfile)
                    for barcode in l.barcodes:
                        this_rg = readgroups[(machine,run,str(lane),barcode)]
                        barcode_dir = os.path.join(l.lane_qseqs_demux_dir, barcode)
                        file1_demux = os.path.join(l.lane_qseqs_dir, my.qseq_strip(os.path.basename(file1))+'.'+barcode+'_qseq.txt.gz') if file1 else None
                        if file1_demux:
                            this_rg.readgroup_qseqs += [file1_demux]
                            this_rg.readgroup_qseqs_read_lengths += [length1]
                        file2_demux = os.path.join(l.lane_qseqs_dir, my.qseq_strip(os.path.basename(file2))+'.'+barcode+'_qseq.txt.gz') if file2 else None
                        if file2_demux:
                            this_rg.readgroup_qseqs += [file2_demux]
                            this_rg.readgroup_qseqs_read_lengths += [length2]
                        qseqtagfile_demux = os.path.join(l.lane_qseqs_dir, my.qseq_strip(os.path.basename(qseqtagfile))+'.'+barcode+'_qseq.txt.gz')
                        this_rg.tile_qseq_created_date = tile_qseq_created_date
                        l.lane_qseqs_demux_move_commands += \
                            ['for qseq in {}/*_qseq.txt.gz; do base=`basename $qseq`; mv -a $qseq {}/${{base%_qseq.txt.gz}}.{}_qseq.txt.gz; done'.format(barcode_dir, l.lane_qseqs_dir, barcode)]
                    #l.lane_qseqs_remove_intermediates_command = \
                    #    'rm -rf {} {} {} {}'.format(l.lane_qseqs_demux_dir, file1, file2, qseqtagfile)
                    #keep multiplexed merged qseq files, so we can auto skip merging on rerun
                    #TODO this wouldn't be necessary if we were smarter about autoskiping novobarcode
                    l.lane_qseqs_remove_intermediates_command = \
                        'rm -rf {}'.format(l.lane_qseqs_demux_dir)
    return lane_info

def configure_readgroups(dirs, readgroups, novoalign, index, samtools):

    for readgroup in [r for r in readgroups if readgroups[r].readgroup_qseqs]:
        rg = readgroups[readgroup]
        if rg.readgroup_qseqs and len(rg.readgroup_qseqs) > 0:
            
            #output files that will be created
            rg.readgroup_bam = os.path.join(
                dirs['readgroup_bams'], '{}.{}.{}{}.novoalign.bam'.format(rg.machine, rg.run, rg.lane, '.'+rg.barcode if rg.barcode else ''))
            rg.novoalign_bam_validate = rg.readgroup_bam + '.validate'
            rg.fixmate_bam = my.bam_strip(rg.readgroup_bam)+'.fixmate.bam'
            rg.fixmate_bai = my.bam_strip(rg.fixmate_bam) + '.bai'
            rg.fixmate_bam_validate = rg.fixmate_bam + '.validate'
            rg.library_bam = os.path.join(dirs['library_bams'], '{}.library.bam'.format(rg.library))
            rg.library_markdup_bam = os.path.join(dirs['library_bams'], '{}.library.markdup.bam'.format(rg.library))
            rg.sample_bam = os.path.join(dirs['sample_bams'], '{}.sample.bam'.format(rg.sample))
            rg.study_bam = os.path.join(dirs['study_bams'], '{}.study.bam'.format(rg.study))
            
            #determine read length from metadata if possible, otherwise from first row of qseq file
            #to calculate novoalign parameter for min quality bases
            read_lengths = [r for r in rg.readgroup_qseqs_read_lengths if my.is_number(r) and r != 0]
            read_length = min(read_lengths) if read_lengths else 50
            min_quality_bases = int(round(read_length/2.0)) #novoalign default is ~20; this yields 25 for 50 base reads and 50 for 100 base reads
           
            #compose RG string for bam header
            if not rg.run_date:
                rg.run_date = rg.tile_qseq_created_date
            if not rg.platform:
                rg.platform = 'ILLUMINA'
            if not rg.instrument_model:
                rg.instrument_model = 'Illumina_HiSeq_2000'
            if not rg.sequencing_center:
                rg.sequencing_center = 'UCLA'
            if not my.is_number(rg.predicted_median_insert_size):
                rg.predicted_median_insert_size = 200
            rg.platform_unit = '{}.{}{}.{}{}'.format(rg.machine, rg.run, '.' + rg.flowcell if rg.flowcell else '', rg.lane, '.' + rg.barcode if rg.barcode else '')
            rg.readgroup_id = '{}.{}.{}{}'.format(rg.machine, rg.run, rg.lane, '.' + rg.barcode if rg.barcode else '')
            rg.readgroup_description = 'sample={};tumor_status={};library={};library_protocol={};sequencing_library={};sequencing_center={};platform={};instrument_model={};platform_unit={};run_date={};run_folder={};machine={};run={};flowcell={};lane={};barcode={};read_length={};predicted_median_insert_size={};adapters={};'.format(
                rg.sample, rg.tumor_status, rg.library, rg.library_protocol, rg.sequencing_library, rg.sequencing_center, rg.platform, rg.instrument_model, rg.platform_unit, rg.run_date, rg.run_folder, rg.machine, rg.run, rg.flowcell, rg.lane, rg.barcode, ','.join([str(r) for r in read_lengths]), rg.predicted_median_insert_size, ','.join((rg.adapter_1, rg.adapter_2)))
            if rg.additional_readgroup_description:
                rg.readgroup_description = rg.readgroup_description + rg.addtional_readgroup_description + '' if rg.addtional_readgroup_description.endswith(';') else ';'
            rg.RG = """@RG\\tID:{}\\tCN:{}\\tDT:{}\\tLB:{}\\tPI:{}\\tPL:{}\\tPU:{}\\tSM:{}\\tDS:{}""".format(
                rg.readgroup_id, rg.sequencing_center, rg.run_date, rg.library, rg.predicted_median_insert_size, rg.platform, rg.platform_unit, rg.sample, rg.readgroup_description)
                
            #novoalign command
            rg.novoalign_cmd = "{} -k -o SAM '{}' -d {} -a {} {} -F QSEQ -l {} -H --hdrhd 1 -f {} --mmapoff | {} view -S1  -  > {};".format(
                novoalign, rg.RG, index, rg.adapter_1, rg.adapter_2, min_quality_bases, ' '.join(rg.readgroup_qseqs), samtools, rg.readgroup_bam)


def main():

    #command line arguments
    parser = argparse.ArgumentParser(parents=[my.default_parser()],
        description = 'align fastq files to readgroup bams and merge into library and sample bams',
        epilog = 'pypeline.fastq2bam version 1.0β1 ©2011-2012 Michael Yourshaw all rights reserved')
    parser.add_argument('--metadata', '-m', required=True,
        help='metadata file to map machine, run, lane, barcode to adapters, library, sample and other readgroup info',)
    parser.add_argument('--fastqs', '-q', nargs='+', required=True,
        help='list of fastq files')
    parser.add_argument('--python',
        help='python executable (default: from config->DEFAULT.python)')
    parser.add_argument('--reference', '-r',
        help='path to reference sequence fasta file (default: from config->reference.default_reference)')
    parser.add_argument('--novoalign',
        help='path to novoalign executable (default: from config->novocraft.novoalign)')
    parser.add_argument('--index',
        help='path to novoindex-produced index (default: from config->novocraft.default_index)')
    parser.add_argument('--samtools',
        help='path to samtools executable (default: from config->samtools.samtools)')
    parser.add_argument('--overwrite_all_output_files', action='store_true', default=False,
        help='overwriting any existing valid files')
    parser.add_argument('--overwrite_merged_tile_qseqs', action='store_true', default=False,
        help='redo merge tile qseqs, overwriting valid merged qseq files')
    parser.add_argument('--overwrite_novoalign_fixmate_bams', action='store_true', default=False,
        help='redo novoalign and FixMateInformation, overwriting existing valid novolaign and fixmate bam files')
    parser.add_argument('--overwrite_merge_readgroup_bams', action='store_true', default=False,
        help='redo merge readgroup bams, overwriting existing valid library bam files')
    parser.add_argument('--overwrite_markdup_library_bams', action='store_true', default=False,
        help='redo MarkDuplicates, overwriting existing valid markdup library bam files')
    parser.add_argument('--overwrite_merge_library_bams', action='store_true', default=False,
        help='redo merge library bams, overwriting existing valid sample bam files')
    parser.add_argument('--NoQseq2SampleBams', action='store_true', default=False,
        help='omit all Qseq2SampleBams steps; just do GATK steps')
    parser.add_argument('--NoProcessQseqs', action='store_true', default=False,
        help='omit ProcessQseqs steps of merging tile and nemultiplexing')
    parser.add_argument('--NoNovoalign', action='store_true', default=False,
        help='omit Novoalign step')
    parser.add_argument('--NoFixMateInformation', action='store_true', default=False,
        help='omit FixMateInformation step')
    parser.add_argument('--NoMergeReadgroups2Libraries', action='store_true', default=False,
        help='omit MergeReadgroups2Libraries step')
    parser.add_argument('--NoMergeLibraries2Samples', action='store_true', default=False,
        help='omit MergeLibraries2Samples step')
    parser.add_argument('--NoGATK', action='store_true', default=False,
        help='omit all GATK steps')
    parser.add_argument('--NoQualityScoreRecalibration', action='store_true', default=False,
        help='omit QualityScoreRecalibration step')
    parser.add_argument('--NoUnifiedGenotyper', action='store_true', default=False,
        help='omit UnifiedGenotyper step')
    parser.add_argument('--NoVariantRecalibration', action='store_true', default=False,
        help='omit VariantRecalibration step')
    parser.add_argument('--NoVariantAnnotation', action='store_true', default=False,
        help='omit VariantAnnotation step')
    #parser.add_argument('--No', action='store_true', default=False,
    #    help='omit  step')
    parser.add_argument('--NoTileQseqMetrics', action='store_true', default=False,
        help='omit TileQseqMetrics')
    parser.add_argument('--NoReadgroupQseqMetrics', action='store_true', default=False,
        help='omit ReadgroupQseqMetrics')
    parser.add_argument('--NoReadgroupBamMetrics', action='store_true', default=False,
        help='omit ReadgroupBamMetrics')
    parser.add_argument('--NoLibraryBamMetrics', action='store_true', default=False,
        help='omit LibraryBamMetrics')
    parser.add_argument('--NoLibraryMarkDupBamMetrics', action='store_true', default=False,
        help='omit LibraryMarkDupBamMetrics')
    parser.add_argument('--NoSampleBamMetrics', action='store_true', default=False,
        help='omit SampleBamMetrics')
    parser.add_argument('--NoVariantMetrics', action='store_true', default=False,
        help='omit VariantMetrics')

    args = parser.parse_args()

    dirs = my.create_default_directory_structure(args.top_dir)
    config = my.get_config(args)
    if not args.python:
        args.python = config.get('DEFAULT','python')
        
    regex = config.get('picard', 'READ_NAME_REGEX')

    my_name = 'fastq2bam'
    my_job_dir =  my.unique_dir('fastq2bam_'+my.localtime_squish(), dirs['jobs'])
    
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
    
    #get readgroup info from metadata file
    print 'Getting metadata ...'
    readgroups = my.get_readgroups_from_metadata(args.metadata)
 
    #populate collection of Lane objects with directories, files, commands
    print 'Getting information from fastq files ...'
    lane_info = get_lane_info_from_fastqs(fastqs=args.fastqs, readgroups=readgroups)
 
    #command for tile qseq metrics
    tile_qseq_metrics_job_dir = my.unique_dir('tile_qseq_metrics', my_job_dir)
    all_tile_qseqs = set(my.flatten([glob(q) for q in args.qseqs]))
    tile_qseq_metrics_cmd = "{} {} -d {} -j {} -q {} -o {} -b '{}' --tiles".format(
       args.python, qseq_metrics.__file__, dirs['top'], tile_qseq_metrics_job_dir, ' '.join(all_tile_qseqs), os.path.join(dirs['qseq_metrics'],'tile_qseq_metrics_'+my.localtime_squish()+'.txt'), args.barcode_re)

    #command for readgroup qseq metrics
    rgs_qseq_metrics_job_dir = my.unique_dir('readgroup_qseq_metrics', my_job_dir)
    all_rg_qseqs = my.flatten([readgroups[r].readgroup_qseqs for r in readgroups if readgroups[r].readgroup_qseqs])
    rgs_qseq_metrics_cmd = "{} {} -d {} -j {} -q {} -o {} -b '{}'".format(
        args.python, qseq_metrics.__file__, dirs['top'], rgs_qseq_metrics_job_dir, ' '.join(all_rg_qseqs), os.path.join(dirs['qseq_metrics'],'readgroup_qseq_metrics_'+my.localtime_squish()+'.txt'), args.barcode_re)

    print "Configuring readgroup information ..."
    configure_readgroups(dirs, readgroups, args.novoalign, args.index, args.samtools)
    
    #verify that all barcodes are in the novobarcode tag file
    all_barcodes = set(my.flatten([lane_info[l].barcodes for l in lane_info if lane_info[l].barcodes != ['']]))
    if all_barcodes:
        print 'Checking barcodes ...'
        if not all_barcodes.issubset(tag_file_barcodes):
            raise Fastq2BamError('Qseq barcodes {} are not in tag file {} barcode set {}'
                .format(all_barcodes-tag_file_barcodes, barcodefile, tag_file_barcodes))

    #create bam files from qseq files
    if not args.NoQseq2SampleBams:
        hold_lane_jids = []
        if not args.NoProcessQseqs:
            
            #qseq metrics by tile
            if not args.NoTileQseqMetrics:
                print 'Submitting jobs to calculate metrics of tile qseq files ...'
                job_dir = tile_qseq_metrics_job_dir
                cmd = tile_qseq_metrics_cmd
                job_name = my_name+'_tile_qseq_metrics'
                job = my.run_job(cmd, job_name, job_dir)
                tile_metrics_jobid = job.jobId
                
            print 'Submitting jobs to merge tile qseqs into lane qseqs ...'
            job_dir = my.unique_dir('merge_tile_qseqs', my_job_dir)
            hold_novoalign_jids = []
            for k in sorted(lane_info):
                key = '{}_{}_{}'.format(*k)
                l = lane_info[k]
                
                #merge tile qseqs into lane qseqs
                for read in l.reads:
                    cmd = l.reads[read].merge_cmd
                    #print cmd
                    job_name = my_name+'_merge_tile_qseqs_'+key
                    job = my.run_job(cmd, job_name, job_dir)
                    merge_tile_qseqs_jobid = job.jobId
                    hold_lane_jids += [job.jobId]
                    
                #demultiplex barcoded lanes
                if l.barcodes != ['']:
                    #print 'Submitting jobs to demultiplex barcoded lanes ...'
                    cmd = l.demultiplex_command
                    #print cmd
                    job_name = my_name+'_demux_lane_qseqs_'+key
                    job = my.run_job(cmd, job_name, job_dir, hold_jid=merge_tile_qseqs_jobid)
                    demux_lane_qseqs_jobid = job.jobId
                    hold_lane_jids += [job.jobId]
                    
                    #move demultiplexed files into lane directory
                    cmds = l.lane_qseqs_demux_move_commands
                    #print '\n'.join(cmds)
                    job_name = my_name+'_compress_demux_lane_qseqs_'+key
                    job = my.run_job(cmds, job_name, job_dir, hold_jid=demux_lane_qseqs_jobid)
                    compress_demux_lane_qseqs_jobid = job.jobId
                    hold_lane_jids += [job.jobId]
                    
                    #delete barcode directories and multiplexed lane qseq file
                    cmd = l.lane_qseqs_remove_intermediates_command
                    #print cmd
                    job_name = my_name+'_remove_demux_intermediates_'+key
                    job = my.run_job(cmd, job_name, job_dir, hold_jid=compress_demux_lane_qseqs_jobid)
                    remove_demux_intermediates_jobid = job.jobId
                    hold_lane_jids += [job.jobId]
                    
            if not args.NoReadgroupQseqMetrics:
                #qseq metrics by readgroup
                print 'Submitting jobs to calculate metrics of readgroup qseq files ...'
                job_dir = rgs_qseq_metrics_job_dir
                cmd = rgs_qseq_metrics_cmd
                #print cmd
                job_name = my_name+'_readgroup_qseq_metrics'
                job = my.run_job(cmd, job_name, job_dir, hold_jid=hold_lane_jids)
                qseq_metrics_jobid = job.jobId
                
        if not hold_lane_jids:
            hold_lane_jids = None
            
        #qseqs to readgroup bams
        hold_rg_jids = []
        novoalign_jobid = None
        job_dir = my.unique_dir('novoalign', my_job_dir)
        rgs = [readgroups[r] for r in readgroups if readgroups[r].novoalign_cmd]
        print 'Submitting jobs to align readgroup qseqs > readgroup bams and fix mate information, '
        for rg in rgs:
            if args.overwrite_all_output_files or args.overwrite_novoalign_fixmate_bams or not my.check_files([rg.fixmate_bam], [(rg.fixmate_bam_validate,'No errors found\n')], [(rg.fixmate_bam,rg.fixmate_bai), (rg.fixmate_bam,rg.fixmate_bam_validate)]):
                if not args.NoNovoalign:
                    #novoalign by readgroup
                    cmd = rg.novoalign_cmd
                    job_name = my_name+'_novoalign_'+rg.readgroup_id
                    #print cmd
                    job = my.run_job(cmd, job_name, job_dir, hold_jid=hold_lane_jids, memory = '24G' if himem else '7G')
                    novoalign_jobid = job.jobId
                    hold_rg_jids += [job.jobId]
                    
                #fix mate information with picard, overwriting the bam file created by novoalign
                #this step sorts the bam file, creates an index,
                #and ensures that all mate-pair information is in sync between each read and its mate pair.
                if not args.NoFixMateInformation:
                    #print 'Submitting jobs to fix mate information ...'
                    job_name = my_name+'_FixMateInformation_'+rg.readgroup_id
                    picard_args = {'CREATE_INDEX': 'true', 'CREATE_MD5_FILE': 'true', 'SORT_ORDER': 'coordinate', 'INPUT': rg.readgroup_bam, 'OUTPUT': rg.fixmate_bam}
                    job = picard.run(tool='FixMateInformation', picard_args=picard_args, job_dir=job_dir, hold_jid=novoalign_jobid)
                    fixmateinformation_jobid = job.jobId
                    hold_rg_jids += [job.jobId]
                    
                    #validate bam file with picard
                    job_name = my_name+'_ValidateFixMateBamFile_'+rg.readgroup_id
                    picard_args = {'INPUT': rg.fixmate_bam, 'OUTPUT': rg.fixmate_bam_validate}
                    job = picard.run(tool='ValidateSamFile', picard_args=picard_args, job_dir=job_dir, reference=args.reference, hold_jid=fixmateinformation_jobid)
                    validatefixmatebamfile_jobid = job.jobId
                    hold_rg_jids += [job.jobId]
                    
                    #verify validation
                    job_name = my_name+'_VerifyValidateFixMateBamFile_'+rg.readgroup_id
                    cmd = '{} {} -v {}'.format(args.python, verify_sam_file_validation.__file__, rg.fixmate_bam_validate)
                    job = my.run_job(cmd, job_name, job_dir, hold_jid=validatefixmatebamfile_jobid, email=args.email)
                    verifyvalidatefixmatebamfile_jobid = job.jobId
                    hold_rg_jids += [job.jobId]
                    
        if not hold_rg_jids:
            hold_rg_jids = None
            
        #readgroup bam metrics
        if not args.NoReadgroupBamMetrics:
            print 'Submitting jobs to calculate metrics of readgroup bam files ...'
            readgroup_bams = [rg.fixmate_bam for rg in rgs]
            readgroup_bam_metrics_job_dir = my.unique_dir('readgroup_bam_metrics', my_job_dir)
            job_dir = readgroup_bam_metrics_job_dir
            cmd = '{} {} -d {} --output_prefix {} -i {} -o {}'.format(
                args.python, bam_metrics.__file__, dirs['top'], 'readgroup_metrics/',' '.join(readgroup_bams) ,dirs['readgroup_bam_metrics'])
            job_name = my_name+'_readgroup_bam_metrics'
            job = my.run_job(cmd, job_name, job_dir, hold_jid=hold_rg_jids)
            readgroup_bam_metrics_jobid = job.jobId
            
        merge_rmdup_job_dir = my.unique_dir('merge_rmdup', my_job_dir)
        job_dir = merge_rmdup_job_dir
        hold_mergereadgroups2library_jids = []
        library_readgroup_bams = {b: sorted(set(my.flatten([rg.fixmate_bam for rg in rgs if rg.library_bam == b]))) for b in set(my.flatten([rg.library_bam for rg in rgs]))}
        library_markdup_bams = {b: set(my.flatten([rg.library_markdup_bam for rg in rgs if rg.library_bam == b])).pop() for b in set(my.flatten([rg.library_bam for rg in rgs]))}
        sample_library_bams = {b: sorted(set(my.flatten([rg.library_markdup_bam for rg in rgs if rg.sample_bam == b]))) for b in set(my.flatten([rg.sample_bam for rg in rgs]))}
        study_sample_bams = {b: sorted(set(my.flatten([rg.sample_bam for rg in rgs if rg.study_bam == b]))) for b in set(my.flatten([rg.study_bam for rg in rgs]))}
        
        #merge readgroups by library
        hold_merge_jobids = []
        
        if not args.NoMergeReadgroups2Libraries:
            print 'Submitting jobs to merge readgroup bam files by library and mark duplicates ...'
            for library_bam in library_readgroup_bams:
                
                library_bai = my.bam_strip(library_bam)+'.bai'
                library_validate = library_bam+'.validate'
                bams2merge = library_readgroup_bams[library_bam]
                
                if len(bams2merge) > 0:
                    if args.overwrite_all_output_files or args.overwrite_merge_readgroup_bams or not my.check_files(None, [(library_validate,'No errors found\n')], [(library_bam,library_bai), (library_bam,library_validate)]):
                        if my.file_exists(library_bam):
                            os.remove(library_bam)
                        if my.file_exists(library_bai):
                            os.remove(library_bai)
                        if my.file_exists(library_validate):
                            os.remove(library_validate)
                            
                        job_name = my_name+'_MergeReadgroups_'+os.path.basename(library_bam)
                        
                        if len(bams2merge) == 1:
                            rg_bam = bams2merge[0]
                            rg_bai = my.bam_strip(rg_bam)+'.bai'
                            rg_validate = rg_bam+'.validate'
                            cmd = 'mv {0} {1}; ln -fs {1} {0}; mv {2} {3}; ln -fs {3} {2}; mv {4} {5}; ln -fs {5} {4};'.format(
                                rg_bam, library_bam, rg_bai, library_bai, rg_validate, library_validate)
                            job = my.run_job(cmd, job_name, job_dir, hold_jid=hold_rg_jids)
                            mergereadgroups2library_jobid = job.jobId
                            hold_merge_jobids += [job.jobId]
                            hold_mergereadgroups2library_jids += [job.jobId]
                        else:
                            multi_valued_arg_strings = 'INPUT='+' INPUT='.join(bams2merge)
                            picard_args = {'CREATE_INDEX': 'true', 'CREATE_MD5_FILE': 'true', 'SORT_ORDER': 'coordinate', 'ASSUME_SORTED': 'false', 'MERGE_SEQUENCE_DICTIONARIES': 'false', 'USE_THREADING': 'true', 'OUTPUT': library_bam}
                            job = picard.run(tool='MergeSamFiles', picard_args=picard_args, multi_valued_arg_strings=multi_valued_arg_strings, job_dir=job_dir, hold_jid=hold_rg_jids)
                            mergereadgroups2library_jobid = job.jobId
                            hold_merge_jobids += [job.jobId]
                            hold_mergereadgroups2library_jids += [job.jobId]
                            
                            #validate bam file with picard
                            job_name = my_name+'_ValidateLibraryBamFile_'+os.path.basename(library_bam)
                            picard_args = {'INPUT': library_bam, 'OUTPUT': library_validate}
                            job = picard.run(tool='ValidateSamFile', picard_args=picard_args, job_dir=job_dir, reference=args.reference, hold_jid=mergereadgroups2library_jobid)
                            validatelibrarybamfile_jobid = job.jobId
                            hold_merge_jobids += [job.jobId]
                            hold_mergereadgroups2library_jids += [job.jobId]
                            
                            #verify validation
                            job_name = my_name+'_VerifyValidateLibraryBamFile_'+os.path.basename(library_bam)
                            cmd = '{} {} -v {}'.format(args.python, verify_sam_file_validation.__file__, library_validate)
                            job = my.run_job(cmd, job_name, job_dir, hold_jid=validatelibrarybamfile_jobid, email=args.email)
                            verifyvalidatelibrarybamfile_jobid = job.jobId
                            hold_merge_jobids += [job.jobId]
                            hold_mergereadgroups2library_jids += [job.jobId]
                            
                    #MarkDuplicates
                    library_markdup_bam = library_markdup_bams[library_bam]
                    library_markdup_bai = my.bam_strip(library_markdup_bam)+'.bai'
                    library_markdup_validate = library_markdup_bam+'.validate'
                    metrics_file = os.path.join(dirs['library_bam_metrics'], '{}.markdup.metrics'.format(os.path.basename(library_markdup_bam)))
                    if args.overwrite_all_output_files or args.overwrite_markdup_library_bams or not my.check_files(None, [(library_markdup_validate,'No errors found\n')], [(library_markdup_bam,library_markdup_bai), (library_markdup_bam,library_markdup_validate)]):
                        job_name = my_name+'_MarkDuplicates_'+os.path.basename(library_markdup_bam)
                        picard_args = {'CREATE_INDEX': 'true', 'CREATE_MD5_FILE': 'true', 'REMOVE_DUPLICATES': 'false', 'MAX_FILE_HANDLES_FOR_READ_ENDS_MAP': 8000, 'READ_NAME_REGEX': '"'+regex+'"', 'OPTICAL_DUPLICATE_PIXEL_DISTANCE': 100, 'METRICS_FILE': metrics_file, 'INPUT': library_bam, 'OUTPUT': library_markdup_bam}
                        job = picard.run(tool='MarkDuplicates', picard_args=picard_args, job_dir=job_dir, hold_jid=hold_merge_jobids)
                        markdup_jobid = job.jobId
                        hold_mergereadgroups2library_jids += [job.jobId]
                        
                        #validate bam file with picard
                        job_name = my_name+'_ValidateLibraryMarkDupBamFile_'+os.path.basename(library_markdup_bam)
                        picard_args = {'INPUT': library_markdup_bam, 'OUTPUT': library_markdup_validate}
                        job = picard.run(tool='ValidateSamFile', picard_args=picard_args, job_dir=job_dir, reference=args.reference, hold_jid=markdup_jobid)
                        validatemarkdup_jobid = job.jobId
                        hold_mergereadgroups2library_jids += [job.jobId]
                        
                        #verify validation
                        job_name = my_name+'_VerifyValidateLibraryMarkDupBamFile_'+os.path.basename(library_markdup_bam)
                        cmd = '{} {} -v {}'.format(args.python, verify_sam_file_validation.__file__, library_markdup_validate)
                        job = my.run_job(cmd, job_name, job_dir, hold_jid=validatemarkdup_jobid, email=args.email)
                        verifyvalidatemarkdup_jobid = job.jobId
                        hold_mergereadgroups2library_jids += [job.jobId]
                        
        if not hold_mergereadgroups2library_jids:
            hold_mergereadgroups2library_jids = None
            
        if not args.NoLibraryBamMetrics:
            print 'Submitting jobs to calculate metrics of library bam files before removing duplicates ...'
            library_bams = [l for l in library_readgroup_bams]
            library_bam_metrics_job_dir = my.unique_dir('library_pre-markdup_metrics', my_job_dir)
            job_dir = library_bam_metrics_job_dir
            cmd = '{} {} -d {} --output_prefix {} -i {} -o {}'.format(
                args.python, bam_metrics.__file__, dirs['top'], 'library_pre-markdup_metrics/',' '.join(library_bams) ,dirs['library_bam_metrics'])
            job_name = my_name+'_library_pre-markdup_metrics'
            job = my.run_job(cmd, job_name, job_dir, hold_jid=hold_mergereadgroups2library_jids)
            library_pre_markdup_metrics_jobid = job.jobId
            
        if not args.NoLibraryMarkDupBamMetrics:
            print 'Submitting jobs to calculate metrics of library bam files after removing duplicates ...'
            library_markdup_bams = [library_markdup_bams[l] for l in library_markdup_bams]
            markdup_library_bam_metrics_job_dir = my.unique_dir('library_post-markdup_metrics', my_job_dir)
            job_dir = markdup_library_bam_metrics_job_dir
            cmd = '{} {} -d {} --output_prefix {} -i {} -o {}'.format(
                args.python, bam_metrics.__file__, dirs['top'], 'library_post-markdup_metrics/',' '.join(library_markdup_bams) ,dirs['library_bam_metrics'])
            job_name = my_name+'_library_post-markdup_metrics'
            job = my.run_job(cmd, job_name, job_dir, hold_jid=hold_mergereadgroups2library_jids)
            library_post_markdup_metrics_jobid = job.jobId
            
        #merge libraries by sample
        hold_mergelibraries2sample_jids = []
        
        if not args.NoMergeLibraries2Samples:
            print 'Submitting jobs to merge library bam files by sample ...'
            for sample_bam in sample_library_bams:
                
                sample_bai = my.bam_strip(sample_bam)+'.bai'
                sample_validate = sample_bam+'.validate'
                bams2merge = sample_library_bams[sample_bam]
                
                if len(bams2merge) > 0:
                    if args.overwrite_all_output_files or args.overwrite_merge_library_bams or not my.check_files(None, [(sample_validate,'No errors found\n')], [(sample_bam,sample_bai), (sample_bam,sample_validate)]):
                        if my.file_exists(sample_bam):
                            os.remove(sample_bam)
                        if my.file_exists(sample_bai):
                            os.remove(sample_bai)
                        if my.file_exists(sample_validate):
                            os.remove(sample_validate)
                        job_name = my_name+'_MergeLibraries_'+os.path.basename(sample_bam)
                        if len(bams2merge) == 1:
                            library_bam = bams2merge[0]
                            library_bai = my.bam_strip(library_bam)+'.bai'
                            library_validate = library_bam+'.validate'
                            cmd = 'mv {0} {1}; ln -fs {1} {0}; mv {2} {3}; ln -fs {3} {2}; mv {4} {5}; ln -fs {5} {4};'.format(
                                library_bam, sample_bam, library_bai, sample_bai, library_validate, sample_validate)
                            job = my.run_job(cmd, job_name, job_dir, hold_jid=hold_mergereadgroups2library_jids)
                            mergelibraries2sample_jobid = job.jobId
                            hold_mergelibraries2sample_jids += [job.jobId]
                        else:
                            multi_valued_arg_strings = 'INPUT='+' INPUT='.join(bams2merge)
                            picard_args = {'CREATE_INDEX': 'true', 'CREATE_MD5_FILE': 'true', 'SORT_ORDER': 'coordinate', 'ASSUME_SORTED': 'false', 'MERGE_SEQUENCE_DICTIONARIES': 'false', 'USE_THREADING': 'true', 'OUTPUT': sample_bam}
                            job = picard.run(tool='MergeSamFiles', picard_args=picard_args, multi_valued_arg_strings=multi_valued_arg_strings, job_dir=job_dir, hold_jid=hold_mergereadgroups2library_jids)
                            mergelibraries2sample_jobid = job.jobId
                            hold_mergelibraries2sample_jids += [job.jobId]
                            
                            #validate bam file with picard
                            job_name = my_name+'_ValidateSampleBamFile_'+os.path.basename(sample_bam)
                            picard_args = {'INPUT': sample_bam, 'OUTPUT': sample_validate}
                            job = picard.run(tool='ValidateSamFile', picard_args=picard_args, job_dir=job_dir, reference=args.reference, hold_jid=mergelibraries2sample_jobid)
                            validatesamplebamfile_jobid = job.jobId
                            hold_mergelibraries2sample_jids += [job.jobId]
                            
                            #verify validation
                            job_name = my_name+'_VerifyValidateSampleBamFile_'+os.path.basename(sample_bam)
                            cmd = '{} {} -v {}'.format(args.python, verify_sam_file_validation.__file__, sample_validate)
                            job = my.run_job(cmd, job_name, job_dir, hold_jid=validatesamplebamfile_jobid, email=args.email)
                            verifyvalidatesamplebamfile_jobid = job.jobId
                            hold_mergelibraries2sample_jids += [job.jobId]
                            
        if not hold_mergelibraries2sample_jids:
            hold_mergelibraries2sample_jids = None
            
        if not args.NoSampleBamMetrics:
            print 'Submitting jobs to calculate metrics of sample bam files ...'
            sample_bams = [l for l in sample_library_bams]
            sample_bam_metrics_job_dir = my.unique_dir('sample_metrics', my_job_dir)
            job_dir = sample_bam_metrics_job_dir
            cmd = '{} {} -d {} --output_prefix {} -i {} -o {}'.format(
                args.python, bam_metrics.__file__, dirs['top'], 'sample_metrics/',' '.join(sample_bams) ,dirs['sample_bam_metrics'])
            job_name = my_name+'_sample_metrics'
            job = my.run_job(cmd, job_name, job_dir, hold_jid=hold_mergelibraries2sample_jids)
            sample_metrics_jobid = job.jobId
            
    print 'All pypeline jobs submitted. Results will be in subdirectories of {}'.format(dirs['top'])
    
    return 0


if __name__ == "__main__": sys.exit(main())
