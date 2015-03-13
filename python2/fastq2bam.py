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
import logging
from warnings import warn
import my
import job
import qseq_metrics
import picard
import verify_sam_file_validation
import bam_metrics
import metrics_consolidate
import merge_validate_qseqs
import validate_brlmm_file
import verify_OK_file
import validate_vcf

class Fastq2BamError(Exception): pass

name = 'pypeline.fastq2bam'
version = '1.0β1'
copyright = '©2011-2014 Michael Yourshaw all rights reserved'
run_time = my.localtime_stamp()
time_squish = my.localtime_squish()

#-d /scratch1/tmp/myourshaw/gmdjen -m /scratch1/tmp/myourshaw/gmdjen/gmdjen_metadata.txt -i /data/storage-1-04/archive/myourshaw/HiSeq2k.Geschwindlab/120224_SN1083_0094_B_2012-045_C08EVACXX/Data/Intensities/BaseCalls/Unaligned/Project_*/Sample_*/*_*_L*_R[12]_*.fastq.gz
#-d /scratch1/tmp/myourshaw/gmdjen -m /scratch1/tmp/myourshaw/gmdjen/gmdjen_metadata.txt --fastqs /data/storage-1-02/archive/myourshaw/HiSeq2k.Pathology/120417_SN973_0083_AC0P8TACXX/Data/Intensities/BaseCalls/Unaligned/Project_GMD/Sample_GMD73/GMD73_CGATGT_L004_R[12]_001.fastq.gz
#-d /scratch1/tmp/myourshaw/gmdjen -m /scratch1/tmp/myourshaw/gmdjen/GMD_JEN_metadata.txt --fastqs /data/storage-1-04/archive/myourshaw/HiSeq2k.Freimerlab/120424_SN860_0135_A_2012-100_D0WCLACXX/Data/Intensities/BaseCalls/Unaligned/Project_GMD/Sample_GMD12A/GMD12A_CAGATC_L005_R[12]_001.fastq.gz

#one sample
#--NoLinkdatagen -d /scratch1/tmp/myourshaw/gmdjen -m /scratch1/tmp/myourshaw/gmdjen/GMD_JEN_metadata.txt --fastqs /data/storage-1-04/archive/myourshaw/HiSeq2k.Freimerlab/120424_SN860_0135_A_2012-100_D0WCLACXX/Data/Intensities/BaseCalls/Unaligned/Project_GMD/Sample_GMD12A/*_*_L*_R[12]_*.fastq.gz

#all samples
#--NoLinkdatagen -d /scratch1/tmp/myourshaw/gmdjen -m /scratch1/tmp/myourshaw/gmdjen/GMD_JEN_metadata.txt --email myourshaw@ucla.edu --fastqs /data/storage-1-04/archive/myourshaw/HiSeq2k.Freimerlab/120424_SN860_0135_A_2012-100_D0WCLACXX/Data/Intensities/BaseCalls/Unaligned/Project_*/Sample_*/*_*_L*_R[12]_*.fastq.gz /data/storage-1-01/archive/myourshaw/HiSeq2k.Freimerlab/120424_SN860_0136_B_2012-101_D0WTRACXX/Data/Intensities/BaseCalls/Unaligned/Project_*/Sample_*/*_*_L00?_R[12]_001.fastq.gz /data/storage-1-04/archive/myourshaw/HiSeq2k.Geschwindlab/120224_SN1083_0094_B_2012-045_C08EVACXX/Data/Intensities/BaseCalls/Unaligned/Project_*/Sample_*/*_*_L00?_R[12]_001.fastq.gz /data/storage-1-02/archive/myourshaw/HiSeq2k.Pathology/120417_SN973_0083_AC0P8TACXX/Data/Intensities/BaseCalls/Unaligned/Project_*/Sample_*/*_*_L00?_R[12]_001.fastq.gz

#mm_48_20130208 flowcell 1
#-d /scratch1/tmp/myourshaw/mm_48_20130208 -m /scratch1/tmp/myourshaw/mm_48_20130208/mm_48_20130208_metadata.txt --fastqs /scratch1/tmp/myourshaw/mm_48_20130208/fastqs/130308_SN608_0315_AC1BWKACXX_VAP020/Unaligned/Project_MM/Sample_*/*_*_L*_R[12]_*.fastq.gz

#mmjj_20130514
#--NoNovoalign --NoFixMateInformation -d /scratch1/tmp/myourshaw/Ymmjj_20130514 -m /scratch1/tmp/myourshaw/mmjj_20130514/mmjj_20130514_metadata.txt --fastqs /data/storage-1-04/archive/myourshaw/illumina/110616_SN860_0065_2011-101_B817EAABXX/Unaligned/Project_*/Sample_*/*_*_L*_R[12]_*.fastq.gz /data/storage-1-02/archive/myourshaw/illumina/110623_SN860_0067_2011-100R_A81MVKABXX/Unaligned/Project_*/Sample_*/*_*_L*_R[12]_*.fastq.gz /data/storage-1-00/archive/myourshaw/HiSeq2k.Freimerlab/110812_SN860_0079_2011-143_AC03YLACXX/Unaligned/Project_*/Sample_*/*_*_L*_R[12]_*.fastq.gz /data/storage-1-00/archive/myourshaw/HiSeq2k.Freimerlab/110929_SN860_0095_B_2011-173_D0A1YACXX/Unaligned/Project_*/Sample_*/*_*_L*_R[12]_*.fastq.gz /data/storage-1-02/archive/myourshaw/solexa/110511_SN430_0243_B817FLABXX/Unaligned/Project_*/Sample_*/*_*_L*_R[12]_*.fastq.gz /data/storage-1-02/archive/myourshaw/solexa/110715_SN973_0041_AC041EACXX/Unaligned/Project_*/Sample_*/*_*_L*_R[12]_*.fastq.gz /data/storage-1-04/archive/myourshaw/HiSeq2k.Freimerlab/120424_SN860_0135_A_2012-100_D0WCLACXX/Unaligned/Project_*/Sample_*/*_*_L*_R[12]_*.fastq.gz /data/storage-1-01/archive/myourshaw/HiSeq2k.Freimerlab/120424_SN860_0136_B_2012-101_D0WTRACXX/Unaligned/Project_*/Sample_*/*_*_L00?_R[12]_001.fastq.gz /data/storage-1-04/archive/myourshaw/HiSeq2k.Geschwindlab/120224_SN1083_0094_B_2012-045_C08EVACXX/Unaligned/Project_*/Sample_*/*_*_L00?_R[12]_001.fastq.gz /data/storage-1-02/archive/myourshaw/HiSeq2k.Pathology/120417_SN973_0083_AC0P8TACXX/Unaligned/Project_*/Sample_*/*_*_L00?_R[12]_001.fastq.gz /data/storage-1-02/archive/myourshaw/HiSeq2k.Pathology/120522_SN973_0090_BD120TACXX/Unaligned/Project_*/Sample_*/*_*_L00?_R[12]_001.fastq.gz /data/storage-1-02/archive/myourshaw/HiSeq2k.Pathology/120605_SN973_0093_BC0V4JACXX/Unaligned/Project_*/Sample_*/*_*_L00?_R[12]_001.fastq.gz /data/storage-1-03/archive/myourshaw/HiSeq2k.Pelligrini/130308_SN608_0315_AC1BWKACXX_VAP020/Unaligned/Project_MM/Sample_*/*_*_L*_R[12]_*.fastq.gz /data/storage-1-03/archive/myourshaw/HiSeq2k.Pelligrini/130405_SN608_0324_AC1DRMACXX_VAP021/Unaligned/Project_MM/Sample_*/*_*_L*_R[12]_*.fastq.gz /data/storage-1-03/archive/myourshaw/HiSeq2k.Pelligrini/130405_SN608_0325_BC1TAYACXX_VB013/Unaligned/Project_MM/Sample_*/*_*_L*_R[12]_*.fastq.gz /data/storage-1-03/archive/myourshaw/HiSeq2k.Pathology/130501_SN973_0182_AH0L2WADXX/Unaligned/Project_MM/Sample_*/*_*_L*_R[12]_*.fastq.gz /data/storage-1-03/archive/myourshaw/HiSeq2k.Pathology/130501_SN973_0183_BH0JE0ADXX/Unaligned/Project_MM/Sample_*/*_*_L*_R[12]_*.fastq.gz

#fhf_20140217 (multiple files per readgroup)
#--overwrite_all_output_files -d /scratch1/tmp/myourshaw/fhf_20140217 -m /scratch1/tmp/myourshaw/fhf_20140217/fhf_20140217_metadata.txt --targets /scratch1/tmp/myourshaw/resources/intervals/Nimblegen/nimblegen_HGSC_VCRome_b37_liftover.interval_list /scratch1/tmp/myourshaw/resources/intervals/EnsemblIntervals/current-ensembl-intervals/protein_coding_genes/ensembl_72_protein_coding_known_CDS_ess_5sr4_3sr13.interval_list --fastqs /scratch1/tmp/myourshaw/fhf_20140217/fastqs/Sample_13_97467/13_97467_GTCCGC_L001_R[12]_*.fastq.gz

#--email myourshaw@ucla.edu --overwrite_all_output_files -d /scratch1/tmp/myourshaw/sarc -m /scratch1/tmp/myourshaw/sarc/sarc_06052014_metadata.txt --bait /scratch1/tmp/myourshaw/resources/intervals/Nimblegen/SeqCapEZ_Exome_v3.0_Design_Annotation_files/SeqCap_EZ_Exome_v3_capture.interval_list --targets /scratch1/tmp/myourshaw/resources/intervals/Nimblegen/SeqCapEZ_Exome_v3.0_Design_Annotation_files/SeqCap_EZ_Exome_v3_capture.interval_list /scratch1/tmp/myourshaw/resources/intervals/EnsemblIntervals/current/protein_coding_genes/ensembl_74_protein_coding_known_CDS_ess_5sr4_3sr13.interval_list --fastqs /scratch1/tmp/myourshaw/sarc/fastqs/*.fastq.gz

#valley fever
#--email myourshaw@ucla.edu -d /scratch1/tmp/myourshaw/vf_20140828 -m /scratch1/tmp/myourshaw/vf_20140828/vf_20140828_metadata.txt --bait /scratch1/tmp/myourshaw/resources/intervals/Nimblegen/SeqCapEZ_Exome_v3.0_Design_Annotation_files/SeqCap_EZ_Exome_v3_capture.interval_list --fastqs /scratch1/tmp/myourshaw/vf_20140828/fastqs/Sample_KMC-*/*.fastq.gz

def gatk(dirs, tool, jar, params, job_name_prefix='', job_dir=None, java_mem='5g', processors=8, memory='6G', queue='all.q', tmp_dir=None, email=None, synchronous=False, hold_jid=None):
    if not job_dir:
        job_dir = dirs['gatk_jobs']
    if not tmp_dir:
        tmp_dir = dirs['tmp']
    if isinstance(params, (list, tuple)):
        params = ' '.join(params)
    cmd = 'java -Xmx{} -Djava.io.tmpdir="{}" -jar "{}" {}{}'.format(
        java_mem, tmp_dir, jar, params, ' -T '+tool if tool else '')
    job_name = '{}_{}'.format(tool if tool else re.sub(r'\.jar$','',os.path.basename(jar)), job_name_prefix)
    job = my.run_job(cmd, job_name, job_dir, email=email, processors=processors, memory=memory, queue=queue, synchronous=synchronous, hold_jid=hold_jid)
    #DEBUG
    #print cmd
    return job

def main():

    #command line arguments
    parser = argparse.ArgumentParser(parents=[my.default_parser()],
        description = 'align fastq files to readgroup bams and merge into library and sample bams',
        epilog = '{} version {} {}'.format(name, version, copyright))
    parser.add_argument('--metadata', '-m', required=True,
        help='metadata file to map machine, run, lane, barcode to adapters, library, sample and other readgroup info',)
    parser.add_argument('--fastqs', '-i', nargs='+', required=True,
        help='list of fastq files')
    parser.add_argument('--python',
        help='python executable (default: from config->DEFAULT.python)')
    parser.add_argument('--perl',
        help='perl executable (default: from config->DEFAULT.perl)')
    parser.add_argument('--reference', '-r',
        help='path to reference sequence fasta file (default: from config->reference.default_reference)')
    parser.add_argument('--novoalign',
        help='path to novoalign executable (default: from config->novocraft.novoalign)')
    parser.add_argument('--index',
        help='path to novoindex-produced index (default: from config->novocraft.default_index)')
    parser.add_argument('--gatk_jar',
        help='path to GenomeAnalysisTK jar')
    parser.add_argument('--bait','-B', metavar='BAIT_INTERVAL_LIST', nargs='?',
        help='bait interval list (default: config->interval.nimblegen_SeqCapEZ_Exome_v30)')
    parser.add_argument('--targets','-T', metavar='TARGET_INTERVAL_LIST', nargs='+',
        help='list of target interval_list files (default: [bait, config->interval.ensembl_protein_coding_genes, config->interval.refgene_interval_list])')
    parser.add_argument('--depth_of_coverage_intervals', nargs='*',
        help='list of interval_list files for depth of coverage (default: targets)')
    parser.add_argument('--depth_of_coverage_gene_list',
        help='list of  gene list file for depth of coverage (default: config->interval.refgene_gatk)')
    parser.add_argument('--samtools',
        help='path to samtools executable (default: from config->samtools.samtools)')
    parser.add_argument('--mpileup_m', type=int, default = 0.002,
        help='mpileup m parameter, minimum gapped reads for indel candidates (default: 1)')
    parser.add_argument('--mpileup_F', type=float, default = 1,
        help='mpileup m parameter, minimum fraction of gapped reads for candidates (default: 0.002)')
    parser.add_argument('--bcftools',
        help='path to bcftools executable (default: from config->samtools.bcftools)')
    parser.add_argument('--annotHapMap2L_b37',
        help='path to linkdatagen annotHapMap2L_b37 file (default: from config->linkdatagen.annotHapMap2L_b37)')
    parser.add_argument('--annotHapMap2',
        help='path to linkdatagen annotHapMap2 file (default: from config->linkdatagen.annotHapMap2)')
    parser.add_argument('--vcf2linkdatagen',
        help='path to vcf2linkdatagen.pl executable (default: from config->linkdatagen.vcf2linkdatagen)')
    parser.add_argument('--linkdatagen_mps',
        help='path to linkdatagen_mps.pl executable (default: from config->linkdatagen.linkdatagen_mps)')
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
    parser.add_argument('--overwrite_mpileups', action='store_true', default=False,
        help='redo mpileups, overwriting existing valid mpileup files')
    parser.add_argument('--overwrite_linkdatagen_pileups', action='store_true', default=False,
        help='redo sample pileups, overwriting existing valid sample pileup files')
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
    parser.add_argument('--Mpileups', action='store_true', default=False,
        help='Create pileups for samtools variant calling')
    parser.add_argument('--NoLinkdatagen', action='store_true', default=False,
        help='omit pileups and linkdatagen for homozygosity mapping')
    parser.add_argument('--NoFastQC', action='store_true', default=False,
        help='Do not run FastQC on fastq files')
    parser.add_argument('--fastqc_executable',
        help='path to FastQC executable (default: from config->fastqc.fastqc_executable)')

    args = parser.parse_args()

    dirs = my.create_default_directory_structure(args.top_dir)
    config = my.get_config(args)

    args.python = args.python if args.python else config.get('DEFAULT','python')
    if not my.is_exe(args.python):
        raise Fastq2BamError('{} is not the python executable'.format(args.python))
    args.perl = args.perl if args.perl else config.get('DEFAULT','perl')
    if not my.is_exe(args.perl):
        raise Fastq2BamError('{} is not the perl executable'.format(args.perl))
    args.novoalign = args.novoalign if args.novoalign else config.get('novocraft','novoalign')
    if not my.is_exe(args.novoalign):
        raise Fastq2BamError('{} is not a novoalign executable'.format(args.novoalign))
    args.index = args.index if args.index else config.get('novocraft','default_index')
    if not my.file_exists(args.index):
        raise Fastq2BamError('cannot find novolaign genome index file {}'.format(args.index))
    himem = bool(os.path.getsize(args.index) > 6500000000)
    args.reference = args.reference if args.reference else config.get('reference','default_decoy_reference')
    if not my.file_exists(args.reference):
        raise Fastq2BamError('cannot find reference sequence fasta file {}'.format(args.reference))
    args.gatk_jar = args.gatk_jar if args.gatk_jar else config.get('gatk','GenomeAnalysisTK')
    if not my.file_exists(args.gatk_jar):
        raise Fastq2BamError('cannot find GATK jar file {}'.format(args.gatk_jar))
    if not args.bait:
        #args.bait = config.get('interval','sureselect_50mb_interval_list')
        args.bait = config.get('interval','nimblegen_SeqCapEZ_Exome_v30')
    if not my.file_exists(args.bait):
        raise Fastq2BamError('cannot find bait file {}'.format(args.bait))
    if not args.targets:
        args.targets = [args.bait, config.get('interval','ensembl_protein_coding_genes'), config.get('interval','refgene_interval_list'),]
    for i in args.targets:
        if not my.file_exists(i):
            raise Fastq2BamError('cannot find target file {}'.format(i))
    args.depth_of_coverage_intervals = args.depth_of_coverage_intervals if args.depth_of_coverage_intervals else args.targets
    for i in args.depth_of_coverage_intervals:
        if not my.file_exists(i):
            raise Fastq2BamError('cannot find depth of coverage interval file {}'.format(i))
    args.depth_of_coverage_gene_list = args.depth_of_coverage_gene_list if args.depth_of_coverage_gene_list else config.get('interval','refgene_gatk')
    if not my.file_exists(args.depth_of_coverage_gene_list):
        raise Fastq2BamError('cannot find depth of coverage gene list file {}'.format(args.depth_of_coverage_gene_list))
    metrics_intervals_parameters = '--bait {} --targets {} --depth_of_coverage_intervals {} --depth_of_coverage_gene_list {}'.format(
        args.bait, ' '.join(args.targets), ' '.join(args.depth_of_coverage_intervals), args.depth_of_coverage_gene_list)
    args.samtools = args.samtools if args.samtools else config.get('samtools','samtools')
    if not my.is_exe(args.samtools):
        raise Fastq2BamError('{} is not a samtools executable'.format(args.samtools))
    args.bcftools = args.bcftools if args.bcftools else config.get('samtools','bcftools')
    if not my.is_exe(args.bcftools):
        raise Fastq2BamError('{} is not a bcftools executable'.format(args.bcftools))
    args.annotHapMap2L_b37 = args.annotHapMap2L_b37 if args.annotHapMap2L_b37 else config.get('linkdatagen','annotHapMap2L_b37')
    if not my.file_exists(args.annotHapMap2L_b37):
        raise Fastq2BamError('{} is not the annotHapMap2L_b37 file'.format(args.annotHapMap2L_b37))
    args.annotHapMap2 = args.annotHapMap2 if args.annotHapMap2 else config.get('linkdatagen','annotHapMap2')
    if not my.file_exists(args.annotHapMap2):
        raise Fastq2BamError('{} is not the annotHapMap2 file'.format(args.annotHapMap2))
    annotHapMap2_dir = os.path.dirname(args.annotHapMap2)
    annotHapMap2_file = os.path.basename(args.annotHapMap2)
    args.vcf2linkdatagen = args.vcf2linkdatagen if args.vcf2linkdatagen else config.get('linkdatagen','vcf2linkdatagen')
    if not my.file_exists(args.vcf2linkdatagen):
        raise Fastq2BamError('{} is not the vcf2linkdatagen perl script'.format(args.vcf2linkdatagen))
    args.linkdatagen_mps = args.linkdatagen_mps if args.linkdatagen_mps else config.get('linkdatagen','linkdatagen_mps')
    if not my.file_exists(args.linkdatagen_mps):
        raise Fastq2BamError('{} is not the linkdatagen_mps perl script'.format(args.linkdatagen_mps))
    
    regex = config.get('picard', 'READ_NAME_REGEX')

    my_name = 'fastq2bam'
    my_job_dir =  my.unique_dir('fastq2bam_'+time_squish, dirs['jobs'])
    
    this_log = os.path.join(my_job_dir,my_name+'.'+time_squish+'.log')
    logging.basicConfig(filename=this_log, format='%(asctime)s\t%(message)s', level=logging.DEBUG)
    #header for log
    logging.info('jobID\tcommand')
    #effective args
    logging.info('{}\t{}'.format(0, args))

    #get readgroup info from metadata file and first lines of fastq files
    #align with novoalign
    #fix mate information with picard
    print 'Processing metadata ...'
    readgroups = my.get_fastq_readgroups(args.fastqs, args.metadata)
    #hold library merge until novoalign and part merge complete
    hold_rg_jids = []
    
    #FastQC metrics
    if not args.NoFastQC:
        args.fastqc_executable = args.fastqc_executable if args.fastqc_executable else config.get('fastqc','fastqc_executable')
        if not my.file_exists(args.fastqc_executable):
            warn ('Cannot find FastQC executable. Skipping FastQC.')
        else:
            print 'Submitting FastQC job'
            fastqc_jobid = None
            job_dir = my.unique_dir('fastqc', my_job_dir)
            cmd = "{} --outdir {} {};".format(
                args.fastqc_executable, dirs['fastq_metrics'], ' '.join(args.fastqs))
            job_name = my_name+'_fastqc_'
            #print cmd
            job = my.run_job(cmd, job_name, job_dir)
            fastqc_jobid = job.jobId
            logging.info('{}\t{}'.format(job.jobId, cmd))


    print 'Submitting jobs to merge readgroup bam part files by readgroup ...'
    
    #align, merge, validate all readgroups
    if not readgroups:
        raise Fastq2BamError('No input files that match metadata.')
    for readgroup in readgroups:
        rg = readgroups[readgroup]
        if rg.do_not_align:
            continue
        #fastq_reads is a dict of read files
        #fastq_reads.keys = ('read', 'path')
        #fastq_reads.values = fastq peek data from first record of file
        if len(rg.fastq_reads) == 0:
            raise Fastq2BamError('No read files for {}'.format(rg))
        read_keys = rg.fastq_reads.keys()
        read_keys.sort()
        #a list of all files (both reads) of the readgroup
        fastq_read_files = [r[1] for r in read_keys]
        #expect one (single end) or two (paired end) reads, but the reads might not be numbered with id '1' and '2'
        #for example with long Illumina indexes, the data reads will be '1' and '3'
        #there can be multiple paths for each read if fastq was split
        fastq_read_ids = list(set([k[0] for k in read_keys]))
        fastq_read_ids.sort()
        #there must be 1 (single end) or 2 (paired end) read keys
        #and each must have the same number of files
        if len(fastq_read_ids) > 2:
            raise Fastq2BamError('More than two instrument read files for {}'.format(rg))
        #dict key = read id, value = [list of paths]
        read_paths = { k: [r[1] for r in read_keys if r[0] == k] for k in fastq_read_ids}
        #require same number of files for each of two paired end reads
        if len(fastq_read_ids) == 2:
             if len(read_paths[fastq_read_ids[0]]) != len(read_paths[fastq_read_ids[1]]):
                raise Fastq2BamError('Unequal number of read files for read {} and read {} reads for {}'.format(fastq_read_ids[0], fastq_read_ids[1], rg))
        #sort by read,path
        #assumes that both reads of paired end run will sort to correspond with each other 
        for i in fastq_read_ids:
            read_paths[i].sort
        #ordered tuples of paths (read1, [read2])
        readgroup_paths = [(read_paths[fastq_read_ids[0]][i], read_paths[fastq_read_ids[1]][i] if len(fastq_read_ids) == 2 else '') for i in range(len(read_paths[fastq_read_ids[0]]))]
        #setup novoalign parameters for each read file or pair of read filess in the readgroup
        #dict key=(fastq1,fastq2) value=bam_part
        rg.readgroup_bam_parts = {}
        rg.novoalign_cmds = []
        #compose RG string for bam header
        #run date from one of the fastq files if not specified in metadata
        if not rg.run_date:
            rg.run_date = strftime('%Y-%m-%d', localtime(os.path.getctime(fastq_read_files[0])))
        #defaults
        if not rg.platform:
            rg.platform = 'ILLUMINA'
        if not rg.instrument_model:
            rg.instrument_model = 'Illumina_HiSeq_2500'
        if not rg.sequencing_center:
            rg.sequencing_center = 'UCLA'
        if not my.is_number(rg.predicted_median_insert_size):
            rg.predicted_median_insert_size = 200
        rg.platform_unit = '{}.{}{}.{}{}'.format(rg.machine, rg.run, '.' + rg.flowcell if rg.flowcell else '', rg.lane, '.' + rg.barcode if rg.barcode else '')
        rg.readgroup_id = '{}.{}.{}{}'.format(rg.machine, rg.run, rg.lane, '.' + rg.barcode if rg.barcode else '')
        #determine read length from first read of fastq file
        #to calculate novoalign parameter for min quality bases
        #first read, first file
        read1_length = len(rg.fastq_reads[(fastq_read_ids[0], read_paths[fastq_read_ids[0]][0])][0]['sequence'])
        #second read, first file
        read2_length = len(rg.fastq_reads[(fastq_read_ids[1], read_paths[fastq_read_ids[1]][0])][0]['sequence']) if len(fastq_read_ids) == 2 else 0
        min_quality_bases = int(round((read1_length if read2_length == 0 else min(read1_length, read2_length))/2.0)) #novoalign default is ~20; this yields 25 for 50 base reads and 50 for 100 base reads
        rg.readgroup_description = 'sample={};original_sample={};tumor_status={};library={};library_protocol={};sequencing_library={};sequencing_center={};platform={};instrument_model={};platform_unit={};run_date={};run_folder={};machine={};run={};flowcell={};lane={};barcode={};read_length={};predicted_median_insert_size={};adapters={};fastq_files={}'.format(
            rg.sample, rg.original_sample, rg.tumor_status, rg.library, rg.library_protocol, rg.sequencing_library, rg.sequencing_center, rg.platform, rg.instrument_model, rg.platform_unit, rg.run_date, rg.run_folder, rg.machine, rg.run, rg.flowcell, rg.lane, rg.barcode, ','.join([str(read1_length), str(read2_length)]), rg.predicted_median_insert_size, ','.join((rg.adapter_1, rg.adapter_2)), ','.join(fastq_read_files))
        if rg.additional_readgroup_description:
            rg.readgroup_description = rg.readgroup_description + rg.addtional_readgroup_description + '' if rg.addtional_readgroup_description.endswith(';') else ';'
        #RG is the readgroup specification used in the Novoalign command line
        rg.RG = """@RG\\tID:{}\\tCN:{}\\tDT:{}\\tLB:{}\\tPI:{}\\tPL:{}\\tPU:{}\\tSM:{}\\tDS:{}""".format(
            rg.readgroup_id, rg.sequencing_center, rg.run_date, rg.library, rg.predicted_median_insert_size, rg.platform, rg.platform_unit, rg.sample, rg.readgroup_description)

        #novoalign output bam (merged if the readgroup fastqs/bams were split)
        rg.readgroup_bam = os.path.join(dirs['readgroup_bams'], '{}.{}.{}{}.novoalign.bam'.format(
            rg.machine, rg.run, rg.lane, '.'+rg.barcode if rg.barcode else ''))
        rg.fixmate_bam = my.bam_strip(rg.readgroup_bam)+'.fixmate.bam'
        rg.fixmate_bams += [rg.fixmate_bam]
        rg.fixmate_baits += [rg.capture_bait_interval_list]
        rg.fixmate_bai = my.bam_strip(rg.fixmate_bam) + '.bai'
        rg.fixmate_bam_validate = rg.fixmate_bam + '.validate'
            
        #hold merging or moving readgroup parts to readgroup until novoalign and post-processing complete
        hold_rg_part_jids = []
        for rpx in xrange(len(readgroup_paths)):
            rpx_str = '{:0=3}'.format(rpx)
            fastq1, fastq2 = readgroup_paths[rpx]
            
            #novoalign output bam (may be for part of a readgroup if fastq files were split)
            readgroup_bam_part = os.path.join(dirs['readgroup_bam_parts'], '{}.{}.{}{}.{}.novoalign.bam'.format(
                rg.machine, rg.run, rg.lane, '.'+rg.barcode if rg.barcode else '', rpx_str))
            rg.readgroup_bam_parts[(fastq1,fastq2)] = readgroup_bam_part
            #picard fixmate output files
            #and fastq-keyed dict of the parts
            fixmate_bam_part = my.bam_strip(readgroup_bam_part)+'.fixmate.bam'
            rg.fixmate_bam_parts[(fastq1,fastq2)] = fixmate_bam_part
            fixmate_bai_part = my.bam_strip(fixmate_bam_part) + '.bai'
            rg.fixmate_bai_parts[(fastq1,fastq2)] = fixmate_bai_part
            fixmate_bam_validate_part = fixmate_bam_part + '.validate'
            rg.fixmate_bam_validate_parts[(fastq1,fastq2)] = fixmate_bam_validate_part
            
            #novoalign fastq parts files to readgroup bam parts
            if args.overwrite_all_output_files or args.overwrite_novoalign_fixmate_bams or not my.check_files([rg.fixmate_bam], [(rg.fixmate_bam_validate,'No errors found\n')], [(rg.fixmate_bam,rg.fixmate_bai), (rg.fixmate_bam,rg.fixmate_bam_validate)]):
                if not args.NoNovoalign:
                    novoalign_jobid = None
                    job_dir = my.unique_dir('novoalign', my_job_dir)
                    cmd = "{} -d {} --mmapoff -f {} {} --hdrhd 1 -F ILM1.8 --ILQ_SKIP -l {} -H -t 17,3 --Q2Off -a {} {} -o SAM '{}' -k | {} view -S1  -  > {};".format(
                        args.novoalign, args.index, fastq1, fastq2, min_quality_bases, rg.adapter_1, rg.adapter_2, rg.RG, args.samtools, readgroup_bam_part)
                    job_name = my_name+'_novoalign_'+rg.readgroup_id+'.'+rpx_str
                    #print cmd
                    job = my.run_job(cmd, job_name, job_dir, memory = '24G' if himem else '7G')
                    novoalign_jobid = job.jobId
                    hold_rg_part_jids += [job.jobId]
                    hold_rg_jids += [job.jobId]
                    logging.info('{}\t{}'.format(job.jobId, cmd))
                    
                    #fix mate information with picard, overwriting the bam file created by novoalign
                    #this step sorts the bam file, creates an index,
                    #and ensures that all mate-pair information is in sync between each read and its mate pair.
                    job_name = my_name+'_FixMateInformation_'+rg.readgroup_id+'.'+rpx_str
                    picard_args = {'CREATE_INDEX': 'true', 'CREATE_MD5_FILE': 'true', 'SORT_ORDER': 'coordinate', 'INPUT': readgroup_bam_part, 'OUTPUT': fixmate_bam_part}
                    job = picard.run(tool='FixMateInformation', picard_args=picard_args, job_dir=job_dir, hold_jid=novoalign_jobid)
                    fixmateinformation_jobid = job.jobId
                    hold_rg_part_jids += [job.jobId]
                    hold_rg_jids += [job.jobId]
                    #logging.info('{}\tFixMateInformation {}'.format(job.jobId, picard_args))
                    
                    #validate bam file with picard
                    job_name = my_name+'_ValidateFixMateBamFile_'+rg.readgroup_id+'.'+rpx_str
                    picard_args = {'INPUT': fixmate_bam_part, 'OUTPUT': fixmate_bam_validate_part}
                    job = picard.run(tool='ValidateSamFile', picard_args=picard_args, job_dir=job_dir, reference=args.reference, hold_jid=fixmateinformation_jobid)
                    validatefixmatebamfile_jobid = job.jobId
                    hold_rg_part_jids += [job.jobId]
                    hold_rg_jids += [job.jobId]
                    #logging.info('{}\tValidateFixMateBamFile {}'.format(job.jobId, picard_args))
                    
                    #verify validation
                    job_name = my_name+'_VerifyValidateFixMateBamFile_'+rg.readgroup_id+'.'+rpx_str
                    cmd = '{} {} -v {}'.format(args.python, verify_sam_file_validation.__file__, fixmate_bam_validate_part)
                    job = my.run_job(cmd, job_name, job_dir, hold_jid=validatefixmatebamfile_jobid, email=args.email)
                    verifyvalidatefixmatebamfile_jobid = job.jobId
                    hold_rg_part_jids += [job.jobId]
                    hold_rg_jids += [job.jobId]
                    logging.info('{}\t{}'.format(job.jobId, cmd))
        #/rpx
        
        #merge readgroup parts by readgroup
        hold_merge_readgroup_parts2readgroup_jids = []
        if len(rg.fixmate_bam_parts) > 0 and len(hold_rg_part_jids) > 0:
            hold_merge_readgroup_bam_parts_jobids = []
            if my.file_exists(rg.fixmate_bam):
                os.remove(rg.fixmate_bam)
            if my.file_exists(rg.fixmate_bai):
                os.remove(rg.fixmate_bai)
            if my.file_exists(rg.fixmate_bam_validate):
                os.remove(rg.fixmate_bam_validate)
                
            job_name = my_name+'_MergeReadgroupParts_'+os.path.basename(rg.fixmate_bam)
            
            #move or merge
            if len(rg.fixmate_bam_parts) == 1:
                #move sole bam part to readgroup bam
                rg_part_bam = rg.fixmate_bam_parts.values()[0]
                rg_part_bai = my.bam_strip(rg_part_bam)+'.bai'
                rg_part_validate = rg_part_bam+'.validate'
                cmd = 'mv {0} {1}; ln -fs {1} {0}; mv {2} {3}; ln -fs {3} {2}; mv {4} {5}; ln -fs {5} {4};'.format(
                    rg_part_bam, rg.fixmate_bam, rg_part_bai, rg.fixmate_bai, rg_part_validate, rg.fixmate_bam_validate)
                job = my.run_job(cmd, job_name, job_dir, hold_jid=hold_rg_part_jids)
                merge_readgroup_parts2readgroup_jobid = job.jobId
                hold_merge_jobids += [job.jobId]
                hold_merge_readgroup_parts2readgroup_jids += [job.jobId]
                logging.info('{}\t{}'.format(job.jobId, cmd))
            else:
                #merge multiple bam parts to readgroup bam
                bamparts2merge = rg.fixmate_bam_parts.values()
                multi_valued_arg_strings = 'INPUT='+' INPUT='.join(bamparts2merge)
                picard_args = {'CREATE_INDEX': 'true', 'CREATE_MD5_FILE': 'true', 'SORT_ORDER': 'coordinate', 'ASSUME_SORTED': 'false', 'MERGE_SEQUENCE_DICTIONARIES': 'false', 'USE_THREADING': 'true', 'OUTPUT': rg.fixmate_bam}
                job = picard.run(tool='MergeSamFiles', picard_args=picard_args, multi_valued_arg_strings=multi_valued_arg_strings, job_dir=job_dir, hold_jid=hold_rg_part_jids)
                merge_readgroup_parts2readgroup_jobid = job.jobId
                hold_rg_jids += [job.jobId]
                hold_merge_readgroup_parts2readgroup_jids += [job.jobId]
                #logging.info('{}\tMergeSamFiles {}'.format(job.jobId, picard_args))
                
                #validate bam file with picard
                job_name = my_name+'_ValidateMergedReadgroupParts_'+os.path.basename(rg.fixmate_bam)
                picard_args = {'INPUT': rg.fixmate_bam, 'OUTPUT': rg.fixmate_bam_validate}
                job = picard.run(tool='ValidateSamFile', picard_args=picard_args, job_dir=job_dir, reference=args.reference, hold_jid=merge_readgroup_parts2readgroup_jobid)
                validate_merged_readgroup_parts_jobid = job.jobId
                hold_rg_jids += [job.jobId]
                hold_merge_readgroup_parts2readgroup_jids += [job.jobId]
                #logging.info('{}\tValidateSamFile {}'.format(job.jobId, picard_args))
                
                #verify validation
                job_name = my_name+'_VerifyValidateMergedReadgroupParts_'+os.path.basename(rg.fixmate_bam)
                cmd = '{} {} -v {}'.format(args.python, verify_sam_file_validation.__file__, rg.fixmate_bam_validate)
                job = my.run_job(cmd, job_name, job_dir, hold_jid=validate_merged_readgroup_parts_jobid, email=args.email)
                verify_validate_merged_readgroup_parts_jobid = job.jobId
                hold_rg_jids += [job.jobId]
                hold_merge_readgroup_parts2readgroup_jids += [job.jobId]
                logging.info('{}\t{}'.format(job.jobId, cmd))
                            
            if not hold_rg_jids:
                hold_rg_jids = None
        #/merge
        
        #output files that will be created after novoalign and merge
        rg.library_bam = os.path.join(dirs['library_bams'], '{}.library.bam'.format(rg.library))
        rg.library_bai = my.bam_strip(rg.library_bam)+'.bai'
        rg.library_bam_validate = rg.library_bam+'.validate'
        rg.library_markdup_bam = os.path.join(dirs['library_bams'], '{}.library.markdup.bam'.format(rg.library))
        rg.library_markdup_bai = my.bam_strip(rg.library_markdup_bam)+'.bai'
        rg.library_markdup_bam_validate = rg.library_markdup_bam+'.validate'
        rg.sample_bam = os.path.join(dirs['sample_bams'], '{}.sample.bam'.format(rg.sample))
        rg.sample_bai = my.bam_strip(rg.sample_bam)+'.bai'
        rg.sample_bam_validate = rg.sample_bam+'.validate'
        rg.study_bam = os.path.join(dirs['study_bams'], '{}.study.bam'.format(rg.study))
        rg.study_bai = my.bam_strip(rg.study_bam)+'.bai'
        rg.study_bam_validate = rg.study_bam+'.validate'
    #/readgroup
    
    #readgroup bam metrics
    if not args.NoReadgroupBamMetrics:
        print 'Submitting jobs to calculate metrics of readgroup bam files ...'
        readgroup_bam_metrics_job_dir = my.unique_dir('readgroup_metrics', my_job_dir)
        job_dir = readgroup_bam_metrics_job_dir
        cmd = '{} {} -d {} --output_prefix {} -i {} -o {} --NoDepthOfCoverage {} --bam_baits {}'.format(
            args.python, bam_metrics.__file__, dirs['top'], '', ' '.join(rg.fixmate_bams), dirs['readgroup_metrics'], metrics_intervals_parameters, ' '.join(list(rg.fixmate_baits)))
        job_name = my_name+'_readgroup_metrics'
        job = my.run_job(cmd, job_name, job_dir, hold_jid=hold_rg_jids)
        readgroup_bam_metrics_jobid = job.jobId
        logging.info('{}\t{}'.format(job.jobId, cmd))
    
    merge_rmdup_job_dir = my.unique_dir('merge_rmdup', my_job_dir)
    job_dir = merge_rmdup_job_dir
    hold_mergereadgroups2library_jids = []
    library_readgroup_bams = {b: sorted(set(my.flatten([rg.fixmate_bam for rg in readgroups.values() if rg.library_bam == b]))) for b in set(my.flatten([rg.library_bam for rg in readgroups.values()]))}
    library_markdup_bams = {b: set(my.flatten([rg.library_markdup_bam for rg in readgroups.values() if rg.library_bam == b])).pop() for b in set(my.flatten([rg.library_bam for rg in readgroups.values()]))}
    sample_library_bams = {b: sorted(set(my.flatten([rg.library_markdup_bam for rg in readgroups.values() if rg.sample_bam == b]))) for b in set(my.flatten([rg.sample_bam for rg in readgroups.values()]))}
    study_sample_bams = {b: sorted(set(my.flatten([rg.sample_bam for rg in readgroups.values() if rg.study_bam == b]))) for b in set(my.flatten([rg.study_bam for rg in readgroups.values()]))}
    
    #merge readgroups by library
    hold_merge_jobids = []
    duplication_metrics_files = []
    
    if not args.NoMergeReadgroups2Libraries:
        print 'Submitting jobs to merge readgroup bam files by library and mark duplicates ...'
        for library_bam in [b for b in library_readgroup_bams if b]:
            
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
                        logging.info('{}\t{}'.format(job.jobId, cmd))
                    else:
                        multi_valued_arg_strings = 'INPUT='+' INPUT='.join(bams2merge)
                        picard_args = {'CREATE_INDEX': 'true', 'CREATE_MD5_FILE': 'true', 'SORT_ORDER': 'coordinate', 'ASSUME_SORTED': 'false', 'MERGE_SEQUENCE_DICTIONARIES': 'false', 'USE_THREADING': 'true', 'OUTPUT': library_bam}
                        job = picard.run(tool='MergeSamFiles', picard_args=picard_args, multi_valued_arg_strings=multi_valued_arg_strings, job_dir=job_dir, hold_jid=hold_rg_jids)
                        mergereadgroups2library_jobid = job.jobId
                        hold_merge_jobids += [job.jobId]
                        hold_mergereadgroups2library_jids += [job.jobId]
                        #logging.info('{}\tMergeSamFiles {}'.format(job.jobId, picard_args))
                        
                        #validate bam file with picard
                        job_name = my_name+'_ValidateLibraryBamFile_'+os.path.basename(library_bam)
                        picard_args = {'INPUT': library_bam, 'OUTPUT': library_validate}
                        job = picard.run(tool='ValidateSamFile', picard_args=picard_args, job_dir=job_dir, reference=args.reference, hold_jid=mergereadgroups2library_jobid)
                        validatelibrarybamfile_jobid = job.jobId
                        hold_merge_jobids += [job.jobId]
                        hold_mergereadgroups2library_jids += [job.jobId]
                        #logging.info('{}\tValidateSamFile {}'.format(job.jobId, picard_args))
                        
                        #verify validation
                        job_name = my_name+'_VerifyValidateLibraryBamFile_'+os.path.basename(library_bam)
                        cmd = '{} {} -v {}'.format(args.python, verify_sam_file_validation.__file__, library_validate)
                        job = my.run_job(cmd, job_name, job_dir, hold_jid=validatelibrarybamfile_jobid, email=args.email)
                        verifyvalidatelibrarybamfile_jobid = job.jobId
                        hold_merge_jobids += [job.jobId]
                        hold_mergereadgroups2library_jids += [job.jobId]
                        logging.info('{}\t{}'.format(job.jobId, cmd))
                        
                #MarkDuplicates
                library_markdup_bam = library_markdup_bams[library_bam]
                library_markdup_bai = my.bam_strip(library_markdup_bam)+'.bai'
                library_markdup_validate = library_markdup_bam+'.validate'
                metrics_file = os.path.join(dirs['library_duplication_metrics'], '{}.duplication.metrics'.format(os.path.basename(library_markdup_bam)))
                duplication_metrics_files.append(metrics_file)
                if args.overwrite_all_output_files or args.overwrite_markdup_library_bams or not my.check_files(None, [(library_markdup_validate,'No errors found\n')], [(library_markdup_bam,library_markdup_bai), (library_markdup_bam,library_markdup_validate)]):
                    job_name = my_name+'_MarkDuplicates_'+os.path.basename(library_markdup_bam)
                    picard_args = {'CREATE_INDEX': 'true', 'CREATE_MD5_FILE': 'true', 'REMOVE_DUPLICATES': 'false', 'MAX_FILE_HANDLES_FOR_READ_ENDS_MAP': 8000, 'READ_NAME_REGEX': '"'+regex+'"', 'OPTICAL_DUPLICATE_PIXEL_DISTANCE': 100, 'METRICS_FILE': metrics_file, 'INPUT': library_bam, 'OUTPUT': library_markdup_bam}
                    job = picard.run(tool='MarkDuplicates', picard_args=picard_args, job_dir=job_dir, hold_jid=hold_merge_jobids)
                    markdup_jobid = job.jobId
                    hold_mergereadgroups2library_jids += [job.jobId]
                    #logging.info('{}\tMarkDuplicates {}'.format(job.jobId, picard_args))
                    
                    #validate bam file with picard
                    job_name = my_name+'_ValidateLibraryMarkDupBamFile_'+os.path.basename(library_markdup_bam)
                    picard_args = {'INPUT': library_markdup_bam, 'OUTPUT': library_markdup_validate}
                    job = picard.run(tool='ValidateSamFile', picard_args=picard_args, job_dir=job_dir, reference=args.reference, hold_jid=markdup_jobid)
                    validatemarkdup_jobid = job.jobId
                    hold_mergereadgroups2library_jids += [job.jobId]
                    #logging.info('{}\tValidateSamFile {}'.format(job.jobId, picard_args))
                    
                    #verify validation
                    job_name = my_name+'_VerifyValidateLibraryMarkDupBamFile_'+os.path.basename(library_markdup_bam)
                    cmd = '{} {} -v {}'.format(args.python, verify_sam_file_validation.__file__, library_markdup_validate)
                    job = my.run_job(cmd, job_name, job_dir, hold_jid=validatemarkdup_jobid, email=args.email)
                    verifyvalidatemarkdup_jobid = job.jobId
                    hold_mergereadgroups2library_jids += [job.jobId]
                    logging.info('{}\t{}'.format(job.jobId, cmd))
                    
    if not hold_mergereadgroups2library_jids:
        hold_mergereadgroups2library_jids = None
        
    if not args.NoLibraryBamMetrics:
        print 'Submitting jobs to calculate metrics of library bam files before removing duplicates ...'
        library_bams = [l for l in library_readgroup_bams if l]
        library_bam_metrics_job_dir = my.unique_dir('library_pre-markdup_metrics', my_job_dir)
        job_dir = library_bam_metrics_job_dir
        cmd = '{} {} -d {} --output_prefix {} -i {} -o {} --NoDepthOfCoverage {}'.format(
            args.python, bam_metrics.__file__, dirs['top'], '',' '.join(library_bams), dirs['library_pre-markdup_metrics'], metrics_intervals_parameters)
        job_name = my_name+'_library_pre-markdup_metrics'
        job = my.run_job(cmd, job_name, job_dir, hold_jid=hold_mergereadgroups2library_jids)
        library_pre_markdup_metrics_jobid = job.jobId
        logging.info('{}\t{}'.format(job.jobId, cmd))
        
    if not args.NoLibraryMarkDupBamMetrics:
        print 'Submitting jobs to calculate metrics of library bam files after removing duplicates ...'
        library_markdup_bams = [library_markdup_bams[l] for l in library_markdup_bams if l]
        markdup_library_bam_metrics_job_dir = my.unique_dir('library_post-markdup_metrics', my_job_dir)
        job_dir = markdup_library_bam_metrics_job_dir
        cmd = '{} {} -d {} --output_prefix {} -i {} -o {} --NoDepthOfCoverage {}'.format(
            args.python, bam_metrics.__file__, dirs['top'], '',' '.join(library_markdup_bams), dirs['library_post-markdup_metrics'], metrics_intervals_parameters)
        job_name = my_name+'_library_post-markdup_metrics'
        job = my.run_job(cmd, job_name, job_dir, hold_jid=hold_mergereadgroups2library_jids)
        library_post_markdup_metrics_jobid = job.jobId
        logging.info('{}\t{}'.format(job.jobId, cmd))
        
        #consolidated_inputs = duplication_metrics_files
        consolidated_inputs = os.path.join(dirs['library_duplication_metrics'], '*')
        consolidated_dir = os.path.join(dirs['library_duplication_metrics'],'consolidated')
        my.makedir(consolidated_dir)
        job_name = my_name+'_consolidate'
        consolidated_job_dir = my.unique_dir('library_duplication_metrics', my_job_dir)
        cmd = '{} {} --input {} --output_dir {}'.format(
            args.python, metrics_consolidate.__file__, consolidated_inputs, consolidated_dir)
        job = my.run_job(cmd, job_name, job_dir, synchronous=False, hold_jid = hold_mergereadgroups2library_jids)
        consolidated_jobid = job.jobId
        logging.info('{}\t{}'.format(job.jobId, cmd))

    #merge libraries by sample
    hold_mergelibraries2sample_jids = []
    
    if not args.NoMergeLibraries2Samples:
        print 'Submitting jobs to merge library bam files by sample ...'
        for sample_bam in [b for b in sample_library_bams if b]:
            
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
                        logging.info('{}\t{}'.format(job.jobId, cmd))
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
                        logging.info('{}\t{}'.format(job.jobId, cmd))
                        
    if not hold_mergelibraries2sample_jids:
        hold_mergelibraries2sample_jids = None
        
    sample_bams = [l for l in sample_library_bams if l]
        
    if not args.NoSampleBamMetrics:
        print 'Submitting jobs to calculate metrics for sample bam files ...'
        #picard metrics and depth of coverage
        sample_bam_metrics_job_dir = my.unique_dir('sample_metrics', my_job_dir)
        job_dir = sample_bam_metrics_job_dir
        cmd = '{} {} -d {} --output_prefix {} -i {} -o {} {}'.format(
            args.python, bam_metrics.__file__, dirs['top'], '',' '.join(sample_bams), dirs['sample_metrics'], metrics_intervals_parameters)
        job_name = my_name+'_sample_metrics'
        job = my.run_job(cmd, job_name, job_dir, hold_jid=hold_mergelibraries2sample_jids)
        sample_metrics_jobid = job.jobId
        logging.info('{}\t{}'.format(job.jobId, cmd))
        
    #mpileup for improved indel calling per novoalign 3 release notes
    #From Colin Hercus 2014-02-25
    #With regard -m & -F, the problem with samtools is that the default -F is too low
    #and as we mentioned they may report the wrong indel. By setting the fraction
    #higher we avoid the problem. I think the -F 0.07 is about right for any read
    #length. We didn't try mixing this with realignment but I expect if you use an
    #indel realigner the -F will be less critical. We kept -m at 1 as one read with
    #indel would call the indel even if read depth was only 1 read and this boosted
    #the TP rate. The issue with long indels as that unless the indel is near the
    #centre of the read it will get soft clipped so alignments that include indel are
    #reduced as indel length increases. If you do use -m 1 and do find indel with
    #read depth 1 then you need to look at soft clipped near the indel. You might use
    #an indel caller that looks for clipped alignments such PINDEL. 

    #process all samples together
    hold_mpileup_jids = []
    if args.Mpileups:
        
            print 'Submitting jobs to call variants with mpileup ...'
            
            mpileup_job_dir = my.unique_dir('mpileup', my_job_dir)
            job_dir = mpileup_job_dir
            sample_bams_list = ' '.join([l for l in sample_bams if l])
            pileup_vcf = os.path.join(dirs['mpileups'], os.path.basename(args.top_dir)+'.mpileup.vcf')
            pileup_vcf_validate = pileup_vcf+'.validate'
            
            #mpileup
            job_name = my_name+'_mpileup_'+os.path.basename(pileup_vcf)
            cmd = '{} mpileup -C50 -D -S -d10000 -m {} -F {} -q13 -guf {} {} | {} view -bvcg - > {};'.format(
                args.samtools, args.mpileup_m, args.mpileup_F, args.reference, sample_bams_list, args.bcftools, pileup_vcf)
            job = my.run_job(cmd, job_name, job_dir, hold_jid=hold_mergelibraries2sample_jids)
            mpileup_jobid = job.jobId
            hold_mpileup_jids += [job.jobId]
            logging.info('{}\t{}'.format(job.jobId, cmd))
            
            ##bcftools
            #job_name = my_name+'_bcftools_'+os.path.basename(pileup_bcf_vcf)
            #cmd = '{} view -cg -t0.5 {} > {};'.format(args.bcftools, pileup_bcf, pileup_bcf_vcf)
            #job = my.run_job(cmd, job_name, job_dir, hold_jid=mpileup_jobid)
            #bcftools_jobid = job.jobId
            #hold_mpileup_jids += [job.jobId]
            #logging.info('{}\t{}'.format(job.jobId, cmd))
            
            #validate pileup.vcf
            job_name = my_name+'_ValidateVcfFile_'+os.path.basename(pileup_vcf)
            cmd = '{} {} -i {} -o {}'.format(args.python, validate_vcf.__file__, pileup_vcf, pileup_vcf_validate)
            job = my.run_job(cmd, job_name, job_dir, hold_jid=mpileup_jobid)
            validatesamplevcffile_jobid = job.jobId
            hold_mpileup_jids += [job.jobId]
            logging.info('{}\t{}'.format(job.jobId, cmd))
            
            #verify pileup.vcf validation
            job_name = my_name+'_VerifyValidateSampleVcfFile_'+os.path.basename(pileup_vcf)
            cmd = '{} {} -v {}'.format(args.python, verify_OK_file.__file__, pileup_vcf_validate)
            job = my.run_job(cmd, job_name, job_dir, hold_jid=validatesamplevcffile_jobid)
            verifyvalidatesamplevcffile_jobid = job.jobId
            hold_mpileup_jids += [job.jobId]
            logging.info('{}\t{}'.format(job.jobId, cmd))
                
    #pileups, linkdatagen, and plink for homozygosity mapping
    hold_linkdatagen_jids = []
    if not args.NoLinkdatagen:
        print 'Submitting jobs to find homozygosity blocks with linkdatagen and plink ...'
        linkdatagen_job_dir = my.unique_dir('linkdatagen', my_job_dir)
        job_dir = linkdatagen_job_dir
        #file to tell linkdatagen to use first (only) genotype in sample_bcf_vcf
        linkdatagen_whichsamplesfile = os.path.join(dirs['linkdatagen'], '1.ws')
        with open(linkdatagen_whichsamplesfile, 'w') as ws:
            ws.write('1')
        for sample_bam in sample_bams:
            sample_bcf = os.path.join(dirs['linkdatagen'], os.path.basename(sample_bam)+'.bcf')
            sample_bcf_vcf = sample_bcf+'.vcf'
            sample_bcf_vcf_validate = sample_bcf_vcf+'.validate'
            sample_brlmm = sample_bcf_vcf+'.brlmm'
            sample_brlmm_validate = sample_brlmm+'.validate'
            sample_prefix = my.r_strip(my.bam_strip(sample_bam), '.sample')
            sample_name = os.path.basename(sample_prefix)
            linkdatagen_mps_stdout = os.path.join(dirs['linkdatagen'], sample_name+'.Ped_HapMap2_pl.out')
            sample_plink_dir = os.path.join(dirs['plink'], sample_name)
            plink_makebed_stdout = os.path.join(dirs['plink'], sample_name+'.plink_makebed.out')
            plink_homozyg_stdout = os.path.join(dirs['plink'], sample_name+'.plink_homozyg.out')
            #fake ped file for linkdatagen (all female to get X)
            linkdatagen_pedfile = os.path.join(dirs['linkdatagen'], sample_name+'.ped')
            with open(linkdatagen_pedfile,'w') as ped:
                ped.write('0001\t001\t0\t0\t2\t2')
            if args.overwrite_all_output_files or args.overwrite_linkdatagen_pileups or not my.check_files(files_exist=[sample_bcf,sample_bcf_vcf,sample_brlmm], files_contain=[(sample_bcf_vcf_validate,'OK'),(sample_brlmm_validate,'OK')], files_not_before=[(sample_bcf,sample_bam),(sample_bcf_vcf,sample_bcf),(sample_brlmm,sample_bcf_vcf)]):
                
                #mpileup for linkdatagen
                job_name = my_name+'_mpileup_'+os.path.basename(sample_bcf)
                cmd = '{} mpileup -C50 -d10000 -q13 -g -f {} -l {} {} > {};'.format(args.samtools, args.reference, args.annotHapMap2L_b37, sample_bam, sample_bcf)
                job = my.run_job(cmd, job_name, job_dir, hold_jid=hold_mergelibraries2sample_jids)
                mpileup_jobid = job.jobId
                hold_linkdatagen_jids += [job.jobId]
                logging.info('{}\t{}'.format(job.jobId, cmd))
                
                #bcftools for linkdatagen
                job_name = my_name+'_bcftools_'+os.path.basename(sample_bcf_vcf)
                cmd = '{} view -cg -t0.5 {} > {};'.format(args.bcftools, sample_bcf, sample_bcf_vcf)
                job = my.run_job(cmd, job_name, job_dir, hold_jid=mpileup_jobid)
                bcftools_jobid = job.jobId
                hold_linkdatagen_jids += [job.jobId]
                logging.info('{}\t{}'.format(job.jobId, cmd))
                
                #validate sample_bcf_vcf for linkdatagen
                job_name = my_name+'_ValidateBcfVcfFile_'+os.path.basename(sample_bcf_vcf)
                cmd = '{} {} -i {} -o {}'.format(args.python, validate_vcf.__file__, sample_bcf_vcf, sample_bcf_vcf_validate)
                job = my.run_job(cmd, job_name, job_dir, hold_jid=bcftools_jobid)
                validatesamplebcfvcffile_jobid = job.jobId
                hold_linkdatagen_jids += [job.jobId]
                logging.info('{}\t{}'.format(job.jobId, cmd))
                
                #verify sample_bcf_vcf validation for linkdatagen
                job_name = my_name+'_VerifyValidateSampleBcfVcfFile_'+os.path.basename(sample_bcf_vcf)
                cmd = '{} {} -v {}'.format(args.python, verify_OK_file.__file__, sample_bcf_vcf_validate)
                job = my.run_job(cmd, job_name, job_dir, hold_jid=validatesamplebcfvcffile_jobid)
                verifyvalidatesamplebcfvcffile_jobid = job.jobId
                hold_linkdatagen_jids += [job.jobId]
                logging.info('{}\t{}'.format(job.jobId, cmd))
                
                #vcf2linkdatagen
                job_name = my_name+'_vcf2linkdatagen_'+os.path.basename(sample_brlmm)
                cmd = '{} {} -annotfile {} -mindepth 5 -missingness 0 {} > {};'.format(args.perl, args.vcf2linkdatagen, args.annotHapMap2, sample_bcf_vcf, sample_brlmm)
                job = my.run_job(cmd, job_name, job_dir, hold_jid=bcftools_jobid)
                vcf2linkdatagen_jobid = job.jobId
                hold_linkdatagen_jids += [job.jobId]
                logging.info('{}\t{}'.format(job.jobId, cmd))
                
                #validate sample_brlmm
                job_name = my_name+'_ValidateBrlmmFile_'+os.path.basename(sample_brlmm)
                cmd = '{} {} -i {} -o {}'.format(args.python, validate_brlmm_file.__file__, sample_brlmm, sample_brlmm_validate)
                job = my.run_job(cmd, job_name, job_dir, hold_jid=vcf2linkdatagen_jobid)
                validatesamplebrlmmfile_jobid = job.jobId
                hold_linkdatagen_jids += [job.jobId]
                logging.info('{}\t{}'.format(job.jobId, cmd))
                
                #verify sample_brlmm validation
                job_name = my_name+'_VerifyValidateSampleBrlmmFile_'+os.path.basename(sample_brlmm)
                cmd = '{} {} -v {}'.format(args.python, verify_OK_file.__file__, sample_brlmm_validate)
                job = my.run_job(cmd, job_name, job_dir, hold_jid=validatesamplebrlmmfile_jobid)
                verifyvalidatesamplebrlmmfile_jobid = job.jobId
                hold_linkdatagen_jids += [job.jobId]
                logging.info('{}\t{}'.format(job.jobId, cmd))
                
                #linkdatagen_mps
                job_name = my_name+'_linkdatagen_mps_'+os.path.basename(sample_brlmm)
                cmd = '{} {} -pedfile {} -whichsamplesfile {} -callfile {} -annot_dir {} -annotfile {} -binsize 0 -merr -prog pl -uninf -outputdir {} > {};'.format(
                    args.perl, args.linkdatagen_mps, linkdatagen_pedfile, linkdatagen_whichsamplesfile, sample_brlmm, annotHapMap2_dir, annotHapMap2_file, sample_plink_dir, linkdatagen_mps_stdout)
                job = my.run_job(cmd, job_name, job_dir, hold_jid=vcf2linkdatagen_jobid)
                linkdatagen_mps_jobid = job.jobId
                hold_linkdatagen_jids += [job.jobId]
                logging.info('{}\t{}'.format(job.jobId, cmd))
                
                #plink makebed
                job_name = my_name+'_plink_makebed_'+os.path.basename(sample_name)
                cmd = 'cd {}_plink/; plink --file plink --make-bed --out {} > {};'.format(sample_plink_dir, sample_name, plink_makebed_stdout)
                job = my.run_job(cmd, job_name, job_dir, hold_jid=linkdatagen_mps_jobid)
                plink_makebed_jobid = job.jobId
                hold_linkdatagen_jids += [job.jobId]
                logging.info('{}\t{}'.format(job.jobId, cmd))
                
                #plink homozyg
                job_name = my_name+'_plink_makebed_'+os.path.basename(sample_name)
                cmd = 'cd {}_plink/; plink --bfile {} --homozyg > {};'.format(sample_plink_dir, sample_name, plink_homozyg_stdout)
                job = my.run_job(cmd, job_name, job_dir, hold_jid=plink_makebed_jobid)
                plink_homozyg_jobid = job.jobId
                hold_linkdatagen_jids += [job.jobId]
                logging.info('{}\t{}'.format(job.jobId, cmd))
                
    print 'All pypeline jobs submitted. Results will be in subdirectories of {}'.format(dirs['top'])
    print 'Log of commands is in {}'.format(this_log)
    
    return 0


if __name__ == "__main__": sys.exit(main())
