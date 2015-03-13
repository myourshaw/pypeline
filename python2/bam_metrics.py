#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import ConfigParser #configparser in python 3
import glob
import tempfile
from time import localtime, strftime
import my
import job
import picard
import metrics_consolidate
import gatk

class BamMetricsError(Exception): pass

name = 'pypeline.bam_metrics'
version = '1.0β1'
copyright = '©2011-2013 Michael Yourshaw all rights reserved'
run_time = my.localtime_stamp()

#gmd
#--NoDepthOfCoverage -d /scratch0/tmp/myourshaw/gmd1 --output_prefix test -i /scratch0/tmp/myourshaw/gmd1/bams/readgroup_bams/*.novoalign.fixmate.bam

#1 small gmd
#--NoDepthOfCoverage -d /scratch0/tmp/myourshaw/gmd1 --output_prefix test -i /scratch0/tmp/myourshaw/gmd1/bams/readgroup_bams/HWI-ST430.243.6.ATCACGA.novoalign.fixmate.bam

#gmd35
#-d /scratch0/tmp/myourshaw/gmdmms -i /scratch0/tmp/myourshaw/gmd/bams/sample_bams/*.sample.recalibrated.bam /scratch0/tmp/myourshaw/mms/bams/sample_bams/*.sample.recalibrated.bam

#truseq aliz bipolar
#-d /scratch1/tmp/myourshaw/hlee --output_prefix sample_metrics/ --bait /scratch1/tmp/myourshaw/resources/intervals/TruSeq/TruSeq_exome_targeted_regions.txt.interval_list -i /scratch0/tmp/aliz/bipolar/data/truseqABCD_myourshaw/bams/sample_bams/*.sample.recalibrated.bam

#-d /scratch1/tmp/myourshaw/gmdjen --output_prefix  -i /scratch1/tmp/myourshaw/gmdjen/bams/readgroup_bams/HWI-ST0860.135.5.CAGATC.novoalign.fixmate.bam -o /scratch1/tmp/myourshaw/gmdjen/metrics/readgroup_bam_metrics

#-d /scratch1/tmp/myourshaw/test -i /scratch1/tmp/myourshaw/gmdjen/bams/readgroup_bams/HWI-ST973.83.5.CGATGT.novoalign.fixmate.bam -o /scratch1/tmp/myourshaw/test/metrics/readgroup_bam_metrics

#-d /scratch1/tmp/myourshaw/gmdjen --output_prefix  -i /scratch1/tmp/myourshaw/gmdjen/bams/readgroup_bams/HWI-ST0860.135.5.CAGATC.novoalign.fixmate.bam -o /scratch1/tmp/myourshaw/gmdjen/metrics/picard_bam_metrics/readgroup_metrics --bait /share/apps/myourshaw/resources/intervals/TruSeq/TruSeq_exome_targeted_regions.txt.interval_list --targets /share/apps/myourshaw/resources/intervals/TruSeq/TruSeq_exome_targeted_regions.txt.interval_list /share/apps/myourshaw/resources/intervals/EnsemblIntervals/intervals_Ensembl_genes_protein_coding_cds.interval_list /share/apps/myourshaw/resources/intervals/RefGeneIntervals/refGene.ucsc.GRCh37.cds_ess.interval_list --depth_of_coverage_intervals /share/apps/myourshaw/resources/intervals/TruSeq/TruSeq_exome_targeted_regions.txt.interval_list /share/apps/myourshaw/resources/intervals/EnsemblIntervals/intervals_Ensembl_genes_protein_coding_cds.interval_list /share/apps/myourshaw/resources/intervals/RefGeneIntervals/refGene.ucsc.GRCh37.cds_ess.interval_list --depth_of_coverage_gene_list /share/apps/myourshaw/resources/refgene.b37.sorted.txt

#-d /scratch1/tmp/myourshaw/mmjj_20130514 --output_prefix  -i /scratch1/tmp/myourshaw/mmjj_20130514/bams/readgroup_bams/HWI-ST0860.79.3.CGATGT.novoalign.fixmate.bam /scratch1/tmp/myourshaw/mmjj_20130514/bams/readgroup_bams/HWI-ST973.183.1.GGACTCCTCTCTCTAT.novoalign.fixmate.bam /scratch1/tmp/myourshaw/mmjj_20130514/bams/readgroup_bams/HWI-ST608.325.5.GCTACGCTCTCTCTAT.novoalign.fixmate.bam -o /scratch1/tmp/myourshaw/mmjj_20130514/metrics/picard_bam_metrics/readgroup_metrics --NoDepthOfCoverage --bait /share/apps/myourshaw/resources/intervals/TruSeq/TruSeq_exome_targeted_regions.txt.interval_list --targets /share/apps/myourshaw/resources/intervals/TruSeq/TruSeq_exome_targeted_regions.txt.interval_list /share/apps/myourshaw/resources/intervals/EnsemblIntervals/intervals_Ensembl_genes_protein_coding_cds.interval_list /share/apps/myourshaw/resources/intervals/RefGeneIntervals/refGene.ucsc.GRCh37.cds_ess.interval_list --depth_of_coverage_intervals /share/apps/myourshaw/resources/intervals/TruSeq/TruSeq_exome_targeted_regions.txt.interval_list /share/apps/myourshaw/resources/intervals/EnsemblIntervals/intervals_Ensembl_genes_protein_coding_cds.interval_list /share/apps/myourshaw/resources/intervals/RefGeneIntervals/refGene.ucsc.GRCh37.cds_ess.interval_list --depth_of_coverage_gene_list /share/apps/myourshaw/resources/refgene.b37.sorted.txt --bam_baits /scratch1/tmp/myourshaw/resources/intervals/SureSelect/SureSelect_All_Exon_50mb_with_annotation_b37_sorted.interval_list /scratch1/tmp/myourshaw/resources/intervals/Nextera/NexteraRapidCapture_ExpandedExome_TargetedRegions.interval_list /scratch1/tmp/myourshaw/resources/intervals/Nextera/NexteraRapidCapture_ExpandedExome_TargetedRegions.interval_list 

#--NoDepthOfCoverage -d /scratch1/tmp/myourshaw/mmjj_tmp --output_prefix   -i /scratch1/tmp/myourshaw/mmjj_20130514/bams/library_bams/GMD181a_TCCTGAGCTAGATCGC_Illumina_Nextera_ExpandedExome_Enrichment.library.bam -o /scratch1/tmp/myourshaw/mmjj_tmp/metrics/picard_bam_metrics/library_metrics/library_pre-markdup_metrics

#--NoPicardMetrics -d /scratch1/tmp/myourshaw/doc_test --output_prefix  -i /scratch1/tmp/myourshaw/vf_20140828/bams/readgroup_bams/HWI-ST1148.219.8.TGACCA.novoalign.fixmate.bam

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

def run(top_dir, input, out_dir=None, output_prefix=None, job_dir=None, reference=None,
        depth_of_coverage_intervals=None, depth_of_coverage_gene_list=None,
        bait=None, bam_baits = None, targets=None, regex=None, pixels=100,
        gatk_jar=None, samtools=None, addtional_metrics_to_consolidate=None,
        NoPicardMetrics=False, NoDepthOfCoverage=False):

    dirs = my.create_default_directory_structure(top_dir)
    config = my.get_config()

    python = config.get('DEFAULT','python')
    
    if not regex:
        regex = config.get('picard', 'READ_NAME_REGEX')
     
    my_name = 'bam_metrics'
    my_job_dir = my.unique_dir(my_name, job_dir) if job_dir else my.unique_dir(my_name, dirs['jobs'])
    
    if not out_dir:
        out_dir = os.path.join(dirs['metrics'],'bam_metrics')
    my.makedir(out_dir)

    if not output_prefix:
        output_prefix = ''
        
    #get defaults from configuration file
    reference = reference if reference else config.get('reference','default_decoy_reference')
    if not my.file_exists(reference):
        raise BamMetricsError('cannot find reference sequence fasta file {}'.format(reference))
    gatk_jar = gatk_jar if gatk_jar else config.get('gatk','GenomeAnalysisTK')
    if not my.file_exists(gatk_jar):
        raise BamMetricsError('cannot find GATK jar file {}'.format(gatk_jar))
    if not bait:
        bait = config.get('interval','nimblegen_SeqCapEZ_Exome_v30')
    if not my.file_exists(bait):
        raise BamMetricsError('cannot find bait file {}'.format(bait))
    if not targets:
        targets = [bait, config.get('interval','ensembl_protein_coding_genes'), config.get('interval','refgene_interval_list'),]
    for i in targets:
        if not my.file_exists(i):
            raise BamMetricsError('cannot find target file {}'.format(i))
    depth_of_coverage_intervals = depth_of_coverage_intervals if depth_of_coverage_intervals else targets
    for i in depth_of_coverage_intervals:
        if not my.file_exists(i):
            raise BamMetricsError('cannot find depth of coverage interval file {}'.format(i))
    depth_of_coverage_gene_list = depth_of_coverage_gene_list if depth_of_coverage_gene_list else config.get('interval','refgene_gatk')
    if not my.file_exists(depth_of_coverage_gene_list):
        raise BamMetricsError('cannot find depth of coverage gene list file {}'.format(depth_of_coverage_gene_list))
    if not samtools:
        samtools = config.get('samtools','samtools')

    #unglob and sort bams unless bam_baits are present
    bams = [] #list of bams or (bam,bam_bait) tuples
    if bam_baits:
        if len(input) != len(bam_baits):
            raise BamMetricsError('number of bam files ({}) is not equal to number of corresponding capture bait interval list files ({})'.format(len(input), len(bam_baits)))
        for bam_ix in range(len(input)):
            unglobbed_bams = my.unglob(input[bam_ix])
            for unglobbed_bam in unglobbed_bams:
                bams.append ((unglobbed_bam,bam_baits[bam_ix]))
    else:
        bams = my.unglob(input)
    
    cmds = []
    hold_jids = []
    bam_metrics_files = []
    depthofcoverage_cmds = []

    for b in bams:
        if isinstance(b,basestring):
            bam = b
            bam_bait = bait
        else:
            bam, bam_bait = b
            if not bam_bait or not my.file_exists(bam_bait):
                bam_bait = bait
        if not bam or not my.file_exists(bam):
            raise BamMetricsError('bam file {} does not exist'.format(bam))
        #output files
        metricsbase = os.path.join(out_dir, output_prefix+os.path.basename(bam)+'.metrics')
        validate = metricsbase + '.validate'
        idxstats = metricsbase + '.idxstats'
        indexstats = metricsbase + '.indexstats'
        gcbiasmetrics = metricsbase + '.gcbias.table'
        gcbiasmetricschart = metricsbase + '.gcbias.pdf'
        gcbiasmetricssummary = metricsbase + '.gcbias.summary'
        librarycomplexity = metricsbase + '.librarycomplexity'
        coverage_output_prefix = os.path.join(out_dir, output_prefix+os.path.basename(bam))

        if not NoPicardMetrics:
            #idxstats is samtools, not picard and requires *.bam.bai for index name
            #cmds += ['{0} idxstats {1} > {2};'.format (samtools, bam, idxstats)]
            #bam_metrics_files += [idxstats]

            #BamIndexStats
            picard_args = 'INPUT={}'.format(bam)
            cmd = 'python {} -t {} --job_dir {} --picard_args {} --redirect_stdout {}'.format(
                picard.__file__, 'BamIndexStats', my_job_dir, picard_args, indexstats)
            cmds += [cmd]
            bam_metrics_files += [indexstats]
        
            #CalculateHsMetrics
            for target in targets:
                hsoutput = '{0}.target_{1}.hsmetrics'.format(metricsbase,os.path.basename(target))
                picard_args = 'INPUT={} OUTPUT={} BAIT_INTERVALS={} TARGET_INTERVALS={}'.format(
                    bam, hsoutput, bam_bait, target)
                cmd = 'python {} -t {} --job_dir {} --picard_args {}'.format(
                    picard.__file__, 'CalculateHsMetrics', my_job_dir, picard_args)
                cmds += [cmd]
                bam_metrics_files += [hsoutput]
        
            #CollectMultipleMetrics
            picard_args = 'INPUT={} OUTPUT={} REFERENCE_SEQUENCE={}'.format(
                bam, metricsbase, reference)
            multi_valued_arg_strings = 'PROGRAM='+' PROGRAM='.join((
                'CollectAlignmentSummaryMetrics','CollectInsertSizeMetrics','QualityScoreDistribution','MeanQualityByCycle'))
            cmd = 'python {} -t {} --job_dir {} --multi_valued_arg_strings {} --picard_args {}'.format(
                picard.__file__, 'CollectMultipleMetrics', my_job_dir, multi_valued_arg_strings, picard_args)
            cmds += [cmd]
            bam_metrics_files += [metricsbase+'.*_metrics']
        
            #CollectGcBiasMetrics
            picard_args = 'INPUT={} OUTPUT={} REFERENCE_SEQUENCE={} CHART_OUTPUT={} SUMMARY_OUTPUT={}'.format(
                bam, gcbiasmetrics, reference, gcbiasmetricschart, gcbiasmetricssummary)
            cmd = 'python {} -t {} --job_dir {} --picard_args {}'.format(
                picard.__file__, 'CollectGcBiasMetrics', my_job_dir, picard_args)
            cmds += [cmd]
            bam_metrics_files += [gcbiasmetrics]
            bam_metrics_files += [gcbiasmetricssummary]
        
            #EstimateLibraryComplexity READ_NAME_REGEX={3} OPTICAL_DUPLICATE_PIXEL_DISTANCE
            picard_args = 'INPUT={} OUTPUT={} READ_NAME_REGEX="\'{}\'" OPTICAL_DUPLICATE_PIXEL_DISTANCE={}'.format(
                bam, librarycomplexity, regex, pixels)
            cmd = 'python {} -t {} --job_dir {} --picard_args {}'.format(
                picard.__file__, 'EstimateLibraryComplexity', my_job_dir, picard_args)
            cmds += [cmd]
            bam_metrics_files += [librarycomplexity]

        #DepthOfCoverage
        if not NoDepthOfCoverage:
            for interval in depth_of_coverage_intervals:
                doc_prefix = os.path.join(dirs['depth_of_coverage'], '{}.{}.depth'.format(os.path.basename(bam), os.path.basename(interval)))
                params = '-I {} -o {} -R {} -L {} -geneList {} -dcov 1000 -l INFO -pt sample --omitDepthOutputAtEachBase --omitLocusTable -ct 1 -ct 2 -ct 3 -ct 4 -ct 5 -ct 6 -ct 7 -ct 8 -ct 9 -ct 10 -ct 11 -ct 12 -ct 13 -ct 14 -ct 15 -ct 16 -ct 17 -ct 18 -ct 19 -ct 20 -ct 21 -ct 22 -ct 23 -ct 24 -ct 25 -ct 26 -ct 27 -ct 28 -ct 29 -ct 30 -ct 40 -ct 50 -ct 60 -ct 70 -ct 80 -ct 90 -ct 100 -ct 110 -ct 120 -ct 130 -ct 140 -ct 150 -ct 160 -ct 170 -ct 180 -ct 190 -ct 200 -ct 300 -ct 400 -ct 500 -ct 600 -ct 700 -ct 800 -ct 900 -ct 1000'.format(
                    bam, doc_prefix, reference, interval, depth_of_coverage_gene_list)
                job = gatk(dirs, 'DepthOfCoverage', gatk_jar, params, job_name_prefix='depth_of_coverage', processors=8, java_mem='30g', memory='30G', queue='all.q@compute-[235]*')
                print 'scheduled GATK DepthOfCoverage job {}.\njob_args = {}'.format(job.jobId, ' '.join(job.args))
                doc_jobid = job.jobId
                hold_jids += [job.jobId]
    if len(cmds) > 0:
        job_name = my_name+'_picard'
        job_dir = my_job_dir
        job = my.run_job(cmds, job_name, job_dir, synchronous=False, processors=8)
        picard_metrics_jobid = job.jobId
        hold_jids += [job.jobId]
        print 'bam metrics job {}, commands\n{}'.format(job.jobId, '\n'.join(cmds))
    else:
        print 'bam metrics nothing to do'
    
    return
   

def main():

    #command line arguments
    parser = argparse.ArgumentParser(parents=[my.default_parser()],
        description = 'Calculate bam file metrics using picard and gatk DepthOfCoverage',
        epilog = '{} version {} {}'.format(name, version, copyright))
    parser.add_argument('--input', '-i', metavar='BAM', nargs='+',
        help='input list of bam files')
    parser.add_argument('--out_dir', '-o',
        help='output directory (default: dirs["bam_metrics"]')
    parser.add_argument('--output_prefix', nargs='?',
        help='prefix for output file names (default: bam_metrics_<timestamp>')
    parser.add_argument('--job_dir', '-j',
        help='parent job directory (default: dirs[jobs]')
    parser.add_argument('--reference', '-r',
        help='path to reference sequence fasta file (default: from config->reference.default_reference)')
    parser.add_argument('--gatk_jar',
        help='path to GenomeAnalysisTK jar')
    parser.add_argument('--bait','-B', metavar='BAIT_INTERVAL_LIST', nargs='?',
        help='bait interval list (default: config->interval.nimblegen_SeqCapEZ_Exome_v30)')
    parser.add_argument('--bam_baits', metavar='BAM_BAITS', nargs='*',
        help='list of bait interval list files corresponding to bam file globs')
    parser.add_argument('--targets','-T', metavar='TARGET_INTERVAL_LIST', nargs='+',
        help='list of target interval_list files (default: [bait, config->interval.ensembl_protein_coding_genes, config->interval.refgene_interval_list])')
    parser.add_argument('--depth_of_coverage_intervals', nargs='*',
        help='list of interval_list files for depth of coverage (default: targets)')
    parser.add_argument('--depth_of_coverage_gene_list',
        help='list of  gene list file for depth of coverage (default: config->interval.refgene_list)')
    parser.add_argument('--regex','-X', nargs='?',
        help='read name regex (default for illumina/novoalign: config->picard.READ_NAME_REGEX)')
    parser.add_argument('--pixels','-D', nargs='?', type=int, default=100,
        help='optical duplicate pixel distance (default for illumina: 100)')
    parser.add_argument('--samtools',
        help='path to samtools (default: config->samtools.samtools)')
    parser.add_argument('--addtional_metrics_to_consolidate', metavar='METRIC', nargs='*',
        help='additional files to include in consolidation (e.g., *.markdup.metrics')
    parser.add_argument('--NoPicardMetrics', action='store_true', default=False, help='do not calculate picard metrics')
    parser.add_argument('--NoDepthOfCoverage', action='store_true', default=False, help='do not calculate DepthOfCoverage')
    args = parser.parse_args()

    dirs = my.create_default_directory_structure(args.top_dir)
    config = my.get_config(args)

    run(top_dir=args.top_dir, input=args.input, out_dir=args.out_dir, output_prefix=args.output_prefix,
        job_dir=args.job_dir, reference=args.reference, depth_of_coverage_intervals=args.depth_of_coverage_intervals,
        depth_of_coverage_gene_list=args.depth_of_coverage_gene_list,
        bait=args.bait, bam_baits = args.bam_baits, targets=args.targets, regex=args.regex, pixels=args.pixels,
        gatk_jar=args.gatk_jar, samtools=args.samtools, addtional_metrics_to_consolidate=args.addtional_metrics_to_consolidate,
        NoPicardMetrics=args.NoPicardMetrics, NoDepthOfCoverage=args.NoDepthOfCoverage)

if __name__ == "__main__": sys.exit(main())



