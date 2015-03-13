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

class BamMetricsError(Exception): pass

#gmd
#--NoDepthOfCoverage -d /scratch0/tmp/myourshaw/gmd1 --output_prefix test -i /scratch0/tmp/myourshaw/gmd1/bams/readgroup_bams/*.novoalign.fixmate.bam

#1 small gmd
#--NoDepthOfCoverage -d /scratch0/tmp/myourshaw/gmd1 --output_prefix test -i /scratch0/tmp/myourshaw/gmd1/bams/readgroup_bams/HWI-ST430.243.6.ATCACGA.novoalign.fixmate.bam

#gmd35
#-d /scratch0/tmp/myourshaw/gmdmms -i /scratch0/tmp/myourshaw/gmd/bams/sample_bams/*.sample.recalibrated.bam /scratch0/tmp/myourshaw/mms/bams/sample_bams/*.sample.recalibrated.bam

def run(top_dir, input, out_dir=None, output_prefix=None, job_dir=None, ref=None, genelist=None,
        targetlimit=None, bait=None, targets=None, regex=r'"[^_]+_[^_]+_[0-9]+_([0-9]+)_([0-9.]+)_([0-9.]+).*"',
        pixels=100, gatk=None, samtools=None, consolidate=None, NoPicardMetrics=False,
        NoDepthOfCoverage=False):

    dirs = my.create_default_directory_structure(top_dir)
    config = my.get_config()
     
    my_name = 'bam_metrics'
    my_job_dir = my.unique_dir(my_name, job_dir) if job_dir else my.unique_dir(my_name, dirs['jobs'])
    
    if not out_dir:
        out_dir = os.path.join(dirs['metrics'],'bam_metrics')
    my.makedir(out_dir)

    if not output_prefix:
        output_prefix = 'bam_metrics_'+my.localtime_squish()
        
    #get defaults from configuration file
    if not ref:
        ref = config.get('reference','b37_fasta')
    if not genelist:
        genelist = config.get('gene','refgene_b37_sorted')
    if not bait:
        bait = config.get('interval','sureselect_50mb_interval_list')
    if not targetlimit:
        targetlimit = bait
    if not targets:
        targets = [config.get('interval','ensembl_protein_coding_genes'), config.get('interval','ensembl_all_genes'), config.get('interval','ensembl_all_genes_utr')]
    targets += [bait]
    if not gatk:
        gatk = config.get('gatk','GenomeAnalysisTK')
    if not samtools:
        samtools = config.get('samtools','samtools')

    bams = my.unglob(input)
    
    cmds = []
    hold_jids = []
    bam_metrics_files = []
    depthofcoverage_cmds = []

    for bam in bams:
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
                    bam, hsoutput, bait, target)
                cmd = 'python {} -t {} --job_dir {} --picard_args {}'.format(
                    picard.__file__, 'CalculateHsMetrics', my_job_dir, picard_args)
                cmds += [cmd]
                bam_metrics_files += [hsoutput]
        
            #CollectMultipleMetrics
            picard_args = 'INPUT={} OUTPUT={} REFERENCE_SEQUENCE={}'.format(
                bam, metricsbase, ref)
            multi_valued_arg_strings = 'PROGRAM='+' PROGRAM='.join((
                'CollectAlignmentSummaryMetrics','CollectInsertSizeMetrics','QualityScoreDistribution','MeanQualityByCycle'))
            cmd = 'python {} -t {} --job_dir {} --multi_valued_arg_strings {} --picard_args {}'.format(
                picard.__file__, 'CollectMultipleMetrics', my_job_dir, multi_valued_arg_strings, picard_args)
            cmds += [cmd]
            bam_metrics_files += [metricsbase+'.*_metrics']
        
            #CollectGcBiasMetrics
            picard_args = 'INPUT={} OUTPUT={} REFERENCE_SEQUENCE={} CHART_OUTPUT={} SUMMARY_OUTPUT={}'.format(
                bam, gcbiasmetrics, ref, gcbiasmetricschart, gcbiasmetricssummary)
            cmd = 'python {} -t {} --job_dir {} --picard_args {}'.format(
                picard.__file__, 'CollectGcBiasMetrics', my_job_dir, picard_args)
            cmds += [cmd]
            bam_metrics_files += [gcbiasmetrics]
            bam_metrics_files += [gcbiasmetricssummary]
        
            #EstimateLibraryComplexity READ_NAME_REGEX={3} OPTICAL_DUPLICATE_PIXEL_DISTANCE
            picard_args = 'INPUT={} OUTPUT={} READ_NAME_REGEX="{}" OPTICAL_DUPLICATE_PIXEL_DISTANCE={}'.format(
                bam, librarycomplexity, regex, pixels)
            cmd = 'python {} -t {} --job_dir {} --picard_args {}'.format(
                picard.__file__, 'EstimateLibraryComplexity', my_job_dir, picard_args)
            cmds += [cmd]
            bam_metrics_files += [librarycomplexity]

        #DepthOfCoverage
        if not NoDepthOfCoverage:
            cmd='\
$(which java) -Xmx16g -Djava.io.tmpdir={0} \
-jar {1} \
-T DepthOfCoverage \
-l INFO \
-R {2} \
-I {3} \
-o {4} \
-L {5} \
-geneList {6} \
--omitDepthOutputAtEachBase \
-dcov 1000 \
-pt sample \
-ct 1 -ct 2 -ct 3 -ct 4 -ct 5 -ct 6 -ct 7 -ct 8 -ct 9 -ct 10 -ct 20 -ct 30 -ct 40 -ct 50 -ct 60 -ct 70 -ct 80 -ct 90 -ct 100 -ct 110 -ct 120 -ct 130 -ct 140 -ct 150 -ct 160 -ct 170 -ct 180 -ct 190 -ct 200 -ct 300 \
;'.format(my.unique_dir('depthofcoverage_tmp', my_job_dir), gatk, ref, bam, coverage_output_prefix, targetlimit, genelist)
            depthofcoverage_cmds += [cmd]
            bam_metrics_files += [coverage_output_prefix+'*']
    print '\n'.join(cmds)
    if len(cmds) > 0:
        job_name = my_name+'_picard'
        job_dir = my_job_dir
        job = my.run_job(cmds, job_name, job_dir, synchronous=False, processors=8)
        picard_metrics_jobid = job.jobId
        hold_jids += [job.jobId]
    
    if len(depthofcoverage_cmds) > 0:
        job_name = my_name+'_depthofcoverage'
        job_dir = my_job_dir
        job = my.run_job(depthofcoverage_cmds, job_name, job_dir, synchronous=False, processors=8, memory = '24G')
        depthofcoverage_jobid = job.jobId
        hold_jids += [job.jobId]

    return
   

def main():

    #command line arguments
    parser = argparse.ArgumentParser(parents=[my.default_parser()],
        description = 'Calculate bam file metrics using picard and gatk DepthOfCoverage',
        epilog = 'pypeline.bam_metrics version 1.0β1 ©2011-2012 Michael Yourshaw all rights reserved')
    parser.add_argument('--input', '-i', metavar='BAM', nargs='+',
        help='input list of bam files')
    parser.add_argument('--out_dir', '-o',
        help='output directory (default: dirs["bam_metrics"]')
    parser.add_argument('--output_prefix',
        help='prefix for output file names (default: bam_metrics_timestamp')
    parser.add_argument('--job_dir', '-j',
        help='parent job directory (default: dirs[jobs]')
    parser.add_argument('--ref','-R', nargs='?',
        help='reference genome fasta (default: config->b37_fasta)')
    parser.add_argument('--genelist','-G', nargs='?',
        help='GATK DepthOfCoverage geneList (default: config->refgene_b37_sorted)')
    parser.add_argument('--targetlimit','-L', nargs='?',
        help='GATK DepthOfCoverage target limit (default: bait interval list)')
    parser.add_argument('--bait','-B', metavar='BAIT_INTERVAL_LIST', nargs='?',
        help='bait interval list (default: config->sureselect_50mb_interval_list)')
    parser.add_argument('--targets','-T', metavar='TARGET_INTERVAL_LIST', nargs='*',
        help='list of target interval_list files (default: <bait interval list> config->ensembl_all_genes config->ensembl_protein_coding_genes)')
    parser.add_argument('--regex','-X', nargs='?',default=r'"[^_]+_[^_]+_[0-9]+_([0-9]+)_([0-9.]+)_([0-9.]+).*"',
        help='read name regex (default for illumina/novoalign: "[^_]+_[^_]+_[0-9]+_([0-9]+)_([0-9.]+)_([0-9.]+).*")')
    parser.add_argument('--pixels','-D', nargs='?', type=int, default=100,
        help='optical duplicate pixel distance (default for illumina: 100)')
    parser.add_argument('--gatk',
        help='path to GATK (default: config->gatk.GenomeAnalysisTK)')
    parser.add_argument('--samtools',
        help='path to samtools (default: config->samtools.samtools)')
    parser.add_argument('--consolidate', metavar='METRIC', nargs='*',
        help='additional files to include in consolidation (e.g., *hahahaaha.markdup.metrics')
    parser.add_argument('--NoPicardMetrics', action='store_true', default=False, help='do not calculate picard metrics')
    parser.add_argument('--NoDepthOfCoverage', action='store_true', default=False, help='do not calculate DepthOfCoverage')
    args = parser.parse_args()

    dirs = my.create_default_directory_structure(args.top_dir)
    config = my.get_config(args)

    run(top_dir=args.top_dir, input=args.input, out_dir=args.out_dir, output_prefix=args.output_prefix,
        job_dir=args.job_dir, ref=args.ref, genelist=args.genelist, targetlimit=args.targetlimit,
        bait=args.bait, targets=args.targets, regex=args.regex, pixels=args.pixels,
        gatk=args.gatk, samtools=args.samtools, consolidate=args.consolidate,
        NoPicardMetrics=args.NoPicardMetrics, NoDepthOfCoverage=args.NoDepthOfCoverage)

if __name__ == "__main__": sys.exit(main())



