#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import stat
import argparse
from ConfigParser import SafeConfigParser #configparser in python 3
import glob
import tempfile
from time import localtime, strftime
import re
from warnings import warn
import shutil
import my
import job
import picard
import validate_RealignerTargetCreator_intervals
import validate_recal_csv
import validate_vcf
import verify_sam_file_validation
import verify_OK_file

class Bam2VcfError(Exception): pass

#--NoUnifiedGenotyper --NoVariantRecalibrator --NoAnnotation -d /scratch0/tmp/myourshaw/mms --bams /scratch0/tmp/myourshaw/mms/bams/sample_bams/*.sample.bam --vcf_prefix /scratch0/tmp/myourshaw/mms/vcfs/mms4
#--NoRealignRecalibrate -d /scratch0/tmp/myourshaw/gmd --bams /scratch0/tmp/myourshaw/gmd/bams/sample_bams/*.sample.bam --vcf_prefix /scratch0/tmp/myourshaw/gmd/vcfs/gmd32
#-d /scratch0/tmp/myourshaw/gmd --bams /scratch0/tmp/myourshaw/gmd/bams/sample_bams/gmd19.sample.bam --vcf_prefix /scratch0/tmp/myourshaw/gmd/vcfs/gmd32
#--NoUnifiedGenotyper --NoVariantRecalibrator --NoAnnotation -d /scratch0/tmp/myourshaw/gmd --bams /scratch0/tmp/myourshaw/gmd/bams/sample_bams/GMD145.sample.bam /scratch0/tmp/myourshaw/gmd/bams/sample_bams/GMD149.sample.bam /scratch0/tmp/myourshaw/gmd/bams/sample_bams/GMD165.sample.bam --vcf_prefix /scratch0/tmp/myourshaw/gmd/vcfs/gmd3
#--NoRealignRecalibrate --NoUnifiedGenotyper --NoVariantRecalibrator -d /scratch0/tmp/myourshaw/gmd --bams /scratch0/tmp/myourshaw/gmd/bams/sample_bams/*.sample.bam --vcf_prefix /scratch0/tmp/myourshaw/gmd/vcfs/gmd32

#gmd35mms
#--NoRealignRecalibrate -d /scratch0/tmp/myourshaw/gmdmms --bams /scratch0/tmp/myourshaw/gmd/bams/sample_bams/*.sample.bam /scratch0/tmp/myourshaw/mms/bams/sample_bams/*.sample.bam --vcf_prefix /scratch0/tmp/myourshaw/gmdmms/vcfs/gmd35mms

#gmd35MMS
#--NoRealignRecalibrate -d /scratch0/tmp/myourshaw/gmdMMS --bams /scratch0/tmp/myourshaw/gmd/bams/sample_bams/*.sample.bam /scratch0/tmp/myourshaw/mms/bams/sample_bams/*.sample.bam --vcf_prefix /scratch0/tmp/myourshaw/gmdMMS/vcfs/gmd35mms

#gmd28mms
#--NoRealignRecalibrate -d /scratch0/tmp/myourshaw/gmd28mms --bams /scratch0/tmp/myourshaw/gmd/bams/sample_bams/GMD*.sample.bam /scratch0/tmp/myourshaw/mms/bams/sample_bams/*.sample.bam --vcf_prefix /scratch0/tmp/myourshaw/gmd28mms/vcfs/gmd28mms

#gmd28mms all loci
#--NoRealignRecalibrate --output_mode EMIT_ALL_CONFIDENT_SITES -d /scratch0/tmp/myourshaw/gmd28mms_all --bams /scratch0/tmp/myourshaw/gmd/bams/sample_bams/GMD*.sample.bam /scratch0/tmp/myourshaw/mms/bams/sample_bams/*.sample.bam --vcf_prefix /scratch0/tmp/myourshaw/gmd28mms_all/vcfs/gmd28mms_all

#run all of bipolar bams with GMD controls, but do RealignRecalibrate step only on bipolar samples
#-d /scratch0/tmp/aliz/bipolar/data/truseqABCD_myourshaw --bams /scratch0/tmp/aliz/bipolar/data/truseqABCD_myourshaw/bams/sample_bams/*.sample.bam --recalibrated_bams /scratch0/tmp/myourshaw/gmd/bams/sample_bams/zila_temp/GMD*.sample.recalibrated.bam --vcf_prefix /scratch0/tmp/aliz/bipolar/data/truseqABCD_myourshaw/vcfs/bp20

#only do RealignRecalibrate on 3 samples that need redoing, then use all samples for UnifiedGenotyper and beyond
#-d /scratch0/tmp/aliz/ocd_tic/data/pipeline1 --bams /scratch0/tmp/myourshaw/gmd/bams/sample_bams/zila_temp/GMD72.sample.bam /scratch0/tmp/myourshaw/gmd/bams/sample_bams/zila_temp/GMD108.sample.bam /scratch0/tmp/myourshaw/gmd/bams/sample_bams/zila_temp/GMD165.sample.bam --recalibrated_bams /scratch0/tmp/myourshaw/gmd/bams/sample_bams/zila_temp/GMD*.sample.recalibrated.bam /scratch0/tmp/aliz/ocd_tic/data/pipeline1/bams/sample_bams/*.sample.recalibrated.bam --vcf_prefix /scratch0/tmp/aliz/ocd_tic/data/pipeline1/vcfs/ocd18

#provide list of samples to run VariantEval selectively on those only, as well
#--NoUnifiedGenotyper -d /scratch0/tmp/aliz/bipolar/data/truseqABCD_myourshaw --vcf_prefix /scratch0/tmp/aliz/bipolar/data/truseqABCD_myourshaw/vcfs/bp20 --samples "8000350017 8000350018 8000350019 8000350031 8000350041 8000350042 8000350043 8000350053 8000350054 8000350055 8000350065 8000350066 8000350077 8000350078 8000350080 8000350089 8000350090 8000350101 8000350102 8000368254"

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
    print cmd
    return job
    
def run(top_dir, bams, recalibrated_bams, vcf_prefix, parent_job_dir=None, targets=None, gatk_jar=None,
        analyzecovariates_jar=None, ref=None, known_indels_vcfs=None, dbsnp=None,
        hapmap=None, omni=None, rscript=None, gatkrscripts=None, output_mode='EMIT_VARIANTS_ONLY',
        email=None, synchronous=False, hold_jids=[], NoRealignRecalibrate=False,
        NoUnifiedGenotyper=False, NoVariantRecalibrator=False, SnpEffAnnotation=False,
        NoVariantAnnotator=False, args=None, snpEff_jar=None,
        snpEff_database='GRCh37.64', snpEff_config=None, additional_VQSR_vcfs=None, samples=None):

    dirs = my.create_default_directory_structure(top_dir)
    if not parent_job_dir or not os.path.isdir(parent_job_dir):
        parent_job_dir = dirs['jobs']
    config = my.get_config()
    python = config.get('DEFAULT','python')

    analyze_covariates_files = ('.CycleCovariate.dat','.CycleCovariate.dat.Cycle_hist.pdf','.CycleCovariate.dat.qual_diff_v_Cycle.pdf','.CycleCovariate.dat.reported_qual_v_Cycle.pdf','.DinucCovariate.dat','.DinucCovariate.dat.Dinuc_hist.pdf','.DinucCovariate.dat.qual_diff_v_Dinuc.pdf','.DinucCovariate.dat.reported_qual_v_Dinuc.pdf','.QualityScoreCovariate.dat','.QualityScoreCovariate.dat.quality_emp_hist.pdf','.QualityScoreCovariate.dat.quality_emp_v_stated.pdf','.QualityScoreCovariate.dat.quality_rep_hist.pdf')
    #return value: a list of job ids that can be used by downstream hold_jid
    job_ids = []
    
    if not hold_jids:
        hold_jids = []
    if not targets:
        targets = config.get('interval','ensembl_all_genes_utr')
    if not my.file_exists(targets):
        raise Bam2VcfError('cannot find target intervals')
    if not gatk_jar:
        gatk_jar = config.get('gatk','GenomeAnalysisTK')
    if not my.file_exists(gatk_jar):
        raise Bam2VcfError('cannot find GenomeAnalysisTK jar')
    if not analyzecovariates_jar:
        analyzecovariates_jar = config.get('gatk','AnalyzeCovariates')
    if not my.file_exists(analyzecovariates_jar):
        raise Bam2VcfError('cannot find AnalyzeCovariates jar')
    if not ref:
        ref = config.get('gatk','human_g1k_v37')
    if not my.file_exists(ref):
        raise Bam2VcfError('cannot find reference genome')
    if not known_indels_vcfs:
        known_indels_vcfs = [config.get('gatk','mills_devine_indels'), config.get('gatk','low_coverage_vqsr_indels')]
    if not isinstance(known_indels_vcfs, (list, tuple)):
        known_indels_vcfs = [known_indels_vcfs]
    if not dbsnp:
        dbsnp = config.get('variant','dbsnp')
    if not my.file_exists(dbsnp):
        raise Bam2VcfError('cannot find dbsnp')
    if not hapmap:
        hapmap = config.get('gatk','hapmap')
    if not my.file_exists(hapmap):
        raise Bam2VcfError('cannot find hapmap')
    if not omni:
        omni = config.get('gatk','omni')
    if not my.file_exists(hapmap):
        raise Bam2VcfError('cannot find hapmap')
    if not rscript:
        rscript = config.get('DEFAULT','Rscript')
    if not my.file_exists(rscript):
        raise Bam2VcfError('cannot find Rscript')
    if not gatkrscripts:
        gatkrscripts = config.get('gatk','gatk_R_scripts')
    if not os.path.isdir(gatkrscripts):
        raise Bam2VcfError('cannot find gatk R script resources directory')
    if SnpEffAnnotation:
        if not snpEff_jar:
            snpEff_jar = config.get('gatk','snpEff_jar')
        if not my.file_exists(snpEff_jar):
            raise Bam2VcfError('cannot find snpEff.jar')
        if not snpEff_config:
            snpEff_config = os.path.join(os.path.dirname(snpEff_jar),'snpEff.config')
        if not my.file_exists(snpEff_config):
            raise Bam2VcfError('cannot find snpEff.config')
        
    job_name_prefix = 'bams2vcf_{}'.format(my.localtime_squish())
    if not os.path.isdir(parent_job_dir):
        parent_job_dir = dirs['job']
    job_dir = my.unique_dir(job_name_prefix, parent_job_dir)
    
    VerifyValidateCountCovariates2_jobid = None
    verifyvalidaterecalibratedbamfile_jobids = []
    neededforunifiedgenotyper_jobids = []
    recalibrated_bams = my.unglob(recalibrated_bams) if recalibrated_bams else []
    bams = my.unglob(bams)
    
    if recalibrated_bams:
        for bam in bams:
            if bam in recalibrated_bams:
                raise Bam2VcfError('bam file (' + bam + ') listed under both bams and recalibrated_bams parameters')
    ug_holds =  list(hold_jids)
    
    for bam in bams:
        bam_header = my.bam_peek(bam)
        rg_ids = [rg['ID'] for rg in bam_header['RG']]
        bam_holds = list(hold_jids)
        #files and directories
        bam_base = my.bam_strip(bam)
        RealignerTargetCreator_intervals = bam_base+'.realignertargetcreator.intervals'
        RealignerTargetCreator_validate = RealignerTargetCreator_intervals+'.validate'
        realigned_bam = bam_base+'.realigned.bam'
        realigned_bai = bam_base+'.realigned.bai'
        realigned_bam_validate = realigned_bam+'.validate'
        recal_csv = bam_base+'.recal.csv'
        recal_csv_validate = recal_csv+'.validate'
        analyze_covariates1_dir = os.path.join(dirs['quality_score_recalibration_metrics'], os.path.basename(bam_base), 'pre-recalibration')
        my.makedir(analyze_covariates1_dir)
        analyze_covariates1_files = [os.path.join(analyze_covariates1_dir,i+f) for f in analyze_covariates_files for i in rg_ids]
        recalibrated_bam = bam_base+'.recalibrated.bam'
        recalibrated_bai = bam_base+'.recalibrated.bai'
        recalibrated_bam_validate = recalibrated_bam+'.validate'
        recal2_csv = bam_base+'.recal2.csv'
        recal2_csv_validate = recal2_csv+'.validate'
        analyze_covariates2_dir = os.path.join(dirs['quality_score_recalibration_metrics'], os.path.basename(bam_base), 'post-recalibration')
        my.makedir(analyze_covariates2_dir)
        analyze_covariates2_files = [os.path.join(analyze_covariates2_dir,i+f) for f in analyze_covariates_files for i in rg_ids]
        bam_job_name = job_name_prefix+'_'+os.path.basename(bam_base)
        
        if recalibrated_bam not in recalibrated_bams:
            recalibrated_bams += [recalibrated_bam]
        else:
            warn('recalibrated bam file (' + recalibrated_bam + ') already exists and will be overwritten')
            
        if not NoRealignRecalibrate:
            
            #RealignerTargetCreator
            if args.overwrite_all_output_files or args.overwrite_RealignerTargetCreator_intervals or not my.check_files(RealignerTargetCreator_intervals, [(RealignerTargetCreator_validate,'OK')], [(RealignerTargetCreator_intervals, RealignerTargetCreator_validate)]):
                params = '-R "{}" -o {} -I "{}" --known "{}"'.format(
                     ref, RealignerTargetCreator_intervals, bam, '" --known "'.join(known_indels_vcfs))
                job = gatk(dirs, 'RealignerTargetCreator', gatk_jar, params, job_name_prefix=bam_job_name, hold_jid=bam_holds)
                RealignerTargetCreator_jobid = job.jobId
                job_ids += [job.jobId]
                bam_holds += [job.jobId]
                
                #validate format of interval file
                job_name = 'ValidateRealignerTargetCreatorIntervalFile_'+bam_job_name
                cmd = '{} {} --input {} --output {}'.format(python, validate_RealignerTargetCreator_intervals.__file__, RealignerTargetCreator_intervals, RealignerTargetCreator_intervals+'.validate')
                job = my.run_job(cmd, job_name, job_dir, hold_jid=bam_holds[-1] if len(bam_holds)>0 else None)
                ValidateRealignerTargetCreator_jobid = job.jobId
                job_ids += [job.jobId]
                bam_holds += [job.jobId]
                
                #verify validation
                job_name = 'VerifyValidateRealignerTargetCreatorIntervalFile_'+bam_job_name
                cmd = '{} {} --input {}'.format(python, verify_OK_file.__file__, RealignerTargetCreator_intervals+'.validate')
                job = my.run_job(cmd, job_name, job_dir, hold_jid=bam_holds[-1] if len(bam_holds)>0 else None)
                VerifyValidateRealignerTargetCreator_jobid = job.jobId
                job_ids += [job.jobId]
                bam_holds += [job.jobId]
                neededforunifiedgenotyper_jobids += [job.jobId]
                    
            #IndelRealigner
            if args.overwrite_all_output_files or args.overwrite_realigned_bams or not my.check_files(realigned_bam, [(realigned_bam_validate,'No errors found\n')], [(realigned_bam, realigned_bai), (realigned_bam, realigned_bam_validate)]):
                params = '-I "{}" -R "{}" -targetIntervals "{}" -o "{}" -known "{}" -model USE_READS'.format(
                    bam, ref, RealignerTargetCreator_intervals, realigned_bam, '" -known "'.join(known_indels_vcfs))
                #job = gatk(dirs, 'IndelRealigner', gatk_jar, params, job_name_prefix=bam_job_name, hold_jid=RealignerTargetCreator_jobid)
                job = gatk(dirs, 'IndelRealigner', gatk_jar, params, job_name_prefix=bam_job_name, hold_jid=bam_holds[-1] if len(bam_holds)>0 else None)
                IndelRealigner_jobid = job.jobId
                job_ids += [job.jobId]
                bam_holds += [job.jobId]
                
                #validate bam file with picard
                picard_args = {'INPUT': realigned_bam, 'OUTPUT': realigned_bam_validate}
                job = picard.run(tool='ValidateSamFile', picard_args=picard_args, job_dir=job_dir, reference=ref, hold_jid=bam_holds[-1] if len(bam_holds)>0 else None)
                ValidateRealignedBamFile_jobid = job.jobId
                job_ids += [job.jobId]
                bam_holds += [job.jobId]
                
                #verify validation
                job_name = 'VerifyValidateRealignedBamFile_'+bam_job_name
                cmd = '{} {} -v {}'.format(python, verify_sam_file_validation.__file__, realigned_bam_validate)
                job = my.run_job(cmd, job_name, job_dir, hold_jid=bam_holds[-1] if len(bam_holds)>0 else None)
                VerifyValidateRealignedBamFile_jobid = job.jobId
                job_ids += [job.jobId]
                bam_holds += [job.jobId]
                neededforunifiedgenotyper_jobids += [job.jobId]
                
            #CountCovariates1
            #SOLiD needs  --solid_nocall_strategy PURGE_READ
            if args.overwrite_all_output_files or args.overwrite_recal_csv or not my.check_files(recal_csv, [(recal_csv_validate,'OK')], [(recal_csv, recal_csv_validate)]):
                params = '-R "{}" -knownSites "{}" -I "{}" -recalFile "{}" -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate -nt 8 --solid_nocall_strategy PURGE_READ'.format(
                    ref, dbsnp, realigned_bam, recal_csv)
                job = gatk(dirs, 'CountCovariates', gatk_jar, params, job_name_prefix=bam_job_name, processors=8, hold_jid=bam_holds[-1] if len(bam_holds)>0 else None)
                CountCovariates1_jobid = job.jobId
                job_ids += [job.jobId]
                bam_holds += [job.jobId]
                
                #validate format of recal_csv file
                job_name = 'ValidateCountCovariates1File_'+bam_job_name
                cmd = '{} {} --input {} --output {}'.format(python, validate_recal_csv.__file__, recal_csv, recal_csv_validate)
                job = my.run_job(cmd, job_name, job_dir, hold_jid=bam_holds[-1] if len(bam_holds)>0 else None)
                ValidateCountCovariates1_jobid = job.jobId
                job_ids += [job.jobId]
                
                #verify validation
                job_name = 'VerifyValidateCountCovariates1File_'+bam_job_name
                cmd = '{} {} --input {}'.format(python, verify_OK_file.__file__, recal_csv_validate)
                job = my.run_job(cmd, job_name, job_dir, hold_jid=bam_holds[-1] if len(bam_holds)>0 else None)
                VerifyValidateCountCovariates1_jobid = job.jobId
                job_ids += [job.jobId]
                bam_holds += [job.jobId]
                
            #AnalyzeCovariates1
            if args.overwrite_all_output_files or args.overwrite_analyze_covariates_1 or not my.check_files(analyze_covariates1_files):
                params = '-recalFile "{}" -outputDir "{}" -ignoreQ 5'.format(
                    recal_csv, analyze_covariates1_dir)
                job = gatk(dirs, None, analyzecovariates_jar, params, job_name_prefix=bam_job_name, hold_jid=bam_holds[-1] if len(bam_holds)>0 else None)
                AnalyzeCovariates1_jobid = job.jobId
                job_ids += [job.jobId]
            
            #TableRecalibration
            #SOLiD needs  --solid_nocall_strategy PURGE_READ
            if args.overwrite_all_output_files or args.overwrite_recalibrated_bams or not my.check_files(recalibrated_bam, [(recalibrated_bam_validate,'No errors found\n')], [(recalibrated_bam, recalibrated_bai), (recalibrated_bam, recalibrated_bam_validate)]):
                params = '-R "{}" -I "{}" -o "{}" -recalFile "{}" -baq RECALCULATE --solid_nocall_strategy PURGE_READ'.format(
                    ref, realigned_bam, recalibrated_bam, recal_csv)
                job = gatk(dirs, 'TableRecalibration', gatk_jar, params, job_name_prefix=bam_job_name, hold_jid=bam_holds[-1] if len(bam_holds)>0 else None)
                TableRecalibration_jobid = job.jobId
                job_ids += [job.jobId]
                bam_holds += [job.jobId]
                
                #validate bam file with picard
                picard_args = {'INPUT': recalibrated_bam, 'OUTPUT': recalibrated_bam_validate}
                job = picard.run(tool='ValidateSamFile', picard_args=picard_args, job_dir=job_dir, reference=ref, hold_jid=bam_holds[-1] if len(bam_holds)>0 else None)
                ValidateRecalibratedBamFile_jobid = job.jobId
                job_ids += [job.jobId]
                bam_holds += [job.jobId]
                
                #verify validation
                job_name = 'VerifyValidateRecalibratedBamFile_'+bam_job_name
                cmd = '{} {} -v {}'.format(python, verify_sam_file_validation.__file__, recalibrated_bam_validate)
                job = my.run_job(cmd, job_name, job_dir, hold_jid=bam_holds[-1] if len(bam_holds)>0 else None, email=args.email)
                VerifyValidateRecalibratedBamFile_jobid = job.jobId
                job_ids += [job.jobId]
                bam_holds += [job.jobId]
                verifyvalidaterecalibratedbamfile_jobids += [job.jobId]
                neededforunifiedgenotyper_jobids += [job.jobId]
                
            #CountCovariates2
            #SOLiD needs  --solid_nocall_strategy PURGE_READ
            if args.overwrite_all_output_files or args.overwrite_recal2_csv or not my.check_files(recal2_csv, [(recal2_csv_validate,'OK')], [(recal2_csv, recal2_csv_validate)]):
                params = '-R "{}" -knownSites "{}" -I "{}" -recalFile "{}" -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate -nt 24 --solid_nocall_strategy PURGE_READ'.format(
                    ref, dbsnp, recalibrated_bam, recal2_csv)
                job = gatk(dirs, 'CountCovariates', gatk_jar, params, job_name_prefix=bam_job_name, processors=24, queue='all.q@compute-6*', hold_jid=bam_holds[-1] if len(bam_holds)>0 else None)
                CountCovariates2_jobid = job.jobId
                #bam_holds += [job.jobId]
                job_ids += [job.jobId]
                
                #validate format of recal_csv file
                job_name = 'ValidateCountCovariates2File_'+bam_job_name
                cmd = '{} {} --input {} --output {}'.format(python, validate_recal_csv.__file__, recal2_csv, recal2_csv_validate)
                job = my.run_job(cmd, job_name, job_dir, hold_jid=CountCovariates2_jobid)
                ValidateCountCovariates2_jobid = job.jobId
                #bam_holds += [job.jobId]
                job_ids += [job.jobId]
                
                #verify validation
                job_name = 'VerifyValidateCountCovariates2File_'+bam_job_name
                cmd = '{} {} --input {}'.format(python, verify_OK_file.__file__, recal2_csv_validate)
                job = my.run_job(cmd, job_name, job_dir, hold_jid=ValidateCountCovariates2_jobid)
                VerifyValidateCountCovariates2_jobid = job.jobId
                #bam_holds += [job.jobId]
                job_ids += [job.jobId]
                
            #AnalyzeCovariates2
            if args.overwrite_all_output_files or args.overwrite_analyze_covariates_2 or not my.check_files(analyze_covariates2_files):
                params = '-recalFile "{}" -outputDir "{}" -ignoreQ 5'.format(
                    recal2_csv, analyze_covariates2_dir)
                job = gatk(dirs, None, analyzecovariates_jar, params, job_name_prefix=bam_job_name, hold_jid=VerifyValidateCountCovariates2_jobid)
                AnalyzeCovariates2_jobid = job.jobId
                job_ids += [job.jobId]
        ug_holds += bam_holds
        
    #files and directories
    ug_raw_vcf = vcf_prefix+'.ug.raw.vcf'
    ug_raw_vcf_validate = ug_raw_vcf+'.validate'
    unifiedgenotyper_metrics = os.path.join(dirs['unified_genotyper_metrics'], job_name_prefix+'ug.raw.metrics')
    ug_raw_snvs_vcf = vcf_prefix+'.ug.raw.snvs.vcf'
    ug_raw_snvs_vcf_validate = ug_raw_snvs_vcf+'.validate'
    ug_raw_indels_vcf = vcf_prefix+'.ug.raw.indels.vcf'
    ug_raw_indels_vcf_validate = ug_raw_indels_vcf+'.validate'
    recal = ug_raw_snvs_vcf+'.recal'
    tranches = ug_raw_snvs_vcf+'.tranches'
    plots = ug_raw_snvs_vcf+'.plots.R'
    recalibrated_snvs_vcf = vcf_prefix+'.recalibrated.snvs.vcf'
    recalibrated_snvs_vcf_validate = recalibrated_snvs_vcf+'.validate'
    filtered_indels_vcf = vcf_prefix+'.filtered.indels.vcf'
    filtered_indels_vcf_validate = filtered_indels_vcf+'.validate'
    analysis_ready_vcf = vcf_prefix+'.analysis_ready.vcf'
    analysis_ready_vcf_validate = analysis_ready_vcf+'.validate'
    snpEff_vcf = re.sub(r'\.vcf$','',analysis_ready_vcf)+'.snpEff.vcf'
    snpEff_vcf_validate = snpEff_vcf+'.validate'
    #TODO:specify snpEff_genes path in command
    snpEff_genes = re.sub(r'\.vcf$','',analysis_ready_vcf)+'.snpEff_genes.txt'
    snpEff_summary = re.sub(r'\.vcf$','',analysis_ready_vcf)+'.snpEff_summary.html'
    annotated_vcf = vcf_prefix+'.annotated.vcf'
    annotated_vcf_validate = annotated_vcf+'.validate'
    vcf_job_name = os.path.basename(vcf_prefix)
    
    verifyvalidatevcf_jobids = []
    neededforselectvariants_jobids = []

    if not NoUnifiedGenotyper:
        #UnifiedGenotyper
        VerifyValidateUnifiedGenotyperVcfFile_jobid = None
        if args.overwrite_all_output_files or args.overwrite_unified_genotyper or not my.check_files(ug_raw_vcf_validate, [(ug_raw_vcf_validate,'OK')], [(ug_raw_vcf, ug_raw_vcf_validate)]):
            params = '-R "{}" -I "{}" --dbsnp "{}" -o "{}" -L "{}" -metrics "{}" --output_mode {} -glm BOTH -stand_call_conf 30.0 -stand_emit_conf 10.0 -dcov 2500 -G Standard -nt 16 -baq CALCULATE_AS_NECESSARY'.format(
                ref, '" -I "'.join(recalibrated_bams), dbsnp, ug_raw_vcf, targets, unifiedgenotyper_metrics, output_mode)
            #job = gatk(dirs, 'UnifiedGenotyper', gatk_jar, params, job_name_prefix=bam_job_name)
            #TODO: chage hold_jid below
            job = gatk(dirs, 'UnifiedGenotyper', gatk_jar, params, job_name_prefix=vcf_job_name, processors=8, hold_jid=ug_holds)
            UnifiedGenotyper_jobid = job.jobId
            job_ids += [job.jobId]
            ug_holds += [job.jobId]
            
            #validate format of vcf file
            job_name = 'ValidateUnifiedGenotyperVcfFile_'+vcf_job_name
            cmd = '{} {} --input {} --output {}'.format(python, validate_vcf.__file__, ug_raw_vcf, ug_raw_vcf_validate)
            job = my.run_job(cmd, job_name, job_dir, hold_jid=ug_holds)
            ValidateUnifiedGenotyperVcfFile_jobid = job.jobId
            job_ids += [job.jobId]
            ug_holds += [job.jobId]
            
            #verify validation
            job_name = 'VerifyValidateUnifiedGenotyperVcfFile_'+vcf_job_name
            cmd = '{} {} --input {}'.format(python, verify_OK_file.__file__, ug_raw_vcf_validate)
            job = my.run_job(cmd, job_name, job_dir, hold_jid=ug_holds)
            VerifyValidateUnifiedGenotyperVcfFile_jobid = job.jobId
            job_ids += [job.jobId]
            ug_holds += [job.jobId]
            verifyvalidatevcf_jobids += [job.jobId]
            neededforselectvariants_jobids += [job.jobId]
            
        #VariantEval
        params = '-R "{}" --eval {} --dbsnp "{}" -o "{}" -ST Sample'.format(
            ref, ug_raw_vcf, dbsnp, ug_raw_vcf+'.VariantEval')
        job = gatk(dirs, 'VariantEval', gatk_jar, params, job_name_prefix=vcf_job_name, hold_jid=ug_holds)
        VariantEvalUGRaw_jobid = job.jobId
        job_ids += [job.jobId]
        
    if not NoVariantRecalibrator:
        
        #merge additional_VQSR_vcfs into ug_raw_vcf
        if additional_VQSR_vcfs:
            pass
            
        snvs_holds = list(ug_holds)
        indels_holds = list(ug_holds)
        #SelectVariants for SNVs
        VerifyValidateSelectVariants_snvs_jobid = None
        if args.overwrite_all_output_files or args.overwrite_select_variants_snvs or not my.check_files(ug_raw_snvs_vcf_validate, [(ug_raw_snvs_vcf_validate,'OK')], [(ug_raw_snvs_vcf, ug_raw_snvs_vcf_validate)]):
            params = '-R "{}" --variant "{}" -o "{}" -selectType SNP'.format(
                ref, ug_raw_vcf, ug_raw_snvs_vcf)
            job = gatk(dirs, 'SelectVariants', gatk_jar, params, job_name_prefix=vcf_job_name, hold_jid=snvs_holds)
            SelectVariants_snvs_jobid = job.jobId
            job_ids += [job.jobId]
            ug_holds += [job.jobId]
            snvs_holds += [job.jobId]
            
            #validate format of vcf file
            job_name = 'ValidateSelectVariantsSnvsVcfFile_'+vcf_job_name
            cmd = '{} {} --input {} --output {}'.format(python, validate_vcf.__file__, ug_raw_snvs_vcf, ug_raw_snvs_vcf_validate)
            job = my.run_job(cmd, job_name, job_dir, hold_jid=snvs_holds)
            ValidateSelectVariants_snvs_jobid = job.jobId
            job_ids += [job.jobId]
            ug_holds += [job.jobId]
            snvs_holds += [job.jobId]
            
            #verify validation
            job_name = 'VerifyValidateSelectVariantsSnvsVcfFile_'+vcf_job_name
            cmd = '{} {} --input {}'.format(python, verify_OK_file.__file__, ug_raw_snvs_vcf_validate)
            job = my.run_job(cmd, job_name, job_dir, hold_jid=snvs_holds)
            VerifyValidateSelectVariants_snvs_jobid = job.jobId
            job_ids += [job.jobId]
            ug_holds += [job.jobId]
            snvs_holds += [job.jobId]
            verifyvalidatevcf_jobids += [job.jobId]
            
        #SelectVariants for indels
        VerifyValidateSelectVariants_indels_jobid = None
        if args.overwrite_all_output_files or args.overwrite_select_variants_indels or not my.check_files(ug_raw_indels_vcf_validate, [(ug_raw_indels_vcf_validate,'OK')], [(ug_raw_indels_vcf, ug_raw_indels_vcf_validate)]):
            params = '-R "{}" --variant "{}" -o "{}" -selectType INDEL'.format(
                ref, ug_raw_vcf, ug_raw_indels_vcf)
            job = gatk(dirs, 'SelectVariants', gatk_jar, params, job_name_prefix=vcf_job_name, hold_jid=indels_holds)
            SelectVariants_indels_jobid = job.jobId
            job_ids += [job.jobId]
            ug_holds += [job.jobId]
            indels_holds += [job.jobId]
            
            #validate format of vcf file
            job_name = 'ValidateSelectVariantsIndelsVcfFile_'+vcf_job_name
            cmd = '{} {} --input {} --output {}'.format(python, validate_vcf.__file__, ug_raw_indels_vcf, ug_raw_indels_vcf_validate)
            job = my.run_job(cmd, job_name, job_dir, hold_jid=indels_holds)
            ValidateSelectVariants_indels_jobid = job.jobId
            job_ids += [job.jobId]
            ug_holds += [job.jobId]
            indels_holds += [job.jobId]
            
            #verify validation
            job_name = 'VerifyValidateSelectVariantsIndelsVcfFile_'+vcf_job_name
            cmd = '{} {} --input {}'.format(python, verify_OK_file.__file__, ug_raw_indels_vcf_validate)
            job = my.run_job(cmd, job_name, job_dir, hold_jid=indels_holds)
            VerifyValidateSelectVariants_indels_jobid = job.jobId
            job_ids += [job.jobId]
            ug_holds += [job.jobId]
            indels_holds += [job.jobId]
            verifyvalidatevcf_jobids += [job.jobId]
         
        #VariantRecalibrator
        params = '-R "{}" -input "{}" --maxGaussians 6 -resource:hapmap,VCF,known=false,training=true,truth=true,prior=15.0 "{}" -resource:omni,VCF,known=false,training=true,truth=false,prior=12.0 "{}" -resource:dbsnp,VCF,known=true,training=false,truth=false,prior=8.0 "{}" -an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an InbreedingCoeff -mode SNP -recalFile "{}" -tranchesFile "{}" -rscriptFile "{}"'.format(
            ref, ug_raw_snvs_vcf, hapmap, omni, dbsnp, recal, tranches, plots)
        #job = gatk(dirs, 'VariantRecalibrator', gatk_jar, params, job_name_prefix=bam_job_name)
        job = gatk(dirs, 'VariantRecalibrator', gatk_jar, params, job_name_prefix=vcf_job_name, hold_jid=snvs_holds)
        VariantRecalibrator_jobid = job.jobId
        job_ids += [job.jobId]
        ug_holds += [job.jobId]
        
        #ApplyRecalibration
        VerifyValidateApplyRecalibration_jobid = None
        if args.overwrite_all_output_files or args.overwrite_recalibrated_snvs or not my.check_files(recalibrated_snvs_vcf_validate, [(recalibrated_snvs_vcf_validate,'OK')], [(recalibrated_snvs_vcf, recalibrated_snvs_vcf_validate)]):
            params = '-R "{}" -input "{}" --ts_filter_level 99.0 -tranchesFile "{}" -recalFile "{}" -o "{}" --mode SNP'.format(
                ref, ug_raw_snvs_vcf, tranches, recal, recalibrated_snvs_vcf)
            job = gatk(dirs, 'ApplyRecalibration', gatk_jar, params, job_name_prefix=vcf_job_name, hold_jid=snvs_holds)
            ApplyRecalibration_jobid = job.jobId
            job_ids += [job.jobId]
            ug_holds += [job.jobId]
            snvs_holds += [job.jobId]
            
            #validate format of vcf file
            job_name = 'ValidateApplyRecalibrationVcfFile_'+vcf_job_name
            cmd = '{} {} --input {} --output {}'.format(python, validate_vcf.__file__, recalibrated_snvs_vcf, recalibrated_snvs_vcf_validate)
            job = my.run_job(cmd, job_name, job_dir, hold_jid=ug_holds)
            ValidateApplyRecalibration_jobid = job.jobId
            job_ids += [job.jobId]
            ug_holds += [job.jobId]
            snvs_holds += [job.jobId]
            
            #verify validation
            job_name = 'VerifyValidateApplyRecalibrationVcfFile_'+vcf_job_name
            cmd = '{} {} --input {}'.format(python, verify_OK_file.__file__, recalibrated_snvs_vcf_validate)
            job = my.run_job(cmd, job_name, job_dir, hold_jid=ug_holds)
            VerifyValidateApplyRecalibration_jobid = job.jobId
            job_ids += [job.jobId]
            ug_holds += [job.jobId]
            snvs_holds += [job.jobId]
            verifyvalidatevcf_jobids += [job.jobId]
            
        #VariantFiltration indels
        VerifyValidateVariantFiltrationIndels_jobid = None
        if args.overwrite_all_output_files or args.overwrite_filtered_indels or not my.check_files(filtered_indels_vcf_validate, [(filtered_indels_vcf_validate,'OK')], [(filtered_indels_vcf, filtered_indels_vcf_validate)]):
            params = '-R "{}" --variant "{}" -o "{}" --filterExpression "QD < 2.0 || ReadPosRankSum < -20.0 || InbreedingCoeff < -0.8 || FS > 200.0"  --filterName GATKStandard'.format(
                ref, ug_raw_indels_vcf, filtered_indels_vcf)
            job = gatk(dirs, 'VariantFiltration', gatk_jar, params, job_name_prefix=vcf_job_name, hold_jid=indels_holds)
            VariantFiltrationIndels_jobid = job.jobId
            job_ids += [job.jobId]
            ug_holds += [job.jobId]
            indels_holds += [job.jobId]
            
            #validate format of vcf file
            job_name = 'ValidateVariantFiltrationVcfFile_'+vcf_job_name
            cmd = '{} {} --input {} --output {}'.format(python, validate_vcf.__file__, filtered_indels_vcf, filtered_indels_vcf_validate)
            job = my.run_job(cmd, job_name, job_dir, hold_jid=indels_holds)
            ValidateVariantFiltrationIndels_jobid = job.jobId
            job_ids += [job.jobId]
            ug_holds += [job.jobId]
            indels_holds += [job.jobId]
            
            #verify validation
            job_name = 'VerifyValidateVariantFiltrationVcfFile_'+vcf_job_name
            cmd = '{} {} --input {}'.format(python, verify_OK_file.__file__, filtered_indels_vcf_validate)
            job = my.run_job(cmd, job_name, job_dir, hold_jid=indels_holds)
            VerifyValidateVariantFiltrationIndels_jobid = job.jobId
            job_ids += [job.jobId]
            ug_holds += [job.jobId]
            verifyvalidatevcf_jobids += [job.jobId]
            indels_holds += [job.jobId]
        
        #CombineVariants
        CombineVariants_jobid = None
        if args.overwrite_all_output_files or args.overwrite_combine_variants or not my.check_files(analysis_ready_vcf_validate, [(analysis_ready_vcf_validate,'OK')], [(analysis_ready_vcf, analysis_ready_vcf_validate)]):
            params = '-R "{}" --variant:snv "{}" --variant:indel "{}" -o "{}" -assumeIdenticalSamples'.format(
                ref, recalibrated_snvs_vcf, filtered_indels_vcf, analysis_ready_vcf)
            #job = gatk(dirs, 'CombineVariants', gatk_jar, params, job_name_prefix=bam_job_name, hold_jid=ApplyRecalibration_jobid)
            job = gatk(dirs, 'CombineVariants', gatk_jar, params, job_name_prefix=vcf_job_name, hold_jid=ug_holds)
            CombineVariants_jobid = job.jobId
            job_ids += [job.jobId]
            ug_holds += [job.jobId]
            
            #validate format of vcf file
            job_name = 'ValidateCombineVariantsVcfFile_'+vcf_job_name
            cmd = '{} {} --input {} --output {}'.format(python, validate_vcf.__file__, analysis_ready_vcf, analysis_ready_vcf_validate)
            job = my.run_job(cmd, job_name, job_dir, hold_jid=ug_holds)
            ValidateCombineVariants_jobid = job.jobId
            job_ids += [job.jobId]
            ug_holds += [job.jobId]
        
            #verify validation
            job_name = 'VerifyValidateCombineVariantsVcfFile_'+vcf_job_name
            cmd = '{} {} --input {}'.format(python, verify_OK_file.__file__, analysis_ready_vcf_validate)
            job = my.run_job(cmd, job_name, job_dir, hold_jid=ug_holds)
            VerifyValidateCombineVariants_jobid = job.jobId
            job_ids += [job.jobId]
            ug_holds += [job.jobId]
            verifyvalidatevcf_jobids += [job.jobId]
        
        #VariantEval
        params = '-R "{}" --eval {} --dbsnp "{}" -o "{}" -ST Sample'.format(
            ref, analysis_ready_vcf, dbsnp, analysis_ready_vcf+'.VariantEval')
        job = gatk(dirs, 'VariantEval', gatk_jar, params, job_name_prefix=vcf_job_name, hold_jid=ug_holds)
        VariantEvalAnalysisReady_jobid = job.jobId
        job_ids += [job.jobId]
        
        #VariantEval only on select samples
        if samples:
            params = '-R "{}" --eval {} --dbsnp "{}" -o "{}" -ST Sample --sample "{}"'.format(
                ref, analysis_ready_vcf, dbsnp, analysis_ready_vcf+'.samples.VariantEval', ' --sample '.join(samples.split(' ')))
            job = gatk(dirs, 'VariantEval', gatk_jar, params, job_name_prefix=vcf_job_name, hold_jid=ug_holds)
            VariantEvalAnalysisReady_jobid = job.jobId
            job_ids += [job.jobId]            

    snpEff_jobid = None
    if SnpEffAnnotation:
        #snpEff
        job_name = 'snpEff_'+os.path.basename(snpEff_vcf)
        cmd = 'java -Xmx4G -jar {} eff -c {} -s {} -v -i vcf -o vcf {} {} > {};'.format(
            snpEff_jar, snpEff_config, snpEff_summary, snpEff_database, analysis_ready_vcf, snpEff_vcf)
        job = my.run_job(cmd, job_name, job_dir, hold_jid=None if NoVariantRecalibrator else CombineVariants_jobid)
        print cmd
        snpEff_jobid = job.jobId
        job_ids += [job.jobId]
    
    VariantAnnotator_jobid = None
    if not NoVariantAnnotator:
        #VariantAnnotator
        #TODO: doesn't work with snpEff 2.0.3
        #params = '-R "{}" -A SnpEff --variant {} --snpEffFile {} -o {}'.format(
        #    ref, analysis_ready_vcf, snpEff_vcf, annotated_vcf)
        params = '-R "{}" --variant {} -o {}'.format(
            ref, analysis_ready_vcf, annotated_vcf)
        job = gatk(dirs, 'VariantAnnotator', gatk_jar, params, job_name_prefix=vcf_job_name, hold_jid=ug_holds)
        VariantAnnotator_jobid = job.jobId
        job_ids += [job.jobId]
    
    return job_ids

def main():

    #command line arguments
    parser = argparse.ArgumentParser(parents=[my.default_parser()],
        description = 'Starting with one or a cohort of markdup bam files, realign indels, recalibrate base quality scores, genotype, and filter/recalibrate variants',
        epilog = 'pypeline.bam2vcf version 1.0β1 ©2011 Michael Yourshaw all rights reserved')
    parser.add_argument('--bams', nargs='+',
        help='input bam file(s)')
    parser.add_argument('--recalibrated_bams', nargs='*',
        help='input bam file(s) that have already been realigned and recalibrated')
    parser.add_argument('--vcf_prefix', required=True,
        help='/path/to/vcfs/file_prefix')
    parser.add_argument('--known_indels_vcfs', nargs='*',
        help='Input VCF file(s) with known indels')
    parser.add_argument('--targets',
        help='targeted intervals file')
    parser.add_argument('--gatk_jar',
        help='path to GenomeAnalysisTK jar')
    parser.add_argument('--analyzecovariates_jar',
        help='path to AnalyzeCovariates jar')
    parser.add_argument('--ref',
        help='path to reference genome fasta file')
    parser.add_argument('--dbsnp',
        help='path to dbsnp vcf')
    parser.add_argument('--hapmap',
        help='path to hapmap vcf')
    parser.add_argument('--omni',
        help='path to 1000 genomes omni vcf')
    parser.add_argument('--rscript',
        help='path to Rscript executable')
    parser.add_argument('--gatkrscripts',
        help='path to gatk R scripts resources directory')
    parser.add_argument('--output_mode', choices=('EMIT_VARIANTS_ONLY','EMIT_ALL_CONFIDENT_SITES','EMIT_ALL_SITES'), default='EMIT_VARIANTS_ONLY',
        help='Should unified genotyper output confident genotypes (i.e. including ref calls) or just the variants')
    parser.add_argument('--overwrite_all_output_files', action='store_true', default=False,
        help='overwriting any existing valid files')
    parser.add_argument('--overwrite_RealignerTargetCreator_intervals', action='store_true', default=False,
        help='redo RealignerTargetCreator, overwriting existing valid interval file')
    parser.add_argument('--overwrite_realigned_bams', action='store_true', default=False,
        help='redo IndelRealigner, overwriting existing valid realigned bam files')
    parser.add_argument('--overwrite_recal_csv', action='store_true', default=False,
        help='redo CountCovariates on realigned bams, overwriting existing valid .recal.csv files')
    parser.add_argument('--overwrite_recalibrated_bams', action='store_true', default=False,
        help='redo TableRecalibration, overwriting existing valid recalibrated bam files')
    parser.add_argument('--overwrite_recal2_csv', action='store_true', default=False,
        help='redo CountCovariates on recalibrated bams, overwriting existing valid .recal2.csv files')
    parser.add_argument('--overwrite_analyze_covariates_1', action='store_true', default=False,
        help='redo AnalyzeCovariates1, overwriting existing files')
    parser.add_argument('--overwrite_analyze_covariates_2', action='store_true', default=False,
        help='redo AnalyzeCovariates2, overwriting existing files')
    parser.add_argument('--overwrite_unified_genotyper', action='store_true', default=False,
        help='redo UnifiedGenotyper, overwriting existing valid vcf file')
    parser.add_argument('--overwrite_select_variants_snvs', action='store_true', default=False,
        help='redo SelectVariants for SNPs, overwriting existing valid vcf file')
    parser.add_argument('--overwrite_select_variants_indels', action='store_true', default=False,
        help='redo SelectVariants for indels, overwriting existing valid vcf file')
    parser.add_argument('--overwrite_recalibrated_snvs', action='store_true', default=False,
        help='redo ApplyRecalibration for vcf with SNPs, overwriting existing valid vcf file')
    parser.add_argument('--overwrite_filtered_indels', action='store_true', default=False,
        help='redo VariantFiltration for vcf with indels, overwriting existing valid vcf file')
    parser.add_argument('--overwrite_combine_variants', action='store_true', default=False,
        help='redo CombineVariants, overwriting existing valid vcf file')
    parser.add_argument('--NoRealignRecalibrate', action='store_true', default=False,
        help='omit RealignerTargetCreator, IndelRealigner, TableRecalibration, and realted steps')
    parser.add_argument('--NoUnifiedGenotyper', action='store_true', default=False,
        help='omit UnifiedGenotyper')
    parser.add_argument('--NoVariantRecalibrator', action='store_true', default=False,
        help='omit VariantRecalibrator and realted steps')
    parser.add_argument('--NoVariantAnnotator', action='store_true', default=False,
        help='omit VariantAnnotator step')
    parser.add_argument('--SnpEffAnnotation', action='store_true', default=False,
        help='run snpEff step')
    parser.add_argument('--snpEff_jar',
        help='path to snpEff.jar')
    parser.add_argument('--snpEff_config',
        help='path to snpEff.config (default: dirname(snpEff_jar)/snpEff.config)')
    parser.add_argument('--snpEff_database', default='GRCh37.64',
        help='snpEff database name (default: GRCh37.64)')
    parser.add_argument('--additional_VQSR_vcfs',
        help='additional vcf files to use with VQSR to get adequate statistical power (recommend >=30 samples for exomes)')
    parser.add_argument('--samples',
        help='names of samples in study to display QC metrics plots for') #Set[String]
    args = parser.parse_args()

    dirs = my.create_default_directory_structure(args.top_dir)
    config = my.get_config(args)

    job_ids = run(top_dir=args.top_dir, bams=args.bams, recalibrated_bams=args.recalibrated_bams, vcf_prefix=args.vcf_prefix,
        parent_job_dir=args.parent_job_dir, targets=args.targets, gatk_jar=args.gatk_jar,
        analyzecovariates_jar=args.analyzecovariates_jar, ref=args.ref, known_indels_vcfs=args.known_indels_vcfs,
        dbsnp=args.dbsnp, hapmap=args.hapmap, omni=args.omni, rscript=args.rscript,
        gatkrscripts=args.gatkrscripts, output_mode=args.output_mode,
        email=args.email, synchronous=args.synchronous, hold_jids=args.hold_jid if args.hold_jid else [],
        NoRealignRecalibrate=args.NoRealignRecalibrate, NoUnifiedGenotyper=args.NoUnifiedGenotyper,
        NoVariantRecalibrator=args.NoVariantRecalibrator, SnpEffAnnotation=args.SnpEffAnnotation,
        NoVariantAnnotator=args.NoVariantAnnotator, args=args,
        snpEff_jar=args.snpEff_jar, snpEff_database=args.snpEff_database, snpEff_config=args.snpEff_config,
        additional_VQSR_vcfs=args.additional_VQSR_vcfs, samples=args.samples)
    
    return job_ids

if __name__ == "__main__": sys.exit(main())
