#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-d /Users/myourshaw/lab/GMD/2011-04_sequencing/sequencing/metrics/pre-rmdup_metrics/summary /Users/myourshaw/lab/GMD/2011-04_sequencing/sequencing/metrics/pre-rmdup_metrics/* 
#-d /Users/myourshaw/lab/GMD/2011-04_sequencing/sequencing/metrics/post-rmdup_metrics/summary /Users/myourshaw/lab/GMD/2011-04_sequencing/sequencing/metrics/post-rmdup_metrics/* 
#-d /Users/myourshaw/lab/GMD/2011-04_sequencing/metrics_test /Users/myourshaw/lab/GMD/2011-04_sequencing/sequencing/metrics/post-rmdup_metrics/GMD* /Users/myourshaw/lab/GMD/2011-04_sequencing/sequencing/metrics/rmdup_metrics/GMD*

#-d /scratch1/tmp/myourshaw/hlee/metrics/bam_metrics/sample_metrics/consolidated /scratch1/tmp/myourshaw/hlee/metrics/bam_metrics/sample_metrics/*

#-d /scratch1/tmp/myourshaw/hlee/metrics/consolidated -i /scratch1/tmp/myourshaw/hlee/metrics/sample_bam_metrics/* /scratch1/tmp/myourshaw/hlee/metrics/library_bam_metrics/*

#-d /scratch1/tmp/myourshaw/hlee/metrics/consolidated -i /scratch1/tmp/myourshaw/hlee/metrics/sample_bam_metrics/*.librarycomplexity

#-d /scratch1/tmp/myourshaw/gmdjen/metrics/readgroup_bam_metrics/readgroup_metrics/consolidated -i /scratch1/tmp/myourshaw/gmdjen/metrics/readgroup_bam_metrics/readgroup_metrics/HWI-ST1070.94.*.*.novoalign.fixmate.bam.*
#-d /scratch1/tmp/myourshaw/gmdjen/metrics/sample_bam_metrics/sample_metrics/consolidated -i /scratch1/tmp/myourshaw/gmdjen/metrics/sample_bam_metrics/sample_metrics/*.sample.bam.*

#--input /scratch1/tmp/myourshaw/gmdjen/metrics/readgroup_bam_metrics/HWI-ST0860.135.5.CAGATC.novoalign.fixmate.bam.metrics.indexstats /scratch1/tmp/myourshaw/gmdjen/metrics/readgroup_bam_metrics/HWI-ST0860.135.5.CAGATC.novoalign.fixmate.bam.metrics.target_TruSeq_exome_targeted_regions.txt.interval_list.hsmetrics /scratch1/tmp/myourshaw/gmdjen/metrics/readgroup_bam_metrics/HWI-ST0860.135.5.CAGATC.novoalign.fixmate.bam.metrics.target_intervals_Ensembl_genes_protein_coding_cds.interval_list.hsmetrics /scratch1/tmp/myourshaw/gmdjen/metrics/readgroup_bam_metrics/HWI-ST0860.135.5.CAGATC.novoalign.fixmate.bam.metrics.target_refseq.coding.exon.merged.2bp.extended.hg19.Dec2010.interval_list.hsmetrics /scratch1/tmp/myourshaw/gmdjen/metrics/readgroup_bam_metrics/HWI-ST0860.135.5.CAGATC.novoalign.fixmate.bam.metrics.*_metrics /scratch1/tmp/myourshaw/gmdjen/metrics/readgroup_bam_metrics/HWI-ST0860.135.5.CAGATC.novoalign.fixmate.bam.metrics.gcbias.table /scratch1/tmp/myourshaw/gmdjen/metrics/readgroup_bam_metrics/HWI-ST0860.135.5.CAGATC.novoalign.fixmate.bam.metrics.gcbias.summary /scratch1/tmp/myourshaw/gmdjen/metrics/readgroup_bam_metrics/HWI-ST0860.135.5.CAGATC.novoalign.fixmate.bam.metrics.librarycomplexity --output_dir /scratch1/tmp/myourshaw/gmdjen/metrics/readgroup_bam_metrics/consolidated

#--input /scratch1/tmp/myourshaw/test/metrics/readgroup_bam_metrics/HWI-ST973.83.5.CGATGT.novoalign.fixmate.bam.metrics.indexstats /scratch1/tmp/myourshaw/test/metrics/readgroup_bam_metrics/HWI-ST973.83.5.CGATGT.novoalign.fixmate.bam.metrics.target_TruSeq_exome_targeted_regions.txt.interval_list.hsmetrics /scratch1/tmp/myourshaw/test/metrics/readgroup_bam_metrics/HWI-ST973.83.5.CGATGT.novoalign.fixmate.bam.metrics.target_intervals_Ensembl_genes_protein_coding_cds.interval_list.hsmetrics /scratch1/tmp/myourshaw/test/metrics/readgroup_bam_metrics/HWI-ST973.83.5.CGATGT.novoalign.fixmate.bam.metrics.target_refseq.coding.exon.merged.2bp.extended.hg19.Dec2010.interval_list.hsmetrics /scratch1/tmp/myourshaw/test/metrics/readgroup_bam_metrics/HWI-ST973.83.5.CGATGT.novoalign.fixmate.bam.metrics.*_metrics /scratch1/tmp/myourshaw/test/metrics/readgroup_bam_metrics/HWI-ST973.83.5.CGATGT.novoalign.fixmate.bam.metrics.gcbias.table /scratch1/tmp/myourshaw/test/metrics/readgroup_bam_metrics/HWI-ST973.83.5.CGATGT.novoalign.fixmate.bam.metrics.gcbias.summary /scratch1/tmp/myourshaw/test/metrics/readgroup_bam_metrics/HWI-ST973.83.5.CGATGT.novoalign.fixmate.bam.metrics.librarycomplexity --output_dir /scratch1/tmp/myourshaw/test/metrics/readgroup_bam_metrics/consolidated

import sys,os
import re
import traceback
import argparse
import my

def out(file,string):
	if (not string.endswith('\n')): string += '\n'
	file.write(string)
	
def parse_metrics(output_dir, input_files, metrics_class, header_startswith, output_file_prefix=None):
    
    #http://picard.sourceforge.net/picard-metric-definitions.shtml#DuplicationMetrics
    docs={'AlignmentSummaryMetrics': """##AlignmentSummaryMetrics
##High level metrics about the alignment of reads within a SAM file, produced by the CollectAlignmentSummaryMetrics program and usually stored in a file with the extension ".alignment_summary_metrics".
##CATEGORY: One of either UNPAIRED (for a fragment run), FIRST_OF_PAIR when metrics are for only the first read in a paired run, SECOND_OF_PAIR when the metrics are for only the second read in a paired run or PAIR when the metrics are aggregated for both first and second reads in a pair.
##TOTAL_READS: The total number of reads including all PF and non-PF reads. When CATEGORY equals PAIR this value will be 2x the number of clusters.
##PF_READS: The number of PF reads where PF is defined as passing Illumina's filter.
##PCT_PF_READS: The percentage of reads that are PF (PF_READS / TOTAL_READS)
##PF_NOISE_READS: The number of PF reads that are marked as noise reads. A noise read is one which is composed entirey of A bases and/or N bases. These reads are marked as they are usually artifactual and are of no use in downstream analysis.
##PF_READS_ALIGNED: The number of PF reads that were aligned to the reference sequence. This includes reads that aligned with low quality (i.e. their alignments are ambiguous).
##PCT_PF_READS_ALIGNED: The percentage of PF reads that aligned to the reference sequence. PF_READS_ALIGNED / PF_READS
##PF_ALIGNED_BASES: The total number of aligned bases, in all mapped PF reads, that are aligned to the reference sequence.
##PF_HQ_ALIGNED_READS: The number of PF reads that were aligned to the reference sequence with a mapping quality of Q20 or higher signifying that the aligner estimates a 1/100 (or smaller) chance that the alignment is wrong.
##PF_HQ_ALIGNED_BASES: The number of bases aligned to the reference sequence in reads that were mapped at high quality. Will usually approximate PF_HQ_ALIGNED_READS * READ_LENGTH but may differ when either mixed read lengths are present or many reads are aligned with gaps.
##PF_HQ_ALIGNED_Q20_BASES: The subest of PF_HQ_ALIGNED_BASES where the base call quality was Q20 or higher.
##PF_HQ_MEDIAN_MISMATCHES: The median number of mismatches versus the reference sequence in reads that were aligned to the reference at high quality (i.e. PF_HQ_ALIGNED READS).
##PF_MISMATCH_RATE: The rate of bases mismatching the reference for all bases aligned to the reference sequence.
##PF_HQ_ERROR_RATE: The percentage of bases that mismatch the reference in PF HQ aligned reads.
##PF_INDEL_RATE: The number of insertion and deletion events per 100 aligned bases. Uses the number of events as the numerator, not the number of inserted or deleted bases.
##MEAN_READ_LENGTH: The mean read length of the set of reads examined. When looking at the data for a single lane with equal length reads this number is just the read length. When looking at data for merged lanes with differing read lengths this is the mean read length of all reads.
##READS_ALIGNED_IN_PAIRS: The number of aligned reads whose mate pair was also aligned to the reference.
##PCT_READS_ALIGNED_IN_PAIRS: The percentage of reads whose mate pair was also aligned to the reference. READS_ALIGNED_IN_PAIRS / PF_READS_ALIGNED
##BAD_CYCLES: The number of instrument cycles in which 80% or more of base calls were no-calls.
##STRAND_BALANCE: The number of PF reads aligned to the positive strand of the genome divided by the number of PF reads aligned to the genome.
##PCT_CHIMERAS: The percentage of reads that map outside of a maximum insert size (usually 100kb) or that have the two ends mapping to different chromosomes.
##PCT_ADAPTER: The percentage of PF reads that are unaligned and match to a known adapter sequence right from the start of the read.""",

'DuplicationMetrics': """##DuplicationMetrics
##Metrics that are calculated during the process of marking duplicates within a stream of SAMRecords.
##LIBRARY: The library on which the duplicate marking was performed.
##UNPAIRED_READS_EXAMINED: The number of mapped reads examined which did not have a mapped mate pair, either because the read is unpaired, or the read is paired to an unmapped mate.
##READ_PAIRS_EXAMINED: The number of mapped read pairs examined.
##UNMAPPED_READS: The total number of unmapped reads examined.
##UNPAIRED_READ_DUPLICATES: The number of fragments that were marked as duplicates.
##READ_PAIR_DUPLICATES: The number of read pairs that were marked as duplicates.
##READ_PAIR_OPTICAL_DUPLICATES: The number of read pairs duplicates that were caused by optical duplication. Value is always < READ_PAIR_DUPLICATES, which counts all duplicates regardless of source.
##PERCENT_DUPLICATION: The percentage of mapped sequence that is marked as duplicate.
##ESTIMATED_LIBRARY_SIZE: The estimated number of unique molecules in the library based on PE duplication.""",

'ExtractIlluminaBarcodes.BarcodeMetric': """##ExtractIlluminaBarcodes.BarcodeMetric
##Metrics produced by the ExtractIlluminaBarcodes program that is used to parse data in the basecalls directory and determine to which barcode each read should be assigned.
##BARCODE: The barcode (from the set of expected barcodes) for which the following metrics apply. Note that the "symbolic" barcode of NNNNNN is used to report metrics for all reads that do not match a barcode.
##BARCODE_NAME:
##LIBRARY_NAME:
##READS: The total number of reads matching the barcode.
##PF_READS: The number of PF reads matching this barcode (always less than or equal to READS).
##PERFECT_MATCHES: The number of all reads matching this barcode that matched with 0 errors or no-calls.
##PF_PERFECT_MATCHES: The number of PF reads matching this barcode that matched with 0 errors or no-calls.
##ONE_MISMATCH_MATCHES: The number of all reads matching this barcode that matched with 1 error or no-call.
##PF_ONE_MISMATCH_MATCHES: The number of PF reads matching this barcode that matched with 1 error or no-call.
##PCT_MATCHES: The percentage of all reads in the lane that matched to this barcode.
##RATIO_THIS_BARCODE_TO_BEST_BARCODE_PCT: The rate of all reads matching this barcode to all reads matching the most prevelant barcode. For the most prevelant barcode this will be 1, for all others it will be less than 1. One over the lowest number in this column gives you the fold-difference in representation between barcodes.
##PF_PCT_MATCHES: The percentage of PF reads in the lane that matched to this barcode.
##PF_RATIO_THIS_BARCODE_TO_BEST_BARCODE_PCT: The rate of PF reads matching this barcode to PF reads matching the most prevelant barcode. For the most prevelant barcode this will be 1, for all others it will be less than 1. One over the lowest number in this column gives you the fold-difference in representation of PF reads between barcodes.
##PF_NORMALIZED_MATCHES: The "normalized" matches to each barcode. This is calculated as the number of pf reads matching this barcode over the sum of all pf reads matching any barcode (excluding orphans). If all barcodes are represented equally this will be 1.""",

'GcBiasDetailMetrics': """##GcBiasDetailMetrics
##Class that holds detailed metrics about reads that fall within windows of a certain GC bin on the reference genome.
##GC: The G+C content of the reference sequence represented by this bin. Values are from 0% to 100%
##WINDOWS: The number of windows on the reference genome that have this G+C content.
##READ_STARTS: The number of reads whose start position is at the start of a window of this GC.
##MEAN_BASE_QUALITY: The mean quality (determined via the error rate) of all bases of all reads that are assigned to windows of this GC.
##NORMALIZED_COVERAGE: The ration of "coverage" in this GC bin vs. the mean coverage of all GC bins. A number of 1 represents mean coverage, a number less than one represents lower than mean coverage (e.g. 0.5 means half as much coverage as average) while a number greater than one represents higher than mean coverage (e.g. 3.1 means this GC bin has 3.1 times more reads per window than average).
##ERROR_BAR_WIDTH: The radius of error bars in this bin based on the number of observations made. For example if the normalized coverage is 0.75 and the error bar width is 0.1 then the error bars would be drawn from 0.65 to 0.85.""",

'GcBiasSummaryMetrics': """##GcBiasSummaryMetrics
##High level metrics that capture how biased the coverage in a certain lane is.
##WINDOW_SIZE: The window size on the genome used to calculate the GC of the sequence.
##TOTAL_CLUSTERS: The total number of clusters that were seen in the gc bias calculation.
##ALIGNED_READS: The total number of aligned reads used to compute the gc bias metrics.
##AT_DROPOUT: Illumina-style AT dropout metric. Calculated by taking each GC bin independently and calculating (%ref_at_gc - %reads_at_gc) and summing all positive values for GC=[0..50].
##GC_DROPOUT: Illumina-style GC dropout metric. Calculated by taking each GC bin independently and calculating (%ref_at_gc - %reads_at_gc) and summing all positive values for GC=[50..100].""",

'HsMetrics': """##HsMetrics
##The set of metrics captured that are specific to a hybrid selection analysis.
##BAIT_SET: The name of the bait set used in the hybrid selection.
##GENOME_SIZE: The number of bases in the reference genome used for alignment.
##BAIT_TERRITORY: The number of bases which have one or more baits on top of them.
##TARGET_TERRITORY: The unique number of target bases in the experiment where target is usually exons etc.
##BAIT_DESIGN_EFFICIENCY: Target terrirtoy / bait territory. 1 == perfectly efficient, 0.5 = half of baited bases are not target.
##TOTAL_READS: The total number of reads in the SAM or BAM file examine.
##PF_READS: The number of reads that pass the vendor's filter.
##PF_UNIQUE_READS: The number of PF reads that are not marked as duplicates.
##PCT_PF_READS: PF reads / total reads. The percent of reads passing filter.
##PCT_PF_UQ_READS: PF Unique Reads / Total Reads.
##PF_UQ_READS_ALIGNED: The number of PF unique reads that are aligned with mapping score > 0 to the reference genome.
##PCT_PF_UQ_READS_ALIGNED: PF Reads Aligned / PF Reads.
##PF_UQ_BASES_ALIGNED: The number of bases in the PF aligned reads that are mapped to a reference base. Accounts for clipping and gaps.
##ON_BAIT_BASES: The number of PF aligned bases that mapped to a baited region of the genome.
##NEAR_BAIT_BASES: The number of PF aligned bases that mapped to within a fixed interval of a baited region, but not on a baited region.
##OFF_BAIT_BASES: The number of PF aligned bases that mapped to neither on or near a bait.
##ON_TARGET_BASES: The number of PF aligned bases that mapped to a targetted region of the genome.
##PCT_SELECTED_BASES: On+Near Bait Bases / PF Bases Aligned.
##PCT_OFF_BAIT: The percentage of aligned PF bases that mapped neither on or near a bait.
##ON_BAIT_VS_SELECTED: The percentage of on+near bait bases that are on as opposed to near.
##MEAN_BAIT_COVERAGE: The mean coverage of all baits in the experiment.
##MEAN_TARGET_COVERAGE: The mean coverage of targets that recieved at least coverage depth = 2 at one base.
##PCT_USABLE_BASES_ON_BAIT: The number of aligned, de-duped, on-bait bases out of the PF bases available.
##PCT_USABLE_BASES_ON_TARGET: The number of aligned, de-duped, on-target bases out of the PF bases available.
##FOLD_ENRICHMENT: The fold by which the baited region has been amplified above genomic background.
##ZERO_CVG_TARGETS_PCT: The number of targets that did not reach coverage=2 over any base.
##FOLD_80_BASE_PENALTY: The fold over-coverage necessary to raise 80% of bases in "non-zero-cvg" targets to the mean coverage level in those targets.
##PCT_TARGET_BASES_2X: The percentage of ALL target bases acheiving 2X or greater coverage.
##PCT_TARGET_BASES_10X: The percentage of ALL target bases acheiving 10X or greater coverage.
##PCT_TARGET_BASES_20X: The percentage of ALL target bases acheiving 20X or greater coverage.
##PCT_TARGET_BASES_30X: The percentage of ALL target bases acheiving 30X or greater coverage.
##HS_LIBRARY_SIZE: The estimated number of unique molecules in the selected part of the library.
##HS_PENALTY_10X: The "hybrid selection penalty" incurred to get 80% of target bases to 10X. This metric should be interpreted as: if I have a design with 10 megabases of target, and want to get 10X coverage I need to sequence until PF_ALIGNED_BASES = 10^7 * 10 * HS_PENALTY_10X.
##HS_PENALTY_20X: The "hybrid selection penalty" incurred to get 80% of target bases to 20X. This metric should be interpreted as: if I have a design with 10 megabases of target, and want to get 20X coverage I need to sequence until PF_ALIGNED_BASES = 10^7 * 20 * HS_PENALTY_20X.
##HS_PENALTY_30X: The "hybrid selection penalty" incurred to get 80% of target bases to 10X. This metric should be interpreted as: if I have a design with 10 megabases of target, and want to get 30X coverage I need to sequence until PF_ALIGNED_BASES = 10^7 * 30 * HS_PENALTY_30X.
##AT_DROPOUT: A measure of how undercovered <= 50% GC regions are relative to the mean. For each GC bin [0..50] we calculate a = % of target territory, and b = % of aligned reads aligned to these targets. AT DROPOUT is then abs(sum(a-b when a-b < 0)). E.g. if the value is 5% this implies that 5% of total reads that should have mapped to GC<=50% regions mapped elsewhere.
##GC_DROPOUT: A measure of how undercovered >= 50% GC regions are relative to the mean. For each GC bin [50..100] we calculate a = % of target territory, and b = % of aligned reads aligned to these targets. GC DROPOUT is then abs(sum(a-b when a-b < 0)). E.g. if the value is 5% this implies that 5% of total reads that should have mapped to GC>=50% regions mapped elsewhere.""",

'InsertSizeMetrics': """##InsertSizeMetrics
##Metrics about the insert size distribution of a paired-end library, created by the CollectInsertSizeMetrics program and usually written to a file with the extension ".insert_size_metrics". In addition the insert size distribution is plotted to a file with the extension ".insert_size_histogram.pdf".
##MEDIAN_INSERT_SIZE: The MEDIAN insert size of all paired end reads where both ends mapped to the same chromosome.
##MEDIAN_ABSOLUTE_DEVIATION: The median absolute deviation of the distribution. If the distribution is essentially normal then the standard deviation can be estimated as ~1.4826 * MAD.
##MIN_INSERT_SIZE: The minimum measured insert size. This is usually 1 and not very useful as it is likely artifactual.
##MAX_INSERT_SIZE: The maximum measure insert size by alignment. This is usually very high representing either an artifact or possibly the presence of a structural re-arrangement.
##MEAN_INSERT_SIZE: The mean insert size of the "core" of the distribution. Artefactual outliers in the distribution often cause calculation of nonsensical mean and stdev values. To avoid this the distribution is first trimmed to a "core" distribution of +/- N median absolute deviations around the median insert size. By default N=10, but this is configurable.
##STANDARD_DEVIATION: Standard deviation of insert sizes over the "core" of the distrubution.
##READ_PAIRS: The total number of read pairs that were examined in the entire distribution.
##PAIR_ORIENTATION: The pair orientation of the reads in this data category.
##WIDTH_OF_10_PERCENT: The "width" of the bins, centered around the median, that encompass 10% of all read pairs.
##WIDTH_OF_20_PERCENT: The "width" of the bins, centered around the median, that encompass 20% of all read pairs.
##WIDTH_OF_30_PERCENT: The "width" of the bins, centered around the median, that encompass 30% of all read pairs.
##WIDTH_OF_40_PERCENT: The "width" of the bins, centered around the median, that encompass 40% of all read pairs.
##WIDTH_OF_50_PERCENT: The "width" of the bins, centered around the median, that encompass 50% of all read pairs.
##WIDTH_OF_60_PERCENT: The "width" of the bins, centered around the median, that encompass 60% of all read pairs.
##WIDTH_OF_70_PERCENT: The "width" of the bins, centered around the median, that encompass 70% of all read pairs. This metric divided by 2 should approximate the standard deviation when the insert size distribution is a normal distribution.
##WIDTH_OF_80_PERCENT: The "width" of the bins, centered around the median, that encompass 80% of all read pairs.
##WIDTH_OF_90_PERCENT: The "width" of the bins, centered around the median, that encompass 90% of all read pairs.
##WIDTH_OF_99_PERCENT: The "width" of the bins, centered around the median, that encompass 100% of all read pairs.""",

'MultilevelMetrics': """##MultilevelMetrics
##SAMPLE: The sample to which these metrics apply. If null, it means they apply to all reads in the file.
##LIBRARY: The library to which these metrics apply. If null, it means that the metrics were accumulated at the sample level.
##READ_GROUP: The read group to which these metrics apply. If null, it means that the metrics were accumulated at the library or sample level.""",

'RnaSeqMetrics': """##RnaSeqMetrics
##Metrics about the alignment of RNA-seq reads within a SAM file to genes, produced by the CollectRnaSeqMetrics program and usually stored in a file with the extension ".rna_metrics".
##PF_BASES: The total number of PF bases including non-aligned reads.
##PF_ALIGNED_BASES: The total number of aligned PF bases. Non-primary alignments are not counted. Bases in aligned reads that do not correspond to reference (e.g. soft clips, insertions) are not counted.
##RIBOSOMAL_BASES: Number of bases in primary aligments that align to ribosomal sequence.
##CODING_BASES: Number of bases in primary aligments that align to a non-UTR coding base for some gene, and not ribosomal sequence.
##UTR_BASES: Number of bases in primary aligments that align to a UTR base for some gene, and not a coding base.
##INTRONIC_BASES: Number of bases in primary aligments that align to an intronic base for some gene, and not a coding or UTR base.
##INTERGENIC_BASES: Number of bases in primary aligments that do not align to any gene.
##IGNORED_READS: Number of primary alignments that map to a sequence specified on command-line as IGNORED_SEQUENCE. These are not counted in PF_ALIGNED_BASES, CORRECT_STRAND_READS, INCORRECT_STRAND_READS, or any of the base-counting metrics. These reads are counted in PF_BASES.
##CORRECT_STRAND_READS: Number of aligned reads that map to the correct strand. 0 if library is not strand-specific.
##INCORRECT_STRAND_READS: Number of aligned reads that map to the incorrect strand. 0 if library is not strand-specific.
##PCT_RIBOSOMAL_BASES: RIBOSOMAL_BASES / PF_ALIGNED_BASES
##PCT_CODING_BASES: CODING_BASES / PF_ALIGNED_BASES
##PCT_UTR_BASES: UTR_BASES / PF_ALIGNED_BASES
##PCT_INTRONIC_BASES: INTRONIC_BASES / PF_ALIGNED_BASES
##PCT_INTERGENIC_BASES: INTERGENIC_BASES / PF_ALIGNED_BASES
##PCT_MRNA_BASES: PCT_UTR_BASES + PCT_CODING_BASES
##PCT_USABLE_BASES: The percentage of bases mapping to mRNA divided by the total number of PF bases.
##PCT_CORRECT_STRAND_READS: CORRECT_STRAND_READS/(CORRECT_STRAND_READS + INCORRECT_STRAND_READS). 0 if library is not strand-specific.
##MEDIAN_CV_COVERAGE: The median CV of coverage of the 1000 most highly expressed transcripts.
##MEDIAN_5PRIME_BIAS: The median 5 prime bias of the 1000 most highly expressed transcripts, where 5 prime bias is calculated per transcript as: mean coverage of the 5' most 100 bases divided by the mean coverage of the whole transcript.
##MEDIAN_3PRIME_BIAS: The median 3 prime bias of the 1000 most highly expressed transcripts, where 3 prime bias is calculated per transcript as: mean coverage of the 3' most 100 bases divided by the mean coverage of the whole transcript.
##MEDIAN_5PRIME_TO_3PRIME_BIAS: The ratio of coverage at the 5' end of to the 3' end based on the 1000 most highly expressed transcripts.""",
}

    if not output_file_prefix:
        output_file_prefix = metrics_class
    if (len(input_files) > 0):
        output_file_name = os.path.join(output_dir,output_file_prefix+'.txt')
        output_file = open(output_file_name,'w')
        if docs.get(metrics_class):
            out(output_file, docs[metrics_class])
        out(output_file, '##Consolidated Files:')
        for file in input_files:
            out (output_file,'## %s' % (file))
        header_written = False
        for file in input_files:
            file_id = os.path.join(os.path.basename(os.path.dirname(file)),os.path.basename(file))
            linecount = 0
            metrics_class_found = False
            for line in open(file):
                linecount += 1
                if (line.strip() == ''): continue
                if (line.startswith('## HISTOGRAM')):
                    break
                if (not metrics_class_found and line.startswith('## METRICS CLASS') and line.rstrip().endswith('.'+metrics_class)):
                    metrics_class_found = True
                    continue
                if(metrics_class_found):
                    if(not header_written):
                        if(line.startswith(header_startswith)):
                            out(output_file,'%s\t%s\t%s' % ('#METRICS_CLASS','FILE',line))
                            header_written = True
                            continue
                    else:
                        fields = line.rstrip('\n').split('\t')
                        if (len(fields) > 0 and fields[0] != header_startswith):
                            out(output_file,'%s\t%s\t%s' % (metrics_class,file_id,line))

def parse_histograms(output_dir, input_files, metrics_class, header_startswith):
    if (len(input_files) > 0):
        output_file_name = os.path.join(output_dir,metrics_class+'.txt')
        output_file = open(output_file_name,'w')
        parallel_output_file_name = os.path.join(output_dir,metrics_class+'.parallel_histograms.txt')
        parallel_output_file = open(parallel_output_file_name,'w')
        out(output_file, '##Consolidated Files:')
        for file in input_files:
            out (output_file,'## %s' % (file))
            out (parallel_output_file,'## %s' % (file))
        header_written = False
        out_list=[]
        header = ()
        min_x = None
        max_x = None
        x_set = set()
        for file in input_files:
            file_id = os.path.join(os.path.basename(os.path.dirname(file)),os.path.basename(file))
            linecount = 0
            histogram_found = False
            for line in open(file):
                linecount += 1
                if (line.strip() == ''): continue
                if (not histogram_found and line.startswith('## HISTOGRAM')):
                    histogram_found = True
                    continue
                if(histogram_found):
                    if(not header_written):
                        if(line.startswith(header_startswith)):
                            out(output_file,'%s\t%s\t%s' % ('#METRICS_CLASS','FILE',line))
                            header_fields = line.rstrip('\n').split('\t')
                            header = (header_fields[0],header_fields[1])
                            header_written = True
                            out(parallel_output_file,"## %s for each file" % header_fields[1])
                            out(parallel_output_file,"#%s\t%s" % (header_fields[0],'\t'.join(input_files)))
                            continue
                    else:
                        fields = line.rstrip('\n').split('\t')
                        if (len(fields) > 0 and fields[0] != header_startswith):
                            out(output_file,'%s\t%s\t%s' % (metrics_class,file_id,line))
                            x = float(fields[0])
                            y = float(fields[1])
                            out_list.append((file,x,y))
                            x_set.add(x)
                            if (min_x == None or x < min_x): min_x = x
                            if (max_x == None or x > max_x): max_x = x
        x_list = list(x_set)
        x_list.sort()
        for i in x_list:
            line_list = ['0']*len(input_files)
            these_x = [o for o in out_list if o[1] == i]
            for fxy in these_x:
                pos = input_files.index(fxy[0])
                line_list[pos] = str(fxy[2])
            out(parallel_output_file,"%s\t%s" % (i,'\t'.join(line_list)) )

def main():
    #command line arguments
    parser = argparse.ArgumentParser(
        description = 'Consolidates multiple picard metrics files into a nice format for a spreadsheet',
        epilog = 'pypeline.picard_metrics_consolidate version 1.0β1 ©2011-2012 Michael Yourshaw all rights reserved')
    parser.add_argument('--input', '-i', required=True, nargs="+",
        help='input files from output of picard MarkDuplicates (output must end with .markdup.metrics or .duplication.metrics), BamIndexStats, CalculateHsMetrics, CollectMultipleMetrics, CollectGcBiasMetrics, and/or EstimateLibraryComplexity')
    parser.add_argument('--output_dir', '-d', '-o', required=True,
        help='output directory')

    args = parser.parse_args()

    output_dir = args.output_dir
    my.makedir(output_dir)

    files_alignment_summary_metrics = []
    files_gcbias_summary = []
    files_gcbias_table = []
    files_hsmetrics = []
    files_idxstats = []
    files_indexstats = []
    files_insert_size_metrics = []
    files_insert_size_histograms = []
    files_librarycomplexity = []
    files_librarycomplexityhistograms = []
    files_quality_by_cycle_metrics = []
    files_quality_distribution_metrics = []
    files_duplication_metrics = []
    files_duplication_histograms = []
    files = [f for f in my.unglob(args.input) if
            f.endswith('.alignment_summary_metrics') or f.endswith('.gcbias.summary') or f.endswith('.gcbias.table') 
            or f.endswith('.hsmetrics') or f.endswith('.idxstats') or f.endswith('.indexstats')
            or f.endswith('.insert_size_metrics') or f.endswith('.librarycomplexity')
            or f.endswith('.quality_by_cycle_metrics') or f.endswith('.quality_distribution_metrics')
            or f.endswith('.markdup.metrics') or f.endswith('.duplication.metrics')]
    for file in files:
        if file.endswith('.pdf'): continue
        elif file.endswith('alignment_summary_metrics'): files_alignment_summary_metrics.append(file)
        elif file.endswith('.gcbias.summary'): files_gcbias_summary.append(file)
        elif file.endswith('.gcbias.table'): files_gcbias_table.append(file)
        elif file.endswith('.hsmetrics'): files_hsmetrics.append(file)
        elif file.endswith('.idxstats'): files_idxstats.append(file)
        elif file.endswith('.indexstats'): files_indexstats.append(file)
        elif file.endswith('.insert_size_metrics'):
            files_insert_size_metrics.append(file)
            files_insert_size_histograms.append(file)
        elif file.endswith('.librarycomplexity'):
            files_librarycomplexity.append(file)
            files_librarycomplexityhistograms.append(file)
        elif file.endswith('.quality_by_cycle_metrics'): files_quality_by_cycle_metrics.append(file)
        elif file.endswith('.quality_distribution_metrics'): files_quality_distribution_metrics.append(file)
        elif file.endswith('.duplication.metrics') or file.endswith('.markdup.metrics'):
            files_duplication_metrics.append(file)
            files_duplication_histograms.append(file)
            
    parse_metrics(output_dir,files_alignment_summary_metrics,'AlignmentSummaryMetrics','CATEGORY')
    parse_metrics(output_dir,files_gcbias_summary,'GcBiasSummaryMetrics','WINDOW_SIZE')
    parse_metrics(output_dir,files_gcbias_table,'GcBiasDetailMetrics','GC')
    parse_metrics(output_dir,files_hsmetrics,'HsMetrics','BAIT_SET')
    #parse_metrics_class_files(output_dir,files_idxstats,'','')
    #parse_metrics_class_files(output_dir,files_indexstats,'','')
    parse_metrics(output_dir,files_insert_size_metrics,'InsertSizeMetrics','MEDIAN_INSERT_SIZE')
    parse_histograms(output_dir,files_insert_size_histograms,'InsertSizeHistograms','insert_size')
    parse_metrics(output_dir,files_librarycomplexity,'DuplicationMetrics','LIBRARY', 'LibraryComplexity')
    parse_histograms(output_dir,files_librarycomplexityhistograms,'LibraryComplexityHistograms','duplication_group_count')
    parse_histograms(output_dir,files_quality_by_cycle_metrics,'QualityByCycleHistograms','CYCLE')
    parse_histograms(output_dir,files_quality_distribution_metrics,'QualityDistributionsHistograms','QUALITY')
    parse_metrics(output_dir,files_duplication_metrics,'DuplicationMetrics','LIBRARY')
    parse_histograms(output_dir,files_duplication_histograms,'DuplicationHistograms','BIN')


if __name__ == "__main__": sys.exit(main())
