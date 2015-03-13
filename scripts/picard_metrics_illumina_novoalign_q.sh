#!/bin/bash

#picard_metrics_illumina_q.sh

usage="usage: $0 [<email address> (default: myourshaw@ucla.edu)] [<ref.fasta> | <ref.fa> (default human_g1k_v37.fasta)] <directory for output> <list of sam or bam files>";
#get path and parent directory of this script
#http://hintsforums.macworld.com/archive/index.php/t-73839.html
# if the call came via symlink, then use its target instead:
arg=$0; [[ -L $0 ]] && arg=$(stat -L -c'%n' "$0")
script_path=$(2>/dev/null cd "${arg%/*}" >&2; echo "`pwd -P`/${arg##*/}")
script_dir=$(dirname "$script_path")

script_dir=$(dirname "$script_path")

HERE=${script_dir};

EMAIL="myourshaw@ucla.edu";
if [[ "$1" =~ '\S+@\S+\.\S+'  ]]; then EMAIL="$1"; shift; fi

source ${HERE}/b37.variables.sh $EMAIL;

ref="$B37_FASTA";
if [[ ${1##*.} == "fasta" || ${1##*.} == "fa" ]]; then ref="$1"; shift; fi
if (( $# < 1 )); then echo $usage; exit; fi

output_dir="$1"; shift;

bait="$SURESELECT_50MB_INTERVAL_LIST";
target_all_genes="$ENSEMBL_ALL_GENES";
target_protein_coding_genes="$ENSEMBL_PROTEIN_CODING_GENES"

for b in $@; do
 if [[ ! -e $b ]]; then echo "file does not exist [$b]"; exit; fi
done

wd=`pwd`;

for bam in $@; do
	#output files
	statsbase="${output_dir}/`basename ${bam}`.picard_metrics";
	validate="$statsbase.validate";
	idxstats="$statsbase.idxstats";
	hsmetrics_all="$statsbase.all_genes.hsmetrics";
	hsmetrics_protein="$statsbase.protein_coding_genes.hsmetrics";
	indexstats="$statsbase.indexstats";
	gcbiasmetrics="$statsbase.gcbias.table";
	gcbiasmetricschart="$statsbase.gcbias.pdf";
	gcbiasmetricssummary="$statsbase.gcbias.summary";
	librarycomplexity="$statsbase.librarycomplexity";
	#deprecated:
	alignmentsummarymetrics="$bam.alignmentsummarymetrics";
	insertsizemetrics="$bam.insertsize.table";
	insertsizehistogram="$bam.insertsize.histogram";
	meanqualitybycycle="$bam.meanqualitybycycle.table";
	meanqualitybycyclechart="$bam.meanqualitybycycle.pdf";
	qualityscoredistribution="$bam.qualityscoredistribution.table";
	qualityscoredistributionchart="$bam.qualityscoredistribution.pdf";
	
	qout="$output_dir/qout";
	mkdir -p "$qout";
	qsub_allq="$QSUB_ALLQ -e $qout -o $qout";
	qsub_lomem="$QSUB_LOMEM -e $qout -o $qout";
	qsub_himem="$QSUB_HIMEM -e $qout -o $qout";

	#idxstats
	name="idxstats_`basename $idxstats`";
	cmd=" \
	$SAMTOOLS idxstats $bam > $idxstats \
	;";
	idxstatsId=`echo "$cmd" | $qsub_allq -N $name`;
	echo "$idxstatsId $name $cmd";

	#CalculateHsMetrics.all_genes
	name="CalculateHsMetrics.all_genes_`basename $hsmetrics_all`";
	cmd=" \
	java -Xmx5g -Djava.io.tmpdir=$TMPSCRATCH \
	 -jar $PICARD/CalculateHsMetrics.jar \
	 TMP_DIR=$TMPSCRATCH \
	 VERBOSITY=INFO \
	 QUIET=false \
	 VALIDATION_STRINGENCY=SILENT \
	 COMPRESSION_LEVEL=5 \
	 MAX_RECORDS_IN_RAM=1000000 \
	 BAIT_INTERVALS=$bait \
	 TARGET_INTERVALS=$target_all_genes \
	 INPUT=$bam \
	 OUTPUT=$hsmetrics_all \
	;";
	CalculateHsMetricsAllId=`echo "$cmd" | $qsub_lomem -N $name`;
	echo "$CalculateHsMetricsAllId $name $cmd";

	#CalculateHsMetrics.protein_coding_genes
	name="CalculateHsMetrics.protein_coding_genes_`basename $hsmetrics_protein`";
	cmd=" \
	java -Xmx5g -Djava.io.tmpdir=$TMPSCRATCH \
	 -jar $PICARD/CalculateHsMetrics.jar \
	 TMP_DIR=$TMPSCRATCH \
	 VERBOSITY=INFO \
	 QUIET=false \
	 VALIDATION_STRINGENCY=SILENT \
	 COMPRESSION_LEVEL=5 \
	 MAX_RECORDS_IN_RAM=1000000 \
	 BAIT_INTERVALS=$bait \
	 TARGET_INTERVALS=$target_protein_coding_genes \
	 INPUT=$bam \
	 OUTPUT=$hsmetrics_protein \
	;";
	CalculateHsMetricsProteinId=`echo "$cmd" | $qsub_lomem -N $name`;
	echo "$CalculateHsMetricsProteinId $name $cmd";

	#BamIndexStats
	name="BamIndexStats_`basename $indexstats`";
	cmd=" \
	java -Xmx5g -Djava.io.tmpdir=$TMPSCRATCH \
	 -jar $PICARD/BamIndexStats.jar \
	 TMP_DIR=$TMPSCRATCH \
	 VERBOSITY=INFO \
	 QUIET=false \
	 VALIDATION_STRINGENCY=SILENT \
	 COMPRESSION_LEVEL=5 \
	 MAX_RECORDS_IN_RAM=1000000 \
	 INPUT=$bam \
	 > $indexstats \
	;";
	BamIndexStatsId=`echo "$cmd" | $qsub_lomem -N $name`;
	echo "$BamIndexStatsId $name $cmd";
	
	#CollectMultipleMetrics
	name="CollectMultipleMetrics_`basename $statsbase`";
	cmd=" \
	java -Xmx5g -Djava.io.tmpdir=$TMPSCRATCH \
	 -jar $PICARD/CollectMultipleMetrics.jar \
	 TMP_DIR=$TMPSCRATCH \
	 VERBOSITY=INFO \
	 QUIET=false \
	 VALIDATION_STRINGENCY=SILENT \
	 COMPRESSION_LEVEL=5 \
	 MAX_RECORDS_IN_RAM=1000000 \
	 REFERENCE_SEQUENCE=$ref \
	 INPUT=$bam \
	 OUTPUT=$statsbase \
	 PROGRAM=CollectAlignmentSummaryMetrics \
	 PROGRAM=CollectInsertSizeMetrics \
	 PROGRAM=QualityScoreDistribution \
	 PROGRAM=MeanQualityByCycle \
	;";
	CollectMultipleMetricsId=`echo "$cmd" | $qsub_lomem -N $name`;
	echo "$CollectMultipleMetricsId $name $cmd";

	#CollectGcBiasMetrics
	name="CollectGcBiasMetrics_`basename $gcbiasmetrics`";
	cmd=" \
	java -Xmx5g -Djava.io.tmpdir=$TMPSCRATCH \
	 -jar $PICARD/CollectGcBiasMetrics.jar \
	 TMP_DIR=$TMPSCRATCH \
	 VERBOSITY=INFO \
	 QUIET=false \
	 VALIDATION_STRINGENCY=SILENT \
	 COMPRESSION_LEVEL=5 \
	 MAX_RECORDS_IN_RAM=1000000 \
	 REFERENCE_SEQUENCE=$ref \
	 INPUT=$bam \
	 OUTPUT=$gcbiasmetrics \
	 CHART_OUTPUT=$gcbiasmetricschart \
	 SUMMARY_OUTPUT=$gcbiasmetricssummary \
	;";
	CollectGcBiasMetricsId=`echo "$cmd" | $qsub_lomem -N $name`;
	echo "$CollectGcBiasMetricsId $name $cmd";

	#EstimateLibraryComplexity
	#SOLiD PE READ_NAME_REGEX=([0-9]+)_([0-9]+)_([0-9]+).*
	#SOLiD PE OPTICAL_DUPLICATE_PIXEL_DISTANCE=10
	#Illumina PE READ_NAME_REGEX=.+:.+:[0-9]+:[0-9]+:([0-9]+):([0-9]+):([0-9]+).*
	
	#Illumina QSEQ.novoalign PE READ_NAME_REGEX=
	#HWI-ST430_243_6_7_10665_92077_0 163
	#[^_]+_[^_]+_[0-9]+_([0-9]+)_([0-9]+)_([0-9]+).*
	
	#Illumina PE OPTICAL_DUPLICATE_PIXEL_DISTANCE=100
	name="EstimateLibraryComplexity_`basename $librarycomplexity`";
	cmd=" \
	java -Xmx5g -Djava.io.tmpdir=$TMPSCRATCH \
	 -jar $PICARD/EstimateLibraryComplexity.jar \
	 TMP_DIR=$TMPSCRATCH \
	 VERBOSITY=INFO \
	 QUIET=false \
	 VALIDATION_STRINGENCY=SILENT \
	 COMPRESSION_LEVEL=5 \
	 MAX_RECORDS_IN_RAM=1000000 \
	 INPUT=$bam \
	 OUTPUT=$librarycomplexity \
	 READ_NAME_REGEX=\"[^_]+_[^_]+_[0-9]+_([0-9]+)_([0-9]+)_([0-9]+).*\" \
	 OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 \
	;";
	EstimateLibraryComplexityId=`echo "$cmd" | $qsub_lomem -N $name`;
	echo "$EstimateLibraryComplexityId $name $cmd";

# included in CollectMultipleMetrics
	#CollectAlignmentSummaryMetrics
	name="CollectAlignmentSummaryMetrics_`basename $alignmentsummarymetrics`";
	cmd=" \
	java -Xmx5g -Djava.io.tmpdir=\"$TMPSCRATCH\" \
	 -jar \"$PICARD/CollectAlignmentSummaryMetrics.jar\" \
	 TMP_DIR=\"$TMPSCRATCH\" \
	 VERBOSITY=INFO \
	 QUIET=false \
	 VALIDATION_STRINGENCY=SILENT \
	 COMPRESSION_LEVEL=5 \
	 MAX_RECORDS_IN_RAM=1000000 \
	 REFERENCE_SEQUENCE=\"$ref\" \
	 INPUT=\"$bam\" \
	 OUTPUT=\"$alignmentsummarymetrics\" \
	 ADAPTER_SEQUENCE=null \
	;";
	#CollectAlignmentSummaryMetricsId=`echo "$cmd" | $qsub_lomem -N $name`;
	#echo "$CollectAlignmentSummaryMetricsId $name $cmd";

	#CollectInsertSizeMetrics
	name="CollectInsertSizeMetrics_`basename $insertsizemetrics`";
	cmd=" \
	java -Xmx5g -Djava.io.tmpdir=\"$TMPSCRATCH\" \
	 -jar \"$PICARD/CollectInsertSizeMetrics.jar\" \
	 TMP_DIR=\"$TMPSCRATCH\" \
	 VERBOSITY=INFO \
	 QUIET=false \
	 VALIDATION_STRINGENCY=SILENT \
	 COMPRESSION_LEVEL=5 \
	 MAX_RECORDS_IN_RAM=1000000 \
	 INPUT=\"$bam\" \
	 OUTPUT=\"$insertsizemetrics\" \
	 HISTOGRAM_FILE=\"$insertsizehistogram\" \
	;";
	#CollectInsertSizeMetricsId=`echo "$cmd" | $qsub_lomem -N $name`;
	#echo "$CollectInsertSizeMetricsId $name $cmd";

	#MeanQualityByCycle
	name="MeanQualityByCycle_`basename $meanqualitybycycle`";
	cmd=" \
	java -Xmx5g -Djava.io.tmpdir=\"$TMPSCRATCH\" \
	 -jar \"$PICARD/MeanQualityByCycle.jar\" \
	 TMP_DIR=\"$TMPSCRATCH\" \
	 VERBOSITY=INFO \
	 QUIET=false \
	 VALIDATION_STRINGENCY=SILENT \
	 COMPRESSION_LEVEL=5 \
	 MAX_RECORDS_IN_RAM=1000000 \
	 INPUT=\"$bam\" \
	 OUTPUT=\"$meanqualitybycycle\" \
	 CHART_OUTPUT=\"$meanqualitybycyclechart\" \
	;";
	#MeanQualityByCycleId=`echo "$cmd" | $qsub_lomem -N $name`;
	#echo "$MeanQualityByCycleId $name $cmd";

	#QualityScoreDistribution
	name="QualityScoreDistribution_`basename $qualityscoredistribution`";
	cmd=" \
	java -Xmx5g -Djava.io.tmpdir=\"$TMPSCRATCH\" \
	 -jar \"$PICARD/QualityScoreDistribution.jar\" \
	 TMP_DIR=\"$TMPSCRATCH\" \
	 VERBOSITY=INFO \
	 QUIET=false \
	 VALIDATION_STRINGENCY=SILENT \
	 COMPRESSION_LEVEL=5 \
	 MAX_RECORDS_IN_RAM=1000000 \
	 INPUT=\"$bam\" \
	 OUTPUT=\"$qualityscoredistribution\" \
	 CHART_OUTPUT=\"$qualityscoredistributionchart\" \
	;";
	#QualityScoreDistributionId=`echo "$cmd" | $qsub_lomem -N $name`;
	#echo "$QualityScoreDistributionId $name $cmd";
done
cd $wd;

