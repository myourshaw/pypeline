#!/bin/bash

usage='usage: '$0' \\\n
[-m <email address> (default: '`whoami`'@ucla.edu)] \\\n
[-R <fasta reference genome> (default: b37.variables.sh->B37_FASTA)] \\\n
[-G <GATK DepthOfCoverage geneList> (default: b37.variables.sh->REFGENE_B37_SORTED)] \\\n
[-L <GATK DepthOfCoverage target interval list> (default: bait interval list)] \\\n
[-B <bait interval list> (default: b37.variables.sh->SURESELECT_50MB_INTERVAL_LIST)] \\\n
[-A <target interval list 1> (default: all genes from b37.variables.sh->ENSEMBL_ALL_GENES)] \\\n
[-P <target interval list 2> (default: protein coding genes from b37.variables.sh->ENSEMBL_PROTEIN_CODING_GENES)] \\\n
[-T <target interval list 3> (default: none)] \\\n
[-X <read name regex> (default for illumina/novoalign: "[^_]+_[^_]+_[0-9]+_([0-9]+)_([0-9]+)_([0-9]+).*")] \\\n
[-D <optical duplicate pixel distance> (default for illumina: 100)] \\\n
[-O <directory for output>(default: <directory of first input file>/picard_metrics)] \\\n
<input list of bam or sam files> \n'

wd=`pwd`;
#get path and parent directory of this script
#http://hintsforums.macworld.com/archive/index.php/t-73839.html
# if the call came via symlink, then use its target instead:
arg=$0; [[ -L $0 ]] && arg=$(stat -L -c'%n' "$0")
script_path=$(2>/dev/null cd "${arg%/*}" >&2; echo "`pwd -P`/${arg##*/}")
script_dir=$(dirname "$script_path")
HERE=${script_dir};

#user options
while getopts  "m:R:G:L:B:A:P:T:X:D:O:" o
do	case "$o" in
	m)	EMAIL="$OPTARG";;
	R)	REF="$OPTARG";;
	G)	GENELIST="$OPTARG";;
	L)	TARGETLIMIT="$OPTARG";;
	B)	BAIT="$OPTARG";;
	A)	TARGET1="$OPTARG";;
	P)	TARGET2="$OPTARG";;
	T)	TARGET3="$OPTARG";;
	X)	READREGEX="$OPTARG";;
	D)	PIXELDISTANCE="$OPTARG";;
	O)	METRICSDIRECTORY="$OPTARG";;
	[?])	echo -e $usage
		exit 1;;
	esac
done
shift $(( $OPTIND - 1 ));

#required parameter(s) and supporting file(s)
if (( $# < 1 )); then echo -e $usage; exit; fi
B37_VARIABLES="${HERE}/b37.variables.sh";
if [[ ! -e "${B37_VARIABLES}" ]]; then echo "this script requires ${B37_VARIABLES}"; exit; fi
if [[ $EMAIL == '' ]]; then EMAIL="`whoami`@ucla.edu"; fi
source ${B37_VARIABLES} $EMAIL;

#user options | defaults
if [[ ! -e "${REF}" ]]; then REF="${B37_FASTA}"; fi
if [[ ! -e "${BAIT}" ]]; then BAIT="${SURESELECT_50MB_INTERVAL_LIST}"; fi
if [[ ! -e "${GENELIST}" ]]; then GENELIST="${REFGENE_B37_SORTED}"; fi
if [[ ! -e "${TARGETLIMIT}" ]]; then TARGETLIMIT="${BAIT}"; fi
if [[ ! -e "${TARGET1}" ]]; then TARGET1="${ENSEMBL_ALL_GENES}"; fi
if [[ ! -e "${TARGET2}" ]]; then TARGET2="${ENSEMBL_PROTEIN_CODING_GENES}"; fi
if [[ "${READREGEX}" == '' ]]; then READREGEX='"[^_]+_[^_]+_[0-9]+_([0-9]+)_([0-9]+)_([0-9]+).*"'; fi
if [[ "${PIXELDISTANCE}" == '' ]]; then PIXELDISTANCE=100; fi
if [[ "${METRICSDIRECTORY}" == '' ]]; then METRICSDIRECTORY="`dirname $1`/picard_metrics"; fi

#check existence of required files
for e in \
 "${REF}" \
 "${BAIT}" \
 "${GENELIST}" \
 "${TARGETLIMIT}" \
 "${SAMTOOLS}" \
 "${PICARD}/CalculateHsMetrics.jar" \
 "${PICARD}/BamIndexStats.jar" \
 "${PICARD}/CollectMultipleMetrics.jar" \
 "${PICARD}/CollectGcBiasMetrics.jar" \
 "${PICARD}/EstimateLibraryComplexity.jar" \
 "${PICARD_METRICS_CONSOLIDATE}" \
 ; do
	if [[ ! -e "${e}" ]]; then echo "this script requires ${e}"; exit; fi
done

#check existence of bam/sam files
for b in $@; do
 if [[ ! -e $b ]]; then echo "file does not exist [$b]"; exit; fi
done

mkdir -p "${METRICSDIRECTORY}";
cd "${METRICSDIRECTORY}";
qout="${METRICSDIRECTORY}/qout";
mkdir -p "${qout}";
consolidated_dir="${METRICSDIRECTORY}/consolidated_metrics";
coverage_output_file="${METRICSDIRECTORY}/DepthOfCoverage"
mkdir -p "${consolidated_dir}";
qsub_allq="${QSUB_ALLQ} -e ${qout} -o ${qout}";
qsub_lomem="${QSUB_LOMEM} -e ${qout} -o ${qout}";
qsub_himem="${QSUB_HIMEM} -e ${qout} -o ${qout}";

#job control lists
holdall='';
files_to_consolidate='';

for bam in $@; do
	#output files
	statsbase="${METRICSDIRECTORY}/`basename ${bam}`.picard_metrics";
	validate="${statsbase}.validate";
	idxstats="${statsbase}.idxstats";
	if [[ $TARGET1 != '' ]]; then hsmetrics_target1="${statsbase}.target_`basename ${TARGET1}`.hsmetrics"; fi
	if [[ $TARGET2 != '' ]]; then hsmetrics_target2="${statsbase}.target_`basename ${TARGET2}`.hsmetrics"; fi
	if [[ $TARGET3 != '' ]]; then hsmetrics_target3="${statsbase}.target_`basename ${TARGET3}`.hsmetrics"; fi
	indexstats="${statsbase}.indexstats";
	gcbiasmetrics="${statsbase}.gcbias.table";
	gcbiasmetricschart="${statsbase}.gcbias.pdf";
	gcbiasmetricssummary="${statsbase}.gcbias.summary";
	librarycomplexity="${statsbase}.librarycomplexity";
	
	#idxstats
	name="idxstats_`basename ${idxstats}`";
	cmd=" \
	${SAMTOOLS} idxstats ${bam} > ${idxstats} \
	;";
	idxstatsId=`echo "${cmd}" | ${qsub_allq} -N ${name}`;
	echo "${idxstatsId} ${name} ${cmd}";
	if [[ ${holdall} == '' ]]; then holdall=${idxstatsId}; else holdall="${holdall},${idxstatsId}"; fi

	#CalculateHsMetrics.target1
	if [[ -e "${BAIT}" && -e "${TARGET1}" ]]; then
		name="CalculateHsMetrics_`basename ${hsmetrics_target1}`";
		cmd=" \
		java -Xmx5g -Djava.io.tmpdir=${TMPSCRATCH} \
		 -jar ${PICARD}/CalculateHsMetrics.jar \
		 TMP_DIR=${TMPSCRATCH} \
		 VERBOSITY=INFO \
		 QUIET=false \
		 VALIDATION_STRINGENCY=SILENT \
		 COMPRESSION_LEVEL=5 \
		 MAX_RECORDS_IN_RAM=1000000 \
		 BAIT_INTERVALS=${BAIT} \
		 TARGET_INTERVALS=${TARGET1} \
		 INPUT=${bam} \
		 OUTPUT=${hsmetrics_target1} \
		;";
		CalculateHsMetricsAllId=`echo "${cmd}" | ${qsub_allq} -N ${name}`;
		echo "$CalculateHsMetricsAllId ${name} ${cmd}";
		if [[ ${holdall} == '' ]]; then holdall=${CalculateHsMetricsAllId}; else holdall="${holdall},${CalculateHsMetricsAllId}"; fi
		if [[ ${files_to_consolidate} == '' ]]; then files_to_consolidate="${hsmetrics_target1}"; else files_to_consolidate="${files_to_consolidate} ${hsmetrics_target1}"; fi
	fi
	
	#CalculateHsMetrics.target2
	if [[ -e "${BAIT}" && -e "${TARGET2}" ]]; then
		name="CalculateHsMetrics_`basename ${hsmetrics_target2}`";
		cmd=" \
		java -Xmx5g -Djava.io.tmpdir=${TMPSCRATCH} \
		 -jar ${PICARD}/CalculateHsMetrics.jar \
		 TMP_DIR=${TMPSCRATCH} \
		 VERBOSITY=INFO \
		 QUIET=false \
		 VALIDATION_STRINGENCY=SILENT \
		 COMPRESSION_LEVEL=5 \
		 MAX_RECORDS_IN_RAM=1000000 \
		 BAIT_INTERVALS=${BAIT} \
		 TARGET_INTERVALS=${TARGET2} \
		 INPUT=${bam} \
		 OUTPUT=${hsmetrics_target2} \
		;";
		CalculateHsMetricsProteinId=`echo "${cmd}" | ${qsub_allq} -N ${name}`;
		echo "$CalculateHsMetricsProteinId ${name} ${cmd}";
		if [[ ${holdall} == '' ]]; then holdall=${CalculateHsMetricsProteinId}; else holdall="${holdall},${CalculateHsMetricsProteinId}"; fi
		if [[ ${files_to_consolidate} == '' ]]; then files_to_consolidate="${hsmetrics_target2}"; else files_to_consolidate="${files_to_consolidate} ${hsmetrics_target2}"; fi
	fi
	
	#CalculateHsMetrics.target3
	if [[ -e "${BAIT}" && -e "${TARGET3}" ]]; then
		name="CalculateHsMetrics_`basename ${hsmetrics_target3}`";
		cmd=" \
		java -Xmx5g -Djava.io.tmpdir=${TMPSCRATCH} \
		 -jar ${PICARD}/CalculateHsMetrics.jar \
		 TMP_DIR=${TMPSCRATCH} \
		 VERBOSITY=INFO \
		 QUIET=false \
		 VALIDATION_STRINGENCY=SILENT \
		 COMPRESSION_LEVEL=5 \
		 MAX_RECORDS_IN_RAM=1000000 \
		 BAIT_INTERVALS=${BAIT} \
		 TARGET_INTERVALS=${TARGET3} \
		 INPUT=${bam} \
		 OUTPUT=${hsmetrics_target3} \
		;";
		CalculateHsMetricsProteinId=`echo "${cmd}" | ${qsub_allq} -N ${name}`;
		echo "$CalculateHsMetricsProteinId ${name} ${cmd}";
		if [[ ${holdall} == '' ]]; then holdall=${CalculateHsMetricsProteinId}; else holdall="${holdall},${CalculateHsMetricsProteinId}"; fi
		if [[ ${files_to_consolidate} == '' ]]; then files_to_consolidate="${hsmetrics_target3}"; else files_to_consolidate="${files_to_consolidate} ${hsmetrics_target3}"; fi
	fi
	
	#BamIndexStats
	name="BamIndexStats_`basename ${indexstats}`";
	cmd=" \
	java -Xmx5g -Djava.io.tmpdir=${TMPSCRATCH} \
	 -jar ${PICARD}/BamIndexStats.jar \
	 TMP_DIR=${TMPSCRATCH} \
	 VERBOSITY=INFO \
	 QUIET=false \
	 VALIDATION_STRINGENCY=SILENT \
	 COMPRESSION_LEVEL=5 \
	 MAX_RECORDS_IN_RAM=1000000 \
	 INPUT=${bam} \
	 > ${indexstats} \
	;";
	BamIndexStatsId=`echo "${cmd}" | ${qsub_allq} -N ${name}`;
	echo "${BamIndexStatsId} ${name} ${cmd}";
	if [[ ${holdall} == '' ]]; then holdall=${BamIndexStatsId}; else holdall="${holdall},${BamIndexStatsId}"; fi
	
	#CollectMultipleMetrics
	name="CollectMultipleMetrics_`basename ${statsbase}`";
	cmd=" \
	java -Xmx5g -Djava.io.tmpdir=${TMPSCRATCH} \
	 -jar ${PICARD}/CollectMultipleMetrics.jar \
	 TMP_DIR=${TMPSCRATCH} \
	 VERBOSITY=INFO \
	 QUIET=false \
	 VALIDATION_STRINGENCY=SILENT \
	 COMPRESSION_LEVEL=5 \
	 MAX_RECORDS_IN_RAM=1000000 \
	 REFERENCE_SEQUENCE=${REF} \
	 INPUT=${bam} \
	 OUTPUT=${statsbase} \
	 PROGRAM=CollectAlignmentSummaryMetrics \
	 PROGRAM=CollectInsertSizeMetrics \
	 PROGRAM=QualityScoreDistribution \
	 PROGRAM=MeanQualityByCycle \
	;";
	CollectMultipleMetricsId=`echo "${cmd}" | ${qsub_allq} -N ${name}`;
	echo "$CollectMultipleMetricsId ${name} ${cmd}";
	if [[ ${holdall} == '' ]]; then holdall=${CollectMultipleMetricsId}; else holdall="${holdall},${CollectMultipleMetricsId}"; fi
	if [[ ${files_to_consolidate} == '' ]]; then files_to_consolidate="${statsbase}.*_metrics"; else files_to_consolidate="${files_to_consolidate} ${statsbase}.*_metrics"; fi

	#CollectGcBiasMetrics
	name="CollectGcBiasMetrics_`basename ${gcbiasmetrics}`";
	cmd=" \
	java -Xmx5g -Djava.io.tmpdir=${TMPSCRATCH} \
	 -jar ${PICARD}/CollectGcBiasMetrics.jar \
	 TMP_DIR=${TMPSCRATCH} \
	 VERBOSITY=INFO \
	 QUIET=false \
	 VALIDATION_STRINGENCY=SILENT \
	 COMPRESSION_LEVEL=5 \
	 MAX_RECORDS_IN_RAM=1000000 \
	 REFERENCE_SEQUENCE=${REF} \
	 INPUT=${bam} \
	 OUTPUT=${gcbiasmetrics} \
	 CHART_OUTPUT=${gcbiasmetricschart} \
	 SUMMARY_OUTPUT=${gcbiasmetricssummary} \
	;";
	CollectGcBiasMetricsId=`echo "${cmd}" | ${qsub_allq} -N ${name}`;
	echo "$CollectGcBiasMetricsId ${name} ${cmd}";
	if [[ ${holdall} == '' ]]; then holdall=${CollectGcBiasMetricsId}; else holdall="${holdall},${CollectGcBiasMetricsId}"; fi
	if [[ ${files_to_consolidate} == '' ]]; then files_to_consolidate="${gcbiasmetrics}"; else files_to_consolidate="${files_to_consolidate} ${gcbiasmetrics}"; fi
	if [[ ${files_to_consolidate} == '' ]]; then files_to_consolidate="${gcbiasmetricssummary}"; else files_to_consolidate="${files_to_consolidate} ${gcbiasmetricssummary}"; fi

	#EstimateLibraryComplexity
	
	#SOLiD PE READ_NAME_REGEX=([0-9]+)_([0-9]+)_([0-9]+).*
	#SOLiD PE OPTICAL_DUPLICATE_PIXEL_DISTANCE=10
	#Illumina PE READ_NAME_REGEX=.+:.+:[0-9]+:[0-9]+:([0-9]+):([0-9]+):([0-9]+).*
	
	#Illumina QSEQ.novoalign PE READ_NAME_REGEX=
	#HWI-ST430_243_6_7_10665_92077_0 163
	#[^_]+_[^_]+_[0-9]+_([0-9]+)_([0-9]+)_([0-9]+).*
	
	#Illumina PE OPTICAL_DUPLICATE_PIXEL_DISTANCE=100
	
	name="EstimateLibraryComplexity_`basename ${librarycomplexity}`";
	cmd=" \
	java -Xmx5g -Djava.io.tmpdir=${TMPSCRATCH} \
	 -jar ${PICARD}/EstimateLibraryComplexity.jar \
	 TMP_DIR=${TMPSCRATCH} \
	 VERBOSITY=INFO \
	 QUIET=false \
	 VALIDATION_STRINGENCY=SILENT \
	 COMPRESSION_LEVEL=5 \
	 MAX_RECORDS_IN_RAM=1000000 \
	 INPUT=${bam} \
	 OUTPUT=${librarycomplexity} \
	 READ_NAME_REGEX=${READREGEX} \
	 OPTICAL_DUPLICATE_PIXEL_DISTANCE=${PIXELDISTANCE} \
	;";
	EstimateLibraryComplexityId=`echo "${cmd}" | ${qsub_allq} -N ${name}`;
	echo "$EstimateLibraryComplexityId ${name} ${cmd}";
	if [[ ${holdall} == '' ]]; then holdall=${EstimateLibraryComplexityId}; else holdall="${holdall},${EstimateLibraryComplexityId}"; fi
	if [[ ${files_to_consolidate} == '' ]]; then files_to_consolidate="${librarycomplexity}"; else files_to_consolidate="${files_to_consolidate} ${librarycomplexity}"; fi
done

#DepthOfCoverage
bamlist='';
for bam in $@; do
	bamlist="${bamlist} -I ${bam}";
done

name="DepthOfCoverage";
cmd=" \
java -Xmx16g -Djava.io.tmpdir=\"${TMPSCRATCH}\" \
 -jar \"${GATK}/GenomeAnalysisTK.jar\" \
 -T DepthOfCoverage \
 -l INFO \
 -R \"${REF}\" \
 ${bamlist} \
 -o \"$coverage_output_file\" \
 -L \"${TARGETLIMIT}\" \
 -geneList \"${GENELIST}\" \
 --omitDepthOutputAtEachBase \
 -dcov 1000 \
 -pt sample -pt readgroup -pt library \
 -ct 1 -ct 2 -ct 3 -ct 4 -ct 5 -ct 6 -ct 7 -ct 8 -ct 9 -ct 10 -ct 20 -ct 30 -ct 40 -ct 50 -ct 60 -ct 70 -ct 80 -ct 90 -ct 100 -ct 110 -ct 120 -ct 130 -ct 140 -ct 150 -ct 160 -ct 170 -ct 180 -ct 190 -ct 200 -ct 300 \
;";
#DepthOfCoverageId=`echo "${cmd}" | ${qsub_himem} -N ${name}`;
#echo "${DepthOfCoverageId} ${name} ${cmd}";

#ConsolidateMetrics
name="ConsolidateMetrics";
cmd=" \
${PYTHON} ${PICARD_METRICS_CONSOLIDATE} \
-d ${consolidated_dir} \
${files_to_consolidate} \
;";
ConsolidateMetricsId=`echo "${cmd}" | ${qsub_allq} -N ${name} -hold_jid ${holdall}`;
echo "${ConsolidateMetricsId} ${name} ${cmd}";

cd $wd;
