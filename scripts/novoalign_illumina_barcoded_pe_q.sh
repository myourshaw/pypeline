#!/bin/bash

usage='usage: '$0' \\\n
[-m <email address> (default: '`whoami`'@ucla.edu)] \\\n
[-r <read group identier> (default: from file name, example: "GMD108_XT_Illumina.2011-05-11.HWI-ST430.243.6.ACTTGAA"] \\\n
[-d <date the run was produced (ISO8601 date or date/time)> (default: from file name, example "2011-05-11")] \\\n
[-l <library> (default: from file name, example: "GMD108_XT")] \\\n
[-u <platform unit> (default: from file name, example: "HWI-ST430.243.6.ACTTGAA")] \\\n
[-s <sample> (default: from file name, example: "GMD108")] \\\n
[-1 <multiplexing adapter 1> (default: Agilent sureSelect XT="AGATCGGAAGAGCACACGTCT")] \\\n
[-2 <reverse complement of multiplexing adapter 2 + A 2> (default: Agilent sureSelect XT="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA")] \\\n
[-x <novoindex> (default: if -h not specified then b37.variables.sh->NOVOALIGN_LOMEM_INDEX else b37.variables.sh->NOVOALIGN_MAX_INDEX)] \\\n
[-h (if present, use high memory node)] \\\n
[-R <fasta reference genome> (default: b37.variables.sh->37_FASTA)] \\\n
[-B <bait interval list> (default: b37.variables.sh->SURESELECT_50MB_INTERVAL_LIST)] \\\n
[-A <target interval list 1> (default: all genes from b37.variables.sh->ENSEMBL_ALL_GENES)] \\\n
[-P <target interval list 2> (default: protein coding genes from b37.variables.sh->ENSEMBL_PROTEIN_CODING_GENES)] \\\n
[-T <target interval list 3> (default: none)] \\\n
[-X <read name regex> (default for illumina/novoalign: "[^_]+_[^_]+_[0-9]+_([0-9]+)_([0-9]+)_([0-9]+).*")] \\\n
[-D <optical duplicate pixel distance> (default for illumina: 100)] \\\n
[-o <directory for output of all aligned, sorted, validated bam files> (default: directory of each input file)] \\\n
[-O <directory for picard metrics> (default: <output directory>/picard_metrics)] \\\n
<list of qseq files for paired-end 1 named as "RG_ID.pe{1,2}.qseq.txt"> \n
where RG_ID = <experiment>.<date>.<machine>.<run>.<lane>.<barcode> \n
and experiment = <sample>_<library>_<platform> \n
example of filename "GMD108_XT_Illumina.2011-05-11.HWI-ST430.243.6.ACTTGAA.pe1.qseq.txt" \n'

wd=`pwd`;
#get path and parent directory of this script
#http://hintsforums.macworld.com/archive/index.php/t-73839.html
# if the call came via symlink, then use its target instead:
arg=$0; [[ -L $0 ]] && arg=$(stat -L -c'%n' "$0")
script_path=$(2>/dev/null cd "${arg%/*}" >&2; echo "`pwd -P`/${arg##*/}")
script_dir=$(dirname "$script_path")
HERE=${script_dir};

#user options
HIMEM=0;
while getopts  "m:r:d:l:u:s:1:2:x:hR:B:A:P:T:X:D:o:O:" o
do	case "$o" in
	m)	EMAIL="$OPTARG";;
	r)	READGROUP="$OPTARG";;
	d)	DATE="$OPTARG";;
	l)	LIBRARY="$OPTARG";;
	u)	PLATFORMUNIT="$OPTARG";;
	s)	SAMPLE="$OPTARG";;
	1)	ADAPTER1="$OPTARG";;
	2)	ADAPTER2="$OPTARG";;
	x)	NOVOINDEX="$OPTARG";;
	h)	HIMEM=1;;
	R)	REF="$OPTARG";;
	B)	BAIT="$OPTARG";;
	A)	TARGET1="$OPTARG";;
	P)	TARGET2="$OPTARG";;
	X)	READREGEX="$OPTARG";;
	D)	PIXELDISTANCE="$OPTARG";;
	o)	OUTPUTDIRECTORY="$OPTARG";;
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

#user options | defaults
if [[ ${EMAIL == '' ]]; then EMAIL="`whoami`@ucla.edu"; fi
source ${B37_VARIABLES} $EMAIL;
if [[ ! -x "${NOVOALIGN}" ]]; then echo "this script requires ${NOVOALIGN}"; exit; fi
if [[ ! -x "${PICARD_METRICS}" ]]; then echo "this script requires ${PICARD_METRICS}"; exit; fi
if [[ ! -e "${PICARD_METRICS_CONSOLIDATE}" ]]; then echo "this script requires ${PICARD_METRICS_CONSOLIDATE}"; exit; fi
if [[ "${ADAPTER1}" == '' ]]; then ADAPTER1="AGATCGGAAGAGCACACGTCT"; fi
if [[ "${ADAPTER2}" == '' ]]; then ADAPTER2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA"; fi
if [[ "${NOVOINDEX}" == '' ]]; then if [[ HIMEM == 0 ]] then NOVOINDEX="${NOVOALIGN_LOMEM_INDEX}" else NOVOINDEX="${NOVOALIGN_MAX_INDEX}"; fi
if [[ "${REF}" == '' ]]; then REF="${B37_FASTA}"; fi
if [[ "${BAIT}" == '' ]]; then BAIT="${SURESELECT_50MB_INTERVAL_LIST}"; fi
if [[ "${TARGET1}" == '' ]]; then TARGET1="${ENSEMBL_ALL_GENES}"; TARGET2="${ENSEMBL_PROTEIN_CODING_GENES}"; fi
if [[ "${TARGET2}" == '' ]]; then TARGET2="${ENSEMBL_PROTEIN_CODING_GENES}"; fi
if [[ "${READREGEX}" == '' ]]; then READREGEX='"[^_]+_[^_]+_[0-9]+_([0-9]+)_([0-9]+)_([0-9]+).*"'; fi
if [[ "${PIXELDISTANCE}" == '' ]]; then PIXELDISTANCE=100; fi
if [[ "${OUTPUTDIRECTORY}" == '' ]]; then OUTPUTDIRECTORY="`dirname $1`/picard_metrics"; fi
if [[ "${METRICSDIRECTORY}" == '' ]]; then METRICSDIRECTORY="`dirname $1`/picard_metrics"; fi

for f in $@; do
 if [[ ! -e $f ]]; then echo "file does not exist [$f]"; exit; fi
done

#########################################################################################################

#job control lists
holdall='';
bam_list='';

for f in $@; do

	if [[ "${OUTPUTDIRECTORY}" == '' ]]; then output_dir="`dirname ${f}`"; else output_dir="${OUTPUTDIRECTORY}"; fi
	if [[ "${METRICSDIRECTORY}" == '' ]]; then metrics_dir="${OUTPUTDIRECTORY}/picard_metrics"; else metrics_dir="${METRICSDIRECTORY}"; fi
	mkdir -p "${output_dir}";
	mkdir -p "${metrics_dir}";
	qout="${output_dir}/qout";
	mkdir -p "${qout}";
	qsub_allq="${QSUB_ALLQ} -e ${qout} -o ${qout}";
	qsub_allq_2gb="${QSUB_ALLQ_2GB} -e ${qout} -o ${qout}";
	qsub_allq_4gb="${QSUB_ALLQ_2GB} -e ${qout} -o ${qout}";
	qsub_lomem="${QSUB_LOMEM} -e ${qout} -o ${qout}";
	qsub_himem="${QSUB_HIMEM} -e ${qout} -o ${qout}";

	f2="${f%.pe1.qseq.txt}.pe2.qseq.txt";
	if [[ ! -e "${f2}" ]]; then f2=""; fi
	base=`basename ${f}`; #GMD108_XT_Illumina.2011-05-11.HWI-ST430.243.6.ACTTGAA.pe1.qseq.txt

	prefix="${f%.pe1.qseq.txt}";
	
	#files created by this script
	sam="${prefix}.novoalign.sam";
	bam="$[prefix}.novoalign.bam";
	validate="${bam}.validate";
	
	#Novoalign
	if [[ "${READGROUP}" == '' ]]; then rgid="${base%.pe1.qseq.txt}"; else rgid="${READGROUP}"; fi #GMD108_XT_Illumina.2011-05-11.HWI-ST430.243.6.ACTTGAA
	experiment="${rgid%%.*}"; #GMD108_XT_Illumina
	if [[ "${SAMPLE}" == '' ]]; then sample="${experiment%%_*}"; else sample="${SAMPLE}"; fi #GMD108
	if [[ "${LIBRARY}" == '' ]]; then library="${experiment%_*}"; else library="${LIBRARY}"; #GMD108_XT
	x="${rgid#*.}"; #2011-05-11.HWI-ST430.243.6.ACTTGAA
	if [[ "${DATE}" == '' ]]; then dt="${x%%.*}"; else dt="${DATE}"; #2011-05-11
	if [[ "${PLATFORMUNIT}" == '' ]]; then pu="${x#*.}"; else pu="${PLATFORMUNIT}"; #HWI-ST430.243.6.ACTTGAA

	rg='@RG\tID:'${rgid}'\tCN:UCLA\tDT:'${dt}'\tLB:'${library}'\tPL:ILLUMINA\tPU:'${pu}'\tSM:'${sample};
	#@RG\tID:GMD108_XT_Illumina.2011-05-11.HWI-ST430.243.6.ACTTGAA\tCN:UCLA\tDT:2011-05-11\tLB:GMD108_XT\tPL:ILLUMINA\tPU:HWI-ST430.243.6.ACTTGAA\tSM:GMD108
	
	name=novoalign_`basename ${sam}`;
	cmd="${NOVOALIGN} -k -o SAM \"$rg\" -d ${NOVOINDEX} -a ${ADAPTER1} ${ADAPTER2} -F QSEQ -f ${f} ${f2} > ${sam};";
	if [[ ${HIMEM} == 0 ]]; then NovoalignId=`echo "${cmd}" | ${qsub_lomem} -N ${name}` else NovoalignId=`echo "${cmd}" | ${qsub_himem} -N ${name}`; fi
	echo "${NovoalignId} ${name} $[cmd}";

	#FixMateInformation
	name=gatkFixMateInformation_`basename ${bam}`;
	cmd=" \
	java -Xmx5g -Djava.io.tmpdir=${TMPSCRATCH} \
	 -jar ${PICARD}/FixMateInformation.jar \
	 TMP_DIR=${TMPSCRATCH} \
	 VERBOSITY=INFO \
	 QUIET=false \
	 VALIDATION_STRINGENCY=STRICT \
	 COMPRESSION_LEVEL=5 \
	 MAX_RECORDS_IN_RAM=1000000 \
	 CREATE_INDEX=true \
	 CREATE_MD5_FILE=true \
	 INPUT=${sam} \
	 OUTPUT=${bam} \
	 SORT_ORDER=coordinate \
	;";
	FixMateInformationId=`echo "${cmd}" | ${qsub_lomem} -N ${name} -hold_jid ${NovoalignId}`;
	echo "${FixMateInformationId} ${name} ${cmd}";
	
	#ValidateSamFile
	name=gatkValidateSamFile_`basename ${validate}`;
	cmd=" \
	java -Xmx5g  -Djava.io.tmpdir=${TMPSCRATCH} \
	 -jar ${PICARD}/ValidateSamFile.jar \
 	 TMP_DIR=${TMPSCRATCH} \
	 VERBOSITY=INFO \
	 QUIET=false \
	 VALIDATION_STRINGENCY=STRICT \
	 COMPRESSION_LEVEL=5 \
	 MAX_RECORDS_IN_RAM=1000000 \
	 INPUT=${bam} \
	 OUTPUT=${validate} \
	 MODE=VERBOSE \
	 REFERENCE_SEQUENCE=${B37_FASTA} \
	 ;";
	ValidateSamFileId=`echo "${cmd}" | ${qsub_lomem} -N ${name} -hold_jid ${FixMateInformationId}`;
	echo "${ValidateSamFileId} ${name} ${cmd}";
	if [[ ${holdall} == '' ]]; then holdall=${ValidateSamFileId}; else holdall="${holdall},${ValidateSamFileId}"; fi
	if [[ ${bam_list} == '' ]]; then bam_list="${bam}"; else bam_list="${bam_list} ${bam}"; fi

done

#PicardMetricsPreRmdup
name="PicardMetricsPreRmdup_pre-rmdup" 
cmd=" \
${PICARD_METRICS} \
-m ${EMAIL} \
-R ${REF} \
-B ${BAIT} \
-A ${TARGET1} \
-P ${TARGET2} \
-T ${TARGET3} \
-X ${READREGEX} \
-D ${PIXELDISTANCE} \
-O ${METRICSDIRECTORY}/pre-rmdup \
${bam_list} \
;";
ConsolidateMetricsId=`echo "${cmd}" | ${qsub_allq} -N ${name} -hold_jid ${holdall}`;
echo "$ConsolidateMetricsId ${name} ${cmd}";

#PicardMetricsPostRmdup
cmd=" \
${PICARD_METRICS} \
-m ${EMAIL} \
-R ${REF} \
-B ${BAIT} \
-A ${TARGET1} \
-P ${TARGET2} \
-T ${TARGET3} \
-X ${READREGEX} \
-D ${PIXELDISTANCE} \
-O ${METRICSDIRECTORY}/post-rmdup \
${rmdup_list} \
;";

cd $wd;

