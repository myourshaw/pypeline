#!/bin/bash

#gatk_DepthOfCoverage_q.sh

usage="usage: $0 [<email address> (default: myourshaw@ucla.edu)] [<ref.fasta> | <ref.fa> (default human_g1k_v37.fasta)] [<target intervals> (default SureSelect_All_Exon_50mb_with_annotation_b37_sorted.interval_list)] <path to output file> <list of sorted bam filest>";
#get path and parent directory of this script
#http://hintsforums.macworld.com/archive/index.php/t-73839.html
# if the call came via symlink, then use its target instead:
arg=$0; [[ -L $0 ]] && arg=$(stat -L -c'%n' "$0")
script_path=$(2>/dev/null cd "${arg%/*}" >&2; echo "`pwd -P`/${arg##*/}")
script_dir=$(dirname "$script_path")

EMAIL="myourshaw@ucla.edu";
if [[ "$1" =~ '\S+@\S+\.\S+'  ]]; then EMAIL="$1"; shift; fi

source ${script_dir}/b37.variables.sh $EMAIL;

ref="$B37_FASTA";
if [[ ${1##*.} == "fasta" || ${1##*.} == "fa" ]]; then ref="$1"; shift; fi
targets="$SURESELECT_50MB_INTERVAL_LIST";
if [[ ${1##*.} == "bed" || ${1##*.} == "interval_list" || ${1##*.} == "intervals" ]]; then targets="$1"; shift; fi
if (( $# < 2 )); then echo $usage; exit; fi

coverage_output_file="$1"; shift;
output_dir="`dirname $coverage_output_file`";
mkdir -p "$output_dir";
qout="$output_dir/qout";
mkdir -p "$qout";
qsub_allq="$QSUB_ALLQ -e $qout -o $qout";
qsub_lomem="$QSUB_LOMEM -e $qout -o $qout";
qsub_himem="$QSUB_HIMEM -e $qout -o $qout";

bamlist='';
for bam in $@; do
	if [[ ! -e $bam ]]; then echo "bam file does not exist [$bam]"; exit; fi
	bamlist="${bamlist} -I ${bam}";
done

wd="`pwd`";

#resources
dbsnp="$DBSNP_132_B37_VCF";
refgene="$REFGENE_B37_ROD";

#applications and resources
sting="$STING";
gatk="$GATK";
gatkrscripts="$GATK_R_SCRIPTS";
perl="$GATK_PERL";
python="$GATK_PYTHON";
rscript="$RSCRIPT";
picard="$PICARD";
tmplocal="$TMPLOCAL";
tmpscratch="$TMPSCRATCH";
mkdir -p "$tmpscratch";

#DepthOfCoverage
name="gatkDepthOfCoverage_`basename $coverage_output_file`";
cmd=" \
java -Xmx16g -Djava.io.tmpdir=\"$TMPSCRATCH\" \
 -jar \"$GATK/GenomeAnalysisTK.jar\" \
 -T DepthOfCoverage \
 -l INFO \
 -R \"$ref\" \
 $bamlist \
 -o \"$coverage_output_file\" \
 -L \"$targets\" \
 -geneList \"${REFGENE_B37_SORTED}\" \
 --omitDepthOutputAtEachBase \
 -pt sample -pt readgroup -pt library \
 -ct 1 -ct 2 -ct 3 -ct 4 -ct 5 -ct 6 -ct 7 -ct 8 -ct 9 -ct 10 -ct 20 -ct 30 -ct 40 -ct 50 -ct 60 -ct 70 -ct 80 -ct 90 -ct 100 -ct 110 -ct 120 -ct 130 -ct 140 -ct 150 -ct 160 -ct 170 -ct 180 -ct 190 -ct 200 -ct 300 \
;";
DepthOfCoverageId=`echo "$cmd" | $qsub_himem -N $name`;
echo "$DepthOfCoverageId $name $cmd";

cd "$wd";

