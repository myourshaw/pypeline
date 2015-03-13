#!/bin/bash
#$ -cwd
#$ -V
#$ -M myourshaw@ucla.edu
#$ -m eas
#$ -terse

#gatk_ReadBackedPhasing_q.sh

usage="usage: $0 [<email address> (default: myourshaw@ucla.edu)] [<ref.fasta> | <ref.fa> (default human_g1k_v37.fasta)] <vcf file> <list of bam files>";
here=$(dirname "$0");

email="myourshaw@ucla.edu";
if [[ "$1" =~ '\S+@\S+\.\S+'  ]]; then email="$1"; shift; fi

source "$here/b37.variables.sh" $email;
ref="$B37_FASTA";
if [[ ${1##*.} == "fasta" || ${1##*.} == "fa" ]]; then ref="$1"; shift; fi
if (( $# < 2 )); then echo $usage; exit; fi
vcf="$1"; shift;
if [[ ! -e "$vcf" ]]; then echo "vcf file does not exist [$vcf]"; exit; fi

bamlist='';
for bam in $@; do
	if [[ ! -e "$bam" ]]; then echo "bam file does not exist [$bam]"; exit; fi
	bamlist="${bamlist} -I ${bam}";
done

wd="`pwd`";

#files created by this script
output_dir="`dirname $bam`";
mkdir -p "$output_dir";
qout="$output_dir/qout";
mkdir -p "$qout";
qsub_allq="$QSUB_ALLQ -e $qout -o $qout";
qsub_lomem="$QSUB_LOMEM -e $qout -o $qout";
qsub_himem="$QSUB_HIMEM -e $qout -o $qout";
phasedvcf="${vcf%.vcf}.phased.vcf";

#ReadBackedPhasing
name="gatkReadBackedPhasing_`basename $phasedvcf`";
cmd=" \
java -Xmx5g -Djava.io.tmpdir=\"$TMPSCRATCH\" \
 -jar \"$GATK/GenomeAnalysisTK.jar\" \
 -T ReadBackedPhasing \
 -R \"$ref\" \
 $bamlist \
 -B:variant,VCF \"$vcf\" \
 --phaseQualityThresh 10 \
 -o \"$phasedvcf\" \
 ;";
ReadBackedPhasingId=`echo "$cmd" | $qsub_lomem -N $name`;
echo "$ReadBackedPhasingId $name $cmd";
	
cd "$wd";

