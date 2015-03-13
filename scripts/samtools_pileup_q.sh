#!/bin/bash
#$ -cwd
#$ -V
#$ -M myourshaw@ucla.edu
#$ -m eas
#$ -terse

#runs pileup on all bases, emits extra consensus columns
usage="usage: /home/myourshaw/local/bin/samtools_pileup_q.sh [<reference.fasta | reference.fa>] <list of bam files>";
if (( $# < 1 )); then echo $usage; exit; fi

source ~/local/bin/b37.variables.sh;
ref=$B37_FASTA;
if [[ ${1##*.} == "fasta" || ${1##*.} == "fa" ]]; then ref=$1; shift; fi

wd=`pwd`;

for b in $@; do
 if [[ ! -e $b ]]; then echo "bam file does not exist [$b]"; exit; fi
done

for bam in $@; do
	#output files
	output_dir=`dirname $bam`;
	mkdir -p $output_dir;
	pileup=${bam%.bam}.pileup;
	qout="$output_dir/qout";
	mkdir -p $qout;
	qsub_allq="$QSUB_ALLQ -e $qout -o $qout";
	qsub_lomem="$QSUB_LOMEM -e $qout -o $qout";
	qsub_himem="$QSUB_HIMEM -e $qout -o $qout";

	#pileup
	name=pileup_`basename $pileup`;
	cmd=" \
	 samtools pileup -c -a -B -f $ref $bam > $pileup \
	 ;";
	mpileupId=`echo $cmd | $qsub_lomem -N $name`;
	echo "$mpileupId $name";
done
cd $wd;

