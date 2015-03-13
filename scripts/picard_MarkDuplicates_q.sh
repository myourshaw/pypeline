#!/bin/bash
usage="usage: picard_MarkDuplicates_q.sh <list of sorted bam files";
if (( $# < 1 )); then echo $usage; exit; fi
wd=`pwd`;
for i in $@;
do
	picard=/share/apps/picard-tools-current;
	resources=/scratch0/tmp/myourshaw/genomes/resources;
	ref=$resources/Homo_sapiens_assembly18.fasta;
	tmp=/state/partition1/tmp/`whoami`;
	if [[ -d `dirname $i` ]]; then cd `dirname $i`; fi
	base=`basename $i`;
	prefix=${i%.bam};
	o=$prefix.rmdup.bam;
	m=$prefix.rmdup.metrics;
	cmd="java -Xmx6g -jar $picard/MarkDuplicates.jar I=$i O=$o M=$m VALIDATION_STRINGENCY=SILENT TMP_DIR=$tmp MAX_RECORDS_IN_RAM=3000000 REMOVE_DUPLICATES=true ASSUME_SORTED=true";
	echo $cmd | qsub -N MarkDuplicates_$base -q all.q@compute-[23]* -pe serial 8 -l vf=3G -cwd -v PATH -M myourshaw@ucla.edu -m eas;
done
cd $wd;

