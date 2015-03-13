#!/bin/bash
#$ -cwd
#$ -V
#$ -M myourshaw@ucla.edu
#$ -m eas
#$ -terse

usage="usage: /home/myourshaw/local/bin/samtools_idxstats_q.sh <list of bam files>";

wd=`pwd`;

#applications and resources
picard="/share/apps/picard-tools-current";
samtools="/share/apps/samtools-current/samtools";
resources="/scratch1/tmp/myourshaw/resources";
ref="$resources/human_g1k_v37.fasta";
tmplocal="/state/partition1/tmp/`whoami`";
tmpscratch="/scratch1/tmp/`whoami`/tmp";
mkdir -p $tmpscratch;

#qsub options
email="myourshaw@ucla.edu";

for bam in $@; do
	#output files
	output_dir=`dirname $bam`;
	mkdir -p $output_dir;
	idxstats=$bam.idxstats;
	
	qout="$output_dir/qout";
	mkdir -p $qout;
	qsub_allq='qsub -q all.q -cwd -V -e '$qout' -o '$qout' -M '$email' -m eas -terse';
	qsub_lomem='qsub -q all.q@compute-4* -cwd -V -e '$qout' -o '$qout' -M '$email' -m eas -terse -hard -pe serial 8 -l mem_free=6G';
	qsub_himem='qsub -q all.q@compute-[23]* -cwd -V -e '$qout' -o '$qout' -M '$email' -m eas -terse -hard -pe serial 8 -l mem_free=24G';

	#idxstats
	name=idxstats_`basename $idxstats`;
	cmd=" \
	$samtools idxstats $bam > $idxstats \
	;";
	idxstatsId=`echo "$cmd" | $qsub_allq -N $name`;

done

cd $wd;

