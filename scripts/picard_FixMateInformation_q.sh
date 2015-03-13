#!/bin/bash
#$ -cwd
#$ -V
#$ -M myourshaw@ucla.edu
#$ -m eas
#$ -terse

usage="usage: /home/myourshaw/local/bin/picard_FixMateInformation_q.sh [<reference.fasta | reference.fa>] <list of sam or bam files>";
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
	validate=$bam.validate;
	
	qout="$output_dir/qout";
	mkdir -p $qout;
	qsub_allq="$QSUB_ALLQ -e $qout -o $qout";
	qsub_lomem="$QSUB_LOMEM -e $qout -o $qout";
	qsub_himem="$QSUB_HIMEM -e $qout -o $qout";

	#FixMateInformation
	name=FixMateInformation_`basename $bam`;
	cmd=" \
	java -Xmx5g -Djava.io.tmpdir=$TMPSCRATCH \
	 -jar $PICARD/FixMateInformation.jar \
	 TMP_DIR=$TMPLOCAL \
	 VERBOSITY=INFO \
	 QUIET=false \
	 VALIDATION_STRINGENCY=STRICT \
	 COMPRESSION_LEVEL=5 \
	 MAX_RECORDS_IN_RAM=1000000 \
	 CREATE_INDEX=true \
	 CREATE_MD5_FILE=true \
	 INPUT=$bam \
	 SORT_ORDER=coordinate \
	;";
	FixMateInformationId=`echo "$cmd" | $qsub_lomem -N $name`;
	echo "$FixMateInformationId $name";
		
	#ValidateSamFile
	name=ValidateSamFile_`basename $validate`;
	cmd=" \
	java -Xmx5g -Djava.io.tmpdir=$TMPSCRATCH \
	 -jar $PICARD/ValidateSamFile.jar \
	 TMP_DIR=$TMPLOCAL \
	 VERBOSITY=INFO \
	 QUIET=false \
	 VALIDATION_STRINGENCY=STRICT \
	 COMPRESSION_LEVEL=5 \
	 MAX_RECORDS_IN_RAM=1000000 \
	 INPUT=$bam \
	 OUTPUT=$validate \
	 MODE=VERBOSE \
	 REFERENCE_SEQUENCE=$ref \
	 ;";
	ValidateSamFileId=`echo "$cmd" | $qsub_lomem -N $name -hold_jid $FixMateInformationId`;
	echo "$ValidateSamFileId $name";
done

cd $wd;

