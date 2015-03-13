#!/bin/bash
#$ -cwd
#$ -V
#$ -M myourshaw@ucla.edu
#$ -m eas
#$ -terse

#picard_ValidateSamFile_q.sh

usage="usage: $0 [<email address> (default: myourshaw@ucla.edu)] [<reference.fasta | reference.fa>] <list of sam or bam files>";
here=$(dirname "$0");

email="myourshaw@ucla.edu";
if [[ "$1" =~ '\S+@\S+\.\S+'  ]]; then email="$1"; shift; fi

source "$here/b37.variables.sh" $email;
ref="$B37_FASTA";
if [[ ${1##*.} == "fasta" || ${1##*.} == "fa" ]]; then ref="$1"; shift; fi
if (( $# < 1 )); then echo $usage; exit; fi

wd="`pwd`";

for b in $@; do
 if [[ ! -e "$b" ]]; then echo "bam file does not exist [$b]"; exit; fi
done

for bam in $@; do
	#output files
	output_dir="`dirname $bam`";
	mkdir -p "$output_dir";
	validate="$bam.validate";
	qout="$output_dir/qout";
	mkdir -p "$qout";
	qsub_allq="$QSUB_ALLQ -e $qout -o $qout";
	qsub_lomem="$QSUB_LOMEM -e $qout -o $qout";
	qsub_himem="$QSUB_HIMEM -e $qout -o $qout";

	#ValidateSamFile
	name="ValidateSamFile_`basename $validate`";
	cmd=" \
	java -Xmx5g -Djava.io.tmpdir=\"$TMPSCRATCH\" \
	-jar \"$PICARD/ValidateSamFile.jar\" \
	 TMP_DIR=\"$TMPSCRATCH\" \
	 VERBOSITY=INFO \
	 QUIET=false \
	 VALIDATION_STRINGENCY=STRICT \
	 COMPRESSION_LEVEL=5 \
	 MAX_RECORDS_IN_RAM=1000000 \
	 INPUT=\"$bam\" \
	 OUTPUT=\"$validate\" \
	 MODE=VERBOSE \
	 REFERENCE_SEQUENCE=\"$ref\" \
	 IGNORE=MISSING_TAG_NM \
	 ;";
	ValidateSamFileId=`echo $cmd | $qsub_lomem -N $name`;
	echo "$ValidateSamFileId $name $cmd";
done

cd "$wd";

