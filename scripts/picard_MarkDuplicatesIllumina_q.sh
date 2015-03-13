#!/bin/bash

#picard_MarkDuplicatesIllumina_q.sh

usage="usage: $0 [<email address> (default: myourshaw@ucla.edu)] [<reference.fasta | reference.fa>] <list of sorted bam files";
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
if (( $# < 1 )); then echo $usage; exit; fi

wd=`pwd`;

for b in $@; do
 if [[ ! -e "$b" ]]; then echo "bam file does not exist [$b]"; exit; fi
done

for bam in $@; do
	#output files
	output_dir="`dirname $bam`";
	mkdir -p "$output_dir";
	markdup="${bam%.bam}.markdup.bam";
	markdupmetrics="${bam%.bam}.markdup.metrics";
	qout="$output_dir/qout";
	mkdir -p "$qout";
	qsub_allq="$QSUB_ALLQ -e $qout -o $qout";
	qsub_lomem="$QSUB_LOMEM -e $qout -o $qout";
	qsub_himem="$QSUB_HIMEM -e $qout -o $qout";

	#MarkDuplicates
	#SOLiD PE READ_NAME_REGEX=([0-9]+)_([0-9]+)_([0-9]+).*
	#SOLiD PE OPTICAL_DUPLICATE_PIXEL_DISTANCE=10
	#Illumina PE READ_NAME_REGEX=.+:.+:[0-9]+:[0-9]+:([0-9]+):([0-9]+):([0-9]+).*
	#Illumina PE OPTICAL_DUPLICATE_PIXEL_DISTANCE=100

	#Illumina QSEQ.novoalign PE READ_NAME_REGEX=
	#HWI-ST430_243_6_7_10665_92077_0 163
	#[^_]+_[^_]+_[0-9]+_([0-9]+)_([0-9]+)_([0-9]+).*
	
	name="MarkDuplicates_`basename $markdup`";
	cmd=" \
	java -Xmx24g -Djava.io.tmpdir=\"$TMPSCRATCH\" \
	-jar \"$PICARD/MarkDuplicates.jar\" \
	 TMP_DIR=$TMPSCRATCH \
	 VERBOSITY=INFO \
	 QUIET=false \
	 VALIDATION_STRINGENCY=STRICT \
	 COMPRESSION_LEVEL=5 \
	 MAX_RECORDS_IN_RAM=1000000 \
	 CREATE_INDEX=true \
	 CREATE_MD5_FILE=true \
	 INPUT=\"$bam\" \
	 OUTPUT=\"$markdup\" \
	 METRICS_FILE=\"$markdupmetrics\" \
	 REMOVE_DUPLICATES=false \
	 MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000 \
	 READ_NAME_REGEX=\"[^_]+_[^_]+_[0-9]+_([0-9]+)_([0-9]+)_([0-9]+).*\" \
	 OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 \
	;";
	MarkDuplicatesId=`echo "$cmd" | $qsub_himem -N $name`;
	echo "$MarkDuplicatesId $name $cmd";
done
cd $wd;

