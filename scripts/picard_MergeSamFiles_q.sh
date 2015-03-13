#!/bin/bash

#picard_MergeSamFiles_q.sh

usage="usage: $0 [<reference.fasta | reference.fa>] <merged output bam file> <list of sam or bam files to be merged>";

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

if (( $# < 3 )); then echo $usage; exit; fi

wd=`pwd`;

bam=$1; shift;if (( $# < 2 )); then echo $usage; exit; fi

bamlist='';
for b in $@; do
 if [[ ! -e $b ]]; then echo "bam file does not exist [$b]"; exit; fi
 bamlist="${bamlist} INPUT=${b}";
done

#output files
output_dir=`dirname $bam`;
mkdir -p $output_dir;
validate=$bam.validate;
qout="$output_dir/qout";
mkdir -p $qout;
qsub_allq="$QSUB_ALLQ -e $qout -o $qout";
qsub_lomem="$QSUB_LOMEM -e $qout -o $qout";
qsub_himem="$QSUB_HIMEM -e $qout -o $qout";

#MergeSamFiles
name=MergeSamFiles_`basename $bam`;
cmd=" \
java -Xmx5g -Djava.io.tmpdir=\"$TMPSCRATCH\" \
 -jar \"$PICARD/MergeSamFiles.jar\" \
 TMP_DIR=\"$TMPSCRATCH\" \
 VERBOSITY=INFO \
 QUIET=false \
 VALIDATION_STRINGENCY=STRICT \
 COMPRESSION_LEVEL=5 \
 MAX_RECORDS_IN_RAM=1000000 \
 CREATE_INDEX=true \
 CREATE_MD5_FILE=true \
 $bamlist \
 OUTPUT=\"$bam\" \
 SORT_ORDER=coordinate \
 ASSUME_SORTED=false \
 MERGE_SEQUENCE_DICTIONARIES=false \
 USE_THREADING=true \
;";
MergeSamFilesID=`echo "$cmd" | $qsub_lomem -N $name`;
echo "$MergeSamFilesID $name $cmd";

#FixMateInformation
name=FixMateInformation_`basename $bam`;
cmd=" \
java -Xmx5g -Djava.io.tmpdir=\"$TMPSCRATCH\" \
 -jar \"$PICARD/FixMateInformation.jar\" \
 TMP_DIR=\"$TMPSCRATCH\" \
 VERBOSITY=INFO \
 QUIET=false \
 VALIDATION_STRINGENCY=STRICT \
 COMPRESSION_LEVEL=5 \
 MAX_RECORDS_IN_RAM=1000000 \
 CREATE_INDEX=true \
 CREATE_MD5_FILE=true \
 INPUT=\"$bam\" \
 SORT_ORDER=coordinate \
;";
FixMateInformationId=`echo "$cmd" | $qsub_lomem -N $name -hold_jid $MergeSamFilesID`;
echo "$FixMateInformationId $name $cmd";
	
#ValidateSamFile
name=ValidateSamFile_`basename $validate`;
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
 ;";
ValidateSamFileId=`echo "$cmd" | $qsub_lomem -N $name -hold_jid $FixMateInformationId`;
echo "$ValidateSamFileId $name $cmd";

cd $wd;

