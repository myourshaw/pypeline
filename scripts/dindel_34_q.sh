#!/bin/bash
#$ -cwd
#$ -V
#$ -M myourshaw@ucla.edu
#$ -m eas
#$ -terse

#dindel_34_q.sh
#called from dindel_q.sh in the same directory
#requires b37.variables.sh in the same directory

#dindel stages 3 and 4
usage="usage: $0 <email address> <reference fasta> <bam file>";
if (( $# < 2 )); then echo $usage; exit; fi
here=$(dirname "$0");
email="$1"; shift;
ref="$1"; shift;
bamFile="$1"; shift;

source "$here/b37.variables.sh" $email;

wd=`pwd`;

#job control lists
holdindels='';

#output files
sample="`basename $bamFile`";
output_dir="`dirname $bamFile`";
mkdir -p "$output_dir";
qout="$output_dir/qout";
mkdir -p "$qout";
qsub_allq="$QSUB_ALLQ -e $qout -o $qout";
qsub_lomem="$QSUB_LOMEM -e $qout -o $qout";
qsub_himem="$QSUB_HIMEM -e $qout -o $qout";

bamBase="`basename $bamFile`";
sampleID="${bamBase%%.*}";
outputFileStage1Prefix="$bamFile.dindel_output"; # > $bamFile.dindel_output.{variants,libraries}.txt
varFile="$outputFileStage1Prefix.variants.txt";
libFile="$outputFileStage1Prefix.libraries.txt";
windowfiledir="`dirname $bamFile`/dindel_realign_windows"; #realignment window files
mkdir -p "$windowfiledir";
windowFilePrefix="$windowfiledir/$bamBase.realign_windows"; # > $windowfiledir/$sample.realign_windows.*.txt
windowFilePattern="$windowFilePrefix'.*.txt'";
outputFilesStage2Prefix="$windowfiledir/$bamBase.dindel_stage2_output_windows"; # > $windowfiledir/$sample.dindel_stage2_output_windows.*.glf.txt
outputFilesStage2Pattern="$outputFilesStage2Prefix'.*.glf.txt'";
stage2OutputFiles="$bamFile.dindel_stage2_outputfiles.txt";
#if [[ -e "$stage2OutputFiles" ]]; then rm "$stage2OutputFiles"; fi
variantCalls="$bamFile.variantCalls.vcf";

#dindel stage 3 indels
for windowFile in "$windowFilePattern"; do
	name="dindelindels_`basename $windowFile`";
	outputFilePrefix="${windowFile/.realign_windows./.dindel_stage2_output_windows.}";
	outputFilePrefix="${outputFilePrefix%.txt}";
	cmd=" \
	\"$DINDEL\" \
	--analysis indels \
	--doDiploid \
	--bamFile \"$bamFile\" \
	--ref \"$ref\" \
	--varFile \"$windowFile\" \
	--libFile \"$libFile\" \
	--outputFile \"$outputFilePrefix\" \
	;";
	indelsId=`echo "$cmd" | $qsub_lomem -N $name`;
	echo "$indelsId $name $cmd";
	if [[ $holdindels == '' ]]; then holdindels=${indelsId}; else holdindels="${holdindels},${indelsId}"; fi
done

#dindel stage 4 mergeOutputDiploid
name="dindelmergeOutputDiploid_`basename $variantCalls`";
cmd=" \
for f in "$outputFilesStage2Pattern"; do \
echo \$f >> \"$stage2OutputFiles\"; \
done; \
python \"$DINDEL_PYTHON/mergeOutputDiploid.py\" \
--inputFiles \"$stage2OutputFiles\" \
--outputFile \"$variantCalls\" \
--ref \"$ref\" \
--sampleID \"$sampleID\" \
;";
mergeOutputDiploidId=`echo "$cmd" | $qsub_lomem -N $name -hold_jid $holdindels`;
echo "$mergeOutputDiploidId $name $cmd";

cd $wd;

