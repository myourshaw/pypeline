#!/bin/bash
#$ -cwd
#$ -V
#$ -M myourshaw@ucla.edu
#$ -m eas
#$ -terse

#dindel_q.sh
#requires dindel_34_q.sh in the same directory
#requires b37.variables.sh in the same directory

usage="usage: $0 [email address] [<ref.fasta> | <ref.fa>] <list of sorted bam files>";
here=$(dirname "$0");

email="myourshaw@ucla.edu";
if [[ "$1" =~ '\S+@\S+\.\S+'  ]]; then email="$1"; shift; fi
source "$here/b37.variables.sh" $email;
ref="$B37_FASTA";
if [[ ${1##*.} == "fasta" || ${1##*.} == "fa" ]]; then ref="$1"; shift; fi
if (( $# < 1 )); then echo $usage; exit; fi

wd=`pwd`;

for bam in $@; do
	#output files
	sample="`basename $bam`";
	output_dir="`dirname $bam`";
	mkdir -p "$output_dir";
	qout="$output_dir/qout";
	mkdir -p "$qout";
	qsub_allq="$QSUB_ALLQ -e $qout -o $qout";
	qsub_lomem="$QSUB_LOMEM -e $qout -o $qout";
	qsub_himem="$QSUB_HIMEM -e $qout -o $qout";
	
	bamBase="`basename $bam`";
	sampleID="${bamBase%%.*}";
	outputFileStage1Prefix="$bam.dindel_output"; # > $bam.dindel_output.{variants,libraries}.txt
	varFile="$outputFileStage1Prefix.variants.txt";
	libFile="$outputFileStage1Prefix.libraries.txt";
	windowfiledir="`dirname $bam`/dindel_realign_windows"; #realignment window files
	mkdir -p "$windowfiledir";
	windowFilePrefix="$windowfiledir/$bamBase.realign_windows"; # > $windowfiledir/$sample.realign_windows.*.txt
	windowFilePattern="$windowFilePrefix'.*.txt'";
	outputFilesStage2Prefix="$windowfiledir/$bamBase.dindel_stage2_output_windows"; # > $windowfiledir/$sample.dindel_stage2_output_windows.*.glf.txt
	outputFilesStage2Pattern="$outputFilesStage2Prefix'.*.glf.txt'";
	stage2OutputFiles="$bam.dindel_stage2_outputfiles.txt";
	#if [[ -e "$stage2OutputFiles" ]]; then rm "$stage2OutputFiles"; fi
	variantCalls="$bam.variantCalls.vcf";

	#dindel stage 1 getCIGARindels
	name="dindelgetCIGARindels_`basename $outputFileStage1Prefix`";
	cmd=" \
	\"$DINDEL\" \
	--analysis getCIGARindels \
	--bamFile \"$bam\" \
	--outputFile \"$outputFileStage1Prefix\" \
	--ref \"$ref\" \
	;";
	getCIGARindelsId=`echo "$cmd" | $qsub_lomem -N $name`;
	echo "$getCIGARindelsId $name $cmd";
	
	#dindel stage 2 makeWindows
	name="dindelmakeWindows_`basename $windowFilePrefix`";
	cmd=" \
	python \"$DINDEL_PYTHON/makeWindows.py\" \
	--inputVarFile \"$varFile\" \
	--windowFilePrefix \"$windowFilePrefix\" \
	--numWindowsPerFile 1000 \
	;";
	makeWindowsId=`echo "$cmd" | $qsub_lomem -N $name -hold_jid $getCIGARindelsId`;
	echo "$makeWindowsId $name $cmd";
	
	#dindel stages 3 and 4
	name="dindel34_`basename $windowFilePrefix`";
	cmd=" \
	\"$here/dindel_34_q.sh\" \
	$email \
	\"$ref\" \
	\"$bam\" \
	;";
	dindel34Id=`echo "$cmd" | $qsub_lomem -N $name -hold_jid $makeWindowsId`;
	echo "$dindel34Id $name $cmd";
done

cd $wd;

