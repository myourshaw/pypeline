#!/bin/bash

#novoalign_illumina_q.sh

#input exome illumina qseq files in <RG_ID>/<RG_ID>.pe{1,2}.qseq.txt
#RG_ID=<experiment>.<date>.<machine>.<run>.<lane>.<barcode>
#experiment=<sample>_<library>_<platform>
#/scratch1/tmp/myourshaw/gmd/new/qseq/rg/GMD108_XT_Illumina.2011-05-11.HWI-ST430.243.6.ACTTGAA/GMD108_XT_Illumina.2011-05-11.HWI-ST430.243.6.ACTTGAA.pe1.qseq.txt

#output: aligned, sorted, validated bam files

usage="usage: $0 [<email address> (default: myourshaw@ucla.edu)] <.pe1.qseq.txt file list>";
if (( $# < 1 )); then echo $usage; exit; fi

#get path and parent directory of this script
#http://hintsforums.macworld.com/archive/index.php/t-73839.html
# if the call came via symlink, then use its target instead:
arg=$0; [[ -L $0 ]] && arg=$(stat -L -c'%n' "$0")
script_path=$(2>/dev/null cd "${arg%/*}" >&2; echo "`pwd -P`/${arg##*/}")
script_dir=$(dirname "$script_path")

EMAIL="myourshaw@ucla.edu";
if [[ "$1" =~ '\S+@\S+\.\S+'  ]]; then EMAIL="$1"; shift; fi

source ${script_dir}/b37.variables.sh $EMAIL;

if (( $# < 1 )); then echo $usage; exit; fi
for f in $@; do #/scratch1/tmp/myourshaw/tibhan/SRX*/SRR*/SRR*.fastq
	if [[ ! -e $f ]]; then echo "qseq file does not exist [$f]"; exit; fi
done

wd=`pwd`;

for f in $@; do
	f2=${f%.pe1.qseq.txt}.pe2.qseq.txt;
	if [[ ! -e $f2 ]]; then f2=""; fi
	base=`basename $f`; #GMD108_XT_Illumina.2011-05-11.HWI-ST430.243.6.ACTTGAA.pe1.qseq.txt
	output_dir=`dirname $f`;
	mkdir -p $output_dir;
	qout="$output_dir/qout";
	mkdir -p $qout;
	qsub_allq="$QSUB_ALLQ -e $qout -o $qout";
	qsub_lomem="$QSUB_LOMEM -e $qout -o $qout";
	qsub_himem="$QSUB_HIMEM -e $qout -o $qout";

	prefix=${f%.pe1.qseq.txt};
	
	#files created by this script
	sam=$prefix.novoalign.sam;
	bam=$prefix.novoalign.bam;
	validate=$bam.validate;
	
	#Novoalign
	rgid=${base%.pe1.qseq.txt}; #GMD108_XT_Illumina.2011-05-11.HWI-ST430.243.6.ACTTGAA
	experiment=${rgid%%.*}; #GMD108_XT_Illumina
	sample=${experiment%%_*}; #GMD108
	library=${experiment%_*}; #GMD108_XT
	x=${rgid#*.}; #2011-05-11.HWI-ST430.243.6.ACTTGAA
	dt=${x%%.*}; #2011-05-11
	pu=${x#*.}; #HWI-ST430.243.6.ACTTGAA

	rg='@RG\tID:'$rgid'\tCN:UCLA\tDT:'$dt'\tLB:'$library'\tPL:ILLUMINA\tPU:'$pu'\tSM:'$sample;
	#@RG\tID:GMD108_XT_Illumina.2011-05-11.HWI-ST430.243.6.ACTTGAA\tCN:UCLA\tDT:2011-05-11\tLB:GMD108_XT\tPL:ILLUMINA\tPU:HWI-ST430.243.6.ACTTGAA\tSM:GMD108
	
	adapter1="AGATCGGAAGAGCACACGTCT"; #A + Multiplexing Adapter 1
	adapter2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA"; #reverse complement of Multiplexing Adapter 2 + A
	
	name=novoalign_`basename $sam`;
	cmd="$NOVOALIGN -k -o SAM \"$rg\" -d $NOVOALIGN_HIMEM_INDEX -a $adapter1 $adapter2 -F QSEQ -f $f $f2 > $sam;";
	#NovoalignId=`echo "$cmd" | $qsub_himem -N $name`;
	#echo "$NovoalignId $name $cmd";

	#FixMateInformation
	name=gatkFixMateInformation_`basename $bam`;
	cmd=" \
	java -Xmx5g -Djava.io.tmpdir=$TMPSCRATCH \
	 -jar $PICARD/FixMateInformation.jar \
	 TMP_DIR=$TMPSCRATCH \
	 VERBOSITY=INFO \
	 QUIET=false \
	 VALIDATION_STRINGENCY=STRICT \
	 COMPRESSION_LEVEL=5 \
	 MAX_RECORDS_IN_RAM=1000000 \
	 CREATE_INDEX=true \
	 CREATE_MD5_FILE=true \
	 INPUT=$sam \
	 OUTPUT=$bam \
	 SORT_ORDER=coordinate \
	;";
	#FixMateInformationId=`echo "$cmd" | $qsub_lomem -N $name -hold_jid $NovoalignId`;
	#echo "$FixMateInformationId $name $cmd";
	
	#ValidateSamFile
	name=gatkValidateSamFile_`basename $validate`;
	cmd=" \
	java -Xmx5g  -Djava.io.tmpdir=$TMPSCRATCH \
	 -jar $PICARD/ValidateSamFile.jar \
 	 TMP_DIR=$TMPSCRATCH \
	 VERBOSITY=INFO \
	 QUIET=false \
	 VALIDATION_STRINGENCY=STRICT \
	 COMPRESSION_LEVEL=5 \
	 MAX_RECORDS_IN_RAM=1000000 \
	 INPUT=$bam \
	 OUTPUT=$validate \
	 MODE=VERBOSE \
	 REFERENCE_SEQUENCE=$B37_FASTA \
	 ;";
echo $qsub_lomem;
	ValidateSamFileId=`echo "$cmd" | $qsub_lomem -N $name`; # -hold_jid $FixMateInformationId`;
	echo "$ValidateSamFileId $name $cmd";

done
cd $wd;

