#!/bin/bash
#$ -cwd
#$ -V
#$ -M myourshaw@ucla.edu
#$ -m eas
#$ -terse

#gatk_UnifiedGenotyper_SingleSample_q.sh

usage="usage: $0 [<email address> (default: myourshaw@ucla.edu)] [<ref.fasta> | <ref.fa> (default human_g1k_v37.fasta)] <list of sorted, recalibrated, realigned bam files>";
here=$(dirname "$0");

email="myourshaw@ucla.edu";
if [[ "$1" =~ '\S+@\S+\.\S+'  ]]; then email="$1"; shift; fi

source "$here/b37.variables.sh" $email;
ref="$B37_FASTA";
if [[ ${1##*.} == "fasta" || ${1##*.} == "fa" ]]; then ref="$1"; shift; fi
if (( $# < 1 )); then echo $usage; exit; fi

for bam in $@; do
	if [[ ! -e "$bam" ]]; then echo "bam file does not exist [$bam]"; exit; fi
done

wd="`pwd`";

#resources
dbsnp="$DBSNP_132_B37_VCF";
refgene="$REFGENE_B37_ROD";
hapmap="$HAPMAP_R27_B37_VCF";
kg="$ONEKG_AUTOSOME_B37_VCF";

for bam in $@; do
	#files created by this script
	output_dir="`dirname $bam`";
	mkdir -p "$output_dir";
	qout="$output_dir/qout";
	mkdir -p "$qout";
	qsub_allq="$QSUB_ALLQ -e $qout -o $qout";
	qsub_lomem="$QSUB_LOMEM -e $qout -o $qout";
	qsub_himem="$QSUB_HIMEM -e $qout -o $qout";
	snpvcf="${prefix}.gatk_ug_snps.raw.vcf";
	indelvcf="${prefix}.gatk_ug_indels.raw.vcf";
	indelPrefix="${bam%.bam}.gatk_ig_indels";
	indelGenotyperBed="$indelPrefix.bed";
	indelGenotyperVerbose="$indelPrefix.txt";
	indelGenotyperVcf="$indelPrefix.vcf";
	indelMask="$indelPrefix.mask.bed";
	
	#UnifiedGenotyper-SNPs
	name="gatkUnifiedGenotyperSnp_`basename $snpvcf`";
	cmd=" \
	java -Xmx5g -Djava.io.tmpdir=\"$TMPSCRATCH\" \
	 -jar \"$GATK/GenomeAnalysisTK.jar\" \
	 -R \"$ref\" \
	 -T UnifiedGenotyper \
	 -I \"$bam\" \
	 -B:dbsnp,VCF \"$dbsnp\" \
	 -o \"$snpvcf\" \
	 -stand_call_conf 50.0 \
	 -stand_emit_conf 10.0 \
	 -dcov 1000 \
	 -G Standard \
	 -G Experimental \
	 -G WorkInProgress \
	 -A SBByDepth \
	 -nt 8 \
	 -baq CALCULATE_AS_NECESSARY \
	 ;";
	UnifiedGenotyperSnpId=`echo "$cmd" | $qsub_lomem -N $name`;
	echo "$UnifiedGenotyperSnpId $name $cmd";
	
	#UnifiedGenotyper-indels
	name="gatkUnifiedGenotyperIndel_`basename $indelvcf`";
	cmd=" \
	java -Xmx5g -Djava.io.tmpdir=\"$TMPSCRATCH\" \
	 -jar \"$GATK/GenomeAnalysisTK.jar\" \
	 -R \"$ref\" \
	 -T UnifiedGenotyper \
	 -I \"$bam\" \
	 -B:dbsnp,VCF \"$dbsnp\" \
	 -o \"$indelvcf\" \
	 -stand_call_conf 50.0 \
	 -stand_emit_conf 10.0 \
	 -dcov 1000 \
	 -G Standard \
	 -G Experimental \
	 -G WorkInProgress \
	 -A SBByDepth \
	 -nt 8 \
	 -baq CALCULATE_AS_NECESSARY \
	 -glm DINDEL \
	 ;";
	UnifiedGenotyperIndelId=`echo "$cmd" | $qsub_lomem -N $name`;
	echo "$UnifiedGenotyperIndelId $name $cmd";

	#IndelGenotyperV2
	name="gatkIndelGenotyperV2_`basename $indelGenotyperVcf`";
	cmd=" \
	java -Xmx5g -Djava.io.tmpdir=\"$TMPSCRATCH\" \
	 -jar \"$GATK/GenomeAnalysisTK.jar\" \
	 -T IndelGenotyperV2 \
	 -l INFO \
	 -R \"$ref\" \
	 -I \"$bam\" \
	 -bed \"$indelGenotyperBed\" \
	 -verbose \"$indelGenotyperVerbose\" \
	 -o \"$indelGenotyperVcf\" \
	 --refseq \"$refgene\" \
	 -ws 600 \
	 ;";
	IndelGenotyperV2Id=`echo "$cmd" | $qsub_lomem -N $name`;
	echo "$IndelGenotyperV2Id $name $cmd";

	#MakeIndelMask
	name="gatkMakeIndelMask_`basename $indelMask`";
	cmd="python \"$GATK_PYTHON/makeIndelMask.py\" \"$indelGenotyperBed\" 10 \"$indelMask\";
	MakeIndelMaskId=`echo "$cmd" | $qsub_lomem -N $name -hold_jid $IndelGenotyperV2Id`;
	echo "$MakeIndelMaskId $name $cmd";
done

cd "$wd";

