#!/bin/bash
#$ -cwd
#$ -V
#$ -M myourshaw@ucla.edu
#$ -m eas
#$ -terse

#gatk_VariantEval_q.sh

usage="usage: $0 [<email address> (default: myourshaw@ucla.edu)] [<ref.fasta> | <ref.fa> (default human_g1k_v37.fasta)] <list of vcf files>";
here=$(dirname "$0");

email="myourshaw@ucla.edu";
if [[ "$1" =~ '\S+@\S+\.\S+'  ]]; then email="$1"; shift; fi

source "$here/b37.variables.sh" $email;
ref="$B37_FASTA";
if [[ ${1##*.} == "fasta" || ${1##*.} == "fa" ]]; then ref="$1"; shift; fi
if (( $# < 1 )); then echo $usage; exit; fi
for vcf in $@; do
 if [[ ! -e "$vcf" ]]; then echo "vcf file does not exist [$vcf]"; exit; fi
done

wd="`pwd`";

#resources
dbsnp="$DBSNP_132_B37_VCF";
refgene="$REFGENE_B37_ROD";
hapmap="$HAPMAP_R27_B37_VCF";
kg="$ONEKG_AUTOSOME_B37_VCF";

for vcf in $@; do
	#files created by this script
	output_dir="`dirname $vcf`";
	mkdir -p "$output_dir";
	qout="$output_dir/qout";
	mkdir -p "$qout";
	qsub_allq="$QSUB_ALLQ -e $qout -o $qout";
	qsub_lomem="$QSUB_LOMEM -e $qout -o $qout";
	qsub_himem="$QSUB_HIMEM -e $qout -o $qout";
	variantEval="${vcf}.eval";
	
	#VariantEval
	name="gatkVariantEval_`basename $variantEval`";
	cmd=" \
	java -Xmx5g -Djava.io.tmpdir=\"$TMPSCRATCH\" \
	 -jar \"$GATK/GenomeAnalysisTK.jar\" \
	 -T VariantEval \
	 -R \"$ref\" \
	 -l INFO \
	 -B:eval,VCF \"$vcf\" \
	 -B:dbsnp,VCF \"$dbsnp\" \
	 -o \"$variantEval\" \
	 -ST Sample \
	 -L /scratch0/tmp/bakeoff/capture.interval/SureSelect_Sanger_hg19_100bp_merged_intersected_RefSeqG37_extended2bp.interval_list \
	 ;";
	VariantEvalId=`echo "$cmd" | $qsub_lomem -N $name`;
	echo "$VariantEvalId $name $cmd";
done

cd "$wd";

