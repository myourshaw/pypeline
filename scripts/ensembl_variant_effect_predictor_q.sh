#!/bin/bash
#$ -cwd
#$ -V
#$ -M myourshaw@ucla.edu
#$ -m eas
#$ -terse

usage="usage: $0 [<email address>] <list of vcf files>";

email="myourshaw@ucla.edu";
if [[ "$1" =~ '\S+@\S+\.\S+'  ]]; then email="$1"; shift; fi

if (( $# < 1 )); then echo $usage; exit; fi

source `dirname $0`/b37.variables.sh $email;

for vcf in $@; do
	if [[ ! -e $vcf ]]; then echo "vcf file does not exist [$vcf]"; exit; fi
done

wd=`pwd`;

for vcf in $@; do
	#files created by this script
	output_dir=`dirname $vcf`;
	mkdir -p $output_dir;
	qout="$output_dir/qout";
	mkdir -p $qout;
	qsub_allq="$QSUB_ALLQ -e $qout -o $qout";
	qsub_lomem="$QSUB_LOMEM -e $qout -o $qout";
	qsub_himem="$QSUB_HIMEM -e $qout -o $qout";
	
	txt=${vcf}.ensembl_variant_effect.txt;
	
	#VariantEffectPredictor
	name=ensemblVariantEffectPredictor_`basename $txt`;
	cmd=" \
	perl $VARIANT_EFFECT_PREDICTOR \
	--input_file $vcf \
	--format vcf \
	--output_file $txt \
	--sift=b \
	--polyphen=b \
	--condel=b \
	--hgnc \
	--hgvs \
	--species human \
	--buffer_size 10000 \
	--host cortex.local \
	--user ensembl \
	--port 3306 \
	--password ensembl \
	;";
	VariantEffectPredictorId=`echo "$cmd" | $qsub_allq -N $name`;
	echo "$VariantEffectPredictorId $name";
done

cd $wd;

