#!/bin/bash

usage="usage: $0 <list of vcf files>";

if (( $# < 1 )); then echo $usage; exit; fi

for vcf in $@; do
	if [[ ! -e $vcf ]]; then echo "vcf file does not exist [$vcf]"; exit; fi
done

VARIANT_EFFECT_PREDICTOR="/Users/myourshaw/apps/variant_effect_predictor/variant_effect_predictor.pl";

wd=`pwd`;

for vcf in $@; do
	#files created by this script
	output_dir=`dirname $vcf`;
	mkdir -p $output_dir;
	
	txt=${vcf}.ensembl_variant_effect.txt;
	
	#VariantEffectPredictor
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
	"$cmd";
done

cd $wd;

