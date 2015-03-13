#!/bin/bash
#$ -cwd
#$ -V
#$ -M sptaylor@ucla.edu
#$ -m eas
#$ -terse

usage="usage: /home/tpaige/local/bin/variant_effect_predictor.tpaige_q.sh  [<email address>] <list of vcf files>";

email="sptaylor@ucla.edu";
if [[ "$1" =~ '\S+@\S+\.\S+'  ]]; then email="$1"; shift; fi

if (( $# < 1 )); then echo $usage; exit; fi

here=$(dirname "$0");
source "$here/b37.variables.sh" $email;

for vcf in $@; do
	if [[ ! -e $vcf ]]; then echo "vcf file does not exist [$vcf]"; exit; fi
done

wd=`pwd`;

for vcf in $@; do
	#files created by this script
	output_dir=`dirname $vcf`;
	mkdir -p $output_dir;
	qout="$output_dir/qout/vep";
	tempdir="$output_dir/tmp";
	mkdir -p $tempdir;
	mkdir -p $qout;
	qsub_allq="$QSUB_ALLQ -o $qout -j y";
	qsub_lomem="qsub -q all.q@compute-4*,all.q@compute-[23]* -cwd -V -M sptaylor@ucla.edu -m eas -terse -hard -pe serial 2 -l vf=800m -o $qout -j y";
	
	bgzip -c $vcf > $tempdir/$vcf.bgz;
	tabix -f -p vcf $tempdir/$vcf.bgz;

	job_ids=();
#	for i in $(cut -f1 $B37_FASTA_FAI ); do 
#		tabix -p vcf -f $tempdir/$vcf.bgz $i > $tempdir/chr$i.$vcf;
#		echo "creating $vcf temp file for chr$i" ; 
	for i in $(vcfutils.pl splitchr -l 12500000 $B37_FASTA_FAI ); do 
		tabix -p vcf -f $tempdir/$vcf.bgz $i > $tempdir/chr$i.$vcf;
		echo "creating $vcf temp file for chr$i" ; 

		txt=${vcf}.chr${i/:/_}.VEP.txt;
	
	#VariantEffectPredictor
		name=vep_`basename $txt`;
		cmd=" \
		perl /home/tpaige/lab/git-tpaige/variant_effect_predictor.gt.cons.phenotype.faster.pl -i $tempdir/chr$i.$vcf \
		-o $tempdir/$txt \
		--use_ncbi \
		--check_existing=1 \
		--failed=0 \
		--use_uniprot \
		-format vcf \
		;";
		VariantEffectPredictorId=`echo "$cmd" | $qsub_lomem -N $name`;
		echo "$VariantEffectPredictorId $name";
		job_ids+=($VariantEffectPredictorId);

	done
	
	jid_list=$(printf ",%s" "${job_ids[@]}");
	jid_list=${jid_list:1};
	
	#Cat temp files
	name=cat_tempfiles_`basename $vcf`;
	header=/scratch1/tmp/tpaige/resources/vep.header.txt;
	cmd=" \
	cat $tempdir/${vcf}.chr*.VEP.txt | grep -v '#' | sort-alt -N +1.0 | cat $header - > $output_dir/${vcf}.VEP.txt
	;";
       CatTempFilesId=`echo "$cmd" | $qsub_lomem -N $name -hold_jid $jid_list`;
       echo "$CatTempFilesId $name";

       #Remove tempfiles
       name=cleanup_`basename $vcf`;
       cmd="\
       rm $tempdir/${vcf}.chr*.VEP.txt $tempdir/chr*.$vcf
       ;";
       CleanupFilesID=`echo "$cmd" | $qsub_lomem -N $name -hold_jid $CatTempFilesId`;
       echo "$CleanupFilesID $name";



done

cd $wd
