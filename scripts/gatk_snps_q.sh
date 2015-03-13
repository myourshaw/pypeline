#!/bin/bash
#$ -cwd
#$ -V
#$ -M myourshaw@ucla.edu
#$ -m eas
#$ -terse

#gatk_snps_q.sh
#run after UnifiedGenotyper

usage="usage: $0 [<email address> (default: myourshaw@ucla.edu)] [<ref.fasta> | <ref.fa> (default human_g1k_v37.fasta)] <sample(s).gatk_ug_snps.raw.vcf> <sample(s).gatk_ug_indels.vcf> [<bam file(s)>]";
here=$(dirname "$0");

email="myourshaw@ucla.edu";
if [[ "$1" =~ '\S+@\S+\.\S+'  ]]; then email="$1"; shift; fi

source "$here/b37.variables.sh" $email;
ref="$B37_FASTA";
if [[ ${1##*.} == "fasta" || ${1##*.} == "fa" ]]; then ref="$1"; shift; fi
if (( $# < 2 )); then echo $usage; exit; fi
vcf="$1"; shift;
mask="$1"; shift;
bamlist='';
for bam in $@; do
 if [[ ! -e "$bam" ]]; then echo "bam file does not exist [$bam]"; exit; fi
 bamlist="${bamlist} -I ${bam}";
done

wd="`pwd`";

#resources
dbsnp="$DBSNP_132_B37_VCF";
refgene="$REFGENE_B37_ROD";
hapmap="$HAPMAP_R27_B37_VCF";
kg="$ONEKG_AUTOSOME_B37_VCF";

#files created by this script
output_dir="`dirname $vcf`";
mkdir -p "$output_dir";
qout="$output_dir/qout";
mkdir -p "$qout";
qsub_allq="$QSUB_ALLQ -e $qout -o $qout";
qsub_lomem="$QSUB_LOMEM -e $qout -o $qout";
qsub_himem="$QSUB_HIMEM -e $qout -o $qout";
prefilteredVcf="${vcf%.vcf}.prefiltered.vcf";
prefilteredVcfTable="$prefilteredVcf.table";
variantRecalibratorReportPrefix="${vcf%.vcf}.`basename $hapmap`-`basename $dbsnp`-`basename $kg`";
clusterFile="$variantRecalibratorReportPrefix.cluster";
variantRecalibratedVcf="${vcf%.vcf}.recalibrated.vcf";
tranchesFile="$variantRecalibratedVcf.tranches";
variantRecalibratedFilteredVcf="${variantRecalibratedVcf%.vcf}.filtered.vcf";
variantEval="${variantRecalibratedFilteredVcf}.eval.csv";
variantRecalibratedFilteredAnnotatedVcf="${variantRecalibratedFilteredVcf%.vcf}.annotated.vcf";

#VariantFiltrationSnp
name="gatkVariantFiltrationSnp_`basename $prefilteredVcf`";
cmd=" \
java -Xmx5g -Djava.io.tmpdir=\"$TMPSCRATCH\" \
 -jar \"$GATK/GenomeAnalysisTK.jar\" \
 -T VariantFiltration \
 -l INFO \
 -R \"$ref\" \
 -o \"$prefilteredVcf\" \
 -B:variant,VCF \"$vcf\" \
 -B:mask,VCF \"$mask\" \
 --maskName InDel \
  --clusterWindowSize 10 \
 --filterExpression 'MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)' \
  --filterName HARD_TO_VALIDATE \
 --filterExpression 'QUAL < 30.0 || QD < 5.0 || HRun > 5 || SB > -0.10' \
  --filterName GATKStandard \
 ;";
#note that the * requires quoting $cmd
VariantFiltrationId=`echo "$cmd" | $qsub_lomem -N $name`; # -hold_jid $UnifiedGenotyperId`;
echo "$VariantFiltrationId $name $cmd";

#GenerateVariantClusters
name="gatkGenerateVariantClusters_`basename $clusterFile`";
cmd=" \
java -Xmx5g -Djava.io.tmpdir=\"$TMPSCRATCH\" \
 -jar \"$GATK/GenomeAnalysisTK.jar\" \
 -R \"$ref\" \
 -B:input,VCF \"$prefilteredVcf\" \
 -B:hapmap,VCF \"$hapmap\" \
 -B:1kg,VCF \"$kg\" \
 -B:dbsnp,VCF \"$dbsnp\" \
 -l INFO \
 -an QD -an SB -an HaplotypeScore -an HRun \
 -clusterFile \"$clusterFile\" \
 -T GenerateVariantClusters \
 --maxGaussians 8 \
;";
GenerateVariantClustersId=`echo "$cmd" | $qsub_lomem -N $name -hold_jid $VariantFiltrationId`;
echo "$GenerateVariantClustersId $name $cmd";

#Vcf2table
name="gatkVcf2table_`basename $prefilteredVcfTable`";
cmd=" \
python \"$GATK_PYTHON/vcf2table.py\" \
 -f CHROM,POS,ID,AC,AF,AN,DB,DP,HRun,MQ,MQ0,HaplotypeScore,QD,SB \
 \"$prefilteredVcf\" \
 > \"$prefilteredVcfTable\" \
 ;";
Vcf2tableId=`echo "$cmd" | $qsub_lomem -N $name -hold_jid $VariantFiltrationId`;
echo "$Vcf2tableId $name $cmd";

#VariantRecalibratorReport
name="gatkVariantRecalibratorReport_`basename $variantRecalibratorReportPrefix`";
cmd=" \
\"$RSCRIPT\" \"$GATK_R/VariantRecalibratorReport/VariantRecalibratorReport.R\" \
 \"$variantRecalibratorReportPrefix\" \
 \"$clusterFile\" \
 \"$prefilteredVcfTable\" \
;";
VariantRecalibratorReportId=`echo "$cmd" | $qsub_lomem -N $name -hold_jid $GenerateVariantClustersId,$Vcf2tableId`;
echo "$VariantRecalibratorReportId $name $cmd";

#VariantRecalibrator
name="gatkVariantRecalibrator_`basename $variantRecalibratedVcf`";
cmd=" \
java -Xmx5g -Djava.io.tmpdir=\"$TMPSCRATCH\" \
 -jar \"$GATK/GenomeAnalysisTK.jar\" \
 -R \"$ref\" \
 -B:input,VCF \"$prefilteredVcf\" \
 -B:hapmap,VCF \"$hapmap\" \
 -B:truth,VCF \"$hapmap\" \
 -B:1kg,VCF \"$kg\" \
 -B:dbsnp,VCF \"$dbsnp\" \
 -l INFO \
 --ignore_filter HARD_TO_VALIDATE \
 --ignore_filter LowQual \
 -clusterFile \"$clusterFile\" \
 -o \"$variantRecalibratedVcf\" \
 -tranchesFile \"$tranchesFile\" \
 --target_titv 3.0 \
 -sm TRUTH_SENSITIVITY \
 -tranche 50 -tranche 10  -tranche 1 -tranche 0.1 \
 -resources \"$GATK_R\" \
 -Rscript \"$RSCRIPT\" \
 -T VariantRecalibrator \
 -baq CALCULATE_AS_NECESSARY \
;";
VariantRecalibratorId=`echo "$cmd" | $qsub_lomem -N $name -hold_jid $GenerateVariantClustersId`;
echo "$VariantRecalibratorId $name $cmd";

#ApplyVariantCuts
name="gatkApplyVariantCuts_`basename $variantRecalibratedFilteredVcf`";
cmd=" \
java -Xmx5g -Djava.io.tmpdir=\"$TMPSCRATCH\" \
 -jar \"$GATK/GenomeAnalysisTK.jar\" \
 -R \"$ref\" \
 -B:input,VCF \"$variantRecalibratedVcf\" \
 -B:dbsnp,VCF \"$dbsnp\" \
 -l INFO \
 --fdr_filter_level 1.0 \
 -tranchesFile \"$tranchesFile\" \
 -o \"$variantRecalibratedFilteredVcf\" \
 -T ApplyVariantCuts \
;";
ApplyVariantCutsId=`echo "$cmd" | $qsub_lomem -N $name -hold_jid $VariantRecalibratorId`;
echo "$ApplyVariantCutsId $name $cmd";

#VariantEval
name="gatkVariantEval_`basename $variantEval`";
cmd=" \
java -Xmx5g -Djava.io.tmpdir=\"$TMPSCRATCH\" \
 -jar \"$GATK/GenomeAnalysisTK.jar\" \
 -T VariantEval \
 -R \"$ref\" \
 -l INFO \
 -B:eval,VCF \"$variantRecalibratedFilteredVcf\" \
 -B:dbsnp,VCF \"$dbsnp\" \
 -o \"$variantEval\" \
 -ST Sample \
 ;";
VariantEvalId=`echo "$cmd" | $qsub_lomem -N $name -hold_jid $ApplyVariantCutsId`;
echo "$VariantEvalId $name $cmd";

#VariantAnnotator
if (( ${#bamlist} != 0 )); then
	name="gatkVariantAnnotator_`basename $variantRecalibratedFilteredAnnotatedVcf`";
	cmd=" \
	java -Xmx5g -Djava.io.tmpdir=\"$TMPSCRATCH\" \
	 -jar \"$GATK/GenomeAnalysisTK.jar\" \
	 -T VariantAnnotator \
	 -l INFO \
	 -R \"$ref\" \
	 $bamlist \
	 -o \"$variantRecalibratedFilteredAnnotatedVcf\" \
	 --useAllAnnotations \
	 -B:variant,VCF \"$variantRecalibratedFilteredVcf\" \
	 -B:compHapMap,VCF \"$hapmap\" \
	 -B:comp1KG,VCF \"$kg\" \
	;";
	#VariantAnnotatorId=`echo "$cmd" | $qsub_lomem -N $name -hold_jid $ApplyVariantCutsId`;
	#echo "$VariantAnnotatorId $name $cmd";
fi

cd "$wd";

