#!/bin/bash
#$ -cwd
#$ -V
#$ -M myourshaw@ucla.edu
#$ -m eas
#$ -terse

#gatk_recalibrate_Illumina_srma_q.sh

usage="usage: $0 [<email address> (default: myourshaw@ucla.edu)] [<ref.fasta> | <ref.fa> (default human_g1k_v37.fasta)] <output dir for coverage> <list of sorted bam filest>";
here=$(dirname "$0");

email="myourshaw@ucla.edu";
if [[ "$1" =~ '\S+@\S+\.\S+'  ]]; then email="$1"; shift; fi

source "$here/b37.variables.sh" $email;
ref=$B37_FASTA;
if [[ ${1##*.} == "fasta" || ${1##*.} == "fa" ]]; then ref="$1"; shift; fi
if (( $# < 2 )); then echo $usage; exit; fi

coveragedir="$1"; shift;
mkdir -p "$coveragedir";
for bam in $@; do
	if [[ ! -e "$bam" ]]; then echo "bam file does not exist [$bam]"; exit; fi
done

wd="`pwd`";

#resources
dbsnp="$DBSNP_132_B37_VCF";
refgene="$REFGENE_B37_ROD";
targets="$SURESELECT_50MB_BED";

#applications and resources
sting="$STING";
gatk="$GATK";
gatkrscripts="$GATK_R_SCRIPTS";
perl="$GATK_PERL";
python="$GATK_PYTHON";
rscript="$RSCRIPT";
picard="$PICARD";
tmplocal="$TMPLOCAL";
tmpscratch="$TMPSCRATCH";
mkdir -p "$tmpscratch";

#job control lists
holdall='';
recalibratedlist='';
validatelist='';

#analyze bam files
for bam in $@; do
	#output files
	output_dir="`dirname "$bam"`";
	mkdir -p "$output_dir";
	qout="$output_dir/qout";
	mkdir -p "$qout";
	qsub_allq="$QSUB_ALLQ -e $qout -o $qout";
	qsub_lomem="$QSUB_LOMEM -e $qout -o $qout";
	qsub_himem="$QSUB_HIMEM -e $qout -o $qout";
	intervals="${bam}.intervals";
	realigned="${bam%.bam}.realigned.bam";
	fixmate="${realigned%.bam}.fixmate.bam";
	markdup="${fixmate%.bam}.markdup.bam";
	markdupmetrics="$markdup.metrics";
	recal1csv="${markdup%.bam}.recal.csv";
	recalibrated="${markdup%.bam}.recalibrated.bam";
	recal2csv="${recalibrated%.bam}.recal.csv";
	validate="${recalibrated}.validate";
	recalibratedlist="${recalibratedlist} -I ${recalibrated}";
	validatelist="${validatelist} $validate";
	analyzeCovariatesOutputDirPre="`dirname "$bam"`/gatkRecalibrationCovariateAnalysis/pre-recalibration";
	mkdir -p "$analyzeCovariatesOutputDirPre";
	analyzeCovariatesOutputDirPost="`dirname "$bam"`/gatkRecalibrationCovariateAnalysis/post-recalibration";
	mkdir -p "$analyzeCovariatesOutputDirPost";	
	
	#RealignerTargetCreator
	name="gatkRealignerTargetCreator_`basename $intervals`";
	cmd=" \
	java -Xmx5g  -Djava.io.tmpdir=\"$TMPSCRATCH\" \
	 -jar \"$GATK/GenomeAnalysisTK.jar\" \
	 -T RealignerTargetCreator \
	 -R \"$ref\" \
	 -o \"$intervals\" \
	 -I \"$bam\" \
	 -B:dbsnp,VCF \"$dbsnp\" \
	 --mismatchFraction 0.10 \
	 --realignReadsWithBadMates \
	 ;";
	RealignerTargetCreatorId=`echo "$cmd" | $qsub_lomem -N $name`;
	echo "$RealignerTargetCreatorId $name $cmd";
	
	#IndelRealigner
	name="gatkIndelRealigner_`basename $realigned`";
	cmd=" \
	java -Xmx5g -Djava.io.tmpdir=\"$TMPSCRATCH\" \
	 -jar \"$GATK/GenomeAnalysisTK.jar\" \
	 -I \"$bam\" \
	 -R \"$ref\" \
	 -T IndelRealigner \
	 -targetIntervals \"$intervals\" \
	 -o \"$realigned\" \
	 -B:dbsnp,VCF \"$dbsnp\" \
	 -compress 0 \
	 --entropyThreshold 0.10 \
	 --maxReadsInRam 1000000 \
	 --realignReadsWithBadMates \
	 ;";
	IndelRealignerId=`echo "$cmd" | $qsub_lomem -N $name -hold_jid $RealignerTargetCreatorId`;
	echo "$IndelRealignerId $name $cmd";
	
	#FixMateInformation
	name="gatkFixMateInformation_`basename $fixmate`";
	cmd=" \
	java -Xmx5g -Djava.io.tmpdir=\"$TMPSCRATCH\" \
	 -jar \"$picard/FixMateInformation.jar\" \
	 TMP_DIR=\"$TMPSCRATCH\" \
	 VERBOSITY=INFO \
	 QUIET=false \
	 VALIDATION_STRINGENCY=SILENT \
	 COMPRESSION_LEVEL=5 \
	 MAX_RECORDS_IN_RAM=1000000 \
	 CREATE_INDEX=true \
	 CREATE_MD5_FILE=true \
	 INPUT=\"$realigned\" \
	 OUTPUT=\"$fixmate\" \
	 SORT_ORDER=coordinate \
	;";
	FixMateInformationId=`echo "$cmd" | $qsub_lomem -N $name -hold_jid $IndelRealignerId`;
	echo "$FixMateInformationId $name $cmd";

	#MarkDuplicates
	#SOLiD PE READ_NAME_REGEX=([0-9]+)_([0-9]+)_([0-9]+).*
	#SOLiD PE OPTICAL_DUPLICATE_PIXEL_DISTANCE=10
	#Illumina PE READ_NAME_REGEX=.+:.+:[0-9]+:[0-9]+:([0-9]+):([0-9]+):([0-9]+).*
	#Illumina PE OPTICAL_DUPLICATE_PIXEL_DISTANCE=100
	name="gatkMarkDuplicates_`basename $markdup`";
	cmd=" \
	java -Xmx24g -Djava.io.tmpdir=\"$TMPSCRATCH\" \
	-jar \"$PICARD/MarkDuplicates.jar\" \
 	 TMP_DIR=\"$TMPSCRATCH\" \
	 VERBOSITY=INFO \
	 QUIET=false \
	 VALIDATION_STRINGENCY=SILENT \
	 COMPRESSION_LEVEL=5 \
	 MAX_RECORDS_IN_RAM=1000000 \
	 CREATE_INDEX=true \
	 CREATE_MD5_FILE=true \
	 INPUT=\"$fixmate\" \
	 OUTPUT=\"$markdup\" \
	 METRICS_FILE=\"$markdupmetrics\" \
	 REMOVE_DUPLICATES=false \
	 READ_NAME_REGEX=\".+:.+:[0-9]+:[0-9]+:([0-9]+):([0-9]+):([0-9]+).*\" \
	 OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 \
	 ;";
	MarkDuplicatesId=`echo "$cmd" | $qsub_himem -N $name -hold_jid $FixMateInformationId`;
	echo "$MarkDuplicatesId $name $cmd";
	
	#CountCovariates1 before recalibration
	name="gatkCountCovariates1_`basename $recal1csv`";
	cmd=" \
	java -Xmx5g -Djava.io.tmpdir=\"$TMPSCRATCH\" \
	-jar \"$GATK/GenomeAnalysisTK.jar\" \
	 -l INFO \
	 -R \"$ref\" \
	 -B:dbsnp,VCF \"$dbsnp\" \
	 -I \"$markdup\" \
	 -T CountCovariates \
	  -cov ReadGroupCovariate \
	  -cov QualityScoreCovariate \
	  -cov CycleCovariate \
	  -cov DinucCovariate \
	 -recalFile \"$recal1csv\" \
	 -nt 8 \
	 ;";
	CountCovariates1Id=`echo "$cmd" | $qsub_lomem -N $name -hold_jid $MarkDuplicatesId`;
	echo "$CountCovariates1Id $name $cmd";
	
	#TableRecalibration
	name="gatkTableRecalibration_`basename $recalibrated`";
	cmd=" \
	java -Xmx5g -Djava.io.tmpdir=\"$TMPSCRATCH\" \
	-jar \"$GATK/GenomeAnalysisTK.jar\" \
	 -l INFO \
	 -R \"$ref\" \
	 -I \"$markdup\" \
	 -T TableRecalibration \
	 -o \"$recalibrated\" \
	 -recalFile \"$recal1csv\" \
	 -compress 0 \
	 -baq RECALCULATE \
	 ;";
	TableRecalibrationId=`echo "$cmd" | $qsub_lomem -N $name -hold_jid $CountCovariates1Id`;
	echo "$TableRecalibrationId $name $cmd";
	
	#BuildRecalibratedBamIndex
	name="gatkBuildRecalibratedBamIndex_`basename $recalibrated`";
	cmd=" \
	java -Xmx5g -Djava.io.tmpdir=\"$TMPSCRATCH\" \
	-jar \"$picard/BuildBamIndex.jar\" \
	 INPUT=\"$recalibrated\" \
	 TMP_DIR=\"$TMPSCRATCH\" \
	 VERBOSITY=INFO \
	 VALIDATION_STRINGENCY=SILENT \
	 QUIET=false \
	 ;";
	BuildRecalibratedBamIndexId=`echo "$cmd" | $qsub_lomem -N $name -hold_jid $TableRecalibrationId`;
	echo "$BuildRecalibratedBamIndexId $name $cmd";
	
	#CountCovariates2 again on the recalibrated bam file
	name="gatkCountCovariates2_`basename $recal2csv`";
	cmd=" \
	java -Xmx5g -Djava.io.tmpdir=\"$TMPSCRATCH\" \
	 -jar \"$GATK/GenomeAnalysisTK.jar\" \
	 -l INFO \
	 -R \"$ref\" \
	 -B:dbsnp,VCF \"$dbsnp\" \
	 -I \"$recalibrated\" \
	 -T CountCovariates \
	  -cov ReadGroupCovariate \
	  -cov QualityScoreCovariate \
	  -cov CycleCovariate \
	  -cov DinucCovariate \
	 -recalFile \"$recal2csv\" \
	 -nt 8 \
	 ;";
	CountCovariates2Id=`echo "$cmd" | $qsub_lomem -N $name -hold_jid $BuildRecalibratedBamIndexId`;
	echo "$CountCovariates2Id $name $cmd";
	
	#AnalyzeCovariates1 before recalibration
	name="gatkAnalyzeCovariates1_`basename $recal1csv`";
	cmd=" \
	java -Xmx5g -Djava.io.tmpdir=\"$TMPSCRATCH\" \
	 -jar \"$GATK/AnalyzeCovariates.jar\" \
	 -recalFile \"$recal1csv\" \
	 -Rscript \"$rscript\" \
	 -resources \"$gatkrscripts\" \
	 -outputDir \"$analyzeCovariatesOutputDirPre\" \
	 -ignoreQ 5 \
	 ;";
	AnalyzeCovariates1Id=`echo "$cmd" | $qsub_lomem -N $name -hold_jid $CountCovariates1Id`;
	echo "$AnalyzeCovariates1Id $name $cmd";
	
	#AnalyzeCovariates2 after recalibration
	name="gatkAnalyzeCovariates2_`basename $recal2csv`";
	cmd=" \
	java -Xmx5g -Djava.io.tmpdir=\"$TMPSCRATCH\" \
	 -jar \"$GATK/AnalyzeCovariates.jar\" \
	 -recalFile \"$recal2csv\" \
	 -Rscript \"$rscript\" \
	 -resources \"$gatkrscripts\" \
	 -outputDir \"$analyzeCovariatesOutputDirPost\" \
	 -ignoreQ 5 \
	 ;";
	AnalyzeCovariates2Id=`echo "$cmd" | $qsub_lomem -N $name -hold_jid $CountCovariates2Id`;
	echo "$AnalyzeCovariates2Id $name $cmd";
	
	#ValidateSamFile
	name="gatkValidateSamFile_`basename $recalibrated`";
	cmd=" \
	java -Xmx5g  -Djava.io.tmpdir=\"$TMPSCRATCH\" \
	 -jar \"$picard/ValidateSamFile.jar\" \
 	 TMP_DIR=\"$TMPSCRATCH\" \
	 VERBOSITY=INFO \
	 QUIET=false \
	 VALIDATION_STRINGENCY=SILENT \
	 COMPRESSION_LEVEL=5 \
	 MAX_RECORDS_IN_RAM=1000000 \
	 INPUT=\"$recalibrated\" \
	 OUTPUT=\"$validate\" \
	 MODE=VERBOSE \
	 REFERENCE_SEQUENCE=\"$ref\" \
	 IGNORE=MISSING_TAG_NM \
	 ;";
	ValidateSamFileId=`echo "$cmd" | $qsub_lomem -N $name -hold_jid $BuildRecalibratedBamIndexId`;
	echo "$ValidateSamFileId $name $cmd";

	if [[ $holdall == '' ]]; then holdall=$ValidateSamFileId; else holdall="${holdall},${ValidateSamFileId}"; fi
done

####################################################################################################
# holdall can be used to trigger downstream processing that uses all bam files ($recalibratedlist) as inputs
####################################################################################################
#echo "please check results of bam file validation before proceeding with downstream analysis: $validatelist";

#DepthOfCoverage
name="gatkDepthOfCoverage_`basename $coveragedir/coverage.out`";
cmd=" \
java -Xmx24g -Djava.io.tmpdir=\"$TMPSCRATCH\" \
 -jar \"$GATK/GenomeAnalysisTK.jar\" \
 -T DepthOfCoverage \
 -l INFO \
 -R \"$ref\" \
 $recalibratedlist \
 -o \"$coveragedir/coverage.out\" \
 -L \"$targets\" \
 -geneList \"$refgene\" \
;";
#DepthOfCoverageId=`echo "$cmd" | $qsub_himem -N $name -hold_jid $holdall`;
#echo "$DepthOfCoverageId $name $cmd";

cd "$wd";

