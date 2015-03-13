#!/bin/bash

#HUMAN BUILD GRCh37

USAGE="usage: source $0 [<email address> (default: myourshaw@ucla.edu)] [<resources directory> (default: $HERE/resources)]"
if [[ $1 == "?" || $1 == "-h" || $1 == "--help" || $1 == "help" ]]; then echo $USAGE; exit; fi

#get path and parent directory of this script
#http://hintsforums.macworld.com/archive/index.php/t-73839.html
# if the call came via symlink, then use its target instead:
arg=$0; [[ -L $0 ]] && arg=$(stat -L -c'%n' "$0")
script_path=$(2>/dev/null cd "${arg%/*}" >&2; echo "`pwd -P`/${arg##*/}")
script_dir=$(dirname "$script_path")

HERE=${script_dir};

EMAIL="myourshaw@ucla.edu";
if [[ "$1" =~ '\S+@\S+\.\S+'  ]]; then EMAIL="$1"; shift; fi

RESOURCES_DIR="${HERE}/resources";
if [[ -d "$1" ]]; then RESOURCES_DIR="$1"; shift; fi

QSUB_ALLQ='qsub -q all.q -cwd -V -M '$EMAIL' -m eas -terse';
QSUB_ALLQ_2GB='qsub -q all.q -cwd -V -M '$EMAIL' -m eas -terse -hard -l mem_free=2G';
QSUB_ALLQ_4GB='qsub -q all.q -cwd -V -M '$EMAIL' -m eas -terse -hard --l mem_free=4G';
QSUB_ALLQ_8='qsub -q all.q -cwd -V -M '$EMAIL' -m eas -terse-hard -pe serial 8';
QSUB_ALLQ_2GB_8='qsub -q all.q -cwd -V -M '$EMAIL' -m eas -terse -hard -pe serial 8 -l mem_free=2G';
QSUB_ALLQ_4GB_8='qsub -q all.q -cwd -V -M '$EMAIL' -m eas -terse -hard -pe serial 8 --l mem_free=4G';
QSUB_LOMEM='qsub -q all.q@compute-4* -cwd -V -M '$EMAIL' -m eas -terse -hard -pe serial 8 -l mem_free=6G';
QSUB_HIMEM='qsub -q all.q@compute-[23]* -cwd -V -M '$EMAIL' -m eas -terse -hard -pe serial 8 -l mem_free=24G';

#APPLICATIONS
APPS="/share/apps/myourshaw";
PYTHON="/share/apps/python/2.7.1/bin/python";
if [[ ! -x ${PYTHON} ]]; then PYTHON="python"; fi

#novoalign
NOVOCRAFT="${APPS}/novocraft-current";
NOVOALIGN="${NOVOCRAFT}/novoalign";
NOVOBARCODE="${NOVOCRAFT}/novobarcode";
NOVOINDEX="${NOVOCRAFT}/novoindex";

#GATK
STING="${APPS}/Sting";
GATK="${STING}/dist";
GATK_PERL="${STING}/perl";
GATK_PYTHON="${STING}/python";
GATK_R="${STING}/R";
GATK_R_SCRIPTS="${STING}/R";

#satools
SAMTOOLS="${APPS}/samtools";
SAMTOOLS_DIR="${APPS}/samtools-current";

#picard
PICARD="${APPS}/picard-current";
PICARD_DIR="${APPS}/picard-current";
PICARD_METRICS="${APPS}/scripts/picard_metrics_q.sh";
PICARD_METRICS_CONSOLIDATE="${APPS}/scripts/python/picard_metrics_consolidate.py";

#R
RSCRIPT=`which Rscript 2>/dev/null`;
if [[ ! -x $RSCRIPT ]]; then RSCRIPT=/usr/bin/Rscript; fi

#ensembl
ENSEMBL_DB="";
VARIANT_EFFECT_PREDICTOR="${APPS}/variant_effect_predictor/variant_effect_predictor.pl";

#dindel
DINDEL=`which dindel 2>/dev/null`;
if [[ ! -x $DINDEL ]]; then DINDEL="/share/apps/myourshaw/dindel-current/binaries/dindel-linux-x64-novoalign"; fi
DINDEL_PYTHON="/share/apps/myourshaw/dindel-python";

#TEMPORARY
TMPLOCAL="/state/partition1/tmp/`whoami`";
TMPSCRATCH="/scratch1/tmp/`whoami`/tmp";
mkdir -p $TMPSCRATCH;

#RESOURCES

#REFERENCE GENOME
B37_FASTA="${RESOURCES_DIR}/human_g1k_v37.fasta";
B37_FASTA_FAI="${RESOURCES_DIR}/human_g1k_v37.fasta.fai";
NOVOALIGN_MAX_INDEX="${RESOURCES_DIR}/ref.b37.k15s1.nix";
NOVOALIGN_HIMEM_INDEX="${RESOURCES_DIR}/ref.b37.himem.nix";
NOVOALIGN_LOMEM_INDEX="${RESOURCES_DIR}/ref.b37.lomem.nix";
#NOVOALIGNCS_MAX_INDEX="${RESOURCES_DIR}/ref.b37.k15s1.ncx";
#NOVOALIGNCS_HIMEM_INDEX="${RESOURCES_DIR}/ref.b37.himem.ncx";
#NOVOALIGNCS_LOMEM_INDEX="${RESOURCES_DIR}/ref.b37.lomem.ncx";

#GATK GENE MODELS
REFGENE_B37_SORTED="${RESOURCES_DIR}/refgene.b37.sorted.txt"
REFGENE_B37_ROD="${RESOURCES_DIR}/refGene_b37.rod";
#CCDS_BIG_TABLE="${RESOURCES_DIR}/";
KNOWNGENE_BIG_TABLE="${RESOURCES_DIR}/knownGene-big-table-b37.txt";
REFGENE_BIG_TABLE="${RESOURCES_DIR}/refGene-big-table-b37.txt";

#VARIANTS
DBSNP_129_B37_ROD="${RESOURCES_DIR}/var.dbsnp129.b37.rod";
#DBSNP_132_B37_ROD="${RESOURCES_DIR}/";
DBSNP_132_B37_VCF="${RESOURCES_DIR}/var.dbsnp132.b37.vcf";
DBSNP_132_EXCLUDING_SITES_AFTER_129_VCF="${RESOURCES_DIR}/dbsnp_132.excluding_sites_after_129.vcf";
HAPMAP_R27_B37_VCF="${RESOURCES_DIR}/var.hapmap_r27.b37.vcf";
ONEKG_AUTOSOME_B37_VCF="${RESOURCES_DIR}/var.1kg.autosome.b37.vcf";
ONEKG_CHRX_B37_VCF="${RESOURCES_DIR}/var.1kg.chrX.b37.vcf";
BUSHMAN_B37_VCF="${RESOURCES_DIR}/var.bushman.b37.vcf";
DANE_B37_VCF="${RESOURCES_DIR}/var.dane.b37.vcf";
TIBHAN_B37_VCF="${RESOURCES_DIR}/var.tibhan.b37.vcf";

#intervals
INTERVALS_DIR="${RESOURCES_DIR}/intervals";
SURESELECT_INTERVALS_DIR="${INTERVALS_DIR}/SureSelect"; #agilent
TRUSEQ_INTERVALS_DIR="${INTERVALS_DIR}/"; #illumina
SEQCAP_INTERVALS_DIR="${INTERVALS_DIR}/"; #nimblegen
ENSEMBL_INTERVALS_DIR="${INTERVALS_DIR}/EnsemblIntervals";
GENE_INTERVALS_DIR="${INTERVALS_DIR}/GeneIntervals":

#BAIT CAPTURE REGIONS
SURESELECT_38MB_BED="${SURESELECT_INTERVALS_DIR}/SureSelect_All_Exon_G3362_with_annotation_b37_sorted.bed";
SURESELECT_38MB_INTERVAL_LIST="${SURESELECT_INTERVALS_DIR}/SureSelect_All_Exon_G3362_with_annotation_b37_sorted.interval_list";
SURESELECT_38MB_V2_BED="${SURESELECT_INTERVALS_DIR}/SureSelect_All_Exon_V2_with_annotation_b37_sorted.bed";
SURESELECT_38MB_V2_INTERVAL_LIST="${SURESELECT_INTERVALS_DIR}/SureSelect_All_Exon_V2_with_annotation_b37_sorted.interval_list";
SURESELECT_50MB_BED="SURESELECT_INTERVALS_DIR/SureSelect_All_Exon_50mb_with_annotation_b37_sorted.bed";
SURESELECT_50MB_INTERVAL_LIST="${SURESELECT_INTERVALS_DIR}/SureSelect_All_Exon_50mb_with_annotation_b37_sorted.interval_list";
SURESELECT_50MB_REFSEQ_INTERSECTED_EXTENDED2BP_INTERVAL_LIST="${SURESELECT_INTERVALS_DIR}/SureSelect_Sanger_b37_100bp_merged_intersected_RefSeqG37_extended2bp.interval_list";

#TARGET WISH LIST GENE MODEL REGIONS
ENSEMBL_ALL_GENES="${ENSEMBL_INTERVALS_DIR}/all_genes.interval_list";
ENSEMBL_CLINICAL_GENES="${ENSEMBL_INTERVALS_DIR}/clinical_genes.interval_list";
ENSEMBL_PROTEIN_CODING_GENES="${ENSEMBL_INTERVALS_DIR}/protein_coding_genes.interval_list";
ALLGENES_14MODEL_CDS_BED="${GENE_INTERVALS_DIR}/intervals_translated_cds/b37.allGenes14ModelsMerged.cds.bed";
ALLGENES_14MODEL_CDS_INTERVAL_LIST="${GENE_INTERVALS_DIR}/intervals_translated_cds/b37.allGenes14ModelsMerged.cds.interval_list";

