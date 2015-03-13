#!/bin/bash

#download_clinvar_data.sh
USAGE="usage: $0 [/top/level/output/directory]"
if [[ $1 == "?" || $1 == "-h" || $1 == "--help" || $1 == "help" ]]; then echo ${USAGE}; exit; fi

#get path and parent directory of this script
#http://hintsforums.macworld.com/archive/index.php/t-73839.html
# if the call came via symlink, then use its target instead:
arg=$0; [[ -L $0 ]] && arg=$(stat -L -c'%n' "$0");
script_path=$(2>/dev/null cd "${arg%/*}" >&2; echo "`pwd -P`/${arg##*/}");
script_dir=$(dirname "$script_path");

HERE=${script_dir};

#output directory
if (( $# > 0 )); then d=$1; else d=${wd}; fi

DIR=${d};
mkdir -p $DIR;

curl "ftp://ftp.ncbi.nih.gov/pub/clinvar/gene_condition_source_id" -o ${DIR}/gene_condition_source_id.txt \

python /home/myourshaw/lab/pypeline/python2/sql_columns.py \
-db v -s clinvar -ci GeneSymbol \
-i ${DIR}/gene_condition_source_id.txt

curl "ftp://ftp.ncbi.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz" -o - \
| gunzip > ${DIR}/variant_summary.txt;

python /home/myourshaw/lab/pypeline/python2/sql_columns.py \
-db v -s clinvar -ci Chromosome Start Stop -x GeneSymbol \
-i ${DIR}/variant_summary.txt

