#!/bin/bash

#download_hgnc_data.sh
USAGE="usage: $0 [/top/level/output/directory]"
if [[ $1 == "?" || $1 == "-h" || $1 == "--help" || $1 == "help" ]]; then echo ${USAGE}; exit; fi

#get path and parent directory of this script
#http://hintsforums.macworld.com/archive/index.php/t-73839.html
# if the call came via symlink, then use its target instead:
arg=$0; [[ -L $0 ]] && arg=$(stat -L -c'%n' "$0");
script_path=$(2>/dev/null cd "${arg%/*}" >&2; echo "`pwd -P`/${arg##*/}");
script_dir=$(dirname "$script_path");

HERE=${script_dir};
wd=`pwd`;

#output directory
if (( $# > 0 )); then d=$1; else d=${wd}; fi

DIR=$d;
mkdir -p ${DIR};

curl "ftp://ftp.ebi.ac.uk/pub/databases/genenames/hgnc_complete_set.txt.gz" -o - \
| gunzip > ${DIR}/hgnc.txt;


python ${HERE}/mysql_table.py \
--header_line 1 -t hgnc -k 'Approved Symbol' -x 'Approved Name' 'Ensembl ID (supplied by Ensembl)' \
-i ${DIR}/hgnc.txt;

#create mysql database table and import data
#MYSQL="$(which mysql) -h ${HOST} -P ${PORT} --user=${USER} -p${PASSWORD}";
#MYSQLIMPORT="$(which mysqlimport) -h ${HOST} -P ${PORT} --user=${USER} -p${PASSWORD}";
#eval "${MYSQL} ${DATABASE} < ${DIR}/hgnc.mysql";
#eval "${MYSQLIMPORT} --delete --local --ignore-lines=1 ${DATABASE} ${DIR}/hgnc.txt";



