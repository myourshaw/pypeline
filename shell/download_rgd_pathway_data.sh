#!/bin/bash

#download_rat_data.sh
#human, rat, and mouse pathways from RGD Pathways - Rat Genome Database
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

curl "ftp://rgd.mcw.edu/pub/data_release/annotated_rgd_objects_by_ontology/with_terms/homo_terms_pw" -o ${DIR}/rgd_homo_terms.txt;
curl "ftp://rgd.mcw.edu/pub/data_release/annotated_rgd_objects_by_ontology/with_terms/rattus_terms_pw" -o ${DIR}/rgd_rattus_terms.txt;
curl "ftp://rgd.mcw.edu/pub/data_release/annotated_rgd_objects_by_ontology/with_terms/mus_terms_pw" -o ${DIR}/rgd_mus_terms.txt;


for f in  ${DIR}/rgd_*_terms.txt; do
python /home/myourshaw/lab/pypeline/github-repositories/vax/VAX_DIR/mysql_table.py \
--header_line 28 -x 'OBJECT_SYMBOL' 'TERM_NAME' 'TERM_ACC_ID' \
-i ${f};
done

#create mysql database table and import data
#MYSQL="$(which mysql) -h ${HOST} -P ${PORT} --user=${USER} -p${PASSWORD}";
#MYSQLIMPORT="$(which mysqlimport) -h ${HOST} -P ${PORT} --user=${USER} -p${PASSWORD}";
#for f in ${DIR}/rgd_*_terms.mysql; do
#eval "${MYSQL} ${DATABASE} < ${f}";
#done
#for f in ${DIR}/rgd_*_terms.txt; do
#eval "${MYSQLIMPORT} --delete --local --ignore-lines=29 ${DATABASE} ${f}";
#done


