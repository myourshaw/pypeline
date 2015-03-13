#!/bin/sh

#import_uniprot_mysql.sh 
#assumes data directory populated by download_uniprot_data.sh and uniprot2db.pl
USAGE="usage: $0 /uniprot/download/directory host user password database"
if [[ $1 == "?" || $1 == "-h" || $1 == "--help" || $1 == "help" ]]; then echo ${USAGE}; exit; fi

#get path and parent directory of this script
#http://hintsforums.macworld.com/archive/index.php/t-73839.html
# if the call came via symlink, then use its target instead:
arg=$0; [[ -L $0 ]] && arg=$(stat -L -c'%n' "$0");
script_path=$(2>/dev/null cd "${arg%/*}" >&2; echo "`pwd -P`/${arg##*/}");
script_dir=$(dirname "$script_path");

HERE=${script_dir};
wd=`pwd`;

if (( $# < 5 )) || [[ ! -d ${1} ]]; then echo ${USAGE}; exit; fi

UNIPROT_DIR=${1}; shift;
HOST=${1}; shift;
USER=${1}; shift;
PASSWORD=${1}; shift;
DATABASE=${1}; shift;
MYSQL="mysql -h ${HOST} --user=${USER} -p${PASSWORD} --execute";
MYSQLIMPORT="mysqlimport -h ${HOST} --user=${USER} -p${PASSWORD}";

table=uniprot_human_protein;
COLLIST=`head -1 ${UNIPROT_DIR}/${table}.txt`
${MYSQL} "TRUNCATE TABLE ${DATABASE}.${table}";
${MYSQL} "LOAD DATA LOCAL INFILE '${UNIPROT_DIR}/${table}.txt' INTO TABLE ${DATABASE}.${table} IGNORE 1 LINES; SHOW WARNINGS;" > ${UNIPROT_DIR}/${table}.output;
table=uniprot_human_xref;
COLLIST=`head -1 ${UNIPROT_DIR}/${table}.txt`
${MYSQL} --execute "TRUNCATE TABLE ${DATABASE}.${table}";
${MYSQL} "LOAD DATA LOCAL INFILE '${UNIPROT_DIR}/${table}.txt' INTO TABLE ${DATABASE}.${table} IGNORE 1 LINES; SHOW WARNINGS;" > ${UNIPROT_DIR}/${table}.output;
table=uniprot_human_feature;
COLLIST=`head -1 ${UNIPROT_DIR}/${table}.txt`
${MYSQL} "TRUNCATE TABLE ${DATABASE}.${table}";
${MYSQL} "LOAD DATA LOCAL INFILE '${UNIPROT_DIR}/${table}.txt' INTO TABLE ${DATABASE}.${table} IGNORE 1 LINES; SHOW WARNINGS;" > ${UNIPROT_DIR}/${table}.output;
