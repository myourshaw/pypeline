#!/bin/bash

#ensembl_table_reload.sh
#assumes databases and tables exist
USAGE="usage: $0 <ensembl current_mysql directory> <host> <user> <password> <list of database.table names>"
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

ENSEMBL_DIR=${1}; shift;
HOST=${1}; shift;
USER=${1}; shift;
PASSWORD=${1}; shift;
TABLES=$@;
MYSQL="/usr/bin/mysql -h ${HOST} --user=${USER} -p${PASSWORD}";
MYSQLIMPORT="/usr/bin/mysqlimport -h ${HOST} --user=${USER} -p${PASSWORD}";

SUFFIX="`date +%Y%m%d%H%M%S`-$RANDOM"
qout=${wd}/qout_${SUFFIX};
mkdir -p ${qout};
QSUB='qsub -q all.q@compute-4*'" -cwd -V -e ${qout} -o ${qout} -l excl=true -N ensembl_table_reload_${SUFFIX}";
cmds="${qout}/ensembl_table_reload_${SUFFIX}.cmds";
if [[ -f "$cmds" ]]; then rm -f "$cmds"; fi

#truncate tables
for dt in ${TABLES}; do
 db=${dt%.*};
 d=${ENSEMBL_DIR}/${db};
 t=${dt#*.};
 if [[ -d ${d} ]] && [[ -e ${d}/${db}.sql ]]; then
  echo "truncating table ${db}.${t}";
  echo "TRUNCATE TABLE ${db}.${t};" | ${MYSQL};
 else
  echo "${d}/${db}.sql does not exist";
 fi
done

echo "creating job array to import data";

#import data
for dt in ${TABLES}; do
 db=${dt%.*};
 d=${ENSEMBL_DIR}/${db};
 t=${dt#*.};
 f=${d}/${t}.txt
 if [[ -d ${d} && -e ${f} ]]; then
  N="${db}.`basename ${f%.txt}`";
  cmd="${MYSQLIMPORT} -L --fields_escaped_by=\\ ${db} ${f}";
  echo ${cmd} >> ${cmds};
 else
  echo "${f} does not exist";
 fi
done

if [[ -f "$cmds" ]]; then
 echo '$(head -n $SGE_TASK_ID' ${cmds} '| tail -n 1)' | ${QSUB} -t 1-$(cat ${cmds} | wc -l);
fi

echo "running cluster jobs to import data";
echo "job array is in ${cmds}";
echo "STDOUT and STDERR will be in ${qout}";

cd ${wd};

