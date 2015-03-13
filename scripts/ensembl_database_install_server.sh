#!/bin/bash

#ensembl_database_install_server.sh
#runs on mysql server
#assumes mysql files downloaded and unzipped
USAGE="usage: $0 <ensembl current_mysql directory> <host> <user> <password>"
if [[ $1 == "?" || $1 == "-h" || $1 == "--help" || $1 == "help" ]]; then echo ${USAGE}; exit; fi

#get path and parent directory of this script
#http://hintsforums.macworld.com/archive/index.php/t-73839.html
# if the call came via symlink, then use its target instead:
arg=$0; [[ -L $0 ]] && arg=$(stat -L -c'%n' "$0");
script_path=$(2>/dev/null cd "${arg%/*}" >&2; echo "`pwd -P`/${arg##*/}");
script_dir=$(dirname "$script_path");

HERE=${script_dir};
wd=`pwd`;

if (( $# < 4 )) || [[ ! -d ${1} ]]; then echo ${USAGE}; exit; fi

ENSEMBL_DIR=${1}; shift;
HOST=${1}; shift;
USER=${1}; shift;
PASSWORD=${1}; shift;
MYSQL="mysql -h ${HOST} --user=${USER} -p${PASSWORD}";
MYSQLIMPORT="mysqlimport -h ${HOST} --user=${USER} -p${PASSWORD}";

#drop and re-create databases
echo 'dropping and recreating databases';
for d in ${ENSEMBL_DIR}/*; do
 if [[ ${d} == 'mysql' ]]; then
  echo "ERROR: accidentally tried to drop mysql database";
  exit;
 fi
 if [[ -d ${d} ]]; then
  db=${d##*/};
  echo "dropping database ${db}";
  echo "DROP DATABASE IF EXISTS ${db};" | ${MYSQL};
  echo "creating database ${db}";
  echo "CREATE DATABASE ${db};" | ${MYSQL};
  echo "granting permissions on ${db}";
  echo "GRANT SELECT ON ${db}.* TO 'ensembl'@'%';" | ${MYSQL};
 fi
done

#create tables
echo 'creating tables';
for d in ${ENSEMBL_DIR}/*; do
 db=${d##*/};
 if [[ -d ${d} && -e ${d}/${db}.sql ]]; then
  echo "creating tables in database ${db}";
  ${MYSQL} ${db} < ${d}/${db}.sql;
 else
  echo "${d}/${db}.sql does not exist";
 fi
done

echo "importing data";

#import data
for d in ${ENSEMBL_DIR}/*; do
 if [[ -d ${d} ]]; then
  db=${d##*/};
  filecount=$(find ${d} -name '*.txt' | wc -l);
  if (( ${filecount} > 0 )); then
   cmd="${MYSQLIMPORT} --local --fields_escaped_by=\\ ${db} ${d}/*.txt";
   echo $cmd;
   ${cmd};
  fi
 fi
done

echo 'ensembl database import done';

cd ${wd};

