#!/bin/bash

usage="usage: $0 [<email address> (default: myourshaw@ucla.edu)] <novobarcode executable> <qseq tag file per novobarcode spec> [options for novobarcode] <list of first paired-end qseq files [the program will match the second paired-end file and the barcode file]>";
if (( $# < 2 )); then echo $usage; exit; fi

#get path and parent directory of this script
#http://hintsforums.macworld.com/archive/index.php/t-73839.html
# if the call came via symlink, then use its target instead:
arg=$0; [[ -L $0 ]] && arg=$(stat -L -c'%n' "$0")
script_path=$(2>/dev/null cd "${arg%/*}" >&2; echo "`pwd -P`/${arg##*/}")
script_dir=$(dirname "$script_path")

EMAIL="myourshaw@ucla.edu";
if [[ "$1" =~ '\S+@\S+\.\S+'  ]]; then EMAIL="$1"; shift; fi

QSUB_ALLQ='qsub -q all.q -cwd -V -M '$EMAIL' -m eas -terse';
QSUB_LOMEM='qsub -q all.q@compute-4*,all.q@compute-[23]* -cwd -V -M '$EMAIL' -m eas -terse -hard -pe serial 8 -l mem_free=6G';
QSUB_HIMEM='qsub -q all.q@compute-[23]* -cwd -V -M '$EMAIL' -m eas -terse -hard -pe serial 8 -l mem_free=24G';

wd=`pwd`;

novobarcode=$1; shift;
if [[ ! -x "$novobarcode" ]]; then echo "novobarcode program $novobarcode doesn't exist or isn't executable"; echo $usage; exit; fi

tags=$1; shift;
if [[ ! -e "$tags" ]]; then echo "qseq tag file $tags doesn't exist"; echo $usage; exit; fi

options="";
while [[ $1 =~ ^- ]]; do
 options="$options $1"; shift;
done

if  (( $# < 1 )); then echo "input must include one or more qseq files containing first paired end reads"; echo $usage; exit; fi

for pe1 in $@; do
	if [[ ! -e "$pe1" ]]; then echo "first qseq file $pe1 doesn't exist"; exit; fi
	pe2=${pe1%1_*_qseq.txt}3${pe1##*s_?_1};
	if [[ ! -e "$pe2" ]]; then echo "second qseq file $pe2 doesn't exist"; exit; fi
	bc=${pe1%1_*_qseq.txt}2${pe1##*s_?_1};
	if [[ ! -e "$bc" ]]; then echo "barcode qseq file $bc doesn't exist"; exit; fi
done

for pe1 in $@; do
pe2=${pe1%1_*_qseq.txt}3${pe1##*s_?_1};
bc=${pe1%1_*_qseq.txt}2${pe1##*s_?_1};
name="unbarcode_`basename $pe1`";
cmd="
 $novobarcode \
  -b $tags \
  -F QSEQ \
  -f $pe1 $pe2 \
  -i $bc \
  --ILQ_SKIP \
  $options \
 ;";

 UnbarcodeId=`echo "$cmd" | $QSUB_ALLQ -N $name`;
 echo "$UnbarcodeId $name $cmd";

done
cd $wd;
