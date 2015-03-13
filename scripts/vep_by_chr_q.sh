#!/bin/bash

usage="usage: $0  [<email address>] <list of vcf files>";

arg=$0; [[ -L $0 ]] && arg=$(stat -L -c'%n' "$0")
script_path=$(2>/dev/null cd "${arg%/*}" >&2; echo "`pwd -P`/${arg##*/}")
script_dir=$(dirname "$script_path")

script_dir=$(dirname "$script_path")

HERE=${script_dir};

EMAIL="myourshaw@ucla.edu";
if [[ "$1" =~ '\S+@\S+\.\S+'  ]]; then EMAIL="$1"; shift; fi

source ${HERE}/b37.variables.sh $EMAIL;

for vcf in $@; do
	if [[ ! -e $vcf ]]; then echo "vcf file does not exist [$vcf]"; exit; fi
done

wd=`pwd`;

for vcf in $@; do
	#files created by this script
	output_dir=`dirname $vcf`;
	vcf_base=`basename ${vcf}`;
	mkdir -p $output_dir;
	tempdir="$output_dir/tmp/vep";
	qout="$tempdir/qout";
	mkdir -p $tempdir;
	mkdir -p $qout;
	
	job_ids=();
	split -a 6 -d -l 1000 ${vcf} ${tempdir}/${vcf_base}.part.
	for f in ${tempdir}/${vcf_base}.part.*; do
		vep=${f}.vep
		name=vep_`basename ${vep}`;
		cmd="\
		perl /share/apps/myourshaw/variant_effect_predictor-current/variant_effect_predictor.pl \
		--no_progress \
		--species human \
		--input_file "${f}" \
		--format vcf \
		--output_file "${vep}" \
		--force_overwrite \
		--host cortex.local \
		--user ensembl \
		--password ensembl \
		--port 3306 \
		--terms so \
		--sift=b \
		--polyphen=b \
		--condel=b \
		--regulatory \
		--hgvs \
		--gene \
		--protein \
		--hgnc \
		--dir ${HOME}/.vep \
		--buffer_size 10000 \
		";
		VEPId=`echo "${cmd}" | qsub -q all.q -cwd -V -e $qout -o $qout -M ${EMAIL} -m eas -terse -hard -pe serial 8 -l mem_free=6G -N ${name}`;
		echo "${VEPId} ${name}";
		job_ids+=(${VEPId});
	done
	
	jid_list=$(printf ",%s" "${job_ids[@]}");
	jid_list=${jid_list:1};
	
	#Cat temp files
	name=cat_tempfiles_`basename ${vcf}`;
	header=/scratch1/tmp/tpaige/resources/vep.header.txt;
	cmd=" \
	cat ${tempdir}/${vcf_base}.part.*.vep | grep -v '#' | cat $header - > $output_dir/${vcf_base}.vep
	;";
	CatTempFilesId=`echo "${cmd}" | qsub -q all.q -cwd -V -e $qout -o $qout -M ${EMAIL} -m eas -terse -N ${name} -hold_jid ${jid_list}`;
	echo "${CatTempFilesId} ${name}";
	
	#Remove tempfiles
	name=cleanup_`basename ${vcf}`;
	cmd="\
	rm ${tempdir}/${vcf_base}.part.*
	;";
	#CleanupFilesID=`echo "$cmd" | qsub -q all.q -cwd -V -e $qout -o $qout -M ${EMAIL} -m eas -terse -N ${name} -hold_jid ${CatTempFilesId}`;
	#echo "${CleanupFilesID} ${name}";



done

cd $wd
