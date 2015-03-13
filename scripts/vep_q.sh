#!/bin/bash

for f in $@; do
cmd="\
perl /share/apps/myourshaw/variant_effect_predictor-current/variant_effect_predictor.pl \
--no_progress \
--species human \
--input_file "${f}" \
--format vcf \
--output_file "${f}.vep" \
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
echo "${cmd}" | qsub -q all.q@compute-[23]* -cwd -V -M myourshaw@ucla.edu -m eas -terse -hard -pe serial 8 -l mem_free=24G -N vep_`basename ${f}`;
done

