#!/bin/sh

#All possible SNVs of GRCh37/hg19
#http://krishna.gs.washington.edu/download/CADD/v1.0/whole_genome_SNVs.tsv.gz ~2 hrs
#http://krishna.gs.washington.edu/download/CADD/v1.0/whole_genome_SNVs.tsv.gz.tbi

#All possible SNVs of GRCh37/hg19 incl. all annotations
#http://krishna.gs.washington.edu/download/CADD/v1.0/whole_genome_SNVs_inclAnno.tsv.gz ~5 hrs
#http://krishna.gs.washington.edu/download/CADD/v1.0/whole_genome_SNVs_inclAnno.tsv.gz.tbi

#12.3M InDels to initiate a local setup 
#http://krishna.gs.washington.edu/download/CADD/v1.0/InDels.tsv.gz
#http://krishna.gs.washington.edu/download/CADD/v1.0/InDels.tsv.gz.tbi

#12.3M InDels incl. all annotations to initiate a local setup
#http://krishna.gs.washington.edu/download/CADD/v1.0/InDels_inclAnno.tsv.gz
#http://krishna.gs.washington.edu/download/CADD/v1.0/InDels_inclAnno.tsv.gz.tbi

u="http://krishna.gs.washington.edu/download/CADD/v1.0/"

dir=/scratch1/vax/75/cadd;
mkdir -p ${dir};
qout=${dir}/qout;
mkdir -p ${qout};

for url in "${u}/whole_genome_SNVs.tsv.gz" "${u}/whole_genome_SNVs.tsv.gz.tbi" "${u}/whole_genome_SNVs_inclAnno.tsv.gz" "${u}/whole_genome_SNVs_inclAnno.tsv.gz.tbi" "${u}/InDels.tsv.gz" "${u}/InDels.tsv.gz.tbi" "${u}/InDels_inclAnno.tsv.gz" "${u}/InDels_inclAnno.tsv.gz.tbi"; do
cmd="wget -N -c "${url}" -P ${dir}";
echo "${cmd}" | qsub -q all.q -cwd -V -M myourshaw@ucla.edu -m eas -e ${qout} -o ${qout} -l excl=true -N $(basename ${url});
done

