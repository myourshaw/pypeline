#!/usr/bin/env python
# -*- coding: utf-8 -*-

import my

file = '/data/storage-1-01/archive/myourshaw/scratch0/1000genomes/phase1_analysis_results_integrated_call_sets_20121012/ALL.wgs.integrated_phase1_v3.20101123.snps_indels_sv.sites.vcf.gz'
with my.open_gz_or_text(file) as vcf:
    that_line = ''
    that_chrom = ''
    that_pos = ''
    that_ref = ''
    that_alt = ''
    for line in vcf:
        if line.startswith('#'):
            continue
        else:
            fields = line.split('\t')
            chrom = fields[0]
            pos = fields[1]
            ref = fields[3]
            alt = fields[4]
            if chrom == that_chrom and pos == that_pos and ref == that_ref and len(alt) == 1 and len(that_alt) == 1:
                print that_line
                print line
            that_line = line
            that_chrom = chrom
            that_pos = pos
            that_ref = ref
            that_alt = alt
