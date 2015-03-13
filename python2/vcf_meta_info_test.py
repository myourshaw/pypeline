#!/usr/bin/env python

import my

vcf_file = '/scratch1/vax/75/dbsnp/00-All.vcf.gz'
p = my.parse_vcf_meta_info(vcf_file)
v = p.get('dbSNP_BUILD_ID')
print p
