#!/usr/bin/env python

import MySQLdb
database = MySQLdb.connect (host = 'myourshaw-dev.genome.ucla.edu', db = 'variants')
cursor =database.cursor (MySQLdb.cursors.DictCursor)

cursor.execute (
"""SELECT sample, chrom, pos, ref, alt
FROM genotypes
WHERE genotype = '1/1' and rs_id is null
ORDER BY sample, chrom, pos'""")

genotypes = cursor.fetchall()

for genotype in genotypes:
	print genotype

for sample in samples:
	picard.run('MergeSamFiles', INPUT=sample.readgroup_bams, OUTPUT=sample.id+'.bam')
