#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import argparse
from ConfigParser import SafeConfigParser
import re
from glob import iglob
import my
import job

#GMD4A.110623_SN860_0067_2011-100R_A81MVKABXX.3.1.TTAGGCA.novoalign.sam

for sam in iglob('/scratch0/tmp/myourshaw/gmd/bam_test/*.sam'):
	items = os.path.basename(sam).split('.')
	sample = items[0]
	exptrun = items[1]
	lane = items[2]
	barcode = items[4]
	rgid = '.'.join([exptrun,lane,barcode])
	metadata = my.get_metadata('/home/myourshaw/lab/git-myourshaw/python/metadata.txt', column_to_search='readgroup_id', value=rgid)
	if len(metadata) != 1:
		print metadata
	library = metadata[0]['library']
	newbase = '.'.join([sample,library,rgid,'novoalign','sam'])
	print 'os.rename({}, {})'.format(sam, os.path.join(os.path.dirname(sam), newbase))
	#os.rename(sam, os.path.join(os.path.dirname(sam), newbase))
