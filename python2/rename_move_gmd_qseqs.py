#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import argparse
from ConfigParser import SafeConfigParser
import re
from glob import glob, iglob
import my
import job
import qseq2rg_novobarcode_and_cat
import qseq_metrics

in_glob = ['/scratch0/tmp/myourshaw/gmd/reads/readgroup_qseqs/*_qseq.txt.gz']
qseqs = [my.qseq_peek(q) for q in set(my.flatten([glob(q) for q in in_glob]))]
top_dir = '/scratch0/tmp/myourshaw/gmd/qseqs'
for q in qseqs:
    barcode = q['path'].rstrip('_qseq.txt.gz').split('.')[-1]
    new_name = '{}.{}.{}.{}.{}_qseq.txt.gz'.format(q['machine'], q['run'],q['lane'],q['read'], barcode)
    new_dir = os.path.join(top_dir, str(q['machine']), str(q['run']), str(q['lane']))
    new_file = os.path.join(new_dir, new_name)
    print 'mkdir -p {}; mv {} {}'.format(new_dir, q['path'], new_file)




