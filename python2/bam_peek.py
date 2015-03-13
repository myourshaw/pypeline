#!/usr/bin/env python

#!/usr/bin/env python

import sys
import os
import tempfile
import glob
import gzip
from time import localtime, strftime
import math
import re
import collections
import argparse
from ConfigParser import SafeConfigParser
import csv
import my
import pysam


b = '/scratch0/tmp/myourshaw/gmd/reads/library_bams/pre-markdup/GMD125.GMD125_XT_10.bam'
bp = my.bam_peek(b)
for i in bp['RG']:
	for j in i:
		print '{}: {}\n'.format(j,i[j]),


