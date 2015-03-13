#!/usr/bin/env python

import sys
import os
import re
import traceback
from optparse import OptionParser, OptionGroup
import glob

parser = OptionParser()
(options, args) = parser.parse_args()
sep = "\t"

for file in args:
    file = os.path.realpath(file)
    trimmed = open(file + ".trimmed","w")
    for line in open(file):
        line = line.strip()
        if len(line) == 0:
            continue
        fields = line.split(sep)
        for i in range(len(fields)):
            fields[i] = fields[i].strip()
        trimmed.write(sep.join(fields)+"\n")
