#!/usr/bin/env python

import sys
import os
import re
import traceback
from optparse import OptionParser, OptionGroup
import glob

parser = OptionParser()
(options, args) = parser.parse_args()

def get_files(names):
	return (file for file in names if os.path.isfile(file))

#in_glob= "/Volumes/raid/rutgers/rutgers_map_b36/chr*.map2"
#chr_re = re.compile(r"chr(?P<chr>.+).map2",re.IGNORECASE)
#out = open("/Volumes/raid/rutgers/rutgers_map_b36/rutgers_map_b36.txt","w")
#header = "chr\tMarkers_name\tType\tPrimer/SNP_ref_name\tInformative_meioses\tHeterozygosity\tBuild36_map_physical_position\tSex-averaged_map_position\tFemale_map_position\tMale_map_position\n"
#sep = "\t"

#in_glob= "/Volumes/raid/rutgers/smooth_map_b36/chr*.sm.map2"
#chr_re = re.compile(r"chr(?P<chr>.+).sm.map2",re.IGNORECASE)
#out = open("/Volumes/raid/rutgers/smooth_map_b36/smooth_map_b36.txt","w")
#header = "chr\tMarkers_name\tBuild36_map_physical_position\tSex.averaged_map_position\tSmoothed_sex.averaged_map_position\tFemale_map_position\tSmoothed_Female_map_position\tMale_map_position\tSmoothed_male_map_position\n"
#sep = "\t"

in_glob = "/Users/myourshaw/Downloads/rutgers_map_b37/chr*.txt"
chr_re = re.compile(r"chr(?P<chr>\d{1,2}).txt",re.IGNORECASE)
out = open("/Users/myourshaw/Downloads/rutgers_map_b37/all.txt","w")
sep = "\t"
header = """chr	Marker_name	Type	Alias	Informative_meioses	Heterozygosity	Build37_map_physical_position	Sex-averaged_map_position	Female_map_position	Male_map_position	Build36_map_physical_position"""
header_fields = header.split(sep)


out.write(header+'\n')
seps_required = header.count("\t")-1
header_written = False

for file in get_files(glob.iglob(in_glob)):
	chr = chr_re.match(os.path.basename(file)).group("chr")
	line_number = 0
	for line in open(file):
		line_number += 1
		if line_number == 1:
			file_header_fields = line.rstrip('\n').split(sep)
			continue
		line_fields=line.rstrip('\n').split(sep)
		output_fields = []
		for i in range(1,len(header_fields)):
			if header_fields[i] in file_header_fields:
				output_fields.append(line_fields[file_header_fields.index(header_fields[i])])
			else:
				output_fields.append('')
		out.write("%s\t%s\n" % (chr,sep.join(output_fields)))
	