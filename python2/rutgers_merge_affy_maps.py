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

in_glob = "/Volumes/raid/rutgers/Affymetrix_50/Affy_50_*M.csv"
chr_re = re.compile(r"Affy_50_(?P<chr>.+)M.csv",re.IGNORECASE)
out = open("/Volumes/raid/rutgers//Affymetrix_50//Affymetrix_50.txt","w")
header = "chr\tMarkers_name\tBuild36_map_physical_position\tSex-averaged_map_position_(Kosambi_cM)\tFemale_map_position_(Kosambi_cM)\tMale_map_position_(Kosambi_cM)\n"
sep = ","

#in_glob = "/Volumes/raid/rutgers/Affymetrix_60/Affy_60_*M.csv"
#chr_re = re.compile(r"Affy_60_(?P<chr>.+)M.csv",re.IGNORECASE)
#out = open("/Volumes/raid/rutgers//Affymetrix_60//Affymetrix_60.txt","w")
#header = "chr\tMarkers_name\tBuild36_map_physical_position\tSex-averaged_map_position_(Kosambi_cM)\tFemale_map_position_(Kosambi_cM)\tMale_map_position_(Kosambi_cM)\n"
#sep = ","

out.write(header)
seps_required = header.count("\t")-1
header_written = False

for file in get_files(glob.iglob(in_glob)):
	chr = chr_re.match(os.path.basename(file)).group("chr")
	line_number = 0
	for line in open(file):
		line_number += 1
		if line_number == 1: #and not header_written:
			#out.write(line)
			#header_written = True
			continue
		line=line.strip()
		seps = line.count(sep)
		if seps < seps_required:
			line += sep*(seps_required-seps)
		out.write("%s\t%s\n" % (chr,"\t".join(line.split(sep))))
	