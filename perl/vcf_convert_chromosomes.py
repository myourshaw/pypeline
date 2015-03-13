#!/usr/bin/env python

import sys
import os
import re
import traceback
from optparse import OptionParser, OptionGroup
from IndentedHelpFormatterWithNL import *
import gzip

description = """
This script merges, converts contig names, and sorts vcf files.
"""
parser = OptionParser( description=description, usage="usage: %prog [options] INPUT-FILE", formatter=IndentedHelpFormatterWithNL())
parser.add_option("-g", "--genome-translate-file", help="Genome translate file path used for filtering chromosomes (e.g., chr6_apd_hap1), converting chromosomes (e.g., chr1 -> 1) and determining chromosome sort order [#genome_in\tchromosome_in\tsort_order_in\tgenome_out\tchromosome_out(case-sensitive)\tsort_order_out")
parser.add_option("-y", "--genome-in", help="genome id for chromosomes in input file; must be in column 0 of genome translate file")
parser.add_option("-z", "--genome-out", help="genome id for chromosomes in output file; must be in column 3 of genome translate file and refers to the reference file that will be used when running the GenomicAnnotator")
parser.add_option("-o", "--output-filename", help="Output file path [Default: %default]", default="stdout")
parser.add_option("-v", "--verbose", action="store_true", default=False, help="Verbose.")
parser.add_option("-d", "--delimiter", help="The delimiter that separates values in a line of INPUT-FILE. Set to 'tab' to make it use tab [Default: whitespace].")

(options, args) = parser.parse_args()

verbose = options.verbose

delimiter = options.delimiter
if delimiter and delimiter.lower() == "tab":
    delimiter = "\t"

if len(args) < 1 or not os.access(args[0], os.R_OK):
    error("Requires a valid INPUT-FILE")
input_filename = args[0]

output_filename = options.output_filename
if output_filename and output_filename != "stdout" and output_filename != "-" and os.access(output_filename, os.F_OK) and not os.access(output_filename, os.W_OK):
    error("Unable to write to: %s" % str(options.output_filename))

OUTPUT_FORMAT_DELIMITER = "\t"

header_lines = []
data_lines = []
chrom_col = 0
pos_col = 1
counter = 0
skipped_lines_counter = 0
previous_n = -1 # Checks whether data is in order
need_to_sort = False

def error(msg):
    print("ERROR: %s.        (Rerun with -h to print help info) \n" % msg)
    #parser.print_help()
    sys.exit(-1)

def warn(msg):
    print("WARNING: %s" % msg)

def fatal(msg):
    print(msg)
    sys.exit(-1)


def join_fields(fields):
    return OUTPUT_FORMAT_DELIMITER.join(fields)

def split_line(line):
    if delimiter:
        return line.split(delimiter)
    else:
        return line.split()

def line_key(line):
	line_split = split_line(line)
	return chrpos_to_n(line_split[0], line_split[1])

#creates a table to translate input chromosome to output chromosome and sort order
#for the specified combination of genome-in and genome-out options
#input format is: genome_in\tchromosome_in\tsort_order_in\tgenome_out\tchromosome_out\tsort_order_out
#chromosome_out is case-sensitive and must match alignment in bam files
# ignores comment lines starting with #
# genome names are arbitrary
# translations from a genome to itself should be provided
# chromosome_translate_table is a dictionary with key = chromosome_in and value = chromosome_out
# chromosome_sort_order_table is a dictionary with key = chromosome_out and value = sort_order_out)
def load_genome_translate_tables(genome_translate_file, genome_in, genome_out, chromosome_translate_table,chromosome_sort_order_table):
    genome_in = genome_in.lower().strip()
    genome_out = genome_out.lower().strip()
    for line in open(genome_translate_file):
        line = line.strip()
        if line.startswith("#") or line == '':
            continue
        line_fields = line.split("\t")
        if len(line_fields) < 6:
            error("Found only %d fields in genome translate table." % len(line_fields))
        genome_in_line = line_fields[0].lower().strip()
        genome_out_line = line_fields[3].lower().strip()
        if genome_in_line == genome_in and genome_out_line == genome_out:
            chromosome_in = line_fields[1].lower().strip()
            chromosome_out = line_fields[4].strip()
            try: sort_order_out = int(line_fields[5].strip())
            except: error("Non-integer sort_order %s in genome-translate-file: [%s]" % (line_fields[5].strip(), line))
            chromosome_translate_table[chromosome_in] = chromosome_out
            chromosome_sort_order_table[chromosome_out] = sort_order_out

# Computes an integer key for this line. These keys can be used to sort the lines by reference
def chrpos_to_n(chr,pos):
	chr_n = int(chromosome_sort_order_table[chr.strip()]) #chromosome name is case-sensitive
	start_n = long(pos.strip().lower())
	end_n = 0
	N = (chr_n * 10L**23) + (start_n * 10L**11) + end_n # Combine chr, start, stop into a single numeric key for sorting
	return N

use_genome_translate_file = bool(options.genome_translate_file and options.genome_in and options.genome_out)
if use_genome_translate_file:
    chromosome_translate_table = {}
    chromosome_sort_order_table = {}
    load_genome_translate_tables(options.genome_translate_file,options.genome_in,options.genome_out,chromosome_translate_table,chromosome_sort_order_table)
else:
	error("Options genome-translate-file, genome-in, and genome-out are required.")
	
# begin processing
GZIP_MAGIC = b"\x1F\x8B"
input_file = open(input_filename, "rb")
magic = input_file.read(len(GZIP_MAGIC))
input_file.close()
if magic == GZIP_MAGIC:
	input_file = gzip.open(input_filename, "rb")
else:
		input_file = open(input_filename, "r")
for line in input_file:
		counter+=1
		line = line.strip()
		if line.startswith("#") or line == "":
				#this is a header line
				if line.startswith("##reference="): ##reference=1000genomes_NCBI36
						line = "##reference=" + options.genome_out
				header_lines += [line]
		else:
				# This is a data line
				line_fields = split_line(line)
				if (line_fields[chrom_col].lower()) not in chromosome_translate_table: #chromosome name input is not case sensitive
					if verbose:
						print "Chromosome %s not in genome_translate_table: [%s]. Skipping..." % (line_fields[chrom_col].lower(),join(split_line(line)))
					skipped_lines_counter += 1
					continue           
				line_fields[chrom_col] = chromosome_translate_table[line_fields[chrom_col].lower()] #chromosome name output is case sensitive
				try:
					n = chrpos_to_n(line_fields[chrom_col],line_fields[pos_col])
					if not need_to_sort and n < previous_n:
						need_to_sort = True
						warn("Line %d is out of order. Will need to sort all lines." % counter)
					previous_n = n
				except Exception, e:
					warn("Couldn't parse line: " + "  ".join(line_fields) + ". " +str(e) + ". Skipping...")
					if verbose: traceback.print_exc()
					skipped_lines_counter += 1
					continue
				data_lines += [ join_fields(line_fields) ]


if verbose and skipped_lines_counter:
    print("Skipped %d / %d lines. (%f%%)" % (skipped_lines_counter, counter, skipped_lines_counter/float(counter)))
if need_to_sort:
  if verbose:
    print("Sorting %d lines..." % len(data_lines))
  data_lines.sort(key=line_key)

if verbose:
  print("Writing data to: " + output_filename)

# Write output file
if output_filename == "stdout" or output_filename == "-":
  output_file = sys.stdout
elif output_filename.endswith("gz"):
  output_file = gzip.open(output_filename, "wb")
else:
  output_file = open(output_filename, "w")

for line in header_lines:
	output_file.write(line + "\n")

for line in data_lines:
	output_file.write(line + "\n")

output_file.close()
