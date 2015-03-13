#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import argparse
import my

#-i /scratch1/tmp/myourshaw/gmd/runs_tmp/110511_SN430_0243_B817FLABXX.6.1.ACAGTGA_qseq.txt.gz /scratch1/tmp/myourshaw/gmd/runs_tmp/110511_SN430_0243_B817FLABXX.6.2.ACAGTGA_qseq.txt.gz /scratch1/tmp/myourshaw/gmd/runs_tmp/110511_SN430_0243_B817FLABXX.6.3.ACAGTGA_qseq.txt.gz -o /scratch1/tmp/myourshaw/gmd/runs_tmp/qseq_validate.txt
#-i /scratch1/tmp/myourshaw/gmd/runs_tmp/tmp_novobarcode_*/novobarcode_lane_6*/ACAGTGA/* -o /scratch1/tmp/myourshaw/gmd/runs_tmp/qseq_validate_tiles_6.ACAGTGA.txt

def main():
	parser = argparse.ArgumentParser(
	description = 'Compare pairs of qseq files by id and calculate pass/fail metrics',
	epilog = 'pypeline.qseq_validat version 1.0Î²1 Â©2011-2012 Michael Yourshaw all rights reserved')
	parser.add_argument('--input', '-i', nargs='+',
											help='list of qseq files')
	parser.add_argument('--output', '-o', required=True,
											help='output report')
	args = parser.parse_args()

	input_files = sorted(set(my.list_files(args.input)))
	with open(args.output,'w') as o:
		for file in input_files:
			print file
			count = 0
			try:
				i = my.open_gz_or_text(file)
			except Exception as e:
				msg = 'open_exception\t{}\t{}\t{}\n'.format(count, file, e)
				print msg
				o.write(msg)
			else:
				while True:
					try:
						line = i.readline()
					except Exception as e:
						msg = 'read_exception\t{}\t{}\t{}\n'.format(count, file, e)
						print msg
						o.write(msg)
					else:
						if not line:
							break
						else:
							count += 1
							if my.qseq_line_valid(line):
								continue
							else:
								msg = 'invalid_line\t{}\t{}\t{}\n'.format(count, file, line)
								print msg
								o.write(msg)


if __name__ == "__main__": sys.exit(main())
