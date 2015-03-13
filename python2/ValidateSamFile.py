#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import my
import picard

class ValidateSamFileError(Exception): pass

#--email myourshaw@ucla.edu --debug -i /scratch0/tmp/myourshaw/gmd/reads/library_bams/GMD5A1.GMD5A1_XT_4.gatk.realigned.recalibrated.bam /scratch0/tmp/myourshaw/gmd/reads/library_bams/GMD5A1.GMD5A1_XT_4.markdup.bam /scratch0/tmp/myourshaw/gmd/reads/library_bams/GMD5A2.GMD5A2_XT_5-2.gatk.realigned.recalibrated.bam /scratch0/tmp/myourshaw/gmd/reads/library_bams/GMD5A2.GMD5A2_XT_5-2.markdup.bam /scratch0/tmp/myourshaw/gmd/reads/readgroup_bams/GMD5A1.GMD5A1_XT_4.110511_SN430_0243_B817FLABXX.6.TGACCAA.novoalign.bam /scratch0/tmp/myourshaw/gmd/reads/readgroup_bams/GMD5A1.GMD5A1_XT_4.110623_SN860_0067_2011-100R_A81MVKABXX.4.TGACCAA.novoalign.bam /scratch0/tmp/myourshaw/gmd/reads/readgroup_bams/GMD5A2.GMD5A2_XT_5-2.110511_SN430_0243_B817FLABXX.6.ACAGTGA.novoalign.bam /scratch0/tmp/myourshaw/gmd/reads/readgroup_bams/GMD5A2.GMD5A2_XT_5-2.110623_SN860_0067_2011-100R_A81MVKABXX.5.ACAGTGA.novoalign.bam

def main():
	parser = argparse.ArgumentParser(
		description = 'wrapper for picard ValidateSamFile',
		epilog = 'pypeline.ValidateSamFile version 1.0β1 ©2011-2012 Michael Yourshaw all rights reserved')
	parser.add_argument('--input', '-i', nargs='+',
											help='list of sam/bam files to validate; results -> <input file>.validate',)
	parser.add_argument('--email', nargs='*',
											help='email address(*) for job completion notices',)
	parser.add_argument('--debug', action='store_true', default=False,
											help='do not delete temporary files')
	args = parser.parse_args()
	
	for i in args.input:
		picard.run('ValidateSamFile', email=args.email, debug=args.debug, synchronous=False, INPUT=i, OUTPUT=i+'.validate')

if __name__ == "__main__": sys.exit(main())
