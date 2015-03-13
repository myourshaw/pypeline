#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import pysam

bam = '/scratch1/tmp/myourshaw/gmdjen/bams/sample_bams/GMD154A.sample.realigned.recalibrated.bam'
samfile = pysam.Samfile( bam, "rb" )
print
