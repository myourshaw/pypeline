#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import argparse
from ConfigParser import SafeConfigParser #configparser in python 3
import shutil
import my
import pysam

class BamPairedEndsError(Exception): pass

#-i /scratch1/tmp/myourshaw/mm_jen_20130505/bams/readgroup_bams/HWI-ST973.93.5.CAGATC.novoalign.fixmate.bam -o /scratch1/tmp/myourshaw/mm_jen_20130505/bams/readgroup_bams/chromosomes_not_same.txt

name = 'pypeline.bam_paired_ends'
version = 0.1
copyright = '©2011-2013 Michael Yourshaw all rights reserved'
run_time = my.localtime_stamp()

def main():

    #command line arguments
    parser = argparse.ArgumentParser(parents=[my.default_parser()],
        description = 'bam paired ends',
        epilog = '{} version {} {}'.format(name, version, copyright))
    parser.add_argument('--input', '-i', nargs='+', required=True,
        help='list of bam files')
    parser.add_argument('--output', '-o', required=True,
        help='list of bam files')

    args = parser.parse_args()

    bams = my.unglob(args.input)
    
    my.makedir(os.path.dirname(args.output))
    with open(args.output,'w') as o:
        o.write('#file\tread1Record\tqname\treadChrom\tread2Chrom\n')
        #https://groups.google.com/forum/?fromgroups=#!topic/pysam-user-group/zvj-REqZodg
        for bam in bams:
            print bam
        #returns dict of header; h['RG'] will be one or more dicts of RG keys and values
            with my.open_bam(bam) as sam:
                header = sam.header
                counter = 0
                record_number = 0
                for read in sam.fetch():
                    record_number += 1
                    if read.is_read1 and not read.is_unmapped and read.is_paired and read.is_proper_pair and not read.mate_is_unmapped:
                        pointer = sam.tell() # pointer to the current position in the BAM file
                        counter += 1
                        if (counter % 100000) == 0:
                            print 'Processed {} pairs'.format(counter)
                        try: 
                            mate = sam.mate(read) 
                        except ValueError: 
                            # Invalid mate (usually post-filtered) 
                            continue 
                        finally: 
                            sam.seek(pointer) # Return the BAM file to the position of read1 in the pair
                        if read.tid != mate.tid:
                            o.write('{}\t{}\t{}\t{}\t{}\n'.format(bam, record_number, read.qname, sam.getrname(read.tid), sam.getrname(mate.tid)))
                print '{} has {} paired reads'.format(bam, counter)



if __name__ == "__main__": sys.exit(main())
