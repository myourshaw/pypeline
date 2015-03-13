#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import argparse
from ConfigParser import SafeConfigParser
from urllib import urlretrieve

#these files must be in same directory as this script
import my

def run(snp_file_url, dir):
    done_file = my.done_file(os.path.join(dir,'cadd_download'))
    if my.file_exists(done_file):
        my.print_log('CADD files already downloaded. To reinstall "rm {}"'.format(done_file))
    else:
        my.print_log('Installing CADD.')
        my.print_log('Downloading CADD files.')
        #download selected version
        snp_file = os.path.join(dir, os.path.basename(snp_file_url))
        done_file = my.done_file(snp_file)
        tbi_file_url = snp_file_url+'.tbi'
        tbi_file = snp_file+'.tbi'
        my.print_log('Downloading CADD file {} to {}. This will take >2 hours.'.format(snp_file_url, snp_file))
        urlretrieve(snp_file_url, snp_file)
        my.print_log('Downloaded CADD file {} to {}. This will take >2 hours.'.format(snp_file_url, snp_file))
        my.print_log('Downloading CADD.tbi file {} to {}. This will take >2 hours.'.format(tbi_file_url, tbi_file))
        urlretrieve(tbi_file_url, tbi_file)
        my.print_log('Downloaded CADD.tbi file {} to {}. This will take >2 hours.'.format(tbi_file_url, tbi_file))
        
        with open(done_file, 'w'):
            pass
        my.print_log('Downloaded CADD files.')



def main():

    #command line arguments
    parser = argparse.ArgumentParser(
        description = 'download CADD',
        epilog = """vax.vax_CADD_downloader 1.0β1 ©2011-2014 Michael Yourshaw all rights reserved.""")
    parser.add_argument('--snp_file_url', '-s', type=str, metavar='SNP_FILE_URL', required=True,
        help='CADD remote SNP file url')
    parser.add_argument('--dir', '-d', type=str, metavar='CADD_LOCAL_DIRECTORY', required=True,
        help='CADD local directory')
    #parse args
    args = parser.parse_args()

    run(snp_file_url=args.snp_file_url, dir=args.dir,)

if __name__ == "__main__": sys.exit(main())

