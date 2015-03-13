#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import argparse
from ConfigParser import SafeConfigParser
try:
    import sh
except ImportError:
    raise ImportError,"The Python sh module is required to run this program. Try 'pip install sh'."

#these files must be in same directory as this script
import my

#requires bgzip and tabix

def run(dbnsfp_zip_url, dir):
    dbnsfp_file = os.path.basename(dbnsfp_zip_url)
    dbnsfp_zip_file = os.path.join(dir, dbnsfp_file)
    dbnsfp_gz = os.path.join(dir, 'dbNSFP.gz')
    done_file = my.done_file(dbnsfp_gz)
    if my.file_exists(done_file):
        my.print_log('dbSNP file already downloaded and installed. To reinstall "rm {}"'.format(done_file))
    else:
        download_done_file = my.done_file(dbnsfp_zip_file)
        if my.file_exists(download_done_file):
            my.print_log('dbSNP file already downloaded. To redownload "rm {}"'.format(download_done_file))
        else:
            my.print_log('Downloading {}.'.format(dbnsfp_zip_url))
            my.getzip(dbnsfp_zip_url, dbnsfp_zip_file, dir)
            my.print_log('Downloaded {}.'.format(dbnsfp_zip_url))
            with open(download_done_file, 'w'):
                pass
        my.print_log('Merging and indexing dbNSFP as {}. This will take ~4 hrs'.format(dbnsfp_gz))
        dbnsfp_glob = os.path.join(dir, 'dbNSFP*_variant.chr*')
        #cat dbNSFP*_variant.chr* | bgzip -c > dbNSFP.gz
        sh.bgzip(sh.cat(sh.glob(dbnsfp_glob), _piped=True), '-c', _out=dbnsfp_gz)
        #tabix -s 1 -b 2 -e 2 dbnsfp_gz
        sh.tabix('-s', '1', '-b', '2', '-e', '2', dbnsfp_gz)
        my.print_log('Merged and indexed dbNSFP as {}.'.format(dbnsfp_gz))
        
        with open(done_file, 'w'):
            pass
        my.print_log('Downloaded and extracted dbNSFP files.')



def main():

    #command line arguments
    parser = argparse.ArgumentParser(
        description = 'download, extract, and index dbSNFP',
        epilog = """vax.vax_dbNSFP_downloader 1.0β1 ©2011-2014 Michael Yourshaw all rights reserved.""")
    parser.add_argument('--dbnsfp_zip_url', '-u', type=str, metavar='URL', required=True,
        help='dbSNFP remote file url')
    parser.add_argument('--dir', '-d', type=str, metavar='DIRECTORY', required=True,
        help='dbSNFP local directory')
    #parse args
    args = parser.parse_args()

    run(dbnsfp_zip_url=args.dbnsfp_zip_url, dir=args.dir,)

if __name__ == "__main__": sys.exit(main())

