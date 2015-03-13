#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import argparse
from ConfigParser import SafeConfigParser
from contextlib import closing
import copy
from getpass import getpass
import gzip
import re
import shutil
import subprocess
import tarfile
from warnings import warn
from zipfile import ZipFile
try:
    import MySQLdb
    import MySQLdb.cursors
except ImportError:
    raise ImportError,"The Python MySQLdb module is required to run this program. Try 'pip install MySQL-python'."
try:
    from requests import get
except ImportError:
    raise ImportError,"The Python requests module is required to run this program. Try 'pip install requests'."
import my

def getzip(self, url, zipfile, unzipdir): 
    """Download zipfile from url, extract contents to unzipdir, and remove zipfile."""
    done_file = os.path.join(unzipdir, '.'+os.path.basename(zipfile)+'.done')
    if my.file_exists(done_file):
        self.print_log('{} already extracted; skipping. To reinstall "rm {}"'.format(os.path.basename(zipfile), done_file))
    else:
        self.print_log('Downloading  {} as .'.format(url, zipfile))
        with closing(get(url, stream=True)) as r:
            with open(zipfile, 'wb') as fd:
                for chunk in r.iter_content():
                    fd.write(chunk)
        self.print_log('Extracting into {}.'.format(unzipdir))
        with ZipFile(zipfile, 'r') as zip:
            zip.extractall(unzipdir)
        os.remove(zipfile)
        with open(done_file, 'w'):
            pass

url = 'http://dbnsfp.houstonbioinformatics.org/dbNSFPzip/'
with closing(get(url, stream=True)) as r:
