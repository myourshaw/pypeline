#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import argparse
from contextlib import closing
import re
from urllib import urlretrieve
import my
try:
    from bs4 import BeautifulSoup
except ImportError:
    raise ImportError,"The Python BeautifulSoup module is required to run this program. Try 'pip install beautifulsoup4'."

#-u http://www.gardnermuseum.org/music/listen/music_library -d /Users/myourshaw/Downloads/gardnermuseum_music -p '.+\.mp3'

def list_hrefs(url, pattern=None, case_insensitive=True):
    """Scrape a web page for hrefs, filter by regex pattern"""
    if pattern:
        rx = re.compile(pattern, re.I) if case_insensitive else re.compile(pattern)
    with closing(requests.get(url)) as r:
        html = r.text
    soup = BeautifulSoup(html)
    hrefs = []
    for link in soup.find_all('a'):
        href =  link.get('href')
        if href and (pattern == None or rx.match(href)):
            hrefs.append(href)
    return hrefs

def run(url, download_directory, pattern=None, case_insensitive=True):
    hrefs = list(set(list_hrefs(url, pattern, case_insensitive)))
    my.makedir(download_directory)
    for href in hrefs:
        print href
        try:
            urlretrieve(href, os.path.join(download_directory, os.path.basename(href)))
        except IOError, e:
            print 'Problem downloading {}. {}'.format(href, e)

def main():

    #command line arguments
    parser = argparse.ArgumentParser(
        description = 'download all or selected hrefs from a url',
        epilog = 'pypeline.download_url version 1.0β1 ©2011-2014 Michael Yourshaw all rights reserved')
    parser.add_argument('--url', '-u', required=True,
        help='url from which to download')
    parser.add_argument('--download_directory', '-d', required=True,
        help='directory for downloaded files')
    parser.add_argument('--pattern', '-p', default=None,
        help='regular expression filter to select files to download (default: None)')
    parser.add_argument('--case_insensitive', '-i', action='store_false', default=True,
        help='whether pattern is case insensitive (default: True)')
    args = parser.parse_args()
    
    run(url=args.url, download_directory=args.download_directory, pattern=args.pattern, case_insensitive=args.case_insensitive)


if __name__ == "__main__": sys.exit(main())
