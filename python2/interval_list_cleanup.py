#!/usr/bin/env python
# -*- coding: utf-8 -*-
__author__ = 'myourshaw'

import sys
import os
import argparse
import my

#-i /share/apps/myourshaw/resources/intervals/RefGeneIntervals/RefGene.cds.interval_list


def run(input, output=None, chroms=[str(c) for c in range(1,23)]+['X','Y','MT'],):
    """remove non standard contigs from a picard-style intervals_list file"""
    if not output:
        output = my.swap_ext(input, 'interval_list', '1-22XYMT.interval_list')


    my.makedir(os.path.dirname(output))
    with open(input, 'r') as in_f, open(output, 'w') as out_f:
        for line in in_f:
            fields = line.rstrip('\n').upper().split('\t')
            if fields[0].startswith('@HD'):
                out_f.write(line)
            elif fields[0].startswith('@SQ'):
                if fields[1][3:] in chroms:
                    out_f.write(line)
            elif fields[0] in chroms:
                out_f.write(line)

def main():
    """command line arguments"""
    parser = argparse.ArgumentParser(
        description = 'remove non standard contigs from a picard-style intervals_list file',
        epilog = 'interval_list_cleanup 1.0β1 ©2014 Michael Yourshaw all rights reserved')

    parser.add_argument('--input', '-i', type=str, required=True,
        help='input interval_list file',)
    parser.add_argument('--output', '-o', type=str,
        help='output cleaned interval_list file (default: <input>.1-22XYMT.interval_list)')
    parser.add_argument('--chroms', '-c', nargs='*', default=[str(c) for c in range(1,23)]+['X','Y','MT'],
        help='valid contig names (default: {1..22} X Y MT)')

    args = parser.parse_args()

    run(input=args.input, output=args.output, chroms=args.chroms,)


if __name__ == "__main__": sys.exit(main())
