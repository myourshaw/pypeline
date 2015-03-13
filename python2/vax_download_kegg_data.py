#!/usr/bin/env python
# -*- coding: utf-8 -*-
__author__ = 'myourshaw'
import os
import sys
import contextlib
import argparse
import urllib2
from warnings import warn
import my
try:
    from requests import get
except ImportError:
    raise ImportError,"The Python requests module is required to run this program. Try 'pip install requests'."

#-o /scratch1/tmp/myourshaw/rodent/kegg_rodent.txt -s mmu rno


@contextlib.contextmanager
def file_writer(file_name = None):
    # Create writer object based on file_name
    writer = open(file_name, 'w') if file_name else sys.stdout
    # yield the writer object for the actual use
    yield writer
    # If it is file, then close the writer object
    if file_name:
        writer.close()


def get_kegg(operation, argument, argument2):
    return get('http://rest.kegg.jp/{}/{}/{}'.format(operation, argument, argument2)).content


def run(output=None, species=['hsa', 'mmu']):
    if output:
        my.makedir(os.path.dirname(output))
    with file_writer(output) as o:
        o.write('#gene_id\tpath_id\tpathway\n')
        for s in species:
            print('Downloading species {} pathways and genes.'.format(s))
            pathways = get_kegg('list', 'pathway', s).strip()
            if not pathways:
                warn('Pathways not found for species {}'.format(s))
            else:
                for p in pathways.split('\n'):
                    if p:
                        path_id, pathway = p.split('\t')
                        path_id = path_id[len('path:'):]
                        #pathway = pathway.split(' - {}'.format(s[0].upper()))[0]
                        genes = get_kegg('link', 'genes', path_id).strip().split('\n')
                        for gene in genes:
                            if gene:
                                path_id, gene_id = gene.split('\t')
                                o.write('\t'.join((gene_id, path_id, pathway))+'\n')
    if output:
        with open(my.done_file(output), 'w'):
            pass

    print('vax_download_kegg_data finished')


def main():
    parser = argparse.ArgumentParser(
        description = 'Download KEGG genes and pathways',
        epilog = 'vax.download_kegg_data 1.0β1 ©2014 Michael Yourshaw all rights reserved.')
    parser.add_argument('--output', '-o', type=str, metavar='FILE',
        help='output file; (default: STDOUT)')
    parser.add_argument('--species', '-s', nargs='*', default=['hsa', 'mmu'], metavar='SPECIES',
        help='List of species for which to get genes (default: hsa mmu')

    #parse args
    args = parser.parse_args()

    run(output=args.output, species=args.species)

if __name__ == "__main__": sys.exit(main())
