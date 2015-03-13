#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
import logging
try:
    from gtf_to_genes import *
except:
    raise ImportError("The gtf_to_genes module is required to run this program. Try 'pip install gtf_to_genes'.")


#--index_file /scratch1/tmp/myourshaw/resources/Illumina_iGenomes/gtf.index --search_path_root /scratch1/tmp/myourshaw/resources/Illumina_iGenomes --regex_input (.+Illumina_iGenomes/([^/]+)/(Ensembl)/([^/]+)/Annotation/Genes)/genes\.gtf$ --cache_file_pattern \1/\2.\3.\4.cache --identifier_pattern 2:\3:\4

def run(index_file, search_path_root, regex_input, cache_file_pattern, identifier_pattern):
    ignore_cache = False
    logger = logging.getLogger("test")

#/scratch1/tmp/myourshaw/resources/Illumina_iGenomes/Mus_musculus/Ensembl/GRCm38/Annotation/Genes/genes.gtf
#/scratch1/tmp/myourshaw/resources/Illumina_iGenomes/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf
#Illumina's NCBI gtf's fail for lack of a gene_biotype attribues
#\1 path root
#\2 organism
#\3 database
#\4 genome build
#or  r"{INDEX_FILE_PATH}/\2.cache"

    # test(index_file, logger, 'Mus_musculus:Ensembl:GRCm38')
    # return
    index_gtf_files(index_file,
                    search_path_root,
                    regex_input,
                    cache_file_pattern,
                    identifier_pattern,
                    ignore_cache,
                    logger)

    print 'Index done. Testing.'


def test(index_file, logger, identifier):
    logger = logging.getLogger("test")
    species, gtf_file_name, genes = get_indexed_genes_for_identifier(index_file, logger, identifier)
    print species
    if genes:
        print genes.keys()
        print "# of protein coding genes = ", len(genes['protein_coding'])


def main():

    #command line arguments
    parser = argparse.ArgumentParser(
        description = 'Create index and cache for gtf_to_gene library. See http://gtf-to-genes.googlecode.com/hg/doc/build/html/index.html',
        epilog = 'gtf_to_gene_indexer 1.0β1 ©2014 Michael Yourshaw all rights reserved')
    parser.add_argument('--index_file', '-x', type=str, required=True,
        help='Index file to hold list of parsed GTF files')
    parser.add_argument('--search_path_root', '-s', type=str, required=True,
        help='Root directory to start scanning for GTF files')
    parser.add_argument('--regex_input', '-r', type=str, required=True,
        help='Regular expression to match GTF or gzipped GTF files. Brackets can be used to construct the corresponding cache file name. E.g. r"(.+\/)(([^.]+)\..+\.(.+)\.gtf(?:\.gz)?)$"')
    parser.add_argument('--cache_file_pattern', '-c', type=str, required=True,
        help='Pattern used to construct cache file name after regular expression substitution. Brackets can be used to construct the corresponding cache file name E.g. r"\1\2.cache" or r"{INDEX_FILE_PATH}/\2.cache" might give “/path/to/gtf/homo_sapiens.cache”')
    parser.add_argument('--identifier_pattern', '-i', type=str, required=True,
        help='Pattern used to constuct the (e.g. species) identifier for this GTF file E.g. r"\1:\2" might give "homo_sapiens:47"')

    #parse args
    args = parser.parse_args()


    run(index_file=args.index_file, search_path_root=args.search_path_root, regex_input=args.regex_input, cache_file_pattern=args.cache_file_pattern, identifier_pattern=args.identifier_pattern)

if __name__ == "__main__": sys.exit(main())
