#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import argparse
from collections import Counter
try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET
import my

#-i /scratch1/tmp/myourshaw/resources/metabolome/hmdb_3.6_proteins/HMDBP*.xml -o /scratch1/tmp/myourshaw/resources/metabolome/hmdb_3.6_proteins.txt.foo

def run (input, output):
    files = my.unglob(input)
    my.makedir(os.path.dirname(output))
    genes = Counter()
    with open (output, 'w') as out:
        out.write('#gene\tfunction\tmetabolites\n')
        for f in files:
            tree = ET.ElementTree(file=f)
            root = tree.getroot()
            #for child_of_root in root:
            #    print child_of_root.tag, child_of_root.attrib
            if tree.find('gene_name') == None and tree.find('protein_associations/protein/gene_name') == None:
                continue
            elif tree.find('gene_name') != None:
                gene = '' if not tree.find('gene_name').text else tree.find('gene_name').text.strip()
            else:
                gene = '' if not tree.find('protein_associations/protein/gene_name').text else tree.find('protein_associations/protein/gene_name').text.strip()
            if not gene or genes[gene.upper()]:
                continue
            genes[gene.upper()] += 1
            metabolites = [m.text.strip() for m in tree.iterfind('metabolite_associations/metabolite/name')]
            if not metabolites and  tree.find('metabolite/name') != None and tree.find('metabolite/name').text:
                metabolites = [tree.find('metabolite/name').text.strip()]
            function = None
            #metabolites has description
            if tree.find('description') != None:
                function = '' if not tree.find('description').text else tree.find('description').text.strip().encode('utf-8')
            #proteins has general/specific function
            if not function and tree.find('general_function') != None:
                function = '' if not tree.find('general_function').text else my.l_strip(tree.find('general_function').text.strip(), 'Involved in ').capitalize()
            if function and metabolites:
                function += '; metabolites: {}'.format('|'.join(metabolites))
            if not function and metabolites:
                function = '|'.join(metabolites)
            #specific function is the same as UniProt function
            if not function and tree.find('specific_function') != None:
                function = '' if not tree.find('specific_function').text else tree.find('specific_function').text.strip()
            if not function:
                function = gene
            #print '{}: {}'.format(sub_pathway, gene)
            line_out = '{}\t{}\t{}'.format(
                gene.replace('\t', ' '), function.replace('\t', ' '), '|'.join(metabolites).replace('\t', ' ')).replace('\n', ' ')
            out.write(line_out+'\n')


def main():

    #command line arguments
    parser = argparse.ArgumentParser(
        description = 'generate CREATE TABLE and import statements',
        epilog = 'pypeline.column_widths version 1.0β1 ©2011-2014 Michael Yourshaw all rights reserved')
    parser.add_argument('--input', '-i', required=True, nargs='+',
        help='metabolome xml files')
    parser.add_argument('--output', '-o', required=True,
        help='output text file')
    args = parser.parse_args()

    run(input=args.input, output=args.output)

if __name__ == "__main__": sys.exit(main())


