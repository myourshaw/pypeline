#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import argparse
import re
import my

class RefGene2DbError(): pass

#http://biostar.stackexchange.com/questions/2145/bulk-download-of-ncbi-gene-summary-field

def run(ref=None, input=None, output=None, dir=None):

    if not dir and not (my.file_exists(ref) and my.file_exists(input)):
        raise RefGene2DbError('No inputs.')
    if dir and not my.file_exists(ref):
        ref = os.path.join(dir, 'gene_RefSeqGene')
    if dir and not input:
        input = os.path.join(dir, 'refseqgene*.genomic.gbff.gz')
    if not output:
            output = os.path.join(dir if dir else os.path,dirname(ref), 'refseq_summary.txt')
            
    locus2gene = {}
    with open(ref) as ref:
        for r in ref:
            if r.startswith('#'):
                continue
            r = r.rstrip('\n').split('\t')
            locus2gene[r[3].split('.')[0]] = r[2]
    
    inputs = my.unglob(input)
    locus2comment = {}
    locus_count = 0
    refseq_re = re.compile(r"\s*\[provided by RefSeq.*\]",re.I)
    for input in inputs:
        f = my.open_gz_or_text(input)
        in_comment=False
        for line in f:
            if line[0:5] == "LOCUS":
                locus = line.split()[1]
                comment = ""
                locus_count += 1
            elif line[0:7] == "COMMENT":
                in_comment=True
                comment += line.split("    ")[1].replace("\n", " ")
            elif line[0:7] == "PRIMARY":
                in_comment = False
                try:
                    summary = comment.split("Summary:")[1]#.strip().split('[provided by RefSeq].')[0].rstrip()
                except:
                    summary = comment#.strip().split('[provided by RefSeq].')[0].rstrip()
                locus2comment[locus] = refseq_re.split(summary)[0]
            elif in_comment:
                comment += line.split("            ")[1].replace("\n", " ")
    with open(output,'w') as output:
        output.write('#NG_ID\tGeneSymbol\tSUMMARY\n')
        for locus in sorted(locus2comment):
            output.write('{}\t{}\t{}\n'.format(locus, locus2gene.get(locus,''), locus2comment[locus]))


def main():

    #command line arguments
    parser = argparse.ArgumentParser(
        description = 'Convert RefGene Genbank protein files to tab-delimited database files',
        epilog = 'pypeline.refgene2db version 1.0β1 ©2011-2014 Michael Yourshaw all rights reserved')
    parser.add_argument('--dir', '-d',
        help='directory that contains gene_RefSeqGene and refseqgene*.genomic.gbff.gz files')
    parser.add_argument('--ref', '-r',
        help='path to gene_RefSeqGene file (overrides ')
    parser.add_argument('--input', '-i', nargs='*',
        help='downloaded refseqgene*.genomic.gbff.gz files')
    parser.add_argument('--output', '-o',
        help='output file')
    args = parser.parse_args()
    
    run(ref=args.ref, input=args.input, output=args.output, dir=args.dir)

if __name__ == "__main__": sys.exit(main())
