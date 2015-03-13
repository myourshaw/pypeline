#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import argparse
import my
import convert_contigs

class ConvertContigsException(Exception): pass

#-i /scratch1/vax/75/refgene/refGene_ucsc_hg19.txt -o /scratch1/vax/75/refgene/refGene_ucsc_b37.txt --chrom_col 2 --pos_col 4 --faidx /share/apps/myourshaw/vax/resources/faidx/faidx_b37_hg19_decoy.txt --input_genome_id hg19 --output_genome_id b37

def run(input, faidx, input_genome_id, output_genome_id, output, chrom_col=0, pos_col = 1, header_line_count=1, header_lines_startwith=['#','@']):
    
    done_file = my.done_file(output)
    
    genomes = {}
    converted = {}
    
    with open(faidx) as fx:
        faidx_columns = fx.readline().lstrip('#').rstrip('\n').split('\t')
        try:
            contig_col = faidx_columns.index('contig')
        except:
            raise ConvertContigsException('faidx file must have contig, sortOrder, genome, genomeAlias, contigAlias, and sortOrderAlias columns.')
        try:
            sortOrder_col = faidx_columns.index('sortOrder')
        except:
            raise ConvertContigsException('faidx file must have contig, sortOrder, genome, genomeAlias, contigAlias, and sortOrderAlias columns.')
        try:
            genome_col = faidx_columns.index('genome')
        except:
            raise ConvertContigsException('faidx file must have contig, sortOrder, genome, genomeAlias, contigAlias, and sortOrderAlias columns.')
        try:
            genomeAlias_col = faidx_columns.index('genomeAlias')
        except:
            raise ConvertContigsException('faidx file must have contig, sortOrder, genome, genomeAlias, contigAlias, and sortOrderAlias columns.')
        try:
            contigAlias_col = faidx_columns.index('contigAlias')
        except:
            raise ConvertContigsException('faidx file must have contig, sortOrder, genome, genomeAlias, contigAlias, and sortOrderAlias columns.')
        try:
            sortOrderAlias_col = faidx_columns.index('sortOrderAlias')
        except:
            raise ConvertContigsException('faidx file must have contig, sortOrder, genome, genomeAlias, contigAlias, and sortOrderAlias columns.')
        for line in fx:
            if line == '\n' or line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            genomes[(fields[genome_col], fields[contig_col])] = (fields[sortOrder_col], fields[genomeAlias_col], fields[contigAlias_col], fields[sortOrderAlias_col], )
            
    
    with open(input) as inputh, open(output, 'w') as outputh:
        outputs = []
        for line in inputh:
            if header_line_count>0 or True in [line.startswith(s) for s in header_lines_startwith]:
                outputh.write(line)
                header_line_count -= 1
                continue
            line = line.rstrip('\n')
            data = line.split('\t')
            chrom = data[chrom_col]
            pos = data[pos_col]
            #new chrom
            data[chrom_col] = genomes[(input_genome_id, chrom)][2]
            #new sort order in data[0]
            data.insert(0, int(genomes[(input_genome_id, chrom)][3]))
            outputs.append(data)
        #sort by new sort_order, pos
        outputs.sort(key=lambda data: (data[0],int(data[pos_col+1])) )
        for o in outputs:
            output_line = '\t'.join(o[1:])+'\n'
            outputh.write('\t'.join(o[1:])+'\n')
    with open(done_file, 'w'):
        pass

def main():

    #command line arguments
    parser = argparse.ArgumentParser(
        description = 'convert the chromosome field to a different genome build and sort',
        epilog = 'convert_contigs 1.0β1 ©2011-2014 Michael Yourshaw all rights reserved')

    parser = argparse.ArgumentParser(
        description = 'create a GRCh-sytle file and an exome cds + essential splice site interval_list from ucsc refGene table',
        epilog = 'pypeline.refGene2interval_list version 1.0β1 ©2011-2014 Michael Yourshaw all rights reserved')
    parser.add_argument('--input', '-i', required=True,
        help='input UCSC refGene file (default: file will be downloaded from UCSC mysql server',)
    parser.add_argument('--faidx', required=True,
        help='faidx file relating hg to GRCh contigs, with GATK sort orders',)
    parser.add_argument('--input_genome_id', required=True,
        help='faidx genome of input file (example: hg19',)
    parser.add_argument('--output_genome_id', required=True,
        help='faidx genome of output file (example: b37',)
    parser.add_argument('--output', '-o', required=True,
        help='output refGene file')
    parser.add_argument('--chrom_col', default=0, type=int,
        help='zero-based column number of chromosome column to be converted (default: 0)')
    parser.add_argument('--pos_col', default=1, type=int,
        help='zero-based column number of position column for sorting (default: 1)')
    parser.add_argument('--header_line_count', default=1, type=int,
        help='number of header lines in input file to be copied to output file (default: 1)')
    parser.add_argument('--header_lines_startwith', default=['#','@'], nargs='*',
        help='list of strings that header lines start with(default: # @)')
    
    args = parser.parse_args()
    
    run(input=args.input, faidx=args.faidx, input_genome_id=args.input_genome_id,
        output_genome_id= args.output_genome_id, output=args.output,
        chrom_col=args.chrom_col, pos_col=args.pos_col, header_line_count=args.header_line_count, header_lines_startwith=args.header_lines_startwith)

if __name__ == "__main__": sys.exit(main())
