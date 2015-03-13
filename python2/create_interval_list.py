#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import argparse
import my
import interval_list_cleanup

#-i /scratch1/vax/75/refgene/refGene_b37_ucsc.tmp.txt --faidx /share/apps/myourshaw/vax/resources/faidx/faidx_b37_hg19_decoy.txt --interval_list_header /share/apps/myourshaw/vax/resources/faidx/interval_list_b37_decoy_header.txt
#-i /scratch1/tmp/myourshaw/resources/intervals/RefGeneIntervals/current/refGene_b37_ucsc.txt --faidx /share/apps/myourshaw/vax/resources/faidx/faidx_b37_hg19_decoy.txt --interval_list_header /share/apps/myourshaw/vax/resources/faidx/interval_list_b37_decoy_header.txt


class ConvertContigsException(Exception): pass


def run(input, cds_interval_list, exon_interval_list, gene_interval_list, faidx, interval_list_header, genome_id='b37',
        zero_based=True, header_line_count=1, header_lines_startwith=['#','@'],
        chrom_col=2, strand_col=3, name_cols=[12,1], splice_bases=2,
        tx_start_col=4, tx_end_col=5, cds_start_col=6, cds_end_col=7, exon_starts_col=9, exon_ends_col=10,):
    
    print('Starting create_interval_list')
    if not cds_interval_list:
        cds_interval_list=my.swap_ext(input, '.txt', '.cds.all_contigs.interval_list')
    if not exon_interval_list:
        exon_interval_list=my.swap_ext(input, '.txt', '.exon.all_contigs.interval_list')
    if not gene_interval_list:
        gene_interval_list=my.swap_ext(input, '.txt', '.gene.all_contigs.interval_list')
        
    cds_interval_list_done_file = my.done_file(cds_interval_list)
    exon_interval_list_done_file = my.done_file(exon_interval_list)
    gene_interval_list_done_file = my.done_file(gene_interval_list)
    
    genomes = {}
    
    Interval = namedtuple('Interval', 'sort_order, chrom,  strand, name, tx_start, tx_end, cds_start, cds_end, exon_starts, exon_ends')
    
    print('Building genome cross-reference.')
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

    print('Writing header metadata.')
    with open(cds_interval_list, 'w') as cds_interval_listh, open(exon_interval_list, 'w') as exon_interval_listh, open(gene_interval_list, 'w') as gene_interval_listh:
        with open(interval_list_header) as interval_list_headerh:
            for line in interval_list_headerh:
                cds_interval_listh.write(line)
                exon_interval_listh.write(line)
                gene_interval_listh.write(line)
    
        print('Reading {} and building interval lists.'.format(input))
        with open(input) as inputh:
            cds_intervals = []
            exon_intervals = []
            gene_intervals = []
            #create arrays of all intervals
            for line in inputh:
                if header_line_count>0 or True in [line.startswith(s) for s in header_lines_startwith]:
                    header_line_count -= 1
                    continue
                line = line.rstrip('\n')
                data = line.split('\t')
                chrom = data[chrom_col]
                sort_order = int(genomes[(genome_id, chrom)][0])
                #convert to 1-based
                tx_start = int(data[tx_start_col]) + (1 if zero_based else 0)
                tx_end = int(data[tx_end_col])
                cds_start = int(data[cds_start_col]) + (1 if zero_based else 0)
                cds_end = int(data[cds_end_col])
                exon_starts = [int(s) + (1 if zero_based else 0) for s in data[exon_starts_col].rstrip(',').split(',')]
                exon_ends = [int(s) for s in data[exon_ends_col].rstrip(',').split(',')]
                if strand_col == None:
                    strand = '+'
                else:
                    strand = '-' if data[strand_col] in ('-' ,'-1') else '+'
                if name_cols == None:
                    name = '{}:{}-{}'.format(data[chrom_col], data[tx_start_col], data[tx_end_col])
                else:
                    name = '|'.join([data[int(n)] for n in name_cols])
                #save gene interval
                gene_intervals.append(Interval(sort_order=sort_order, chrom=chrom, strand=strand, name=name,
                                          tx_start=tx_start, tx_end=tx_end, cds_start=cds_start, cds_end=cds_end, exon_starts=exon_starts, exon_ends=exon_ends,))
                exons = zip(exon_starts, exon_ends)
                for i,e in enumerate(exons):
                    exon_start, exon_end = e
                    #add splice bases
                    if len(exons) == 1:
                        splice_exon_start, splice_exon_end = (exon_start, exon_end)
                    elif i == 0:
                        splice_exon_start, splice_exon_end = (exon_start, exon_end+splice_bases)
                    elif i == len(exons)-1:
                        splice_exon_start, splice_exon_end = (exon_start-splice_bases, exon_end)
                    else:
                        splice_exon_start, splice_exon_end = (exon_start-splice_bases,exon_end+splice_bases)
                    #save exon (with utr) interval
                    exon_intervals.append(Interval(sort_order=sort_order, chrom=chrom, strand=strand, name=name+'[e{}]'.format(i+1),
                                          tx_start=tx_start, tx_end=tx_end, cds_start=cds_start, cds_end=cds_end, exon_starts=splice_exon_start, exon_ends=splice_exon_end,))
                    #extract translated region
                    # ---uuu---
                    if cds_start > exon_end or cds_end < exon_start:
                        continue
                    # ---uuuCCCuuu---
                    elif cds_start > exon_start and cds_end < exon_end:
                        splice_cds_start, splice_cds_end = cds_start, cds_end
                    # ---uuuCCC---
                    elif cds_start > exon_start and cds_end >= exon_end:
                        splice_cds_start, splice_cds_end = cds_start, splice_exon_end
                    # ---CCCuuu---
                    elif cds_start <= exon_start and cds_end < exon_end:
                        splice_cds_start, splice_cds_end = splice_exon_start, cds_end
                    # ---CCC---
                    else:
                        splice_cds_start, splice_cds_end = splice_exon_start, splice_exon_end
                    #save cds (no utr) interval
                    cds_intervals.append(Interval(sort_order=sort_order, chrom=chrom, strand=strand, name=name+'[e{}]'.format(i+1),
                                              tx_start=tx_start, tx_end=tx_end, cds_start=cds_start, cds_end=cds_end, exon_starts=splice_cds_start, exon_ends=splice_cds_end,))
                        
        #sort by sort_order(chrom), start, end
        cds_intervals.sort( key=lambda data: (data.sort_order, data.exon_starts, data.exon_ends) )
        exon_intervals.sort( key=lambda data: (data.sort_order, data.exon_starts, data.exon_ends) )
        gene_intervals.sort( key=lambda data: (data.sort_order, data.tx_start, data.tx_end) )
        
        #write cds interval interval_list
        print('Writing {}.'.format(cds_interval_list))
        that_chrom = None
        that_start = None
        that_end = None
        that_strand = None
        names = set()
        for o in cds_intervals:
            #chrom, start, end, strand, name = o[1:6]
            if o.chrom == that_chrom and o.strand == that_strand:
                if o.exon_starts > that_end:
                    cds_interval_listh.write('{}\t{}\t{}\t{}\t{}\n'.format(that_chrom, that_start, that_end, that_strand, ','.join(names)))
                    #print('{}\t{}\t{}\t{}\t{}'.format(that_chrom, that_start, that_end, that_strand, ','.join(names)))
                    that_start, that_end = o.exon_starts, o.exon_ends
                    names = {o.name}
                else:
                    that_end = max(o.exon_ends, that_end)
                    names.add(o.name)
            else:
                if that_chrom:
                    cds_interval_listh.write('{}\t{}\t{}\t{}\t{}\n'.format(that_chrom, that_start, that_end, that_strand, ','.join(names)))
                    #print('{}\t{}\t{}\t{}\t{}'.format(that_chrom, that_start, that_end, that_strand, ','.join(names)))
                that_chrom, that_start, that_end , that_strand = o.chrom, o.exon_starts, o.exon_ends, o.strand
                names = {o.name}
        with open(cds_interval_list_done_file, 'w'):
            pass
        
        #write exon interval interval_list
        print('Writing {}.'.format(exon_interval_list))
        that_chrom = None
        that_start = None
        that_end = None
        that_strand = None
        names = set()
        for o in exon_intervals:
            #chrom, start, end, strand, name = o[1:6]
            if o.chrom == that_chrom and o.strand == that_strand:
                if o.exon_starts > that_end:
                    exon_interval_listh.write('{}\t{}\t{}\t{}\t{}\n'.format(that_chrom, that_start, that_end, that_strand, ','.join(names)))
                    #print('{}\t{}\t{}\t{}\t{}'.format(that_chrom, that_start, that_end, that_strand, ','.join(names)))
                    that_start, that_end = o.exon_starts, o.exon_ends
                    names = {o.name}
                else:
                    that_end = max(o.exon_ends, that_end)
                    names.add(o.name)
            else:
                if that_chrom:
                    exon_interval_listh.write('{}\t{}\t{}\t{}\t{}\n'.format(that_chrom, that_start, that_end, that_strand, ','.join(names)))
                    #print('{}\t{}\t{}\t{}\t{}'.format(that_chrom, that_start, that_end, that_strand, ','.join(names)))
                that_chrom, that_start, that_end , that_strand = o.chrom, o.exon_starts, o.exon_ends, o.strand
                names = {o.name}
        with open(exon_interval_list_done_file, 'w'):
            pass

        #write gene interval interval_list
        print('Writing {}.'.format(gene_interval_list))
        that_chrom = None
        that_start = None
        that_end = None
        that_strand = None
        names = set()
        for o in gene_intervals:
            #chrom, start, end, strand, name = o[1:6]
            if o.chrom == that_chrom and o.strand == that_strand:
                if o.tx_start > that_end:
                    gene_interval_listh.write('{}\t{}\t{}\t{}\t{}\n'.format(that_chrom, that_start, that_end, that_strand, ','.join(names)))
                    #print('{}\t{}\t{}\t{}\t{}'.format(that_chrom, that_start, that_end, that_strand, ','.join(names)))
                    that_start, that_end = o.tx_start, o.tx_end
                    names = {o.name}
                else:
                    that_end = max(o.tx_end, that_end)
                    names.add(o.name)
            else:
                if that_chrom:
                    gene_interval_listh.write('{}\t{}\t{}\t{}\t{}\n'.format(that_chrom, that_start, that_end, that_strand, ','.join(names)))
                    #print('{}\t{}\t{}\t{}\t{}'.format(that_chrom, that_start, that_end, that_strand, ','.join(names)))
                that_chrom, that_start, that_end , that_strand = o.chrom, o.tx_start, o.tx_end, o.strand
                names = {o.name}
        with open(gene_interval_list_done_file, 'w'):
            pass
    print('Interval lists created.')
    for l in (cds_interval_list, exon_interval_list, gene_interval_list):
        print('Creating interval list with only contigs 1-22, X, Y, M {}'.format(l))
        clean_file = my.swap_ext(l, '.all_contigs.interval_list', '.1-22XYMT.interval_list')
        done_file = my.done_file(clean_file)
        interval_list_cleanup.run(l, clean_file)
        print('Created interval list with only contigs 1-22, X, Y, M {}'.format(l))
    print('create_interval_list job done.')


def main():
    """command line arguments"""
    parser = argparse.ArgumentParser(
        description = 'create sorted interval_list file with header',
        epilog = 'refGene2Interval_list 1.0β1 ©2011-2014 Michael Yourshaw all rights reserved')

    parser.add_argument('--input', '-i', required=True, type=str,
        help='input source file',)
    parser.add_argument('--zero_based', action='store_false', default=True,
        help='input file start positions are zero-based (default: True)')
    parser.add_argument('--header_line_count', default=1, type=int,
        help='number of header lines in input file to be skipped (default: 1)')
    parser.add_argument('--header_lines_startwith', default=['#','@'], nargs='*',
        help='list of strings that header lines start with(default: # @)')
    parser.add_argument('--chrom_col', default=2, type=int,
        help='zero-based column number of chromosome column to be converted (default: 2)')
    parser.add_argument('--tx_start_col', default=4, type=int,
        help='zero-based column number of txStart position in input file (default: 4)')
    parser.add_argument('--tx_end_col', default=5, type=int,
        help='zero-based column number of txEnd position in input file (default: 5)')
    parser.add_argument('--cds_start_col', default=6, type=int,
        help='zero-based column number of cdsStart position in input file (default: 6)')
    parser.add_argument('--cds_end_col', default=7, type=int,
        help='zero-based column number of cdsEnd position in input file (default: 7)')
    parser.add_argument('--exon_starts_col', default=9, type=int,
        help='zero-based column number of exonStarts position in input file (default: 9)')
    parser.add_argument('--exon_ends_col', default=10, type=int,
        help='zero-based column number of exonEnds position in input file (default: 10)')
    parser.add_argument('--strand_col', default=3, type=int,
        help='zero-based column number of strand in input file (default: 3; None = all intervals will have strand +)')
    parser.add_argument('--name_cols', default=[12,1], nargs='*',
        help='zero-based column numbers of names (e.g., NR_* and HGNC) in input file (default: 12 1; None = all intervals will have name <chrom>:<start>-<end>)')
    parser.add_argument('--splice_bases', default=2, type=int,
        help='number of additional splice-site base to add  in exon interval_list (default: 2)')
    parser.add_argument('--faidx', required=True, type=str,
        help='faidx file relating hg to GRCh contigs, with GATK sort orders',)
    parser.add_argument('--interval_list_header', required=True, type=str,
        help='picard-style interval_list header file, compatible with contigs in faidx',)
    parser.add_argument('--genome_id', default='b37', type=str,
        help='faidx genome of output file (default: b37',)
    parser.add_argument('--cds_interval_list', '-c', type=str,
        help="output interval_list file for cds (excluding 5' and 3' utr) and splice sites (default=<input>.exon.interval_list")
    parser.add_argument('--exon_interval_list', '-e', type=str,
        help="output interval_list file for exons (including 5' and 3' utr) and splice sites (default=<input>.exon.interval_list")
    parser.add_argument('--gene_interval_list', '-g', type=str,
        help='output interval_list file for entire genes (default=<input>.gene.interval_list')
    
    args = parser.parse_args()
    
    run(input=args.input, cds_interval_list=args.cds_interval_list, exon_interval_list=args.exon_interval_list, gene_interval_list=args.gene_interval_list,
        faidx=args.faidx, interval_list_header=args.interval_list_header, genome_id= args.genome_id, 
        zero_based=args.zero_based, header_line_count=args.header_line_count, header_lines_startwith=args.header_lines_startwith,
        chrom_col=args.chrom_col, strand_col=args.strand_col, name_cols=args.name_cols,splice_bases=args.splice_bases,
        tx_start_col=args.tx_start_col, tx_end_col=args.tx_end_col, cds_start_col=args.cds_start_col, cds_end_col=args.cds_end_col, exon_starts_col=args.exon_starts_col, exon_ends_col=args.exon_ends_col,)

if __name__ == "__main__": sys.exit(main())

