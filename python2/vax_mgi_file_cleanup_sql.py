#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import argparse
from warnings import warn
import my
import sql_columns

class MGIMousePhenotypeFilesError(Exception): pass

#-H /scratch1/tmp/myourshaw/mouse_phenoypes/HMD_HumanPhenotype.rpt -V /scratch1/tmp/myourshaw/mouse_phenoypes/VOC_MammalianPhenotype.rmwhitespace.rpt

name = 'pypeline.MGI_mouse_phenotype_files'
version = 1.1
copyright = '©2013-2014 Michael Yourshaw all rights reserved'
run_time = my.localtime_stamp()

def run(input, database=None, schema=None):

    if not (database or schema):
        raise MGIMousePhenotypeFilesError('Database and/or schema parameters required.')
    if not database:
        database = schema
    if not schema:
        schema = database

    tables = {
        'MGI_AllGenes': {'header': ('MGI_Gene_ID', 'Gene_Symbol', 'Gene_Name', 'KOMP_Target', 'Regeneron', 'CSD', 'Other_Knockout_Source', 'Chr', 'Start_Coordinate', 'End_Coordinate', 'Strand', 'NCBI_ID', 'Ensembl_ID', 'VEGA_ID', 'cCDS'),
                         'primary_key': ['Gene_Symbol',],
                         'indexes': ['Ensembl_ID',],
                         },
        'HMD_HumanPhenotype': {'header': ('Human_Gene', 'Human_Entrez_Gene_ID', 'HomoloGene_ID', 'Mouse_Marker_Symbol', 'MGI_Marker_ID', 'Mammalian_Phenotype_ID'),
                         'primary_key': ['Human_Gene', 'Mouse_Marker_Symbol',],
                         'indexes': ['Mouse_Marker_Symbol', 'HomoloGene_ID'],
                         },
        'HOM_MouseHumanSequence': {'header': ('HomoloGene_ID', 'Common_Organism_Name', 'NCBI_Taxon_ID', 'Symbol', 'EntrezGene_ID', 'Mouse_MGI_ID', 'HGNC_ID', 'OMIM_Gene_ID', 'Genetic_Location', 'Genomic_Coordinates_mouse_human', 'Nucleotide_RefSeq IDs', 'Protein_RefSeq_IDs', 'SWISS_PROT_IDs'),
                         'primary_key': ['HomoloGene_ID', 'NCBI_Taxon_ID',],
                         'indexes': ['Symbol',],
                         },
        'VOC_MammalianPhenotype': {'header': ('Mammalian_Phenotype_ID', 'Mammalian_Phenotype_Name', 'Description'),
                         'primary_key': ['Mammalian_Phenotype_ID',],
                         'indexes': [],
                         },
    }

    for rpt in input:
        table = os.path.splitext(os.path.basename(rpt))[0]
        #add header, remove empty last field, trim fields, flatten phenotypes
        header = tables[table]['header']
        line_count = 0
        data_input_line_count = 0
        data_output_line_count = 0
        print 'Processing {}.'.format(rpt)
        fixed_rpt = my.swap_ext(rpt, '.rpt', '.txt')
        with open(rpt) as rpt, open(fixed_rpt, 'w') as fixed:
            fixed.write('\t'.join(header)+'\n')
            for line in rpt:
                line_count += 1
                line = line.rstrip('\n')
                if not bool(line.strip()) or line.startswith('#'):
                    continue
                else: #data line
                    data_input_line_count += 1
                    fields = [f.strip() for f in line.split('\t')]
                    if table == 'HMD_HumanPhenotype':
                        #expand multiple phenotypes to one per output record
                        phenotypes = [p.strip() for p in fields[5].split()]
                        if not phenotypes:
                            phenotypes = ['']
                        for phenotype in phenotypes:
                            output_line = '\t'.join(fields[:5]) + '\t' + phenotype + '\n'
                            fixed.write(output_line)
                            data_output_line_count += 1
                    else:
                        output_line = '\t'.join(fields) + '\n'
                        fixed.write(output_line)
                        data_output_line_count += 1
        print 'Read {} input lines, {} data input lines.\nWrote {} data output lines to {}'.format(line_count, data_input_line_count, data_output_line_count, fixed_rpt)
        print('Making CREATE TABLE `{}` script.'.format(table))
        sql_columns.run(fixed_rpt, database=database, schema=database, table=table, primary_key=tables[table]['primary_key'], indexes=tables[table]['indexes'], header_line=1,)
        print('Made CREATE TABLE `{}` script.'.format(table))

def main():

    #command line arguments
    parser = argparse.ArgumentParser(
        description = 'cleanup files from MGI-Mouse Genome Informatics',
        epilog = '{} version {} {}'.format(name, version, copyright))
    parser.add_argument('--input', '-i', nargs='+',
        help='input .rpt files')
    parser.add_argument('--database', '-d', type=str,
        help='database')
    parser.add_argument('--schema', '-s', type=str,
        help='database schema')
    args = parser.parse_args()

    run(input=args.input, database=args.database, schema=args.schema,)

if __name__ == "__main__": sys.exit(main())


