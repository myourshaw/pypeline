#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = 'myourshaw'

import sys
import os
import argparse
import sqlite3
import csv
from my import makedir, done_file, open_gz_or_text

#-a /Volumes/scratch/broad/ExAC_release/release0.2/ExAC.r0.2.sites.vep.vcf.alleles.gz -g /Volumes/scratch/broad/ExAC_release/release0.2/ExAC.r0.2.sites.vep.vcf.genotypes.gz -d /Volumes/scratch/broad/ExAC_release/release0.2 -n exac


def run(allele_file, genotype_file, database_dir, database_name):
    dir = '/Volumes/scratch/broad/ExAC_release/release0.2/'

    makedir(dir)
    database_file= os.path.join(database_dir, database_name)
    database_file_done = done_file(database_file)
    if all((map(os.path.exists, (database_file, database_file_done)))):
        print('Database {} already created. To redo:\nrm {}'.format(os.path.basename(database_name), database_file_done))
        return

    print('Creating ExAC SQLite database.')
    try:
        os.remove(database_file)
    except OSError:
        pass

    with sqlite3.connect(database_file) as c:
        c.execute('''DROP TABLE IF EXISTS alleles;''')
        c.execute('''DROP INDEX IF EXISTS IX_alleles;''')
        c.execute('''CREATE TABLE alleles (
                        CHROM TEXT,
                        POS INTEGER,
                        REF TEXT,
                        allele TEXT,
                        AN INTEGER,
                        AN_AFR INTEGER,
                        AN_AMR INTEGER,
                        AN_EAS INTEGER,
                        AN_FIN INTEGER,
                        AN_NFE INTEGER,
                        AN_OTH INTEGER,
                        AN_SAS INTEGER,
                        AN_Adj INTEGER,
                        AC INTEGER,
                        AC_AFR INTEGER,
                        AC_AMR INTEGER,
                        AC_EAS INTEGER,
                        AC_FIN INTEGER,
                        AC_NFE INTEGER,
                        AC_OTH INTEGER,
                        AC_SAS INTEGER,
                        AC_Adj INTEGER,
                        AF REAL,
                        AF_AFR REAL,
                        AF_AMR REAL,
                        AF_EAS REAL,
                        AF_FIN REAL,
                        AF_NFE REAL,
                        AF_OTH REAL,
                        AF_SAS REAL
                        );'''
        )

        print('Importing ExAC SQLite alleles table.')
        with open_gz_or_text(allele_file) as fh:
            r =csv.DictReader(fh, dialect=csv.excel_tab)
            sql = '''INSERT INTO alleles ({}) values ({})'''.format(
                ','.join(r.fieldnames), ','.join([':'+f for f in r.fieldnames]))
            c.executemany(sql, r)
        c.execute('''CREATE INDEX IX_alleles ON alleles (CHROM, POS, REF, allele);''')

        c.execute('''DROP TABLE IF EXISTS genotypes;''')
        c.execute('''DROP INDEX IF EXISTS IX_genotypes;''')
        c.execute('''CREATE TABLE genotypes (
                        CHROM TEXT,
                        POS INTEGER,
                        REF TEXT,
                        allele1 TEXT,
                        allele2 TEXT,
                        GC_All INTEGER,
                        GC INTEGER,
                        GC_AFR INTEGER,
                        GC_AMR INTEGER,
                        GC_EAS INTEGER,
                        GC_FIN INTEGER,
                        GC_NFE INTEGER,
                        GC_OTH INTEGER,
                        GC_SAS INTEGER
                        );'''
        )

        print('Importing ExAC SQLite genotypes table.')
        with open_gz_or_text(genotype_file) as fh:
            r =csv.DictReader(fh, dialect=csv.excel_tab)
            sql = '''INSERT INTO genotypes ({}) values ({})'''.format(
                ','.join(r.fieldnames), ','.join([':'+f for f in r.fieldnames]))
            c.executemany(sql, r)
        c.execute('''CREATE INDEX IX_genotypes ON genotypes (CHROM, POS, REF, allele1, allele2);''')

    with open(database_file_done, 'w'):
        pass
    print('Created ExAC SQLite database.')


def main():

    parser = argparse.ArgumentParser(
        description = 'create a sqlite database for ExAC allele and genotype data',
        epilog = 'pypeline.exac2db version 1.0β1 ©2011-2014 Michael Yourshaw all rights reserved')
    parser.add_argument('--allele_file', '-a', required=True,
        help='allele file')
    parser.add_argument('--genotype_file', '-g', required=True,
        help='allele file')
    parser.add_argument('--database_dir', '-d', required=True,
        help='database directory')
    parser.add_argument('--database_name', '-n', required=True,
        help='database name')
    args = parser.parse_args()

    run(allele_file=args.allele_file, genotype_file=args.genotype_file,
        database_dir=args.database_dir, database_name=args.database_name)


if __name__ == "__main__": sys.exit(main())
