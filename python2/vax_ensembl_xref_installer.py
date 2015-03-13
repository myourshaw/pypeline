#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import argparse
from ConfigParser import SafeConfigParser
from contextlib import closing
import copy
import fileinput
from ftplib import FTP
import ftputil
from getpass import getpass
import gzip
import re
from shlex import split
import shutil
import subprocess
import tarfile
from time import localtime, strftime
from urllib import urlretrieve
from warnings import warn
from zipfile import ZipFile
try:
    from bs4 import BeautifulSoup
except ImportError:
    raise ImportError,"The Python BeautifulSoup module is required to run this program. Try 'pip install beautifulsoup4'."
try:
    from ftputil import FTPHost
except ImportError:
    raise ImportError,"The Python ftputil module is required to run this program. Try 'pip install ftputil'."
try:
    import MySQLdb
    import MySQLdb.cursors
except ImportError:
    raise ImportError,"The Python MySQLdb module is required to run this program. Try 'pip install MySQL-python'."
try:
    import requests
except ImportError:
    raise ImportError,"The Python requests module is required to run this program. Try 'pip install requests'."
try:
    import sh
    from sh import git
except ImportError:
    raise ImportError,"The Python sh module is required to run this program. Try 'pip install sh'."

#these files must be in same directory as this script
import csv2tab
import my
import job
import omim2db
import refgene2db
import sql_columns
import MGI_mouse_phenotype_files
#requires uniprot2db.pl in this directory and SWISS::Entry in PERL5LIB
#requires download_kegg_data.pl

class EnsemblXrefInstallerError(Exception):
    """Log vax_installer-generated exceptions"""
    def __init__(self,msg, log=True):
        self.msg = msg
        if log:
            my.print_log(msg)
        
class EnsemblXrefInstaller:
    """Install Ensembl MySQL databases"""
    log_file = ''
    
    def print_log(self, s):
        my.print_log(s, self.log_file)
        
    def mysql(self, cursor, query):
        """Execute each MySQL statement in query and print the statement(s) to the log"""
        if isinstance(query, basestring):
            query = [query]
        for q in query:
            cursor.execute(q)
            q.replace('\n', ' ')
            self.print_log('Executed: {}'.format(q))
    
    def run(self, ensembl_xref_with_dups_file, database, password=None, config=None, host=None, port=None, user=None,):
        
        if not my.file_exists(ensembl_xref_with_dups_file):
            raise EnsemblXrefInstallerError('Ensembl_xref file not specified', log=False)
        ensembl_xref_dir = os.path.dirname(ensembl_xref_with_dups_file)
        if not config:
            config = my.get_config('vax.cfg')

        #log
        self.log_file = my.swap_ext(ensembl_xref_with_dups_file, '.txt', '.log')
        try:
            os.remove(self.log_file)
        except OSError:
            pass
        self.print_log(['Starting Ensembl_xref table installer',
            'Ensembl_xref file = {}'.format(ensembl_xref_with_dups_file),
            'Installation log file = {}'.format(self.log_file),
            ])
        
        #MySQL install user and password
        if not user:
                user = config.get('db','USER')
        if not password:
            raise EnsemblXrefInstallerError('A password is required. Installer aborted.')
        if not host:
            host = config.get('db','HOST')
        if not port:
            port = int(config.get('db','PORT'))
                        
        #test MySQL connection and create vax database
        self.print_log('Cheking database connection')
        with MySQLdb.connect(host,user,password,port=port) as cursor:
            try:
                query = 'SELECT VERSION()'
                self.mysql(cursor, query)
                mysql_version = cursor.fetchone()
            except MySQLdb.Error, e:
                msg = ''
                try:
                    msg = 'MySQL Error [{}]: {}'.format(e.args[0], e.args[1])
                except IndexError:
                    msg = 'MySQL Error {}'.format(e)
                raise EnsemblXrefInstallerError('Database connection failed. {}'.format(msg))
            else:
                if not mysql_version:
                    raise EnsemblXrefInstallerError('Database connected but SELECT VERSION() query failed.')
                self.print_log('Database connected. MySQL version: {}'.format(mysql_version))
        self.print_log('Database connection OK')

################################################################
#                     Ensembl xref table                       #
################################################################

        done_file = my.done_file(ensembl_xref_with_dups_file+'.install_table')
        if my.file_exists(done_file):
            self.print_log('Ensembl_xref table already installed. To reinstall "rm {}"'.format(done_file))
        else:
            text_file = ensembl_xref_with_dups_file
            sql_file = text_file+'.mysql'
            table = my.r_strip(os.path.basename(text_file), '.txt')
            unique_table = my.r_strip(table, '_with_dups')
            self.print_log('Installing Ensembl xref table `{}`.`{}` from {}'.format(database, unique_table, ensembl_xref_with_dups_file))
            
            self.print_log('Making CREATE TABLE `{}` script.'.format(table))
            sql_columns.run(text_file, database=database, schema=database, table=table, primary_key=['ENSG', 'ENST', 'ENSP', 'dbname', 'primary_id', 'display_id'], indexes=['ENSG', 'ENST', 'ENSP', 'dbname,primary_id', 'dbname,display_id'])
            
            self.print_log('Creating `{}`.`{}` table.'.format(database, table))
            with open(sql_file, 'r') as sql:
                sh.mysql('-h', host, '-P', port, '-u', user, '-p{}'.format(password), database, _in=sql)
                
            self.print_log('Importing data into `{}`.`{}`.'.format(database, table))
            sh.mysqlimport('-h', host, '-P', port, '-u', user, '-p{}'.format(password), '--delete', '--local', '--ignore-lines=1', database, text_file)
            self.print_log('Imported data into `{}`.`{}`.'.format(database, table))
            
            self.print_log('Selecting unique rows from `{0}`.`{1}` into `{0}`.`{2}`.'.format(database, table, unique_table))
            with MySQLdb.connect(host,user,password,port=port) as cursor:
                query = ['DROP TABLE IF EXISTS `{}`.`{}`'.format(database, unique_table),

"""CREATE TABLE `{0}`.`{1}` AS
SELECT DISTINCT ENSG,ENST,ENSP,dbname,primary_id,display_id,chrom,geneStart,geneEnd,strand,txStart,txEnd,cdsStart,cdsEnd,biotype,is_canonical
FROM `{0}`.`{2}`""".format(database, unique_table, table),

"UPDATE `{}`.`{}` SET `ENST` = '' WHERE `ENST` IS NULL".format(database, unique_table),

"UPDATE `{}`.`{}` SET `ENSP` = '' WHERE `ENSP` IS NULL".format(database, unique_table),

"""ALTER TABLE `{}`.`{}`
ADD PRIMARY KEY (ENSG, ENST, ENSP, dbname, primary_id, display_id)".format(database, unique_table)
ADD INDEX `IX_ENST` (`ENST` ASC),
ADD INDEX `IX_ENSP` (`ENSP` ASC),
ADD INDEX `IX_db_name` (`dbname` ASC),
ADD INDEX `IX_display_id` (`display_id` ASC)""".format(database, unique_table),

"DROP TABLE `{}`.`{}`".format(database, table),
]
                self.mysql(cursor, query)
            with open(done_file, 'w'):
                pass
            self.print_log('Installed Ensembl xref table `{}`.`{}` from {}'.format(database, unique_table, ensembl_xref_with_dups_file))

def main():

    #command line arguments
    parser = argparse.ArgumentParser(
        description = 'install Ensembl xref table on server',
        epilog = 'vax.ensebml_xref_installer 1.0β1 ©2011-2014 Michael Yourshaw all rights reserved')
    #required
    parser.add_argument('--ensembl_xref_with_dups_file', '-e', type=str, required=True,
        help='ensembl_xref_<version>_with_dups.txt file (required; example: /scratch1/vax/75/ensembl/ensembl_xref/ensembl_xref_75_with_dups.txt)')
    parser.add_argument('--database', '-d', type=str, required=True,
        help='MySQL database')
    #databases
    parser.add_argument('--host', '-H', type=str,
        help='MySQL database server hostname or ip address (default: config.db.HOST, example: cortex.local)')
    parser.add_argument('--port', '-P', type=int,
        help='MySQL database (default: config.db.PORT, example: 3306)')
    parser.add_argument('--user', '-u', type=str,
        help='MySQL database user with priveleges to install Ensembl (default: config.db.USER, example: sa)')
    parser.add_argument('--password', '-p', type=str,
        help='MySQL password for installing user (default: config.db.USER, example: None = enter silently at prompt)')
    #parse args
    args = parser.parse_args()
    config = my.get_config(args, 'vax.cfg')

    #MySQL install user and password
    if not args.user:
        try:
            args.user = config.get('db','USER')
        except:
            pass
        else:
            if not args.user:
                try:
                    args.user = raw_input('MySQL user with install permissions: ')
                except:
                    raise EnsemblXrefInstallerError('A user is required. Installer aborted.')
                else:
                    #empty user
                    if not args.user:
                        raise EnsemblXrefInstallerError('A user is required. Installer aborted.')
    if not args.password:
        try:
            args.password = config.get('db','PASSWORD')
        except:
            pass
        else:
            if not args.password:
                try:
                    args.password = getpass('MySQL password with install permissions: ')
                except:
                    raise EnsemblXrefInstallerError('A password is required. Installer aborted.')
                else:
                    #empty password
                    if not args.password:
                        raise EnsemblXrefInstallerError('A password is required. Installer aborted.')
    if not args.host:
        try:
            args.host = config.get('db','HOST')
        except:
            pass
        else:
            if not args.host:
                try:
                    args.host = raw_input('MySQL host (e.g., 127.0.0.1, localhost, or cortex.local): ')
                except:
                    raise EnsemblXrefInstallerError('A host is required. Installer aborted.')
                else:
                    #empty host
                    if not args.host:
                        raise EnsemblXrefInstallerError('A host is required. Installer aborted.')
    if not args.port:
        try:
            args.port = int(config.get('db','PORT'))
        except:
            pass
        else:
            if not args.port:
                try:
                    args.port = int(raw_input('MySQL port (e.g., 3306: '))
                except:
                    raise EnsemblXrefInstallerError('A port is required. Installer aborted.')
                else:
                    #empty port
                    if not args.port:
                        raise EnsemblXrefInstallerError('A port is required. Installer aborted.')    

    EXI = EnsemblXrefInstaller()
    EXI.run(ensembl_xref_with_dups_file=args.ensembl_xref_with_dups_file, database=args.database,
        password=args.password, config=config, host=args.host, port=args.port, user=args.user,)

if __name__ == "__main__": sys.exit(main())

