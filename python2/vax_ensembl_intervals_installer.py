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

class EnsemblIntervalsInstallerError(Exception):
    """Log vax_installer-generated exceptions"""
    def __init__(self,msg, log=True):
        self.msg = msg
        if log:
            print_log(msg)
        
class EnsemblIntervalsInstaller:
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

#--ensembl_intervals_dir /share/apps/myourshaw/vax/75/ensembl_intervals/all_genes -d vax_test -H cortex.local -P 3306 -u sa -p 

    def run(self, ensembl_intervals_dir, database, password=None, config=None, host=None, port=None, user=None,):
        
        if not os.path.isdir(ensembl_intervals_dir):
            raise EnsemblIntervalsInstallerError('Ensembl intervals directory not specified', log=False)

        if not config:
            config = my.get_config('vax.cfg')

        #executable sh Command objects
        perl = sh.Command(config.get('DEFAULT','perl'))
        python = sh.Command(config.get('DEFAULT','python'))
        
        #log
        self.log_file = os.path.join(ensembl_intervals_dir, 'ensembl_intervals_installer.log')
        try:
            os.remove(self.log_file)
        except OSError:
            pass
        self.print_log(['Starting Ensembl intervals table installer',
            'Ensembl intervals output directory = {}'.format(ensembl_intervals_dir),
            'Installation log file = {}'.format(self.log_file),
            ])
        
        #MySQL install user and password
        if not user:
                user = config.get('db','USER')
        if not password:
            raise EnsemblIntervalsInstallerError('A password is required. Installer aborted.')
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
                raise EnsemblIntervalsInstallerError('Database connection failed. {}'.format(msg))
            else:
                if not mysql_version:
                    raise EnsemblIntervalsInstallerError('Database connected but SELECT VERSION() query failed.')
                self.print_log('Database connected. MySQL version: {}'.format(mysql_version))
        self.print_log('Database connection OK')

################################################################
#                   Ensembl intervals table                    #
################################################################

        done_file = os.path.join(ensembl_intervals_dir,'.ensembl_intervals_table_install.done')
        if my.file_exists(done_file):
            self.print_log('Ensembl intervals table already installed. To reinstall "rm {}"'.format(done_file))
        else:
            for w in os.walk(ensembl_intervals_dir):
                intervals_file = os.path.join(ensembl_intervals_dir, [f for f in w[2] if f.endswith('.genomic_region_details.txt')][0])
                break
            if not intervals_file:
                raise EnsemblIntervalsInstallerError('Ensembl .genomic_region_details.txt intervals file not found')
            
            text_file = os.path.join(ensembl_intervals_dir, my.r_strip(os.path.basename(intervals_file), '.genomic_region_details.txt')[:64]) +'.txt'
            #grep out comments and get proper table name in file
            grep_done_file = os.path.join(os.path.dirname(text_file), '.'+os.path.basename(text_file)+'.done')
            if my.file_exists(grep_done_file):
                self.print_log('Ensembl intervals grep already done. To redo "rm {}"'.format(grep_done_file))
            else:
                self.print_log('Greping {} to {}.'.format(intervals_file, text_file))
                sh.grep('-e', '^##', '-v', intervals_file, _out=text_file)
                with open(grep_done_file, 'w'):
                    pass
            
            sql_file = text_file+'.mysql'
            table = my.r_strip(os.path.basename(text_file), '.txt')
            self.print_log('Installing Ensembl intervals table `{}`.`{}` from {}'.format(database, table, intervals_file))
            
            sql_done_file = os.path.join(os.path.dirname(sql_file), '.'+os.path.basename(sql_file)+'.done')
            if my.file_exists(sql_done_file):
                self.print_log('Ensembl intervals CREATE TABLE already created. To recreate "rm {}"'.format(sql_done_file))
            else:
                self.print_log('Making CREATE TABLE `{}` script.'.format(table))
                sql_columns.run(text_file, database=database, schema=database, table=table, clustered_index=['Chromosome','Start','End'], indexes=['GeneID','TranscriptID','GeneSymbol'])
                with open(sql_done_file, 'w'):
                    pass

            import_done_file = os.path.join(os.path.dirname(sql_file), '.import.'+os.path.basename(sql_file)+'.done')
            if my.file_exists(import_done_file):
                self.print_log('Ensembl intervals table already imported. To reimport "rm {}"'.format(import_done_file))
            else:
                self.print_log('Creating `{}`.`{}` table.'.format(database, table))
                with open(sql_file, 'r') as sql:
                    sh.mysql('-h', host, '-P', port, '-u', user, '-p{}'.format(password), database, _in=sql)
                
                self.print_log('Importing data into `{}`.`{}`.'.format(database, table))
                sh.mysqlimport('-h', host, '-P', port, '-u', user, '-p{}'.format(password), '--delete', '--local', '--ignore-lines=1', database, text_file)
                self.print_log('Imported data into `{}`.`{}`.'.format(database, table))
                with open(import_done_file, 'w'):
                    pass
                
            with open(done_file, 'w'):
                pass
            self.print_log('Installed Ensembl intervals table {} from {}'.format(table, intervals_file))


def main():

    #command line arguments
    parser = argparse.ArgumentParser(
        description = 'install Ensembl intervals table on server',
        epilog = 'vax.ensebml_intervals_installer 1.0β1 ©2011-2014 Michael Yourshaw all rights reserved')
    #required
    parser.add_argument('--ensembl_intervals_dir', '-e', type=str, required=True,
        help='ensembl intervals directory (required; example: /scratch1/vax/75/ensembl/ensembl_intervals)')
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
                    raise EnsemblIntervalsInstallerError('A user is required. Installer aborted.')
                else:
                    #empty user
                    if not args.user:
                        raise EnsemblIntervalsInstallerError('A user is required. Installer aborted.')
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
                    raise EnsemblIntervalsInstallerError('A password is required. Installer aborted.')
                else:
                    #empty password
                    if not args.password:
                        raise EnsemblIntervalsInstallerError('A password is required. Installer aborted.')
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
                    raise EnsemblIntervalsInstallerError('A host is required. Installer aborted.')
                else:
                    #empty host
                    if not args.host:
                        raise EnsemblIntervalsInstallerError('A host is required. Installer aborted.')
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
                    raise EnsemblIntervalsInstallerError('A port is required. Installer aborted.')
                else:
                    #empty port
                    if not args.port:
                        raise EnsemblIntervalsInstallerError('A port is required. Installer aborted.')    

    EXI = EnsemblIntervalsInstaller()
    EXI.run(ensembl_intervals_dir=args.ensembl_intervals_dir, database=args.database,
        password=args.password, config=config, host=args.host, port=args.port, user=args.user,)

if __name__ == "__main__": sys.exit(main())

