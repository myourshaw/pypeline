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
import my

#--ensembl_db_dir /scratch1/vax/75/ensembl --host cortex.local --port 3306 --user sa --password PASSWORD --ensembl_user ensembl --ensembl_password ensembl
#--ensembl_db_dir /scratch1/vax/75/ensembl --host cortex.local --port 3306 --user sa --password m.cha3ly --ensembl_user ensembl --ensembl_password ensembl

class EnsemblInstallerError(Exception):
    """Log vax_installer-generated exceptions"""
    def __init__(self,msg, log=True):
        self.msg = msg
        if log:
            print_log(msg)
        
class EnsemblInstaller:
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
    
    def run(self, ensembl_db_dir, password=None, config=None, host=None, port=None, user=None, ensembl_user=None, ensembl_password=None,):
        
        if not os.path.isdir(ensembl_db_dir):
            raise EnsemblInstallerError('Ensembl directory not specified', log=False)
        if not config:
            config = my.get_config('vax.cfg')

        #executable sh Command objects
        perl = sh.Command(config.get('DEFAULT','perl'))
        python = sh.Command(config.get('DEFAULT','python'))
        
        #log
        self.log_file = os.path.join(ensembl_db_dir, 'ensembl_db_installer.log')
        try:
            os.remove(self.log_file)
        except OSError:
            pass
        self.print_log(['Starting Ensembl database installer',
            'Ensembl  directory = {}'.format(ensembl_db_dir),
            'Installation log file = {}'.format(self.log_file),
            ])
        
        #MySQL install user and password
        if not user:
                user = config.get('db','USER')
        if not password:
            raise EnsemblInstallerError('A password is required. Installer aborted.')
        if not host:
            host = config.get('db','HOST')
        if not port:
            port = int(config.get('db','PORT'))
                        
        #runtime databases and passwords
        if not ensembl_user:
            ensembl_user = config.get('ensembl','ENSEMBL_USER')
        if not ensembl_password:
            ensembl_password = config.get('ensembl','ENSEMBL_PW')

        self.print_log(['Ensembl database will be installed with these parameters:',
                        'user = {}'.format(user),
                        'host = {}'.format(host),
                        'port = {}'.format(port),
                        'ensembl runtime user = {}'.format(ensembl_user),
                        'ensembl runtime password = {}'.format(ensembl_password),
                        ])
        
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
                raise EnsemblInstallerError('Database connection failed. {}'.format(msg))
            else:
                if not mysql_version:
                    raise EnsemblInstallerError('Database connected but SELECT VERSION() query failed.')
                self.print_log('Database connected. MySQL version: {}'.format(mysql_version))
        self.print_log('Database connection OK')

################################################################
#                   Ensembl MySQL database                     #
################################################################

        done_file = os.path.join(ensembl_db_dir,'.ensembl_db_install.done')
        if my.file_exists(done_file):
            self.print_log('Ensembl MySQL database already installed. To reinstall "rm {}"'.format(done_file))
        else:
            for root,dirs,files in os.walk(ensembl_db_dir):
                break
            database_dirnames = sorted([d for d in dirs
                                        if d != 'mysql' and not d.startswith('ensembl_mart_')
                                        and my.file_exists(os.path.join(ensembl_db_dir, d, d+'.sql'))
                                        ])
            if not database_dirnames:
                raise EnsemblInstallerError('There are no Ensembl database data subdirectories in {}'.format(ensembl_db_dir))
            self.print_log('Installing Ensembl MySQL databases from {}'.format(ensembl_db_dir))
            database_dirs = [os.path.join(ensembl_db_dir, d) for d in database_dirnames]
            for database in database_dirnames:
                db_done_file =my.done_file(os.path.join(ensembl_db_dir, database))
                if my.file_exists(db_done_file):
                    self.print_log('Ensembl MySQL database {} already installed. To reinstall "rm {}"'.format(database, db_done_file))
                else:
                    sql_file = os.path.join(ensembl_db_dir, database, database+'.sql')
                    if my.file_exists(sql_file):
                        with MySQLdb.connect(host,user,password,port=port) as cursor:
                            self.print_log('Creating database `{}`'.format(database))
                            query = ['DROP DATABASE IF EXISTS `{}`'.format(database),
                                     'CREATE DATABASE  `{}`'.format(database),
                                     "GRANT SELECT ON {}.* TO '{}'@'%' IDENTIFIED BY '{}'".format(database, ensembl_user, ensembl_password)]
                            self.mysql(cursor, query)
                            self.print_log('Created database `{}`'.format(database))
                        self.print_log('Creating tables and other objects in database `{}` by executing {}'.format(database, sql_file))
                        with open(sql_file, 'r') as sql:
                            sh.mysql('-h', host, '-P', port, '-u', user, '-p{}'.format(password), database, _in=sql)
                        self.print_log('Created tables and other objects in database `{}` by executing {}'.format(database, sql_file))
                        self.print_log('Importing data. This will take ~20 hours')
                        self.print_log('Populating tables of database `{}`'.format(database))
                        for root, dirs, files in os.walk(os.path.join(ensembl_db_dir, database)):
                            break
                        text_files = [os.path.join(ensembl_db_dir, database, f) for f in files if f.endswith('.txt')]
                        for text_file in text_files:
                            table = my.r_strip(os.path.basename(text_file), '.txt')
                            self.print_log('Importing `{}`.`{}` table.'.format(database, table))
                            sh.mysqlimport('-h', host, '-P', port, '-u', user, '-p{}'.format(password), '--delete', '--local', database, text_file)
                            self.print_log('Imported `{}`.`{}` table.'.format(database, table))
                        with open(db_done_file, 'w'):
                            pass
                        self.print_log('Populated tables of database `{}`'.format(database))
            with MySQLdb.connect(host,user,password,port=port) as cursor:
                self.print_log('Purging MySQL logs.')
                query = "PURGE BINARY LOGS BEFORE '{}'".format(strftime("%Y-%m-%d %H:%M:%S", localtime()))           
                self.mysql(cursor, query)
                self.print_log('Purged MySQL logs.')
            with open(done_file, 'w'):
                pass
            self.print_log('Installed Ensembl MySQL databases from {}'.format(ensembl_db_dir))



def main():

    #command line arguments
    parser = argparse.ArgumentParser(
        description = 'install Ensembl MySQL databases on server',
        epilog = 'vax.ensebml_database_installer 1.0β1 ©2011-2014 Michael Yourshaw all rights reserved')
    #required
    parser.add_argument('--ensembl_db_dir', '-e', type=str, required=True,
        help='directory that contains unzipped subdirectories with Ensembl data and SL scripts (required; example: /scratch1/vax/75/ensembl)')
    #databases
    parser.add_argument('--host', '-H', type=str,
        help='MySQL database server hostname or ip address (default: config.db.HOST, example: cortex.local)')
    parser.add_argument('--port', '-P', type=int,
        help='MySQL database (default: config.db.PORT, example: 3306)')
    parser.add_argument('--user', '-u', type=str,
        help='MySQL database user with priveleges to install Ensembl (default: config.db.USER, example: sa)')
    parser.add_argument('--password', '-p', type=str,
        help='MySQL password for installing user (default: config.db.USER, example: None = enter silently at prompt)')
    parser.add_argument('--ensembl_user', '-U', type=str,
        help='user for runtime access to ensembl databases (default: config.db.ENSEMBL_USER, example: ensembl)')
    parser.add_argument('--ensembl_password', '-W', type=str,
        help='UNSECURE password for runtime access to ensembl databases (default: config.db.ENSEMBL_PW, example: ensembl)')
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
                    raise EnsemblInstallerError('A user is required. Installer aborted.')
                else:
                    #empty user
                    if not args.user:
                        raise EnsemblInstallerError('A user is required. Installer aborted.')
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
                    raise EnsemblInstallerError('A password is required. Installer aborted.')
                else:
                    #empty password
                    if not args.password:
                        raise EnsemblInstallerError('A password is required. Installer aborted.')
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
                    raise EnsemblInstallerError('A host is required. Installer aborted.')
                else:
                    #empty host
                    if not args.host:
                        raise EnsemblInstallerError('A host is required. Installer aborted.')
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
                    raise EnsemblInstallerError('A port is required. Installer aborted.')
                else:
                    #empty port
                    if not args.port:
                        raise EnsemblInstallerError('A port is required. Installer aborted.')    

    EI = EnsemblInstaller()
    EI.run(ensembl_db_dir=args.ensembl_db_dir, password=args.password, config=config,
        host=args.host, port=args.port, user=args.user, ensembl_user=args.ensembl_user, ensembl_password=args.ensembl_password,)

if __name__ == "__main__": sys.exit(main())

