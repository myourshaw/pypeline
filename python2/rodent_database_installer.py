#!/usr/bin/env python
# -*- coding: utf-8 -*-
__author__ = 'myourshaw'

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
try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET

#these files must be in same directory as this script
import csv2tab
import my
import job
import sql_columns
import vax_mgi_file_cleanup_sql
import vax_download_kegg_data
#vax_ensembl_database_installer.py
#vax_dbNSFP_downloader.py
#vax_ensembl_genomic_regions.pl
#vax_ensembl_intervals_installer.py
#vax_ensembl_xref2db.pl
#vax_ensembl_xref_installer.py
#vax_ensembl_version.pl
#vax_uniprot2db.pl (requires SWISS::Entry in PERL5LIB)

required_programs = ['vax_ensembl_database_installer.py', 'vax_ensembl_version.pl',
                     'vax_uniprot2db.pl','vax_download_kegg_data.py',]

for p in required_programs:
    if not my.file_exists(os.path.join(os.path.dirname(__file__), p)):
        raise Exception('Cannot find {}. The following programs must be in the same directory as this installer: {}'.format(p, ', '.join(prequired_programs)))

config_file = 'rodent.cfg'

install_choices = ['all', 'mouse_ensembl_database', 'rat_ensembl_database', 'kegg', 'mgi', 'rgd', 'uniprot', ]


class RodentInstallerError(Exception):
    """Log rodent_installer-generated exceptions"""
    def __init__(self,msg, log=True):
        self.msg = msg
        if log:
            my.print_log(msg)


class cd():
    """Context manager for changing the current working directory

    Usage:
        with cd(new_dir):
            do something in new_dir, return to current directory when done
    """
    def __init__(self, newPath):
        self.newPath = newPath

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)
    #
    #def __exit__(self, etype, value, traceback):
    #    self.connection(close)


class RodentInstaller:
    """Install Variant Effect Predictor, Ensembl API, data, and plugins for rodent"""

    log_file = ''
    ensembl_jobs = []
    gunzip_jobs = []
    all_jobs = []
    #illegal characters for mysql identifiers
    forbidden = re.compile(r'[^0-9,a-z,A-Z$_]')

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

    def download_ftp_dir(self, host, remote_directory, local_directory, pattern=None, case_insensitive=True):
        """Download all files from an FTP directory; return list of downloaded files"""
        my.makedir(local_directory)
        self.print_log('Downloading {} to {}.'.format(os.path.join(host, remote_directory), local_directory))
        with FTPHost(host, 'anonymous', '') as h:
            h.chdir(remote_directory)
            names = h.listdir(h.curdir)
            if pattern:
                rx = re.compile(pattern, re.I) if case_insensitive else re.compile(pattern)
                names = [n for n in names if rx.match(n)]
            if not names:
                self.print_log('WARNING: nothing to download from {}.'.format(os.path.join(host, remote_directory)))
                return []
            else:
                for name in names:
                    if h.path.isfile(name):
                        name_done_file = os.path.join(local_directory, '.'+name+'.done')
                        if my.file_exists(name_done_file):
                            self.print_log('{} {} already downloaded; skipping. To reinstall "rm {}"'.format(host, name, name_done_file))
                        else:
                            self.print_log('Downloading {} to {}.'.format(os.path.join(host, remote_directory, name), os.path.join(local_directory, name)))
                            h.download(name, os.path.join(local_directory, name))
                            with open(name_done_file, 'w'):
                                pass
                            self.print_log('Downloaded {} to {}.'.format(os.path.join(host, remote_directory, name), os.path.join(local_directory, name)))
                self.print_log('Downloaded {} to {}.'.format(os.path.join(host, remote_directory), local_directory))
                return [os.path.join(local_directory, n) for n in names]

    def getzip(self, url, zipfile, unzipdir):
        """Download zipfile from url, extract contents to unzipdir, and remove zipfile."""
        done_file = os.path.join(unzipdir, '.'+os.path.basename(zipfile)+'.done')
        if my.file_exists(done_file):
            self.print_log('{} already downloaded and extracted; skipping. To reinstall "rm {}"'.format(os.path.basename(zipfile), done_file))
        else:
            self.print_log('Downloading  {} as {}.'.format(url, zipfile))
            urlretrieve(url, zipfile)
            self.print_log('Extracting {} into {}.'.format(zipfile, unzipdir))
            with ZipFile(zipfile, 'r') as zip:
                zip.extractall(unzipdir)
            os.remove(zipfile)
            with open(done_file, 'w'):
                pass

    def run(self, install=['all'], ensembl_version=None, delete_ensembl_version=None, config=None,
            host=None, port=None, user=None, password=None, mouse_database=None, rat_database=None,
            rodent_user=None, rodent_password=None, ensembl_user=None, ensembl_password=None,
            rodent_dir=None,
            overwrite_rodent_database=False, mouse=False, rat=False):

################################################################
#                            setup                             #
################################################################

        for c in install:
            if c not in install_choices:
                raise RodentInstallerError("argument --install: invalid choice: '{}' (choose one or more from '{}')".format(c, "', '".join(install_choices)))
        if not config:
            config = my.get_config(config_file)

        #executable sh Command objects
        perl = sh.Command(config.get('rodent','perl'))
        python = sh.Command(config.get('rodent','python'))

        if not rodent_dir:
            rodent_dir = config.get('rodent', 'RODENT_DIR')
        my.makedir(rodent_dir)

        if mouse and not mouse_database:
            mouse_database = 'mus'
        if rat and not rat_database:
            rat_database = 'rattus'

        log_dir = my.makedir(os.path.join(rodent_dir,'logs'))
        self.log_file = os.path.join(log_dir, 'rodent_installer_{}.log'.format(strftime("%Y%m%d%H%M%S", localtime())))
        try:
            os.remove(self.log_file)
        except OSError:
            pass
        self.print_log(['Starting RODENT installer',
            'Installation log file = {}'.format(self.log_file),
            ])

        #versions of modules
        versions = {}

        species_dict = {'mouse': 'mus_musculus',
                        'rat': 'rattus_norvegicus'}
        #Ensembl version
        #DEBUG:
        #self.print_log('Debugging with Ensembl database/API version 75.')
        #ensembl_current_version = 75
        #ensembl_local_current_version = 75
        self.print_log('Getting current Ensembl database/API version. This will take a minute.')
        vax_ensembl_version_pl = os.path.join(os.path.dirname(__file__), 'vax_ensembl_version.pl')
        for s in (mouse, rat):
            if s:
                species = species_dict['mouse'] if mouse else species_dict['rat']
                pipe = subprocess.Popen([config.get('DEFAULT','perl'), vax_ensembl_version_pl, '--host', 'ensembldb.ensembl.org', '--user', 'anonymous', '--password', '', '--species', species], stdout=subprocess.PIPE)
                v = pipe.stdout.read()
                ensembl_current_version = int(v) if my.is_int(v) else None
                pipe = subprocess.Popen([config.get('DEFAULT','perl'), vax_ensembl_version_pl, '--host', 'cortex.local', '--user', 'ensembl', '--password', 'ensembl', '--species', species], stdout=subprocess.PIPE)
                v = pipe.stdout.read()
                ensembl_local_current_version = int(v) if my.is_int(v) else None
                self.print_log('Local version of ensembl for {} is {}. Current version {} will be installed.'.format(species, ensembl_local_current_version, ensembl_current_version))

                if mouse:
                    if not mouse_ensembl_version or (my.is_int(ensembl_current_version) and my.is_int(ensembl_local_current_version) and mouse_ensembl_version > ensembl_current_version):
                        mouse_ensembl_version = ensembl_current_version
                    else:
                        #TODO: support installing older versions
                        mouse_ensembl_version = ensembl_current_version

                    #older version of Ensembl databases to delete
                    if not delete_mouse_ensembl_version:
                        delete_mouse_ensembl_version = str(int(mouse_ensembl_version)-2) if my.is_int(mouse_ensembl_version) else None
                    if not delete_mouse_ensembl_version:
                        self.print_log('Could not figure out which Ensembl version is two older than {}. Delete older databases manually to free up disk space'.format(mouse_ensembl_version))
                    else:
                        #if not my.query_yes_no('Ensembl database version {} will be deleted from the local MySQL server. Is this OK?'.format(delete_mouse_ensembl_version)):
                            delete_mouse_ensembl_version = None
                if rat:
                    if not rat_ensembl_version or (my.is_int(ensembl_current_version) and my.is_int(ensembl_local_current_version) and rat_ensembl_version > ensembl_current_version):
                        rat_ensembl_version = ensembl_current_version
                    else:
                        #TODO: support installing older versions
                        rat_ensembl_version = ensembl_current_version

                    #older version of Ensembl databases to delete
                    if not delete_rat_ensembl_version:
                        delete_rat_ensembl_version = str(int(rat_ensembl_version)-2) if my.is_int(rat_ensembl_version) else None
                    if not delete_rat_ensembl_version:
                        self.print_log('Could not figure out which Ensembl version is two older than {}. Delete older databases manually to free up disk space'.format(rat_ensembl_version))
                    else:
                        #if not my.query_yes_no('Ensembl database version {} will be deleted from the local MySQL server. Is this OK?'.format(delete_rat_ensembl_version)):
                            delete_rat_ensembl_version = None

        if mouse:
            mouse_version_dir = os.path.join(rodent_dir, 'mouse',mouse_ensembl_version)
            my.makedir(mouse_version_dir)
            mouse_dir_link = os.path.join(mouse_version_dir, 'mouse-current')
            try:
                os.remove(mouse_dir_link)
            except OSError:
                pass
            with cd(rodent_dir):
                os.symlink(mouse_ensembl_version, 'mouse-current')

            self.print_log('mouse installation directory = {}'.format(mouse_version_dir),)

            mouse_ensembl_db_dir = my.makedir(os.path.join(mouse_version_dir,'ensembl'))

            mouse_ensembl_intervals_dir = my.makedir(os.path.join(mouse_version_dir,'ensembl_intervals', 'all_genes'))

            mouse_ensembl_protein_coding_intervals_dir = my.makedir(os.path.join(mouse_version_dir,'ensembl_intervals', 'protein_coding_genes'))

            mouse_ensembl_xref_dir = my.makedir(os.path.join(mouse_version_dir,'ensembl_xref'))

            mouse_kegg_dir = my.makedir(os.path.join(mouse_version_dir,'kegg'))

            mouse_mgi_dir = my.makedir(os.path.join(mouse_version_dir,'mgi'))

            mouse_rgd_dir = my.makedir(os.path.join(mouse_version_dir,'rgd'))

        if rat:
            rat_version_dir = os.path.join(rodent_dir, 'rat',rat_ensembl_version)
            my.makedir(rat_version_dir)
            rat_dir_link = os.path.join(rat_version_dir, 'rat-current')
            try:
                os.remove(rat_dir_link)
            except OSError:
                pass
            with cd(rodent_dir):
                os.symlink(rat_ensembl_version, 'rat-current')

            self.print_log('rat installation directory = {}'.format(rat_version_dir),)

            rat_ensembl_db_dir = my.makedir(os.path.join(rat_version_dir,'ensembl'))

            rat_ensembl_intervals_dir = my.makedir(os.path.join(rat_version_dir,'ensembl_intervals', 'all_genes'))

            rat_ensembl_protein_coding_intervals_dir = my.makedir(os.path.join(rat_version_dir,'ensembl_intervals', 'protein_coding_genes'))

            rat_ensembl_xref_dir = my.makedir(os.path.join(rat_version_dir,'ensembl_xref'))

            rat_kegg_dir = my.makedir(os.path.join(rat_version_dir,'kegg'))

            rat_mgi_dir = my.makedir(os.path.join(rat_version_dir,'mgi'))

            rat_rgd_dir = my.makedir(os.path.join(rat_version_dir,'rgd'))

        job_dir = my.makedir(os.path.join(rodent_dir, 'jobs'))

        tmp = my.makedir(os.path.join(rodent_dir, 'tmp'))

        uniprot_dir = my.makedir(os.path.join(rodent_dir,'uniprot'))

        #MySQL install user and password
        if not user:
            try:
                user = config.get('db','USER')
            except:
                pass
            else:
                if not user:
                    try:
                        user = raw_input('MySQL user with install permissions: ')
                    except:
                        raise RodentInstallerError('A user is required. Installer aborted.')
                    else:
                        #empty user
                        if not user:
                            raise RodentInstallerError('A user is required. Installer aborted.')
        if not password:
            try:
                password = config.get('db','PASSWORD')
            except:
                pass
            else:
                if not password:
                    try:
                        password = getpass('MySQL password with install permissions: ')
                    except:
                        raise RodentInstallerError('A password is required. Installer aborted.')
                    else:
                        #empty password
                        if not password:
                            raise RodentInstallerError('A password is required. Installer aborted.')
        if not host:
            try:
                host = config.get('db','HOST')
            except:
                pass
            else:
                if not host:
                    try:
                        host = raw_input('MySQL host (e.g., 127.0.0.1, localhost, or cortex.local): ')
                    except:
                        raise RodentInstallerError('A host is required. Installer aborted.')
                    else:
                        #empty host
                        if not host:
                            raise RodentInstallerError('A host is required. Installer aborted.')
        if not port:
            try:
                port = int(config.get('db','PORT'))
            except:
                pass
            else:
                if not port:
                    try:
                        port = int(raw_input('MySQL port (e.g., 3306: '))
                    except:
                        raise RodentInstallerError('A port is required. Installer aborted.')
                    else:
                        #empty port
                        if not port:
                            raise RodentInstallerError('A port is required. Installer aborted.')
        if rat and not mouse_database:
            try:
                database = config.get('rodent','MOUSE_DATABASE_FORMAT_STRING').format(ensembl_version)
            except:
                pass
            else:
                if not mouse_database:
                    try:
                        mouse_database = raw_input('MySQL rodent database (e.g., mus_75: ')
                    except:
                        raise RodentInstallerError('A database is required. Installer aborted.')
                    else:
                        #empty database
                        if not mouse_database:
                            raise RodentInstallerError('A database is required. Installer aborted.')

        if rat and not rat_database:
            try:
                rat_database = config.get('rodent','MOUSE_DATABASE_FORMAT_STRING').format(ensembl_version)
            except:
                pass
            else:
                if not rat_database:
                    try:
                        rat_database = raw_input('MySQL rodent database (e.g., mus_75: ')
                    except:
                        raise RodentInstallerError('A database is required. Installer aborted.')
                    else:
                        #empty database
                        if not rat_database:
                            raise RodentInstallerError('A database is required. Installer aborted.')

        #runtime databases and passwords
        if not rodent_user:
            rodent_user = config.get('rodent','RODENT_USER')
        if not rodent_password:
            rodent_password = config.get('rodent','RODENT_PW')
        if not ensembl_user:
            ensembl_user = config.get('ensembl','ENSEMBL_USER')
        if not ensembl_password:
            ensembl_password = config.get('ensembl','ENSEMBL_PW')

        self.print_log(['Database will be installed with these parameters:',
            'user = {}'.format(user),
            'host = {}'.format(host),
            'port = {}'.format(port),
            'Ensembl version = {}'.format(ensembl_version),
            'rodent runtime user = {}'.format(rodent_user),
            'rodent runtime password = {}'.format(rodent_password),
            ])

        if mouse:
            self.print_log(['mouse version = {}'.format(mouse_ensembl_version),
                'mouse database name = {}'.format(mouse_database),
                ])
        if rat:
            self.print_log(['rat version = {}'.format(rat_ensembl_version),
                'rat database name = {}'.format(rat_database),
                ])


################################################################
# test MySQL connection and create mouse and/or rat databases  #
################################################################

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
                raise RodentInstallerError('Database connection failed. {}'.format(msg))
            else:
                if not mysql_version:
                    raise RodentInstallerError('Database connected but SELECT VERSION() query failed.')
                self.print_log('Database connected. MySQL version: {}'.format(mysql_version))
                query = "SELECT SCHEMA_NAME FROM INFORMATION_SCHEMA.SCHEMATA WHERE SCHEMA_NAME = '{}'".format(database)
                self.mysql(cursor, query)
                schema = cursor.fetchone()
                self.print_log('Database `{}` {}.'.format(database, "doesn't exist" if schema == None else "exists"))

                self.print_log('Purging MySQL logs.')
                with MySQLdb.connect(host,user,password,database,port=port) as cursor:
                    query = "PURGE BINARY LOGS BEFORE '{}'".format(strftime("%Y-%m-%d %H:%M:%S", localtime()))
                    self.mysql(cursor, query)

                if mouse:
                    #get existing database's versions table
                    versions = my.get_key_value_dict(host, user, password, port, mouse_database, table='versions', key_col='module', value_col='version',)
                    versions['mouse'] = str(mouse_ensembl_version)
                    #create new database with empty versions table
                    if schema == None or overwrite_rodent_database:
                        self.print_log('Creating new database `{}`. Any existing data and procedures will be deleted.'.format(mouse_database))
                        query = ['DROP DATABASE IF EXISTS `{}`'.format(mouse_database),
                                 'CREATE DATABASE `{}`'.format(mouse_database),
                                 """CREATE TABLE `{}`.`versions` (
`module` VARCHAR(255) NOT NULL,
`version` VARCHAR(255) NULL,
PRIMARY KEY (`module`))""".format(mouse_database,),]
                        self.mysql(cursor, query)
                        self.print_log('Created new database `{}`.'.format(mouse_database))
                    else:
                        self.print_log('Retaining database `{}`, but overwriting tables and procedures for selected modules.'.format(mouse_database))
                    query = "GRANT SELECT, EXECUTE ON {}.* TO '{}'@'%' IDENTIFIED BY '{}'".format(mouse_database, rodent_user, rodent_password)
                    self.mysql(cursor, query)
                if rat:
                    #get existing database's versions table
                    versions = my.get_key_value_dict(host, user, password, port, rat_database, table='versions', key_col='module', value_col='version',)
                    versions['rat'] = str(rat_ensembl_version)
                    #create new database with empty versions table
                    if schema == None or overwrite_rodent_database:
                        self.print_log('Creating new database `{}`. Any existing data and procedures will be deleted.'.format(rat_database))
                        query = ['DROP DATABASE IF EXISTS `{}`'.format(rat_database),
                                 'CREATE DATABASE `{}`'.format(rat_database),
                                 """CREATE TABLE `{}`.`versions` (
`module` VARCHAR(255) NOT NULL,
`version` VARCHAR(255) NULL,
PRIMARY KEY (`module`))""".format(rat_database,),]
                        self.mysql(cursor, query)
                        self.print_log('Created new database `{}`.'.format(rat_database))
                    else:
                        self.print_log('Retaining database `{}`, but overwriting tables and procedures for selected modules.'.format(rat_database))
                    query = "GRANT SELECT, EXECUTE ON {}.* TO '{}'@'%' IDENTIFIED BY '{}'".format(rat_database, rodent_user, rodent_password)
                    self.mysql(cursor, query)

        self.print_log('Database connection OK')


################################################################
#                   Ensembl MySQL databases                    #
################################################################


        if mouse and ('all' in install or 'mouse_ensembl_database' in install):
            done_file = os.path.join(mouse_ensembl_db_dir,'.schedule_mouse_ensembl_db_download.done')
            if my.file_exists(done_file):
                self.print_log('Download Ensembl mouse MySQL databases cluster jobs already scheduled. To reschedule "rm {}"'.format(done_file))
            else:
                ensembl_host = config.get('ensembl', 'ENSEMBL_HOST')
                ensembl_ftp_dir = config.get('ensembl', 'ENSEMBL_FTP_DIR')
                ensembl_ftp_pattern = config.get('ensembl', 'ENSEMBL_MOUSE_FTP_PATTERN')
                ensembl_rsync_url =  config.get('ensembl', 'ENSEMBL_RSYNC_URL')
                #filter out unnecessary files/databases
                ensembl_blacklist_patterns =  config.get('ensembl', 'ENSEMBL_BLACKLIST_PATTERNS').split()
                rx = re.compile('({})'.format(')|('.join(ensembl_blacklist_patterns)))
                ftp_files = [f for f in my.list_ftp_files(ensembl_host, ensembl_ftp_dir, ensembl_ftp_pattern) if not rx.match(f)]
                #rsync a few very large databases separately, the rest in a single rsync to avoid overloading the server
                ensembl_rsync_parallel_patterns =  config.get('ensembl', 'ENSEMBL_RSYNC_PARALLEL_PATTERNS').split()
                rx = re.compile('({})'.format(')|('.join(ensembl_rsync_parallel_patterns)))
                #list of rsync batches; batches will be scheduled in parallel jobs]
                rsyncs = [[os.path.join(ensembl_rsync_url,f)] for f in ftp_files if rx.match(f)]
                rsyncs.append([os.path.join(ensembl_rsync_url,f) for f in ftp_files if not rx.match(f)])
                self.print_log('Deleting older Ensembl databases.')
                if delete_ensembl_version:
                    for file in ftp_files:
                        old_file = file.replace('_'+ensembl_version, '_'+delete_ensembl_version)
                        with MySQLdb.connect(host,user,password,port=port) as cursor:
                            self.print_log('Deleting old database `{}`'.format(old_file))
                            query = 'DROP DATABASE IF EXISTS `{}`'.format(old_file)
                            self.mysql(cursor, query)
                            self.print_log('Deleted old database `{}`'.format(old_file))
                self.print_log('Deleted older Ensembl databases.')
                self.print_log('Scheduling cluster jobs to download and unzip Ensembl MySQL database files.')
                for batch in rsyncs:
                    cmd = 'rsync -av "{}" "{}/"'.format('" "'.join(batch), mouse_ensembl_db_dir)
                    job_name = 'rsync_ensembl_{}_db.format(ensembl_version)'
                    job = my.run_job(cmd, job_name, job_dir)
                    ensembl_job_id = job.jobId
                    self.ensembl_jobs.append(job.jobId)
                    self.all_jobs.append(job.jobId)
                    self.print_log('Scheduled job {}: {}.'.format(job.jobId, cmd))
                    for file in batch:
                        ensembl_subdir = os.path.join(mouse_ensembl_db_dir, os.path.basename(file))
                        cmd = 'gunzip -rf {}/*'.format(ensembl_subdir)
                        job_name = 'gunzip_ensembl_{}_db'.format(mouse_ensembl_version)
                        job = my.run_job(cmd, job_name, job_dir, hold_jid = ensembl_job_id)
                        gunzip_job_id = job.jobId
                        self.gunzip_jobs.append(job.jobId)
                        self.ensembl_jobs.append(job.jobId)
                        self.all_jobs.append(gunzip_job_id)
                        self.print_log('Scheduled job {}: {}.'.format(job.jobId, cmd))
                with open(done_file, 'w'):
                    pass
                self.print_log('Scheduled cluster jobs to download and unzip Ensembl MySQL database files.')

            done_file = os.path.join(ensembl_db_dir,'.schedule_ensembl_db_import.done')
            if my.file_exists(done_file):
                self.print_log('Import Ensembl MySQL databases cluster jobs already scheduled. To reschedule "rm {}"'.format(done_file))
            else:
                self.print_log('Scheduling job to create new  Ensembl database version {} and import data.'.format(ensembl_version))
                rodent_ensembl_database_installer_py = os.path.join(os.path.dirname(__file__), 'rodent_ensembl_database_installer.py')
                cmd= '{} {} --ensembl_db_dir {} --host {} --port {} --user {} --password {} --ensembl_user {} --ensembl_password {} '.format(
                    config.get('DEFAULT','python'), rodent_ensembl_database_installer_py, ensembl_db_dir, host, port, user, password, ensembl_user, ensembl_password)
                job_name = 'ensembl_{}_db_import'.format(ensembl_version)
                job = my.run_job(cmd, job_name, job_dir, hold_jid=self.gunzip_jobs)
                ensembl_import_job_id = job.jobId
                self.ensembl_jobs.append(job.jobId)
                self.all_jobs.append(job.jobId)
                self.print_log('Scheduled job {}: {}.'.format(job.jobId, cmd))

                with open(done_file, 'w'):
                    pass
                self.print_log('Scheduled cluster jobs {} to install Ensembl MySQL database.'.format(','.join(self.ensembl_jobs)))

        if rat and ('all' in install or 'rat_ensembl_database' in install):
            done_file = os.path.join(rat_ensembl_db_dir,'.schedule_rat_ensembl_db_download.done')
            if my.file_exists(done_file):
                self.print_log('Download Ensembl rat MySQL databases cluster jobs already scheduled. To reschedule "rm {}"'.format(done_file))
            else:
                ensembl_host = config.get('ensembl', 'ENSEMBL_HOST')
                ensembl_ftp_dir = config.get('ensembl', 'ENSEMBL_FTP_DIR')
                ensembl_ftp_pattern = config.get('ensembl', 'ENSEMBL_RAT_FTP_PATTERN')
                ensembl_rsync_url =  config.get('ensembl', 'ENSEMBL_RSYNC_URL')
                #filter out unnecessary files/databases
                ensembl_blacklist_patterns =  config.get('ensembl', 'ENSEMBL_BLACKLIST_PATTERNS').split()
                rx = re.compile('({})'.format(')|('.join(ensembl_blacklist_patterns)))
                ftp_files = [f for f in my.list_ftp_files(ensembl_host, ensembl_ftp_dir, ensembl_ftp_pattern) if not rx.match(f)]
                #rsync a few very large databases separately, the rest in a single rsync to avoid overloading the server
                ensembl_rsync_parallel_patterns =  config.get('ensembl', 'ENSEMBL_RSYNC_PARALLEL_PATTERNS').split()
                rx = re.compile('({})'.format(')|('.join(ensembl_rsync_parallel_patterns)))
                #list of rsync batches; batches will be scheduled in parallel jobs]
                rsyncs = [[os.path.join(ensembl_rsync_url,f)] for f in ftp_files if rx.match(f)]
                rsyncs.append([os.path.join(ensembl_rsync_url,f) for f in ftp_files if not rx.match(f)])
                self.print_log('Deleting older Ensembl databases.')
                if delete_ensembl_version:
                    for file in ftp_files:
                        old_file = file.replace('_'+ensembl_version, '_'+delete_ensembl_version)
                        with MySQLdb.connect(host,user,password,port=port) as cursor:
                            self.print_log('Deleting old database `{}`'.format(old_file))
                            query = 'DROP DATABASE IF EXISTS `{}`'.format(old_file)
                            self.mysql(cursor, query)
                            self.print_log('Deleted old database `{}`'.format(old_file))
                self.print_log('Deleted older Ensembl databases.')
                self.print_log('Scheduling cluster jobs to download and unzip Ensembl MySQL database files.')
                for batch in rsyncs:
                    cmd = 'rsync -av "{}" "{}/"'.format('" "'.join(batch), rat_ensembl_db_dir)
                    job_name = 'rsync_ensembl_{}_db.format(ensembl_version)'
                    job = my.run_job(cmd, job_name, job_dir)
                    ensembl_job_id = job.jobId
                    self.ensembl_jobs.append(job.jobId)
                    self.all_jobs.append(job.jobId)
                    self.print_log('Scheduled job {}: {}.'.format(job.jobId, cmd))
                    for file in batch:
                        ensembl_subdir = os.path.join(rat_ensembl_db_dir, os.path.basename(file))
                        cmd = 'gunzip -rf {}/*'.format(ensembl_subdir)
                        job_name = 'gunzip_ensembl_{}_db'.format(rat_ensembl_version)
                        job = my.run_job(cmd, job_name, job_dir, hold_jid = ensembl_job_id)
                        gunzip_job_id = job.jobId
                        self.gunzip_jobs.append(job.jobId)
                        self.ensembl_jobs.append(job.jobId)
                        self.all_jobs.append(gunzip_job_id)
                        self.print_log('Scheduled job {}: {}.'.format(job.jobId, cmd))
                with open(done_file, 'w'):
                    pass
                self.print_log('Scheduled cluster jobs to download and unzip Ensembl MySQL database files.')

            done_file = os.path.join(ensembl_db_dir,'.schedule_ensembl_db_import.done')
            if my.file_exists(done_file):
                self.print_log('Import Ensembl MySQL databases cluster jobs already scheduled. To reschedule "rm {}"'.format(done_file))
            else:
                self.print_log('Scheduling job to create new  Ensembl database version {} and import data.'.format(ensembl_version))
                rodent_ensembl_database_installer_py = os.path.join(os.path.dirname(__file__), 'rodent_ensembl_database_installer.py')
                cmd= '{} {} --ensembl_db_dir {} --host {} --port {} --user {} --password {} --ensembl_user {} --ensembl_password {} '.format(
                    config.get('DEFAULT','python'), rodent_ensembl_database_installer_py, ensembl_db_dir, host, port, user, password, ensembl_user, ensembl_password)
                job_name = 'ensembl_{}_db_import'.format(ensembl_version)
                job = my.run_job(cmd, job_name, job_dir, hold_jid=self.gunzip_jobs)
                ensembl_import_job_id = job.jobId
                self.ensembl_jobs.append(job.jobId)
                self.all_jobs.append(job.jobId)
                self.print_log('Scheduled job {}: {}.'.format(job.jobId, cmd))

                with open(done_file, 'w'):
                    pass
                self.print_log('Scheduled cluster jobs {} to install Ensembl MySQL database.'.format(','.join(self.ensembl_jobs)))

################################################################
#                    Mouse Phenotypes (MGI)                    #
################################################################

        if mouse and ('all' in install or 'mouse' in install):
            done_file = os.path.join(mgi_dir,'.mgi.done')
            if my.file_exists(done_file):
                self.print_log('Mouse Phenotypes already installed. To reinstall "rm {}"'.format(done_file))
            else:
                self.print_log('Installing Mouse Phenotypes.')

                versions['mgi'] = 'accessed {}'.format(strftime("%Y-%m-%d", localtime()))
                human_url = config.get('mgi', 'MGI_HUMAN_URL')
                human_table = os.path.splitext(os.path.basename(human_url))[0]
                human_rpt_file = os.path.join(mgi_dir, os.path.basename(human_url))
                human_text_file = os.path.splitext(human_rpt_file)[0]+'.txt'
                human_sql_file = human_text_file+'.mysql'
                self.print_log('Downloading {} to {}.'.format(human_url, human_rpt_file))
                urlretrieve(human_url, human_rpt_file)
                self.print_log('Downloaded {} to {}.'.format(human_url, human_rpt_file))

                mouse_url = config.get('mgi', 'MGI_MOUSE_URL')
                mouse_table = os.path.splitext(os.path.basename(mouse_url))[0]
                mouse_rpt_file = os.path.join(mgi_dir, os.path.basename(mouse_url))
                mouse_text_file = os.path.splitext(mouse_rpt_file)[0]+'.txt'
                mouse_sql_file = mouse_text_file+'.mysql'
                self.print_log('Downloading {} to {}.'.format(mouse_url, mouse_rpt_file))
                urlretrieve(mouse_url, mouse_rpt_file)
                self.print_log('Downloaded {} to {}.'.format(mouse_url, mouse_rpt_file))

                self.print_log('Converting {} and {} to {} and {}.'.format(human_rpt_file, mouse_rpt_file, human_text_file, mouse_text_file))
                rodent_MGI_mouse_phenotype_files.run(human_rpt_file, mouse_rpt_file, human_text_file, mouse_text_file)
                self.print_log('Converted {} and {} to {} and {}.'.format(human_rpt_file, mouse_rpt_file, human_text_file, mouse_text_file))

                self.print_log('Making CREATE TABLE `{}` script.'.format(human_table))
                sql_columns.run(human_text_file, database=database, schema=database, table=human_table, indexes=['Mammalian_Phenotype_ID', 'Human_Gene'])
                self.print_log('Made CREATE TABLE `{}` script.'.format(human_table))

                self.print_log('Making CREATE TABLE `{}` script.'.format(mouse_table))
                sql_columns.run(mouse_text_file, database=database, schema=database, table=mouse_table, indexes=['Mammalian_Phenotype_ID'])
                self.print_log('Made CREATE TABLE `{}` script.'.format(mouse_table))

                self.print_log('Creating `{}` table.'.format(human_table))
                with open(human_sql_file, 'r') as sql:
                    sh.mysql('-h', host, '-P', port, '-u', user, '-p{}'.format(password), database, _in=sql)
                self.print_log('Importing `{}` data to `{}` database.'.format(human_table, database))
                sh.mysqlimport('-h', host, '-P', port, '-u', user, '-p{}'.format(password), '--delete', '--local', '--ignore-lines=1', database, human_text_file)
                self.print_log('Imported `{}` to `{}` database.'.format(human_table, database))


                self.print_log('Creating `{}` table.'.format(mouse_table))
                with open(mouse_sql_file, 'r') as sql:
                    sh.mysql('-h', host, '-P', port, '-u', user, '-p{}'.format(password), database, _in=sql)
                self.print_log('Importing `{}` data to `{}` database.'.format(mouse_table, database))
                sh.mysqlimport('-h', host, '-P', port, '-u', user, '-p{}'.format(password), '--delete', '--local', '--ignore-lines=1', database, mouse_text_file)
                self.print_log('Imported `{}` to `{}` database.'.format(mouse_table, database))

                proc = 'hgnc2mouse_phenotype'
                self.print_log('Creating {} stored procedure'.format(proc))
                query = [
'DROP PROCEDURE IF EXISTS {}'.format(proc),
"""CREATE DEFINER=CURRENT_USER PROCEDURE {}(hgnc VARCHAR(15))
BEGIN
    SELECT DISTINCT v.`Mammalian_Phenotype_Name`
    from vw.`VOC_MammalianPhenotype` as v
    join vw.`HMD_HumanPhenotype` as h
    on v.`Mammalian_Phenotype_ID` = h.`Mammalian_Phenotype_ID`
    WHERE h.`Human_Gene` = hgnc
    ORDER BY v.`Mammalian_Phenotype_Name`;
END
""".format(proc)]
                with MySQLdb.connect(host,user,password,database,port=port) as cursor:
                    self.mysql(cursor, query)
                    self.print_log('Created {} stored procedure'.format(proc))

                with open(done_file, 'w'):
                    pass
                self.print_log('Installed Mouse Phenotypes.')


################################################################
#            RGD (Rat Genome Database) ontologies              #
################################################################

        if 'all' in install or 'ontologies' in install:
            done_file = os.path.join(rgd_dir,'.rgd.done')
            if my.file_exists(done_file):
                self.print_log('RGD ontologies already installed. To reinstall "rm {}"'.format(done_file))
            else:
                self.print_log('Installing RGD ontologies. This will take ~25 minutes.')

                #download RGD directory
                rgd_ftp = config.get('rgd', 'RGD_FTP')
                rgd_ftp_dir = config.get('rgd', 'RGD_FTP_DIR')
                rgd_pattern = config.get('rgd', 'RGD_PATTERN')
                files = self.download_ftp_dir(rgd_ftp, rgd_ftp_dir, rgd_dir, pattern=rgd_pattern)
                rgd_version = ''
                if files:
                    for f in files:
                        if not rgd_version:
                            with open(f) as fh:
                                for line in fh:
                                    if line.startswith('#'):
                                        g = '# GENERATED-ON: '
                                        if line.startswith(g):
                                            rgd_version = 'RGD '+line.rstrip('/n')[len(g):]
                                            versions['rgd'] = 'RGD accessed {}'.format(rgd_version)
                                    else:
                                        break

                        table = 'rgd_'+os.path.basename(f)
                        text_file = os.path.join(os.path.dirname(f), table+'.txt')
                        sql_file = text_file+'.mysql'
                        #grep out comments and get proper table name in file
                        grep_done_file = os.path.join(os.path.dirname(text_file), '.'+os.path.basename(text_file)+'.done')
                        if my.file_exists(grep_done_file):
                            self.print_log('RGD ontologies grep already done. To redo "rm {}"'.format(grep_done_file))
                        else:
                            self.print_log('Greping {} to {}.'.format(f, text_file))
                            sh.grep('-e', '^#', '-v', f, _out=text_file)
                            with open(grep_done_file, 'w'):
                                pass
                            self.print_log('Created {}.'.format(text_file))

                        sql_done_file = os.path.join(os.path.dirname(sql_file), '.'+os.path.basename(sql_file)+'.done')
                        if my.file_exists(sql_done_file):
                            self.print_log('{} already created. To redo "rm {}"'.format(sql_file, sql_done_file))
                        else:
                            self.print_log('Making CREATE TABLE `{}` script.'.format(table))
                            sql_columns.run(text_file, database=database, schema=database, table=table, indexes=['OBJECT_SYMBOL'], header_line=1)
                            with open(sql_done_file, 'w'):
                                pass
                            self.print_log('Made CREATE TABLE `{}` script.'.format(table))

                        self.print_log('Creating `{}` table.'.format(table))
                        with open(sql_file, 'r') as sql:
                            sh.mysql('-h', host, '-P', port, '-u', user, '-p{}'.format(password), database, _in=sql)

                        self.print_log('Importing data to `{}`.`{}` table.'.format(database, table))
                        sh.mysqlimport('-h', host, '-P', port, '-u', user, '-p{}'.format(password), '--delete', '--local', '--ignore-lines=1', database, text_file)
                        self.print_log('Imported data to `{}`.`{}` table.'.format(database, table))

                    proc = 'hgnc2rgd_ontologies'
                    self.print_log('Creating `{0}`.`{1}` stored procedure'.format(database, proc))
                    query = [
'DROP PROCEDURE IF EXISTS `{0}`.`{1}`'.format(database, proc),
"""CREATE DEFINER=CURRENT_USER PROCEDURE `{0}`.`{1}`(hgnc VARCHAR(255))
BEGIN
SELECT
(SELECT GROUP_CONCAT(DISTINCT TERM_NAME SEPARATOR '|') FROM `rgd_homo_terms_bp` WHERE OBJECT_SYMBOL = hgnc) AS biological_process,
(SELECT GROUP_CONCAT(DISTINCT TERM_NAME SEPARATOR '|') FROM `rgd_homo_terms_cc` WHERE OBJECT_SYMBOL = hgnc) AS cellular_component,
(SELECT GROUP_CONCAT(DISTINCT TERM_NAME SEPARATOR '|') FROM `rgd_homo_terms_chebi` WHERE OBJECT_SYMBOL = hgnc) AS ChEBI_ontology,
(SELECT GROUP_CONCAT(DISTINCT TERM_NAME SEPARATOR '|') FROM `rgd_homo_terms_hp` WHERE OBJECT_SYMBOL = hgnc) AS human_phenotype,
(SELECT GROUP_CONCAT(DISTINCT TERM_NAME SEPARATOR '|') FROM `rgd_homo_terms_mf` WHERE OBJECT_SYMBOL = hgnc) AS molecular_function,
(SELECT GROUP_CONCAT(DISTINCT TERM_NAME SEPARATOR '|') FROM `rgd_homo_terms_mp` WHERE OBJECT_SYMBOL = hgnc) AS mammalian_phenotype,
(SELECT GROUP_CONCAT(DISTINCT TERM_NAME SEPARATOR '|') FROM `rgd_homo_terms_nbo` WHERE OBJECT_SYMBOL = hgnc) AS neuro_behavioral_ontology,
(SELECT GROUP_CONCAT(DISTINCT TERM_NAME SEPARATOR '|') FROM `rgd_homo_terms_pw` WHERE OBJECT_SYMBOL = hgnc) AS pathway_ontology,
(SELECT GROUP_CONCAT(DISTINCT TERM_NAME SEPARATOR '|') FROM `rgd_homo_terms_rdo` WHERE OBJECT_SYMBOL = hgnc) AS RGD_disease_ontology
;
END
""".format(database, proc)
]
                    with MySQLdb.connect(host,user,password,database,port=port) as cursor:
                        self.mysql(cursor, query)

                    self.print_log('Created `{0}`.`{1}` stored procedure'.format(database, proc))

                with open(done_file, 'w'):
                    pass
                self.print_log('Installed RGD ontologies.')


################################################################
#                           UniProt                            #
################################################################

        #required by KEGG
        if 'all' in install or 'uniprot' in install:
            done_file = os.path.join(uniprot_dir,'.uniprot.done')
            if my.file_exists(done_file):
                self.print_log('UniProt already installed. To reinstall "rm {}"'.format(done_file))
            else:
                self.print_log('Installing UniProt.')
                #download UniProt data
                versions['uniprot'] = 'UniProt accessed {}'.format(strftime("%Y-%m-%d", localtime()))
                files = []
                if mouse or rat:
                    files.append(
                        config.get('uniprot', 'UNIPROT_SPROT_RODENTS'),
                        config.get('uniprot', 'UNIPROT_TREMBL_RODENTS'),
                    )
                if mouse:
                    files.append(
                        config.get('uniprot', 'UNIPROT_IDMAPPING_DAT_MOUSE'),
                        config.get('uniprot', 'UNIPROT_IDMAPPING_Selected_MOUSE')
                    )
                if rat:
                    files.append(
                        config.get('uniprot', 'UNIPROT_IDMAPPING_DAT_RAT'),
                        config.get('uniprot', 'UNIPROT_IDMAPPING_Selected_RAT')
                    )
                    for file in files:
                        with FTPHost(config.get('uniprot', 'UNIPROT_FTP'), 'anonymous', '') as host:
                            self.print_log('Downloading {}.'.format(file))
                            host.chdir(os.path.dirname(file))
                            host.download_if_newer(os.path.basename(file), os.path.join(uniprot_dir, os.path.basename(file)))
                            with open(my.done_file(file), 'w'):
                                pass
                            self.print_log('Downloaded {}.'.format(file))

                #create tab-delimited UniProt files
                self.print_log('Create Uniprot tab-delimited files. This will take ~15 min')
                vax_uniprot2db_pl = os.path.join(os.path.dirname(__file__), 'rodent_uniprot2db.pl')
                perl(vax_uniprot2db_pl, '-dir', uniprot_dir)
                self.print_log('Uniprot tab-delimited files created.')

                #create Uniprot table definitions and tables, and import UniProt data
                self.print_log('Creating Uniprot MySQL tables. This will take ~1 hr.')
                for table in['uniprot_rodents_protein',
                             'uniprot_rodents_protein',
                             'uniprot_rodents_protein',]:
                    text_file = os.path.join(uniprot_dir, table+'.txt')
                    sql_file = text_file+'.mysql'
                    self.print_log('Making CREATE TABLE `{}` script.'.format(table))
                    sql_columns.run(text_file, database=database, schema=database, table=table, indexes=['UniProtKB_AC'])
                    self.print_log('Creating `{}` table.'.format(table))
                    with open(sql_file, 'r') as sql:
                        sh.mysql('-h', host, '-P', port, '-u', user, '-p{}'.format(password), database, _in=sql)
                    self.print_log('Importing data to `{}`.`{}` table.'.format(database, table))
                    sh.mysqlimport('-h', host, '-P', port, '-u', user, '-p{}'.format(password), '--delete', '--local', '--ignore-lines=1', database, text_file)
                    self.print_log('Imported data to `{}`.`{}` table.'.format(database, table))

                with MySQLdb.connect(host,user,password,database,port=port) as cursor:

                    #create enst_uniprot tables for rodent access
                    self.print_log('Creating enst_uniprot table')
                    query = [
'DROP TABLE IF EXISTS `enst_uniprot`',
"""CREATE TABLE `enst_uniprot`
(KEY IX_ENST_UniProtKB_AC_topic(`ENST`,`UniProtKB_AC`,`topic`),
KEY IX_ENST_topic (`ENST`,`topic`),
KEY IX_UniProtKB_AC_topic (`UniProtKB_AC`,`topic`))
ENGINE=MyISAM DEFAULT CHARSET=latin1
SELECT CAST(xe.`RESOURCE_IDENTIFIER` AS CHAR(15)) ENST, u.*
FROM `uniprot_human_protein` u
JOIN `uniprot_human_xref` xe
ON u.`UniProtKB_AC` = xe.`UniProtKB_AC` AND xe.`RESOURCE_ABBREVIATION` = 'Ensembl' AND xe.`RESOURCE_IDENTIFIER` LIKE 'ENST%'
"""]
                    self.mysql(cursor, query)

                    self.print_log('Creating enst_uniprot_feature table')
                    query = [
'DROP TABLE IF EXISTS `enst_uniprot_feature`',
"""CREATE TABLE enst_uniprot_feature
(KEY IX_ENST_feature_aaStart_aaEnd (ENST,feature,aaStart,aaEnd),
KEY IX_UniProtKB_AC_feature_aaStart_aaEnd (UniProtKB_AC,feature,aaStart,aaEnd))
ENGINE=MyISAM DEFAULT CHARSET=latin1
SELECT DISTINCT CAST(xe.RESOURCE_IDENTIFIER AS CHAR(15)) ENST, u.*
FROM uniprot_human_protein_feature u
JOIN uniprot_human_xref xe
ON u.UniProtKB_AC = xe.UniProtKB_AC AND xe.RESOURCE_ABBREVIATION = 'Ensembl' AND xe.RESOURCE_IDENTIFIER LIKE 'ENST%'
"""]
                    self.mysql(cursor, query)

                    #create UniProt stored procedures for rodent access
                    self.print_log('Creating enst2uniprot stored procedure')
                    query = [
'DROP PROCEDURE IF EXISTS enst2uniprot',
"""CREATE DEFINER=CURRENT_USER PROCEDURE enst2uniprot(enst CHAR(15))
BEGIN
SELECT u.topic, u.value
FROM enst_uniprot AS u
WHERE u.ENST = enst
AND COALESCE(u.value, '') <> ''
ORDER BY u.topic;
END
"""]
                    self.mysql(cursor, query)

                    self.print_log('Creating enst2uniprot_feature stored procedure')
                    query = [
'DROP PROCEDURE IF EXISTS enst2uniprot_feature',
"""CREATE DEFINER=CURRENT_USER PROCEDURE enst2uniprot_feature(enst CHAR(15), aaStart INT, aaEnd INT)
BEGIN
SELECT u.feature, u.aaStart, u.aaEnd, u.description
FROM enst_uniprot_feature AS u
WHERE u.ENST = enst
AND ((aaStart BETWEEN u.aaStart and u.aaEnd)
OR  (aaEnd BETWEEN u.aaStart and u.aaEnd))
ORDER BY u.feature, u.aaStart, u.aaEnd;
END
"""]
                    self.mysql(cursor, query)

                with open(done_file, 'w'):
                    pass
                self.print_log('Installed UniProt.')


################################################################
#                             KEGG                             #
################################################################

        #requires prior UniProt installation

        if 'all' in install or 'kegg' in install:
            done_file = os.path.join(kegg_dir,'.kegg.done')
            if my.file_exists(done_file):
                self.print_log('KEGG already installed. To reinstall "rm {}"'.format(done_file))
            else:
                self.print_log('Installing KEGG. This will take ~10 mins.')

                #download KEGG pathways
                versions['kegg'] = 'KEGG accessed {}'.format(strftime("%Y-%m-%d", localtime()))
                rodent_download_kegg_data_pl = os.path.join(os.path.dirname(__file__), 'rodent_download_kegg_data.pl')
                kegg_file = os.path.join(kegg_dir, 'kegg_gene_pathways.txt')
                self.print_log('Downloading KEGG data to {}.'.format(kegg_file))
                perl(rodent_download_kegg_data_pl, '-o', kegg_file)
                self.print_log('Downloaded KEGG data.')

                #create KEGG table definitions and tables, and import UniProt data
                table = 'kegg_gene_pathways'
                text_file = os.path.join(kegg_dir, table+'.txt')
                sql_file = text_file+'.mysql'
                self.print_log('Making CREATE TABLE `{}` script.'.format(table))
                sql_columns.run(text_file, database=database, schema=database, table=table, indexes=['gene_id'])
                self.print_log('Creating `{}` table.'.format(table))
                with open(sql_file, 'r') as sql:
                    sh.mysql('-h', host, '-P', port, '-u', user, '-p{}'.format(password), database, _in=sql)
                self.print_log('Importing data to `{}`.`{}` table.'.format(database, table))
                sh.mysqlimport('-h', host, '-P', port, '-u', user, '-p{}'.format(password), '--delete', '--local', '--ignore-lines=1', database, text_file)
                self.print_log('Imported data to `{}`.`{}` table.'.format(database, table))

                with MySQLdb.connect(host,user,password,database,port=port) as cursor:

                    #create KEGG table for rodent access
                    table = 'ensg_kegg'
                    self.print_log('Creating {} table'.format(table))
                    query = [
'DROP TABLE IF EXISTS `{}`'.format(table),
"""CREATE TABLE `{}`
(PRIMARY KEY (ENSG,gene_id,path_id))
ENGINE=MyISAM DEFAULT CHARSET=latin1
SELECT DISTINCT CAST(xe.`OPTIONAL_INFORMATION_2` AS CHAR(15)) ENSG, k.`gene_id`, k.`path_id`, LTRIM(RTRIM(k.pathway)) pathway
FROM `kegg_gene_pathways` AS k
JOIN `uniprot_human_xref` AS x
ON k.`gene_id` = x.`RESOURCE_IDENTIFIER` AND x.`RESOURCE_ABBREVIATION` = 'KEGG'
JOIN `uniprot_human_xref` AS xe
ON x.`UniProtKB_AC` = xe.`UniProtKB_AC` AND xe.`RESOURCE_ABBREVIATION` = 'Ensembl' AND xe.`OPTIONAL_INFORMATION_2` LIKE 'ENSG%'
WHERE IFNULL(LTRIM(RTRIM(k.`pathway`)), '') <> ''
""".format(table)]
                    self.mysql(cursor, query)

                    #create KEGG stored procedure for rodent access
                    proc = 'ensg2kegg'
                    self.print_log('Creating {} stored procedure'.format(proc))
                    query = [
'DROP PROCEDURE IF EXISTS {}'.format(proc),
"""CREATE DEFINER=CURRENT_USER PROCEDURE {}(ensg CHAR(15))
BEGIN
SELECT DISTINCT k.`pathway`
FROM `ensg_kegg` AS k
WHERE k.`ENSG` = ensg
ORDER BY k.`path_id`;
END
""".format(proc)]
                    self.mysql(cursor, query)

                with open(done_file, 'w'):
                    pass
                self.print_log('Installed KEGG.')


################################################################
#                                                              #
################################################################


################################################################
#                            Finis                             #
################################################################


        self.print_log('Rodent database installer done.')


def main():

    #command line arguments
    parser = argparse.ArgumentParser(
        description = 'Installs rodent (mouse, actually) database from Ensembl MGI, RGD, UniProt, and KEGG on mysql server',
        epilog = 'rodent.rodent_database_installer 1.0β1 ©2014 Michael Yourshaw all rights reserved.')
    parser.add_argument('--install', nargs='*', default=['all'], type=str, metavar='MODULE', choices = install_choices,
        help='Modules to install. Choose one or more from "{}" (default: all)'.format(' '.join(install_choices)))
    #Ensembl versions
    parser.add_argument('--ensembl_version', '-e', type=int,
        help='Version of Ensembl API/VEP/cache to install. At present, the latest version will be installed regardless of this parameter. (default: current version; example: 75)')
    parser.add_argument('--delete_ensembl_version', type=str,
        help='old version of Ensembl databases to delete (default: ensembl_version - 2)')
    #databases
    parser.add_argument('--host', '-H', type=str,
        help='MySQL database server hostname or ip address (default: config.db.HOST, example: cortex.local)')
    parser.add_argument('--port', '-P', type=int,
        help='MySQL database (default: config.db.PORT, example: 3306)')
    parser.add_argument('--user', '-u', type=str,
        help='MySQL database user with priveleges to installt (default: config.db.USER, example: sa)')
    parser.add_argument('--password', '-p', type=str,
        help='MySQL password for installing user (default: config.db.USER, example: None = enter silently at prompt)')
    parser.add_argument('--mouse_database', '-m', type=str,
        help='mouse MySQL database; EXISTING DATABASE WILL BE OVERWRITTEN! (default: config.rodent.MOUSE_DATABASE_FORMAT_STRING_{ensembl_version}, example: mus)')
    parser.add_argument('--rat_database', '-r', type=str,
        help='rat MySQL database; EXISTING DATABASE WILL BE OVERWRITTEN! (default: config.rodent.RAT_DATABASE_FORMAT_STRING_{ensembl_version}, example: rattus)')
    parser.add_argument('--rodent_user', '-U', type=str,
        help='user for runtime access to rodent database (default: config.rodent.RODENT_USER, example: rodent)')
    parser.add_argument('--rodent_password', '-W', type=str,
        help='UNSECURE password for runtime access to rodent database (default: config.rodent.RODENT_PW, example: rodent)')
    parser.add_argument('--ensembl_user', type=str,
        help='user for runtime access to ensembl databases (default: config.ensembl.ENSEMBL_USER, example: ensembl)')
    parser.add_argument('--ensembl_password', type=str,
        help='UNSECURE password for runtime access to ensembl databases (default: config.ensembl.ENSEMBL_PW, example: ensembl)')
    #top level directory
    parser.add_argument('--rodent_dir', '-d', type=str,
        help='top-level directory for rodent installation; parent of mus and rattus (default: config.rodent.RODENT_DIR)')
    #versions
    #options
    parser.add_argument('--overwrite_rodent_database', action='store_true', default=False,
        help='overwrite mouse and rat database,s deleting all data and procedures (default: if mus or rattus databases exist, just overwrite specific tables and procedures installed by this application')
    parser.add_argument('--mouse', action='store_true', default=False,
        help='install mouse data (default: False)')
    parser.add_argument('--rat', action='store_true', default=False,
        help='install rat data NOT IMPLEMENTED (default: False)')

    #parse args
    args = parser.parse_args()
    print """These modules will be installed or updated:
    {}""".format(', '.join(args.install))

    config = my.get_config(args, config_file)

    R = RodentInstaller()
    R.run(install=args.install, ensembl_version=args.ensembl_version, delete_ensembl_version=args.delete_ensembl_version, config=config,
        host=args.host, port=args.port, user=args.user, password=args.password, mouse_database=args.mouse_database, rat_database=args.rat_database,
        rodent_user=args.rodent_user, rodent_password=args.rodent_password, ensembl_user=args.ensembl_user, ensembl_password=args.ensembl_password,
        rodent_dir=args.rodent_dir,
        overwrite_rodent_database=args.overwrite_rodent_database, mouse=args.mouse, rat=args.rat)

if __name__ == "__main__": sys.exit(main())

