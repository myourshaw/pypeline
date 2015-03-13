#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-d vax_test -e 75 -p XXXXX
#-d vax_test -e 75 -p XXXXX --install ensembl_database kegg mgi rgd uniprot
#-p XXXXX


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
import json
import re
from shlex import split
import shutil
import subprocess
import tarfile
from time import localtime, strftime
from urllib import urlretrieve
from urllib2 import urlopen
from warnings import warn
from zipfile import ZipFile
try:
    from bs4 import BeautifulSoup
except ImportError:
    raise ImportError("The Python BeautifulSoup module is required to run this program. Try 'pip install beautifulsoup4'.")
try:
    from ftputil import FTPHost
except ImportError:
    raise ImportError("The Python ftputil module is required to run this program. Try 'pip install ftputil'.")
try:
    import MySQLdb
    import MySQLdb.cursors
except ImportError:
    raise ImportError("The Python MySQLdb module is required to run this program. Try 'pip install MySQL-python'.")
try:
    import requests
except ImportError:
    raise ImportError("The Python requests module is required to run this program. Try 'pip install requests'.")
try:
    import sh
    from sh import git
except ImportError:
    raise ImportError("The Python sh module is required to run this program. Try 'pip install sh'.")
try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET

#these files must be in same directory as this script
import csv2tab
import my
import job
import omim2db
import refgene2db
import sql_columns
import convert_contigs_sort
import create_interval_list
import vax_dbsnp_vcfs2db
import vax_download_kegg_data
import vax_metabolome_xml2db
import vax_mgi_file_cleanup_sql
#vax_CADD_downloader.py
#vax_dbsnp_installer.py
#vax_ensembl_database_installer.py
#vax_dbNSFP_downloader.py
#vax_ensembl_genomic_regions.pl
#vax_ensembl_intervals_installer.py
#vax_ensembl_xref2db.pl
#vax_ensembl_xref_installer.py
#vax_ensembl_version.pl
#vax_nhlbi_vcfs2db.py
#vax_uniprot2db.pl (requires SWISS::Entry in PERL5LIB)

required_programs = ['vax_dbsnp_installer.py', 'vax_download_kegg_data.py', 'vax_ensembl_database_installer.py',
                     'vax_ensembl_genomic_regions.pl', 'vax_ensembl_intervals_installer.py', 'vax_ensembl_xref2db.pl',
                     'vax_ensembl_xref_installer.py', 'vax_ensembl_version.pl', 'vax_nhlbi_vcfs2db.py', 'vax_uniprot2db.pl',]

for p in required_programs:
    if not my.file_exists(os.path.join(os.path.dirname(__file__), p)):
        raise Exception('Cannot find {}. The following programs must be in the same directory as this installer or in PYTHONPATH: {}'.format(p, ', '.join(required_programs)))
    
config_file = 'vax.cfg'

install_choices = ['all', 'bioperl', 'cadd', 'dbnsfp', 'dbsnp', 'ensembl', 'ensembl_api', 'ensembl_database',
                   'ensembl_intervals', 'ensembl_xref', 'faidx', 'fathmm', 'gatk', 'hgmd', 'kegg', 'human_protein_atlas',
                   'hpa',  'metabolome', 'hmd', 'mitocarta', 'mgi', 'nhlbi', 'evs', 'omim', 'ontologies', 'rgd',
                   'refgene', 'refseq', 'uniprot', '1000_genomes', '1kg', 'vax_plugins', 'vep_cache', 'vep_ini',
                   'vep_plugins',]


class VaxInstallerError(Exception):
    """Log vax_installer-generated exceptions"""
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

class VaxInstaller:
    """Install Variant Effect Predictor, Ensembl API, data, and plugins for VAX"""
    
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

    def run(self, config=None, install=['all'],
            ensembl_version=None, delete_ensembl_version=None, ensembl_species=['homo_sapiens', 'mus_musculus'],
            host=None, port=None, user=None, password=None, database=None,
            vax_user=None, vax_password=None, ensembl_user=None, ensembl_password=None,
            hgmd_user=None, hgmd_password=None, hgmd_version=None, hgmd_dir=None, install_hgmd_schemas=None, hgmd_database=None,
            kegg_species=['hsa', 'mmu'], rgd_genera=['homo', 'mus'],
            uniprot_taxonomic_divisions=['homo_sapiens', 'mus_musculus'], uniprot_organisms=['HUMAN_9606', 'MOUSE_10090'],
            vax_dir=None, scratch_dir=None,
            bioperl_version=None, cadd_snp_file_url=None, dbnsfp_file=None, fathmm_file=None, nhlbi_vcf_file_url=None,
            overwrite_vax_database=False,):

################################################################
#                            setup                             #
################################################################

        for c in install:
            if c not in install_choices:
                raise VaxInstallerError("argument --install: invalid choice: '{}' (choose one or more from '{}')".format(c, "', '".join(install_choices)))
        if not config:
            config = my.get_config(config_file)
        
        #executable sh Command objects
        perl = sh.Command(config.get('DEFAULT','perl'))
        python = sh.Command(config.get('DEFAULT','python'))
        
        #directories and log
        apps_dir = config.get('DEFAULT', 'apps')
        my.makedir(apps_dir)
        
        if not vax_dir:
            vax_dir = config.get('DEFAULT', 'VAX_DIR')
        my.makedir(vax_dir)
        
        if not scratch_dir:
            scratch_dir = config.get('vax', 'SCRATCH_DIR')
        my.makedir(scratch_dir)

        log_dir = my.makedir(os.path.join(vax_dir,'logs'))
        self.log_file = os.path.join(log_dir, 'vax_installer_{}.log'.format(strftime("%Y%m%d%H%M%S", localtime())))
        try:
            os.remove(self.log_file)
        except OSError:
            pass
        self.print_log(['Starting VAX installer',
            'Installation log file = {}'.format(self.log_file),
            ])

        self.print_log(['These modules will be installed or updated:']+ install)

        #MySQL install user and password; if not in a parameter or the config file, ask user
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
                        raise VaxInstallerError('A user is required. Installer aborted.')
                    else:
                        #empty user
                        if not user:
                            raise VaxInstallerError('A user is required. Installer aborted.')
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
                        raise VaxInstallerError('A password is required. Installer aborted.')
                    else:
                        #empty password
                        if not password:
                            raise VaxInstallerError('A password is required. Installer aborted.')
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
                        raise VaxInstallerError('A host is required. Installer aborted.')
                    else:
                        #empty host
                        if not host:
                            raise VaxInstallerError('A host is required. Installer aborted.')
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
                        raise VaxInstallerError('A port is required. Installer aborted.')
                    else:
                        #empty port
                        if not port:
                            raise VaxInstallerError('A port is required. Installer aborted.')

        self.print_log('Getting current Ensembl database/API version. This will take a minute.')
       #current ensembl database versions at ensembl.org
        ensembl_software_release_url = config.get('ensembl','ENSEMBL_SOFTWARE_RELEASE_URL')
        ensembl_org_software_release_current_version = int(json.load(urlopen(ensembl_software_release_url))['release'])
        self.print_log('ensembl.org software release current version = {}'.format(ensembl_org_software_release_current_version))
        ensembl_assembly_version_url_pattern = config.get('ensembl','ENSEMBL_ASSEMBLY_VERSION_URL_PATTERN')
        ensembl_org_assembly_versions = {}
        ensembl_installed_api_versions = {}
        vax_ensembl_version_pl = os.path.join(os.path.dirname(__file__), 'vax_ensembl_version.pl')
        for s in ensembl_species:
            ensembl_org_assembly_versions[s] = json.load(urlopen(ensembl_assembly_version_url_pattern.format(species=s)))
            self.print_log('ensembl.org {} assembly default coord system version = {}'.format(s, ensembl_org_assembly_versions[s]['default_coord_system_version']))
        #TODO: A better solution to failure if vax-current link has been updated but perl api yet not installed
        try:
            for s in ensembl_species:
                pipe = subprocess.Popen([config.get('DEFAULT','perl'), vax_ensembl_version_pl, '--host', host, '--port', str(port), '--user', user, '--password', password, '--species', s], stdout=subprocess.PIPE)
                v = pipe.stdout.read()
                if v and '\t' in v:
                    v = v.split('\t')
                    v[0] = int(v[0]) if my.is_int(v[0]) else None
                    ensembl_installed_api_versions[s] = {'api_version': int(v[0]), 'genome_version': v[1]}
                else:
                    ensembl_installed_api_versions[s] = {'api_version': None, 'genome_version': None}
                self.print_log('ensembl installed {} database/api version = {}'.format(s, ensembl_installed_api_versions[s]['api_version']))
                self.print_log('ensembl installed {} database genome version = {}'.format(s, ensembl_installed_api_versions[s]['genome_version']))
        except:
            self.print_log('Failed to get information about current Ensembl Perl API and database, perhaps because they have not yet been installed.')

        if not ensembl_version:
            ensembl_version = ensembl_org_software_release_current_version

        #delete older version of Ensembl database (2 less than current)
        if not delete_ensembl_version:
            delete_ensembl_version = str(min(ensembl_version, ensembl_org_software_release_current_version)-2) if ensembl_version and my.is_int(ensembl_version) else None
        if not delete_ensembl_version:
            self.print_log('Could not figure out which Ensembl version is two older than {}. Delete older databases manually to free up disk space'.format(ensembl_version))
        else: #TODO: uncomment next line for production build to implement deleting
            #if not my.query_yes_no('Ensembl database version {} will be deleted from the local MySQL server. Is this OK?'.format(delete_ensembl_version)):
                delete_ensembl_version = None

        #vax version
        vax_version = str(ensembl_version)
        
        if not database:
            try:
                database = config.get('vax','VAX_DATABASE_FORMAT_STRING').format(ensembl_version)
            except:
                pass
            else:
                if not database:
                    try:
                        database = raw_input('MySQL vax database (e.g., vax_77: ')
                    except:
                        raise VaxInstallerError('A database is required. Installer aborted.')
                    else:
                        #empty database
                        if not database:
                            raise VaxInstallerError('A database is required. Installer aborted.')

        #runtime databases and passwords
        if not vax_user:
            vax_user = config.get('vax','VAX_USER')
        if not vax_password:
            vax_password = config.get('vax','VAX_PW')
        if not ensembl_user:
            ensembl_user = config.get('ensembl','ENSEMBL_USER')
        if not ensembl_password:
            ensembl_password = config.get('ensembl','ENSEMBL_PW')
        if not hgmd_user:
            hgmd_user = config.get('hgmd','HGMD_USER')
        if not hgmd_password:
            hgmd_password = config.get('hgmd','HGMD_PW')        #versions of modules
        versions = {}

        vax_version_dir = os.path.join(vax_dir, vax_version)
        my.makedir(vax_version_dir)
        vax_version_dir_link = os.path.join(vax_dir, 'vax-current')
        with cd(vax_dir):
            try:
                os.remove(os.path.basename(vax_version_dir_link))
            except OSError:
                pass
            os.symlink(vax_version, os.path.basename(vax_version_dir_link))
        
        scratch_version_dir = os.path.join(scratch_dir, vax_version)
        my.makedir(scratch_version_dir)
        scratch_version_dir_link = os.path.join(scratch_dir, 'scratch-current')
        with cd(scratch_dir):
            try:
                os.remove(os.path.basename(scratch_version_dir_link))
            except OSError:
                pass
            os.symlink(vax_version, os.path.basename(scratch_version_dir_link))

        with cd(vax_dir):
            try:
                os.remove(os.path.basename(scratch_version_dir_link))
            except OSError:
                pass
            os.symlink(scratch_version_dir_link, os.path.basename(scratch_version_dir_link))

        #links to modules that must be in vax_version_dir
        vax_modules = config.get('DEFAULT', 'VAX_MODULES').split(';')
        for m in vax_modules:
            module_link = os.path.join(vax_version_dir_link, os.path.basename(m))
            try:
                os.remove(module_link)
            except OSError:
                pass
            with cd(vax_version_dir_link):
                os.symlink(m, os.path.basename(m))

        plugins_version_dir = my.makedir(os.path.join(vax_version_dir_link,'Plugins'))
        plugins_dir_link = os.path.join(vax_dir,'Plugins')
        try:
            os.remove(plugins_dir_link)
        except OSError:
            pass
        with cd(vax_dir):
            os.symlink(os.path.join(os.path.basename(vax_version_dir_link), 'Plugins'), 'Plugins')

        self.print_log(['vax installation directory = {}'.format(vax_version_dir_link),
            'vax installation scratch_directory = {}'.format(scratch_version_dir),
            'Plugins directory = {}.'.format(plugins_version_dir),
            ])
        
        bioperl_dir = my.makedir(os.path.join(vax_version_dir_link,'bioperl'))
        bioperl_live_dir = my.makedir(os.path.join(bioperl_dir,'bioperl-live'))
        bioperl_live_link = os.path.join(vax_dir,'bioperl-live')
        try:
            os.remove(bioperl_live_link)
        except OSError:
            pass
        with cd(vax_dir):
            os.symlink(os.path.join(os.path.basename(vax_version_dir_link), 'bioperl', 'bioperl-live'), 'bioperl-live')

        cadd_dir = my.makedir(os.path.join(scratch_version_dir,'cadd'))
        dbnsfp_dir = my.makedir(os.path.join(scratch_version_dir,'dbnsfp'))
        dbsnp_dir = my.makedir(os.path.join(scratch_version_dir,'dbsnp'))
        
        ensembl_api_dir = my.makedir(os.path.join(vax_version_dir_link,'ensembl'))
        ensembl_api_dir_link = os.path.join(vax_dir,'ensembl')
        try:
            os.remove(ensembl_api_dir_link)
        except OSError:
            pass
        with cd(vax_dir):
            os.symlink(os.path.join(os.path.basename(vax_version_dir_link), 'ensembl'), 'ensembl')
        
        ensembl_db_dir = my.makedir(os.path.join(scratch_version_dir,'ensembl'))
        ensembl_intervals_dir = my.makedir(os.path.join(vax_version_dir_link,'ensembl_intervals', 'all_genes'))
        ensembl_protein_coding_intervals_dir = my.makedir(os.path.join(vax_version_dir_link,'ensembl_intervals', 'protein_coding_genes'))
        ensembl_PERL5LIB = '{PLUGINS_DIR}:{ENSEMBL_DIR}/ensembl/modules:{ENSEMBL_DIR}/ensembl-variation/modules:{ENSEMBL_DIR}/ensembl-funcgen/modules:{ENSEMBL_DIR}/ensembl-compara/modules:{BIOPERL_DIR}'.format(PLUGINS_DIR=plugins_dir_link, ENSEMBL_DIR=ensembl_api_dir_link, BIOPERL_DIR=bioperl_live_link)
        ensembl_xref_dir = my.makedir(os.path.join(vax_version_dir_link,'ensembl_xref'))
        fathmm_dir = my.makedir(os.path.join(scratch_version_dir,'fathmm'))
        fathmm_plugin_dir = my.makedir(os.path.join(plugins_version_dir,'fathmm'))
        gatk_dir = my.makedir(os.path.join(apps_dir,'gatk'))
        hpa_dir = my.makedir(os.path.join(scratch_version_dir,'hpa'))
        job_dir = my.makedir(os.path.join(scratch_version_dir, 'jobs'))
        kegg_dir = my.makedir(os.path.join(scratch_version_dir,'kegg'))
        #TODO: 1000 genomes
        kg_dir = my.makedir(os.path.join(scratch_version_dir,'1000genomes'))
        metabolome_dir = my.makedir(os.path.join(scratch_version_dir,'metabolome'))
        mgi_dir = my.makedir(os.path.join(scratch_version_dir,'mgi'))
        mitocarta_dir = my.makedir(os.path.join(scratch_version_dir,'mitocarta'))
        nhlbi_dir = my.makedir(os.path.join(scratch_version_dir,'nhlbi'))
        omim_dir = my.makedir(os.path.join(scratch_version_dir,'omim'))
        refgene_dir = my.makedir(os.path.join(scratch_version_dir,'refgene'))
        refseq_dir = my.makedir(os.path.join(scratch_version_dir,'refseq'))
        rgd_dir = my.makedir(os.path.join(scratch_version_dir,'rgd'))
        tmp = my.makedir(os.path.join(scratch_version_dir, 'tmp'))
        uniprot_dir = my.makedir(os.path.join(scratch_version_dir,'uniprot'))
        
        vep_cache_dir = scratch_version_dir_link
        with cd(vax_dir):
            try:
                os.remove('vep_cache')
            except OSError:
                pass
            os.symlink(os.path.basename(scratch_version_dir_link), 'vep_cache')
        
        vep_dir = my.makedir(os.path.join(vax_version_dir_link, 'VEP'))
        
        vep_ini = os.path.join(vax_version_dir_link, 'vep.ini')
        vep_ini_link = os.path.join(os.path.basename(vax_version_dir_link), 'vep.ini')
        with cd(vax_dir):
            try:
                os.remove('vep.ini')
            except OSError:
                pass
            os.symlink(vep_ini_link, 'vep.ini')

        vep_pl = os.path.join(ensembl_api_dir_link, config.get('ensembl', 'ENSEMBL_API_VEP_PL_PATH'))
        with cd(vax_dir):
            try:
                os.remove('vep.pl')
            except OSError:
                pass
            with cd(vax_dir):
                os.symlink(vep_pl, 'vep.pl')
        
        vep_plugins_git_dir = my.makedir(os.path.join(vax_version_dir_link,'vep_plugins_git'))

        self.print_log(['Database will be installed with these parameters:',
                        'user = {}'.format(user),
                        'host = {}'.format(host),
                        'port = {}'.format(port),
                        'Ensembl version = {}'.format(ensembl_version),
                        'vax version = {}'.format(vax_version),
                        'vax database name = {}'.format(database),
                        'vax runtime user = {}'.format(vax_user),
                        'vax runtime password = {}'.format(vax_password),
                        'ensembl runtime user = {}'.format(ensembl_user),
                        'ensembl runtime password = {}'.format(ensembl_password),
                        'hgmd runtime user = {}'.format(hgmd_user),
                        'hgmd runtime password = {}'.format(hgmd_password),
                        ])

################################################################
#       test MySQL connection and create VAX database          #
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
                raise VaxInstallerError('Database connection failed. {}'.format(msg))
            else:
                if not mysql_version:
                    raise VaxInstallerError('Database connected but SELECT VERSION() query failed.')
                self.print_log('Database connected. MySQL version: {}'.format(mysql_version))
                query = "SELECT SCHEMA_NAME FROM INFORMATION_SCHEMA.SCHEMATA WHERE SCHEMA_NAME = '{}'".format(database)
                self.mysql(cursor, query)
                schema = cursor.fetchone()
                self.print_log('Database `{}` {}.'.format(database, "doesn't exist" if schema == None else "exists"))
                
                #get existing database's versions table, if this version's database already exists
                try:
                    versions = my.get_key_value_dict(host, user, password, port, database, table='versions', key_col='module', value_col='version',)
                except:
                    versions = {}
                versions['ensembl'] = str(ensembl_version)
                versions['vax'] = str(vax_version)
                
                #create new database with empty versions table
                if schema == None or overwrite_vax_database:
                    self.print_log('Creating new database `{}`. Any existing data and procedures will be deleted.'.format(database))
                    query = ['DROP DATABASE IF EXISTS `{}`'.format(database),
                             'CREATE DATABASE `{}`'.format(database),
                             """CREATE TABLE `{}`.`versions` (
`module` VARCHAR(255) NOT NULL,
`version` VARCHAR(255) NULL,
PRIMARY KEY (`module`))""".format(database,),]
                    self.mysql(cursor, query)
                    self.print_log('Created new database `{}`.'.format(database))
                else:
                    self.print_log('Retaining database `{}`, but overwriting tables and procedures for selected modules.'.format(database))
                query = "GRANT SELECT, EXECUTE ON {}.* TO '{}'@'%' IDENTIFIED BY '{}'".format(database, vax_user, vax_password)
                self.mysql(cursor, query)

        self.print_log('Purging MySQL logs.')
        with MySQLdb.connect(host,user,password,database,port=port) as cursor:
            query = "PURGE BINARY LOGS BEFORE '{}'".format(strftime("%Y-%m-%d %H:%M:%S", localtime()))
            self.mysql(cursor, query)

        self.print_log('Database connection OK')


################################################################
#                   Ensembl MySQL database                     #
################################################################

        if 'all' in install or 'ensembl_database' in install:
            done_file = os.path.join(ensembl_db_dir,'.schedule_ensembl_db_download.done')
            if my.file_exists(done_file):
                self.print_log('Download Ensembl MySQL databases cluster jobs already scheduled. To reschedule "rm {}"'.format(done_file))
            else:
                ensembl_host = config.get('ensembl', 'ENSEMBL_HOST')
                ensembl_ftp_dir = config.get('ensembl', 'ENSEMBL_FTP_DIR_PATTERN').format(ensembl_version)
                ensembl_ftp_pattern = config.get('ensembl', 'ENSEMBL_FTP_PATTERN')
                ensembl_species_ftp_pattern = config.get('ensembl', 'ENSEMBL_SPECIES_FTP_PATTERN')
                ftp_pattern = '|'.join([ensembl_ftp_pattern] + [ensembl_species_ftp_pattern.format(s) for s in ensembl_species])
                # ensembl_rsync_url =  config.get('ensembl', 'ENSEMBL_RSYNC_URL_CUTTENT')
                ensembl_rsync_url =  config.get('ensembl', 'ENSEMBL_RSYNC_URL_PATTERN').format(ensembl_version)
                #filter out unnecessary files/databases
                ensembl_blacklist_patterns =  config.get('ensembl', 'ENSEMBL_BLACKLIST_PATTERNS').split()
                rx = re.compile('({})'.format(')|('.join(ensembl_blacklist_patterns)))
                ftp_files = [f for f in my.list_ftp_files(ensembl_host, ensembl_ftp_dir, ftp_pattern) if not rx.match(f)]
                #rsync a few very large databases separately, the rest in a single rsync to avoid overloading the server
                ensembl_rsync_parallel_patterns =  config.get('ensembl', 'ENSEMBL_RSYNC_PARALLEL_PATTERNS').split()
                ensembl_species_rsync_parallel_pattern =  config.get('ensembl', 'ENSEMBL_SPECIES_RSYNC_PARALLEL_PATTERN')
                rsync_parallel_pattern = '({})'.format(')|('.join(ensembl_rsync_parallel_patterns + [ensembl_species_rsync_parallel_pattern.format(s) for s in ensembl_species]))
                rx = re.compile(rsync_parallel_pattern)
                #list of rsync batches; batches will be scheduled in parallel jobs]
                rsyncs = [[os.path.join(ensembl_rsync_url,f)] for f in ftp_files if rx.match(f)]
                rsyncs.append([os.path.join(ensembl_rsync_url,f) for f in ftp_files if not rx.match(f)]) 
                if delete_ensembl_version:
                    self.print_log('Deleting older Ensembl databases.')
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
                    cmd = 'rsync -av "{}" "{}/"'.format('" "'.join(batch), ensembl_db_dir)
                    job_name = 'rsync_ensembl_{}_db'.format(ensembl_version)
                    job = my.run_job(cmd, job_name, job_dir)
                    ensembl_job_id = job.jobId
                    self.ensembl_jobs.append(job.jobId)
                    self.all_jobs.append(job.jobId)
                    self.print_log('Scheduled job {}: {}.'.format(job.jobId, cmd))
                    for file in batch:
                        ensembl_subdir = os.path.join(ensembl_db_dir, os.path.basename(file))
                        cmd = 'gunzip -rf {}/*'.format(ensembl_subdir)
                        job_name = 'gunzip_ensembl_{}_db'.format(ensembl_version)
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
                self.print_log('Scheduling job to create new Ensembl database version {} and import data.'.format(ensembl_version))
                vax_ensembl_database_installer_py = os.path.join(os.path.dirname(__file__), 'vax_ensembl_database_installer.py')
                cmd= '{} {} --ensembl_db_dir {} --host {} --port {} --user {} --password {} --ensembl_user {} --ensembl_password {} '.format(
                    config.get('DEFAULT','python'), vax_ensembl_database_installer_py, ensembl_db_dir, host, port, user, password, ensembl_user, ensembl_password)
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
#         Ensembl Perl API & Variant Effect Predictor          #
################################################################

        if 'all' in install or 'ensembl_api' in install:
            done_file = os.path.join(ensembl_api_dir,'.ensembl_api.done')
            if my.file_exists(done_file):
                self.print_log('Ensembl Perl API already installed. To reinstall "rm {}"'.format(done_file))
            else:
                self.print_log('Installing Ensembl Perl API.')
                #remove any pre-existing directory
                if os.path.isdir(ensembl_api_dir):
                    shutil.rmtree(ensembl_api_dir)
                #download ensembl git-tools
                git_tools_url = config.get('ensembl','GIT_TOOLS_URL')
                git_tools_dir =  os.path.join(ensembl_api_dir, 'ensembl-git-tools')
                self.print_log('Cloning {} in {}.'.format(git_tools_url, os.path.join(ensembl_api_dir, 'ensembl-git-tools')))
                git.clone(git_tools_url, git_tools_dir)
                #use ensembl-git to download Ensembl API
                self.print_log('Cloning Ensembl Perl API in {}.'.format(ensembl_api_dir))
                git_ensembl = sh.Command(os.path.join(ensembl_api_dir, 'ensembl-git-tools', 'bin', 'git-ensembl'))
                #tools gets api and ensembl-tools
                #use pull instead of clone for preexisting directories
                git_ensembl('--dir', ensembl_api_dir, '--clone', 'tools')
                if not my.file_exists(vep_pl):
                    self.print_log('WARNING: Failed to install Variant Effect Predictor perl script at {}'.format(vep_pl))
                    
                with open(done_file, 'w'):
                    pass
                self.print_log('Installed Ensembl Perl API.')


################################################################
#                           VEP Cache                          #
################################################################

#TODO: #http://www.ensembl.org/info/docs/tools/vep/script/vep_other.html
#
#Variant Effect Predictor  Other information
#
#With the arrival of GRCh38, Ensembl now supports two different assembly versions for the human genome while users
# transition from GRCh37. We provide a VEP cache download on the latest software version (76) for both assembly versions.
#
#The VEP installer will install and set up the correct cache and FASTA file for your assembly of interest. If using
# the --AUTO functionality to install without prompts, remember to add the assembly version required using
# e.g. "--ASSEMBLY GRCh37". It is also possible to have concurrent installations of caches from both assemblies;
# just use the --assembly to select the correct one when you run the VEP script.
#
#Once you have installed the relevant cache and FASTA file, you are then able to use the VEP as normal.
# For those using GRCh37 and requiring database access in addition to the cache (for example, to look up variant
# identifiers using --format id, see cache limitations), the script will warn you that you must change the database
# port in order to connect to the correct database:
# ERROR: Cache assembly version (GRCh37) and database or selected assembly version (GRCh38) do not match
# If using human GRCh37 add "--port 3337" to use the GRCh37 database, or --offline to avoid database connection entirely


#TODO: provide for non-human species
        if 'all' in install or 'vep_cache' in install:
            self.print_log('Installing VEP Cache.')
            url = config.get('vep', 'VEP_DOWNLOAD_URL_FORMAT_STRING').format(ensembl_version)
            zipfile = os.path.join(tmp, config.get('vep', 'ENSEMBL_TOOLS_ZIPFILE_FORMAT_STRING').format(ensembl_version))
            unzipdir = vep_dir
            self.getzip(url, zipfile, unzipdir)
            
            #run VEP installer to download VEP cache
            done_file = os.path.join(vep_dir, '.vep.install.done')
            if my.file_exists(done_file):
                self.print_log('VEP cache already installed. To reinstall "rm {}"'.format(done_file))
            else:
                vep_installer_pl = os.path.join(vep_dir, config.get('vep', 'VEP_INSTALLER_PL_FORMAT_STRING').format(ensembl_version))
                self.print_log('Running VEP installer.')
                vip = perl(vep_installer_pl, '--DESTDIR', vep_dir, '--CACHEDIR', vep_cache_dir,'--SPECIES', config.get('vep', 'VEP_SPECIES'), '--AUTO', config.get('vep', 'VEP_AUTO'))
                
                with open(done_file, 'w'):
                    pass
                self.print_log('Installed VEP Cache.')

################################################################
#                          VEP Plugins                         #
################################################################
        
        if 'all' in install or 'vep_plugins' in install:
            done_file = os.path.join(plugins_version_dir,'.vep_plugins.done')
            if my.file_exists(done_file):
                self.print_log('Ensembl VEP plugins already installed. To reinstall "rm {}"'.format(done_file))
            else:
                self.print_log('Installing VEP Plugins.')
                #remove any pre-existing directory
                if os.path.isdir(vep_plugins_git_dir):
                    shutil.rmtree(vep_plugins_git_dir)
                #download VEP plugins
                vep_plugin_git_url = config.get('vep','VEP_PLUGINS_GIT_URL')
                self.print_log('Cloning {} into {}.'.format(vep_plugin_git_url, vep_plugins_git_dir))
                git.clone(vep_plugin_git_url, vep_plugins_git_dir)
                #use ensembl-git to download Ensembl API
                self.print_log('VEP Plugins cloned.')
                self.print_log('Copying VEP Plugins to {}'.format(plugins_version_dir))
                vep_plugins_glob = vep_plugins_git_dir+'/*'
                sh.cp('-a', sh.glob(vep_plugins_glob), plugins_version_dir)
    
                #fix Condel plugin local path
                self.print_log('Fixing Condel plugin local path')
                condel_conf = os.path.join(plugins_version_dir, 'config', 'Condel', 'config', 'condel_SP.conf')
                for line in fileinput.input(condel_conf, inplace=1):
                    line=line.strip()
                    if(line == "condel.dir='path/to/config/Condel/'"):
                        print("condel.dir='{}/config/Condel/'".format(plugins_version_dir))
                    else:
                        print(line)
                        
                with open(done_file, 'w'):
                    pass
                self.print_log('Installed VEP Plugins.')

################################################################
#                           BioPerl                            #
################################################################

        if 'all' in install or 'bioperl' in install:
            done_file = os.path.join(bioperl_live_dir,'.bioperl-live.done')
            if my.file_exists(done_file):
                self.print_log('Bioperl already installed. To reinstall "rm {}"'.format(done_file))
            else:
                self.print_log('Installing BioPerl.')
                if os.path.isdir(bioperl_live_dir):
                    shutil.rmtree(bioperl_live_dir)
                #git bioperl-live
                self.print_log('Cloning BioPerl.')
                bioperl_url = config.get('bioperl','BIOPERL_URL')
                git.clone(bioperl_url, bioperl_live_dir)
                self.print_log('BioPerl cloned.')
                #determine version if not specified as a run parameter
                if not bioperl_version:
                    default_bioperl_version = None
                    try:
                        default_bioperl_version = config.get('bioperl', 'DEFAULT_BIOPERL_VERSION')
                    except:
                        pass
                    self.print_log('Getting list of BioPerl versions.')
                    with cd(bioperl_live_dir):
                        bioperl_versions = git.tag().strip().split('\n')
                    if not bioperl_versions:
                        self.print_log('No bioperl version git tags in {}. Check the bioperl online documentation.'.format(bioperl_live_dir))
                    else:
                        bioperl_version = my.select_item(bioperl_versions, default_bioperl_version, 'BioPerl versions', 'of BioPerl version to install')
                if not bioperl_version:
                    self.print_log('WARNING: No Bioperl version selected. Rerun or checkout manually.')
                else:
                    #checkout selected version
                    versions['bioperl'] = str(bioperl_version)
                    self.print_log('Checking out BioPerl {} '.format(bioperl_version))
                    with cd(bioperl_live_dir):
                        git.checkout(bioperl_version)
                    self.print_log('BioPerl version {} checked out.'.format(bioperl_version))
                    #remember version as default
                    self.print_log('setting BioPerl version {} as default for next install.'.format(bioperl_version))
                    config.set('bioperl', 'DEFAULT_BIOPERL_VERSION', bioperl_version)
                    
                    with open(done_file, 'w'):
                        pass
                    self.print_log('Installed BioPerl.')

################################################################
#                      Ensembl intervals                       #
################################################################

#requires Ensembl database and Ensembl Perl API to be installed and configured in PERL5LIB

        #TODO: provide for non-human species

        if 'all' in install or 'ensembl_intervals' in install:
            done_file = os.path.join(ensembl_intervals_dir,'.schedule_ensembl_intervals.done')
            if my.file_exists(done_file):
                self.print_log('Ensembl intervals cluster job already scheduled. To reschedule "rm {}"'.format(done_file))
            else:
                self.print_log('Scheduling job to create and install ensembl_intervals table.')
                vax_ensembl_genomic_regions_pl = os.path.join(os.path.dirname(__file__), 'vax_ensembl_genomic_regions.pl')
                vax_ensembl_intervals_installer_py = os.path.join(os.path.dirname(__file__), 'vax_ensembl_intervals_installer.py')
                #make all_genes and protein coding intervals files; install all_genes as a mysql table
                #protein_coding
                cmd = '{} {} -host {} -port {} -user {} -password {} -output_dir {} --all_genes --include_biotype protein_coding --five_prime_UTR --three_prime_UTR --cis_splice_site --five_prime_splice_region_bases 4 --three_prime_splice_region_bases 13'.format(
                    config.get('DEFAULT', 'perl'), vax_ensembl_genomic_regions_pl, host, port, ensembl_user, ensembl_password, ensembl_protein_coding_intervals_dir)
                job_name = 'intervals_ensembl_{}'.format(ensembl_version)
                job = my.run_job(cmd, job_name, job_dir, hold_jid=self.ensembl_jobs, memory = '24G')
                ensembl_intervals_job_id = job.jobId
                self.all_jobs.append(job.jobId)
                self.print_log('Scheduled job {}: {}.'.format(job.jobId, cmd))
                #all_genes
                cmd = '{} {} -host {} -port {} -user {} -password {} -output_dir {} --all_genes --include_putative_genes --include_novel_genes --include_pseudogenes --five_prime_UTR --three_prime_UTR --cis_splice_site --five_prime_splice_region_bases 4 --three_prime_splice_region_bases 13'.format(
                    config.get('DEFAULT', 'perl'), vax_ensembl_genomic_regions_pl, host, port, ensembl_user, ensembl_password, ensembl_intervals_dir)
                job_name = 'intervals_ensembl_{}'.format(ensembl_version)
                job = my.run_job(cmd, job_name, job_dir, hold_jid=self.ensembl_jobs, memory = '24G')
                ensembl_intervals_job_id = job.jobId
                self.all_jobs.append(job.jobId)
                self.print_log('Scheduled job {}: {}.'.format(job.jobId, cmd))
                self.print_log('Scheduling job to create and import ensembl_intervals_{} tables.'.format(ensembl_version))
                cmd= '{} {} --ensembl_intervals_dir {} --host {} --port {} --user {} --password {} --database {}'.format(
                    config.get('DEFAULT','python'), vax_ensembl_intervals_installer_py, ensembl_intervals_dir, host, port, user, password, database)
                job_name = 'intervals_import_ensembl_{}'.format(ensembl_version)
                job = my.run_job(cmd, job_name, job_dir, hold_jid=ensembl_intervals_job_id)
                ensembl_intervals_import_job_id = job.jobId
                self.all_jobs.append(job.jobId)
                self.print_log('Scheduled job {}: {}.'.format(job.jobId, cmd))

################################################################
#                         Ensembl xref                         #
################################################################

#requires Ensembl database and Ensembl Perl API to be installed and configured in PERL5LIB
        if 'all' in install or 'ensembl_xref' in install:
            for s in ensembl_species:
                done_file = os.path.join(ensembl_xref_dir,'.schedule_ensembl_xref_{}_{}.done'.format(ensembl_version, s))
                if my.file_exists(done_file):
                    self.print_log('Ensembl xref cluster job already scheduled. To reschedule "rm {}"'.format(done_file))
                else:
                    self.print_log('Scheduling job to create and install ensembl_xref_{} table.'.format(s))
                    vax_ensembl_xref2db_pl = os.path.join(os.path.dirname(__file__), 'vax_ensembl_xref2db.pl')
                    vax_ensembl_xref_installer_py = os.path.join(os.path.dirname(__file__), 'vax_ensembl_xref_installer.py')
                    ensembl_xref_with_dups_file = os.path.join(ensembl_xref_dir, 'ensembl_xref_{}_{}_with_dups.txt'.format(ensembl_version, s))
                    cmd = '{} {} -host {} -port {} -user {} -pass {} -species {} -output {}'.format(
                        config.get('DEFAULT','perl'), vax_ensembl_xref2db_pl,
                        host, port, ensembl_user, ensembl_password, s, ensembl_xref_with_dups_file)
                    job_name = 'xref_ensembl_{}_{}'.format(ensembl_version, s)
                    job = my.run_job(cmd, job_name, job_dir, hold_jid=self.ensembl_jobs)
                    ensembl_xref_job_id = job.jobId
                    self.all_jobs.append(job.jobId)
                    self.print_log('Scheduled job {}: {}.'.format(job.jobId, cmd))
                    self.print_log('Scheduling job to create and import ensembl_xref_{}_{} table.'.format(ensembl_version, s))
                    cmd= '{} {} --ensembl_xref_with_dups_file {} --host {} --port {} --user {} --password {} --database {}'.format(
                        config.get('DEFAULT','python'), vax_ensembl_xref_installer_py, ensembl_xref_with_dups_file, host, port, user, password, database)
                    job_name = 'xref_import_ensembl_{}_{}'.format(ensembl_version, s)
                    job = my.run_job(cmd, job_name, job_dir, hold_jid=ensembl_xref_job_id)
                    ensembl_xref_import_job_id = job.jobId
                    self.all_jobs.append(job.jobId)
                    self.print_log('Scheduled job {}: {}.'.format(job.jobId, cmd))

################################################################
#                             CADD                             #
################################################################

        if 'all' in install or 'cadd' in install:
            done_file = my.done_file(os.path.join(cadd_dir,'cadd_schedule_download'))
            if my.file_exists(done_file):
                self.print_log('CADD already installed. To reinstall "rm {}"'.format(done_file))
            else:
                self.print_log('Installing CADD.')
                default_cadd_snp_file_url = None
                if not cadd_snp_file_url:
                    try:
                        default_cadd_snp_file_url = config.get('cadd', 'DEFAULT_CADD_SNP_FILE_URL')
                    except:
                        pass
                    #determine version if not specified as a run parameter
                    self.print_log('Getting list of CADD SNP file urls.')
                    cadd_info_url = config.get('cadd', 'CADD_INFO_URL')
                    cadd_snp_file_pattern = config.get('cadd', 'CADD_SNP_FILE_PATTERN')
                    cadd_snp_file_urls = list(set(my.list_hrefs(cadd_info_url, cadd_snp_file_pattern)))
                    if not cadd_snp_file_urls:
                        self.print_log('No urls at {} match the pattern {}. Check the web site and CADD_SNP_FILE_PATTERN in vax.cfg.'.format(cadd_info_url, cadd_snp_file_pattern))
                    else:
                        cadd_snp_file_url = my.select_item(cadd_snp_file_urls, default_cadd_snp_file_url, 'CADD SNP files', 'CADD SNP file to install')
                if not cadd_snp_file_url:
                    self.print_log('WARNING: No CADD SNP file selected. Rerun or checkout manually.')
                else:
                    #download selected version
                    #version in the bottom-level directory that contains the file, e.g., v1.0
                    versions['cadd'] = os.path.basename(os.path.dirname(cadd_snp_file_url))
                    #versions['cadd'] = 'CADD accessed {}'.format(strftime("%Y-%m-%d", localtime()))
                    self.print_log('Scheduling job to download CADD files; ~30 minute run time.')
                    cadd_snp_file = os.path.join(cadd_dir, os.path.basename(cadd_snp_file_url))
                    vax_CADD_downloader_py = os.path.join(os.path.dirname(__file__), 'vax_CADD_downloader.py')
                    cmd= '{} {} --snp_file_url {} --dir {}'.format(
                        config.get('DEFAULT','python'), vax_CADD_downloader_py, cadd_snp_file_url, cadd_dir,)
                    job_name = 'CADD_download'
                    job = my.run_job(cmd, job_name, job_dir)
                    CADD_download_job_id = job.jobId
                    self.all_jobs.append(job.jobId)
                    self.print_log('Scheduled job {}: {}.'.format(job.jobId, cmd))
                    #remember version as default
                    self.print_log('setting CADD SNP file url {} as default for next install.'.format(cadd_snp_file_url))
                    config.set('cadd', 'DEFAULT_CADD_SNP_FILE_URL', cadd_snp_file_url)
                    
                    with open(done_file, 'w'):
                        pass
                    self.print_log('Scheduled job to download CADD files.')

################################################################
#                            dbNSFP                            #
################################################################

        if 'all' in install or 'dbnsfp' in install:
            done_file = my.done_file(os.path.join(dbnsfp_dir,'dbnsfp_download_scheduled'))
            dbnsfp_gz = os.path.join(dbnsfp_dir, 'dbNSFP.gz')
            if my.file_exists(done_file):
                self.print_log('dbNSFP already scheduled for installation. To reinstall "rm {}"'.format(done_file))
            else:
                self.print_log('Installing dbNSFP.')
                default_dbnsfp_file = None
                if not dbnsfp_file:
                    try:
                        default_dbnsfp_file = config.get('dbnsfp', 'DEFAULT_DBNSFP_FILE')
                    except:
                        pass
                    #determine version if not specified as a run parameter
                    self.print_log('Getting list of dbNSFP file urls.')
                    dbnsfp_url = config.get('dbnsfp', 'DBNSFP_URL')
                    dbnsfp_snp_file_pattern = config.get('dbnsfp', 'DBNSFP_FILE_PATTERN')
                    dbnsfp_files = my.list_hrefs(dbnsfp_url, dbnsfp_snp_file_pattern)
                    if not dbnsfp_files:
                        self.print_log('No file names at {} match the pattern {}. Check the web site and DBNSFP_FILE_PATTERN in vax.cfg.'.format(dbnsfp_url, dbnsfp_snp_file_pattern))
                    else:
                        dbnsfp_file = my.select_item(dbnsfp_files, default_dbnsfp_file, 'dbNSFP files', 'dbNSFP file to install')
                if not dbnsfp_file:
                    self.print_log('WARNING: No dbNSFP file selected. Rerun or checkout manually.')
                else:
                    #download selected version
                    #version is file name without extension
                    versions['dbnsfp'] = my.r_strip(my.l_strip(dbnsfp_file, 'dbNSFP'), '.zip')
                    #versions['dbnsfp'] = 'dbNSFP accessed {}'.format(strftime("%Y-%m-%d", localtime()))
                    dbnsfp_zip_url = os.path.join(dbnsfp_url, dbnsfp_file)
                    dbnsfp_zip_file = os.path.join(dbnsfp_dir, dbnsfp_file)
                    
                    self.print_log('Scheduling job to download, merge, and index dbNSFP files; ~5 hour run time.')
                    vax_dbNSFP_downloader_py = os.path.join(os.path.dirname(__file__), 'vax_dbNSFP_downloader.py')
                    cmd= '{} {} --dbnsfp_zip_url {} --dir {}'.format(
                        config.get('DEFAULT','python'), vax_dbNSFP_downloader_py, dbnsfp_zip_url, dbnsfp_dir,)
                    job_name = 'dbNSFP_download'
                    job = my.run_job(cmd, job_name, job_dir)
                    dbNSFP_download_job_id = job.jobId
                    self.all_jobs.append(job.jobId)
                    self.print_log('Scheduled job {}: {}.'.format(job.jobId, cmd))
                    
                    
                    self.print_log('Setting dbNSFP file {} as default for next install.'.format(dbnsfp_file))
                    config.set('dbnsfp', 'DEFAULT_DBNSFP_FILE', dbnsfp_file)
                    
                    with open(done_file, 'w'):
                        pass
                    self.print_log('Scheduled job to download, merge, and index dbNSFP files.')

################################################################
#                            dbSNP                             #
################################################################

        if 'all' in install or 'dbsnp' in install:
            done_file = my.done_file(os.path.join(dbsnp_dir,'dbsnp_schedule_download'))
            if my.file_exists(done_file):
                self.print_log('dbSNP already installed. To reinstall "rm {}"'.format(done_file))
            else:
                self.print_log('Installing dbSNP. This will take >2 days.')
                
                dbsnp_version = config.get('dbsnp', 'DBSNP_VERSION')
                dbsnp_ftp = config.get('dbsnp', 'DBSNP_FTP')
                dbsnp_ftp_dir = config.get('dbsnp', 'DBSNP_FTP_DIR')
                dbsnp_file_pattern = config.get('dbsnp', 'DBSNP_FILE_PATTERN')
                dbsnp_file_rx = re.compile(dbsnp_file_pattern, re.IGNORECASE)
                dbsnp_clinvar_file_pattern = config.get('dbsnp', 'DBSNP_CLINVAR_FILE_PATTERN')
                dbsnp_clinvar_file_rx = re.compile(dbsnp_clinvar_file_pattern, re.IGNORECASE)
                dbsnp_common_no_known_medical_impact_file_pattern = config.get('dbsnp', 'DBSNP_COMMON_NO_KNOWN_MEDICAL_IMPACT_FILE_PATTERN')
                dbsnp_common_no_known_medical_impact_file_rx = re.compile(dbsnp_common_no_known_medical_impact_file_pattern, re.IGNORECASE)
    
                #download dbSNP directory
                dbsnp_files = self.download_ftp_dir(dbsnp_ftp, dbsnp_ftp_dir, dbsnp_dir)
                input_files = []
                m = [f for f in dbsnp_files if dbsnp_file_rx.match(os.path.basename(f))]
                dbsnp_file = m[0] if m else None
                #input_files = [(path,stored_procedure), ...]
                if dbsnp_file:
                    input_files.append((dbsnp_file, 'get_dbsnp'))
                    dbsnp_meta = my.parse_vcf_meta_info(dbsnp_file)
                    dbsnp_version = dbsnp_meta.get('dbSNP_BUILD_ID', '?')
                    versions['dbsnp'] = dbsnp_version
                    self.print_log('dbSNP version is {}'.format(dbsnp_version))
                else:
                    self.print_log('WARNING: could not find a file matching "{}" on the dbSNP ftp site.'.format(dbsnp_file_pattern))
                    
                m = [f for f in dbsnp_files if dbsnp_clinvar_file_rx.match(os.path.basename(f))]
                dbsnp_clinvar_file = m[0] if m else None
                if dbsnp_clinvar_file:
                    input_files.append((dbsnp_clinvar_file, 'get_dbsnp_clinvar'))
                else:
                    self.print_log('WARNING: could not find a file matching "{}" on the dbSNP ftp site.'.format(dbsnp_clinvar_file_pattern))
                m = [f for f in dbsnp_files if dbsnp_common_no_known_medical_impact_file_rx.match(os.path.basename(f))]
                dbsnp_common_no_known_medical_impact_file = m[0] if m else None
                if dbsnp_common_no_known_medical_impact_file:
                    input_files.append((dbsnp_common_no_known_medical_impact_file, 'get_dbsnp_common_no_known_medical_impact'))
                else:
                    self.print_log('WARNING: could not find a file matching "{}" on the dbSNP ftp site.'.format(dbsnp_common_no_known_medical_impact_file_pattern))
                    
                self.print_log('Scheduling jobs to install dbSNP mysql tables. The larger dbSNP files can take up to three days to process')
                vax_dbsnp_installer_py = os.path.join(os.path.dirname(__file__), 'vax_dbsnp_installer.py')
                for i in input_files:
                    this_file, stored_procedure = i
                    table = 'dbsnp_'+dbsnp_version+'_'+my.r_strip(os.path.basename(this_file), '.vcf.gz')
                    text_file = os.path.join(os.path.dirname(this_file), table+'.txt')
                    text_file_done = my.done_file(text_file)
                    if my.file_exists(text_file_done):
                        self.print_log('{} already imported. To reimport "rm {}"'.format(text_file, text_file_done))
                    else:
                        self.print_log('Scheduling job to create and import `{}`.`{}` table from {}.'.format(database, table, this_file))
                        cmd= '{} {} --dbsnp_file {} --dbsnp_version {} --database {} --proc {} --host {} --port {} --user {} --password {}'.format(
                            config.get('DEFAULT','python'), vax_dbsnp_installer_py, this_file, dbsnp_version, database, stored_procedure, host, port, user, password)
                        job_name = '{}_import'.format(table)
                        job = my.run_job(cmd, job_name, job_dir)
                        dbsnp_import_job_id = job.jobId
                        self.all_jobs.append(job.jobId)
                        self.print_log('Scheduled job {}: {}.'.format(job.jobId, cmd))
    
                with open(done_file, 'w'):
                    pass
                self.print_log('Installed dbSNP.')

################################################################
#                            faidx                             #
################################################################

        if 'all' in install or 'faidx' in install:
            # required by RefGene
            faidx_dir = config.get('faidx','FAIDX_DIR')
            faidx_file = config.get('faidx','FAIDX_FILE')
            faidx_decoy_file = config.get('faidx','FAIDX_DECOY_FILE')
            done_file = os.path.join(faidx_dir,'.faidx.done')
            if my.file_exists(done_file):
                self.print_log('Faidx installed. To reinstall "rm {}"'.format(done_file))
            else:
                self.print_log('Installing faidx.')
                for text_file in (faidx_file, faidx_decoy_file):
                    table = my.r_strip(os.path.basename(text_file), '.txt')
                    sql_file = text_file+'.mysql'
                    sql_done_file = os.path.join(os.path.dirname(sql_file), '.'+os.path.basename(sql_file)+'.done')
                    
                    if my.file_exists(sql_done_file):
                        self.print_log('{} already created. To redo "rm {}"'.format(sql_file, sql_done_file))
                    else:
                        self.print_log('Making CREATE TABLE `{}` script.'.format(table))
                        sql_columns.run(text_file, database=database, schema=database, table=table, primary_key=['genome', 'contig'], header_line=1)
                        with open(sql_done_file, 'w'):
                            pass
                        self.print_log('Made CREATE TABLE `{}` script.'.format(table))
                        
                    self.print_log('Creating `{}` table.'.format(table))
                    with open(sql_file, 'r') as sql:
                        sh.mysql('-h', host, '-P', port, '-u', user, '-p{}'.format(password), database, _in=sql)
                        
                    self.print_log('Importing data to `{}`.`{}` table.'.format(database, table))
                    sh.mysqlimport('-h', host, '-P', port, '-u', user, '-p{}'.format(password), '--delete', '--local', '--ignore-lines=1', database, text_file)
                    self.print_log('Imported data to `{}`.`{}` table.'.format(database, table))
                    
                with open(done_file, 'w'):
                    pass
                self.print_log('Installed faidex.')

################################################################
#                           FATHMM                             #
################################################################

        if 'all' in install or 'fathmm' in install:
            done_file = os.path.join(fathmm_dir,'.fathmm.done')
            if my.file_exists(done_file):
                self.print_log('FATHMM already installed. To reinstall "rm {}"'.format(done_file))
            else:
                self.print_log('Installing FATHMM.')
                default_fathmm_file = None
                if not fathmm_file:
                    try:
                        default_fathmm_file = config.get('fathmm', 'DEFAULT_FATHMM_FILE')
                    except:
                        pass
                    #determine version if not specified as a run parameter
                    self.print_log('getting list of FATHMM file names.')
                    fathmm_data_host = config.get('fathmm', 'FATHMM_DATA_FTP_HOST')
                    fathmm_data_dir = config.get('fathmm', 'FATHMM_DATA_FTP_DIR')
                    fathmm_file_pattern = config.get('fathmm', 'FATHMM_FILE_PATTERN')
                    fathmm_files = my.list_ftp_files(fathmm_data_host, fathmm_data_dir, fathmm_file_pattern)
                    if not fathmm_files:
                        self.print_log('No file names at {}/{} match the pattern {}. Check the ftp site and FATHMM_FILE_PATTERN in vax.cfg.'.format(fathmm_data_host, fathmm_data_dir, fathmm_file_pattern))
                    else:
                        fathmm_file = my.select_item(fathmm_files, default_fathmm_file, 'FATHMM files', 'FATHMM file to install')
                if not fathmm_file:
                    self.print_log('WARNING: No FATHMM file selected. Rerun or checkout manually.')
                else:
                    #download selected database version
                    versions['fathmm'] = my.r_strip(my.l_strip(fathmm_file, 'fathmm.'), '.SQL.gz')
                    fathmm_data_url = config.get('fathmm', 'FATHMM_DATA_URL_FORMAT_STRING').format(fathmm_file)
                    fathmm_data_file = os.path.join(fathmm_dir, fathmm_file)
                    self.print_log('Downloading current FATHMM SQL database from {} into {}. >8 Gb (unzipped >21 Gb). This will take ~10 minutes'.format(fathmm_data_url, fathmm_data_file))
                    urlretrieve(fathmm_data_url, fathmm_data_file)
                    self.print_log('FATHMM SQL database downloaded.')
                    #remember version as default
                    self.print_log('setting FATHMM file {} as default for next install.'.format(fathmm_file))
                    config.set('fathmm', 'DEFAULT_FATHMM_FILE', fathmm_file)
                    #create FATHMM MySQL database
                    fathmm_database = config.get('fathmm','FATHMM_DATABASE')
                    self.print_log('Creating FATHMM database `{}`. Any preexisting data will be deleted.'.format(fathmm_database))
                    with MySQLdb.connect(host,user,password,port=port) as cursor:
                        query = 'DROP SCHEMA IF EXISTS `{}`'.format(fathmm_database)
                        self.mysql(cursor, query)
                        query = 'CREATE SCHEMA `{}`'.format(fathmm_database)
                        self.mysql(cursor, query)
                    self.print_log('FATHMM database `{}` created.'.format(fathmm_database))
                    #create tables and import FATHMM mysqldump
                    self.print_log('Importing FATHMM data from {} into database `{}`. This will take ~1 hour'.format(fathmm_data_file, fathmm_database))
                    os.system('gunzip < {} | mysql -h {} -P {} -u {} -p{} {}'.format(fathmm_data_file, host, port, user, password, fathmm_database))
                    self.print_log('FATHMM data imported.')
                    #download Python script
                    fathmm_py_url = config.get('fathmm', 'FATHMM_PY_URL')
                    fathmm_py_file = os.path.join(fathmm_dir, os.path.basename(fathmm_py_url))
                    self.print_log('Downloading current FATHMM Python script from {} into {}.'.format(fathmm_py_url, fathmm_py_file))
                    urlretrieve(fathmm_py_url, fathmm_py_file)
                    self.print_log('FATHMM Python script downloaded.')
                    self.print_log('Copying FATHMM Python script from {} into {}.'.format(fathmm_py_file, fathmm_plugin_dir))
                    shutil.copy(fathmm_py_file, fathmm_plugin_dir)
                    self.print_log('FATHMM Python script copied.')
                    #create FATHMM config.ini
                    fathmm_config_file = os.path.join(fathmm_plugin_dir, 'config.ini')
                    self.print_log('Creating FATHMM config.ini file in {}.'.format(fathmm_plugin_dir))
                    with open(fathmm_config_file , 'w') as fathmm_config_fileh:
                        fathmm_config_fileh.write("""[DATABASE]
HOST = {}
PORT = {}
USER = {}
PASSWD = {}
DB = {}
""".format(host, port, vax_user, vax_password, fathmm_database))
                        
                    with open(done_file, 'w'):
                        pass
                    self.print_log('Installed FATHMM.')

################################################################
#                             GATK                             #
################################################################

        if 'all' in install or 'gatk' in install:
            gatk_current_link = os.path.join(apps_dir,'gatk-current')
            done_file = os.path.join(gatk_current_link,'.gatk.done')
            if my.file_exists(done_file):
                self.print_log('GATK already downloaded. To redownload "rm {}"'.format(done_file))
            else:
                self.print_log('Installing GATK.')
                gatk_date = '{}'.format(strftime("%Y-%m-%d", localtime()))
                gatk_date_dir = os.path.join(gatk_dir, gatk_date)
                if os.path.isdir(gatk_date_dir):
                    shutil.rmtree(gatk_date_dir)
                    
                self.print_log('Downloading GATK to {}.'.format(gatk_date_dir))
                gatk_git_url = config.get('gatk','GATK_GIT_URL')
                git.clone(gatk_git_url, gatk_date_dir)
                self.print_log('Downloaded GATK to {}.'.format(gatk_date_dir))
                
                self.print_log('Determining GATK version.')
                #get GATK version from the pom.xml file
                pom_xml = os.path.join(gatk_date_dir,'pom.xml')
                gatk_version = None
                if my.file_exists(pom_xml):
                    tree = ET.ElementTree(file=pom_xml)
                    root = tree.getroot()
                    #get namespace
                    #http://stackoverflow.com/questions/1953761/accessing-xmlns-attribute-with-python-elementree
                    if root.tag[0] == '{':
                        uri, ignore, tag = root.tag[1:].partition('}')
                    else:
                        uri = None
                        tag = root.tag
                    ns_version = 'parent/version' if uri == None else '{{{0}}}parent/{{{0}}}version'.format(uri)
                    #parent_node = root.find(ns_parent)
                    version_node = root.find(ns_version)
                    if version_node != None:
                        gatk_version = version_node.text
                if gatk_version:
                    self.print_log('GATK version is {}.'.format(gatk_version))                
                else:
                    gatk_version = gatk_date
                    self.print_log('WARNING: Cannot determine version of GATK; using {} as version.'.format(gatk_date))
                versions['gatk'] = gatk_version
    
                self.print_log('Compiling GATK in {}. This will take ~30 minutes'.format(gatk_date_dir))
                with cd(gatk_date_dir):
                    sh.mvn('verify')
                self.print_log('Compiled GATK in {}.'.format(gatk_date_dir))
                
                self.print_log('Linking to GATK and Queue jar files.')
                #link gatk-current -> gatk_version -> gatk_date_dir
                gatk_version_link = os.path.join(gatk_dir, 'gatk_'+gatk_version)
                if os.path.isdir(gatk_version_link):
                    shutil.rmtree(gatk_version_link)
                elif os.path.islink(gatk_version_link):
                    try:
                        os.remove(gatk_version_link)
                    except OSError:
                        pass
                os.symlink(gatk_date_dir, gatk_version_link)
                try:
                    os.remove(gatk_current_link)
                except OSError:
                    pass
                os.symlink(gatk_version_link, gatk_current_link)
                #link GenomeAnalysisTK.jar and Queue.jar
                gatk_genomeanalysistk_jar_subpath = config.get('gatk', 'GATK_GENOMEANALYSISTK_JAR_SUBPATH')
                gatk_genomeanalysistk_jar = os.path.join(gatk_current_link, gatk_genomeanalysistk_jar_subpath)
                gatk_genomeanalysistk_link = os.path.join(apps_dir, 'GenomeAnalysisTK.jar')
                try:
                    os.remove(gatk_genomeanalysistk_link)
                except OSError:
                    pass
                os.symlink(gatk_genomeanalysistk_jar, gatk_genomeanalysistk_link)
                gatk_queue_jar_subpath = config.get('gatk', 'GATK_QUEUE_JAR_SUBPATH')
                gatk_queue_jar = os.path.join(gatk_current_link, gatk_queue_jar_subpath)
                gatk_queue_link = os.path.join(apps_dir, 'Queue.jar')
                try:
                    os.remove(gatk_queue_link)
                except OSError:
                    pass
                os.symlink(gatk_queue_jar, gatk_queue_link)
                self.print_log('Linked to GATK and Queue jar files.')
    
                with open(done_file, 'w'):
                    pass
                self.print_log('Installed GATK.')

################################################################
#             Human Gene Mutation Database (HGMD)              #
################################################################
        
        if 'all' in install or 'hgmd' in install:
            if not hgmd_user:
                hgmd_user = config.get('hgmd', 'HGMD_USER')
            if not hgmd_password:
                hgmd_password = config.get('hgmd', 'HGMD_PW')
            if not hgmd_version:
                hgmd_version = config.get('hgmd', 'HGMD_VERSION')
            versions['hgmd'] = hgmd_version
            if not hgmd_dir:
                hgmd_dir = os.path.join(config.get('hgmd', 'HGMD_DIR'), hgmd_version)
            if not install_hgmd_schemas:
                install_hgmd_schemas = config.get('hgmd', 'INSTALL_HGMD_SCHEMAS').split(',')
            if not hgmd_database:
                hgmd_database = config.get('hgmd', 'HGMD_DATABASE')
    
            done_file = os.path.join(hgmd_dir,'.hgmd.done')
            if my.file_exists(done_file):
                self.print_log('HGMD already installed. To reinstall "rm {}"'.format(done_file))
            else:
                self.print_log('Installing HGMD. This will take ~10 minutes')
                #drop/recreate databases and import data
                #hgmd_views must be last
                for schema in install_hgmd_schemas:
                    dump_file = os.path.join(hgmd_dir, '{}-{}.dump.gz'.format(schema, hgmd_version))
                    if not my.file_exists(dump_file):
                        self.print_log('Cannot find HGMD dump file {}; skipping'.format(dump_file))
                    else:
                        self.print_log('Dropping and recreating HGMD schema `{}`'.format(schema))
                        query = ['DROP DATABASE IF EXISTS {}'.format(schema),
                                 'CREATE DATABASE {}'.format(schema),]
                        with MySQLdb.connect(host,user,password,database,port=port) as cursor:
                            self.mysql(cursor, query)
                        self.print_log('Importing HGMD schema `{}` from {}'.format(schema, dump_file))
                        os.system('gunzip < {} | mysql -h {} -P {} -u {} -p{} {}'.format(dump_file, host, port, user, password, schema))
                        self.print_log('Imported HGMD schema `{}` from {}'.format(schema, dump_file))
                        #add indices and stored procedures (the stored procedures go into the vax database)
                        if schema == 'hgmd_pro':
                            query = ["""ALTER TABLE `hgmd_pro`.`hg19_coords`
ADD INDEX IX_chromosome (chromosome ASC)
, ADD INDEX IX_coordSTART (coordSTART ASC)
, ADD INDEX IX_coordEND (coordEND ASC)""",
'DROP procedure IF EXISTS `{}`.`coord2hgmd`'.format(database),
"""CREATE DEFINER=CURRENT_USER PROCEDURE `{}`.`coord2hgmd`(chromosome VARCHAR(2), coordSTART INT(11), coordEND INT(11))
BEGIN
SELECT DISTINCT
allmut.acc_num,
allmut.disease,
allmut.tag,
allmut.base,
allmut.hgvs,
allmut.codon,
allmut.amino,
allmut.deletion,
allmut.insertion,
allmut.descr,
allmut.pmid
FROM hgmd_pro.hg19_coords
JOIN hgmd_pro.allmut
ON hg19_coords.acc_num = allmut.acc_num
WHERE hg19_coords.chromosome = chromosome
AND (hg19_coords.coordSTART BETWEEN  coordSTART AND coordEND
OR hg19_coords.coordEND BETWEEN  coordSTART AND coordEND);
END""".format(database),
'DROP procedure IF EXISTS `{}`.`gene2hgmd_disease`'.format(database),
"""CREATE DEFINER=CURRENT_USER PROCEDURE `{}`.`gene2hgmd_disease`(gene VARCHAR(10))
BEGIN
SELECT distinct disease
from hgmd_pro.allgenes
where allgenes.gene = gene;
END""".format(database),
]
                            with MySQLdb.connect(host,user,password,database,port=port) as cursor:
                                self.mysql(cursor, query)
                        #permissions
                        query = ["""GRANT SELECT, INSERT, UPDATE, CREATE TEMPORARY TABLES ON {}.* TO 'hgmduser'@'%' IDENTIFIED BY 'mootation'""".format(schema),
                                 """GRANT SELECT, EXECUTE ON {}.* TO '{}'@'%' IDENTIFIED BY '{}'""".format(schema, hgmd_user, hgmd_password),]
                        with MySQLdb.connect(host,user,password,database,port=port) as cursor:
                            self.mysql(cursor, query)
    
                with open(done_file, 'w'):
                    pass
                self.print_log('Installed HGMD.')

################################################################
#                Human Metabolome Database (HMDB)              #
################################################################

        if 'all' in install or 'metabolome' in install or 'hmd' in install:
            done_file = os.path.join(metabolome_dir,'.metabolome.done')
            if my.file_exists(done_file):
                self.print_log('Metabolome data already installed. To reinstall "rm {}"'.format(done_file))
            else:
                self.print_log('Installing metabolome data. This will take ~10 minutes.')
    
                #download metabolome file
                #TODO: scrape url for version
                versions['hmdb'] = '3.6'
                #versions['hmdb'] = 'HMDB accessed {}'.format(strftime("%Y-%m-%d", localtime()))
                url = config.get('metabolome','METABOLOME_URL')
                zipbase = os.path.basename(url)
                table = my.r_strip(zipbase, '.zip')
                zipfile = os.path.join(metabolome_dir, zipbase)
                text_file = my.swap_ext(zipfile, '.zip', '.txt')
                sql_file = text_file+'.mysql'
                unzipdir = metabolome_dir
                
                self.getzip(url, zipfile, unzipdir)
                xml2db_done_file = os.path.join(metabolome_dir, '.'+os.path.basename(text_file)+'.done')
                
                if my.file_exists(xml2db_done_file):
                    self.print_log('Metabolome data already extracted from xml. To extract "rm {}"'.format(xml2db_done_file))
                else:
                    self.print_log('Extracting metabolome data from xml.')
                    vax_metabolome_xml2db.run(input=os.path.join(unzipdir, 'HMDBP*.xml'), output=text_file)
                    with open(xml2db_done_file, 'w'):
                        pass
                    self.print_log('Extracted metabolome data from xml.')
                
                sql_done_file = os.path.join(os.path.dirname(sql_file), '.'+os.path.basename(sql_file)+'.done')
                if my.file_exists(sql_done_file):
                    self.print_log('{} already created. To redo "rm {}"'.format(sql_file, sql_done_file))
                else:
                    self.print_log('Making CREATE TABLE `{}` script.'.format(table))
                    sql_columns.run(text_file, database=database, schema=database, table=table, primary_key=['gene'])
                    with open(sql_done_file, 'w'):
                        pass
                    self.print_log('Made CREATE TABLE `{}` script.'.format(table))
                    
                self.print_log('Creating `{}` table.'.format(table))
                with open(sql_file, 'r') as sql:
                    sh.mysql('-h', host, '-P', port, '-u', user, '-p{}'.format(password), database, _in=sql)
                    
                self.print_log('Importing data to `{}`.`{}` table.'.format(database, table))
                sh.mysqlimport('-h', host, '-P', port, '-u', user, '-p{}'.format(password), '--delete', '--local', '--ignore-lines=1', database, text_file)
                self.print_log('Imported data to `{}`.`{}` table.'.format(database, table))
                                    
                with MySQLdb.connect(host,user,password,database,port=port) as cursor:
                    proc = 'hgnc2metabolome'
                    self.print_log('Creating {} stored procedure'.format(proc))
                    query = [
'DROP PROCEDURE IF EXISTS {}'.format(proc),
"""CREATE DEFINER=CURRENT_USER PROCEDURE {}(hgnc VARCHAR(255))
BEGIN
SELECT `function` AS metabolic_function, `metabolites`
FROM `{}` AS m
WHERE m.`gene` = hgnc;
END
""".format(proc, table)]
                    self.mysql(cursor, query)
                    self.print_log('Created {} stored procedure'.format(proc))
                    
                with open(done_file, 'w'):
                    pass
                self.print_log('Installed Metabolome data.')

################################################################
#                  Human Protein Atlas (HPA)                   #
################################################################

        if 'all' in install or 'human_protein_atlas' in install or 'hpa' in install:
            done_file = os.path.join(hpa_dir,'.hpa.done')
            if my.file_exists(done_file):
                self.print_log('Human Protein Atlas already installed. To reinstall "rm {}"'.format(done_file))
            else:
                self.print_log('Installing Human Protein Atlas. This will take ~10 mins')
    
                #download HPA data and make CREATE TABLE scripts
                #TODO: scrape url for version
                versions['hpa'] = '12'
                #versions['hpa'] = 'HPA accessed {}'.format(strftime("%Y-%m-%d", localtime()))
                url = config.get('hpa', 'HPA_NORMAL_TISSUE_URL')
                zipbase = os.path.basename(url)
                table = 'hpa_'+my.r_strip(zipbase, '.csv.zip')
                zipfile = os.path.join(hpa_dir, zipbase)
                unzipped_file = os.path.join(hpa_dir, my.r_strip(zipbase, '.zip'))
                text_file = os.path.join(hpa_dir, table+'.txt')
                sql_file = text_file+'.mysql'
                self.print_log('Downloading {} to {}.'.format(url, unzipped_file))
                self.getzip(url, zipfile, hpa_dir)
                self.print_log('Downloaded {} to {}.'.format(url, unzipped_file,))
                self.print_log('Converting {} to {}.'.format(unzipped_file, text_file))
                csv2tab.run(unzipped_file, text_file)
                self.print_log('Converted {} to {}.'.format(unzipped_file, text_file))
                self.print_log('Making CREATE TABLE `{}` script.'.format(table))
                sql_columns.run(text_file, database=database, schema=database, table=table, primary_key=['Gene','Tissue','Cell_type'], header_line=1)
                self.print_log('Made CREATE TABLE `{}` script.'.format(table))
                self.print_log('Creating `{}` table.'.format(table))
                with open(sql_file, 'r') as sql:
                    sh.mysql('-h', host, '-P', port, '-u', user, '-p{}'.format(password), database, _in=sql)
                self.print_log('Importing data to `{}`.`{}` table.'.format(database, table))
                sh.mysqlimport('-h', host, '-P', port, '-u', user, '-p{}'.format(password), '--delete', '--local', '--ignore-lines=1', database, text_file)
                self.print_log('Imported data to `{}`.`{}` table.'.format(database, table))
                with MySQLdb.connect(host,user,password,database,port=port) as cursor:

                    proc = 'ensg2hpa_tissue'
                    self.print_log('Creating {} stored procedure'.format(proc))
                    query = [
'DROP PROCEDURE IF EXISTS {}'.format(proc),
"""CREATE DEFINER=CURRENT_USER PROCEDURE {}()
BEGIN
SELECT DISTINCT h.`Tissue`, h.`Cell_type`, h.`Level`
FROM `hpa_normal_tissue` AS h
ORDER BY h.`Tissue`, h.`Cell_type`;
END
""".format(proc)]
                    self.mysql(cursor, query)
                    self.print_log('Created {} stored procedure'.format(proc))
    
                    proc = 'hpa_tissue_cell_type_list'
                    self.print_log('Creating {} stored procedure'.format(proc))
                    query = [
'DROP PROCEDURE IF EXISTS {}'.format(proc),
"""CREATE DEFINER=CURRENT_USER PROCEDURE {}()
BEGIN
SELECT DISTINCT h.`Tissue`, h.`Cell_type`
FROM `hpa_normal_tissue` AS h
ORDER BY h.`Tissue`, h.`Cell_type`;
END
""".format(proc)]
                    self.mysql(cursor, query)
                    self.print_log('Created {} stored procedure'.format(proc))
    
                url = config.get('hpa', 'HPA_SUBCELLULAR_LOCATION_URL')
                zipbase = os.path.basename(url)
                table = 'hpa_'+my.r_strip(zipbase, '.csv.zip')
                zipfile = os.path.join(hpa_dir, zipbase)
                unzipped_file = os.path.join(hpa_dir, my.r_strip(zipbase, '.zip'))
                text_file = os.path.join(hpa_dir, table+'.txt')
                sql_file = text_file+'.mysql'
                self.print_log('Downloading {} to {}.'.format(url, unzipped_file))
                self.getzip(url, zipfile, hpa_dir)
                self.print_log('Downloaded {} to {}.'.format(url, unzipped_file))
                self.print_log('Converting {} to {}.'.format(unzipped_file, text_file))
                csv2tab.run(unzipped_file, text_file)
                self.print_log('Converted {} to {}.'.format(unzipped_file, text_file))
                self.print_log('Making CREATE TABLE `{}` script.'.format(table))
                sql_columns.run(text_file, database=database, schema=database, table=table, primary_key=['Gene'], header_line=1)
                self.print_log('Made CREATE TABLE `{}` script.'.format(table))
                self.print_log('Creating `{}` table.'.format(table))
                with open(sql_file, 'r') as sql:
                    sh.mysql('-h', host, '-P', port, '-u', user, '-p{}'.format(password), database, _in=sql)
                self.print_log('Importing data to `{}`.`{}` table.'.format(database, table))
                sh.mysqlimport('-h', host, '-P', port, '-u', user, '-p{}'.format(password), '--delete', '--local', '--ignore-lines=1', database, text_file)
                self.print_log('Imported data to `{}`.`{}` table.'.format(database, table))
                with MySQLdb.connect(host,user,password,database,port=port) as cursor:
                    proc = 'ensg2hpa_subcellular_location'
                    self.print_log('Creating {} stored procedure'.format(proc))
                    query = [
'DROP PROCEDURE IF EXISTS {}'.format(proc),
"""CREATE DEFINER=CURRENT_USER PROCEDURE {}(ensg CHAR(15))
BEGIN
    SELECT DISTINCT h.`Main_location` AS location
    FROM `hpa_subcellular_location` AS h
    WHERE h.`Gene` = ensg
    UNION
    SELECT DISTINCT h.`Other_location` AS location
    FROM `hpa_subcellular_location` AS h
    WHERE h.`Gene` = ensg
    ORDER BY location;
END
""".format(proc)]
                    self.mysql(cursor, query)
                    self.print_log('Created {} stored procedure'.format(proc))
                    
                with open(done_file, 'w'):
                    pass
                self.print_log('Installed Human Protein Atlas.')

################################################################
#                          MitoCarta                           #
################################################################

        if 'all' in install or 'mitocarta' in install:
            done_file = os.path.join(mitocarta_dir,'.mitocarta.done')
            if my.file_exists(done_file):
                self.print_log('MitoCarta already installed. To reinstall "rm {}"'.format(done_file))
            else:
                self.print_log('Installing MitoCarta.')
    
                #download mitocarta
                versions['mitocarta'] = 'accessed {}'.format(strftime("%Y-%m-%d", localtime()))
                url = config.get('mitocarta', 'MITOCARTA_TXT_URL')
                table = 'HumanMitoCartaAll'
                text_file = os.path.join(mitocarta_dir, table+'.txt')
                sql_file = text_file+'.mysql'
                self.print_log('Downloading MitoCarta data to {}.'.format(text_file))
                urlretrieve(url, text_file)
                self.print_log('Downloaded MitoCarta data.')
                self.print_log('Making CREATE TABLE `{}` script.'.format(table))
                sql_columns.run(text_file, database=database, schema=database, table=table, primary_key=['SYM'], header_line=1)
                self.print_log('Made CREATE TABLE `{}` script.'.format(table))
                self.print_log('Creating `{}` table.'.format(table))
                with open(sql_file, 'r') as sql:
                    sh.mysql('-h', host, '-P', port, '-u', user, '-p{}'.format(password), database, _in=sql)
                self.print_log('Importing data to `{}`.`{}` table.'.format(database, table))
                sh.mysqlimport('-h', host, '-P', port, '-u', user, '-p{}'.format(password), '--delete', '--local', '--ignore-lines=1', database, text_file)
                self.print_log('Imported data to `{}`.`{}` table.'.format(database, table))
                with MySQLdb.connect(host,user,password,database,port=port) as cursor:
                    proc = 'get_mitocarta_gene'
                    self.print_log('Creating {} stored procedure'.format(proc))
                    query = [
'DROP PROCEDURE IF EXISTS {}'.format(proc),
"""CREATE DEFINER=CURRENT_USER PROCEDURE {}(hgnc varchar(50))
BEGIN
SELECT m.`SYM` AS mito_gene
FROM `HumanMitoCartaAll` AS m
WHERE m.`SYM` = hgnc AND m.`MITOCARTA_LIST` = 1;
END
""".format(proc)]
                    self.mysql(cursor, query)
                    self.print_log('Created {} stored procedure'.format(proc))
                    
                with open(done_file, 'w'):
                    pass
                self.print_log('Installed MitoCarta.')

################################################################
#                 Mouse Genes/Phenotypes (MGI)                 #
################################################################

        #This database is used for humans so it will be installed without regard to species parameters.
        if 'all' in install or 'mgi' in install:
            done_file = os.path.join(mgi_dir,'.mgi.done')
            if my.file_exists(done_file):
                self.print_log('Mouse Genes/Phenotypes already installed. To reinstall "rm {}"'.format(done_file))
            else:
                self.print_log('Installing MGI Mouse Genes/Phenotypes.')

                versions['mgi'] = 'accessed {}'.format(strftime("%Y-%m-%d", localtime()))
                mgi_url = config.get('mgi', 'MGI_URL')
                mgi_files = config.get('mgi', 'MGI_FILES')
                for f in mgi_files.split(' '):
                    url = os.path.join(mgi_url, f)
                    table = os.path.splitext(f)[0]
                    rpt_file = os.path.join(mgi_dir, f)
                    text_file = my.swap_ext(rpt_file, '.rpt', '.txt')
                    sql_file = text_file+'.mysql'
                    self.print_log('Downloading {} to {}.'.format(url, rpt_file))
                    urlretrieve(url, rpt_file)
                    self.print_log('Downloaded {} to {}.'.format(url, rpt_file))

                    self.print_log('Converting {} to {}.'.format(rpt_file, text_file))
                    vax_mgi_file_cleanup_sql.run([rpt_file], database, database)
                    self.print_log('Converted {} to {}.'.format(rpt_file, text_file))
    
                    self.print_log('Creating `{}` table.'.format(table))
                    with open(sql_file, 'r') as sql:
                        sh.mysql('-h', host, '-P', port, '-u', user, '-p{}'.format(password), database, _in=sql)
                    self.print_log('Importing `{}` data to `{}` database.'.format(table, database))
                    sh.mysqlimport('-h', host, '-P', port, '-u', user, '-p{}'.format(password), '--delete', '--local', '--ignore-lines=1', database, text_file)
                    self.print_log('Imported `{}` to `{}` database.'.format(table, database))

                #table with human and mouse records joined
                table = 'HOM_MouseHumanSequence_mm_hs'
                self.print_log('Creating {} table'.format(table))
                query = [
'DROP TABLE IF EXISTS {}'.format(table),
"""CREATE TABLE `{}`
(PRIMARY KEY (mm_HomoloGene_ID, hs_HomoloGene_ID), KEY(mm_Symbol), KEY(hs_Symbol))
AS
SELECT
m.HomoloGene_ID AS mm_HomoloGene_ID,
m.Common_Organism_Name AS mm_Common_Organism_Name,
m.NCBI_Taxon_ID AS mm_NCBI_Taxon_ID,
m.Symbol AS mm_Symbol,
m.EntrezGene_ID AS mm_EntrezGene_ID,
m.Mouse_MGI_ID AS mm_Mouse_MGI_ID,
m.HGNC_ID AS mm_HGNC_ID,
m.OMIM_Gene_ID AS mm_OMIM_Gene_ID,
m.Genetic_Location AS mm_Genetic_Location,
m.Genomic_Coordinates_mouse_human AS mm_Genomic_Coordinates_mouse_human,
m.Nucleotide_RefSeq_IDs AS mm_Nucleotide_RefSeq_IDs,
m.Protein_RefSeq_IDs AS mm_Protein_RefSeq_IDs,
m.SWISS_PROT_IDs AS mm_SWISS_PROT_IDs,
h.HomoloGene_ID AS hs_HomoloGene_ID, \
h.Common_Organism_Name AS hs_Common_Organism_Name,
h.NCBI_Taxon_ID AS hs_NCBI_Taxon_ID,
h.Symbol AS hs_Symbol,
h.EntrezGene_ID AS hs_EntrezGene_ID,
h.Mouse_MGI_ID AS hs_Mouse_MGI_ID,
h.HGNC_ID AS hs_HGNC_ID,
h.OMIM_Gene_ID AS hs_OMIM_Gene_ID,
h.Genetic_Location AS hs_Genetic_Location,
h.Genomic_Coordinates_mouse_human AS hs_Genomic_Coordinates_mouse_human,
h.Nucleotide_RefSeq_IDs AS hs_Nucleotide_RefSeq_IDs,
h.Protein_RefSeq_IDs AS hs_Protein_RefSeq_IDs,
h.SWISS_PROT_IDs AS hs_SWISS_PROT_IDs
FROM mus.HOM_MouseHumanSequence AS m
JOIN mus.HOM_MouseHumanSequence AS h
on m.NCBI_Taxon_ID = 10090 and h.NCBI_Taxon_ID = 9606 and m.HomoloGene_ID =h.HomoloGene_ID
""".format(table)]
                with MySQLdb.connect(host,user,password,database,port=port) as cursor:
                    self.mysql(cursor, query)
                    self.print_log('Created {} table'.format(table))

                #vax proc
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
                self.print_log('Installed MGI Mouse Genes/Phenotypes.')

################################################################
#                          NHLBI-EVS                           #
################################################################

        if 'all' in install or 'nhlbi' in install or 'evs' in install:
            done_file = os.path.join(nhlbi_dir,'.nhlbi.done')
            if my.file_exists(done_file):
                self.print_log('NHLBI-EVS already installed. To reinstall "rm {}"'.format(done_file))
            else:
                self.print_log('Installing NHLBI-EVS.')
                default_nhlbi_vcf_file_url = None
                if not nhlbi_vcf_file_url:
                    try:
                        default_nhlbi_vcf_file_url = config.get('nhlbi', 'DEFAULT_NHLBI_VCF_FILE_URL')
                    except:
                        pass
                    #determine version if not specified as a run parameter
                    self.print_log('Getting list of NHLBI-EVS VCF file urls.')
                    nhlbi_info_url = config.get('nhlbi', 'NHLBI_INFO_URL')
                    nhlbi_vcf_file_pattern = config.get('nhlbi', 'NHLBI_VCF_FILE_PATTERN')
                    nhlbi_vcf_file_urls = list(set(my.list_hrefs(nhlbi_info_url, nhlbi_vcf_file_pattern)))
                    if not nhlbi_vcf_file_urls:
                        self.print_log('No urls at {} match the pattern {}. Check the web site and NHLBI_VCF_FILE_PATTERN in vax.cfg.'.format(nhlbi_info_url, nhlbi_vcf_file_pattern))
                    else:
                        if len(nhlbi_vcf_file_urls) == 1:
                            nhlbi_vcf_file_url = nhlbi_vcf_file_urls[0]
                        else:
                            nhlbi_vcf_file_url = my.select_item(nhlbi_vcf_file_urls, default_nhlbi_vcf_file_url, 'NHLBI-EVS VCF files', 'NHLBI-EVS VCF file to install')
                if not nhlbi_vcf_file_url:
                    self.print_log('WARNING: No NHLBI-EVS VCF file selected. Rerun or install manually.')
                else:
                    #download selected version
                    nhlbi_vcf_file = os.path.join(nhlbi_dir, os.path.basename(nhlbi_vcf_file_url))
                    nhlbi_version = os.path.basename(nhlbi_vcf_file_url).split('.')[0]
                    versions['nhlbi'] = nhlbi_version
                    self.print_log('Downloading NHLBI-EVS VCF file {} to {}.'.format(nhlbi_vcf_file_url, nhlbi_vcf_file))
                    urlretrieve(nhlbi_vcf_file_url, nhlbi_vcf_file)
                    self.print_log('NHLBI-EVS VCF downloaded.')
                    #remember version as default
                    self.print_log('setting NHLBI-EVS VCF file url {} as default for next install.'.format(nhlbi_vcf_file_url))
                    config.set('nhlbi', 'DEFAULT_NHLBI_VCF_FILE_URL', nhlbi_vcf_file_url)
                    
                    #process VCF file to create data for mySQL table
                    self.print_log('Scheduling job to create and import `{}`.`{}` table from {}.'.format(database, table, nhlbi_vcf_file))
                    vax_nhlbi_vcfs2db_py = os.path.join(os.path.dirname(__file__), 'vax_nhlbi_vcfs2db.py')
                    forbidden = re.compile(r'[^0-9,a-z,A-Z$_]')
                    table = forbidden.sub('_', my.r_strip(nhlbi_vcf_file, '.vcf.tar.gz'))[:64]
                    output_prefix = os.path.join(os.path.dirname(nhlbi_vcf_file), table+'.txt')
                    stored_procedure = 'get_nhlbi_allele_frequencies'
                    cmd= '{} {} --input {} --output {} --database {} --table {} --proc {} --host {} --port {} --user {} --password {}'.format(
                        config.get('DEFAULT','python'), vax_nhlbi_vcfs2db_py, nhlbi_vcf_file, output_prefix, database, table, stored_procedure, host, port, user, password)
                    job_name = '{}_import'.format(table)
                    job = my.run_job(cmd, job_name, job_dir)
                    nhlbi_import_job_id = job.jobId
                    self.all_jobs.append(job.jobId)
                    self.print_log('Scheduled job {}: {}.'.format(job.jobId, cmd))
                    
                    with open(done_file, 'w'):
                        pass
                    self.print_log('Installed NHLBI-EVS.')


################################################################
#                             OMIM                             #
################################################################

        if 'all' in install or 'omim' in install:
            done_file = os.path.join(omim_dir,'.omim.done')
            if my.file_exists(done_file):
                self.print_log('OMIM already installed. To reinstall "rm {}"'.format(done_file))
            else:
                self.print_log('Installing OMIM. This will take ~10 minutes.')
    
                #download OMIM directory
                versions['omim'] = 'accessed {}'.format(strftime("%Y-%m-%d", localtime()))
                omim_ftp = config.get('omim', 'OMIM_FTP')
                omim_ftp_dir = config.get('omim', 'OMIM_FTP_DIR')
                self.download_ftp_dir(omim_ftp, omim_ftp_dir, omim_dir)
                self.print_log('Creating OMIM tab-delimited files.')
                omim2db.run(omim_dir)
                self.print_log('Created OMIM tab-delimited files.')
    
                table = 'omim_genemap'
                text_file = os.path.join(omim_dir, table+'.txt')
                sql_file = text_file+'.mysql'
                self.print_log('Making CREATE TABLE `{}` script.'.format(table))
                sql_columns.run(text_file, database=database, schema=database, table=table, indexes=['gene'])
                self.print_log('Made CREATE TABLE `{}` script.'.format(table))
                self.print_log('Creating `{}` table.'.format(table))
                with open(sql_file, 'r') as sql:
                    sh.mysql('-h', host, '-P', port, '-u', user, '-p{}'.format(password), database, _in=sql)
                self.print_log('Importing data to `{}`.`{}` table.'.format(database, table))
                sh.mysqlimport('-h', host, '-P', port, '-u', user, '-p{}'.format(password), '--delete', '--local', '--ignore-lines=1', database, text_file)
                self.print_log('Imported data to `{}`.`{}` table.'.format(database, table))
                with MySQLdb.connect(host, user, password, database, port=port) as cursor:
                    proc = 'hgnc2omim'
                    self.print_log('Creating {} stored procedure'.format(proc))
                    query = [
'DROP PROCEDURE IF EXISTS {}'.format(proc),
"""CREATE DEFINER=CURRENT_USER PROCEDURE {}(hgnc VARCHAR(255))
BEGIN
SELECT distinct CAST(CONCAT(o.`disorder`, ' (', o.`mim_locus`, ')') AS CHAR(255)) AS disorder
FROM `omim_genemap` AS o
WHERE o.`gene` = hgnc
ORDER BY o.`disorder`;
END
""".format(proc)]
                    self.mysql(cursor, query)
                    self.print_log('Created {} stored procedure'.format(proc))
                    
                with open(done_file, 'w'):
                    pass
                self.print_log('Installed OMIM.')

################################################################
#                           RefGene                            #
################################################################

    #TODO: provide for non-human species
        #requires faidx
        #creates files for picard/gatk intervals and gene list
        if 'all' in install or 'refgene' in install:
            done_file = os.path.join(refgene_dir, '.refgene.done')
            if my.file_exists(done_file):
                self.print_log('RefGene already installed. To reinstall "rm {}"'.format(done_file))
            else:
                self.print_log('Installing RefGene data and interval_list from ucsc. This will take ~10 minutes.')
                ucsc_host = config.get('refgene', 'UCSC_HOST')
                ucsc_user = config.get('refgene', 'UCSC_USER')
                ucsc_db = config.get('refgene', 'UCSC_DB')
                ucsc_table = config.get('refgene', 'UCSC_TABLE')
                ucsc_genome_id = config.get('refgene', 'UCSC_GENOME_ID')
                genome_id = config.get('refgene', 'GENOME_ID')
                ucsc_genome_table = ucsc_table+'_ucsc_'+ucsc_genome_id
                table = ucsc_table+'_'+genome_id+'_ucsc_'
                ucsc_text_file = os.path.join(refgene_dir, ucsc_genome_table+'.txt')
                text_file = os.path.join(refgene_dir, table+'.txt')
                cds_intervals_file = my.swap_ext(text_file, '.txt', '.cds.interval_list')
                exon_intervals_file = my.swap_ext(text_file, '.txt', '.exon.interval_list')
                gene_intervals_file = my.swap_ext(text_file, '.txt', '.gene.interval_list')
                sql_file = text_file+'.mysql'
                faidx_table = my.r_strip(os.path.basename(config.get('faidx', 'FAIDX_DECOY_FILE')), '.txt')
                interval_list_decoy_header = config.get('faidx', 'INTERVAL_LIST_DECOY_HEADER')
    
                #download UCSC refGene table from UCSC mysql server
                versions['refgene'] = 'accessed {}'.format(strftime("%Y-%m-%d", localtime()))
                ucsc_done_file = my.done_file(ucsc_text_file)
                if my.file_exists(ucsc_done_file):
                    self.print_log('`{}`.`{}` already downloaded. To redownload "rm {}"'.format(ucsc_db, ucsc_table, ucsc_done_file))
                else:
                    #download refGene table from UCSC mysql server
                    self.print_log('Selecting `{}`.`{}` from ucsc into {}.'.format(ucsc_db, ucsc_table, ucsc_text_file))
                    sh.mysql('--host', ucsc_host,'--user', ucsc_user, '--column-names', '--batch', '-e', 'SELECT * FROM {}.{};'.format(ucsc_db, ucsc_table), _out=ucsc_text_file)
                    with open(ucsc_done_file, 'w'):
                        pass
                    self.print_log('Selected `{}`.`{}` from ucsc into {}.'.format(ucsc_db, ucsc_table, ucsc_text_file))
                    
                #convert downloaded UCSC refGene table chroms to Ensembl/GATK style and sort
                convert_done_file = my.done_file(text_file)
                if my.file_exists(convert_done_file):
                    self.print_log('{} already converted and sorted as {}. To reconvert "rm {}"'.format(ucsc_text_file, text_file, convert_done_file))
                else:
                    self.print_log('Converting {} to {}.'.format(ucsc_text_file, text_file))
                    convert_contigs_sort.run(input=ucsc_text_file, output=text_file, chrom_col=2, pos_col=4, faidx=faidx_decoy_file, input_genome_id=ucsc_genome_id, output_genome_id=genome_id)
                    with open(convert_done_file, 'w'):
                        pass
                    self.print_log('Converted {} to {}.'.format(ucsc_text_file, text_file))
                    
                #create picard exon.intervals_list and gene.intervals_list files
                exon_intervals_done_file = my.done_file(exon_intervals_file)
                gene_intervals_done_file = my.done_file(gene_intervals_file)
                if my.file_exists(exon_intervals_done_file) and my.file_exists(gene_intervals_done_file):
                    self.print_log('{} and {} already created. To reconvert "rm {} {}"'.format(exon_intervals_file, gene_intervals_file, exon_intervals_done_file, gene_intervals_done_file))
                else:
                    self.print_log('Creating {} and {} from {}.'.format(exon_intervals_file, gene_intervals_file, text_file))
                    create_interval_list.run(input=text_file, exon_interval_list=exon_intervals_file, gene_interval_list=gene_intervals_file, chrom_col=2, strand_col=3, name_cols=[12,1], splice_bases=2, tx_start_col=4, tx_end_col=5, cds_start_col=6, cds_end_col=7, exon_starts_col=9, exon_ends_col=10, faidx=faidx_decoy_file, interval_list_header=interval_list_decoy_header, genome_id=genome_id, zero_based=True)
                    #create_interval_list will make .done files
                
                #TODO: interval list tables (and ucsc table?) 
                #create refGene mysql table
                sql_done_file = my.done_file(sql_file)
                if my.file_exists(sql_done_file):
                    self.print_log('{} already created. To redo "rm {}"'.format(sql_file, sql_done_file))
                else:
                    self.print_log('Making CREATE TABLE `{}` script.'.format(table))
                    sql_columns.run(text_file, database=database, schema=database, table=table, primary_key=['name'], indexes=['name2', 'chrom,txStart,txEnd'], header_line=1)
                    with open(sql_done_file, 'w'):
                        pass
                    self.print_log('Made CREATE TABLE `{}` script.'.format(table))
                    
                self.print_log('Creating `{}` table.'.format(table))
                with open(sql_file, 'r') as sql:
                    sh.mysql('-h', host, '-P', port, '-u', user, '-p{}'.format(password), database, _in=sql)
                    
                self.print_log('Importing data to `{}`.`{}` table.'.format(database, table))
                sh.mysqlimport('-h', host, '-P', port, '-u', user, '-p{}'.format(password), '--delete', '--local', '--ignore-lines=1', database, text_file)
                self.print_log('Imported data to `{}`.`{}` table.'.format(database, table))
    
                with open(done_file, 'w'):
                    pass
                self.print_log('Installed RefGene.')


################################################################
#                            RefSeq                            #
################################################################

    #TODO: provide for non-human species

        #creates refseq gene summary
        if 'all' in install or 'refseq' in install:
            done_file = os.path.join(refseq_dir,'.refseq.done')
            if my.file_exists(done_file):
                self.print_log('RefSeq already installed. To reinstall "rm {}"'.format(done_file))
            else:
                self.print_log('Installing RefSeq. This will take ~1 hour.')
    
                #download RefSeq directory
                versions['refseq'] = 'accessed {}'.format(strftime("%Y-%m-%d", localtime()))
                refseq_ftp = config.get('refseq', 'REFSEQ_FTP')
                refseqgene_ftp_dir = config.get('refseq', 'REFSEQGENE_FTP_DIR')
                refseqgene_pattern = config.get('refseq', 'REFSEQGENE_PATTERN')
                table = 'refseq_summary'
                text_file = os.path.join(refseq_dir, table+'.txt')
                sql_file = text_file+'.mysql'
                self.download_ftp_dir(refseq_ftp, refseqgene_ftp_dir, refseq_dir, refseqgene_pattern)
                self.print_log('Creating RefSeq tab-delimited file {}.'.format(text_file))
                refgene2db.run(dir=refseq_dir, output=text_file)
                self.print_log('Created RefSeq tab-delimited file {}.'.format(text_file))
                self.print_log('Making CREATE TABLE `{}` script.'.format(table))
                sql_columns.run(text_file, database=database, schema=database, table=table, indexes=['GeneSymbol'])
                self.print_log('Made CREATE TABLE `{}` script.'.format(table))
                self.print_log('Creating `{}` table.'.format(table))
                with open(sql_file, 'r') as sql:
                    sh.mysql('-h', host, '-P', port, '-u', user, '-p{}'.format(password), database, _in=sql)
                self.print_log('Importing data to `{}`.`{}` table.'.format(database, table))
                sh.mysqlimport('-h', host, '-P', port, '-u', user, '-p{}'.format(password), '--delete', '--local', '--ignore-lines=1', database, text_file)
                self.print_log('Imported data to `{}`.`{}` table.'.format(database, table))
                with MySQLdb.connect(host,user,password,database,port=port) as cursor:
                    proc = 'hgnc2refseq_summary'
                    self.print_log('Creating {} stored procedure'.format(proc))
                    query = [
'DROP PROCEDURE IF EXISTS {}'.format(proc),
"""CREATE DEFINER=CURRENT_USER PROCEDURE {}(hgnc VARCHAR(255))
BEGIN
SELECT DISTINCT r.`SUMMARY` AS summary
FROM `refseq_summary` AS r
WHERE r.`GeneSymbol` = hgnc
ORDER BY r.`SUMMARY`;
END
""".format(proc)]
                    self.mysql(cursor, query)
                    self.print_log('Created {} stored procedure'.format(proc))
                    
                with open(done_file, 'w'):
                    pass
                self.print_log('Installed RefSeq.')


################################################################
#            RGD (Rat Genome Database) ontologies              #
################################################################

        if 'all' in install or 'ontologies' in install or 'rgd' in install:
            done_file = os.path.join(rgd_dir,'.rgd.done')
            if my.file_exists(done_file):
                self.print_log('RGD ontologies already installed. To reinstall "rm {}"'.format(done_file))
            else:
                self.print_log('Installing RGD ontologies for {}. This will take ~25 minutes.'.format(' '.join(rgd_genera)))

                #download RGD directory
                rgd_ftp = config.get('rgd', 'RGD_FTP')
                rgd_ftp_dir = config.get('rgd', 'RGD_FTP_DIR')
                rgd_pattern = config.get('rgd', 'RGD_PATTERN')
                rgd_pattern = '|'.join([rgd_pattern.format(s) for s in rgd_genera])
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
                    if 'homo' in rgd_genera:
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
                self.print_log('Installed RGD ontologies for {}.'.format(' '.join(rgd_genera)))


################################################################
#                         VAX Plugins                          #
################################################################

        if 'all' in install or 'vax_plugins' in install:
            done_file = os.path.join(plugins_version_dir, '.vax_plugins.done')
            if my.file_exists(done_file):
                self.print_log('VAX plugins already installed. To reinstall "rm {}"'.format(done_file))
            else:
                #copy vax plugins to current version plugins dir
                vax_plugins_source_dir = config.get('vax', 'VAX_PLUGINS_SOURCE_DIR')
                vax_plugins_glob = vax_plugins_source_dir+'/*'
                self.print_log('Copying vax plugins from {} to {}.'.format(vax_plugins_source_dir, plugins_version_dir))
                sh.cp('-a', sh.glob(vax_plugins_glob), plugins_version_dir)
                
                with open(done_file, 'w'):
                    pass
                self.print_log('Copied vax plugins from {} to {}.'.format(vax_plugins_source_dir, plugins_version_dir))


################################################################
#                         1000 Genomes                         #
################################################################
    #TODO: add this?


################################################################
#                           UniProt                            #
################################################################

    #TODO: Update Swissknife
    #wget -N -P  https://sourceforge.net/projects/swissknife/files/latest/download/Swissknife_<version>.tar.gz
    #tar xf Swissknife_<version>.tar.gz
    ##EITHER (to install in perl libary)
    #cd Swissknife_<version)
    #perl Makefile.PL
    #make install
    #cd /share/apps/myourshaw/vax
    #if [ -L Swissknife ]; then rm Swissknife; fi
    #ln -s Swissknife_<version/SWISS Swissknife
    ##OR (if not doing make install)
    #if [ -L Swissknife ]; then rm Swissknife; fi
    #ln -s Swissknife_<version/lib Swissknife
    #note: make install may be unnecessary as long as PERL5LIB contains /share/apps/myourshaw/vax/Swissknife/lib

#TODO: convert to cluster job
        #required by KEGG
        if 'all' in install or 'uniprot' in install:
            done_file = os.path.join(uniprot_dir, '.uniprot.done')
            if my.file_exists(done_file):
                self.print_log('UniProt already installed. To reinstall "rm {}"'.format(done_file))
            else:
                self.print_log('Installing UniProt.')
                #download UniProt data
                versions['uniprot'] = 'UniProt accessed {}'.format(strftime("%Y-%m-%d", localtime()))
                uniprot_sprot_pattern = config.get('uniprot', 'UNIPROT_SPROT_PATTERN')
                uniprot_trembl_pattern = config.get('uniprot', 'UNIPROT_TREMBL_PATTERN')
                uniprot_idmapping_dat_pattern = config.get('uniprot', 'UNIPROT_IDMAPPING_DAT_PATTERN')
                uniprot_idmapping_selected_pattern = config.get('uniprot', 'UNIPROT_IDMAPPING_SELECTED_PATTERN')
                uniprot_files = []
                for t in uniprot_taxonomic_divisions:
                    uniprot_files += [uniprot_sprot_pattern.format(t), uniprot_trembl_pattern.format(t)]
                for o in uniprot_organisms:
                    uniprot_files += [uniprot_idmapping_dat_pattern.format(o), uniprot_idmapping_selected_pattern.format(o)]
                for file in uniprot_files:
                    with FTPHost(config.get('uniprot', 'UNIPROT_FTP'), 'anonymous', '') as h:
                        self.print_log('Downloading {}.'.format(file))
                        h.chdir(os.path.dirname(file))
                        h.download_if_newer(os.path.basename(file), os.path.join(uniprot_dir, os.path.basename(file)))
                        self.print_log('Downloaded {}.'.format(file))
    
                #create tab-delimited UniProt files
                vax_uniprot2db_done_file = my.done_file(os.path.join(uniprot_dir, 'vax_uniprot2db.pl'))
                if my.file_exists(vax_uniprot2db_done_file):
                    self.print_log('UniProt tab-delimited files already created. To recreate "rm {}"'.format(vax_uniprot2db_done_file))
                else:
                    self.print_log('Create UniProt tab-delimited files. This will take ~{} min'.format(20*len(uniprot_taxonomic_divisions)))
                    vax_uniprot2db_pl = os.path.join(os.path.dirname(__file__), 'vax_uniprot2db.pl')
                    perl(vax_uniprot2db_pl, '-dir', uniprot_dir)
                    with open(vax_uniprot2db_done_file, 'w'):
                        pass
                    self.print_log('Created UniProt tab-delimited files.')

                #create Uniprot table definitions and tables, and import UniProt data
                self.print_log('Creating Uniprot MySQL tables. This will take ~{} hr.'.format(1*len(uniprot_taxonomic_divisions)))
                for t in uniprot_taxonomic_divisions:
                    for table in['uniprot_{}_protein'.format(t),
                                 'uniprot_{}_protein_feature'.format(t),
                                 'uniprot_{}_xref'.format(t)]:
                        text_file = os.path.join(uniprot_dir, table+'.txt')
                        sql_file = text_file+'.mysql'
                        sql_done_file = my.done_file(sql_file)
                        if my.file_exists(sql_done_file):
                            self.print_log('MySQL script {} already created. To recreate "rm {}"'.format(sql_file, sql_done_file))
                        else:
                            self.print_log('Making CREATE TABLE `{}` script.'.format(table))
                            sql_columns.run(text_file, database=database, schema=database, table=table, indexes=['UniProtKB_AC'])
                            self.print_log('Made CREATE TABLE `{}` script.'.format(table))
                            self.print_log('Creating `{}` table on MySQL server.'.format(table))
                            with open(sql_file, 'r') as sql:
                                sh.mysql('-h', host, '-P', port, '-u', user, '-p{}'.format(password), database, _in=sql)
                            self.print_log('Created `{}` table on MySQL server.'.format(table))
                            self.print_log('Importing data to `{}`.`{}` table.'.format(database, table))
                            sh.mysqlimport('-h', host, '-P', port, '-u', user, '-p{}'.format(password), '--delete', '--local', '--ignore-lines=1', database, text_file)
                            self.print_log('Imported data to `{}`.`{}` table.'.format(database, table))

                    with MySQLdb.connect(host,user,password,database,port=port) as cursor:

                        #create enst_uniprot tables for vax access
                        self.print_log('Creating ens_uniprot_{} table. This will take ~10 minutes'.format(t))
                        query = [
'DROP TABLE IF EXISTS `enst_uniprot_{0}`'.format(t),
"""CREATE TABLE `enst_uniprot_{0}`
(ENST VARCHAR(18),
KEY IX_ENST_UniProtKB_AC_topic(`ENST`,`UniProtKB_AC`,`topic`),
KEY IX_ENS_topic (`ENST`,`topic`),
KEY IX_UniProtKB_AC_topic (`UniProtKB_AC`,`topic`))
ENGINE=MyISAM DEFAULT CHARSET=latin1
SELECT xe.`RESOURCE_IDENTIFIER` AS ENST, u.*
FROM `uniprot_{0}_protein` u
JOIN `uniprot_{0}_xref` xe
ON u.`UniProtKB_AC` = xe.`UniProtKB_AC` AND xe.`RESOURCE_ABBREVIATION` = 'Ensembl' AND xe.`RESOURCE_IDENTIFIER` LIKE 'ENS%'
""".format(t)]
                        self.mysql(cursor, query)

                        self.print_log('Creating enst_uniprot_{}_feature table. This will take ~5 min.'.format(t))
                        query = [
'DROP TABLE IF EXISTS `enst_uniprot_{0}_feature`'.format(t),
"""CREATE TABLE enst_uniprot_{0}_feature
(ENST VARCHAR(18),
KEY IX_ENS_feature_aaStart_aaEnd (ENST,feature,aaStart,aaEnd),
KEY IX_UniProtKB_AC_feature_aaStart_aaEnd (UniProtKB_AC,feature,aaStart,aaEnd))
ENGINE=MyISAM DEFAULT CHARSET=latin1
SELECT DISTINCT xe.RESOURCE_IDENTIFIER AS ENST, u.*
FROM uniprot_{0}_protein_feature u
JOIN uniprot_{0}_xref xe
ON u.UniProtKB_AC = xe.UniProtKB_AC AND xe.RESOURCE_ABBREVIATION = 'Ensembl' AND xe.RESOURCE_IDENTIFIER LIKE 'ENS%'
""".format(t)]
                        self.mysql(cursor, query)

                        #create UniProt stored procedures for vax access
                        self.print_log('Creating enst2uniprot_{0} stored procedure'.format(t))
                        query = [
'DROP PROCEDURE IF EXISTS enst2uniprot_{0}'.format(t),
"""CREATE DEFINER=CURRENT_USER PROCEDURE enst2uniprot_{0}(enst VARCHAR(18))
BEGIN
SELECT u.topic, u.value
FROM enst_uniprot_{0} AS u
WHERE u.ENST = enst
AND COALESCE(u.value, '') <> ''
ORDER BY u.topic;
END
""".format(t)]
                        self.mysql(cursor, query)

                        self.print_log('Creating enst2uniprot_{0}_feature stored procedure'.format(t))
                        query = [
'DROP PROCEDURE IF EXISTS enst2uniprot_{0}_feature'.format(t),
"""CREATE DEFINER=CURRENT_USER PROCEDURE enst2uniprot_{0}_feature(enst VARCHAR(18), aaStart INT, aaEnd INT)
BEGIN
SELECT u.feature, u.aaStart, u.aaEnd, u.description
FROM enst_uniprot_{0}_feature AS u
WHERE u.ENST = enst
AND ((aaStart BETWEEN u.aaStart and u.aaEnd)
OR  (aaEnd BETWEEN u.aaStart and u.aaEnd))
ORDER BY u.feature, u.aaStart, u.aaEnd;
END
""".format(t)]
                        self.mysql(cursor, query)
                    
                with open(done_file, 'w'):
                    pass
                self.print_log('Installed UniProt.')

################################################################
#                             KEGG                             #
################################################################

        #requires prior UniProt installation of the same species
        if 'all' in install or 'kegg' in install:
            done_file = os.path.join(kegg_dir,'.kegg.done')
            if my.file_exists(done_file):
                self.print_log('KEGG already installed. To reinstall "rm {}"'.format(done_file))
            else:
                self.print_log('Installing KEGG. This will take ~10 mins.')
    
                #download KEGG pathways
                versions['kegg'] = 'KEGG accessed {}'.format(strftime("%Y-%m-%d", localtime()))
                kegg_file = os.path.join(kegg_dir, 'kegg_gene_pathways.txt')
                self.print_log('Downloading KEGG data to {}.'.format(kegg_file))
                vax_download_kegg_data.run(kegg_file, kegg_species)
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
                    
                    #create KEGG table for vax access
                    if 'hsa' in kegg_species:
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

                        #create KEGG stored procedure for vax access
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
                    if 'mmu' in kegg_species:
                        table = 'ensmusg_kegg'
                        self.print_log('Creating {} table'.format(table))
                        query = [
'DROP TABLE IF EXISTS `{}`'.format(table),
"""CREATE TABLE `{}`
(PRIMARY KEY (ENSG,gene_id,path_id))
ENGINE=MyISAM DEFAULT CHARSET=latin1
SELECT DISTINCT CAST(xe.`OPTIONAL_INFORMATION_2` AS CHAR(18)) ENSG, k.`gene_id`, k.`path_id`, LTRIM(RTRIM(k.pathway)) pathway
FROM `kegg_gene_pathways` AS k
JOIN `uniprot_rodents_xref` AS x
ON k.`gene_id` = x.`RESOURCE_IDENTIFIER` AND x.`RESOURCE_ABBREVIATION` = 'KEGG'
JOIN `uniprot_rodents_xref` AS xe
ON x.`UniProtKB_AC` = xe.`UniProtKB_AC` AND xe.`RESOURCE_ABBREVIATION` = 'Ensembl' AND xe.`OPTIONAL_INFORMATION_2` LIKE 'ENSMUSG%'
WHERE IFNULL(LTRIM(RTRIM(k.`pathway`)), '') <> ''
""".format(table)]
                        self.mysql(cursor, query)

                        #create KEGG stored procedure for vax access
                        proc = 'ensmusg2kegg'
                        self.print_log('Creating {} stored procedure'.format(proc))
                        query = [
'DROP PROCEDURE IF EXISTS {}'.format(proc),
"""CREATE DEFINER=CURRENT_USER PROCEDURE {}(ensg CHAR(18))
BEGIN
SELECT DISTINCT k.`pathway`
FROM `ensmusg_kegg` AS k
WHERE k.`ENSG` = ensg
ORDER BY k.`path_id`;
END
""".format(proc)]
                        self.mysql(cursor, query)

                    if 'rno' in kegg_species:
                        table = 'ensrnog_kegg'
                        self.print_log('Creating {} table'.format(table))
                        query = [
'DROP TABLE IF EXISTS `{}`'.format(table),
"""CREATE TABLE `{}`
(PRIMARY KEY (ENSG,gene_id,path_id))
ENGINE=MyISAM DEFAULT CHARSET=latin1
SELECT DISTINCT CAST(xe.`OPTIONAL_INFORMATION_2` AS CHAR(18)) ENSG, k.`gene_id`, k.`path_id`, LTRIM(RTRIM(k.pathway)) pathway
FROM `kegg_gene_pathways` AS k
JOIN `uniprot_rodents_xref` AS x
ON k.`gene_id` = x.`RESOURCE_IDENTIFIER` AND x.`RESOURCE_ABBREVIATION` = 'KEGG'
JOIN `uniprot_human_xref` AS xe
ON x.`UniProtKB_AC` = xe.`UniProtKB_AC` AND xe.`RESOURCE_ABBREVIATION` = 'Ensembl' AND xe.`OPTIONAL_INFORMATION_2` LIKE 'ENSRNOG%'
WHERE IFNULL(LTRIM(RTRIM(k.`pathway`)), '') <> ''
""".format(table)]
                        self.mysql(cursor, query)

                        #create KEGG stored procedure for vax access
                        proc = 'ensrnog2kegg'
                        self.print_log('Creating {} stored procedure'.format(proc))
                        query = [
'DROP PROCEDURE IF EXISTS {}'.format(proc),
"""CREATE DEFINER=CURRENT_USER PROCEDURE {}(ensg CHAR(18))
BEGIN
SELECT DISTINCT k.`pathway`
FROM `ensrnog_kegg` AS k
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
#                           vep.ini                            #
################################################################

        if 'all' in install or 'vep_ini' in install:
            done_file = my.done_file(vep_ini)
            if my.file_exists(done_file):
                self.print_log('vep.ini already installed. To reinstall "rm {}"'.format(done_file))
            else:
                self.print_log('Installing VEP config as {}.'.format(vep_ini))
                #remove commas for use in plugin parameters
                versions_tmp = {}
                for k,v in versions.items():
                    versions_tmp[k] = v.replace(',', '_').replace(' ','_')
                    
                vep_ini_pattern_file = config.get('vax', 'VEP_INI_PATTERN_FILE')
                with open(vep_ini_pattern_file) as f:
                    vep_ini_pattern = f.read()
                vep_ini_contents = vep_ini_pattern.format(
                    BIOPERL_DIR = bioperl_live_dir,
                    BIOPERL_VERSION = versions_tmp.get('bioperl', 'BioPerl'),
                    CADD_SNP_FILE = os.path.join(cadd_dir, config.get('CADD', 'CADD_SNP_FILE_NAME')),
                    CADD_VERSION = versions_tmp.get('cadd', 'CADD'),
                    DBNSFP_GZ = dbnsfp_gz,
                    DBNSFP_VERSION = versions_tmp.get('dbnsfp', 'dbnsfp'),
                    DBSNP_VERSION = versions_tmp.get('dbsnp', 'dbSNP'),
                    DIR = vax_version_dir_link,
                    DIR_CACHE = vep_cache_dir_link,
                    DIR_PLUGINS = plugins_dir_link,
                    ENSEMBL_DIR = ensembl_api_dir_link,
                    ENSEMBL_HOST = host,
                    ENSEMBL_PASSWORD = ensembl_password,
                    ENSEMBL_PORT = port,
                    ENSEMBL_USER = ensembl_user,
                    ENSEMBL_VERSION = versions_tmp.get('ensembl', 'Ensembl'),
                    FATHMM_VERSION = versions_tmp.get('fathmm', 'FATHMM'),
                    GATK_VERSION = versions_tmp.get('gatk', 'GATK'),
                    HGMD_HOST = host,
                    HGMD_PASSWORD = hgmd_password,
                    HGMD_PLATFORM = 'mysql',
                    HGMD_PORT = port,
                    HGMD_DATABASE = hgmd_database,
                    HGMD_USER = hgmd_user,
                    HGMD_VERSION = versions_tmp.get('hgmd', 'HGMD'),
                    HMDB_VERSION = versions_tmp.get('HMDB', 'hmdb'),
                    HPA_VERSION = versions_tmp.get('hpa', 'HPA'),
                    KEGG_VERSION = versions_tmp.get('kegg', 'KEGG'),
                    MITOCARTA_VERSION = versions_tmp.get('mitocarta', 'MitoCarta'),
                    MGI_VERSION = versions_tmp.get('mgi', 'MGI'),
                    NHLBI_VERSION = versions_tmp.get('nhlbi', 'NHLBI-EVS'),
                    OMIM_VERSION = versions_tmp.get('omim', 'OMIM'),
                    PYTHON_EXECUTABLE = python,
                    RGD_VERSION = versions_tmp.get('rgd', 'RGD'),
                    REFGENE_VERSION = versions_tmp.get('refgene', 'RefGene'),
                    REFSEQ_VERSION = versions_tmp.get('refseq', 'RefSeq'),
                    UNIPROT_VERSION = versions_tmp.get('uniprot', 'UniProt'),
                    VAX_DATABASE = database,
                    VAX_HOST = host,
                    VAX_PASSWORD = vax_password,
                    VAX_PERL5LIB = ensembl_PERL5LIB,
                    VAX_PLATFORM = 'mysql',
                    VAX_PORT = port,
                    VAX_USER = vax_user,
                    VAX_VERSION = versions_tmp.get('vax', 'VAX'),
                )
                with open(vep_ini, 'w') as f:
                    f.write(vep_ini_contents)
                    
                with open(done_file, 'w'):
                    pass
                self.print_log('Installed VEP config as {}.'.format(vep_ini))

################################################################
#                           VAX test                           #
################################################################
        
        
################################################################
#                            Finis                             #
################################################################

        #update versions table
        if versions:
            my.set_key_value_dict(host, user, password, port, database, table='versions', key_col='module', value_col='version', data=versions)

        
        self.print_log('VAX installer done.')
        self.print_log("""To run vax:
    1. The user's PERL5LIB must include '{}'""")
        if self.all_jobs:
            self.print_log("""
    2. These cluster jobs must complete successfully: {}""".format(','.join(self.all_jobs), ensembl_PERL5LIB))



def main():

    #command line arguments
    parser = argparse.ArgumentParser(
        description = 'install or update VAX',
        epilog = """vax.vax_installer 1.0β1 ©2011-2014 Michael Yourshaw all rights reserved.
Critical information about databases, versions, and file locations must be up-to-date in the configutation file:
    {}""".format(os.path.join(os.path.dirname(__file__), config_file)))
    parser.add_argument('--install', nargs='*', default=['all'], type=str, metavar='MODULE', choices = install_choices,
        help='Modules to install. Choose one or more from "{}" (default: all)'.format(' '.join(install_choices)))
    #Ensembl versions
    parser.add_argument('--ensembl_version', '-e', type=int,
        help='Version of Ensembl API/VEP/cache to install. At present, the latest version will be installed regardless of this parameter. (default: current version; example: 77)')
    parser.add_argument('--delete_ensembl_version', type=str,
        help='old version of Ensembl databases to delete (default: ensembl_version - 2)')
    parser.add_argument('--ensembl_species', nargs='*', default=['homo_sapiens', 'mus_musculus'],
        help='list of Ensembl species to install (default: homo_sapiens mus_musculus)')
    #databases
    parser.add_argument('--host', '-H', type=str,
        help='MySQL database server hostname or ip address (default: config.db.HOST, example: cortex.local)')
    parser.add_argument('--port', '-P', type=int,
        help='MySQL database (default: config.db.PORT, example: 3306)')
    parser.add_argument('--user', '-u', type=str,
        help='MySQL database user with priveleges to install VAX (default: config.db.USER, example: sa)')
    parser.add_argument('--password', '-p', type=str,
        help='MySQL password for installing user (default: config.db.USER, example: None = enter silently at prompt)')
    parser.add_argument('--database', '-d', type=str,
        help='VAX MySQL database; EXISTING DATABASE WILL BE OVERWRITTEN! (default: config.vax.VAX_DATABASE_FORMAT_STRING_{ensembl_version}, example: vax_77)')
    parser.add_argument('--vax_user', '-U', type=str,
        help='user for runtime access to vax database (default: config.vax.VAX_USER, example: vax)')
    parser.add_argument('--vax_password', '-W', type=str,
        help='UNSECURE password for runtime access to vax database (default: config.vax.VAX_PW, example: vax)')
    parser.add_argument('--ensembl_user', type=str,
        help='user for runtime access to ensembl databases (default: config.ensembl.ENSEMBL_USER, example: ensembl)')
    parser.add_argument('--ensembl_password', type=str,
        help='UNSECURE password for runtime access to ensembl databases (default: config.ensembl.ENSEMBL_PW, example: ensembl)')
    parser.add_argument('--hgmd_user', '-g', type=str,
        help='user for runtime access to HGMD Pro database (default: config.hgmd.HGMD_USER, example: hgmd)')
    parser.add_argument('--hgmd_password', '-G', type=str,
        help='UNSECURE password for runtime access to HGMD Pro database (default: config.hgmd.HGMD_PW, example: hgmd)')
    parser.add_argument('--hgmd_version', type=str,
        help='HGMD version (example: 2014.1); used as subdirectory if hgmd_directory not specified (default: config.hgmd.HGMD_VERSION)')
    parser.add_argument('--hgmd_dir', type=str,
        help='directory into which HGMD files were previously downloaded (default: config.hgmd.HGMD_DIR/<hgmd_version>, example: /scratch1/vax/hgmd/2014.1)')
    parser.add_argument('--install_hgmd_schemas', nargs='*',
        help='List of all HGMD mySQL database schemas to be installed (default: config.hgmd.INSTALL_HGMD_SCHEMAS, example: hgmd_pro hgmd_snp hgmd_phenbase hgmd_views)')
    parser.add_argument('--hgmd_database', type=str,
        help='HGMD mySQL database schema used by vax (default: config.hgmd.HGMD_DATABASE, example: hgmd_pro)')
    parser.add_argument('--kegg_species', nargs='*', default=['hsa', 'mmu'],
        help='list of KEGG species pathways to install (default: hsa mmu)')
    parser.add_argument('--rgd_genera', nargs='*', default=['homo', 'mus'],
        help='list of Rat Genome Database genera to install (default: homo mus)')
    parser.add_argument('--uniprot_taxonomic_divisions', nargs='*', default=['human', 'rodents'],
        help='list of UniProt taxonomic divisions to install (default: human rodents)')
    parser.add_argument('--uniprot_organisms', nargs='*', default=['HUMAN_9606', 'MOUSE_10090'],
        help='list of UniProt organisms (NAME_number) to download (default: HUMAN_9606 MOUSE_10090)')
    #top level directories
    parser.add_argument('--vax_dir', '-v', type=str,
        help='top-level directory for VAX installation (default: config.vax.VAX_DIR)')
    parser.add_argument('--scratch_dir', '-s', type=str,
        help='scratch directory for temporary install files, cache, and dbnsfp (default: config.vax.SCRATCH_DIR)')
    #versions
    parser.add_argument('--bioperl_version', type=str,
        help='BioPerl version to install (example: release-1-6-923)')
    parser.add_argument('--cadd_snp_file_url', type=str,
        help='url of CADD snp file to install; large ~80Gb; (example: http://cadd.gs.washington.edu/download/v1.0/whole_genome_SNVs.tsv.gz)')
    parser.add_argument('--dbnsfp_file', type=str,
        help='dbNSFP file to install; large ~3Gb; (example: dbNSFPv2.4.zip)')
    parser.add_argument('--fathmm_file', type=str,
        help='FATHMM file to install; large ~80Gb; (example: fathmm.v2.3.SQL.gz)')
    parser.add_argument('--nhlbi_vcf_file_url', type=str,
        help='URL of NHLBI-EVS VCF file to install (example: http://evs.gs.washington.edu/evs_bulk_data/ESP6500SI-V2-SSA137.protein-hgvs-update.snps_indels.vcf.tar.gz)')
    #options
    parser.add_argument('--overwrite_vax_database', action='store_true', default=False,
        help='overwrite VAX database, deleting all data and procedures (default: if VAX database exists, just overwrite specific tables and procedures installed by this application')
    #parse args
    args = parser.parse_args()

    config = my.get_config(args, config_file)

    VI = VaxInstaller()
    VI.run(config=config, install=args.install,
           ensembl_version=args.ensembl_version, delete_ensembl_version=args.delete_ensembl_version, ensembl_species=args.ensembl_species,
           host=args.host, port=args.port, user=args.user, password=args.password, database=args.database,
           vax_user=args.vax_user, vax_password=args.vax_password, ensembl_user=args.ensembl_user, ensembl_password=args.ensembl_password,
           hgmd_user=args.hgmd_user, hgmd_password=args.hgmd_password, hgmd_version=args.hgmd_version, hgmd_dir=args.hgmd_dir, install_hgmd_schemas=args.install_hgmd_schemas, hgmd_database=args.hgmd_database,
           kegg_species=args.kegg_species, rgd_genera=args.rgd_genera,
           uniprot_taxonomic_divisions=args.uniprot_taxonomic_divisions, uniprot_organisms=args.uniprot_organisms,
           vax_dir=args.vax_dir, scratch_dir=args.scratch_dir,
           bioperl_version=args.bioperl_version, cadd_snp_file_url=args.cadd_snp_file_url,
           dbnsfp_file=args.dbnsfp_file, fathmm_file=args.fathmm_file, nhlbi_vcf_file_url=args.nhlbi_vcf_file_url,
           overwrite_vax_database=args.overwrite_vax_database,)

if __name__ == "__main__": sys.exit(main())

