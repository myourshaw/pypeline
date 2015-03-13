#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import argparse
from ConfigParser import SafeConfigParser
from contextlib import closing
from getpass import getpass
import re
try:
    import MySQLdb
    import MySQLdb.cursors
except ImportError:
    raise ImportError,"The Python MySQLdb module is required to run this program. Try 'pip install MySQL-python'."
try:
    import sh
except ImportError:
    raise ImportError,"The Python sh module is required to run this program. Try 'pip install sh'."

#these files must be in same directory as this script
import my
import sql_columns
import vax_dbsnp_vcfs2db

class DbsnpInstallerError(Exception):
    """Log vax_installer-generated exceptions"""
    def __init__(self,msg, log=True):
        self.msg = msg
        if log:
            print_log(msg)
        
class DbsnpInstaller:
    """Install dbSNP MySQL database file"""
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

#--dbsnp_file /scratch1/vax/75/dbsnp/All.vcf.gz --dbsnp_version b141_GRCh37p13 --proc get_dbsnp -d vax_test -H cortex.local -P 3306 -u sa -p 
#--dbsnp_file /scratch1/vax/75/dbsnp/clinvar_20140702.vcf.gz --dbsnp_version b141_GRCh37p13 --proc get_dbsnp_clinvar -d vax_test -H cortex.local -P 3306 -u sa -p 
#--dbsnp_file /scratch1/vax/75/dbsnp/common_no_known_medical_impact_20140702.vcf.gz --dbsnp_version b141_GRCh37p13 --proc get_dbsnp_common_no_medical_impact -d vax_test -H cortex.local -P 3306 -u sa -p 

    def run(self, dbsnp_file, dbsnp_version, database, proc=None, config=None, host=None, port=None, user=None, password=None,):
        
################################################################
#                            dbSNP                             #
################################################################

        #flatten ALT and expand INFO in dbSNP VCF
        table = 'dbsnp_'+dbsnp_version+'_'+my.r_strip(os.path.basename(dbsnp_file), '.vcf.gz')
        done_file = my.done_file(os.path.join(os.path.dirname(dbsnp_file), table))
        if my.file_exists(done_file):
            self.print_log('{} already imported. To reimport "rm {}"'.format(table, done_file))
        else:
            self.print_log('Importing {} ~6 hours for All'.format(table))
            text_file = os.path.join(os.path.dirname(dbsnp_file), table+'.txt')
            text_file_done = my.done_file(text_file)
            if my.file_exists(text_file_done):
                self.print_log('{} already created. To recreate "rm {}"'.format(text_file, text_file_done))
            else:
                self.print_log('Creating {} ~1 hour for All'.format(text_file))
                vax_dbsnp_vcfs2db.run(input=dbsnp_file, output=text_file, uncompressed=True)
                with open(text_file_done, 'w'):
                    pass
                self.print_log('Created {}'.format(text_file))
            
            #create mySQL table
            sql_file = text_file+'.mysql'
            sql_done_file = my.done_file(sql_file)
            if my.file_exists(sql_done_file):
                self.print_log('{} already created. To redo "rm {}"'.format(sql_file, sql_done_file))
            else:
                self.print_log('Making CREATE TABLE `{}` script. ~5 hours for All'.format(table))
                sql_columns.run(text_file, database=database, schema=database, table=table, clustered_index=['CHROM', 'POS'], indexes=['id'])
                with open(sql_done_file, 'w'):
                    pass
                self.print_log('Made CREATE TABLE `{}` script.'.format(table))
                
            self.print_log('Creating `{}` table.'.format(table))
            with open(sql_file, 'r') as sql:
                sh.mysql('-h', host, '-P', port, '-u', user, '-p{}'.format(password), database, _in=sql)
            
            #import table data
            self.print_log('Importing data to `{}`.`{}` table ~30 minutes for All.'.format(database, table))
            sh.mysqlimport('-h', host, '-P', port, '-u', user, '-p{}'.format(password), '--delete', '--local', '--ignore-lines=1', database, text_file)
            self.print_log('Imported data to `{}`.`{}` table.'.format(database, table))
    
            #create stored procedure
            if not proc:
                proc = 'get_'.format(table)
            self.print_log('Creating `{}`.`{}` stored procedure'.format(database, proc))
            query = [
'DROP PROCEDURE IF EXISTS {}'.format(proc),
"""CREATE DEFINER=CURRENT_USER PROCEDURE {}(CHROM VARCHAR(2), POS INT(11), REF VARCHAR(255), ALT VARCHAR(255))
BEGIN
SELECT *
FROM `{}` AS d
WHERE d.`CHROM` = CHROM
AND d.`POS` = POS
AND d.`REF` = REF
AND d.`ALT` = ALT;
END
""".format(proc, table)]
            with MySQLdb.connect(host,user,password,database,port=port) as cursor:
                self.mysql(cursor, query)
            self.print_log('Created `{}`.`{}` stored procedure'.format(database, proc))

        with open(sql_done_file, 'w'):
            pass
        self.print_log('Imported {}'.format(table))


def main():

    #command line arguments
    parser = argparse.ArgumentParser(
        description = 'install dbSNP files table on mysql server',
        epilog = 'vax.dbsnp_installer 1.0β1 ©2011-2014 Michael Yourshaw all rights reserved')
    #required
    parser.add_argument('--dbsnp_file', '-i', type=str, required=True,
        help='dbSNP .vcf.gz file downloaded from ncbi')
    parser.add_argument('--dbsnp_version', '-v', type=str, required=True,
        help='dbSNP version for naming table (example: b141_GRCh37p13)')
    parser.add_argument('--database', '-d', type=str, required=True,
        help='MySQL database')
    parser.add_argument('--proc', type=str,
        help='name of stored procedure that will be used to access this table')
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
                    raise DbsnpInstallerError('A user is required. Installer aborted.')
                else:
                    #empty user
                    if not args.user:
                        raise DbsnpInstallerError('A user is required. Installer aborted.')
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
                    raise DbsnpInstallerError('A password is required. Installer aborted.')
                else:
                    #empty password
                    if not args.password:
                        raise DbsnpInstallerError('A password is required. Installer aborted.')
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
                    raise DbsnpInstallerError('A host is required. Installer aborted.')
                else:
                    #empty host
                    if not args.host:
                        raise DbsnpInstallerError('A host is required. Installer aborted.')
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
                    raise DbsnpInstallerError('A port is required. Installer aborted.')
                else:
                    #empty port
                    if not args.port:
                        raise DbsnpInstallerError('A port is required. Installer aborted.')    

    DbsnpInstaller().run(dbsnp_file=args.dbsnp_file, dbsnp_version=args.dbsnp_version, database=args.database, proc=args.proc,
        password=args.password, config=config, host=args.host, port=args.port, user=args.user,)

if __name__ == "__main__": sys.exit(main())

