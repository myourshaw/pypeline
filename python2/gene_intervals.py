#!/usr/bin/env python

#usage: python /home/myourshaw/local/bin/gene_intervals.py <path/to/results/directory>

import sys
import os
import re
import traceback
from optparse import OptionParser, OptionGroup
import MySQLdb

# Init cmd-line args
description = """
This script creates bed and interval_list files from ucsc gene tables.
Each row represents an exon of a gene transcript. exonStart and exonEnd are
expanded by 2 bases to include essential splice site locations.
"""
parser = OptionParser(description=description)
parser.add_option("-d", "--dir", metavar="DIR", dest="dir", help="Directory for output files")
(options, args) = parser.parse_args()

def tableExists(cursor,database,table):
	cursor.execute("""SELECT COUNT(*) as count FROM information_schema.tables WHERE table_schema = %s AND table_name = %s""",(database,table))
	return cursor.fetchone()['count'] > 0

outputDir = options.dir

dbs=["hg19","hg18"]
tables=['ccdsGene','ensGene','geneid','genscan','knownGene','mgcGenes','orfeomeGenes','refGene','sgpGene','vegaGene','wgEncodeGencode2wayConsPseudoV4','wgEncodeGencodeAutoV4','wgEncodeGencodeManualV4','wgEncodeGencodePolyaV4','xenoRefGene']
#dbs=["hg19"]
#tables=['knownGene']

for db in dbs:
	try:
		conn = MySQLdb.connect (host = "genome-mysql.cse.ucsc.edu", user = "genome", db = db)
		#cursor = conn.cursor ()
		cursor = conn.cursor (MySQLdb.cursors.DictCursor)
		for table in tables:
			if not tableExists(cursor,db,table):
				continue
			table = db+'.'+table
			txBedFile = open(os.path.join(outputDir,table+'.tx.bed'),"w")
			txIntervalFile = open(os.path.join(outputDir,table+'.tx.interval_list'),"w")
			cdsBedFile = open(os.path.join(outputDir,table+'.cds.bed'),"w")
			cdsIntervalFile = open(os.path.join(outputDir,table+'.cds.interval_list'),"w")
			cursor.execute ('SELECT * FROM ' + table + ' ORDER BY chrom,txStart,txEnd,name')
			while (1):
				row = cursor.fetchone ()
				if row == None:
					break
				#(name,chrom,strand,txStart,txEnd,cdsStart,cdsEnd,exonCount,exonStarts,exonEnds) = row
				chrom = row["chrom"]
				txStart = row["txStart"]
				txEnd =  row["txEnd"]
				cdsStart = row["cdsStart"]
				cdsEnd =  row["cdsEnd"]
				exonStarts = row["exonStarts"].strip(',').split(",")
				exonEnds = row["exonEnds"].strip(',').split(",")
				name = row["name"]
				name2 = ':'+row.get("name2") if row.get("name2") is not None and len(row.get("name2").strip()) > 0 else ""
				exonCount = row["exonCount"]
				for i in range(exonCount):
					exonNumber = i+1 if row["strand"] == '+' else row["exonCount"]-i
					exonStart = int(exonStarts[i])
					exonEnd = int(exonEnds[i])
					txFragStart0 = exonStart-2 if i > 0 else exonStart
					txFragStart1 = exonStart-1 if i > 0 else exonStart+1
					txFragEnd = exonEnd+2 if i < exonCount-1 else exonEnd
					txBedFile.write("%s\t%u\t%u\t%s\n" % (chrom,txFragStart0,txFragEnd,name+name2+'|exon'+str(exonNumber)))
					txIntervalFile.write("%s\t%u\t%u\t+\t%s\n" % (chrom,txFragStart1,txFragEnd,name+name2+'|exon'+str(exonNumber)))
					if exonEnd >= cdsStart and exonStart <= cdsEnd:
						if exonStart <= cdsStart and exonEnd >= cdsStart:
							cdsFragStart0 = cdsStart
							cdsFragStart1 = cdsStart+1
						else:
							cdsFragStart0 = txFragStart0
							cdsFragStart1 = txFragStart1
						if exonStart <= cdsEnd and exonEnd >= cdsEnd:
							cdsFragEnd = cdsEnd
						else:
							cdsFragEnd = txFragEnd
						cdsBedFile.write("%s\t%u\t%u\t%s\n" % (chrom,cdsFragStart0,cdsFragEnd,name+name2+'|exon'+str(exonNumber)))
						cdsIntervalFile.write("%s\t%u\t%u\t+\t%s\n" % (chrom,cdsFragStart1,cdsFragEnd,name+name2+'|exon'+str(exonNumber)))
			txBedFile.close()
			txIntervalFile.close()
			cdsBedFile.close()
			cdsIntervalFile.close()
		conn.close ()
	except MySQLdb.Error, e:
		print "Error %d: %s" % (e.args[0], e.args[1])
		sys.exit (1)
