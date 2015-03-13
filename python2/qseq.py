#!/usr/bin/env python

import sys
import os
import re

class QseqError(Exception): pass
class QseqIOError(QseqError): pass
class QseqFileNotFoundError(QseqIOError): pass
class BarcodeError(Exception): pass

qseq_re = re.compile(r".*/(?P<date>\d+)_<?P<instrument>([^_/]+)_(?P<run>\d+)_(?P<flowcell>[^_/]+)/.*/s_(?P<lane>\d)_(?P<read>\d+)_(?P<tile>\d+)_qseq\.txt$")

#a qseq file
class QseqFile(object):
	def __init__(self, path, owner = None, study = None, sample = None, species = None,
							 library_protocol = None, capture_library = None, sequencing_library = None,
							 predicted_median_insert_size = None, experiment = None, platform = None,
							 instrument_model = None, run_date = None, qseq_machine = None,
							 qseq_run = None, qseq_lane = None, barcode = None, index = None,
							 read_group_id = None, sequencing_center = None, read_group_description = None,
							 experiment_run_storage_dir = None, experiment_run = None,
							 qseq_glob = None, qseq_merged_dir = None, qseq_merged = None,
							 bam_read_group_dir = None, bam_read_group = None, bam_sample_dir = None,
							 bam_sample = None):
		self.path = path
		path_match = qseq_re.match(self.path,re.I)
		if path_match:
			path_date = path_match.group('date')
			path_instrument = path_match.group('instrument')
			path_run = path_match.group('run')
			path_flowcell = path_match.group('flowcell')
			path_lane = path_match.group('lane')
			path_tile = path_match.group('tile')
			path_read = path_match.group('read')
		self.owner = owner,
		self.study = study,
		self.sample = sample,
		self.species = species if species else 'Homo sapiens',
		self.library_protocol = library_protocol
		self.capture_library = capture_library
		self.sequencing_library = sequencing_library
		self.predicted_median_insert_size = predicted_median_insert_size
		self.experiment_run = experiment_run if experiment_run else 
		self.experiment = experiment
		self.platform = platform if platform else 'ILLUMINA'
		self.instrument_model = instrument_model if instrument_model else 'Illumina HiSeq 2000'
		self.run_date = run_date
		self.qseq_machine = qseq_machine
		self.qseq_run = qseq_run
		self.qseq_lane = qseq_lane
		self.barcode = barcode
		self.index = index
		self.read_group_id = read_group_id
		self.sequencing_center = sequencing_center
		self.read_group_description = read_group_description
		self.experiment_run_storage_dir = experiment_run_storage_dir
		self.qseq_glob = qseq_glob
		self.qseq_merged_dir = qseq_merged_dir
		self.qseq_merged = qseq_merged
		self.bam_read_group_dir = bam_read_group_dir
		self.bam_read_group = bam_read_group
		self.bam_sample_dir = bam_sample_dir
		self.bam_sample = bam_sample	

	@property
	def basename(self):
		return(os.path.basename(self.path))

class BarcodeList(dict):
	def __init__(self,
							 #default is Illumina TruSeq/Agilent SureSelect XT barcode set
							 barcodes = {'ATCACGA': '01',
							'CGATGTA': '02',
							'TTAGGCA': '03',
							'TGACCAA': '04',
							'ACAGTGA': '05',
							'GCCAATA': '06',
							'CAGATCA': '07',
							'ACTTGAA': '08',
							'GATCAGA': '09',
							'TAGCTTA': '10',
							'GGCTACA': '11',
							'CTTGTAA': '12'
							},
							 **kwds
	):
		self.__barcodes = barcodes
		super(BarcodeList,self).__init__(barcodes, **kwds)
		
	def __getitem__(self, barcode):
		assert isinstance(barcode, basestring), "barcode must be a string"
		return super(BarcodeList,self).__getitem__(barcode)
		
	def __setitem__(self, barcode, library):
		assert isinstance(barcode, basestring), "barcode must be a string"
		assert isinstance(library, basestring), "library must be a string"
		if re.search(r"[^ACGTacgt]", barcode, re.I):
			raise BarcodeError('barcodes must contain only ACGT')
		else:
			super(BarcodeList,self).__setitem__(barcode.upper(), library)
	
	def __str__ (self):
		return '\n'.join(sorted(['{} {}'.format(b[1],b[0]) for b in self.items()]))

	@property
	def barcodes(self):
		return self.items
	
class NovobarcodeList(BarcodeList):
	def __init__(self, distance = 2, format = 'N 5', **kwds):
		self.distance = distance
		self.format = format
		super(NovobarcodeList,self).__init__(**kwds)
	
	def __str__ (self):
		return 'Distance {}\nFormat {}\n{}'.format(self.distance, self.format, str(super(NovobarcodeList,self).__str__()))

#read_1_file[,read_2_file][,read_barcode_file, barcode_list]
#1-3 QseqFile objects that have the same set of (experiment run lane tile x y) in the same order
#with a BarcodeList object if the set contains a qseq
class QseqPairedFiles(object):
	def __init__(self,
							 read_1_file = None,
							 read_2_file = None,
							 read_barcode_file = None,
							 barcode_list = None
							):
		self.read_1_file = None
		self.read_2_file = None
		self.read_barcode_file = None
		self.barcode_list = None
	
	
