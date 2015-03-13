#!/usr/bin/env python

import sys
import os
import collections
import glob
import re
import warnings
import pickle
import my


class QseqError(Exception): pass
class QseqIOError(QseqError): pass
class QseqFileNotFoundError(QseqIOError): pass
class BarcodeError(Exception): pass
class SequencingMetadataError(Exception): pass

OwnerId = collections.namedtuple('OwnerId', 'owner id')
class Owner(object):
	def __init__(self, id, email, samples = None, libraries = None, readgroups = None):
		self.id = id
		self.email = email
		self.samples = samples
		self.libraries = libraries
		self.readgroups = readgroups

class Owners(dict):
	def __init__(self, **kwds):
		super(Owners, self).__init__(**kwds)
		
StudyId = collections.namedtuple('StudyId', 'owner id')
class Study(object):
	def __init__(self, owner, id):
		self.owner = owner
		self.id = id

class Studies(dict):
	def __init__(self, **kwds):
		super(Studies, self).__init__(**kwds)
		
SampleId = collections.namedtuple('SampleId', 'owner id')
class Sample(object):
	def __init__(self, owner, id, libraries = None, readgroups = None, alignments = None, analyses = None, species = 'Homo sapiens'):
		self.owner = owner
		self.id = id
		self.libraries = libraries
		self.readgroups = readgroups
		self.alignments = alignments
		self.analyses = analyses
		self.species = species
		
class Samples(dict):
	def __init__(self, **kwds):
		super(Samples, self).__init__(**kwds)

LibraryId = collections.namedtuple('LibraryId', 'owner id')
class Library(object):
	def __init__(self, owner, id, sample, library_preparation_protocol = None, predicted_median_insert_size = None):
		self.owner = owner
		self.sample = sample
		self.id = id
		self.library_preparation_protocol = library_preparation_protocol
		self.predicted_median_insert_size = predicted_median_insert_size
	
class Libraries(dict):
	def __init__(self, **kwds):
		super(Libraries, self).__init__(**kwds)

ReadGroupId = collections.namedtuple('ReadGroupId', 'owner id')
class ReadGroup(object):
	def __init__(self, owner, id, sequencing_center = None, description = None, library = None,
							 platform = 'ILLUMINA', instrument_model = 'Illumina HiSeq 2000',
							 experiment_run_id = None, run_date = None, run_number = None, instrument = None,
							 lane = None, barcode = None, platform_unit = None, sample = None,
							 qseq_files = None, runs_output_directory = None, analysis_output_directory = None):
		self.owner = owner
		self.id = id
		self.sequencing_center = sequencing_center
		self.description = description
		self.library = library
		self.platform = platform
		self.instrument_model = instrument_model
		self.experiment_run_id = experiment_run_id
		self.run_date = run_date
		self.instrument = instrument
		self.run_number = run_number
		self.lane = lane
		self.barcode = barcode
		self.platform_unit = platform_unit
		self.sample = sample
		self.qseq_files = qseq_files
		self.runs_output_directory = runs_output_directory
		self.analysis_output_directory = analysis_output_directory

class ReadGroups(dict):
	def __init__(self, **kwds):
		super(ReadGroups,self).__init__(**kwds)

AlignmentId = collections.namedtuple('AlignmentId', 'owner id')
class Alignment(object):
	def __init__(self, id, owner):
		self.owner = owner
		self.id = id

class Alignments(dict):
	def __init__(self, **kwds):
		super(Alignments,self).__init__(**kwds)

AnalysisId = collections.namedtuple('AnalysisId', 'owner id')
class Analysis (object):
	def __init__(self, id, owner):
		self.owner = owner
		self.id = id

class Analyses(dict):
	def __init__(self, **kwds):
		super(Analyses,self).__init__(**kwds)

class BarcodeSet(set):
	def __init__(self, **kwds):
		super(BarcodeSet,self).__init__(**kwds)

class MultiplexQseqFiles(dict):
	def __init__(self, **kwds):
		super(MultiplexQseqFiles,self).__init__(**kwds)	
			
def table2metadata(path_to_metadata):
	
	owners = Owners()
	studies = Studies()
	samples = Samples()
	libraries = Libraries()
	readgroups = ReadGroups()
	alignments = Alignments()
	analyses = Analyses()
	multiplexqseqfiles = MultiplexQseqFiles()
	
	with open(path_to_metadata) as metadata:
		line_count = 0
		for line in metadata:
			line_count += 1
			if line.startswith('##') or line.strip() == '':
				continue
			elif line.startswith('#'):
				col_names = line.lstrip('#').rstrip('\n').lower().split('\t')
				#	sequencer_qseq_glob	analysis_output_directory					
				OWNER = col_names.index('owner')
				EMAIL = col_names.index('email')
				STUDY = col_names.index('study')
				SAMPLE = col_names.index('sample')
				SPECIES = col_names.index('species')
				LIBRARY = col_names.index('library')
				LIBRARY_PROTOCOL = col_names.index('library_protocol')
				LIBRARY_PREDICTED_MEDIAN_INSERT_SIZE = col_names.index('library_predicted_median_insert_size')
				SEQUENCING_LIBRARY_MIX = col_names.index('sequencing_library_mix')
				EXPERIMENT = col_names.index('experiment')
				SEQUENCING_CENTER = col_names.index('sequencing_center')
				PLATFORM = col_names.index('platform')
				INSTRUMENT_MODEL = col_names.index('instrument_model')
				RUN_DATE = col_names.index('run_date')
				MACHINE = col_names.index('machine')
				RUN_NUMBER = col_names.index('run_number')
				FLOWCELL = col_names.index('flowcell')
				LANE = col_names.index('lane')
				BARCODE = col_names.index('barcode')
				INDEX = col_names.index('index')
				EXPERIMENT_RUN_ID = col_names.index('experiment_run_id')
				READ_GROUP_ID = col_names.index('read_group_id')
				READ_GROUP_DESCRIPTION = col_names.index('read_group_description')
				SEQUENCER_QSEQ_DIR = col_names.index('sequencer_qseq_dir')
				RUNS_OUTPUT_DIRECTORY = col_names.index('runs_output_directory')
				ANALYSIS_OUTPUT_DIRECTORY = col_names.index('analysis_output_directory')
				continue
			else:
				data = line.rstrip('\n').split('\t')
				if data[EXPERIMENT_RUN_ID].find(data[FLOWCELL]) == -1:
					warnings.warn(\
																				'at line {} flowcell [{}] is not in experiment run id [{}]\n'\
																				.format(line_count, data[FLOWCELL], data[EXPERIMENT_RUN_ID]))
				qseq_dir = data[SEQUENCER_QSEQ_DIR]
				if data[SEQUENCER_QSEQ_DIR].find(data[EXPERIMENT_RUN_ID]) == -1:
					warnings.warn(\
																				'at line {} experiment run id [{}] is not in sequencer qseq dir [{}]\n'\
																				.format(line_count, data[EXPERIMENT_RUN_ID], data[SEQUENCER_QSEQ_GLOB]))
				qseq_lane_glob = sorted(glob.glob(os.path.join(data[SEQUENCER_QSEQ_DIR],'s_{}_?_*_qseq.txt'.format(data[LANE]))))
				if not qseq_lane_glob:
					raise SequencingMetadataError(\
																				'at line {} sequencer qseq dir [{}] has no matching files for lane [{}]\n'\
																				.format(line_count, data[SEQUENCER_QSEQ_DIR], data[LANE]))
				with open(qseq_lane_glob[0]) as qseq:
					qseq_machine, qseq_run, qseq_lane, qseq_tile, qseq_x, qseq_y, qseq_index, qseq_read, qseq_sequence, qseq_quality, qseq_filter = qseq.readline().rstrip('\n').split('\t')
				if qseq_machine != data[MACHINE]:
					raise SequencingMetadataError(\
																				'at line {} qseq machine [{}] is not same as metadata file machine [{}]\n'\
																				.format(line_count, qseq_machine, data[MACHINE]))
				if qseq_run != data[RUN_NUMBER]:
					raise SequencingMetadataError(\
																				'at line {} qseq run number [{}] is not same as metadata file run number [{}]\n'\
																				.format(line_count, qseq_run, data[RUN_NUMBER]))
				if qseq_lane != data[LANE]:
					raise SequencingMetadataError(\
																				'at line {} qseq lane [{}] is not same as metadata file lane [{}]\n'\
																				.format(line_count, qseq_lane, data[LANE]))
				owners[OwnerId(data[OWNER], data[EMAIL])] = Owner(data[OWNER],data[EMAIL])
				studies[StudyId(data[OWNER], data[STUDY])] = Study(data[OWNER],data[STUDY])
				samples[SampleId(data[OWNER], data[SAMPLE])] = Sample(data[OWNER],data[SAMPLE])
				libraryid = '{}.{}'.format(data[SAMPLE],data[LIBRARY])
				libraries[LibraryId(data[OWNER], libraryid)] = Library(data[OWNER],data[SAMPLE],data[LIBRARY])
				readgroupid = '{}.{}.{}.{}.{}{}'.format(data[RUN_DATE],qseq_machine,qseq_run,data[FLOWCELL],qseq_lane,'.'+data[BARCODE] if data[BARCODE] else '')
				readgroupownerid = ReadGroupId(data[OWNER], readgroupid)
				readgroups[readgroupownerid] = ReadGroup(data[OWNER],readgroupid,
					sequencing_center = data[SEQUENCING_CENTER],
					description = data[READ_GROUP_DESCRIPTION],
					library = data[LIBRARY],
					platform = data[PLATFORM],
					instrument_model = data[INSTRUMENT_MODEL],
					run_date = data[RUN_DATE],
					instrument = qseq_machine,
					run_number = qseq_run,
					experiment_run_id = data[EXPERIMENT_RUN_ID],
					lane = qseq_lane,
					barcode = data[BARCODE],
					platform_unit = '{}.{}{}'.format(data[EXPERIMENT_RUN_ID], data[LANE], '.'+data[BARCODE] if data[BARCODE] else ''),
					sample = data[SAMPLE],
					#qseq_files = qseq_lane_glob if not data[BARCODE] else None,
					qseq_files = qseqs2qseqtrios(qseq_lane_glob),
					runs_output_directory = os.path.join(data[RUNS_OUTPUT_DIRECTORY], data[OWNER], data[STUDY], 'runs', data[EXPERIMENT_RUN_ID]),
					analysis_output_directory = os.path.join(data[ANALYSIS_OUTPUT_DIRECTORY], data[OWNER], data[STUDY],'analysis')
					)
	barcodefilegroup = BarcodeFilegroup()
	barcodefilegroupkeys = set([BarcodeFilegroupKey(rg.owner, rg.experiment_run_id, rg.lane) for rg in readgroups.values()])
	for k in barcodefilegroupkeys:
		qseq_trios = [rg.qseq_files for rg in readgroups.values() if rg.owner == k.owner and rg.experiment_run_id == k.experiment_run_id and rg.lane == k.lane]
		s = QseqTrioSet()
		for q in qseq_trios:
			s = s | q
		qseq_trios = s
		barcodes = [(rg.barcode,rg.sample) for rg in readgroups.values() if rg.owner == k.owner and rg.experiment_run_id == k.experiment_run_id and rg.lane == k.lane]
		barcodesample = BarcodeSample()
		for b in barcodes:
			barcodesample[b[0]] = b[1]
		barcodefilegroupitem = BarcodeFilegroupItem(qseq_trios,barcodesample)
		barcodefilegroup[k] = barcodefilegroupitem
	return (owners, studies, samples, libraries, readgroups, barcodefilegroup)
	
def get_reads(read):
	#qseq.txt read -> (read1, read2, read3)
	dirname,basename,flowcell,lane,read,barcode,tile = split_qseq_path(read)
	return (join_qseq_path((dirname,lane,'1',tile)),
					join_qseq_path((dirname,lane,'2',tile)),
					join_qseq_path((dirname,lane,'3',tile)))

flowcell_re = re.compile(r"/?(?P<flowcell>\d{6}_[^_/]+_\d+(?:_[^_/]+)?_[^_/\.]+)/?", re.I)
qseq_re = re.compile(r"s_(?P<lane>\d)_(?P<read>\d+)_(?P<tile>\d+)_qseq\.txt", re.I)
merge_re = re.compile(r"\.(?P<lane>\d)\.(?P<read>\d+).*_qseq\.txt", re.I)
merge_barcode_re = re.compile(r"\.(?P<lane>\d)\.(?P<read>\d+)\.(?P<barcode>[ACGT]+).*_qseq\.txt", re.I)
read1_re = re.compile(r"s_(?P<lane>[1-8])_(?P<read>1)_(?P<tile>\d{4})_qseq\.txt", re.I)
read2_re = re.compile(r"s_(?P<lane>[1-8])_(?P<read>2)_(?P<tile>\d{4})_qseq\.txt", re.I)
read3_re = re.compile(r"s_(?P<lane>[1-8])_(?P<read>3)_(?P<tile>\d{4})_qseq\.txt", re.I)
merge1_re = re.compile(r"\.(?P<lane>[1-8])\.(?P<read>1).*_qseq\.txt", re.I)
merge2_re = re.compile(r"\.(?P<lane>[1-8])\.(?P<read>2).*_qseq\.txt", re.I)
merge3_re = re.compile(r"\.(?P<lane>[1-8])\.(?P<read>3).*_qseq\.txt", re.I)

def split_qseq_path(path):
	#<dirname=*<exptrun>*>/s_<lane>_<read>_<tile>_qseq.txt*-> (dirname,basename,flowcell,lane,read,None,tile)
	#<dirname>/<exptrun>.<lane>.<read>_qseq.txt*-> (dirname,basename,flowcell,lane,read,None,None)
	#<dirname>/<exptrun>.<lane>.<read>.<barcode>_qseq.txt*-> (dirname,basename,flowcell,lane,read,barcode,None)
	match = flowcell_re.search(path)
	dirname,basename = os.path.split(path)
	flowcell = dirname.strip('/').replace('/','.') if match == None else match.group('flowcell')
	match = qseq_re.search(basename)
	if match:
		return (dirname,basename,flowcell,match.group('lane'),match.group('read'),None,match.group('tile'))
	else:
		match = merge_re.search(basename)
		if match:
			return (dirname,basename,flowcell,match.group('lane'),match.group('read'),None,None)
		else:
			match = merge_barcode_re.search(basename)
			if match == None:
				return (dirname,basename,flowcell,match.group('lane'),match.group('read'),match.group('barcode'),None)
			else:
				print '{0} is not a recognized qseq file name'.format(path)
				sys.exit(1)

def join_qseq_path(p):
	#(dirname,lane,read,tile) -> <dirname>/s_<lane>_<read>_<tile>_qseq.txt
	foo = p[0]
	bar = p[1:]
	return os.path.join(p[0],'s_{0}_qseq.txt'.format('_'.join(p[1:])))
	
QseqTrio = collections.namedtuple('QseqTrio', 'read1 read2 read3')
class QseqTrioSet(set):
	def __init__(self, **kwds):
		super(QseqTrioSet,self).__init__(**kwds)
	
class BarcodeSample(dict): #key = barcode, value = sample | index
	def __init__(self, **kwds):
		super(BarcodeSample,self).__init__(**kwds)
	def __str__ (self):
		return '\n'.join(sorted(['{} {}'.format(b[1],b[0]) for b in self.items()]))
	@property
	def barcodes(self):
		return self.keys()
	@property
	def samples(self):
		return self.values()

class NovobarcodeList(BarcodeSample):
	def __init__(self, distance = 2, format = 'N 5', **kwds):
		self.distance = distance
		self.format = format
		super(NovobarcodeList,self).__init__(**kwds)
	
	def __str__ (self):
		return 'Distance {}\nFormat {}\n{}'.format(self.distance, self.format, str(super(NovobarcodeList,self).__str__()))

class BarcodeFilegroupItem(object):
	def __init__(self, qseq_trios = None, barcode_set = None, **kwds):
		self.qseq_trios = qseq_trios
		self.barcode_set = barcode_set

BarcodeFilegroupKey = collections.namedtuple('BarcodeFilegroupKey', 'owner experiment_run_id lane')
class BarcodeFilegroup(dict): #key = BarcodeFilegroupKey, value = BarcodeFilegroupItem
	def __init__(self, **kwds):
		super(BarcodeFilegroup,self).__init__(**kwds)

def main():
	owners, studies, samples, libraries, readgroups, barcodefilegroup = table2metadata('/home/myourshaw/lab/git-myourshaw/python/metadata.txt')
	with open(os.path.splitext('/home/myourshaw/lab/git-myourshaw/python/metadata.txt')[0]+'.pkl','wb') as metadata_pickle:
		pickle.dump(studies, metadata_pickle, -1)
		pickle.dump(samples, metadata_pickle, -1)
		pickle.dump(libraries, metadata_pickle, -1)
		pickle.dump(readgroups, metadata_pickle, -1)
		pickle.dump(barcodefilegroup, metadata_pickle, -1)
	print 'done'
		

if __name__ == "__main__": sys.exit(main())
