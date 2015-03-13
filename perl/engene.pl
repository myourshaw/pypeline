#!/usr/bin/perl -w

=head1 LICENSE

  Copyright (c) 2007-2011 Michael Yourshaw All rights reserved.

=head1 CONTACT

  Please email comments or questions to <myourshaw@ucla.edu>.

=cut

=head1 NAME

Genomic Interval Lister - a script to list genomic intervals of genes nad transcripts

Version 1.0.beta.1

by Michael Yourshaw (myourshaw@ucla.edu)
=cut

use strict;
use warnings;
use Getopt::Long;
use FileHandle;
use File::Spec;
use Date::Calc qw(Today_and_Now Delta_YMDHMS);
use POSIX qw(ceil floor);
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
#BioPerl
use Bio::Seq;
use Bio::SeqIO;
#EnsEMBL Perl API
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor;
use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(MAX_DISTANCE_FROM_TRANSCRIPT);
use Bio::EnsEMBL::Variation::Utils::Sequence qw(unambiguity_code sequence_with_ambiguity);

# configure from command line opts
#my $config = &configure(scalar @ARGV);
my $config;

# run the main sub routine
#&main($config);

# this is the main sub-routine - it needs the configured $config hash
sub main {
	my $config = shift;
	
	debug("Starting...") if defined $config->{verbose};
	
	my $species = $config->{species};
	
	# get adaptors
	my $vfa = $config->{reg}->get_adaptor($species, 'variation', 'variationfeature');
	$config->{tva} = $config->{reg}->get_adaptor($species, 'variation', 'transcriptvariation');
	
	# get fake ones for species with no var DB
	if(!defined($vfa)) {
		$vfa = Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor->new_fake($species);
	}
	else {
		$vfa->db->include_failed_variations($config->{include_failed}) if defined($vfa->db) && $vfa->db->can('include_failed_variations');
	}
	
	
	$config->{sa} = $config->{reg}->get_adaptor($species, 'core', 'slice');
	$config->{ga} = $config->{reg}->get_adaptor($species, 'core', 'gene');
	
	# check we got slice adaptor - can't continue without a core DB
	die("ERROR: Could not connect to core database\n") unless defined $config->{sa} and defined $config->{ga};
	
	
	# create a hash to hold slices so we don't get the same one twice
	my %slice_hash = ();
	my @new_vfs;
	my %vf_hash;
	
	my $transcript_cache;
	
	my $line_number = 0;
	my $vf_count;
	my $in_file_handle = $config->{in_file_handle};
	
	# read the file
	while(<$in_file_handle>) {
	  chomp;
	  
	  $line_number++;
	  
	  # header line?
	  next if /^\#/;

	} # while(<$in_file_handle>)
} # sub main

#####################################################
#change header to geneSymbol|chr01:123456-987654|+|gene_stable_id|e01-10+20,e02-10+10|fragment_types|roi_coordinates|gene_biotype|gene_description
#standard fasta header options
#opional header contents and order
#####################################################
use constant HEADER_FORMAT => "header format:\n"
	. "geneSymbol|gene_stable_id|length|parent gene|Ensembl stable id|region of interest coordinates|region of interest length|padding-upstream+downstream|fragment type(s)|gene biotype|gene description\n"
	. "example:\n"
	. ">SEZ6_e1-35+20|chr17:24356939-24357324[-]|386b|SEZ6|ENSG00000063015|roi:24356939-24357324|386b|pad-0+0|CDS,UTR|protein_coding|seizure related 6 homolog isoform 1 [Source:RefSeq_peptide;Acc:NP_849191]\n"
	. "fragment id = gene name or region coordinates[_up|_e#[_5utr|_3utr|_utr]|_i#|_down ][-additional upstream roi][+additional downstream roi][_fragment part/total fragment parts][_rc]\n"
	. "up=upstream, e=exon, utr=untranslated region, i=intron, down=downstream, rc=reverse complement\n"
	. "fragment types = REGION|INTERGENETIC|GENE|UPSTREAM|CDS[,UTR]|UTR|INTRON|DOWNSTREAM|NONCODING\n"
	. "biotype = C_segment|D_segment|J_segment|miRNA|miRNA_pseudogene|misc_RNA|misc_RNA_pseudogene|Mt_rRNA|Mt_tRNA_pseudogene|protein_coding|pseudogene|repeat|retrotransposed|rRNA|rRNA_pseudogene|scRNA|scRNA_pseudogene|snoRNA|snoRNA_pseudogene|snRNA|snRNA_pseudogene|tRNA_pseudogene|V_segment\n"
;

use constant BIG_NUMBER => 9999999999;
#the following must have unique values
use constant CDS => "cds";
use constant DETAIL => "detail";
use constant DOWNSTREAM => "down";
use constant ERROR => "error";
use constant EXON => "e";
use constant GENE => "gene";
use constant INTERGENETIC => "intergenetic";
use constant INFO => "info";
use constant INTRON => "i";
use constant INPUT => "input";
use constant NONCODING => "noncoding";
use constant REGION => "region";
use constant SUMMARY => "region";
use constant TRANSCRIPT => "tr";
use constant UPSTREAM => "up";
use constant WARNING => "warning";
use constant UPDATE => "update";
use constant UTR3 => "3utr";
use constant UTR5 => "5utr";


my @start_time = Today_and_Now();

#parameter keys must be lower case
my %params = (
	add_bases_downstream_to_region_of_interest_to_get_minimum_length => 0, #bool
	add_bases_upstream_to_region_of_interest_to_get_minimum_length => 0, #bool
	bases_region_of_interest_3_prime_downstream_of_exons => 2,
	bases_region_of_interest_3_prime_downstream_of_introns => 0,
	bases_region_of_interest_5_prime_upstream_of_exons => 2,
	bases_region_of_interest_5_prime_upstream_of_introns => 0,
	bases_region_of_interest_gene_downstream_region => 0,
	bases_region_of_interest_gene_upstream_promoter_region => 0,
	create_log_file => 1, #bool (default = 1)
	create_genomic_intervals_file => 1, #bool (default = 1)
	create_sequence_file => 1, #bool (default = 1)
	create_summary_file => 1, #bool (default = 1)
	debug => 0, #bool (default = 0)
	exclude_region_of_interest_less_than_bases => 0,
	filespec_log_file => '', #default will be <job file>.log
	filespec_genomic_intervals_file => '', #default will be <job file>_intervals.txt
	filespec_sequence_file => '', #default will be <job file>.seq
	filespec_summary_file => '', #default will be <job file>.summary
	hard_mask_repeats => 0, #bool (default = 0)
	include_3_prime_utr_exons => 0, #bool
	include_5_prime_utr_exons => 0, #bool
	include_all_elements_one_record_per_gene => 0, #bool overrides all includes & separates
	include_cds_exons => 1, #bool
	include_gene_downstream_region => 0,
	include_gene_upstream_promoter_region => 0,
	include_intergenetic_sequence_in_genomic_regions => 0, #bool
	include_introns => 0, #bool
	include_only_these_biotypes => "", #list separated by non-alphanumeric character(s)
	include_rna_like_strand => 0, #bool
	include_template_strand => 1, #bool
	input_format => "", # '' (default) | omim_result gene_result
	mask_repeats => 0, #bool
	maximum_record_length => 1000000000000000, #split records if necessary
	maximum_region_of_interest_per_record => 1000000000000000,
	merge_transcripts => 1, #bool default 1
	minimum_record_length => 0, #pad records as necessary
	minimum_region_of_interest_per_record => 0,
	padding_bases_3_prime_downstream_of_each_record => 0,
	padding_bases_3_prime_downstream_of_each_region_of_interest => 0,
	padding_bases_5_prime_upstream_of_each_record => 0,
	padding_bases_5_prime_upstream_of_each_region_of_interest => 0,
	separate_records_for_exons_introns => 1, #bool
	separate_records_for_utr_cds => 0, #bool ignored unless separate_records_for_exons_introns
	separate_records_for_promoter_upstream_downstream_regions => 1, #bool
	sequence_file_format => "fasta", #fasta (default) | piecemaker | Bio::SeqIO formats
	show_variations => 0, #bool
	single_fasta_identifier => 0, #bool
	slice_format => 'chr%s:%u-%u[%s]', #e.g., 'chr%s:%u-%u[%s]'
	stop => 0, #bool forces a premature stop of reading the input file. for debugging.
	use_ensembl_human => 1,
	use_ensembl_vega => 1
	
	#ensembl_version => "" # default "" = latest
	#host => 'ensembldb.ensembl.org' # default = 'ensembldb.ensembl.org'
	#user => 'anonymous' # default = 'anonymous'
	#password => '' # default = no password
	#include_essential_splice_site => 1 #bool
	#regulatory_region_near_gene_start
	#regulatory_region_near_gene_start_or_end
	#create_annotation_of_input_file
	
);

#statistics
my %all_biotypes;# = (C_segment=>0,D_segment=>0,J_segment=>0,miRNA=>0,miRNA_pseudogene=>0,misc_RNA=>0,misc_RNA_pseudogene=>0,Mt_rRNA=>0,Mt_tRNA_pseudogene=>0,protein_coding=>0,pseudogene=>0,repeat=>0,retrotransposed=>0,rRNA=>0,rRNA_pseudogene=>0,scRNA=>0,scRNA_pseudogene=>0,snoRNA=>0,snoRNA_pseudogene=>0,snRNA=>0,snRNA_pseudogene=>0,tRNA_pseudogene=>0,V_segment=>0);
my %all_lengths;
my $max_roi = 0;
my $min_roi = BIG_NUMBER;
#my @roi_lengths;
my %genes;
my %counts;

#globals
my $id;
my $gene_name;
my $that_id;
my @all_coding;
my @all_exons;
my @all_transcripts;
my @all_regions;

use constant USAGE => "engene usage:\n"
	 . "perl engene.pl prints usage.\n"
	 . "perl engene.pl [-S] [-g genomic_intervals_file | -G] [-l log_file | -L] [-d detailed_log_file | -D] [-Q] input_file_list [> sequence_file]\n"
	 . "gets information from Ensembl and creates a sequence file (default format is fasta),\n"
	 . "a genomic intervals file (Agilent eArray genomic tiling format) and a log file.\n"
	 . "input_file_list = - | space-separated list of files; - = input from STDIN, ctrl-D when done.\n"
	 . "Input files consist of a list of gene names (e.g., MYO1D or ENSG00000063015)\n"
	 . "or genomic regions (e.g., chr17:27945334-28228015 +), which will be expanded\n"
	 . "to include all of any genes that overlap the region's end points.\n"
	 . "Comments may follow the gene/region identifier.\n"
	 . "Input files may optionally include program settings (e.g., \$input_format = omim_result).\n"
	 . "Input files are processed sequentially, line-by-line;\n"
	 . "setting changes only apply to subsequent genes/regions.\n"
	 . "Input files may optionally include comments (e.g., # foo bar).\n"
	 . "sequence_file = sequence file path (default is STDOUT; - is input_file1.seq)\n"
	 . "-S = no sequence file output\n"
	 . "-g genomic_intervals_file = genomic intervals file path or STDOUT (default is input_file1.txt with . replaced by _)\n"
	 . "-G = no genomic intervals file output\n"
	 . "-l log_file = specify log file path (default is input_file1.summary and STDERR)\n"
	 . "-L = log to STDERR only\n"
	 . "-d detailed_log_file = specify log file path (default is input_file1.log and STDERR)\n"
	 . "-D = no detailed log file\n"
	 . "-Q = quiet mode, no output to STDERR\n"
	 . "WARNING: sequence, genomic intervals and log files are created by default\n"
	 . "and overwrite existing files of the same name.\n"
;

#options
my %opts; #sequence, genomic intervals, summary, log
getopts("s:Sg:Gm:Ml:LQ", \%opts);
#files
unless (@ARGV) {
	print_usage();
	exit;
}
my $input_file = "";
my $input_file1 = $ARGV[0];
my ($vol,$dir,$file);
if ( $input_file1 eq "-" ) {
	print STDERR "Input genes/genomic regions or program settings, one per line. Ctrl-D when done\n" unless $opts{Q};
	($vol,$dir,$file) = (File::Spec->curdir()."engene_output");
}
else{
	($vol,$dir,$file) = File::Spec->splitpath($input_file1);
}
my $sequence_file;
#remove special characters from file names
(my $clean_file = $file) =~ s/\W/_/;
#redirect STDOUT
my $gt = "&gt;"; #komodo saves command line > as &gt;
if($opts{s} || $params{create_sequence_file} || $params{filespec_sequence_file} || !$opts{S}){
	if($ARGV[$#ARGV] eq "-"){
		$sequence_file = File::Spec->catpath($vol,$dir,"$file.seq");
		pop(@ARGV);
	}
	elsif(index($ARGV[$#ARGV],">") == 0){
		$sequence_file = substr pop(@ARGV),1;
	}
	elsif ($ARGV[$#ARGV-1] eq ">" || $ARGV[$#ARGV-1] eq $gt) {
		$sequence_file = pop(@ARGV);
		pop(@ARGV);
	}
	elsif (rindex($ARGV[$#ARGV-1],">") == length($ARGV[$#ARGV-1])-1) {
		$sequence_file = pop(@ARGV);
		$ARGV[$#ARGV-1] = substr($ARGV[$#ARGV-1],0,rindex($ARGV[$#ARGV-1],">"));
	}
	elsif (rindex($ARGV[$#ARGV-1],$gt) == length($ARGV[$#ARGV-1])-1) {
		$sequence_file = pop(@ARGV);
		$ARGV[$#ARGV-1] = substr($ARGV[$#ARGV-1],0,rindex($ARGV[$#ARGV-1],$gt));
	}
	elsif($opts{s}){
		$sequence_file = $opts{s};
	}
	elsif($params{filespec_sequence_file}){
		$sequence_file = $params{filespec_sequence_file};
	}
	else{
		$sequence_file = File::Spec->catpath($vol,$dir,"$file.seq");
	}
	while($ARGV[$#ARGV] eq ">" || $ARGV[$#ARGV] eq $gt){
		pop(@ARGV);
	}
	open STDOUT, ">", $sequence_file or die "can't open $sequence_file: $!";
}
else{$opts{S}=1;}
my $genomic_intervals_file;
if($opts{g} || $params{create_genomic_intervals_file} || $params{filespec_genomic_intervals_file} || !$opts{G}){
	if($opts{g}){
		$genomic_intervals_file = $opts{g};
	}
	elsif($params{filespec_genomic_intervals_file}){
		$genomic_intervals_file = $params{filespec_genomic_intervals_file};
	}
	else{
		$genomic_intervals_file = File::Spec->catpath($vol,$dir, $clean_file."_intervals.txt");
	}
	open GENOMICINTERVALSFILE, ">", $genomic_intervals_file or die "Can't open $genomic_intervals_file: $!";
}
else{$opts{G}=1;}
my $log_file;
if($opts{d} || $params{create_log_file} || $params{filespec_log_file} || !$opts{D}){
	if($opts{d}){
		$log_file = $opts{d};
	}
	elsif($params{filespec_log_file}){
		$log_file = $params{filespec_log_file};
	}
	else{
		$log_file = File::Spec->catpath($vol,$dir,"$file.log");
	}
	open LOGFILE, ">", $log_file or die "Can't open $log_file: $!";
}
else{$opts{D}=1;}
my $summary_file;
if($opts{l} || $params{create_summary_file} || $params{filespec_summary_file} || !$opts{L}){
	if($opts{l}){
		$log_file = $opts{l};
	}
	if($params{filespec_summary_file}){
		$summary_file = $params{filespec_summary_file};
	}
	else{
		$summary_file = File::Spec->catpath($vol,$dir,"$file.summary");
	}
	open SUMMARYFILE, ">", $summary_file or die "Can't open $summary_file: $!";
}
else{$opts{L}=1;}
#start
print_log(sprintf("START: %s", datetimef(@start_time)),INFO);
print_log("Default settings:",INFO);
print LOGFILE map { "\t\$$_ = $params{$_}\n" } sort keys %params unless $opts{D};
print STDERR map { "\t\$$_ = $params{$_}\n" } sort keys %params unless $opts{Q};
print SUMMARYFILE map { "\t\$$_ = $params{$_}\n" } sort keys %params unless $opts{L};
print_log(sprintf("sequence file: %s",$sequence_file ? $sequence_file : "STDOUT"),INFO);
print_log(sprintf("genomic intervals file: %s",$genomic_intervals_file),INFO);
print_log(sprintf("detailed log file: %s",$log_file),INFO);
print_log(sprintf("summary log file: %s",$summary_file),INFO);
#connect to EnsEMBL database
my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db(
	-host => 'ensembldb.ensembl.org',
	-user => 'anonymous');
my $exon_adaptor = $registry->get_adaptor('human', 'core', 'exon');
my $gene_adaptor = $registry->get_adaptor('human', 'core', 'gene');
my $vega_adaptor = $registry->get_adaptor('human', 'vega', 'gene');
my $repeat_adaptor = $registry->get_adaptor('human', 'core', 'repeatfeature');
my $slice_adaptor = $registry->get_adaptor('human', 'core', 'slice');
my $transcript_adaptor = $registry->get_adaptor('human', 'core', 'transcript');
my @db_adaptors = @{ $registry->get_all_DBAdaptors() };
#read from a list of filenames on the command line
#read from SDTIN if no filenames were given
#"-" indicates STDIN
#"someprogram |" indicates the output of another program (7.14)
#read and parse input file
while (<>) {
	unless ($params{stop}){
		$counts{'input line'}++;
		if ($ARGV ne $input_file) {
			print_log(sprintf("input file: %s", File::Spec->rel2abs($ARGV)),INFO);
			$input_file = $ARGV;
			$counts{'input file'}++;
		}
		chomp;
		if($_){
			print_log($_,INPUT);
			undef $id;
			undef @all_coding;
			undef @all_exons;
			undef @all_transcripts;
			#comment
			if (/^\s*#/) { }
			#parameter
			elsif (/^\s*\$/) {
				parameter();
			}
			#omim_results file
			elsif (lc($params{input_format}) eq "omim_result") {
				if (/\s*[*+#%]?\d{6}.*;\s*(\S+)/) {
					my $this_id = $1;
					if (defined($that_id)) {
						get_by_gene_name($that_id);
						$counts{'input total record'}++;
						$counts{'input gene (OMIM) record'}++;
					}
					$that_id = $this_id;
				}
				elsif ((/\s*Gene map/i || !$_) && defined($that_id)) {
					get_by_gene_name($that_id);
					$counts{'input total record'}++;
					$counts{'input gene (OMIM) record'}++;
					undef $that_id;
				}
				elsif (eof && defined($that_id)) {
					if(defined($that_id)){
						get_by_gene_name($that_id);
						$counts{'input total record'}++;
						$counts{'input gene (OMIM) record'}++;
					}
					undef $that_id;
				}
			}
			#NCBI gene_results file
			elsif (lc($params{input_format}) eq "gene_result") {
				if (/\s*\d+:\s*(\S+)/) {
					my $this_id = $1;
					if (defined($that_id)) {
						get_by_gene_name($that_id);
						$counts{'input total record'}++;
						$counts{'input gene (ncbi) record'}++;
					}
					$that_id = $this_id;
				}
				elsif ((/\s*GeneID/i || !$_) && defined($that_id)) {
					get_by_gene_name($that_id);
					$counts{'input total record'}++;
					$counts{'input gene (ncbi) record'}++;
					undef $that_id;
				}
				elsif (eof && defined($that_id)) {
					get_by_gene_name($that_id);
					$counts{'input total record'}++;
					$counts{'input gene (ncbi) record'}++;
					undef $that_id;
				}
			}
			#Chromosome position
			elsif (/^\s*(?:chr)?[\s]*([0-9XYMT]{1,2})[\s:]+([,\d]+)\D+([,\d]+)(?:\s+([01+-]{1}))?/i){
				if ($1 && $2 && $3) {
					my $chr = $1;
					my $start = $2;
					my $end = $3;
					my $strand = defined($4)
					 && $4 == "+" ? 1 : defined($4)
					 && $4 == "-" ? -1 : $4;
					$start =~ s/,//g;
					$end =~ s/,//g;
					$start = $start <= $end ? $start : $end;
					$end = $end >= $start ? $end : $start;
					get_by_region($chr, $start, $end, $strand);
					$counts{'input total record'}++;
					$counts{'input chromosome position record'}++;
				}
			}
			#Chromosome
			elsif (/^\s*(?:chr)?[\s]*([0-9XYMT]{1,2})[\s:]+/i){
				if ($1) {
					my $chr = $1;
					get_by_region($chr);
					$counts{'input total record'}++;
					$counts{'input chromosome record'}++;
				}
			}
			#Ensmbl gene stable id
			elsif (/^\s*(ENSG\d{11})/) {
				get_by_gene_id($1);
				$counts{'input total record'}++;
				$counts{'input gene record'}++;
				$counts{'input gene (ensembl) record'}++;
			}
			#Ensmbl transcript stable id
			elsif (/^\s*(ENST\d{11})/) {
				get_by_transcript_id($1);
				$counts{'input total record'}++;
				$counts{'input transcript record'}++;
			}
			#Ensmbl exon stable id
			elsif (/^\s*(ENSE\d{11})/) {
				get_by_exon_id($1);
				$counts{'input total record'}++;
				$counts{'input exon record'}++;
			}
			#Gene external name
			elsif (/^\s*([^#\s]+).*#*/) {
				get_by_gene_name($1);
				$counts{'input total record'}++;
				$counts{'input gene record'}++;
				$counts{'input gene (name) record'}++;
			}
			#unrecognized
			else {
				print_log(sprintf("%s has an unexpected format; not processed", $_),WARNING);
			}
		}
	}
}
#all input records done

#genomic interval output
print_all_regions() unless $opts{G};

#print statistics
print_log("lengths of regions of interest (count):",INFO);
print LOGFILE map { "$_ ($all_lengths{$_})," } sort {$a <=> $b} keys %all_lengths unless $opts{D};
print LOGFILE "\n" unless $opts{D};
print SUMMARYFILE map { "$_ ($all_lengths{$_})," } sort {$a <=> $b} keys %all_lengths unless $opts{L};
print SUMMARYFILE "\n" unless $opts{L};
print STDERR map { "$_ ($all_lengths{$_})," } sort {$a <=> $b} keys %all_lengths  unless $opts{Q};
print STDERR "\n" unless $opts{Q};
print_log("gene\tstableID\tbasesROI\tbasesROI+padded",INFO);
print LOGFILE map {"$genes{$_}[2]\t$_\t$genes{$_}[0]\t$genes{$_}[1]\n"} sort {$genes{$a}[2] cmp $genes{$b}[2]} keys %genes unless $opts{D};
print SUMMARYFILE map {"$genes{$_}[2]\t$_\t$genes{$_}[0]\t$genes{$_}[1]\n"} sort {$genes{$a}[2] cmp $genes{$b}[2]} keys %genes unless $opts{L};
print STDERR map {"$genes{$_}[2]\t$_\t$genes{$_}[0]\t$genes{$_}[1]\n"} sort {$genes{$a}[2] cmp $genes{$b}[2]} keys %genes unless $opts{Q};
print_log("gene biotypes:",INFO);
print LOGFILE map { " $all_biotypes{$_} $_\n" } sort keys %all_biotypes unless $opts{D};
print SUMMARYFILE map { " $all_biotypes{$_} $_\n" } sort keys %all_biotypes unless $opts{L};
print STDERR map { " $all_biotypes{$_} $_\n" } sort keys %all_biotypes unless $opts{Q};
$counts{'ROI minimum base'} = $min_roi;
$counts{'ROI maximum base'} = $max_roi;
print_log("statistics:",INFO);
print LOGFILE map { sprintf "%10d %s%s\n",$counts{$_},$_,ss($counts{$_}) } sort keys %counts unless $opts{D};
print SUMMARYFILE map { sprintf "%10d %s%s\n",$counts{$_},$_,ss($counts{$_}) } sort keys %counts unless $opts{L};
print STDERR map { sprintf "%10d %s%s\n",$counts{$_},$_,ss($counts{$_}) } sort keys %counts unless $opts{Q};
my @end_time = Today_and_Now();
print_log(sprintf("job duration %s",durationf(Delta_YMDHMS(@start_time,@end_time))),INFO);
print_log(sprintf("DONE: %s",datetimef(@end_time)),INFO);
exit;

sub get_by_exon_id {
	$id = shift;
	print_log($id,UPDATE);
	my $exon = $exon_adaptor->fetch_by_stable_id($id);
	if ($exon) {
		#coding or noncoding exon
		$counts{'exon'}++;
		print_slice(
			$id,
			$id,
			$id,
			$exon->slice->seq_region_name,
			$exon->start,
			$exon->end,
			$exon->strand,
			0,
			0,
			EXON);
	}
	else {
		print_log(sprintf("%s isn't an exon", $id), ERROR);
	}
} # !get_by_exon_id

sub get_by_gene_id {
	$id = shift;
	undef $gene_name;
	print_log($id,UPDATE);
	my $gene = $gene_adaptor->fetch_by_stable_id($id);
	if ($gene) {
		undef $gene_name;
		gene($gene);
	}
	else {
		print_log(sprintf("%s isn't a gene", $id), ERROR);
	}
} # !get_by_gene_id

sub get_by_gene_name {
	$id = shift;
	$gene_name = $id;
	my @genes;
	if($params{use_ensembl_human}){
		@genes = @{ $gene_adaptor->fetch_all_by_external_name($id) };
	}
	unless (@genes){
		if($params{use_ensembl_vega}){
			@genes = @{ $vega_adaptor->fetch_all_by_external_name($id) };
		}
	}
	if (@genes) {
		my $gene = $genes[0];
		if (@genes > 1) {
			foreach my $g (@genes) {
				if ($g->external_name eq $id) {
					$gene = $g;
					last;
				}
			}
			my $count_temp = @genes;
			my $external_name = $gene->external_name;
			my $stable_id = $gene->stable_id;
			print_log(sprintf("%s identifies %u genes; using %s (%s)",
					$id,
					$count_temp,
					$external_name,
					$stable_id),
				WARNING);
			foreach my $g (@genes) {
				print_log(sprintf("\t%s %s %s (%s)",
					defined($g->external_name) ? $g->external_name : "",
					slicef($gene->slice),
					defined($g->description) ? $g->description : "",
					$g->stable_id),INFO);
			}
		}
		if ($gene->external_name ne $id) {
			print_log(
				sprintf(
					"Ensembl's external name for %s is %s (%s %s %s)",
					$id,
					defined($gene->external_name) ? $gene->external_name : "",
					slicef($gene->slice),
					defined($gene->description) ? $gene->description : "",
					$gene->stable_id
				),
				WARNING
			);
		}
		gene($gene);
	}
	else {
		print_log(sprintf("%s isn't a known gene", $id), ERROR);
	}
	undef $gene_name;
} # !get_by_gene_name

#genomic region (chromosome start end [strand])
sub get_by_region {
	my ($chr, $start, $end, $strand) = @_;
	print_log(sprintf($params{slice_format}, $chr, $start, $end, defined($strand)?$strand:''),UPDATE);
	if (!defined($strand)) {
		$strand = 1;
	}
	my $slice;
	if(defined($start) && defined($end)){
		$slice = $slice_adaptor->fetch_by_region('chromosome', $chr, $start, $end, $strand);
	}
	else{
		$slice = $slice_adaptor->fetch_by_region('chromosome', $chr);
	}
	if ($slice) {
		$id =  slicef($slice);
		my @genes = @{ $slice->get_all_Genes };
		if(@genes || $params{include_intergenetic_sequence_in_genomic_regions}){
			if($params{include_intergenetic_sequence_in_genomic_regions} && ($params{separate_records_for_exons_introns} || @genes == 0)) {
				print_slice(
					$id,
					"",
					"",
					$chr,
					$start,
					$end,
					$strand,
					0,
					0,
					REGION);
			}
			elsif($params{include_intergenetic_sequence_in_genomic_regions}){
				my $ig_start = $start;
				my $last = my $ig = 1;
				foreach my $gene (@genes) {
					#get gene not relative to slice
					$gene = $gene_adaptor->fetch_by_stable_id($gene->stable_id);
					if ($gene->start > $ig_start) {
						print_slice(
							$id . "_".INTERGENETIC.$ig,
							"",
							"",
							$chr,
							$ig_start,
							$gene->start - 1,
							$strand,
							0,
							0,
							INTERGENETIC);
						$ig_start = $gene->end + 1;
						$ig++;
					}
					undef $gene_name;
					gene($gene);
				}
				#intergenetic sequence after last gene, if any
				if ($ig_start <= $end) {
					print_slice(
						$id . "_".INTERGENETIC.$ig,
						"",
						"",
						$chr,
						$ig_start,
						$end,
						$strand,
						0,
						0,
						INTERGENETIC);
				}
			}
			else {
				foreach my $gene (@genes) {
					#get gene not relative to slice
					undef $gene_name;
					gene($gene_adaptor->fetch_by_stable_id($gene->stable_id));
				}
			}
		}
	}
	else {
		print_log(sprintf("%s isn't a region", $id), ERROR);
	}
} # !get_by_region

sub get_by_transcript_id {
	$id = shift;
	print_log($id,UPDATE);
	my $transcript = $transcript_adaptor->fetch_by_stable_id($id);
	if ($transcript) {
		my $gene = $gene_adaptor->fetch_by_transcript_stable_id($id);
		transcript($gene, $transcript);
		print_gene_regions($gene);
	}
	else {
		print_log(sprintf("%s isn't a transcript", $id), ERROR);
	}
} # !get_by_transcript_id

sub gene {
	my $gene = shift;
	unless(defined($gene_name)){
		$gene_name = defined($gene->external_name)?$gene->external_name:$gene->stable_id;
	}
	print_log($gene_name,UPDATE);
    $genes{$gene->stable_id} = [0,0,$gene_name];
    $counts{'output gene'}++;
    $all_biotypes{$gene->biotype}++;
    my @exons = @{ $gene->get_all_Exons };
    print_log(sprintf("%s has no exons", $gene_name), WARNING) unless @exons;
    my $gene_strand = $exons[0]->strand or 1;
    my $gene_extra_up;
    my $gene_extra_down;
    if($params{bases_region_of_interest_gene_upstream_promoter_region}>$params{bases_region_of_interest_5_prime_upstream_of_exons}){
        $gene_extra_up = $params{bases_region_of_interest_gene_upstream_promoter_region};
    }
    else{
        $gene_extra_up = $params{bases_region_of_interest_5_prime_upstream_of_exons};
    }
    if($params{bases_region_of_interest_gene_downstream_region}>$params{bases_region_of_interest_3_prime_downstream_of_exons}
       &&$params{include_gene_downstream_region}){
        $gene_extra_down = $params{bases_region_of_interest_gene_downstream_region};
    }
    else{
        $gene_extra_down = $params{bases_region_of_interest_3_prime_downstream_of_exons};
    }
    if($params{include_gene_upstream_promoter_region}
       &&$params{separate_records_for_promoter_upstream_downstream_regions}&&$gene_extra_up){
		print $gene->start;
        print_slice(
            $gene_name.UPSTREAM,
            $gene_name,
            $gene->stable_id,
            $gene->slice->seq_region_name,
            $gene_strand==1?$gene->start-$gene_extra_up:$gene->end+1,
            $gene_strand==1?$gene->start-1:$gene->end+$gene_extra_up,
            $gene_strand,
            0,
            0,
            UPSTREAM,
            $gene->biotype,
            $gene->description);
    }
    my @transcripts = @{$gene->get_all_Transcripts};
    print_log(sprintf("%s is a pseudogene (no transcripts)", $id), WARNING) unless @transcripts;
    #print gene as one record
    if($params{include_all_elements_one_record_per_gene}
       ||($params{include_3_prime_utr_exons}
       &&$params{include_5_prime_utr_exons}
       &&$params{include_cds_exons}
       &&$params{include_introns}
       &&!$params{separate_records_for_exons_introns})) {
        print_slice(
            $gene_name,
            $gene_name,
            $gene->stable_id,
            $gene->slice->seq_region_name,
            $gene->start,
            $gene->end,
            $gene_strand,
            $params{include_gene_upstream_promoter_region}&&$params{separate_records_for_promoter_upstream_downstream_regions}?0:$gene_extra_up,
            $params{include_gene_downstream_region}&&$params{separate_records_for_promoter_upstream_downstream_regions}?0:$gene_extra_down,
            GENE,
            $gene->biotype,
            $gene->description);
    }
    #print gene as more than one record
    else{
        if (@transcripts) {
            undef(@all_transcripts);
            undef(@all_exons);
            undef(@all_coding);
            #get arrays of transcripts and of their coding regions
            foreach my $transcript (@transcripts) {
                push @all_transcripts, $transcript;
                transcript($gene, $transcript);
            }
            print_gene_regions($gene);
        }
    }
    if($params{include_gene_downstream_region}
       &&$params{separate_records_for_promoter_upstream_downstream_regions}&&$gene_extra_down){
        print_slice(
            $gene_name.DOWNSTREAM,
            DOWNSTREAM,
            $gene->stable_id,
            $gene->slice->seq_region_name,
            $gene_strand==1?$gene->end+1:$gene->start-$gene_extra_down,
            $gene_strand==1?$gene->end+$gene_extra_down:$gene->start-1,
            $gene_strand,
            0,
            0,
            DOWNSTREAM,
            $gene->biotype,
            $gene->description);
    }
} # !gene

sub transcript {
	my ($gene, $transcript) = @_;
	$counts{'analyzed transcript'}++;
	my @exons = @{ $transcript->get_all_Exons() };
	if (@exons) {
		foreach my $exon (@exons) {
			push @all_exons, new Bio::Range(
				-start => $exon->start,
				-end => $exon->end,
				-strand => $exon->strand);
			my $coding_start = $exon->coding_region_start($transcript);
			my $coding_end = $exon->coding_region_end($transcript);
			if (defined($coding_start) && defined($coding_end)) {
				push @all_coding, new Bio::Range(
					-start => $coding_start,
					-end => $coding_end,
					-strand => $exon->strand);
			}
		}
	}
	else {
		print_log(sprintf("%s has no exons", $id), WARNING);
	}
} # !transcript

sub bases {
	my($start,$end) = @_;
	return $end-$start+1;
} # !bases

#returns true if array contains string
sub contains{
	my ($array,$string) = @_;
	for my $array_string(@$array){
		if($string eq $array_string){
			return 1;
		}
	}
	return 0;
}

sub datetimef{
	return sprintf "%04d-%02d-%02d %02d:%02d:%02d", $_[0],$_[1],$_[2],$_[3],$_[4],$_[5];
}

sub durationf{
	my $years = $_[0]>0 ? sprintf("%u year%s ",$_[0],ss($_[0])):"";
	my $months = $_[1]>0 ? sprintf("%u year%s ",$_[1],ss($_[1])):"";
	my $days = $_[2]>0 ? sprintf("%u year%s ",$_[2],ss($_[2])):"";
	return sprintf "%s%s%s%02d:%02d:%02d",$years,$months,$days,$_[3],$_[4],$_[5];
}

sub fetch_repeat_masked_sequence_with_ambiguity {
	my ($slice,$soft) = @_;
	my $aseq = fetch_sequence_with_ambiguity($slice->seq_region_name, $slice->start, $slice->end, $slice->strand);
	my $rseq = $slice->get_repeatmasked_seq()->soft_mask($soft)->seq();
#	my $repeat_masked_slice = $slice->get_repeatmasked_seq();
#	$repeat_masked_slice->soft_mask($soft);
#	my $rseq = $repeat_masked_slice->seq();

} # !fetch_repeat_masked_sequence_with_ambiguity

sub fetch_sequence_with_ambiguity {
	my ($chr, $start, $end, $strand) = @_;
	my $registry = 'Bio::EnsEMBL::Registry';
	$registry->load_registry_from_db(
		-host => 'ensembldb.ensembl.org',
		-user => 'anonymous');
	my $dbCore = $registry->get_DBAdaptor('human', 'core');
	my $dbVar = $registry->get_DBAdaptor('human', 'variation');
	my $ambiguous_slice = sequence_with_ambiguity($dbCore, $dbVar, $chr, $start, $end, $strand);
	my $ambiguous_seq = $ambiguous_slice->seq();
	return $ambiguous_seq;
} # !fetch_sequence_with_ambiguity

sub add_upstream{
	my ($pos,$bases,$strand) = @_;
	if($strand==1){
		return $pos+$bases;
	}
	else{
		return $pos-$bases;
	}
}

sub add_downstream{
	my ($pos,$bases,$strand) = @_;
	if($strand==1){
		return $pos-$bases;
	}
	else{
		return $pos+$bases;
	}
}

sub is_downstream{
	my ($pos1,$pos2,$strand) = @_;
	if($strand==1&&$pos1>$pos2){
		return 1;
	}
	elsif($strand==-1&&$pos1<$pos2){
		return 1;
	}
	else{
		return 0;
	}
}

sub is_upstream{
	my ($pos1,$pos2,$strand) = @_;
	if($strand==1&&$pos1<$pos2){
		return 1;
	}
	elsif($strand==-1&&$pos1>$pos2){
		return 1;
	}
	else{
		return 0;
	}
}

#change parameter ($key value)
sub parameter {
	if(/\s*\$([^\s=>#]+){1}[\s=>#]+([^#]+)#*$/){
		my $key = lc($1);
		my $value = $2;
		$value =~ s/['"]+//g;
		if(defined($value)){
			if ($value =~ /^false$/i || $value =~ /^no$/i) {
				$value = 0;
			}
			elsif ($value =~ /^true$/i || $value =~ /^yes$/i) {
				$value = 1;
			}
		}
		if (exists($params{$key})) {
			print_log(sprintf("%s changed from %s to %s\n",$key,$params{$key},$value),INFO) unless $params{$key} eq $value;
			$params{$key} = $value;
		}
		else {
			print_log(sprintf("%s isn't a setting", $key),WARNING);
		}
	}
} # !parameter

sub print_all_regions {
	my @all_chr_regions;
	@all_regions = sort {$a->{num_chr} <=> $b->{num_chr} || $a->{start} <=> $b->{start} || $a->{end} <=> $b->{end}} @all_regions;
	my $that_chr = 0;
	my $this_chr = 0;
	foreach my $region(@all_regions){
		$this_chr = $region->{num_chr};
		if ($this_chr!=$that_chr){
			if(@all_chr_regions){
				my @union_chr_regions = Bio::Range->unions(@all_chr_regions);
				@union_chr_regions = sort {$a->{start} <=> $b->{start} || $a->{end} <=> $b->{end}} @union_chr_regions;
				foreach my $chr_region(@union_chr_regions) {
					printf GENOMICINTERVALSFILE "chr%s:%u-%u\n", numchr2chr($that_chr), $chr_region->start, $chr_region->end;
				}
				undef @all_chr_regions;
			}
			$that_chr = $this_chr;
		}
		push @all_chr_regions, new Bio::Range(
		-start => $region->{start},
		-end => $region->{end},
		-strand => $region->{strand});
	}
	#last region
	if(@all_chr_regions){
		my @union_chr_regions = Bio::Range->unions(@all_chr_regions);
		@union_chr_regions = sort {$a->{start} <=> $b->{start} || $a->{end} <=> $b->{end}} @union_chr_regions;
		foreach my $chr_region(@union_chr_regions) {
			printf GENOMICINTERVALSFILE "chr%s:%u-%u\n", numchr2chr($that_chr), $chr_region->start, $chr_region->end;
		}
	}
} # !print_all_regions

sub upmost{ #strand,pos
	return pop()==1?min(@_):max(@_);
}
sub downmost{ #strand,pos
	return pop()==1?max(@_):min(@_);
}
sub get_bio_range{
	my ($pos1,$pos2,$strand)=@_;
	return  new Bio::Range(
		-start => $pos1<$pos2?$pos1:$pos2,
		-end => $pos1>$pos2?$pos1:$pos2,
		-strand => $strand)
}
sub get_5prime_3prime{
	my $bio_range = shift;
	return $bio_range->strand ? ($bio_range->start,$bio_range->end) : ($bio_range->end,$bio_range->start);
}
sub print_gene_regions { #exons and introns in separate records
	my $gene = shift;
	unless (defined($gene_name)){
	    $gene_name = defined($gene->external_name)?$gene->external_name:$gene->stable_id;
	}
    #merge intersecting exons and intersecting coding regions
    print_log(sprintf("%s has no exons", $gene_name), WARNING) unless @all_exons;
    print_log(sprintf("%s has no coding regions", $gene_name), WARNING) unless @all_coding;
	my @union_exons = Bio::Range->unions(@all_exons);
    my @union_coding;
	if(@all_coding){@union_coding = Bio::Range->unions(@all_coding);}
	my $strand = $union_exons[0]->strand;
	my $up = $strand == 1 ? -1 : 1;
	my $down = $strand == 1 ? 1 : -1;
	my $gene_5prime;
	my $gene_3prime;
	my $exons_5prime;
	my $exons_3prime;
	my $coding_5prime;
	my $coding_3prime;
    my $utr5_3prime;
    my $utr3_5prime;
	if ($strand==1) {
		if(@union_exons){@union_exons = sort { $a->start <=> $b->start || $a->end <=> $b->end } @union_exons;}
		if(@union_coding){@union_coding = sort { $a->start <=> $b->start || $a->end <=> $b->end } @union_coding;}
		$gene_5prime = $gene->start;
		$gene_3prime = $gene->end;
		if(@union_exons){$exons_5prime = $union_exons[0]->start;}
		if(@union_exons){$exons_3prime = $union_exons[$#union_exons]->end;}
		if(@union_coding){$coding_5prime = $union_coding[0]->start;}
		if(@union_coding){$coding_3prime = $union_coding[$#union_coding]->end;}
	}
	else {
		if(@union_exons){@union_exons = sort { $b->start <=> $a->start || $b->end <=> $a->end } @union_exons;}
		if(@union_coding){@union_coding = sort { $b->start <=> $a->start || $b->end <=> $a->end } @union_coding;}
		$gene_5prime = $gene->end;
		$gene_3prime = $gene->start;
		if(@union_exons){$exons_5prime = $union_exons[0]->end;}
		if(@union_exons){$exons_3prime = $union_exons[$#union_exons]->start;}
		if(@union_coding){$coding_5prime = $union_coding[0]->end;}
		if(@union_coding){$coding_3prime = $union_coding[$#union_coding]->start;}
	}
    if(@union_coding){$utr5_3prime = $coding_5prime+$up*1;}
    if(@union_coding){$utr3_5prime = $coding_3prime+$down*1;}
	my $upstream_region_5prime = upmost($gene_5prime,$exons_5prime+$up*$params{bases_region_of_interest_gene_upstream_promoter_region},$strand);
	my $upstream_region_3prime = $exons_5prime + $up*1;
	my $downstream_region_5prime = $exons_3prime + $down*1;
	my $downstream_region_3prime = downmost($gene_3prime,$exons_3prime+$down*$params{bases_region_of_interest_gene_downstream_region},$strand);
#	#make an array of all the ranges in the gene that might become separate records
#	my @print_infos;
	my $that_exon_3prime;
	my $exon_number=0;
	foreach my $exon(@union_exons){
		my $exon_5prime = upmost($exon->start,$exon->end,$strand);
		my $exon_3prime = downmost($exon->start,$exon->end,$strand);
		if($params{include_introns} && $exon_number>0 && is_downstream($exon_5prime,$that_exon_3prime,$strand)){
            print_slice(
                $gene_name.INTRON.$exon_number,
                $gene_name,
                $gene->stable_id,
                $gene->slice->seq_region_name,
                $that_exon_3prime+$down*1,
                $exon_5prime+$up*1,
                $strand,
                $params{bases_region_of_interest_5_prime_upstream_of_introns},
                $params{bases_region_of_interest_3_prime_downstream_of_introns},
                INTRON,
                $gene->biotype,
                $gene->description);
		}
		$exon_number++;
		$that_exon_3prime = $exon_3prime;
		my @exon_infos;
		my @exon_ranges;
		my $utr5_info;
		my $cds_info;
		my $utr3_info;
        my $this_up = $params{bases_region_of_interest_5_prime_upstream_of_exons};
        my $this_down = $params{bases_region_of_interest_3_prime_downstream_of_exons};
        if($params{separate_records_for_promoter_upstream_downstream_regions}){
            $this_up = 0;
            $this_down = 0;
        }
        elsif($exon_number==1
              &&$params{bases_region_of_interest_gene_upstream_promoter_region}>$params{bases_region_of_interest_5_prime_upstream_of_exons}
              &&$params{include_gene_upstream_promoter_region}){
            $this_up = $params{bases_region_of_interest_gene_upstream_promoter_region};
        }
        if($exon_number==@union_exons
           &&$params{bases_region_of_interest_gene_downstream_region}>$params{bases_region_of_interest_3_prime_downstream_of_exons}
           &&$params{include_gene_downstream_region}){
            $this_down = $params{bases_region_of_interest_gene_downstream_region};
        }
		#gene is non-coding
		unless(@union_coding){
            print_slice(
                $gene_name.EXON.$exon_number,
                $gene_name,
                $gene->stable_id,
                $gene->slice->seq_region_name,
                $exon->start,
                $exon->end,
                $strand,
                $this_up,
                $this_down,
                NONCODING,
                $gene->biotype,
                $gene->description);
		}
		else{
			unless($params{separate_records_for_utr_cds}
				   ||(($params{include_5_prime_utr_exons}&&is_upstream($exon_5prime,$coding_5prime,$strand))
					&&!($params{include_cds_exons} && is_upstream($exon_5prime,$coding_3prime,$strand) && is_downstream($exon_3prime,$coding_5prime,$strand))
					&&($params{include_3_prime_utr_exons} && is_downstream($exon_3prime,$coding_3prime,$strand)))){
				print_slice(
					$gene_name.EXON.$exon_number,
					$gene_name,
					$gene->stable_id,
					$gene->slice->seq_region_name,
					$exon->start,
					$exon->end,
					$strand,
					$this_up,
					$this_down,
					EXON,
					$gene->biotype,
					$gene->description);
			}
			else{
				if($params{include_5_prime_utr_exons} && is_upstream($exon_5prime,$coding_5prime,$strand)){
					my $element_range = get_bio_range($exon_5prime,upmost($exon_3prime,$utr5_3prime,$strand),$strand);
					my $extra_up = $this_up;
					my $extra_down = $exon_3prime==$utr5_3prime?$this_down:0;
					print_slice(
						$gene_name.EXON.$exon_number.UTR5,
						$gene_name,
						$gene->stable_id,
						$gene->slice->seq_region_name,
						$element_range->start,
						$element_range->end,
						$strand,
						$this_up,
						$this_down,
						UTR5,
						$gene->biotype,
						$gene->description);
				}
				if($params{include_cds_exons} && is_upstream($exon_5prime,$coding_3prime,$strand) && is_downstream($exon_3prime,$coding_5prime,$strand)){
					my $element_range = get_bio_range(downmost($coding_5prime,$exon_5prime,$strand),upmost($coding_3prime,$exon_3prime,$strand),$strand);
					my $extra_up = $coding_5prime==$exon_5prime?$this_up:0;
					my $extra_down = $coding_3prime==$exon_3prime?$this_down:0;
					print_slice(
						$gene_name.EXON.$exon_number.CDS,
						$gene_name,
						$gene->stable_id,
						$gene->slice->seq_region_name,
						$element_range->start,
						$element_range->end,
						$strand,
						$this_up,
						$this_down,
						CDS,
						$gene->biotype,
						$gene->description);
				}
				if($params{include_3_prime_utr_exons} && is_downstream($exon_3prime,$coding_3prime,$strand)){
					my $element_range = get_bio_range(downmost($utr3_5prime,$exon_5prime,$strand),$exon_3prime,$strand);
					my $extra_up = $utr3_5prime==$exon_5prime?$this_up:0;
					my $extra_down = $this_down;
					print_slice(
						defined($gene_name)?$gene_name:defined($gene->external_name)?$gene->external_name:$gene->stable_id.EXON.$exon_number.UTR3,
						defined($gene_name)?$gene_name:defined($gene->external_name)?$gene->external_name:$gene->stable_id,
						$gene->stable_id,
						$gene->slice->seq_region_name,
						$element_range->start,
						$element_range->end,
						$strand,
						$this_up,
						$this_down,
						UTR3,
						$gene->biotype,
						$gene->description);
				}
			}
		}
	}
} # !print_gene_regions

sub print_log {
	my ($msg, $type) = @_;
	$msg = "$msg\n";
	if ($type){
#		$type = lc($type);
		if ($type=~/ERROR/i){
			$counts{'! error'}++;
			$msg = "$type: $msg";
		}
		if ($type=~/WARNING/i){
			$counts{'? warning'}++;
			$msg = "$type: $msg";
		}
		print STDERR $msg unless $opts{Q} || $type=~/DETAIL/i;
		print LOGFILE $msg unless $opts{D} || $type=~/UPDATE/i;
		print SUMMARYFILE $msg unless $opts{L} || $type=~/DETAIL/i || $type=~/UPDATE/i;
	}
} # !print_log

sub print_slice {
	my (
		$id,
		$gene_external_name,
		$gene_stable_id,
		$chr,
		$roi_start,
		$roi_end,
		$strand,
		$roi_up, #extra bases upstream of roi
		$roi_down, #extra bases downstream of roi
		$type,
		$biotype,
		$description) = @_;
	#add extra bases to roi
	if(defined($params{include_only_these_biotypes})){
		my @biotypes = split(/[\W]+/,$params{include_only_these_biotypes});
		if(@biotypes && $biotypes[0] && defined($biotype) && !contains(\@biotypes,$biotype)) {
			return;
		}
	}
	if($strand==1){
		$roi_start -= $roi_up;
		$roi_end += $roi_down;
	}
	else{
		$roi_start -= $roi_down;
		$roi_end += $roi_up;
	}
	if(bases($roi_start,$roi_end) < $params{exclude_region_of_interest_less_than_bases}){
		return;
	}
	$roi_up += $params{padding_bases_5_prime_upstream_of_each_region_of_interest};
	$roi_down += $params{padding_bases_3_prime_downstream_of_each_region_of_interest};
	#pad ends of original roi
	if($strand==1){
		$roi_start -= $params{padding_bases_5_prime_upstream_of_each_region_of_interest};
		$roi_end += $params{padding_bases_3_prime_downstream_of_each_region_of_interest};
	}
	else{
		$roi_end += $params{padding_bases_5_prime_upstream_of_each_region_of_interest};
		$roi_start -= $params{padding_bases_3_prime_downstream_of_each_region_of_interest};
	}
	print_sub_slice(
		$id,
		$gene_external_name,
		$gene_stable_id,
		$chr,
		$roi_start,
		$roi_end,
		$strand,
		$roi_up,
		$roi_down,
		$type,
		defined($biotype) ? $biotype : "",
		defined($description) ? $description : "")
} # !print_slice

sub print_sub_slice {
	my (
		$id,
		$gene_external_name,
		$gene_stable_id,
		$chr,
		$roi_start,
		$roi_end,
		$strand,
		$roi_up,
		$roi_down,
		$type,
		$biotype,
		$description) = @_;
	my $length = bases($roi_start,$roi_end);
	#if region of interest is too long, recursively call print_slice with smaller slices
	if ($length > $params{maximum_region_of_interest_per_record}) {
		my $count = ceil($length / $params{maximum_region_of_interest_per_record});
		if ($strand == 1) {
			my $start = $roi_start;
			my $end = $roi_start + $params{maximum_region_of_interest_per_record} - 1;
			for (my $i = 1; $i <= $count; $i++) {
				print_sub_slice(
					$id . "_" . $i . "/" . $count,
					$gene_external_name,
					$gene_stable_id,
					$chr,
					$start,
					$end,
					$strand,
					$i == 1 ? $roi_up : 0,
					$i == $count ? $roi_down : 0,
					$type,
					$biotype,
					$description);
				$start += $params{maximum_region_of_interest_per_record};
				$end = min($end += $params{maximum_region_of_interest_per_record},$roi_end),
			}
			return;
		}
		else {
			my $start = $roi_end;
			my $end = $roi_end - $params{maximum_region_of_interest_per_record} + 1;
			for (my $i = 1 ; $i <= $count ; $i++) {
				print_sub_slice(
					$id . "_" . $i . "/" . $count,
					$gene_external_name,
					$gene_stable_id,
					$chr,
					$start,
					$end,
					$strand,
					$i == 1 ? 0 : $roi_up,
					$i == $count ? 0 : $roi_down,
					$type,
					$biotype,
					$description);
				$start -= $params{maximum_region_of_interest_per_record};
				$end = max($end -= $params{maximum_region_of_interest_per_record},$roi_start),
			}
			return;
		}
	}
	my $extra = $params{minimum_region_of_interest_per_record} - $length;
	if ($extra > 0) {
		$length += $extra;
		if ($params{add_bases_upstream_to_region_of_interest_to_get_minimum_length}&&!$params{add_bases_downstream_to_region_of_interest_to_get_minimum_length}) {
			if ($strand == 1) {
				$roi_start -= $extra;
				$roi_down += $extra;
			}
			else {
				$roi_end += $extra;
				$roi_up += $extra;
			}
		}
		elsif ($params{add_bases_downstream_to_region_of_interest_to_get_minimum_length}&&!$params{add_bases_upstream_to_region_of_interest_to_get_minimum_length}) {
			if ($strand == 1) {
				$roi_end += $extra;
				$roi_up += $extra;
			}
			else {
				$roi_start -= $extra;
				$roi_down += $extra;
			}
		}
		else { #center
			my $extra1 = int($extra / 2);
			my $extra2 = $extra - $extra1;
			if ($strand == 1) {
				$roi_start -= $extra1;
				$roi_down += $extra1;
				$roi_end += $extra2;
				$roi_up += $extra2;
			}
			else {
				$roi_end += $extra1;
				$roi_up += $extra1;
				$roi_start -= $extra2;
				$roi_down += $extra2;
			}
		}
	}
	my $slice_start = $roi_start - $strand * $params{padding_bases_5_prime_upstream_of_each_record};
	my $slice_end = $roi_end + $strand * $params{padding_bases_3_prime_downstream_of_each_record};
	#genomic interval records for Agilent eArray probes
	unless($opts{G}){
		my $id = sprintf "%s%s%s",
			$id,
			$roi_up==0?'':($strand==1?'-':'+').$roi_up,
			$roi_down==0?'':($strand==1?'+':'-').$roi_down;
		my $record_description = sprintf("%s|chr%s:%u-%u|%s|%ub|%s|%s|roi:%u-%u|%ub|pad-%u+%u|%s|%s|%s",
			$id,
			$chr,
			$slice_start,
			$slice_end,
			$strand==1?"+":"-",
			bases($slice_start,$slice_end),
			$gene_external_name,
			$gene_stable_id,
			$roi_start,
			$roi_end,
			bases($roi_start,$roi_end),
			$strand == 1
			? $params{padding_bases_5_prime_upstream_of_each_record}
			: $params{padding_bases_3_prime_downstream_of_each_record},
			$strand == 1
			? $params{padding_bases_3_prime_downstream_of_each_record}
			: $params{padding_bases_5_prime_upstream_of_each_record},
			$type,
			$biotype,
			$description);
		push @all_regions, {num_chr=>chr2numchr($chr), start=>$slice_start, end=>$slice_end, strand=>$strand};
		if ($opts{S}){
			print_log($record_description,DETAIL);
		}
		$counts{'output genomic region record'}++;
	}
	#sequence records
	unless($opts{S}){
		if(!defined($params{sequence_file_format}) || !$params{sequence_file_format}){
			$params{sequence_file_format} = "fasta";
		}
		my $seq;
		if ($params{include_template_strand}) {
			my $slice = $slice_adaptor->fetch_by_region('chromosome', $chr, $slice_start, $slice_end, $strand);
			if ($slice) {
				my $sequence;
				if ($params{show_variations}) {
					$sequence = fetch_sequence_with_ambiguity($chr, $slice_start, $slice_end, $strand);
				}
				elsif($params{mask_repeats}) {
					my $repeat_masked_slice = $slice->get_repeatmasked_seq();
					$repeat_masked_slice->soft_mask($params{hard_mask_repeats}?0:1);
					$sequence = $repeat_masked_slice->seq();
				}
				else {
					$sequence = $slice->seq;
				}
				my $id = sprintf "%s%s%s",
					$id,
                    $roi_up==0?'':($strand==1?'-':'+').$roi_up,
                    $roi_down==0?'':($strand==1?'+':'-').$roi_down;
				my $record_description = sprintf("|%s|%ub|%s|%s|roi:%u-%u|%ub|pad-%u+%u|%s|%s|%s",
					slicef($slice),
					bases($slice_start,$slice_end),
					$gene_external_name,
					$gene_stable_id,
					$roi_start,
					$roi_end,
					bases($roi_start,$roi_end),
					$strand == 1
					? $params{padding_bases_5_prime_upstream_of_each_record}
					: $params{padding_bases_3_prime_downstream_of_each_record},
					$strand == 1
					? $params{padding_bases_3_prime_downstream_of_each_record}
					: $params{padding_bases_5_prime_upstream_of_each_record},
					$type,
					$biotype,
					$description);
				if($params{single_fasta_identifier}){
					$seq = Bio::Seq->new(
						-seq => $sequence,
						-display_id => $id,
						-alphabet => 'dna');
				}
				else{
					$seq = Bio::Seq->new(
						-seq => $sequence,
						-description => $record_description,
						-display_id => $id,
						-alphabet => 'dna');
				}
				my $outseq = Bio::SeqIO->new(
					-fh => \*STDOUT,
					-format => $params{sequence_file_format});
				$outseq->write_seq($seq);
				print_log(
					sprintf("%s %s\n%s...%s",
					$seq->display_id, $seq->description,
					substr($seq->seq, 0, 30),
					substr($seq->seq, -30)),DETAIL);
				$counts{'output sequence record'}++;
			}
			else {
				print_log(sprintf("Can't get %s", $id), ERROR);
			}
		}
		if ($params{include_rna_like_strand}) {
			my $slice = $slice_adaptor->fetch_by_region('chromosome', $chr, $slice_start, $slice_end, -$strand);
			if ($slice) {
				my $sequence;
				if ($params{show_variations}) {
					$sequence = fetch_sequence_with_ambiguity($chr, $slice_start, $slice_end, $strand);
				}
				elsif($params{mask_repeats}) {
					my $repeat_masked_slice = $slice->get_repeatmasked_seq();
					$repeat_masked_slice->soft_mask($params{hard_mask_repeats}?0:1);
					$sequence = $repeat_masked_slice->seq();
				}
				else {
					$sequence = $slice->seq;
				}
				my $id = sprintf "%s%s%s",
					$id,
			$roi_up==0?'':($strand==1?'-':'+').$roi_up,
			$roi_down==0?'':($strand==1?'+':'-').$roi_down;
				my $record_description = sprintf("|%s|%ub|%s|%s|roi:%u-%u|%ub|pad-%u+%u|%s",
					slicef($slice),
					bases($slice_start,$slice_end),
					$gene_external_name,
					$gene_stable_id,
					$roi_start,
					$roi_end,
					bases($roi_start,$roi_end),
					$strand == 1
					? $params{padding_bases_5_prime_upstream_of_each_record}
					: $params{padding_bases_3_prime_downstream_of_each_record},
					$strand == 1
					? $params{padding_bases_3_prime_downstream_of_each_record}
					: $params{padding_bases_5_prime_upstream_of_each_record},
					$type,
					$biotype,
					$description);
				if($params{single_fasta_identifier}){
					$seq = Bio::Seq->new(
						-seq => $sequence,
						-display_id => $id . "_rc",
						-alphabet => 'dna');
				}
				else{
					$seq = Bio::Seq->new(
						-seq => $sequence,
						-description => $record_description,
						-display_id => $id . "_rc",
						-alphabet => 'dna');
				}
				my $outseq = Bio::SeqIO->new(
					-fh => \*STDOUT,
					-format => $params{format});
				$outseq->write_seq($seq);
				print_log(
					sprintf("%s %s\n%s...%s",
						$seq->display_id, $seq->description,
						substr($seq->seq, 0, 30), substr($seq->seq, -30)),DETAIL);
				$counts{'output sequence record'}++;
			}
			else {
				print_log(sprintf("Can't get %s_rc", $id), ERROR);
			}
		}
	}
	if ($length > $max_roi) {
		$max_roi = $length;
	}
	if ($length < $min_roi) {
		$min_roi = $length;
	}
	$genes{$gene_stable_id}[0] += bases($roi_start,$roi_end);
	$genes{$gene_stable_id}[1] += bases($slice_start,$slice_end);
	$counts{'ROI base'} += bases($roi_start,$roi_end);
	$counts{'ROI + padding base'} += bases($slice_start,$slice_end);
	$all_lengths{$length}++;
#	push @roi_lengths, $length;
} # !print_sub_slice

sub print_usage {
	print STDERR USAGE . "Default settings:\n";
	print STDERR map { "\t\$$_ = $params{$_}\n" } sort keys %params;
	print STDERR HEADER_FORMAT;
}

sub chr2numchr{
    my $chr = uc(shift);
	unless ($chr =~ /[XYMT0-9]{1,2}/i){
		print_log(sprintf("%s isn't a chromosome", $chr),ERROR);
	}
    return $chr eq "X" ? 23 : $chr eq "Y" ? 24 : ($chr eq "M" || $chr eq "MT") ? 25 : $chr;
}

sub numchr2chr{
    my $chr = shift;
    return $chr == 23 ? "X" : $chr == 24 ? "Y" : $chr == 25 ? "M" : $chr;
}

sub slicef {
	my $slice = shift;
	return sprintf($params{slice_format},$slice->seq_region_name,$slice->start,$slice->end,$slice->strand==1?"+":"-");
} # !slicef

sub ss {
	return $_[0] == 1 ? "" : "s";
} # !ss

########################################################
sub configure {
	my $args = shift;
	
	my $config = {};
	
	GetOptions(
		$config,
		'help',
		
		# input options,
		'config=s',
		'input_file=s',
		'format=s',
		
		# DB options
		'species=s',
		'registry=s',
		'host=s',
		'user=s',
		'port=s',
		'password=s',
		'db_version=i',
		'genomes',
		
		# runtime options
		'most_severe',
		'buffer_size=i',
		'chunk_size=s',
		'check_ref',
		'check_existing=i',
		'failed=i',
		'whole_genome',
		'tmp_dir=s',
		'gp',
		
		# output options
		'output_file=s',
		'terms=s',
		'verbose',
		'quiet',
		'coding_only',
		'protein',
		'hgnc',
		'hgvs',
		'sift=s',
		'polyphen=s',
		'condel=s',
	);
	
	# print usage message if requested or no args supplied
	if(defined($config->{help}) || !$args) {
		&usage;
		exit(0);
	}
	
	# config file?
	if(defined $config->{config}) {
		
		open CONFIG, $config->{config} or die "ERROR: Could not open config file \"".$config->{config}."\"\n";
		
		while(<CONFIG>) {
			next if /^\#/;
			my ($key, $value) = split /\s+|\=/;
			$key =~ s/^\-//g;
			$config->{$key} = $value unless defined $config->{$key};
		}
		
		close CONFIG;
	}
	
	# check file format
	if(defined $config->{input_format}) {
		die "ERROR: Unrecognised input format specified \"".$config->{input_format}."\"\n" unless $config->{input_format} =~ /pileup|vcf|guess/i;
	}
	
	# output term
	if(defined $config->{terms}) {
		die "ERROR: Unrecognised consequence term type specified \"".$config->{terms}."\" - must be one of ensembl, so, ncbi\n" unless $config->{terms} =~ /ensembl|display|so|ncbi/i;
		if($config->{terms} =~ /ensembl|display/i) {
			$config->{terms} = 'display';
		}
		else {
			$config->{terms} = uc($config->{terms});
		}
	}
	
	# check nsSNP tools
	foreach my $tool(grep {defined $config->{lc($_)}} qw(SIFT PolyPhen Condel)) {
		die "ERROR: Unrecognised option for $tool \"", $config->{lc($tool)}, "\" - must be one of p (prediction), s (score) or b (both)\n" unless $config->{lc($tool)} =~ /^(s|p|b)/;
	}
	
	# summarise options if verbose
	if(defined $config->{verbose}) {
		my $header =<<INTRO;
#----------------------------------#
# ENSEMBL VARIANT EFFECT PREDICTOR #
#----------------------------------#

version 2.0

By Will McLaren (wm2\@ebi.ac.uk)

Configuration options:

INTRO
		print $header;
		
		my $max_length = (sort {$a <=> $b} map {length($_)} keys %$config)[-1];
		
		foreach my $key(sort keys %$config) {
			print $key.(' ' x (($max_length - length($key)) + 4)).$config->{$key}."\n";
		}
		
		print "\n".("-" x 20)."\n\n";
	}
	
	# connection settings for Ensembl Genomes
	if($config->{genomes}) {
		$config->{host} ||= 'mysql.ebi.ac.uk';
		$config->{port} ||= 4157;
	}
	
	# connection settings for main Ensembl
	else {
		$config->{species} ||= "homo_sapiens";
		$config->{host}    ||= 'ensembldb.ensembl.org';
		$config->{port}    ||= 5306;
	}
	
	# set defaults
	$config->{user}         ||= 'anonymous';
	$config->{buffer_size}  ||= 5000;
	$config->{chunk_size}   ||= '50kb';
	$config->{output_file}  ||= "variant_effect_output.txt";
	$config->{tmpdir}       ||= '/tmp';
	$config->{format}       ||= 'guess';
	$config->{terms}        ||= 'display';
	
	$config->{include_failed} = 1 unless defined $config->{include_failed};
	$config->{check_existing} = 1 unless (defined $config->{check_existing} || defined $config->{whole_genome});
	$config->{chunk_size} =~ s/mb?/000000/i;
	$config->{chunk_size} =~ s/kb?/000/i;
	
	# connect to databases
	$config->{reg} = &connect_to_dbs($config);
	
	# get input file handle
	$config->{in_file_handle} = &get_in_file_handle($config);
	
	# configure output file
	$config->{out_file_handle} = &get_out_file_handle($config);
	
	return $config;
}

sub usage {
	my $usage =<<END;
#----------------------------------#
# GENOMIC INTERVAL LISTER #
#----------------------------------#

version 1.0.beta.1

By Michael Yourshaw (myourshaw\@ucla.edu)

Usage:
perl genomic_interval_lister.pl [arguments]

Options
--help                 Display this message and quit
--verbose              Display verbose output as the script runs [default: off]
--quiet                Suppress status and warning messages [default: off]

--config               Load configuration from file. Any command line options
                       specified overwrite those in the file [default: off]

-i | --input_file      Input file - if not specified, reads from STDIN. Files
                       may be gzip compressed.
--format               Alternative input file format - one of "pileup", "vcf"
-o | --output_file     Output file [default: "genomic_interval_list.txt"]

-t | --terms           Type of consequence terms to output - one of "ensembl", "SO",
                       "NCBI" [default: ensembl]
					   
--sift=[p|s|b]         Add SIFT [p]rediction, [s]core or [b]oth [default: off]
--polyphen=[p|s|b]     Add PolyPhen [p]rediction, [s]core or [b]oth [default: off]
--condel=[p|s|b]       Add Condel SIFT/PolyPhen consensus [p]rediction, [s]core or
                       [b]oth [default: off]

NB: SIFT, PolyPhen and Condel predictions are currently available for human only

--hgnc                 If specified, HGNC gene identifiers are output alongside the
                       Ensembl Gene identifier [default: off]
--hgvs                 Output HGVS identifiers (coding and protein) [default: off]
--protein              Output Ensembl protein identifer [default: off]

--coding_only          Only return consequences that fall in the coding region of
                       transcripts [default: off]

--check_ref            If specified, checks supplied reference allele against stored
                       entry in Ensembl Core database [default: off]
--check_existing=[0|1] If set to 1, checks for existing co-located variations in the
                       Ensembl Variation database [default: 1]
--failed=[0|1]         If set to 1, includes variations flagged as failed when checking
                       for co-located variations. Only applies in Ensembl 61 or later
                       [default: 1]
--gp                   If specified, tries to read GRCh37 position from GP field in the
                       INFO column of a VCF file. Only applies when VCF is the input
                       format and human is the species [default: off]

-s | --species         Species to use [default: "human"]
--host                 Manually define database host [default: "ensembldb.ensembl.org"]
-u | --user            Database username [default: "anonymous"]
--port                 Database port [default: 5306]
--password             Database password [default: no password]
--genomes              Sets DB connection params for Ensembl Genomes [default: off]
-r | --registry_file   Registry file to use defines DB connections [default: off]
                       Defining a registry file overrides above connection settings.
--db_version=[number]  Force script to load DBs from a specific Ensembl version. Not
                       advised due to likely incompatibilities between API and DB
					   
-w | --whole_genome    Run in whole genome mode [default: off]
                       Recommended for use with data covering a whole genome,
                       chromosome or gene/set of genes e.g. from resequencing. For
                       better performance, set --buffer_size higher (>10000 if memory
                       allows).
-b | --buffer_size     Sets the number of variants sent in each batch [default: 5000]
                       Increasing buffer size can retrieve results more quickly
                       but requires more memory. Only applies to whole genome mode.
--chunk_size           Sets the chunk size of internal data structure [default: 50kb]
                       Setting this lower may improve speed for variant-dense
                       datasets. Only applies to whole genome mode.
					   
NB: Whole genome mode disables --check_existing option and gene column by default.
END

	print $usage;
}


sub connect_to_dbs {
	my $config = shift;
	
	# get registry
	my $reg = 'Bio::EnsEMBL::Registry';
	
	# load DB options from registry file if given
	if(defined($config->{registry})) {
		debug("Loading DB config from registry file ", $config->{registry});
		$reg->load_all($config->{registry});
	}
	
	# otherwise manually connect to DB server
	else {
		$reg->load_registry_from_db(
			-host       => $config->{host},
			-user       => $config->{user},
			-pass       => $config->{password},
			-port       => $config->{port},
			-db_version => $config->{db_version},
			-species    => $config->{species} =~ /^[a-z]+\_[a-z]+/i ? $config->{species} : undef,
			-verbose    => $config->{verbose},
		);
	}
	
	$reg->set_disconnect_when_inactive();
	
	if($config->{verbose}) {
		# get a meta container adaptors to check version
		my $core_mca = $reg->get_adaptor($config->{species}, 'core', 'metacontainer');
		my $var_mca = $reg->get_adaptor($config->{species}, 'variation', 'metacontainer');
		
		if($core_mca && $var_mca) {
			debug(
				"Connected to core version ", $core_mca->get_schema_version, " database ",
				"and variation version ", $var_mca->get_schema_version, " database"
			);
		}
	}
	
	return $reg;
}


sub get_in_file_handle {
	my $config = shift;

	# define the filehandle to read input from
	my $in_file_handle = new FileHandle;
	
	if(defined($config->{input_file})) {
		
		# check defined input file exists
		die("ERROR: Could not find input file ", $config->{input_file}, "\n") unless -e $config->{input_file};
		
		if($config->{input_file} =~ /\.gz$/){
			$in_file_handle->open("zcat ". $config->{input_file} . " | " ) or die("ERROR: Could not read from input file ", $config->{input_file}, "\n");
		}
		else {
			$in_file_handle->open( $config->{input_file} ) or die("ERROR: Could not read from input file ", $config->{in_file}, "\n");
		}
	}
	
	# no file specified - try to read data off command line
	else {
		$in_file_handle = 'STDIN';
		debug("Reading input from STDIN (or maybe you forgot to specify an input file?)...") unless defined $config->{quiet};
	}
	
	return $in_file_handle;
}


sub get_out_file_handle {
	my $config = shift;
	
	# define filehandle to write to
	my $out_file_handle = new FileHandle;
	$out_file_handle->open(">".$config->{output_file}) or die("ERROR: Could not write to output file ", $config->{output_file}, "\n");
	
	# make header
	my $time = &get_time;
	my $core_mca = $config->{reg}->get_adaptor($config->{species}, 'core', 'metacontainer');
	my $db_string = $core_mca->dbc->dbname." on ".$core_mca->dbc->host if defined $core_mca;
	my $version_string =
		"Using API version ".$config->{reg}->software_version.
		", DB version ".(defined $core_mca && $core_mca->get_schema_version ? $core_mca->get_schema_version : '?');
	
	my $header =<<HEAD;
## Genomic Interval Lister v1.0.beta.1
## Output produced at $time
## Connected to $db_string
## $version_string
## Extra column keys:
## HGNC     : HGNC gene identifier
## ENSP     : Ensembl protein identifer
## HGVSc    : HGVS coding sequence name
## HGVSp    : HGVS protein sequence name
## SIFT     : SIFT prediction
## PolyPhen : PolyPhen prediction
## Condel   : Condel SIFT/PolyPhen consensus prediction
HEAD
	
	# add headers
	print $out_file_handle $header;
	
	# add column headers
	print $out_file_handle join "\t", qw(
		#Uploaded_variation
		Location
		Allele
		Gene
		Transcript
		Consequence
		cDNA_position
		CDS_position
		Protein_position
		Amino_acids
		Codons
		Existing_variation
		Extra
	);
	
	print $out_file_handle "\n";
	
	return $out_file_handle;
}


# gets time
sub get_time() {
	my @time = localtime(time());

	# increment the month (Jan = 0)
	$time[4]++;

	# add leading zeroes as required
	for my $i(0..4) {
		$time[$i] = "0".$time[$i] if $time[$i] < 10;
	}

	# put the components together in a string
	my $time =
 		($time[5] + 1900)."-".
 		$time[4]."-".
 		$time[3]." ".
		$time[2].":".
		$time[1].":".
		$time[0];

	return $time;
}

# prints debug output with time
sub debug {
	my $text = (@_ ? (join "", @_) : "No message");
	my $time = get_time;
	
	print $time." - ".$text.($text =~ /\n$/ ? "" : "\n");
}

sub format_coords {
	my ($start, $end) = @_;
	
	if(!defined($start)) {
		return '-';
	}
	elsif(!defined($end)) {
		return $start;
	}
	elsif($start == $end) {
		return $start;
	}
	elsif($start > $end) {
		return $end.'-'.$start;
	}
	else {
		return $start.'-'.$end;
	}
}

sub find_existing {
	my $new_vf = shift;
	
	if(defined($new_vf->adaptor->db)) {
		
		my $sth = $new_vf->adaptor->db->dbc->prepare(qq{
		  SELECT variation_name, source_id
		  FROM variation_feature
		  WHERE seq_region_id = ?
		  AND seq_region_start = ?
		  AND seq_region_end = ?
		});
		
		$sth->execute($new_vf->slice->get_seq_region_id, $new_vf->start, $new_vf->end);
		
		my ($name, $source);
		$sth->bind_columns(\$name, \$source);
		
		my %by_source;
		
		push @{$by_source{$source}}, $name while $sth->fetch;
		$sth->finish;
		
		if(scalar keys %by_source) {
			foreach my $s(sort {$a <=> $b} keys %by_source) {
				return shift @{$by_source{$s}};
			}
		}
	}
	
	return undef;
}
