#!/usr/bin/perl -w
# ©2007,2008 Michael Yourshaw All Rights Reserved

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
use strict;
use warnings;
use Getopt::Std;
use File::Spec;
use Date::Calc qw(Today_and_Now Delta_YMDHMS);
use POSIX qw(ceil floor);
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
#BioPerl
use Bio::Seq;
use Bio::SeqIO;
#EnsEMBL Perl API
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::Utils::Sequence qw (sequence_with_ambiguity);

use constant BIG_NUMBER => 9999999999;
#the following must have unique values
use constant CDS => "cds";
use constant DOWNSTREAM => "down";
use constant ERROR => "error";
use constant EXON => "e";
use constant GENE => "gene";
use constant INTERGENETIC => "intergenetic";
use constant INTRON => "i";
use constant NONCODING => "noncoding";
use constant REGION => "region";
use constant TRANSCRIPT => "tr";
use constant UPSTREAM => "up";
use constant WARNING => "warning";
use constant UTR3 => "3utr";
use constant UTR5 => "5utr";


my @start_time = Today_and_Now();

#parameter keys must be lower case
my %params = (
	add_bases_downstream_to_region_of_interest_to_get_minimum_length => 0, #bool
	add_bases_upstream_to_region_of_interest_to_get_minimum_length => 0, #bool
	bases_region_of_interest_3_prime_downstream_of_exons => 0,
	bases_region_of_interest_3_prime_downstream_of_introns => 0,
	bases_region_of_interest_5_prime_upstream_of_exons => 0,
	bases_region_of_interest_5_prime_upstream_of_introns => 0,
	bases_region_of_interest_gene_downstream_region => 0,
	bases_region_of_interest_gene_upstream_promoter_region => 0,
	biotypes => "",
	create_log_file => 1, #bool
	create_padded_genomic_intervals_file => 0, #bool (default = 0)
	create_regions_of_interest_genomic_intervals_file => 1, #bool (default = 1)
	create_sequence_file => 0, #bool
	exclude_region_of_interest_less_than_bases => 0,
	hard_mask_repeats => 0, #bool
	include_3_prime_utr_exons => 1, #bool
	include_5_prime_utr_exons => 1, #bool
	include_all_elements_one_record_per_gene => 0, #bool overrides all includes & separates
	include_cds_exons => 1, #bool
	include_gene_downstream_region => 0,
	include_gene_upstream_promoter_region => 0,
	include_intergenetic_sequence_in_genomic_regions => 0, #bool
	include_introns => 0, #bool
	include_rna_like_strand => 0, #bool
	include_template_strand => 1, #bool
	input_format => "", # '' (default) | omim_result gene_result
	mask_repeats => 0, #bool
	maximum_record_length => 1000000000000000, #split records if necessary
	maximum_region_of_interest_per_record => 1000000000000000,
	minimum_record_length => 0, #pad records as necessary
	minimum_region_of_interest_per_record => 0,
	merge_transcripts => 1, #bool default 1
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
my $that_id;
my @all_coding;
my @all_exons;
my @all_transcripts;
my @all_regions;

use constant USAGE => "engene usage:\n"
	 . "perl engene.pl prints usage.\n"
	 . "perl engene.pl [-S] [-g genomic_intervals_file | -G] [-l log_file | -L] input_file_list [> sequence_file]\n"
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
	 . "-G = no genomic intervals output\n"
	 . "-l log_file = specify log file path (default is input_file1.log and STDERR)\n"
	 . "-L = log to STDERR only\n"
	 . "WARNING: sequence, genomic intervals and log files are created by default\n"
	 . "and overwrite existing files of the same name.\n"
;
use constant HEADER_FORMAT => "header format:\n"
	. "fragment id|coordinates[strand]|length|parent gene|Ensembl stable id|region of interest coordinates|region of interest length|padding-upstream+downstream|fragment type(s)|gene biotype|gene description\n"
	. "example:\n"
	. ">SEZ6_e1-35+20|chr17:24356939-24357324[-]|386b|SEZ6|ENSG00000063015|roi:24356939-24357324|386b|pad-0+0|CDS,UTR|protein_coding|seizure related 6 homolog isoform 1 [Source:RefSeq_peptide;Acc:NP_849191]\n"
	. "fragment id = gene name or region coordinates[_up|_e#[_5utr|_3utr|_utr]|_i#|_down ][-additional upstream roi][+additional downstream roi][_fragment part/total fragment parts][_rc]\n"
	. "up=upstream, e=exon, utr=untranslated region, i=intron, down=downstream, rc=reverse complement\n"
	. "fragment types = REGION|INTERGENETIC|GENE|UPSTREAM|CDS[,UTR]|UTR|INTRON|DOWNSTREAM|NONCODING\n"
	. "biotype = C_segment|D_segment|J_segment|miRNA|miRNA_pseudogene|misc_RNA|misc_RNA_pseudogene|Mt_rRNA|Mt_tRNA_pseudogene|protein_coding|pseudogene|repeat|retrotransposed|rRNA|rRNA_pseudogene|scRNA|scRNA_pseudogene|snoRNA|snoRNA_pseudogene|snRNA|snRNA_pseudogene|tRNA_pseudogene|V_segment\n"
;

#options
my %opts; #sequence, genomic intervals, log
getopts("s:Sg:Gl:L", \%opts);
#files
unless (@ARGV) {
	print_usage();
	exit;
}
my $input_file = "";
my $input_file1 = $ARGV[0];
my ($vol,$dir,$file);
if ( $input_file1 eq "-" ) {
	print STDERR "Input genes/genomic regions or program settings, one per line. Ctrl-D when done\n";
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
unless($opts{S}){
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
	else{
		$sequence_file = File::Spec->catpath($vol,$dir,"$file.seq");
	}
	while($ARGV[$#ARGV] eq ">" || $ARGV[$#ARGV] eq $gt){
		pop(@ARGV);
	}
	open STDOUT, ">", $sequence_file or die "can't open $sequence_file: $!";
}

my $genomic_intervals_file;
$genomic_intervals_file = $opts{g} or $genomic_intervals_file = File::Spec->catpath($vol,$dir, $clean_file."_intervals.txt");
open GENOMICINTERVALSFILE, ">", $genomic_intervals_file or die "Can't open $genomic_intervals_file: $!" unless $opts{G};
my $log_file;
$log_file = $opts{l} or $log_file = File::Spec->catpath($vol,$dir,"$file.log");
open LOGFILE, ">", $log_file or die "Can't open $log_file: $!" unless $opts{L};
#start
print_log(sprintf "START: %s", datetimef(@start_time));
print_log(HEADER_FORMAT);
print_log("Default settings:");
print LOGFILE map { "\t\$$_ = $params{$_}\n" } sort keys %params unless $opts{L};
print STDERR map { "\t\$$_ = $params{$_}\n" } sort keys %params;
print_log(sprintf "sequence file: %s",$sequence_file ? $sequence_file : "STDOUT") unless $opts{S};
print_log(sprintf "genomic intervals file: %s",$genomic_intervals_file ? $genomic_intervals_file : "STDOUT") unless $opts{G};
print_log(sprintf "log file: %s",$log_file ? $log_file : "STDOUT") unless $opts{L};
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
foreach my $db_adaptor (@db_adaptors) {
    my $db_connection = $db_adaptor->dbc();

    printf STDERR (
        "species/group\t%s/%s\ndatabase\t%s\nhost:port\t%s:%s\n\n",
        $db_adaptor->species(),   $db_adaptor->group(),
        $db_connection->dbname(), $db_connection->host(),
        $db_connection->port()
    );
}
#read from a list of filenames on the command line
#read from SDTIN if no filenames were given
#"-" indicates STDIN
#"someprogram |" indicates the output of another program (7.14)
#read and parse input file
while (<>) {
	unless ($params{stop}){
		$counts{'input line'}++;
		if ($ARGV ne $input_file) {
			print_log(sprintf("input file: %s", File::Spec->rel2abs($ARGV)));
			$input_file = $ARGV;
			$counts{'input file'}++;
		}
		chomp;
		if($_){
			print_log($_);
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
				printf "%s not processed", $_;
				$counts{'? warning'}++;
			}
		}
	}
}
#all input records done

##genomic interval output
#print_all_regions() unless $opts{G};
#
##print statistics
#print_log("lengths of regions of interest (count):");
#print LOGFILE map { "$_ ($all_lengths{$_})," } sort {$a <=> $b} keys %all_lengths unless $opts{L};
#print LOGFILE "\n" unless $opts{L};
#print STDERR map { "$_ ($all_lengths{$_})," } sort {$a <=> $b} keys %all_lengths;
#print STDERR "\n";
#print_log("gene\tstableID\tbasesROI\tbasesROI+padded");
#print LOGFILE map {"$genes{$_}[2]\t$_\t$genes{$_}[0]\t$genes{$_}[1]\n"} sort {$genes{$a}[2] cmp $genes{$b}[2]} keys %genes unless $opts{L};
#print STDERR map {"$genes{$_}[2]\t$_\t$genes{$_}[0]\t$genes{$_}[1]\n"} sort {$genes{$a}[2] cmp $genes{$b}[2]} keys %genes;
#print_log("gene biotypes:");
#print LOGFILE map { " $all_biotypes{$_} $_\n" } sort keys %all_biotypes unless $opts{L};
#print STDERR map { " $all_biotypes{$_} $_\n" } sort keys %all_biotypes;
#$counts{'ROI minimum base'} = $min_roi;
#$counts{'ROI maximum base'} = $max_roi;
#print_log("statistics:");
#print LOGFILE map { sprintf "%10d %s%s\n",$counts{$_},$_,ss($counts{$_}) } sort keys %counts unless $opts{L};
#print STDERR map { sprintf "%10d %s%s\n",$counts{$_},$_,ss($counts{$_}) } sort keys %counts;
my @end_time = Today_and_Now();
print_log(sprintf "job duration %s",durationf(Delta_YMDHMS(@start_time,@end_time)));
print_log(sprintf "DONE: %s",datetimef(@end_time));
exit;

sub get_by_exon_id {
	$id = shift;
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
	my $gene = $gene_adaptor->fetch_by_stable_id($id);
	if ($gene) {
		gene($gene);
	}
	else {
		print_log(sprintf("%s isn't a gene", $id), ERROR);
	}
} # !get_by_gene_id

sub get_by_gene_name {
	$id = shift;
	my @genes;
	if($params{use_ensembl_human}){
		@genes = @{ $gene_adaptor->fetch_all_by_external_name($id) };
	}
	unless (@genes){
		if($params{use_ensembl_human}){
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
			print_log(sprintf("%s identifies more than one gene (%u genes); using %s (%s)",
					$id,
					$count_temp,
					$external_name,
					$stable_id),
				WARNING);
			foreach my $g (@genes) {
				print_log(sprintf("\t%s %s %s (%s)",
					$g->external_name,
					slicef($gene->slice),
					$g->description,
					$g->stable_id));
			}
		}
		if ($gene->external_name ne $id) {
			print_log(
				sprintf(
					"%s will be named %s %s %s (%s)",
					$id,
					$gene->external_name,
					slicef($gene->slice),
					$gene->description,
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
} # !get_by_gene_name

#genomic region (chromosome start end [strand])
sub get_by_region {
	my ($chr, $start, $end, $strand) = @_;
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

sub strand_symbol{
	return shift==1?'+':'-';
}
sub gene {
	my $gene = shift;
    my $gene_name = defined($gene->external_name)?$gene->external_name:$gene->stable_id;
    $genes{$gene->stable_id} = [0,0,$gene_name,0,0];
    $counts{'output gene'}++;
    $all_biotypes{$gene->biotype}++;
	my @transcripts = @{ $gene->get_all_Transcripts };
	my @exons;
	printf STDERR "%s\t%s\t%s\t$params{slice_format}\t%u bases\t%u transcripts\n",$gene_name,$gene->stable_id,$gene->biotype,$gene->slice->seq_region_name,$gene->start,$gene->end,strand_symbol($gene->strand),bases($gene->start,$gene->end),scalar(@transcripts);
	foreach my $transcript(@transcripts){
		@exons = @{$transcript->get_all_Exons};
		printf STDERR "\t%s\t%u bases\t%u exons\n",$transcript->stable_id,$transcript->length,scalar(@exons);
		foreach my $exon(@exons){
			printf STDERR "\t\t%s\t$params{slice_format}\t%u bases\n",$exon->stable_id,$exon->slice->seq_region_name,$exon->start,$exon->end,strand_symbol($exon->strand),bases($exon->start,$exon->end);
		}
	}
    print_log(sprintf("%s has no exons", $gene_name), WARNING) unless @exons;
	return;
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
        print_slice(
            $gene_name.UPSTREAM,
            $gene_name,
            $gene->stable_id,
            $gene->slice->seq_region_name,
            $gene_strand==1?gene->start-$gene_extra_up:$gene->gene_end+1,
            $gene_strand==1?$gene->start-1:$gene->end+$gene_extra_up,
            $gene_strand,
            0,
            0,
            UPSTREAM,
            $gene->biotype,
            $gene->description);
    }
    @transcripts = @{$gene->get_all_Transcripts};
    $genes{$gene->stable_id}[3] = @transcripts;

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
            $gene_strand==1?$gene->end+$gene_extra_down:gene->start-1,
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
	return abs($end-$start)+1;
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
	if(/s*\$([^\s=>#]+){1}[\s=>#]+([^#]+)#*$/){;
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
		if (exists($params{lc$key})) {
			$params{$key} = $value;
		}
		else {
			printf STDERR "%s isn't a setting", $key;
			$counts{'? warning'}++;
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
    my $gene_name = defined($gene->external_name)?$gene->external_name:$gene->stable_id;
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
    if(@union_coding){@$utr5_3prime = $coding_5prime+$up*1;}
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
						defined($gene->external_name)?$gene->external_name:$gene->stable_id.EXON.$exon_number.UTR3,
						defined($gene->external_name)?$gene->external_name:$gene->stable_id,
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
		if ($type =~ /ERROR/i){
			$counts{'! error'}++;
		}
		if ($type =~ /WARNING/i){
			$counts{'? warning'}++;
		}
		$msg = "$type: $msg";
	}
	print STDERR $msg;
	print LOGFILE $msg unless $opts{L};
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
	if(defined($params{biotypes})){
		my @biotypes = split(/[\D]+/,$params{biotypes});
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
			print_log($record_description);
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
					sprintf "%s %s\n%s...%s",
					$seq->display_id, $seq->description,
					substr($seq->seq, 0, 30),
					substr($seq->seq, -30));
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
						substr($seq->seq, 0, 30), substr($seq->seq, -30)));
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
}

sub chr2numchr{
    my $chr = uc(shift);
	unless ($chr =~ /[XYMT0-9]{1,2}/i){
		print $chr;
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
