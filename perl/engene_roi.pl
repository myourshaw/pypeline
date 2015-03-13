#!/usr/bin/perl -w
# ©2007,2008 Michael Yourshaw All Rights Reserved
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

use constant BIG_NUMBER => 1000000000000000;
#the following must have unique values
use constant UTR5 => "5utr";
use constant UTR3 => "3utr";
use constant CDS => "cds";
use constant INTRON => "i";
use constant EXON => "e";
use constant GENE => "gene";
use constant INTERGENETIC => "intergenetic";
use constant TRANSCRIPT => "tr";
use constant REGION => "region";
use constant UPSTREAM => "up";
use constant DOWNSTREAM => "down";


my @start_time = Today_and_Now();

#SEZ6_e1-35+20|chr17:24356939-24357324[-]|386b|SEZ6|ENSG00000063015|roi:24356939-24357324|386b|pad-0+0|CDS,UTR|protein_coding|seizure related 6 homolog isoform 1 [Source:RefSeq_peptide;Acc:NP_849191]
#fragment id|coordinates[strand]|length|parent gene|Ensembl stable id|region of interest coordinates|region of interest length|padding-upstream+downstream|fragment type(s)|gene biotype|gene description
#fragment id = gene name or region coordinates[_up|_e#[_5utr|_3utr|_utr]|_i#|_down ][-additional upstream roi][+additional downstream roi][_fragment part/total fragment parts][_rc]
#up=upstream, e=exon, utr=untranslated region, i=intron, down=downstream, rc=reverse complement
#fragment types = REGION|INTERGENETIC|GENE|UPSTREAM|CDS[,UTR]|UTR|INTRON|DOWNSTREAM
#biotype = C_segment|D_segment|J_segment|miRNA|miRNA_pseudogene|misc_RNA|misc_RNA_pseudogene|Mt_rRNA|Mt_tRNA_pseudogene|protein_coding|pseudogene|repeat|retrotransposed|rRNA|rRNA_pseudogene|scRNA|scRNA_pseudogene|snoRNA|snoRNA_pseudogene|snRNA|snRNA_pseudogene|tRNA_pseudogene|V_segment

#parameter keys must be lower case
my %params = (
	add_bases_downstream_and_upstream_to_region_of_interest_to_get_minimum_length => 1, #bool (default)
	add_bases_downstream_to_region_of_interest_to_get_minimum_length => 0, #bool
	add_bases_upstream_to_region_of_interest_to_get_minimum_length => 0, #bool
	add_padding_bases_downstream_and_upstream_to_record_to_get_minimum_length => 1, #bool (default)
	add_padding_bases_downstream_to_record_to_get_minimum_length => 0, #bool
	add_padding_bases_upstream_to_record_to_get_minimum_length => 0, #bool
	bases_region_of_interest_3_prime_downstream_of_exons => 20,
	bases_region_of_interest_3_prime_downstream_of_introns => 0,
	bases_region_of_interest_5_prime_upstream_of_exons => 20,
	bases_region_of_interest_5_prime_upstream_of_introns => 0,
	bases_region_of_interest_gene_downstream_region => 0,
	bases_region_of_interest_gene_upstream_promoter_region => 100,
	biotypes => "",
	create_log_file => 1, #bool
	create_padded_genomic_intervals_file => 0, #bool (default = 0)
	create_regions_of_interest_genomic_intervals_file => 1, #bool (default = 1)
	create_sequence_file => 1, #bool
	exclude_region_of_interest_less_than_bases => 0,
	hard_mask_repeats => 0, #bool
	include_3_prime_utr_exons => 1, #bool
	include_5_prime_utr_exons => 1, #bool
	include_all_elements_one_record_per_gene => 0, #bool overrides all includes & separates
	include_cds_exons => 1, #bool
	include_gene_downstream_region => 1,
	include_gene_upstream_promoter_region => 1,
	include_intergenetic_sequence_in_genomic_regions => 0, #bool
	include_introns => 0, #bool
	include_rna_like_strand => 0, #bool
	include_template_strand => 1, #bool
	input_format => "", # '' (default) | omim_result gene_result
	mask_repeats => 1, #bool
	maximum_record_length => 1000000000000000, #split records if necessary
	maximum_region_of_interest_per_record => 1000000000000000,
	minimum_record_length => 0, #pad records as necessary
	minimum_region_of_interest_per_record => 120,
	padding_bases_3_prime_downstream_of_each_record => 0,
	padding_bases_3_prime_downstream_of_each_region_of_interest => 0,
	padding_bases_5_prime_upstream_of_each_record => 0,
	padding_bases_5_prime_upstream_of_each_region_of_interest => 0,
	separate_records_for_exons_introns => 1, #bool
	separate_records_for_utr_cds => 1, #bool
	separate_records_for_promoter_upstream_downstream_regions => 1, #bool
	sequence_file_format => "fasta", #fasta (default) | piecemaker | Bio::SeqIO formats
	show_variations => 1, #bool
	single_fasta_identifier => 0, #bool
	slice_format => 'chr%s:%d-%d[%s]' #e.g., 'chr%s:%d-%d[%s]'
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

my $usage = "engene usage:\n"
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
my $repeat_adaptor = $registry->get_adaptor('human', 'core', 'repeatfeature');
my $slice_adaptor = $registry->get_adaptor('human', 'core', 'slice');
my $transcript_adaptor = $registry->get_adaptor('human', 'core', 'transcript');
#read from a list of filenames on the command line
#read from SDTIN if no filenames were given
#"-" indicates STDIN
#"someprogram |" indicates the output of another program (7.14)
#read and parse input file
while (<>) {
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
#all input records done

#genomic interval output
print_all_regions() unless $opts{G};

#print statistics
print_log("lengths of regions of interest (count):");
print LOGFILE map { "$_ ($all_lengths{$_})," } sort {$a <=> $b} keys %all_lengths unless $opts{L};
print LOGFILE "\n" unless $opts{L};
print STDERR map { "$_ ($all_lengths{$_})," } sort {$a <=> $b} keys %all_lengths;
print STDERR "\n";
print_log("gene\tstableID\tbasesROI\tbasesROI+padded");
print LOGFILE map {"$genes{$_}[2]\t$_\t$genes{$_}[0]\t$genes{$_}[1]\n"} sort {$genes{$a}[2] cmp $genes{$b}[2]} keys %genes unless $opts{L};
print STDERR map {"$genes{$_}[2]\t$_\t$genes{$_}[0]\t$genes{$_}[1]\n"} sort {$genes{$a}[2] cmp $genes{$b}[2]} keys %genes;
print_log("gene biotypes:");
print LOGFILE map { " $all_biotypes{$_} $_\n" } sort keys %all_biotypes unless $opts{L};
print STDERR map { " $all_biotypes{$_} $_\n" } sort keys %all_biotypes;
$counts{'ROI minimum base'} = $min_roi;
$counts{'ROI maximum base'} = $max_roi;
print_log("statistics:");
print LOGFILE map { sprintf "%10d %s%s\n",$counts{$_},$_,ss($counts{$_}) } sort keys %counts unless $opts{L};
print STDERR map { sprintf "%10d %s%s\n",$counts{$_},$_,ss($counts{$_}) } sort keys %counts;
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
		print_log(sprintf("%s isn't an exon", $id), "error");
	}
} # !get_by_exon_id

sub get_by_gene_id {
	$id = shift;
	my $gene = $gene_adaptor->fetch_by_stable_id($id);
	if ($gene) {
		gene($gene);
	}
	else {
		print_log(sprintf("%s isn't a gene", $id), "error");
	}
} # !get_by_gene_id

sub get_by_gene_name {
	$id = shift;
	my @genes = @{ $gene_adaptor->fetch_all_by_external_name($id) };
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
			print_log(sprintf("%s identifies more than one gene (%d genes); using %s (%s)",
					$id,
					$count_temp,
					$external_name,
					$stable_id),
				"warning");
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
				"warning"
			);
		}
		gene($gene);
	}
	else {
		print_log(sprintf("%s isn't a gene", $id), "error");
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
		print_log(sprintf("%s isn't a region", $id), "error");
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
		print_log(sprintf("%s isn't a transcript", $id), "error");
	}
} # !get_by_transcript_id

sub gene {
	my $gene = shift;
	#if(!exists($genes{$gene->stable_id})){
		$genes{$gene->stable_id} = [0,0,defined($gene->external_name) ? $gene->external_name : $gene->stable_id]; #roi bases, bases including padding, external name
		$counts{'output gene'}++;
		$all_biotypes{$gene->biotype}++;
		if ($params{include_all_elements_one_record_per_gene}
			||($params{include_3_prime_utr_exons}
			 &&$params{include_5_prime_utr_exons}
			 &&$params{include_cds_exons}
			 &&$params{include_introns}
			 &&!$params{separate_records_for_exons_introns}
			 &&!$params{separate_records_for_utr_cds})) {
			my @exons = @{ $gene->get_all_Exons };
			my $strand = $exons[0]->strand;
			print_slice(
				defined($gene->external_name) ? $gene->external_name : $gene->stable_id,
				defined($gene->external_name) ? $gene->external_name : $gene->stable_id,
				$gene->stable_id,
				$gene->slice->seq_region_name,
				$gene->start,
				$gene->end,
				$strand,
				$params{additional_region_of_interest_5_prime_upstream_of_first_exon},
				$params{bases_additional_region_of_interest_3_prime_downstream_of_exons},
				GENE,
				$gene->biotype,
				$gene->description);
		}
		else {
			my @transcripts = @{ $gene->get_all_Transcripts };
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
			else {
				print_log(sprintf("%s is a pseudogene", $id), "warning");
			}
		}
	#}
	#else{
		#print_log(sprintf("duplicate gene %s (%s) ignored", defined($gene->external_name) ? $gene->external_name : $gene->stable_id, $gene->stable_id), "warning");
	#}
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
		print_log(sprintf("%s has no exons", $id), "warning");
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
	my $years = $_[0]>0 ? sprintf("%d year%s ",$_[0],ss($_[0])):"";
	my $months = $_[1]>0 ? sprintf("%d year%s ",$_[1],ss($_[1])):"";
	my $days = $_[2]>0 ? sprintf("%d year%s ",$_[2],ss($_[2])):"";
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
					printf GENOMICINTERVALSFILE "chr%s:%d-%d\n", numchr2chr($that_chr), $chr_region->start, $chr_region->end;
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
			printf GENOMICINTERVALSFILE "chr%s:%d-%d\n", numchr2chr($that_chr), $chr_region->start, $chr_region->end;
		}
	}
} # !print_all_regions

sub upmost{
	return pop()==1?min(@_):max(@_);
}
sub downmost{
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
sub print_gene_regions {
	my $gene = shift;
	#merge intersecting exons and intersecting coding regions
	my @union_exons = Bio::Range->unions(@all_exons);
	my @union_coding = Bio::Range->unions(@all_coding);
	my $strand = $union_exons[0]->strand;
	my $up = $strand == 1 ? -1 : 1;
	my $down = $strand == 1 ? 1 : -1;
	my $gene_first;
	my $gene_last;
	my $exons_first;
	my $exons_last;
	my $coding_first;
	my $coding_last;
	my $utr5_last;
	my $utr3_first;
	if ($strand==1) {
		@union_exons = sort { $a->start <=> $b->start || $a->end <=> $b->end } @union_exons;
		@union_coding = sort { $a->start <=> $b->start || $a->end <=> $b->end } @union_coding;
		$gene_first = $gene->start;
		$gene_last = $gene->end;
		$exons_first = $union_exons[0]->start;
		$exons_last = $union_exons[$#union_exons]->end;
		$coding_first = $union_coding[0]->start;
		$coding_last = $union_coding[$#union_coding]->end;
		$utr5_last = $coding_first+$up*1;
		$utr3_first = $coding_last+$down*1;
	}
	else {
		@union_exons = sort { $b->start <=> $a->start || $b->end <=> $a->end } @union_exons;
		@union_coding = sort { $b->start <=> $a->start || $b->end <=> $a->end } @union_coding;
		$gene_first = $gene->end;
		$gene_last = $gene->start;
		$exons_first = $union_exons[0]->end;
		$exons_last = $union_exons[$#union_exons]->start;
		$coding_first = $union_coding[0]->end;
		$coding_last = $union_coding[$#union_coding]->start;
		$utr5_last = $coding_first+$up*1;
		$utr3_first = $coding_last+$down*1;
	}
	#make an array of all the ranges in the gene that might become separate records
	my @print_infos;
	my $upstream_region;
	my $offset = "";
	my $upstream_region_first = upmost($gene_first,$exons_first+$up*$params{bases_region_of_interest_gene_upstream_promoter_region},$strand);
	my $upstream_region_last = $exons_first + $up*1;
	my $element_range = get_bio_range($upstream_region_first,$upstream_region_last,$strand);
	my $roi_range = get_bio_range($upstream_region_first,$upstream_region_last,$strand);
	my $upstream_info;
	if($params{include_gene_upstream_promoter_region} && $upstream_region_first != $exons_first){
		$upstream_info = {element_range=>$element_range, roi_range=>$roi_range, type=>UPSTREAM, up=>0, down=>0, number=>UPSTREAM};
	}
	my $that_exon_last;
	my $exon_number=0;
	foreach my $exon(@union_exons){
		my $exon_first = upmost($exon->start,$exon->end,$strand);
		my $exon_last = downmost($exon->start,$exon->end,$strand);
		if($params{include_introns} && $exon_number>0 && is_downstream($exon_first,$that_exon_last,$strand)){
			my $extra_up = $up*$params{bases_region_of_interest_5_prime_upstream_of_introns};
			my $extra_down = $down*$params{bases_region_of_interest_3_prime_downstream_of_introns};
			my $element_range = get_bio_range($that_exon_last+$down*1,$exon_first+$up*1,$strand);
			my $roi_range = get_bio_range($that_exon_last+$down*1+$extra_up,$exon_first+$up*1+$extra_down,$strand);
			push @print_infos,{element_range=>$element_range, roi_range=>$roi_range, type=>INTRON, up=>$extra_up, down=>$extra_down, number=>INTRON.$exon_number};
		}
		$exon_number++;
		$that_exon_last = $exon_last;
		my @exon_infos;
		my @exon_ranges;
		my $utr5_info;
		my $cds_info;
		my $utr3_info;
		if($params{include_5_prime_utr_exons} && is_upstream($exon_first,$coding_first,$strand)){
			my $extra_up = $up*$params{bases_region_of_interest_5_prime_upstream_of_exons};
			my $extra_down = $exon_last==$utr5_last?$down*$params{bases_region_of_interest_3_prime_downstream_of_exons}:0;
			my $offset = ($extra_up?"+".$extra_up:"").($extra_down?"-".$extra_down:"");
			my $element_range = get_bio_range($exon_first,upmost($exon_last,$utr5_last,$strand),$strand);
			my $roi_range = get_bio_range($exon_first+$extra_up,upmost($exon_last,$utr5_last,$strand)+$extra_down,$strand);
			push @exon_ranges,$roi_range;
			#$utr5_info = {element_range=>$element_range, roi_range=>$roi_range, type=>(UTR5), up=>$extra_up, down=>$extra_down, number=>EXON.$exon_number.UTR5};
			push @print_infos, {element_range=>$element_range, roi_range=>$roi_range, type=>UTR5, up=>$extra_up, down=>$extra_down, number=>EXON.$exon_number.UTR5};
		}
		if($params{include_cds_exons} && is_upstream($exon_first,$coding_last,$strand) && is_downstream($exon_last,$coding_first,$strand)){
			my $extra_up = $coding_first==$exon_first?$up*$params{bases_region_of_interest_5_prime_upstream_of_exons}:0;
			my $extra_down = $coding_last==$exon_last?$down*$params{bases_region_of_interest_3_prime_downstream_of_exons}:0;
			my $offset = ($extra_up?"+".$extra_up:"").($extra_down?"-".$extra_down:"");
			my $element_range = get_bio_range(downmost($coding_first,$exon_first,$strand),upmost($coding_last,$exon_last,$strand),$strand);
			my $roi_range = get_bio_range(downmost($coding_first,$exon_first,$strand)+$extra_up,upmost($coding_last,$exon_last,$strand)+$extra_down,$strand);
			push @exon_ranges,$roi_range;
			#$cds_info = {element_range=>$element_range, roi_range=>$roi_range, type=>(CDS), up=>$extra_up, down=>$extra_down, number=>EXON.$exon_number.CDS};
			push @print_infos, {element_range=>$element_range, roi_range=>$roi_range, type=>CDS, up=>$extra_up, down=>$extra_down, number=>EXON.$exon_number.CDS};
		}
		if($params{include_3_prime_utr_exons} && is_downstream($exon_last,$coding_last,$strand)){
			my $extra_up = $utr3_first==$exon_first?$up*$params{bases_region_of_interest_5_prime_upstream_of_exons}:0;
			my $extra_down = $down*$params{bases_region_of_interest_3_prime_downstream_of_exons};
			my $offset = ($extra_up?"+".$extra_up:"").($extra_down?"-".$extra_down:"");
			my $element_range = get_bio_range(downmost($utr3_first,$exon_first,$strand),$exon_last,$strand);
			my $roi_range = get_bio_range(downmost($utr3_first,$exon_first,$strand)+$extra_up,$exon_last+$extra_down,$strand);
			push @exon_ranges,$roi_range;
			#$utr3_info = {element_range=>$element_range, roi_range=>$roi_range, type=>(UTR3), up=>$extra_up, down=>$extra_down, number=>EXON.$exon_number.UTR3};
			push @print_infos, {element_range=>$element_range, roi_range=>$roi_range, type=>UTR3, up=>$extra_up, down=>$extra_down, number=>EXON.$exon_number.UTR3};
		}
#		my @union_exon_ranges = Bio::Range->unions(@exon_ranges);
		if($params{separate_records_for_utr_cds} || (defined($utr5_info)&&!defined($cds_info)&&defined($utr3_info))){
			push @print_infos, $utr5_info if defined($utr5_info);
			push @print_infos, $cds_info if defined($cds_info);
			push @print_infos, $utr3_info if defined($utr3_info);
		}
		else{
			my $extra_up = $up*$params{bases_region_of_interest_5_prime_upstream_of_exons};
			my $extra_down = $down*$params{bases_region_of_interest_3_prime_downstream_of_exons};
			my $element_range = get_bio_range($exon_first,$exon_last,$strand);
			my $roi_range = get_bio_range($exon_first+$extra_up,$exon_last+$extra_down,$strand);
			my @exon_types;
			push @exon_types, UTR5 if defined($utr5_info);
			push @exon_types, CDS if defined($cds_info);
			push @exon_types, UTR3 if defined($utr3_info);
			my $offset = ($extra_up?"+".$extra_up:"").($extra_down?"-".$extra_down:"");
			push @print_infos, {element_range=>$element_range, roi_range=>$roi_range, type=>\@exon_types, up=>$extra_up, down=>$extra_down, number=>EXON.$exon_number.join("&",@exon_types)};
		}
	}
	$offset = "";
	my $downstream_region_first = $exons_last + $down*1;
	my $downstream_region_last = downmost($gene_last,$exons_last+$down*$params{bases_region_of_interest_gene_downstream_region},$strand);
	$element_range = get_bio_range($downstream_region_first,$downstream_region_last,$strand);
	$roi_range = get_bio_range($downstream_region_first,$downstream_region_last,$strand);
	my $downstream_info;
	if($params{include_gene_downstream_region} && $downstream_region_last != $exons_last){
		$downstream_info = {element_range=>$element_range, roi_range=>$roi_range, type=>DOWNSTREAM, up=>0, down=>0, number=>DOWNSTREAM};
	}
	unless($params{separate_records_for_promoter_upstream_downstream_regions}){
		if($upstream_info){
			my ($roi_first,$roi_last) =  get_5prime_3prime($print_infos[0]->{roi_range});
			unless(upstream($upstream_region_last,add_upstream($roi_first,1,$strand),$strand)){
				my $first_element = shift @print_infos;
				push my @new_roi_range,$first_element->{roi_range};
				push @new_roi_range,$upstream_info->{roi_range};
				@new_roi_range = Bio::Range->unions(@new_roi_range);
				push my @new_element_range,$first_element->{element_range};
				push @new_element_range,$upstream_info->{element_range};
				@new_element_range = Bio::Range->unions(@new_element_range);
			}
			else{
				unshift @print_infos,$upstream_info if($upstream_info);
			}
		}
		if($downstream_info){
			my $last_element = $print_infos[$#print_infos];
		}
	}
	else{
		unshift @print_infos,$upstream_info if($upstream_info);
		push @print_infos,$downstream_info if($downstream_info);
	}
	my @print_slices;
	my @all_roi_ranges;
	my @all_types;
	foreach my $print_info(@print_infos){
		if ($print_info->{roi_range}){
			push @all_roi_ranges, $print_info->{roi_range};
			push @all_types, $print_info->{type}
		}
	}
	#collapse overlapping ranges
	my @union_all_bio_ranges = Bio::Range->unions(@all_roi_ranges);
	unless($params{separate_records_for_exons_introns} || $params{separate_records_for_promoter_upstream_downstream_regions} || @union_all_bio_ranges>1){
		my ($element_first, $foo) = get_5prime_3prime($print_infos[0]->element_range);
		my ($bar, $element_last) = get_5prime_3prime($print_infos[$#print_infos]->element_range);
		my $element_range = get_bio_range($element_first,$element_last,$strand);
		my $extra_up = abs($element_first-$gene_first);
		my $extra_down = abs($element_last-$gene_last);
		my $offset = ($extra_up?"+".$extra_up:"").($extra_down?"-".$extra_down:"");
		push @print_slices, {element_range=>$element_range, roi_range=>$roi_range, type=>@all_types, up=>$extra_up, down=>$extra_down, number=>""};
	}
	#print gene as one record if all ranges are contiguous
	unless($params{separate_records_for_exons_introns}){
		my @bio_ranges;
		my @range_types;
		my $up = 0;
		my $down = 0;
		foreach my $print_info(@print_infos){
			push @bio_ranges, $print_info->{roi_range};
			$up = $up?$up:$print_info->{up};
			$down= $print_info->{down};
			push @range_types, $print_info->{type};
		}
		@bio_ranges = Bio::Range->unions(@bio_ranges);
		if (@bio_ranges == 1 || 1){
			my $range = {range=>shift(@bio_ranges),type=>GENE,up=>$up,down=>$down};
			push @print_slices, {range=>$range,record_number=>""};
		}
	}
	#print gene as multiple records
	unless(@print_slices){
		for(my $exon_number=1;$exon_number<=@print_infos;$exon_number++){
			my $print_info = shift @print_infos;
			if($print_info->{type} eq UPSTREAM){
				push @print_slices, {range=>$print_info,record_number=>"_".EXON.$exon_number.$print_info->{type}};
			}
			elsif($print_info->{type} eq INTRON){
				push @print_slices, {range=>$print_info,record_number=>"_".INTRON.$exon_number-1};
			}
			elsif($print_info->{type} eq UPSTREAM){
				push @print_slices, {range=>$print_info,record_number=>"_".EXON.$exon_number.$print_info->{type}};
			}
			elsif($print_info->{type} eq UTR5 || $print_info->{type} eq CDS || $print_info->{type} eq UTR3){
				my @exon_print_infos;
				if($print_info->{type} eq UTR5){
					if($params{separate_records_for_utr_cds}){
						push @print_slices, {range=>$print_info,record_number=>"_".EXON.$exon_number.$print_info->{type}};
					}
					else{
						push @exon_print_infos,$print_info;
					}
				}
				if($print_infos[0]->{type} eq CDS){
					$print_info = shift @print_infos;
				}
				if($print_info->{type} eq CDS){
					if($params{separate_records_for_utr_cds}){
						push @print_slices, {range=>$print_info,record_number=>"_".EXON.$exon_number.$print_info->{type}};
					}
					else{
						push @exon_print_infos,$print_info;
					}
				}
				if($print_infos[0]->{type} eq UTR3){
					$print_info = shift @print_infos;
				}
				if($print_info->{type} eq UTR3){
					if($params{separate_records_for_utr_cds}){
						push @print_slices, {range=>$print_info,record_number=>"_".EXON.$exon_number.$print_info->{type}};
					}
					else{
						push @exon_print_infos,$print_info;
					}
				}
				if (@exon_print_infos){
					my @bio_ranges;
					my @exon_types;
					my $exon_up = 0;
					my $exon_down = 0;
					foreach my $exon_print_info(@exon_print_infos){
						push @bio_ranges, $exon_print_info->{range};
						$exon_up = $up?$up:$exon_print_info->{up};
						$exon_down= $exon_print_info->{down};
						push @exon_types, $exon_print_info->{type};
					}
					@bio_ranges = Bio::Range->unions(@bio_ranges);
					if(@bio_ranges == 1){
						my $roi_range = {range=>shift(@bio_ranges),type=>join(",",@exon_types),up=>$exon_up,down=>$exon_down};
						push @print_slices, {range=>$roi_range,record_number=>"_".EXON.$exon_number.join("_",@exon_types)};
					}
					else{
						foreach my $roi_range(@bio_ranges){
							push @print_slices, {range=>$roi_range,record_number=>"_".EXON.$exon_number.shift(@exon_types)};
						}
					}
				}
			}
			elsif($print_info->{type} eq DOWNSTREAM){
				push @print_slices, {range=>$print_info,record_number=>$print_info->{type}};
			}
		}
	}
	foreach my $print_slice(@print_slices){
		my $id = (defined($gene->external_name) ? $gene->external_name : $gene->stable_id) . $print_slice->{record_number};
		my $gene_external_name = defined($gene->external_name) ? $gene->external_name : $gene->stable_id;
		my $gene_stable_id = $gene->stable_id;
		my $chr = $gene->slice->seq_region_name;
		my $roi_start = $print_slice->start;
		my $roi_end = $print_slice->end;
		my $strand = $print_slice->slice;
		my $roi_minus = $print_slice->up;#bases in roi 5' (+ strand) or 3' (- strand) of range implied by $id
		my $roi_plus = $print_slice->down;#bases in roi 3' (+ strand) or 5' (- strand) of range implied by $id
		my $type = $print_slice->type;
		my $biotype = $gene->biotype;
		my $description = $gene->description;
		print_slice($id,$gene_external_name,$gene_stable_id,$chr,$roi_start,$roi_end,$strand,$roi_minus,$roi_plus,$type,$biotype,$description);
	}
} # !print_gene_regions

sub print_log {
	my ($msg, $type) = @_;
	$msg = "$msg\n";
	if ($type) {
		$type = uc($type);
		if ( $type =~ /ERR/) {
			$counts{'! error'}++;
		}
		if ( $type =~ /WARN/) {
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
		$roi_minus, #bases in roi 5' (+ strand) or 3' (- strand) of range implied by $id
		$roi_plus, #bases in roi 3' (+ strand) or 5' (- strand) of range implied by $id
		$type,
		$biotype,
		$description) = @_;
	if(bases($roi_start,$roi_end) < $params{exclude_region_of_interest_less_than_bases}){
		return;
	}
	if(defined($params{biotypes})){
		my @biotypes = split(/[\D]+/,$params{biotypes});
		if(@biotypes && $biotypes[0] && defined($biotype) && !contains(\@biotypes,$biotype)) {
			return;
		}
	}
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
		$roi_minus, #bases in roi 5' (+ strand) or 3' (- strand) of range implied by $id
		$roi_plus, #bases in roi 3' (+ strand) or 5' (- strand) of range implied by $id
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
		$roi_minus, #bases in roi 5' (+ strand) or 3' (- strand) of range implied by $id
		$roi_plus, #bases in roi 3' (+ strand) or 5' (- strand) of range implied by $id
		$type,
		$biotype,
		$description) = @_;
	$roi_minus += $params{padding_bases_5_prime_upstream_of_each_region_of_interest};
	$roi_plus += $params{padding_bases_3_prime_downstream_of_each_region_of_interest};
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
					$i == 1 ? $roi_minus : 0,
					$i == $count ? $roi_plus : 0,
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
					$i == 1 ? 0 : $roi_plus,
					$i == $count ? 0 : $roi_minus,
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
		if ($params{add_bases_upstream_to_region_of_interest_to_get_minimum_length}) {
			if ($strand == 1) {
				$roi_start -= $extra;
				$roi_minus += $extra;
			}
			else {
				$roi_end += $extra;
				$roi_plus += $extra;
			}
		}
		elsif ($params{add_bases_downstream_to_region_of_interest_to_get_minimum_length}) {
			if ($strand == 1) {
				$roi_end += $extra;
				$roi_plus += $extra;
			}
			else {
				$roi_start -= $extra;
				$roi_minus += $extra;
			}
		}
		else { #center
			my $extra1 = int($extra / 2);
			my $extra2 = $extra - $extra1;
			if ($strand == 1) {
				$roi_start -= $extra1;
				$roi_minus += $extra1;
				$roi_end += $extra2;
				$roi_plus += $extra2;
			}
			else {
				$roi_end += $extra1;
				$roi_plus += $extra1;
				$roi_start -= $extra2;
				$roi_minus += $extra2;
			}
		}
	}
	my $slice_start = $roi_start - $strand * $params{padding_bases_5_prime_upstream_of_each_record};
	my $slice_end = $roi_end + $strand * $params{padding_bases_3_prime_downstream_of_each_record};
	#genomic interval records for Agilent eArray probes
	unless($opts{G}){
		my $id = sprintf "%s%s%s",
			$id,
			$roi_minus == 0 ? '' : ('-' . $roi_minus),
			$roi_plus == 0 ? '' : ('+' . $roi_plus);
		my $record_description = sprintf("%s|chr%s:%d-%d[%s]|%db|%s|%s|roi:%d-%d|%db|pad-%d+%d|%s|%s|%s",
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
		push @all_regions, {num_chr=>chr2numchr($chr), start=>$slice_start, end=>$slice_end, $strand=>$strand};
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
					$roi_minus == 0 ? '' : ('-' . $roi_minus),
					$roi_plus == 0 ? '' : ('+' . $roi_plus);
				my $record_description = sprintf("|%s|%db|%s|%s|roi:%d-%d|%db|pad-%d+%d|%s|%s|%s",
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
				print_log(sprintf("Can't get %s", $id), "error");
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
					$roi_minus == 0 ? '' : ('-' . $roi_minus),
					$roi_plus == 0 ? '' : ('+' . $roi_plus);
				my $record_description = sprintf("|%s|%db|%s|%s|roi:%d-%d|%db|pad-%d+%d|%s",
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
				print_log(sprintf("Can't get %s_rc", $id), "error");
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
	print STDERR $usage . "Default settings:\n";
	print STDERR map { "\t\$$_ = $params{$_}\n" } sort keys %params;
}

sub chr2numchr{
    my $chr = uc(shift);
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
