#!/usr/bin/perl -w

#"C:\bio\primer3\bin\primer3_core.exe" < "C:\Users\Michael\Documents\Lab\DMD\DMD_roi.bed_primer3in" > "C:\Users\Michael\Documents\Lab\DMD\DMD_roi.bed_primer3out"
#"C:\bio\primer3\bin\primer3_core.exe" -format_output < "C:\Users\Michael\Documents\Lab\DMD\DMD_roi.bed_primer3in" > "C:\Users\Michael\Documents\Lab\DMD\DMD_roi.bed_primer3out"

use strict;
use warnings;
use POSIX qw(ceil floor);
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
#BioPerl
use Bio::Seq;
use Bio::SeqIO;
#EnsEMBL Perl API
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::Utils::Sequence qw (sequence_with_ambiguity);

my $usage = "perl dmd2.pl input_file.bed\nCreate primer3 batch file from positions in a bed file, using ensembl database\n";
unless (@ARGV || $ARGV[0] ne 'help') {
	print $usage;
	exit;
}
my $no_ensembl = 0;
my $no_out_file = 0;
my $no_primer3_format = 0;
my $max_primer = 28;
my $extra_target = 74; #splice sites, etc.
my $min_pairs = 2;
my $guard_band = $max_primer+25; #coverage dip next to twin tower
my $primer_band = ceil($max_primer*2.1); #a place to find a primer
my $max_product = 1000;
my $max_target=1000;
#my $max_target = $max_product-2*$guard_band-2*$primer_band;
#my $min_window = 2*$max_primer;
#my $max_window = $primer_band;
my $id_num = 1;
my$params = 'PRIMER_EXPLAIN_FLAG=1
PRIMER_MISPRIMING_LIBRARY=C:\bio\primer3\bin\humrep_and_simple.txt
PRIMER_PRODUCT_SIZE_RANGE=150-1000
PRIMER_NUM_RETURN=5
PRIMER_MAX_END_STABILITY=9.0
PRIMER_MAX_MISPRIMING=12.00
PRIMER_PAIR_MAX_MISPRIMING=24.00
PRIMER_MAX_TEMPLATE_MISPRIMING=12.00
PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING=24.00
PRIMER_MIN_SIZE=20
PRIMER_OPT_SIZE=22
PRIMER_MAX_SIZE=28
PRIMER_MIN_TM=56.0
PRIMER_OPT_TM=58.0
PRIMER_MAX_TM=60.0
PRIMER_MAX_DIFF_TM=100.0
PRIMER_TM_SANTALUCIA=1
PRIMER_MIN_GC=30.0
PRIMER_MAX_GC=60.0
PRIMER_SELF_ANY=8.00
PRIMER_SELF_END=3.00
PRIMER_NUM_NS_ACCEPTED=0
PRIMER_MAX_POLY_X=5
PRIMER_OUTSIDE_PENALTY=0
PRIMER_GC_CLAMP=1
PRIMER_SALT_CONC=60.0
PRIMER_SALT_CORRECTIONS=1
PRIMER_DIVALENT_CONC=2.0
PRIMER_DNTP_CONC=0.2
PRIMER_DNA_CONC=50.0
PRIMER_LIBERAL_BASE=1
PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS=0
PRIMER_MIN_QUALITY=0
PRIMER_MIN_END_QUALITY=0
PRIMER_QUALITY_RANGE_MIN=0
PRIMER_QUALITY_RANGE_MAX=100
PRIMER_WT_TM_LT=1.0
PRIMER_WT_TM_GT=1.0
PRIMER_WT_SIZE_LT=1.0
PRIMER_WT_SIZE_GT=1.0
PRIMER_WT_GC_PERCENT_LT=0.0
PRIMER_WT_GC_PERCENT_GT=0.0
PRIMER_WT_COMPL_ANY=0.0
PRIMER_WT_COMPL_END=0.0
PRIMER_WT_NUM_NS=0.0
PRIMER_WT_REP_SIM=0.0
PRIMER_WT_SEQ_QUAL=0.0
PRIMER_WT_END_QUAL=0.0
PRIMER_WT_POS_PENALTY=0.0
PRIMER_WT_END_STABILITY=0.0
PRIMER_WT_TEMPLATE_MISPRIMING=0.0
PRIMER_PAIR_WT_PRODUCT_SIZE_LT=0.0
PRIMER_PAIR_WT_PRODUCT_SIZE_GT=0.0
PRIMER_PAIR_WT_PRODUCT_TM_LT=0.0
PRIMER_PAIR_WT_PRODUCT_TM_GT=0.0
PRIMER_PAIR_WT_DIFF_TM=0.0
PRIMER_PAIR_WT_COMPL_ANY=0.0
PRIMER_PAIR_WT_COMPL_END=0.0
PRIMER_PAIR_WT_REP_SIM=0.0
PRIMER_PAIR_WT_PR_PENALTY=1.0
PRIMER_PAIR_WT_IO_PENALTY=0.0
PRIMER_PAIR_WT_TEMPLATE_MISPRIMING=0.0
PRIMER_INTERNAL_OLIGO_MIN_SIZE=18
PRIMER_INTERNAL_OLIGO_OPT_SIZE=20
PRIMER_INTERNAL_OLIGO_MAX_SIZE=27
PRIMER_INTERNAL_OLIGO_MIN_TM=57.0
PRIMER_INTERNAL_OLIGO_OPT_TM=60.0
PRIMER_INTERNAL_OLIGO_MAX_TM=63.0
PRIMER_INTERNAL_OLIGO_MIN_GC=20.0
PRIMER_INTERNAL_OLIGO_MAX_GC=80.0
PRIMER_INTERNAL_OLIGO_SELF_ANY=12.00
PRIMER_INTERNAL_OLIGO_SELF_END=12.00
PRIMER_INTERNAL_OLIGO_NUM_NS=0
PRIMER_INTERNAL_OLIGO_MAX_POLY_X=5
PRIMER_INTERNAL_OLIGO_MISHYB_LIBRARY=
PRIMER_INTERNAL_OLIGO_MAX_MISHYB=12.00
PRIMER_INTERNAL_OLIGO_MIN_QUALITY=0
PRIMER_INTERNAL_OLIGO_SALT_CONC=50.0
PRIMER_INTERNAL_OLIGO_DNA_CONC=50.0
PRIMER_INTERNAL_OLIGO_DIVALENT_CONC=0.0
PRIMER_INTERNAL_OLIGO_DNTP_CONC=0.0
PRIMER_IO_WT_TM_LT=1.0
PRIMER_IO_WT_TM_GT=1.0
PRIMER_IO_WT_SIZE_LT=1.0
PRIMER_IO_WT_SIZE_GT=1.0
PRIMER_IO_WT_GC_PERCENT_LT=0.0
PRIMER_IO_WT_GC_PERCENT_GT=0.0
PRIMER_IO_WT_COMPL_ANY=0.0
PRIMER_IO_WT_NUM_NS=0.0
PRIMER_IO_WT_REP_SIM=0.0
PRIMER_IO_WT_SEQ_QUAL=0.0
PRIMER_TASK=pick_pcr_primers
PRIMER_PICK_ANYWAY=1
PRIMER_FIRST_BASE_INDEX=1
';
unless($no_out_file){
my $out = $ARGV[0]."_primer3in".".txt";
open OUT, ">", $out or die unless $no_out_file;
print $params;
print OUT $params;
}
while(<>){
	if(/chr([0-9MXY]{1,2})\s+(\d+)\s+(\d+)\s+(\S+)\s+\d+\s+([-+])/i){
		my ($chr,$chromStart,$chromEnd,$id,$strand_symbol)=($1,$2,$3,$4,$5);
		my $intron = 1;#substr($id,0,1) eq "I"?1:0;
		my $extra_target = $intron?0:74; #splice sites, etc.
		my $min_pairs = $intron?1:2;
		my $guard_band = $intron?0:$max_primer+25; #coverage dip next to twin tower
		my $primer_band = $intron?ceil($max_primer*3):ceil($max_primer*2.1); #a place to find a primer
		my $target_chromStart = $chromStart-$extra_target;
		my $target_chromEnd = $chromEnd+$extra_target;
		my $target_len = $target_chromEnd-$target_chromStart;
		my $strand = $strand_symbol eq "+" ? 1 : -1;
		my $i=0;
		my $this_target_chromStart = $target_chromStart;
		while($this_target_chromStart<$target_chromEnd){
			my $this_target_len = min($max_product-2*$primer_band-2*$guard_band,$target_chromEnd-$this_target_chromStart,ceil($target_len/$min_pairs),$max_target);
			my $this_target_chromEnd = $this_target_chromStart+$this_target_len;
			my $offset = $max_product;
			my $window_chromStart = $this_target_chromStart-$offset;
			my $window_chromEnd = $this_target_chromEnd+$offset;
			my $window_len = $window_chromEnd-$window_chromStart;
			my $primer_id = sprintf("%04u_%s_%s:%u%s%u[%u%s%u]",$id_num++,$id,$chr,$window_chromStart+1,$strand_symbol,$window_chromEnd,$this_target_chromStart,$strand_symbol,$this_target_chromEnd);
			my $PRIMER_SEQUENCE_ID=sprintf("PRIMER_SEQUENCE_ID=%s",$primer_id);
			my $PRIMER_FIRST_BASE_INDEX=sprintf("PRIMER_FIRST_BASE_INDEX=%u",$strand==1?$window_chromStart+1:$window_chromEnd);
			my $TARGET=sprintf("TARGET=%u,%u",$offset+1,$this_target_len);
			my $EXCLUDED_REGION=sprintf("EXCLUDED_REGION=%u,%u",$offset-$guard_band,$this_target_len+2*$guard_band);
			my $PRIMER_PRODUCT_OPT_SIZE=sprintf("PRIMER_PRODUCT_OPT_SIZE=%u",$max_product);
			my $SEQUENCE=sprintf("SEQUENCE=%s",fetch_sequence_with_ambiguity($chr,$window_chromStart+1,$window_chromEnd,$strand)) unless $no_ensembl;
			my $line;
			if($no_primer3_format){
				$line = sprintf("%s\t%u\t%u\t%u\t%u\t%u\t%u\t\n",$primer_id,$this_target_chromStart,$this_target_chromEnd,$this_target_len,$window_chromStart,$window_chromEnd,$window_len);
			}
			else{
				$line = sprintf("%s\n%s\n%s\n%s\n%s\n=\n",$PRIMER_SEQUENCE_ID,$TARGET,$EXCLUDED_REGION,$PRIMER_PRODUCT_OPT_SIZE,$SEQUENCE);

			}
			print $line;
			print OUT $line unless $no_out_file;
			$this_target_chromStart=$this_target_chromEnd+1;
		}
	}
}
1;
sub fetch_sequence_with_ambiguity {
	my ($chr, $start, $end, $strand) = @_;
	my $registry = 'Bio::EnsEMBL::Registry';
	$registry->load_registry_from_db(
		-host => 'ensembldb.ensembl.org',
		-user => 'anonymous');
	my $dbCore = $registry->get_DBAdaptor('human', 'core');
	my $dbVar = $registry->get_DBAdaptor('human', 'variation');
	my $ambiguous_slice = sequence_with_ambiguity($dbCore, $dbVar, $chr, $start, $end, $strand);
	###############################################################################################
	# modify Bio::EnsEMBL::Variation::Utils::Sequence.pm sub ambiguity_code to always return "N"; #
	###############################################################################################
	my $ambiguous_seq = $ambiguous_slice->seq();
	return $ambiguous_seq;
} # !fetch_sequence_with_ambiguity
