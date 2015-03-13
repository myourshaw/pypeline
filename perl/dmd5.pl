#!/usr/bin/perl -w
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

my $usage = "perl dmd5.pl input_file\ncheck for duplicate output from primer3";
unless (@ARGV || $ARGV[0] ne 'help') {
	print $usage;
	exit;
}
my $gfPcr=1; #in-silico PCR format
my $isPcr=0; #in-silico PCR format

my $id_num = 1;
my ($thatF,$thatR)=('','');
while(<>){
	if(/(\S+)\s+([ACGT]+)\s+([ACGT]+)/){
		if($2 eq $thatF || $3 eq $thatR){
			print;
		}
		$thatF=$2;
		$thatR=$3;
	}
}
