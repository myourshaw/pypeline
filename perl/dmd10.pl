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

my $usage = "perl dmd10.pl input_file\nrenumber raindance primers";
unless (@ARGV || $ARGV[0] ne 'help') {
	print $usage;
	exit;
}
my $i = 1;
while(<>){
	unless(/Amplicon_Name/){
		if(/(\d+)(_.+)/){
			my($num,$rest)=($1,$2);
			printf("%04u%s\n",$i++,$rest)
		}
	}
	else{
		print;
	}
}
