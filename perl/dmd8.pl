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

my $usage = "perl dmd8.pl input_file\ncheck for duplicate or missing raindance primers";
unless (@ARGV || $ARGV[0] ne 'help') {
	print $usage;
	exit;
}
my %both_primers;
my %primer_names;
my %missing_primers;
while(<>){
	unless(/Amplicon_Name/){
		if(/(\S+)\s+([ACGT]+)\s+([ACGT]+)/){
			my($Amplicon_Name,$Fwd_Primer_sequence,$Rev_Primer_sequence)=($1,$2,$3);
			$both_primers{"$Fwd_Primer_sequence|$Rev_Primer_sequence"}++;
			$primer_names{"$Fwd_Primer_sequence|$Rev_Primer_sequence"}.="$Amplicon_Name|";
			unless ($Fwd_Primer_sequence && $Rev_Primer_sequence){
				$missing_primers{$Amplicon_Name}++;
			}
		}
	}
}
if(%missing_primers){
	printf "\n%u MISSING PRIMERS\n",scalar keys %missing_primers ;
	print map {"$_\n"} sort keys %missing_primers;
}
if(%both_primers){
	printf "\n%u DUPLICATE PRIMER PAIRS\n",scalar keys %both_primers;
	print map {$_&&$both_primers{$_}>1?"$primer_names{$_}\t$_\t$both_primers{$_}\n":""} sort keys %both_primers;
}
