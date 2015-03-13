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

my $usage = "perl dmd4.pl input_file\ncheck for missing output from gfPcr bed format";
unless (@ARGV || $ARGV[0] ne 'help') {
	print $usage;
	exit;
}
my $gfPcr=1; #in-silico PCR format
my $isPcr=0; #in-silico PCR format

my $id_num = 1;
while(<>){
	chomp;
	if(/chr.+\s+\d+\s+\d+\s+(\d+)_/){
		while($1!=$id_num++){
			printf("%u\n", $id_num-1);
		}
	}
}
