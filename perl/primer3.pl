#!/usr/bin/perl -w
use strict;
use warnings;

use POSIX qw(ceil floor);
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
#BioPerl
use Bio::Seq;
use Bio::SeqIO;
use Bio::Tools::Primer3;
#EnsEMBL Perl API
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::Utils::Sequence qw (sequence_with_ambiguity);

my $usage = "perl primer3 input_file\n";
unless (@ARGV) {
	print $usage;
	exit;
}
