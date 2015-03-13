#!/usr/bin/perl -w
use strict;
use warnings;

#./gfPcr -maxSize=10000 -out=bed localhost 17779 /data/genomes /home/myourshaw/Documents/Lab/DMD/DMD_roi.bed.primer3in_384N.txt.out.gfPcrin /home/myourshaw/Documents/Lab/DMD/DMD_roi.bed.primer3in_384N.txt.out.gfPcrin.out.bed
#./gfPcr -maxSize=10000 -out=fa localhost 17779 /data/genomes /home/myourshaw/Documents/Lab/DMD/DMD_roi.bed.primer3in_384N.txt.out.gfPcrin /home/myourshaw/Documents/Lab/DMD/DMD_roi.bed.primer3in_384N.txt.out.gfPcrin.out.fa

use POSIX qw(ceil floor);
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
#BioPerl
use Bio::Seq;
use Bio::SeqIO;
#EnsEMBL Perl API
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::Utils::Sequence qw (sequence_with_ambiguity);

my $usage = "perl dmd7.pl input_file\nsubstitute primers from primer3 output";
unless (@ARGV || $ARGV[0] ne 'help') {
	print $usage;
	exit;
}
my $primer="";
while(<>){
	if(/^=$/){
		$primer = "";
		print;
	}
	elsif(/^(\d)$/){
		$primer = $1;
	}
	elsif(/^(PRIMER_\S+)_(\d)(_|=)(.+)/i){
		if($2 eq $primer){
			print "$1$3$4\n";
		}
		else{
			print;
		}
	}
	elsif(/^PRIMER_LEFT|PRIMER_RIGHT|PRIMER_PAIR|PRIMER_PRODUCT_SIZE/){
		if($primer eq ""){
			print;
		}
	}
	else{
		print;
	}
}
