#!/usr/bin/perl -w

#"C:\bio\primer3\bin\primer3_core.exe" < "C:\Users\Michael\Documents\Lab\DMD\DMD_roi.bed.primer3in" > "C:\Users\Michael\Documents\Lab\DMD\DMD_roi.bed.primer3out"
#"C:\bio\primer3\bin\primer3_core.exe" -format_output < "C:\Users\Michael\Documents\Lab\DMD\DMD_roi.bed.primer3in" > "C:\Users\Michael\Documents\Lab\DMD\DMD_roi.bed.primer3out"

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

my $usage = "perl dmd6.pl input_file.bed\nList primer3 batch file positions";
unless (@ARGV || $ARGV[0] ne 'help') {
	print $usage;
	exit;
}
my $max_product = 1000;
my $extra_target = 74; #splice sites, etc.
my $max_primer = 28;
my $guard_band = $max_primer+25; #coverage dip next to twin tower
my $primer_band = 2.1*$max_primer; #a place to find a primer
#my $max_target = $max_product-2*$guard_band-2*$primer_band;
#my $min_window = 2*$max_primer;
#my $max_window = $primer_band;
my $id_num = 1;
while(<>){
	if(/chr([0-9MXY]{1,2})\s+(\d+)\s+(\d+)\s+(\S+)\s+\d+\s+([-+])/i){
		my ($chr,$chromStart,$chromEnd,$id,$strand_symbol)=($1,$2,$3,$4,$5);
		my $target_chromStart = $chromStart-$extra_target;
		my $target_chromEnd = $chromEnd+$extra_target;
		my $target_len = $target_chromEnd-$target_chromStart;
		my $strand = $strand_symbol eq "+" ? 1 : -1;
		my $i=0;
		my $this_target_chromStart = $target_chromStart;
		while($this_target_chromStart<$target_chromEnd){
			my $this_target_len = min($max_product-2*$primer_band-2*$guard_band,$target_chromEnd-$this_target_chromStart,ceil($target_len/2));
			my $this_target_chromEnd = $this_target_chromStart+$this_target_len;
			my $offset = $max_product;#-$this_target_len-2*$guard_band-$max_primer;
			my $window_chromStart = $this_target_chromStart-$offset;
			my $window_chromEnd = $this_target_chromEnd+$offset;
			my $window_len = $window_chromEnd-$window_chromStart;
			my $PRIMER_SEQUENCE_ID=sprintf("%04u_%s_%s:%u%s%u",$id_num++,$id,$chr,$window_chromStart+1,$strand_symbol,$window_chromEnd);
			my $line = sprintf("%s\t%u\t%u\t%u\t%u\t%u\t%u\t\n",$PRIMER_SEQUENCE_ID,$this_target_chromStart,$this_target_chromEnd,$this_target_len,$window_chromStart,$window_chromEnd,$window_len);
			print $line;
			$this_target_chromStart=$this_target_chromEnd+1;
		}
	}
}
1;
