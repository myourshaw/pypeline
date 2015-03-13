#!/usr/bin/perl -w
use strict;
use warnings;
my $ei = "e";
my $i = 11;
my $id;
open OUT, ">", 'C:\Users\Michael\Documents\Lab\DMD\paper\DMD_primers.bed';
print OUT 'track name="PCR primers" description="DMD Raindance PCR primers" group="user" priority=3 visibility="dense" color=255,100,100'."\n";
while(<>){
	my $line = $_;
	if(/>chrX:(\d+)-(\d+) (\d{4})_\S+ \d+bp ([ACTG]+) ([ACTG]+)/){
		printf OUT "chrX\t%u\t%u\tfwd_primer:%s\t1000\t\+\n",$1-1,$1+length($5)-1,$3;
		printf OUT "chrX\t%u\t%u\trev_primer:%s\t1000\t-\n",$2-length($4),$2,$3;
	}
}
