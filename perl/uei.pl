#!/usr/bin/perl -w
use strict;
use warnings;
my $ei = "e";
my $i = 11;
my $id;
open OUT, ">", 'C:\Users\Michael\Documents\Lab\DMD\paper\DMD_004006_ueiN.fasta';
while(<>){
	my $line = $_;
	if(/NM_004006_\d+/){
		$id = "NM_004006_".$ei.$i;
		$line =~s/NM_004006_\d+/$id/;
		print $line;
		if($ei eq "e"){
			$ei = "i";
		}
		else{
			$ei="e";
			$i++;
		}
	}
	print OUT $line;
}
