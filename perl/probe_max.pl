#!/usr/bin/perl -w
use strict;
use warnings;
# "C:\Users\Michael\Documents\Lab\heatmap\lin\x\ASDGENES_CORREL_ALL2.txt" "C:\Users\Michael\Documents\Lab\heatmap\lin\x\gene scoring.txt"
sub trim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}

my ($stan_file,$lin_file) = @ARGV;
open STAN,"<",$stan_file;
open LIN,"<",$lin_file;
my %max;
my %genes;
while(<LIN>){
	my @line = split /\t/;
	unless ($line[0] eq "Gene"){
		$genes{trim($line[0])}=-999;
	}
}
while(<STAN>){
	my @line = split /\t/;
	my @genes = split(/ \/\/\/ /,uc(trim($line[2])));
	my $corr = trim($line[5]);
	unless ($line[1] eq "Alignments"){
		foreach(@genes){
			if(exists($genes{$_}) && $genes{$_}<$corr){
				$genes{$_} = $corr;
		}
		}
	}
}
my $out_file = $lin_file.".txt";
open OUT,">",$out_file;
print map {"$_\t$genes{$_}\n"} sort keys %genes;
