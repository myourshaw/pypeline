#!/usr/bin/perl -w
use strict;
use warnings;

my @ampStarts;
my @ampEnds;
my @ampIds;
my $partStart = 0;
my $partEnd = 0;
my $ampCount = 0;
my @partId;
my $isEnd;
my $pos;
my %ids;
my $i = 0;
print "ampliconPartId\tchromStart\tchromEnd\tampliconCount\tampliconIds\n";
while(<>){
	chomp;
	if(/(\d+)\s(\d+)\s(\d+)/){
		my ($isEnd,$pos,$ampId) = ($1,$2,$3);
		if(keys(%ids)){
			my @keys = keys %ids;
			if($partStart != $pos){
				printf "%u\t%u\t%u\t%u\t%s\n", 633-$i++, $partStart, $pos, scalar keys(%ids), join(",",sort @keys);
			}
			$partStart = $pos;
			if($isEnd){
				delete($ids{$ampId});
			}
			else{
				$ids{$ampId}++;
			}
		}
		else{
			$partStart = $pos;
			$ids{$ampId}++;
		}
	}
}
