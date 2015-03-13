#!/usr/bin/perl -w
use strict;
my $file = shift;
my $wig ="$file.wig";
open IN,"<",$file or die "can't open $file. $!";
open WIG,">",$wig or die "can't open $wig. $!";
print WIG 'track type=wiggle_0 name=coverage description="U87 number of reads"'."\n";
my $chr="";
while(<IN>){
    if(/(chr[0-9XY]{1,2})\s+(\d+)\s+\S+\s+(\d+)/){
        print WIG "variableStep chrom=$1\n" unless $1 eq $chr;
        $chr = $1;
        print WIG "$2\t$3\n";
    }
}

