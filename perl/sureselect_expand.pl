#!/usr/bin/perl -w
use strict;
my $file = shift;
my $expandedfile = "$file.expanded";
open IN,"<",$file or die;
open OUT,">",$expandedfile or die;
while(my $line = <IN>){
    my @pos = split /\t/,$line;
    for ($pos[1] .. $pos[2]){
        print OUT "$pos[0]\t$_\n";
    }
}


