#!/usr/bin/perl -w
use strict;
use warnings;
my (%pct,%cnt);
my $tot = 0;
while(<>){
    if(/([ACGT]{3})\s([0-9.]+)\s(\d+)/i){
        $pct{$1}=$2;
        $cnt{$1}=$3;
        $tot+=$3;
    }
    else{
        if(%pct){
            print map { sprintf("%s\t%8.4f%%\t%d\n",$_,100*$pct{$_},$cnt{$_}) } sort {$cnt{$b} <=> $cnt{$a} || $a cmp $b} (keys %cnt);
            print "total\t$tot\n";
            undef(%pct);
            undef(%cnt);
            $tot = 0;
        }
        print;
    }
}
if(%pct){
    print map { sprintf("%s\t%8.4f%%\t%d\n",$_,100*$pct{$_},$cnt{$_}) } sort {$cnt{$b} <=> $cnt{$a} || $a cmp $b} (keys %cnt);
    print "total\t$tot\n";
}
