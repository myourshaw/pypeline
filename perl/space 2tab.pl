#!/usr/bin/perl -w
use strict;
use warnings;
open STDOUT, ">", $ARGV[0].".tab" or die "can't open $ARGV[0]: $!";
while(<>){
    #chomp;
    tr/ /\t/;
    print ;
}
