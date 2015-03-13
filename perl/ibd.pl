#!/usr/bin/perl -w
my %markers=qw(rs1945085 1 D11S937 1 rs1945465 1 rs7102569 1 rs1247726 1 rs1459937 1 rs1459952 1 D11S901 1 rs2077889 1 rs1278402 1);
use strict;
while (<>){
    if(/(\d+)\s+(3\d+)\s+(3\d+)\s+(\S+)\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)/){
        if($markers{$4} && $2 != $3){
            print "wana3\t$1\t$2\t$3\t$4\t$5\t$6\t$7\n";
        }
    }
}



