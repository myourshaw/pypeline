#!/usr/bin/perl -w
use strict;

my $f = 'SIFT';
my $s = 'deleterious(0.945)';
if ($f eq 'SIFT' && $s =~ /([^(]+)\(([0-9.-]+)/){
    print $1;
    print $2;
}
else{
    print;
}
