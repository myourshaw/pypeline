#!/usr/bin/perl -w

use strict;
use warnings;
my $in = shift;
open IN,"<",$in;
open OUT,">","$in.out";
my $lineno = 0;
while(<IN>){
	print OUT;
}