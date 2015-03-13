#!/usr/bin/perl -w

use strict;
use warnings;
my $in = shift;
open IN,"<",$in;
open OUT,">","$in.sib";
my $lineno = 0;
while(<IN>){
	$lineno++;
	my @foo=split /,/;
	if(($foo[0] eq "famid") || ($foo[1] >= 300)){
		print OUT;
	}
	else{
		#remove founder phenotypes
		printf OUT ("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,\n",$foo[0],$foo[1],$foo[2],$foo[3],$foo[4],$foo[5],$foo[6],$foo[7],$foo[8],$foo[9],$foo[10]);
	}
}