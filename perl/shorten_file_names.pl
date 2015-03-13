#!/usr/bin/perl -w
use strict;
use warnings;
#my$a = shift;
use File::Find;
my $i = 1;
my $op = 's/Section /§/';
finddepth sub{
	if(length($File::Find::name)){
		my $old = $File::Find::name;
		(my $new =$old) =~ s/Section /§/g;
		$new =~ s/ - / /g;
		#rename($old,$new);
		printf ("%d\t%d\t%s\n", $i++,length($File::Find::name),$new);
	}
#	print "$f\n" if defined(-f $f);
}, shift;
