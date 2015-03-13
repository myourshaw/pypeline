#!/usr/bin/perl -w
# ©2009 Michael Yourshaw All Rights Reserved

use strict;
use warnings;

my $that_id;

while (<>) {
	chomp;
	if($_){
		#omim_results file
		if (/\s*[*+#%]?\d{6}.*;\s*(\S+)/) {
			my $this_id = $1;
			if (defined($that_id)) {
				print "$that_id\n";
			}
			$that_id = $this_id;
		}
		elsif ((/\s*Gene map/i || !$_) && defined($that_id)) {
			print "$that_id\n";
			undef $that_id;
		}
		elsif (eof && defined($that_id)) {
			if(defined($that_id)){
				print "$that_id\n";
			}
			undef $that_id;
		}
	}
}
