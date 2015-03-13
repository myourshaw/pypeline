#!/usr/bin/perl -w
use strict;
#!/usr/bin/perl -w
use strict;
use warnings;
use integer;

my %params = (
	input_format => "gene_result",
);

my $input_file  = "";
my $output_file = "";


#redirect STDOUT, which the debugger doesn't do
if ( $ARGV[ $#ARGV - 1 ] eq ">" ) {
	$output_file = pop(@ARGV);
	pop(@ARGV);
	open STDOUT, ">", $output_file or die "can't open $output_file: $!";
}

print_log( "Gene List\n" );

##statistics
my $output_record_count = 0;

#read from a list of filenames on the command line
#read from SDTIN if no filenames were given
#"-" indicates STDIN
#"someprogram |" indicates the output of another program (7.14)
#read and parse input file
while (<>) {
	if ( $ARGV ne $input_file ) {
		print_log( sprintf( "< %s", File::Spec->rel2abs($ARGV) ) );
		$input_file = $ARGV;
	}
	chomp;
	if(!$_){
		next;
	}

	#comment
	if (/^\s*#/) { }

	#omim_results file
	elsif ( lc( $params{input_format} ) eq "omim_result" ) {
		if (/\s*[*+#%]?\d{6}.*;\s*(\S+)/) {
            print_log("$1");
            $output_record_count++;
		}
	}

	#NCBI gene_results file
	elsif ( lc( $params{input_format} ) eq "gene_result" ) {
		if (/\s*\d+:\s*(\S+)/) {
            print_log("$1");
            $output_record_count++;
		}
	}

	#Chromosome position
	elsif (/^\s*(?:chr)?[\s:]*([XYM]|[0-9]{1,2})[\s:]+([,\d]+)\D+([,\d]+)(?:\s+([01+-]{1}))?/i){
		if ( $1 && $2 && $3 ) {
			my $chr  = $1;
			my $start  = $2;
			my $end  = $3;
			my $strand = defined($4)
			 && $4 == "+" ? 1  : defined($4)
			 && $4 == "-" ? -1 : $4;
			$start =~ s/,//g;
			$end  =~ s/,//g;
			$start = $start <= $end ? $start : $end;
			$end  = $end >= $start ? $end  : $start;
			print_log( sprintf("%s:%d-%d%s",$chr, $start, $end, $strand==1?"+":"-" ));
            $output_record_count++;
		}
	}

	#Ensmbl gene stable id
	elsif (/^\s*(ENSG\d{11})/) {
		print_log($1);
		$output_record_count++;
	}

	#Ensmbl transcript stable id
	elsif (/^\s*(ENST\d{11})/) {
		print_log($1);
		$output_record_count++;
	}

	#Ensmbl exon stable id
	elsif (/^\s*(ENSE\d{11})/) {
		print_log($1);
		$output_record_count++;
	}

	#Gene external name
	elsif (/^\s*([^#\s]+).*#*/) {
		print_log($1);
		$output_record_count++;
	}
}

#print statistics
print_log(sprintf("%d record%s",$output_record_count,ss($output_record_count)));
exit;

sub now {
	use Time::localtime;
	my $tm = localtime;
	return sprintf "%04d-%02d-%02d %02d:%02d:%02d", $tm->year + 1900,
	 $tm->mon + 1, $tm->mday, $tm->hour, $tm->min, $tm->sec;
} # !now

sub ss {
	return $_[0] == 1 ? "" : "s";
} # !ss

sub print_log {
	print STDERR $_;
	print STDOUT $_ or die "Can't write STDOUT: $1";
}
