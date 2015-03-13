#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;
use File::Spec;
use Date::Calc qw(Today_and_Now Delta_YMDHMS);
use POSIX qw(ceil floor);
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
#BioPerl
use Bio::Seq;
use Bio::SeqIO;
#EnsEMBL Perl API
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::Utils::Sequence qw (sequence_with_ambiguity);
#options
my %opts; #sequence, genomic intervals, log
getopts("s:Sg:Gl:L", \%opts);
#files
unless (@ARGV) {
	print_usage();
	exit;
}
my $input_file = "";
my $input_file1 = $ARGV[0];
my ($vol,$dir,$file);
if ( $input_file1 eq "-" ) {
	print STDERR "Input line-by line. Ctrl-D when done\n";
	($vol,$dir,$file) = (File::Spec->curdir()."wig_output");
}
else{
	($vol,$dir,$file) = File::Spec->splitpath($input_file1);
}
#redirect STDOUT
my $output_file;
$output_file = $opts{l} or $output_file = File::Spec->catpath($vol,$dir,"$file.density.txt");
open STDOUT, ">", $output_file or die "can't open $output_file: $!";
my $reads = 0;
my $pos = 0;
print "chr\tpos\treads\n";
while (<>) {
    #chr 0basepos 1basepos reads
    if(/chr([0-9MTXY]{1,2})\s+\d+\s+(\d+)\s+(\d+)/i){
            printf "%s\t%u\t%u\n",$1,$2,$3;
    }
 }

sub print_usage{
    use constant USAGE => "nothing to be found here\n";
    print STDERR USAGE;
}
