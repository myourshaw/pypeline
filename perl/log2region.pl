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
getopts("o:i:", \%opts);
my $id = "";
if ($opts{i}){$id = $opts{i};}
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
	($vol,$dir,$file) = (File::Spec->curdir()."region_output");
}
else{
	($vol,$dir,$file) = File::Spec->splitpath($input_file1);
}
#redirect STDOUT
my $output_file;
$output_file = $opts{o} or $output_file = File::Spec->catpath($vol,$dir,"$file.regions.txt");
open STDOUT, ">", $output_file or die "can't open $output_file: $!";
print "id\tgeneSymbol\troi\tchrom\troiStart\troiEnd\tstrand\n";
while (<>) {
    #chr 0basepos 1basepos reads
    if(/(\S+) \|(chr[0-9MTXY]{1,2}):(\d+)-(\d+)\[([+-])\]\|\d+b\|(\S+)\|ENS/i){
        printf"%s\t%s\t%s\t%s\t%u\t%u\t%s\n",$id,$6,substr($1,length($6)+1),$2,$3,$4,$5;
    }
 }


sub print_usage{
    use constant USAGE => "perl log2region.pl [-i id] input_file\n";
    print STDERR USAGE;
}
