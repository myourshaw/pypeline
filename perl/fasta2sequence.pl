#!/usr/bin/perl -w
use strict;
use warnings;
#BioPerl
use Bio::Seq;
use Bio::SeqIO;
#get sequences one per line from fasta files
my $usage = "perl fasta2sequence input fasta file list [> output_file]\n";
unless (@ARGV) {
	print $usage;
	exit;
}
my $input_file;
my $gt = "&gt;"; #komodo saves command line > as &gt;
my $output_file;
if(index($ARGV[$#ARGV],">") == 0){
	$output_file = substr pop(@ARGV),1;
}
elsif ($ARGV[$#ARGV-1] eq ">" || $ARGV[$#ARGV-1] eq $gt) {
	$output_file = pop(@ARGV);
	pop(@ARGV);
}
elsif (rindex($ARGV[$#ARGV-1],">") == length($ARGV[$#ARGV-1])-1) {
	$output_file = pop(@ARGV);
	$ARGV[$#ARGV-1] = substr($ARGV[$#ARGV-1],0,rindex($ARGV[$#ARGV-1],">"));
}
elsif (rindex($ARGV[$#ARGV-1],$gt) == length($ARGV[$#ARGV-1])-1) {
	$output_file = pop(@ARGV);
	$ARGV[$#ARGV-1] = substr($ARGV[$#ARGV-1],0,rindex($ARGV[$#ARGV-1],$gt));
}
if($output_file){
	open STDOUT, ">", $output_file or die "can't open $output_file: $!";
}

while ($input_file = shift @ARGV) {
    my $sequence_file = Bio::SeqIO->new(-file => $input_file, -format => 'fasta' );
		while (my $sequence = $sequence_file->next_seq)
		{
            print STDERR $sequence->seq."\n";
			print $sequence->seq."\n";
		}
}
