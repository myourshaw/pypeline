#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;
my %opts; # -a convert correlations to absolute values
getopts("a", \%opts);

unless (@ARGV) {
	print "usage: perl clean_gene-gene_corr.pl [-a] input_file output_file/n-a convert correlations to absolute values\n";
	exit;
}

#output file name is last input argument
open OUT, ">", pop;

my @gene_column_names;
my @probe_column_names;
my @column_labels;
my $column_labels;
my %correlations;

#read input file line-by-line
while(<>){
	chomp;
	#gene header row
	if(/\t\tGene\.Symbol\t\t\t(.+)/){
		my $genes = $1;
		#convert " /// " to ";" and spaces to "_"
		$genes =~ s/ \/\/\/ /;/g;
		$genes =~ s/ /_/g;
		#put gene column names in an array
		@gene_column_names = split(/\t/, $genes);
	}
	#probe header row
	elsif(/\t\t\tProbe\.Set\.ID\tMean\t(.+)/){
		my $probes = $1;
		$probes =~ s/ /_/g;
		#put probe column names in an array
		@probe_column_names = split(/\t/, $probes);
	}
	#correlation data rows
	elsif(/^.+\t.+\t(.+)\t(\d+.*_at)\t.+?\t(.+)/){
		my ($gene, $probe, $corr) = ($1,$2,$3);
		$gene =~ s/ \/\/\/ /;/g;
		$gene =~ s/ /_/g;
		$probe =~ s/ /_/g;
		#if command line had -a option, convert correlations to absolute values
		if($opts{a}){
			$corr =~ s/-//;
		}
		#put correlations in a hash with "gene|probe" as key
		$correlations{"$gene|$probe"} = $corr;
	}
}
#column labels are "gene|probe"
@column_labels = map{"$_|".shift(@probe_column_names)} @gene_column_names;
$column_labels = join("\t",@column_labels);
print OUT "$column_labels\n";
#rows are "gene|probe	corr.1	corr.2	...	corr.n"
print OUT map{"$_\t$correlations{$_}\n"} sort keys %correlations;
