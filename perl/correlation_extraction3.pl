#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Std;

my %opts;
getopts("astv", \%opts);

unless (@ARGV>1) {
     print STDERR "usage: perl clean_gene-gene_corr.pl [-[a][s][t][v]] input_file output_file\n"
	 ."-a convert correlations to absolute values\n"
	 ."-s sort rows by gene|probe (slow)\n"
	 ."-t insert tab in column 1 row 1 (default is to omit tab for R heatmap format)\n"
	 ."-v verbose output\n";
     exit;
}

#output file name is last input argument
open OUT, ">", pop;

my @gene_column_names;
my @probe_column_names;
my @column_labels;
my $column_labels;
my %correlations;
my $header_printed;;

#read input file line-by-line
while(<>){
     chomp;
	 unless($header_printed){
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
		  if(@probe_column_names && @gene_column_names){
			   #column labels are "gene|probe"
			   @column_labels = map{"$_|".shift(@probe_column_names)} @gene_column_names;
			   $column_labels = join("\t",@column_labels);
			   my $line = ($opts{t}?"\t":"")."$column_labels\n";
			   print STDERR $line if $opts{v};
			   print OUT $line;
			   $header_printed = 1;
		  }
	 }
     #correlation data rowsprint OUT "$column_labels\n";
     elsif(/^.+\t.+\t(.+)\t(\d+.*_at)\t.+?\t(.+)/){
		  my ($gene, $probe, $corr) = ($1,$2,$3);
		  $gene =~ s/ \/\/\/ /;/g;
		  $gene =~ s/ /_/g;
		  $probe =~ s/ /_/g;
		  #if command line had -a option, convert correlations to absolute values
		  if($opts{a}){
		  $corr =~ s/-//;
		  }
		  unless($opts{s}){
			   my $line = sprintf("%s|%s\t%s\n",$gene,$probe,$corr);
			   print STDERR $line if $opts{v};
			   print OUT $line;
		  }
		  else{
			   #put correlations in a hash with "gene|probe" as key
			   $correlations{"$gene|$probe"} = $corr;
			   print STDERR "." if $opts{v};
		  }
     }
}
if($opts{s}){
	 #rows are "gene|probe corr.1     corr.2     ...  corr.n"
	 print STDERR "\n" if $opts{v};
	 print STDERR map{"$_\t$correlations{$_}\n"} sort keys %correlations if $opts{v};
	 print OUT map{"$_\t$correlations{$_}\n"} sort keys %correlations;
}