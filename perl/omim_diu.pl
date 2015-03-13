#!/usr/bin/perl -w
use strict;
use warnings;
#compares any number of omim_result files and prints difference, intersection and/or union files based on MIM numbers
#last record read wins
my $usage = "usage: perl omim_diu.pl [-d difference_file] [-i intersection_file] [-u union_file] input file list\n";
unless (@ARGV){
    print $usage;
    exit;
}
use Getopt::Std;
my %files;
getopt('diu', \%files); #output files: difference, intersection, union
my $omim_file_count = @ARGV;
my (%records,%count);
my ($MIM_number,$record);
while(<>){
	if (/\s*[*+#%]?(\d{6})\s+\S+/) {
        if($MIM_number){
            $records{$MIM_number} = $record;
			$count{$MIM_number}++;
        }
        $MIM_number = $1;
        $record = $_;
    }
    else{
        if($MIM_number){
            $record .= $_;
        }
    }
}
if($MIM_number){
	$records{$MIM_number} = $record;
	$count{$MIM_number}++;
}
#difference (records that are in only one file)
if($files{d}){
	open DIFFERENCE, ">", $files{d} or die "can't open $files{d}: $!";
	while((my $key,my $value) = each(%records)){
		print STDERR $value unless $count{$key}!=1;
		print DIFFERENCE $value unless $count{$key}!=1;
	}
}
#intersection (records that are in all files)
if($files{i}){
	open INTERSECTION, ">", $files{i} or die "can't open $files{i}: $!";
	while((my $key,my $value) = each(%records)){
		print INTERSECTION $value unless $count{$key}<$omim_file_count;
	}
}
#union (records that are in any file)
if($files{u}){
	open UNION, ">", $files{u} or die "can't open $files{u}: $!";
	while((my $key,my $value) = each(%records)){
		print UNION $value;
	}
}
