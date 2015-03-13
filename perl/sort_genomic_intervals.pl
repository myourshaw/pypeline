#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;
#sorts genomic interval file by chromosome (1-22,X,Y,M), start, end
my $usage = "usage: perl [-s] sort_genomic_intervals.pl input file list [> output_file]\n-s = include strand in output\n";
my %opts; #sequence, genomic intervals, log
getopts("s", \%opts); # print strand
unless (@ARGV){
    print $usage;
    exit;
}
#redirect STDOUT
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
my @intervals;
while (<>) {
	chomp;
	if (/^\s*(?:chr)?[\s]*([0-9XYMT]{1,2})[\s:]+([,\d]+)\D+([,\d]+)(?:\s+([01+-]{1}))?/i){
        if ($1 && $2 && $3) {
            my $chr = $1;
			my $num_chr = chr2numchr($chr);
            my $start = $2;
            my $end = $3;
            my $strand = defined($4)
             && $4 == "+" ? 1 : defined($4)
             && $4 == "-" ? -1 : $4;
            $start =~ s/,//g;
            $end =~ s/,//g;
            $start = $start <= $end ? $start : $end;
            $end = $end >= $start ? $end : $start;
		push @intervals, {
            chr=>$chr,
            num_chr=>$num_chr,
            start=>$start,
            end=>$end,
			strand=>$strand};
        }
    }
}

if(@intervals){
	@intervals = sort {$a->{num_chr} <=> $b->{num_chr} || $a->{start} <=> $b->{start} || $a->{end} <=> $b->{end} || $a->{strand} cmp $b->{strand}} @intervals;
	foreach my $region(@intervals){
		printf "chr%s:%d-%d%s\n",$region->{chr},$region->{start},$region->{end},$opts{s}?" $region->{end}":"";
	}
}

sub chr2numchr{
    my $chr = uc(shift);
    return $chr eq "X" ? 23 : $chr eq "Y" ? 24 : ($chr eq "M" || $chr eq "MT") ? 25 : $chr;
}

sub numchr2chr{
    my $chr = shift;
    return $chr == 23 ? "X" : $chr == 24 ? "Y" : $chr == 25 ? "M" : $chr;
}
