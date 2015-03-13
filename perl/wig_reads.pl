#!/usr/bin/perl -w
use strict;
use warnings;
use POSIX qw(ceil floor);
my $bin_width = 1;
my $bin_max = 20;
# [bool include,rangeStart,rangeEnd]
my %chroms = (
    chr1 => [0],
    chr2 => [0],
    chr3 => [0],
    chr4 => [0],
    chr5 => [0],
    chr6 => [0],
    chr7 => [0],
    chr8 => [0],
    chr9 => [0],
    chr10 => [0],
    chr11 => [0],
    chr12 => [0],
    chr13 => [0],
    chr14 => [0],
    chr15 => [0],
    chr16 => [0],
    chr17 => [1,20739925,46995695],
    chr18 => [0],
    chr19 => [0],
    chr20 => [0],
    chr21 => [0],
    chr22 => [0],
    chrM => [0],
    chrX => [0],
    chrY => [0],
    );

my $wigregex = '^(chr\S+)\s+([0-9-]+)\s+([0-9-]+)\s+([0-9-]+)';
my $usage = "perl wig_reads.pl input_file\n";
@ARGV or die $usage;
#my $input_file = shift(@ARGV);
#-f $input_file or die "$input_file not found\n";
my %chr_positions = ();
my %chr_reads = ();
my %chr_histogram = ();
while(<>){
    chomp;
    if(/$wigregex/i){
        my ($chr,$chromStart,$chromEnd,$reads)=($1,$2,$3,$4);
         if ($chroms{$chr}[0]){
            if(!exists($chroms{$chr}[1]) || $chromEnd>=$chroms{$chr}[1]){
                if(!exists($chroms{$chr}[2]) || $chromEnd<=$chroms{$chr}[2]){
                    $chr_positions{$chr}++;
                    $chr_reads{$chr}+=$reads;
                    if($reads > $bin_max){
                        $chr_histogram{"$chr\t>=\t".($bin_max+1)}++;
                    }
                    else{
                        $chr_histogram{"$chr\t<=\t".(ceil(($reads+1)/$bin_width)-1)}++;
                    }
                }
            }
        }
    }
}
print "positions with reads\n";
print map { "$_\t$chr_positions{$_}\n" } sort {$a cmp $b} keys %chr_positions;
print "reads\n";
print map { "$_\t$chr_reads{$_}\n" } sort {$a cmp $b} keys %chr_reads;
print "histogram\n";
print map { "$_\t$chr_histogram{$_}\n" } sort {$a cmp $b} keys %chr_histogram;
