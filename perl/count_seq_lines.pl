#!/usr/bin/perl -w
use strict;
use warnings;

use constant SEQFILEREGEXL => 's_';
use constant SEQFILEREGEXR => '_\d+_seq.txt';
my $seqlineregex = '^[0-9-]+\s+[0-9-]+\s+[0-9-]+\s+[0-9-]+\s+([.ACGTacgt]+)\s*$';
my $usage = "perl count_seq_lines.pl parent_directory [comma_delimited_lane_numbers (default 1,2,3,4,5,6,7,8)]\n";
@ARGV or die $usage;
my $dir = shift(@ARGV);
-d $dir or die "$dir not found\n";
my $lanes = "1|2|3|4|5|6|7|8";
if (@ARGV && $ARGV[0] =~ /^\d+(,\d+)*\d*$/){$lanes = shift(@ARGV);}
print STDERR "$dir\nlanes $lanes\n";
my $seqfileregex = SEQFILEREGEXL.$lanes.SEQFILEREGEXR;
my $files = 0;
my $lines = 0;
my %good_read_lengths;
use File::Find;
find sub{
    if(/$seqfileregex/i){
        my $seqfile = $_;
        open(SEQFILE,"< $seqfile") or die "Can't open $seqfile $!\n";
        $files++;
        while(<SEQFILE>){
            if(/$seqlineregex/i){
                $lines++;
                my $last_read = index($1,'.');
                if ($last_read == -1){$last_read = length($1);}
                $good_read_lengths{$last_read}++;
            }
            else{print STDERR "unrecognized sequence in $seqfile\n\t$_\n";}
        }
    }
}, $dir;
print STDERR "$files files\n$lines sequence lines\n";
print STDERR "good reads\nfirst\tcount\n";
print STDERR map { "$_\t$good_read_lengths{$_}\n" } sort {$b<=>$a} keys %good_read_lengths;
