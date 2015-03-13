#!/usr/bin/perl -w
use strict;
use warnings;

use constant SEQFILEREGEXL => 's_';
use constant SEQFILEREGEXR => '_\d+_seq.txt';
my $seqlineregex = '^([0-9-]+)\s+([0-9-]+)\s+([0-9-]+)\s+([0-9-]+)\s+([.ACGTacgt]+)';
my $usage = "perl bustard2fasta.pl parent_directory [pipe_delimited_lane_numbers (default 1|2|3|4|5|6|7|8)]\n";
@ARGV or die $usage;
my $dir = shift(@ARGV);
-d $dir or die "$dir not found\n";
my $flowcell = "";
if (@ARGV){$flowcell = shift(@ARGV);}
my $lanes = "1|2|3|4|5|6|7|8";
if (@ARGV && $ARGV[0] =~ /^\d+(\|\d+)*$/){$lanes = shift(@ARGV);}
print STDERR "$dir\nlanes $lanes\n";
my $seqfileregex = SEQFILEREGEXL."(".$lanes.")".SEQFILEREGEXR;
my $files = 0;
my $lines = 0;
my $outfile = "";
use File::Find;
find sub{
    if(/$seqfileregex/i){
        my $lane = $1;
        my $seqfile = $_;
        my $new_outfile = "$dir\\$flowcell"."_lane_$lane.txt";
        if ($new_outfile ne $outfile){
            close(OUTFILE) or die "Can't close $outfile $!\n" unless !$outfile;
            $outfile = $new_outfile;
            open(OUTFILE,"> $outfile") or die "Can't open $outfile $!\n";
        }
        open(SEQFILE,"< $seqfile") or die "Can't open $seqfile $!\n";
        $files++;
        while(<SEQFILE>){
            if(/$seqlineregex/i){
                $lines++;
                printf(OUTFILE "%s\t%u\t%u\t%u\t%u\t%s\n",$flowcell,$1,$2,$3,$4,$5);
            }
            else{print STDERR "unrecognized sequence in $seqfile\n\t$_\n";}
        }
    }
}, $dir;
print STDERR "$files files\n$lines sequence lines\n";
