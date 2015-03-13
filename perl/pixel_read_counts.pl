#!/usr/bin/perl -w
use strict;
use warnings;
use POSIX qw(ceil floor);

#camera parameters
use constant TILECOLS => 3;
use constant TILEROWS => 110;
use constant TILEXPIXELS => 1000; #old camera is 1000, new camera is 2000 (there can be slightly out of range readings)
use constant TILEYPIXELS => 1000; #old camera is 1000, new camera is 2000 (there can be slightly out of range readings)

use constant SEQFILEREGEXL => 's_';
use constant SEQFILEREGEXR => '_(\d+)_seq.txt';
my $seqlineregex = '^(\d+)\s+(\d+)\s+-?(\d+)\s+-?(\d+)\s+([.ACGTacgt]+).*$';
my $usage = "perl pixel_read_counts.pl parent_directory [comma_delimited_lane_numbers (default 1,2,3,4,5,6,7,8)]\n";

@ARGV or die $usage;
my $dir = shift(@ARGV);
-d $dir or die "$dir not found\n";
open(OUT,"> $dir\\pixel_read_counts.txt") or die "Can't open $dir \\A_bad.txt $!\n";
print OUT "lane\ttile\ttileCol\ttileRow\ttileX\ttileY\tlaneX\tlaneY\tAs\tCs\tGs\tTs\tNs\n";
my $lanes = "1|2|3|4|5|6|7|8";
if (@ARGV && $ARGV[0] =~ /^\d+(,\d+)*\d*$/){$lanes = shift(@ARGV);}
print STDERR "$dir\nlanes $lanes\n";
my $seqfileregex = SEQFILEREGEXL.$lanes.SEQFILEREGEXR;
my $files = 0;
my $lines = 0;
use File::Find;
find sub{
    if(/$seqfileregex/i){
        my $lane = $1;
        my $tile = $2;
        my $seqfile = $_;
        open(SEQFILE,"< $seqfile") or die "Can't open $seqfile $!\n";
        $files++;
        while(<SEQFILE>){
            chomp;
            if(/$seqlineregex/i){
                $lines++;
                my $input_line = $_;
                my ($lane,$tile,$tileX,$tileY,$seq) = ($1,$2,$3,$4,$5);
                my $As = ($seq =~ tr/A//);
                my $Cs = ($seq =~ tr/C//);
                my $Gs = ($seq =~ tr/G//);
                my $Ts = ($seq =~ tr/T//);
                my $Ns = ($seq =~ tr/\.//);
                #http://www.genographia.org/portal/topics/example-data-sets/illumina-solexa-examples/flow-cell-wall
                my ($tileCol,$tileRow);
                if ($lane%2){
                    #odd lane
                    $tileCol = ceil($tile/TILEROWS);
                }
                else{
                    #even lane
                    $tileCol = TILECOLS - ceilceil($tile/TILEROWS);
                }
                if($lane%2 + $tileCol%2 == 1){
                    #odd lane + even tileCol or even lane + odd tileCol
                    $tileRow = TILEROWS + 1 - ($tile - TILEROWS*($tileCol - 1));
                }
                else{
                   #odd lane + odd tileCol or even lane + even tileCol
                    $tileRow = $tile - TILEROWS*($tileCol - 1);
                }
                my $laneX = $tileX + TILEXPIXELS * ($tileCol - 1);
                my $laneY = $tileY + TILEYPIXELS * ($tileRow - 1);
                printf OUT "%u\t%u\t%d\t%d\t%u\t%u\t%d\t%d\t%u\t%u\t%u\t%u\t%u\t%s\n",$lane,$tile,$tileX,$tileY,$tileCol,$tileRow,$laneX,$laneY,$As,$Cs,$Gs,$Ts,$Ns,$seq;
            }
            else{print STDERR "unrecognized sequence in $seqfile\n\t$_\n";}
        }
    }
}, $dir;
print STDERR "$files files, $lines lines\n";
