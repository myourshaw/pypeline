#!/usr/bin/perl -w
use strict;
#my $wigfile = 'C:\Users\Michael\solexa\solexa_datasets\080620_HWI-EAS172_20E07\Data\C1-36_Firecrest1.9.2_24-06-2008_solexa\Bustard1.9.2_24-06-2008_solexa\pixel_read_counts.txt';
#open (WIG, "< $wigfile");
my $wig_regex = '(chr[0-9XYMT]{1,2}(?:_random)?)\s(\d+)\s(\d+)\s(\d+)';
my $positions = 0;
my $reads = 0;
while(<>){
    if (/$wig_regex/){
        $positions++;
        $reads+=$4;
    }
    my($lane,$tile,$tileCol,$tileRow,$tileX,$tileY,$laneX,$laneY,$As,$Cs,$Gs,$Ts,$Ns) = split /\s/;
    if($As == 36){
        print "$laneX,$laneY,$As\n";
    }
}
