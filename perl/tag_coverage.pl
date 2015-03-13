#!/usr/bin/perl -w
use strict;
use warnings;
use LWP::Simple;
my $usage = "tag_coverage usage\n"
	 . "perl tag_coverage.pl prints usage.\n"
	 . "perl tag_coverage.pl index_url [output_file] lane_list\n"
     . "prints percent and count for each 3-base tag at 3' (ignoring final 3' T)\n"
     . "in a Solexa seq file with a name in the form s_<lane>_<tile>_seq.txt\n"
     . 'e.g., perl tag_coverage.pl "http://genomics.ctrl.ucla.edu/~solexa/solexa_datasets/080108_HWI-EAS172_200G8/Data/C1-36_Firecrest1.8.28_25-01-2008_solexa/Bustard1.8.28_25-01-2008_solexa/" "C:\My Data\tags.txt" 1 2 3'."\n";
unless (@ARGV) {
	print $usage;
	exit;
}
#get url of directory from command line
my $dir = shift(@ARGV);
$dir .= "/" unless rindex($dir,"/") == length($dir)-1;
print "$dir\n";
#get optional output file
if($ARGV[0] !~ /^\d+$/){
    my $outfile = shift(@ARGV);
    open STDOUT, ">", $outfile or die "can't open $outfile: $!";
}
#get directory index
$_ = get($dir) or die $!;
#get file names
my @files = /<a href="(s_\d+_\d+_seq\.txt)">/gi;
#get lanes from command line
while(defined(my $lane = shift(@ARGV))){
    my %tags;
    my $tag_total;
    #get sequence files
    foreach my $file(@files){
        if($file=~/s_($lane)_\d+_seq\.txt/){
            #read file
            $_ = get($dir.$file) or die $!;
            while(/([ACGT]{3})T\s/gi){
                $tags{$1}++;
                $tag_total++;
            }
        }
    }
    if (%tags){
        printf "lane $lane\n";
        print map { sprintf("%s\t%8.4f%%\t%8d\n",$_,100*$tags{$_}/$tag_total,$tags{$_}) } sort {$tags{$b} <=> $tags{$a} || $a cmp $b} (keys %tags);
        print "total\t$tag_total\n";
    }
}
