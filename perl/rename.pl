#!perl -w
use strict;
use warnings;
use File::Glob;

#Cookbook recipe 9-5
# 's/non-white/non_white/' ~/Lab/ADHD/merlin_files/npls/*non-white*

my $op = shift or die "Usage: rename expr [files]\n";
chomp(@ARGV = <STDIN>) unless @ARGV;

for(@ARGV){
    for(glob(shift)){
        my $was = $_;
        eval $op;
        die $@ if $@;
        rename($was,$_) unless $was eq $_;
    }
}


