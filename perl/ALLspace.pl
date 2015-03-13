#!/usr/bin/perl -w
use strict;
use File::Glob;

my $pattern = 'D:\link\merlin_2010-05-19\merlin-00123456\*\*\merlin-00123456-*-*-regress-chrX.tbl';
my @glob = glob($pattern);
foreach my $file(@glob){
    open  IN, "<", $file;
    open OUT,">","$file.temp";
    while(<IN>){
        if(/^(.*)\tModel\s38\t(.*)$/g){
            print OUT "$1\treading\t$2\n";
        }
        else{
            print OUT;
        }
    }
    close IN;
    close OUT;
    unlink $file;
    rename "$file.temp", $file;
}

