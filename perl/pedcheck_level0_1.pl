#!/usr/bin/perl -w
use strict;
my $ dir = "/Volumes/raid/link/pedcheck_2010-01-06";
my @errfiles = qw(autosome-pedcheck-all.err  x-pedcheck-all.err);

my $level0file = "$dir/pedcheck_level0.err";
open L0,">",$level0file or die "can't open $level0file";
print L0 "fam\tind\terr\n";

my $level1file = "$dir/pedcheck_level1.err";
open L1,">",$level1file or die "can't open $level1file";
print L1 "fam\tmarker\terr\n";

foreach my $errfile(@errfiles){
    my $level = -1;
    my $pedigree = -1;
    my $marker = "";
    open ERR,"<","$dir/$errfile" or die "can't open $errfile";
    while(<ERR>){
        if(/\*+\s+LEVEL\s+(\d)\s+ERRORS\s+\*/){
            $level = $1;
        }
        elsif($level == 0){
            if(/ERROR:\s+Ped:\s+(\d+)\s+Per:\s+(\d+)\s+-\s+(.+)/){
                print L0 "$1\t$2\t$3\n"
            }
            elsif(/ERROR/i){
                print STDERR "$_\n";
            }
        }
        elsif($level == 1){
            if(/#+\s+GENOTYPE\s+ERROR:\s+Pedigree\s+(\d+)\s+Locus\s+\d+\s+Name\s+(\S+)\s+#+/){
                $pedigree = $1;
                $marker = $2;
            }
            elsif(/ERROR:\s+(.+)/){
                print L1 "$pedigree\t$marker\t$1\n"
            }
            elsif(/ERROR/i){
                print STDERR "$_\n";
            }
        }
        elsif($level == 2){
             ##### GENOTYPE ERROR: Pedigree 721   Locus  3098    Name D8S1784 #####
            if(/#+\s+GENOTYPE\s+ERROR:\s+Pedigree\s+(\d+)\s+Locus\s+\d+\s+Name\s+(\S+)\s+#+/){
                $pedigree = $1;
                $marker = $2;
                print L1 "$pedigree\t$marker\tlevel 2 error\n"
            }
            elsif(/ERROR/i){
                print STDERR "$_\n";
            }
        }
        elsif($level == 3){
             if(/#+\s+GENOTYPE\s+ERROR:\s+Pedigree\s+(\d+)\s+Locus\s+\d+\s+Name\s+(\S+)\s+#+/){
                $pedigree = $1;
                $marker = $2;
                print L1 "$pedigree\t$marker\tlevel 3 error\n"
            }
            elsif(/ERROR/i){
                print STDERR "$_\n";
            }
       }
        elsif(/ERROR/i){
            print STDERR "$_\n";
        }
    }
}
