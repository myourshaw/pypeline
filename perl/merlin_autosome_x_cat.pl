#!/usr/bin/perl -w
use strict;
use warnings;

use File::Path;
use File::Spec;
use File::Glob;
use Switch;

my %codes = (
    phase1  => "1", #phase1 markers, phase1 sample and part phase2 sample
    phase2  => "2", #phase2 markers, part phase1 and all phase2 sample
    "phase1+2"  => "3", #phase1 and phase2 markers, phase1 and phase2 sample (combined)
    "phase2-1" => "4", #phase2 markers, no phase1 and all phase2 sample
    freq => "f",
    ibd => "i",
    lod => "l",
    out => "o",
    pdf => "p",
    table => "t",
    zscore => "z",
    all_ethnic  => "a",
    non_white  => "n",
    white  => "w",
    all_gender  => "a",
    female_family  => "f",
    female_sibs  => "g",
    male_family  => "m",
    male_sibs  => "n",
    broad  => "b",
    narrow  => "n",
    all_affected_family  => "b",
    all_type  => "a",
    combined  => "c",
    combined_family  => "d",
    combined_inattentive  => "e",
    combined_inattentive_family  => "f",
    hyperactive  => "h",
    hyperactive_family  => "k",
    inattentive  => "i",
    inattentive_family  => "j",
    ALL => "a",
    Pairs => "p",
    QTL => "q",
    regress => "r"
);

#input files from merlin: /Users/myourshaw/lab/ADHD/Paper/link/merlin_2011-03-15/merlin-00123456/{'00000'..'01000'}/{autosome,x}/{*-nonparametric.tbl,*-regress-chr*.tbl}
#output files: /Users/myourshaw/lab/ADHD/Paper/link/merlin_2011-03-15//merlin-00123456
my $dir = shift;
my $header = "sim\tmarker\tanalysis\tvariable\tlod\tpvalue\n";
open OUT,">",File::Spec->catfile($dir,"cat-sim-npl_qtl_regress.txt") or die "can't open output file $!\n";
print OUT $header;
#open NPL,">",File::Spec->catfile($dir,"cat-npl_qtl.txt") or die "can't open output file $!\n";
#print NPL $header;
#open REG,">",File::Spec->catfile($dir,"cat-regress.txt") or die "can't open output file $!\n";
#print REG $header;
for ('00001'..'01000'){
    my $sim=$_;
    my $glob = "$dir/$sim/*/*-nonparametric.tbl";
    for my $file(glob($glob)){
        open IN, "<", $file or die "can't open $file $!\n";
        while(<IN>){
            if(/\s*([0-9XY]{1,3})\s+([Ee\d.-]+)\s+(\S+)\s+(\S+)\s+\[(\S+)\]\s+([Ee\d.-]+)\s+([Ee\d.-]+)\s+([Ee\d.-]+)\s+([Ee\d.-]+)/){
                my ($chr,$pos,$marker,$variable,$analysis,$zscore,$delta,$lod,$pvalue) = ($1,$2,$3,$4,$5,$6,$7,$8,$9);
                print OUT "$sim\t$marker\t$codes{$analysis}\t$variable\t$lod\t$pvalue\n";
            }
            elsif(!/^(na)|(CHR)/){
                print;
            }
        }
    }
    $glob = "$dir/$sim/*/*-regress-*.tbl";
    for my $file(glob($glob)){
        open IN, "<", $file or die "can't open $file $!\n";
        while(<IN>){
            if(/\s*([0-9XY]{1,3})\s+(\S+)\s+(\S+)\s+([Eena\d.-]+)\s+([Eena\d.-]+)\s+([Eena\d.-]+)\s+([Eena\d.-]+)\s+([Eena\d.-]+)/){
                my ($chr,$marker,$variable,$h2,$sd,$info,$lod,$pvalue) = ($1,$2,$3,$4,$5,$6,$7,$8);
                $h2 = $h2 eq "na" ? 0 : $h2;
                $sd = $sd eq "na" ? 0 : $sd;
                $info = $info eq "na" ? 0 : $info;
                $lod = $lod eq "na" ? 0 : $lod;
                $pvalue = $pvalue eq "na" ? 1 : $pvalue;
                print OUT "$sim\t$marker\t$codes{regress}\t$variable\t$lod\t$pvalue\n";
            }
            elsif(!/^(na)|(CHR)/){
                print "$file\n";
                print;
            }
        }
    }
}



