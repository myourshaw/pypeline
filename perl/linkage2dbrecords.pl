#!/usr/bin/perl -w
use strict;
use File::Path;
use File::Spec;

my @chrs = qw(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X XY Y);
my @files = qw(chr chr chr chr fine_chr finemap_chr finemap_chr chr finemap_chr chr finemap_chr chr finemap_chr chr finemap_chr chr finemap_chr chr chr finemap_chr chr chr chr chr chr);
my $out = "/Volumes/raid/link/genotypes.txt";
open OUT,">",$out or die "can't open $out";
print OUT "chr\tfam\tind\tfid\tmid\tsex\tadhd\tmarker\ta\tb\n";
for(my $c = 0;$c<=24;$c++){
    my $chr = $chrs[$c];
    my $file=$files[$c];
    my $folder = ($chr eq "XY" || $chr eq "Y") ? "" : "/wave1_2/new_pheno";
    my $dir = "/Volumes/raid/data/ADHD_data/ADHD\ -\ projects\ 2006-2008/PROJECTS/Wave1_2_3_4_files_for_analyses/chr$chr$folder";
    my $dat = File::Spec->catfile($dir,"$file$chr.dat");
    my $ped = File::Spec->catfile($dir,"$file$chr.ped");
    my @markers;
    open DAT,"<",$dat or die "can't open $dat";
    while(<DAT>){
        if(/M\s+(\S+)/){
            push @markers,$1;
        }
    }
    close DAT;
    open PED,"<",$ped or die "can't open $ped";
    my $markercount = scalar @markers;
    while(<PED>){
        if(/\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)/){
            my @line = split /\s+/, $_;
            for(my $i = 0;$i<@markers;$i++){
                my($fam,$ind,$fid,$mid,$sex,$adhd,$marker,$a,$b)=($1,$2,$3,$4,$5,$6,$markers[$i],$line[$i*2+6],$line[$i*2+7]);
                if($a eq "Y"){$a=$b;}
                print OUT "$chr\t$fam\t$ind\t$fid\t$mid\t$sex\t$adhd\t$marker\t$a\t$b\n";
            }
        }
    }
    @markers = ();
}


