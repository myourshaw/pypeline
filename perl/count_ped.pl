#!perl -w
use strict;
use warnings;

use File::Path;
use File::Spec;

#count ped files

my $dir = "/Volumes/raid/link/merlin_2009-11-07a";
my $peds = 0;
my $nopeds = 0;
my @phases = qw(phase1 phase2 phase1+2 phase2-1);
my @ethnics = qw(all_ethnic non_white white);
my @genders = qw(all_gender female_family female_sibs male_family male_sibs);
my @adhds = qw(broad narrow);
my @types = qw(all_type combined combined_family combined_inattentive combined_inattentive_family hyperactive hyperactive_family inattentive inattentive_family);
my @chromosomes = ("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22","X","XY","Y");

    foreach my $phase(@phases){
        foreach my $ethnic(@ethnics){
            foreach my $gender(@genders){
                foreach my $adhd(@adhds){
                    foreach my $type(@types){
                        foreach my $chr(@chromosomes){
                            my $ped_file = File::Spec->catfile($dir,$ethnic,$gender,$adhd,$type,"$phase-chr$chr","$phase-chr$chr.ped");
                            if (-f $ped_file){
                                $peds++;
                            }
                             else{
                                $nopeds++;
                            }
                        }
                    }
                }
            }
        }
    }



print "peds: $peds\nmissing: $nopeds/n";
exit;
