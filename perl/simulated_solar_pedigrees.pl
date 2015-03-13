#!/usr/bin/perl -w
use strict;
use File::Path;
use File::Spec;
# sibs 75% male; number of sibs in family 85% 2, 12% 3, 3% 4, 0.1% 5, 0.01% 6
my $dir = "/Volumes/raid/link/merlin/simulated_solar_pedigrees";
my @familysizes = (100,500,1000,5000,7000);
my $sibs;
for my $familysize(@familysizes){
    my $outfile = File::Spec->catfile($dir,"solar_$familysize.ped");
    open(OUT,">",$outfile);
    print OUT "famid,id,fa,mo,sex,mztwin,hhid\n";
    my $sibcount=0;
    for(my $i = 1;$i <= $familysize;$i++){
        my $rand = rand();
        if($rand<.85){
            $sibs = 302;
        }
        elsif($rand<.97){
            $sibs = 303;
        }
        elsif($rand<.999){
            $sibs = 304;
        }
        elsif($rand<.9999){
            $sibs = 305;
        }
        else{
            $sibs = 306;
        }
        printf(OUT "%u,201,0,0,1,,%u\n",$i,$i);
        printf(OUT "%u,202,0,0,2,,%u\n",$i,$i);
        for(my $j=301;$j<=$sibs;$j++){
            printf(OUT "%u,%u,201,202,%u,,%u\n",$i,$j,rand()<.75?1:2,$i);
            $sibcount++;
        }
    }
    print "$familysize families, $sibcount sibs\n";
}


