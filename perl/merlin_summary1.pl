#!/usr/bin/perl -w
use strict;
use warnings;

use File::Spec;
my $threshhold = 3;

my $dir = "/Volumes/raid/merlin_files/m-10-20";
my @ethnics = qw(all_ethnic non-white white);
my @genders = qw(all_gender female_family female_sibs male_family male_sibs);
#my @genders = qw(all_gender female_sibs male_sibs);
my @adhd_lifes = qw(adhd_broad_life adhd_narrow_life);
my @adhd_types = qw(all_type combined hyperactive inattentive);
#my @chromosomes = (1 .. 22,"X","XY","Y");
my @chromosomes = (1 .. 22);

my $npl_out_file = File::Spec->catfile($dir,"merlin-npl.tbl");
open OUT,">",$npl_out_file;
my $npl_summary_file = File::Spec->catfile($dir,"merlin-npl-summary.tbl");
open SUM,">",$npl_summary_file;
print "THRESHHOLD\tETHNIC\tGENDER\tADHD_LIFE\tADHD_TYPE\tCHR\tPOS\tMARKER\tLABEL\tANALYSIS\tZSCORE\tDELTA\tLOD\tPVALUE\tExDELTA\tExLOD\tExPVALUE\n";
print SUM "THRESHHOLD\tETHNIC\tGENDER\tADHD_LIFE\tADHD_TYPE\tCHR\tPOS\tMARKER\tLABEL\tANALYSIS\tZSCORE\tDELTA\tLOD\tPVALUE\tExDELTA\tExLOD\tExPVALUE\n";
print OUT "ETHNIC\tGENDER\tADHD_LIFE\tADHD_TYPE\tCHR\tPOS\tMARKER\tLABEL\tANALYSIS\tZSCORE\tDELTA\tLOD\tPVALUE\tExDELTA\tExLOD\tExPVALUE\n";
foreach my $ethnic(@ethnics){
    foreach my $gender(@genders){
        foreach my $adhd_life(@adhd_lifes){
            foreach my $adhd_type(@adhd_types){
                my $data = "$ethnic\t$gender\t$adhd_life\t$adhd_type";
                foreach my $chr(@chromosomes){
                    my $merlin_npl_in_file = File::Spec->catfile($dir,$ethnic,$gender,$adhd_life,$adhd_type,"chr$chr","merlin-npl-chr$chr.tbl");
                    if (!-f $merlin_npl_in_file){
                        print STDERR "\n*** mising $merlin_npl_in_file\n\n";
                    }
                    else{
                        open IN,"<",$merlin_npl_in_file;
                        my ($chr2,$pos,$marker,$label,$analysis,$zscore,$delta,$lod,$pvalue,$exdelta,$exlod,$expvalue);
                        while(<IN>){
                            if(/\s*([0-9XY]{1,2})\s+([e\d.-]+)\s+(\S+)\s+(\S+)\s+\[(\S+)\]\s+([e\d.-]+)\s+([e\d.-]+)\s+([e\d.-]+)\s+([e\d.-]+)\s+([e\d.-]+)\s+([e\d.-]+)\s+([e\d.-]+)/){
                                ($chr,$pos,$marker,$label,$analysis,$zscore,$delta,$lod,$pvalue,$exdelta,$exlod,$expvalue) = ($1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12);
                                printf OUT "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $data,$chr,$pos,$marker,$label,$analysis,$zscore,$delta,$lod,$pvalue,$exdelta,$exlod,$expvalue;
                                if(($lod >= $threshhold && $lod != 9.999) || ($exlod >= $threshhold && $exlod != 9.999)){
                                    printf SUM ">=%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $threshhold,$data,$chr,$pos,$marker,$label,$analysis,$zscore,$delta,$lod,$pvalue,$exdelta,$exlod,$expvalue;
                                    printf ">=%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $threshhold,$data,$chr,$pos,$marker,$label,$analysis,$zscore,$delta,$lod,$pvalue,$exdelta,$exlod,$expvalue;                                    
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
print "done\n";
exit;




