#!/usr/bin/perl -w
use strict;
use warnings;

use File::Path;
use File::Spec;
my $dir = "/Volumes/raid/mega2_files/temp";
#my @merlins = ([qw(merlin m)], [qw(pedwiped p)]);
my @merlins = ([qw(pedwiped p)]);
my @studies = ([qw(all a)], [qw(cidr c)], [qw(phase1 p)]);
my @ethnics = ([qw(all_ethnic a)], [qw(non_white n)], [qw(white w)]);
my @genders = ([qw(all_gender a)], [qw(female_family f)], [qw(female_sibs g)], [qw(male_family m)], [qw(male_sibs n)]);
my @adhd_lifes = ([qw(adhd_broad_life b)], [qw(adhd_narrow_life n)]);
my @adhd_types = ([qw(all_type a)], [qw(combined c)], [qw(hyperactive h)], [qw(inattentive i)]);
my %analyses = (ALL => "a", Pairs => "p", QTL => "q");
my @chromosomes = ("01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20","21","22","X","XY","Y");
#my @chromosomes = (1 .. 22);

my $err_file = File::Spec->catfile($dir,"merlin.err");
open ERR,">",$err_file or die "can't open $err_file ($!)\n";
my $err_header = "fam\tid\tmarker\tratio";
print ERR "$err_header\n";
foreach my $chr(@chromosomes){
    my $in_file = File::Spec->catfile($dir,"merlin_$chr","merlin.err");
    if (!-f $in_file){
        print STDERR "\n*** mising $in_file\n\n";
    }
    else{
        open IN,"<",$in_file or die "cant't open $in_file $!\n";
        while(<IN>){
            my $line = $_;
            chomp;
            if(/\s*(\d+)\s+(\d+)\s+(\S+)\s+([0-9.eE,-]+)/){
                my @values = my ($fam,$ind,$marker,$ratio) = ($1,$2,$3,$4);
                my $format = "%s\t%s\t%s\t%s\n";
                printf ERR $format, $fam,$ind,$marker,$ratio;
            }
        }
    }
}

print "done\n";
exit;
