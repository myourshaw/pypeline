#!/usr/bin/perl -w
use strict;
use warnings;

use File::Spec;
use File::Glob;

my $check_missing_files = 1;
my $glob = "/Volumes/raid/merlin_files/merlin_files_20090927/pedcheck/pedcheck*.comment";
#my $glob = "/Volumes/raid/merlin_files/merlin_files_20090927/pedcheck/pedcheck-all-chr*-all_ethnic-all_gender-adhd_broad_life-all_type.comment";

my $dir = "/Volumes/raid/merlin_files/merlin_files_20090927/pedcheck";
my @ethnics = qw(all_ethnic non_white white);
my @genders = qw(all_gender female_family female_sibs male_family male_sibs);
my @adhd_lifes = qw(adhd_broad_life adhd_narrow_life);
my @adhd_types = qw(all_type combined hyperactive inattentive);
my @chromosomes = (1..22, "X", "XY", "Y");
#my @chromosomes = (1..22);
my @studies = qw(all cidr phase1);
my %chrs;
my %markers;
my %pedigrees;
my %zero;
my $pedigree;
if($check_missing_files){
    foreach my $ethnic(@ethnics){
        foreach my $gender(@genders){
            foreach my $adhd_life(@adhd_lifes){
                foreach my $adhd_type(@adhd_types){
                    foreach my $chr(@chromosomes){
                        foreach my $study(@studies){
                            my $pedcheck_file = File::Spec->catfile($dir,"pedcheck-$study-chr$chr-$ethnic-$gender-$adhd_life-$adhd_type.comment");
                            if (!-f $pedcheck_file){
                                print STDERR "*** mising $pedcheck_file\n";
                            }
                        }
                    }
                }
            }
        }
    }
    exit;            
}
foreach my $file(glob($glob)){
    $file =~ /pedcheck-.+chr([0-9XY]{1,2})\-.+\.comment/;
    my $chr = $1;
    $chrs{$chr}++;
    open IN,"<",$file or die "can't open $file ($!)";
    while(<IN>){
        if(/#\s+Pedigree\s+(\d+)/){
            $pedigree = $1;
        }
        elsif(/#\s+marker\s+(\S+)\s+zeroed out/){
            $markers{"$1"}++; #count pedigrees per zeroed out markers
            $pedigrees{$pedigree}++; #count zeroed out markers per pedigree
            $zero{"$pedigree\t$1"}++; #pedigree,marker to zero out
        }
    }
}
print "chr\tstrata\n";
print map {"$_\t$chrs{$_}\n"} sort keys %chrs;
print "\n\n\n\nmarker\tpedigrees\n";
print map {"$_\t$markers{$_}\n"} sort keys %markers;
print "\n\n\n\npedigree\tmarkers\n";
print map {"$_\t$pedigrees{$_}\n"} sort keys %pedigrees;
print "\n\n\n\npedigree\tmarker\tcount\n";
print map {"$_\t$zero{$_}\n"} sort keys %zero;
