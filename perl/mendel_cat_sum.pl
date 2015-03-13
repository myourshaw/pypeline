#!perl -w
use strict;
use warnings;

use File::Path;
use File::Spec;

#concatenate merlin lod score tables

my $dir = "/Volumes/raid/link/mendel_2009-11-07";

my ($nominal_threshhold, $suggestive_threshhold, $interesting_threshhold, $significant_threshhold, $highly_significant_threshhold, $bizarre, $exclude) = (1.0, 2.2, 3.0, 3.6, 5.4, 9.999, -2);

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
    QTL => "q"
);
my %extensions = (
    freq => "freq",
    ibd => "ibd",
    lod => "lod",
    out => "out",
    pdf => "pdf",
    table => "tbl",
    zscore => "zscore"
);
my @folders = qw(sum);#qw(out sum);
my @phases = qw(phase1 phase2 phase1+2 phase2-1);
my @ethnics = qw(all_ethnic non_white white);
my @genders = qw(all_gender female_family female_sibs male_family male_sibs);
my @adhds = qw(broad narrow);
my @types = qw(all_type combined combined_family combined_inattentive combined_inattentive_family hyperactive hyperactive_family inattentive inattentive_family);
my @chromosomes = ("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22","X","XY","Y");
my @analyses = ("ALL", "Pairs", "QTL");
my %nominals = ();
my %suggestives = ();
my %interestings = ();
my %significants = ();
my %highly_significants = ();

my $cat_dir = File::Spec->catfile($dir,"cat-20091118-part");
File::Path->make_path($cat_dir); 
my $cat_header = "phase\tethnic\tgender\tadhd\ttype\tchr\tmarker\tpos\trecessive_blocks\tadditive_pairs\tadditive_all\tdominant_blocks";
my $format = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%g\t%g\t%g\t%g\t%g\n";
    my $cat_file = File::Spec->catfile($cat_dir,"mendel.sum.cat");
    my $sig_cat_file = File::Spec->catfile($cat_dir,"significant-mendel.sum.cat");
    open CAT,">",$cat_file or die "can't open $cat_file ($!)\n";
    #open SIG, ">",$sig_cat_file or die "can't open $sig_cat_file ($!)\n";
    #print "$cat_header\n";
    print CAT "$cat_header\n";
    #print SIG "$cat_header\n";
    foreach my $phase(@phases){
        foreach my $ethnic(@ethnics){
            foreach my $gender(@genders){
                foreach my $adhd(@adhds){
                    foreach my $type(@types){
                        my @partition = ($codes{$phase},$codes{$ethnic},$codes{$gender},$codes{$adhd},$codes{$type});
                        foreach my $chr(@chromosomes){
                            my $in_file = File::Spec->catfile($dir,"results-20091118-part","results","mendel_results-$phase","sum","$phase-chr$chr-$ethnic-$gender-$adhd-$type.sum");
                            my $ped_file = File::Spec->catfile($dir,$ethnic,$gender,$adhd,$type,"$phase-chr$chr","$phase-chr$chr.ped");
                            if (-f $ped_file && !-f $in_file){
                                print STDERR "mising $in_file\n";
                            }
                            elsif(-f $in_file){
                                open IN,"<",$in_file or die "cant't open $in_file $!\n";
                                while(<IN>){
                                    my $line = $_;
                                    chomp;
                                    if(/\s*(\S+)\s+([e\d.-]+)\s+([e\d.-]+)\s+([e\d.-]+)\s+([e\d.-]+)\s+([e\d.-]+)/i){
                                        my @values = my ($marker,$pos,$recessive,$pairs,$all,$dominant) = ($1 ne '--' ? $1 : "chr$chr"."cm$2",$2,$3,$4,$5,$6);
                                        my @row = (@partition,$chr,@values);
                                        printf CAT $format, @row;
                                    }
                                }
                            }
                            else{
                                print STDERR "no      $in_file\n";
                            }
                        }
                    }
                }
            }
        }
    }



print STDERR "\ndone\n";
exit;
