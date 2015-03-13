#!perl -w
use strict;
use warnings;

use File::Path;
use File::Spec;

#concatenate merlin -regress lod score tables

my $print_missing = 1;
my $print_absent = 0;
my $create_files = 1;
my $doall = 0;

my $dir = "/Volumes/raid/link/pedstatsmerlinall_2010-01-13-partitions";
my $resultsdate = "201001142127";

my $cat_dir = File::Spec->catfile($dir,"cat");
File::Path->make_path($cat_dir); 

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
my @folders = qw(table);#qw(freq ibd lod out pdf table zscore);
my @phases = $doall ? qw(phase1 phase2 phase1+2 phase2-1) : qw(phase1 phase2 phase1+2);
my @ethnics = $doall ? qw(all_ethnic non_white white) : qw(all_ethnic non_white white);
my @genders = $doall ? qw(all_gender female_family female_sibs male_family male_sibs) : qw(all_gender male_family);
my @adhds = $doall ? qw(broad narrow) : qw(broad narrow);
my @types = $doall ? qw(all_type combined combined_family combined_inattentive combined_inattentive_family hyperactive hyperactive_family inattentive inattentive_family) : qw(all_type combined combined_inattentive hyperactive inattentive);
my @chromosomes = ("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22","X","XY","Y"); #("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22","X","XY","Y");
my @chrsets = qw (autosome x);
my @analyses = ("ALL", "Pairs", "QTL");
my %nominals = ();
my %suggestives = ();
my %interestings = ();
my %significants = ();
my %highly_significants = ();

my $cat_header = "phase\tethnic\tgender\tadhd\ttype\tchr\tmarker\tvariable\th2\tsd\tinfo\tlod\tpvalue";
my $format = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%g\t%g\t%g\t%g\t%g\n";
my $tbl_count = 0;
my $missing_count =0;
my $absent_count = 0;
my $qtl_record_count = 0;
my %qtls =();
if($create_files){  
    my $cat_file = File::Spec->catfile($cat_dir,"regress-cat_$resultsdate.txt");
    open CAT,">",$cat_file or die "can't open $cat_file ($!)\n";
    print CAT "$cat_header\n";
}
foreach my $phase(@phases){
    foreach my $ethnic(@ethnics){
        foreach my $gender(@genders){
            foreach my $adhd(@adhds){
                foreach my $type(@types){
                    my @partition = ($codes{$phase},$codes{$ethnic},$codes{$gender},$codes{$adhd},$codes{$type});
                    foreach my $chrset(@chrsets){
                        my $in_file = File::Spec->catfile($dir,"regress-$phase-results_$resultsdate","table","$ethnic-$gender-$adhd-$type-$chrset-regress.tbl");
                        my $ped_file = File::Spec->catfile($dir,$phase,$ethnic,$gender,$adhd,$type,$chrset,"$chrset.ped");
                            if (-s $ped_file && !-s $in_file){
                            $missing_count++;
                            print "mising $in_file\n" unless !$print_missing;
                        }
                        elsif(-f $in_file){
                            $tbl_count++;
                            open IN,"<",$in_file or die "cant't open $in_file $!\n";
                            while(<IN>){
                                my $line = $_;
                                chomp;
                                if(/\s*([0-9XY]{1,3})\s+([\S]+)\s+\S+\s+(\S+)\s+\S+\s+([eE\d.-]+)\s+([eE\d.-]+)\s+([eE\d.-]+)\s+([eE\d.-]+)\s+([eE\d.-]+)/){
                                    my ($chr,$marker,$variable,$h2,$sd,$info,$lod,$pvalue) = ($1 eq "999"?"X":$1,$2,$3,$4,$5,$6,$7,$8);
                                    my @values = ($chr,$marker,$variable,$h2,$sd,$info,$lod,$pvalue);
                                    my @row = (@partition,@values);
                                    $qtl_record_count++;
                                    $qtls{$variable}++;
                                    printf CAT $format, @row unless !$create_files;
                                }
                            }
                        }
                        else{
                            $absent_count++;
                            print "absent  $in_file\n" unless !$print_absent;
                        }
                    }
                }
            }
        }
    }
}
print "\n$tbl_count tbl files\n";
print "$missing_count missing files (ped but no tbl)\n";
print "$absent_count absent files (no ped && no tbl)\n";
print "$qtl_record_count qtl records\n";
print map{"$qtls{$_} qtl $_ records\n"} sort keys %qtls;
print "\n     done\n";
exit;
