#!perl -w
use strict;
use warnings;

use File::Path;
use File::Spec;
use File::Glob;
use Switch;

my $create_files = 1;
#concatenate merlin lod score tables
my $dir = "/Volumes/raid/link/pedstatsmerlinall_2010-01-18-partitions_chrs_age";
my $catdir = File::Spec->catfile($dir,"cat");
File::Path->make_path($catdir);  

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
my $phase='phase1+2';
my $ethnic='all_ethnic';
my $gender='all_gender';
my $adhd = 'broad';
my $type = 'all_type';
my @partition = ($codes{$phase},$codes{$ethnic},$codes{$gender},$codes{$adhd},$codes{$type});
my $nplqtl_header = "phase\tethnic\tgender\tadhd\ttype\tanalysis\tvariable\tmarker\tchr\tpos\tzscore\tdelta\tlod\tpvalue\texdelta\texlod\texpvalue";
my $format = "$codes{$phase}\t$codes{$ethnic}\t$codes{$gender}\t$codes{$adhd}\t$codes{$type}\t%s\t%s\t%s\t%s\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n";
my $tbl_count = 0;
my $missing_count =0;
my $absent_count = 0;
my $npl_all_record_count = 0;
my $npl_pairs_record_count = 0;
my $qtl_record_count = 0;
my %qtls =();
my $nplqtl_file = File::Spec->catfile($catdir,"nplqtlcat.txt");
my $vc_file = File::Spec->catfile($catdir,"vccat.txt");
my $regress_file = File::Spec->catfile($catdir,"regresscat.txt");
my @nplqtlfiles = glob("$dir/chr*/merlin-nonparametric.tbl");
if (@nplqtlfiles){
    open NPLQTL,">",$nplqtl_file or die "can't open $nplqtl_file ($!)\n"  unless !$create_files;
    print NPLQTL "$nplqtl_header\n";
    foreach my $file(@nplqtlfiles){
        $file =~ /\/chr([0-9XY]{1,2})\/merlin-nonparametric\.tbl/g or die "malformed $file\n";
        my $chr = "$1";
        print "$chr\n";
        open IN,"<",$file or die "can't open $file ($!)\n";
        while(<IN>){
            chomp;
            s/\tinf\t/\t0\t/;
            if(/\s*([0-9XY]{1,3})\s+([e\d.-]+)\s+(\S+)\s+(\S+)\s+\[(\S+)\]\s+([e\d.-]+)\s+([e\d.-]+)\s+([e\d.-]+)\s+([e\d.-]+)\s+([e\d.-]+)\s+([e\d.-]+)\s+([e\d.-]+)/){
                my ($analysis,$variable,$marker,$chr,$pos,$zscore,$delta,$lod,$pvalue,$exdelta,$exlod,$expvalue) = ($codes{$5},$4,$3,$chr,$2,$6,$7,$8,$9,$10,$11,$12);
                     my @values = ($analysis,$variable,$marker,$chr,$pos,$zscore,$delta,$lod,$pvalue,$exdelta,$exlod,$expvalue);
                     my @row = (@values);
                     switch ($analysis){
                         case "a"{
                             $npl_all_record_count++;
                         }
                         case "p"{
                             $npl_pairs_record_count++;
                         }
                         case "q"{
                             $qtl_record_count++;
                             $qtls{$variable}++;
                         }
                     }
                     printf NPLQTL $format, @row unless !$create_files;
            }
            else{
                print "$_\n" unless /(^CHR)|(^na)/;
            }
        }
    }
}
print "$tbl_count tbl files\n";
print "$npl_all_record_count npl all records\n";
print "$npl_pairs_record_count npl pairs records\n";
print "$qtl_record_count qtl records\n";
print map{"$qtls{$_} qtl $_ records\n"} sort keys %qtls;

#my @nplqtlfiles = glob("$dir/chr*/merlin-nonparametric.tbl");
#if (@nplqtlfiles){
#    open NPLQTL,">",$nplqtl_file or die "can't open $nplqtl_file ($!)\n";
#    print NPLQTL "$nplqtl_header\n";
#    foreach my $file(@nplqtlfiles){
#        open IN,"<",$file or die "can't open $file ($!)\n";
#        while(<IN>){
#            chomp;
#            s/\tinf\t/\t0\t/;
#            if(/\s*([0-9XY]{1,3})\s+([e\d.-]+)\s+(\S+)\s+(\S+)\s+\[(\S+)\]\s+([e\d.-]+)\s+([e\d.-]+)\s+([e\d.-]+)\s+([e\d.-]+)\s+([e\d.-]+)\s+([e\d.-]+)\s+([e\d.-]+)/){
#                my ($analysis,$variable,$marker,$chr,$pos,$zscore,$delta,$lod,$pvalue,$exdelta,$exlod,$expvalue) = ($codes{$5},$4,$3,$1 eq "999"?"X":$1,$2,$6,$7,$8,$9,$10,$11,$12);
#                     my @values = ($analysis,$variable,$marker,$chr,$pos,$zscore,$delta,$lod,$pvalue,$exdelta,$exlod,$expvalue);
#                     my @row = (@partition,@values);
#                     switch ($analysis){
#                         case "a"{
#                             $npl_all_record_count++;
#                         }
#                         case "p"{
#                             $npl_pairs_record_count++;
#                         }
#                         case "q"{
#                             $qtl_record_count++;
#                             $qtls{$variable}++;
#                         }
#                     }
#                     printf NPLQTL $format, @row unless !$create_files;
#            }
#            else{
#                print "$_\n" unless /(^CHR)|(^na)/;
#            }
#        }
#    }
#}
#print "$tbl_count tbl files\n";
#print "$npl_all_record_count npl all records\n";
#print "$npl_pairs_record_count npl pairs records\n";
#print "$qtl_record_count qtl records\n";
#print map{"$qtls{$_} qtl $_ records\n"} sort keys %qtls;


print "\ndone\n";
exit;
