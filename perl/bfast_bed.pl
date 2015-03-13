#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;
use File::Spec::Functions;
use Cwd;

#C:\Users\Michael\solexa\solexa_datasets\Reports\Michael_Yourshaw_targeted_sequencing\targeted_alignment_many_indexes_chr17_27843500-28132200
my $alignment = "myo1d_targetd_bfast";
my $flowcell_dir_regex = '^(\d{6}_HWI-EAS172_[A-Za-z0-9]+)_s_(\d)_all_seq\.fasta\.bfast\.report$';
my $barcode_dir_regex = "([ACGTacgt]{4})";
my $bed_file_regex = "mismatches_([ACGTacgt]+).bed";
my $bed_record_regex = '(chr[0-9XYMT]{1,2}(?:_random)?)\s(\d+)\s(\d+)\s(INS:|DEL:)?([ACTG-]+)->([ACTG-]+)\((\d+):(\d+):\S+%\[F:(\d+)\S*\|R:(\d+)\S*\]\)\s+(\d+)\s+([+-])';
my $flowcell;
my $lane;
my $barcode;

my $usage = 'perl bfast_bed.pl [input_root_directory[\output_file]]'."\n"
."\tconsolidates bfast produced bed files into one flat file\n";
unless (@ARGV) {
	print $usage;
	exit;
}
my $dir = shift;
my $bed_file = catfile($dir,"bfast_bed.txt");
open BED, ">", $bed_file or die "Can't open $bed_file $!\n";
process_dir($dir);
print "done\n";
exit;

sub process_dir{
    use File::Find;
    find sub{
        my $found = $_;
        my $d = $File::Find::dir;
        if($_ =~ /$flowcell_dir_regex/i){
            $flowcell = $1;
            $lane = $2;
		}
		elsif($_ =~ /$barcode_dir_regex/i){
			my $barcode_dir = $1;
            if($found =~ /$bed_file_regex/i){
				$barcode = $1;
                open(B,"< $_") or die "Can't open $_ $!\n";
                print "processing mismatch file $d/$_\n";
                while(<B>){
					#if(substr($_,0,5) ne "track" && substr($_,0,5) ne "brows"){
						#print;
					#}
                    if(/$bed_record_regex/i){
                        my $chrom = $1;
                        my $chromStart = $2;
                        my $chromEnd = $3; #expect this to be the same as start1base
                        my $class = defined($4)?$4 eq "INS:"?"insertion":"deletion":"single";
                        my $refNCBI = $5;
                        my $observed = $6;
                        my $reads = $7;
                        my $var_reads = $8;
                        my $f_var_reads = $9;
                        my $r_var_reads = $10;
                        my $score = $11;
                        my $strand = $12;
                        printf BED ("%s\t%s\t%s\t%s\t%u\t%u\t%s\t%s\t%s\t%s\t%u\t%u\t%u\t%u\t%u\n")
                        ,$flowcell,$lane,$barcode
                        ,$chrom,$chromStart,$chromEnd,$strand,$refNCBI,$observed,$class,$reads,$var_reads,$f_var_reads,$r_var_reads,$score;
                    }
                }
            }
        }
    }, shift;
}

sub chr2numchr{
    my $chr = uc(shift);
    return $chr eq "X" ? 23 : $chr eq "Y" ? 24 : ($chr eq "M" || $chr eq "MT") ? 25 : $chr;
}

sub numchr2chr{
    my $chr = shift;
    return $chr == 23 ? "X" : $chr == 24 ? "Y" : $chr == 25 ? "M" : $chr;
}
