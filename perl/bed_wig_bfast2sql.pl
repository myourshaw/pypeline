#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;
use File::Spec;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use IO::Compress::Gzip qw(gzip $GzipError);
require File::Temp;
use File::Temp ();
use File::Temp qw/ :seekable /;

#"C:\Users\Michael\solexa\solexa_datasets\Reports\Mike_Yourshaw_targeted_sequencing\09-11-08\09-11-08.job"

#parameter keys must be lower case
my %params = (
    alignment => "BLAT_ID",
	file_regex => '\\\(\d{6}_HWI-EAS\S+_\S+)_s_(\d+)_all_seq\.fasta\.bfast\.report\\\([ACGT]+)\\\((?:mismatches)|(?:seq_density))_([ACGT]+)\.((?:bed)|(?:wig))\.?(gz)?',
    bed_regex => '(chr[0-9XYMT]{1,2}(?:_random)?)\s(\d+)\s(\d+)\s(INS:|DEL:)?([ACTG-]+)->([ACTG-]+)\((\d+):(\d+):\S+%\[F:(\d+)\S*\|R:(\d+)\S*\]\)\s+(\d+)\s+([+-])',
    wig_regex => '(chr[0-9XYMT]{1,2}(?:_random)?)\s(\d+)\s(\d+)\s(\d+)'
);
my $alignment = "BFAST_TARGETED_MYO1D_big_index_20080911";
my $file_regex = '\\\(\d{6}_HWI-EAS\S+_\S+)_s_(\d+)_all_seq\.fasta\.bfast\.report\\\([ACGT]+)\\\((?:mismatches)|(?:seq_density))_([ACGT]+)\.((?:bed)|(?:wig))\.?(gz)?';
#chr17	27929043	27929044	G->T(9:3:33.3%[F:1:11.1%|R:2:22.2%])	333	+	27929043	27929044	80,175,175	1	1	0
my $bed_regex = '(chr[0-9XYMT]{1,2}(?:_random)?)\s(\d+)\s(\d+)\s(INS:|DEL:)?([ACTG-]+)->([ACTG-]+)\((\d+):(\d+):\S+%\[F:(\d+)\S*\|R:(\d+)\S*\]\)\s+(\d+)\s+([+-])';
#chr17 27931425 27931426 1
my $wig_regex = '(chr[0-9XYMT]{1,2}(?:_random)?)\s(\d+)\s(\d+)\s(\d+)';
my $usage =
	"Merge bed and wig files, which may be gzipped, for sql database import. Usage:\n"
	."perl bed_wig_bfast2sql.pl [-b] [-w] [-z] job_control_file [> output_file_root]\n"
    ."-b ignore bed files\n"
    ."-w ignore wig files\n"
	."-z zip output\n"
    ."job_control_file lists fileglobs of wig and/or bed files\n";

unless (@ARGV) {
	print STDERR $usage;
	exit;
}
#options
my %opts; #sequence, genomic intervals, log
getopts("bwz", \%opts);
if ($opts{b} && $opts{w}) {die "Neither bed files nor wig files\n"};
my $input_file = "";
my $input_file1 = $ARGV[0];
my ($vol,$dir,$file) = File::Spec->splitpath($input_file1);
my $default_dir = $vol.$dir;
my $out_file_root;
#replace special characters in file names with _
(my $clean_file = $file) =~ s/\W/_/;
#redirect STDOUT
my $gt = "&gt;"; #komodo saves command line > as &gt;
if(index($ARGV[$#ARGV],">") == 0){
    $out_file_root = substr pop(@ARGV),1;
}
elsif ($ARGV[$#ARGV-1] eq ">" || $ARGV[$#ARGV-1] eq $gt) {
    $out_file_root = pop(@ARGV);
    pop(@ARGV);
}
elsif (rindex($ARGV[$#ARGV-1],">") == length($ARGV[$#ARGV-1])-1) {
    $out_file_root = pop(@ARGV);
    $ARGV[$#ARGV-1] = substr($ARGV[$#ARGV-1],0,rindex($ARGV[$#ARGV-1],">"));
}
elsif (rindex($ARGV[$#ARGV-1],$gt) == length($ARGV[$#ARGV-1])-1) {
    $out_file_root = pop(@ARGV);
    $ARGV[$#ARGV-1] = substr($ARGV[$#ARGV-1],0,rindex($ARGV[$#ARGV-1],$gt));
}
else{
    $out_file_root = File::Spec->catpath($vol,$dir,"$file");
}
while($ARGV[$#ARGV] eq ">" || $ARGV[$#ARGV] eq $gt){
    pop(@ARGV);
}
my $bed_file = $out_file_root."_$params{alignment}.bed.txt";
my $wig_file = $out_file_root."_$params{alignment}.wig.txt";
my $bed;
my $wig;
if($opts{z}){
	unless($opts{b}){
		$bed = new IO::Compress::Gzip "$bed_file.gz" or die "IO::Compress::Gzip failed for $bed_file: $GzipError\n";
	}
	unless($opts{w}){
		$wig = new IO::Compress::Gzip "$wig_file.gz" or die "IO::Compress::Gzip failed for $wig_file: $GzipError\n";
	}
}
else{
	unless ($opts{b}){
		open $bed, ">", $bed_file or die "can't open $bed_file: $!\n";
	}
	unless ($opts{w}){
		open $wig, ">", $wig_file or die "can't open $wig_file: $!\n";
	}
}
print $bed "blat\tflowcell\tlane\ttag\tchrom\tchromStart\tchromEnd\tstrand\trefNCBI\tobserved\tclass\treads\tvar_reads\tf_var_reads\tr_var_reads\tscore\n" unless ($opts{b});
print $wig "blat\tflowcell\tlane\ttag\tchrom\tchromStart\tchromEnd\treads\n" unless ($opts{w});
while(<>){
	#C:\Users\Michael\solexa\solexa_datasets\Reports\Mike_Yourshaw_targeted_sequencing\09-11-08\080???_HWI-EAS172_?????_s_?_all_seq.fasta.bfast.report\????\mismatches_????.bed
	#C:\Users\Michael\solexa\solexa_datasets\Reports\Mike_Yourshaw_targeted_sequencing\09-11-08\080???_HWI-EAS172_?????_s_?_all_seq.fasta.bfast.report\????\seq_density_????.wig.gz
	foreach my $file(glob){
		#if($file=~/\\(\d{6}_HWI-EAS\S+_\S+)_s_(\d+)_all_seq\.fasta\.bfast\.report\\([ACGT]+)\\((?:mismatches)|(?:seq_density))_([ACGT]+)\.((?:bed)|(?:wig))\.?(gz)?/i){
		if($file=~/$params{file_regex}/i){
			my ($flowcell,$lane,$tag,$file_type,$tag2,$bed_wig,$gz) = ($1,$2,$3,$4,$5,$6,$7,$8);
			unless(($bed_wig eq "bed" && $opts{b})||$bed_wig eq "wig" && $opts{w}){
				my $fh;
				my $fname;
				if($gz){
					$fh = new File::Temp(SUFFIX => '.'.$bed_wig) or croak();
					$fname = $fh->filename;
					gunzip $file => $fh or die "gunzip failed for $file: $GunzipError\n";
					$fh->seek(0,SEEK_SET);
				}
				else{
					open $fh,"<",$file or die "Can't open $file: $!\n";
				}
				while(<$fh>){
					if($bed_wig eq "bed"){
						if(/$params{bed_regex}/i){
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
							printf $bed ("%s\t%s\t%s\t%s\t%s\t%u\t%u\t%s\t%s\t%s\t%s\t%u\t%u\t%u\t%u\t%u\n"
							,$params{alignment},$flowcell,$lane,$tag,$chrom,$chromStart,$chromEnd,$strand,
							$refNCBI,$observed,$class,$reads,$var_reads,$f_var_reads,$r_var_reads,$score);
						}
					}
					elsif($bed_wig eq "wig"){
						if(/$wig_regex/i){
							my $chrom = $1;
							my $chromStart = $2; # 0 based position
							my $chromEnd = $3; #1 based position
							my $reads = $4;
							printf $wig ("%s\t%s\t%s\t%s\t%s\t%u\t%u\t%u\n"
							,$params{alignment},$flowcell,$lane,$tag,$chrom,$chromStart,$chromEnd,$reads);
						}
					}
				}
			}
		}
	}
}
print STDERR "done\n";
exit;
