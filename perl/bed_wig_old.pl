#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;
use File::Spec;
#parameter keys must be lower case
my %params = (
	flowcell => "FLOWCELL_ID", #flowcell id
    lanes => "1|2|3|4|5|6|7|8", #pipe-delimited list of lane numbers
    chromosomes => "1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22|X|Y|M", #pipe-delimited list of chromosomes
    dir => "." #parent directory of .bed files
);
my $usage = "perl bed_wig.pl [-b] [-w] job_control_file [> output_file_root]\n"
            ."-b ignore bed files\n"
            ."-w ignore wig files\n";
#options
my %opts; #sequence, genomic intervals, log
getopts("bw", \%opts);
if ($opts{b} && $opts{w}) {die "Neither bed files nor wig files\n"};
unless (@ARGV) {
	print $usage;
	exit;
}
my @DIR;
my $input_file = "";
my $input_file1 = $ARGV[0];
my ($vol,$dir,$file) = File::Spec->splitpath($input_file1);
my $bed_file;
my $wig_file;
#remove special characters from file names
(my $clean_file = $file) =~ s/\W/_/;
#redirect STDOUT
my $gt = "&gt;"; #komodo saves command line > as &gt;
if(index($ARGV[$#ARGV],">") == 0){
    my $out_file_root = substr pop(@ARGV),1;
    $bed_file = $out_file_root."_bed.txt";
    $wig_file = $out_file_root."_wig.txt";
}
elsif ($ARGV[$#ARGV-1] eq ">" || $ARGV[$#ARGV-1] eq $gt) {
    my $out_file_root = pop(@ARGV);
    $bed_file = $out_file_root."_bed.txt";
    $wig_file = $out_file_root."_wig.txt";
    pop(@ARGV);
}
elsif (rindex($ARGV[$#ARGV-1],">") == length($ARGV[$#ARGV-1])-1) {
    my $out_file_root = pop(@ARGV);
    $bed_file = $out_file_root."_bed.txt";
    $wig_file = $out_file_root."_wig.txt";
    $ARGV[$#ARGV-1] = substr($ARGV[$#ARGV-1],0,rindex($ARGV[$#ARGV-1],">"));
}
elsif (rindex($ARGV[$#ARGV-1],$gt) == length($ARGV[$#ARGV-1])-1) {
    my $out_file_root = pop(@ARGV);
    $bed_file = $out_file_root."_bed.txt";
    $wig_file = $out_file_root."_wig.txt";
    $ARGV[$#ARGV-1] = substr($ARGV[$#ARGV-1],0,rindex($ARGV[$#ARGV-1],$gt));
}
else{
    my $out_file_root = File::Spec->catpath($vol,$dir,"$file");
    $bed_file = $out_file_root."_bed.txt";
    $wig_file = $out_file_root."_wig.txt";
}
while($ARGV[$#ARGV] eq ">" || $ARGV[$#ARGV] eq $gt){
    pop(@ARGV);
}
unless ($opts{b}){
    open BED, ">", $bed_file or die "can't open $bed_file: $!";
    print BED "flowcellId\tlaneNumber\ttag\tchrom\tchromStart\tchromEnd\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts\n";
}
unless ($opts{w}){
    open WIG, ">", $wig_file or die "can't open $wig_file: $!";
    print WIG "flowcellId\tlaneNumber\ttag\tchrom\tchromStart\tchromEnd\treads\n";
}
#@ARGV = qw(.) unless @ARGV;
while (<>) {
    chomp;
    if($_){
		#comment
		if (/^\s*#/) { }
		#parameter
		elsif (/^\s*\$/) {
			parameter();
		}
        #directory
        if(/^\s*\$dir/){
            process_dir($params{dir});
        }
    }
}
exit;

sub process_dir{
    my $chromosomes = $params{chromosomes};
    my $mismatches_regex = 'chr('.$params{chromosomes}.')_mismatches_([ACGT]{4})\.bed';
    my $seq_density_regex = 'chr('.$params{chromosomes}.')_seq_density_([ACGT]{4})\.wig';
    my $lanes = $params{lanes};
    my $lanes_regex = 'lane_('.$params{lanes}.')_\d+$';
    my @x = @_;
    use File::Find;
    find sub{
        my $found = $_;
        my $d = $File::Find::dir;
        if($File::Find::dir =~ /$lanes_regex/i){
            my $lane = $1;
            if($found =~ /$mismatches_regex/i && !$opts{b}){
                my $chr = $1;
                my $tag = $2;
                my $regex = '(chr'.$chr.')\s(\d+)\s(\d+)\s(\S+)\s';
                open(B,"< $_") or die "Can't open $_ $!\n";
                while(<B>){
                    if(/$regex/i){
                        print BED "$params{flowcell}\t$lane\t$tag\t$_";
                    }
                }
            }
            if($found =~ /$seq_density_regex/i && !$opts{w}){
                my $chr = $1;
                my $tag = $2;
                my $regex = '(chr'.$chr.')\s(\d+)\s(\d+)\s(\d+)';
                open(W,"< $_") or die "Can't open $_ $!\n";
                while(<W>){
                    if(/$regex/i){
                        print WIG "$params{flowcell}\t$lane\t$tag\t$1\t$2\t$3\t$4\n";
                    }
                }
            }
        }
    }, shift;
}

#change parameter ($key value)
sub parameter {
	if(/s*\$([^\s=>#]+){1}[\s=>#]+([^#]+)#*$/){;
		my $key = lc($1);
		my $value = $2;
		$value =~ s/['"]+//g;
		if(defined($value)){
			if ($value =~ /^false$/i || $value =~ /^no$/i) {
				$value = 0;
			}
			elsif ($value =~ /^true$/i || $value =~ /^yes$/i) {
				$value = 1;
			}
		}
		if (exists($params{lc$key})) {
			$params{$key} = $value;
		}
		else {
			printf STDERR "%s isn't a job control parameter", $key;
		}
	}
} # !parameter

sub print_usage {
	print STDERR $usage . "Job control file parameters:\n";
	print STDERR map { "\t\$$_ = $params{$_}\n" } sort keys %params;
}

sub chr2numchr{
    my $chr = uc(shift);
    return $chr eq "X" ? 23 : $chr eq "Y" ? 24 : ($chr eq "M" || $chr eq "MT") ? 25 : $chr;
}

sub numchr2chr{
    my $chr = shift;
    return $chr == 23 ? "X" : $chr == 24 ? "Y" : $chr == 25 ? "M" : $chr;
}
