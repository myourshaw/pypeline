#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;
use File::Spec;
#parameter keys must be lower case
my %params = (
	flowcell => "FLOWCELL_ID", #flowcell id
    lanes => "1|2|3|4|5|6|7|8", #pipe-delimited list of lane numbers
    seq_file_regex => 's_(\d)_(\d{4})_seq\.txt',
    read_seq_regex => '(\d+)\t(\d+)\t(\d+)\t(\d+)\t([ACGT]+)',
    min_length => 36,
    offset => 0,
    dir => "." #parent directory of .bed files
);
my $usage = "perl seq.pl job_control_file [> output_file_root]\n";
#options
unless (@ARGV) {
	print $usage;
	exit;
}
my @DIR;
my $input_file = "";
my $input_file1 = $ARGV[0];
my ($vol,$dir,$file) = File::Spec->splitpath($input_file1);
#remove special characters from file names
(my $clean_file = $file) =~ s/\W/_/;
#redirect STDOUT
my $gt = "&gt;"; #komodo saves command line > as &gt;
my $out_file_root;
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
my $out_file = $out_file_root.".seq";
open(OUT,">", $out_file) or die "Can't open $out_file!\n";
print OUT "flowcell\tlane\ttile\tx\ty\tseq\n";
while($ARGV[$#ARGV] eq ">" || $ARGV[$#ARGV] eq $gt){
    pop(@ARGV);
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
    my $read_seq_regex = $params{read_seq_regex};
    my $lanes = $params{lanes};
    my $seq_file_regex = $params{seq_file_regex};
    my @x = @_;
    use File::Find;
    find sub{
        my $found = $_;
        my $d = $File::Find::dir;
        if($found =~ /$seq_file_regex/i){
            open(S,"<", $_) or die "Can't open $_ $!\n";
            while(<S>){
                if(/$read_seq_regex/i){
                    my $lane = $1;
                    my $tile = $2;
                    my $x = $3;
                    my $y = $4;
                    my $read = $5;
                    if (length($read)>=$params{min_length}){
                        printf OUT "%s\t%u\t%u\t%u\t%u\t%s\n",
                        $params{flowcell},$lane,$tile,$x,$y,$read;
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
			printf STDERR "%s isn't a job control parameter\n", $key;
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
