#!/usr/bin/perl -w
use strict;
use warnings;
use File::Spec;
#parameter keys must be lower case
my %params = (
    blat => "", #id of BLAT run
	flowcell => "FLOWCELL_ID", #flowcell id
    lanes => "1|2|3|4|5|6|7|8", #pipe-delimited list of lane numbers
	chromosome => "", #specific chromosome, because input file doesn't have chr yet
    chromosomes => "1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22|X|Y|M", #pipe-delimited list of chromosomes
    dir => ".", #parent directory of summary .txt files
    offset => 0 #added to position
);
my $usage = "perl summary.pl job_control_file [> output_file_root]\n";
unless (@ARGV) {
	print $usage;
	exit;
}
my @DIR;
my $input_file = "";
my $input_file1 = $ARGV[0];
my ($vol,$dir,$file) = File::Spec->splitpath($input_file1);
my $summary_file;
#remove special characters from file names
(my $clean_file = $file) =~ s/\W/_/;
#redirect STDOUT
my $gt = "&gt;"; #komodo saves command line > as &gt;
if(index($ARGV[$#ARGV],">") == 0){
    my $out_file_root = substr pop(@ARGV),1;
    $summary_file = $out_file_root."_sum.txt";
}
elsif ($ARGV[$#ARGV-1] eq ">" || $ARGV[$#ARGV-1] eq $gt) {
    my $out_file_root = pop(@ARGV);
    $summary_file = $out_file_root."_sum.txt";
    pop(@ARGV);
}
elsif (rindex($ARGV[$#ARGV-1],">") == length($ARGV[$#ARGV-1])-1) {
    my $out_file_root = pop(@ARGV);
    $summary_file = $out_file_root."_sum.txt";
    $ARGV[$#ARGV-1] = substr($ARGV[$#ARGV-1],0,rindex($ARGV[$#ARGV-1],">"));
}
elsif (rindex($ARGV[$#ARGV-1],$gt) == length($ARGV[$#ARGV-1])-1) {
    my $out_file_root = pop(@ARGV);
    $summary_file = $out_file_root."_sum.txt";
    $ARGV[$#ARGV-1] = substr($ARGV[$#ARGV-1],0,rindex($ARGV[$#ARGV-1],$gt));
}
else{
    my $out_file_root = File::Spec->catpath($vol,$dir,"$file");
    $summary_file = $out_file_root."_sum.txt";
}
while($ARGV[$#ARGV] eq ">" || $ARGV[$#ARGV] eq $gt){
    pop(@ARGV);
}
open OUT, ">", $summary_file or die "can't open $summary_file: $!";
print OUT "blat\tflowcell\tlane\ttag\tchr\tpos\tref\tAfwd\tArev\tCfwd\tCrev\tGfwd\tGrev\tTfwd\tTrev\tinsFwd\tinsRev\tdelFwd\tdelRev\tnoFwd\tnoRev\n";

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
    my $chromosomes = $params{chromosomes}; #eventually we need to filter on lanes and chromosomes
    my $summary_regex = 'summary_lane\.('.$params{lanes}.')\.([ACGT]{4})\.txt';
    my $lanes = $params{lanes};
    my $lanes_regex = '^lane('.$params{lanes}.')$';
    use File::Find;
    my @x = @_;
    find sub{
        my $found = $_;
        my $d = $File::Find::dir;
        print "$d\n";
        if($found =~ /$summary_regex/i){
			my $blat = $params{blat};
			my $flowcell = $params{flowcell};
            my $lane = $1;
            my $tag = $2;
			my $chr = $params{chromosome}; # this should eventually be in the file
            my $regex = '(\d+)\s([ACGT]|NA)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)';
           open(O,"< $_") or die "Can't open $_ $!\n";
            while(<O>){
                if(/$regex/i){
					my $pos = $1+$params{offset};
					my $ref = uc($2)eq"NA"?"N":uc($2);
					my $Afwd = $3;
					my $Tfwd = $4;
					my $Gfwd = $5;
					my $Cfwd = $6;
					my $noFwd = $7;
					my $insFwd = $8;
					my $delFwd = $9;
					my $Arev = $10;
					my $Trev = $11;
					my $Grev = $12;
					my $Crev = $13;
					my $noRev = $14;
					my $insRev = $15;
					my $delRev = $16;
                    print OUT "$blat\t$flowcell\t$lane\t$tag\t$chr\t$pos\t$ref\t$Afwd\t$Arev\t$Cfwd\t$Crev\t$Gfwd\t$Grev\t$Tfwd\t$Trev\t$insFwd\t$insRev\t$delFwd\t$delRev\t$noFwd\t$noRev\n";
                }
            }
        }
    }, shift;
	print STDERR "done\n";
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
