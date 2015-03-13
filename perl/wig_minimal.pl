#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;
use File::Spec;
#parameter keys must be lower case
my %params = (
    alignment => "alignment_ID",
	flowcell => "FLOWCELL_ID", #flowcell id
    lanes => "1|2|3|4|5|6|7|8", #pipe-delimited list of lane numbers
    chromosomes => "1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22|X|Y|M", #pipe-delimited list of chromosomes
    lanes_regex => 'lane_(1|2|3|4|5|6|7|8)_\d+$',
    mismatches_regex => 'chr([0-9XYMT]{1,2})_(\d+)_(\d+)_mismatches_([ACGT]{4})\.bed',
    seq_density_regex => 'chr([0-9XYMT]{1,2})_(\d+)_(\d+)_seq_density_([ACGT]{4})\.wig',
    #mismatches_regex => 'chr([0-9XYMT]{1,2})(_\d+_\d+)?_mismatches_([ACGT]{4})\.bed',
    #seq_density_regex => 'chr([0-9XYMT]{1,2})(_\d+_\d+)?_seq_density_([ACGT]{4})\.wig',
    #mismatches_regex => 'chr([0-9XYMT]{1,2})(_\d+_\d+)?_mismatches_([ACGT]{4}|notag)\.bed',
    #seq_density_regex => 'chr([0-9XYMT]{1,2})(_\d+_\d+)?_seq_density_([ACGT]{4}|notag)\.wig',
    #mismatches_regex => 'chr([0-9XYMT]{1,2})(_\d+_\d+)?_mismatches_(anytag)\.bed',
    #seq_density_regex => 'chr([0-9XYMT]{1,2})(_\d+_\d+)?_seq_density_(anytag)\.wig',
    offset => 0, #
	include_wig_alignment => 0, #include in wig output
	include_wig_flowcell => 0, #
	include_wig_lane => 1, #
	include_wig_barcode => 0, #
	include_wig_chrom => 1, #
	include_wig_chromStart => 0, #
	include_wig_chromEnd => 1, #
	include_wig_reads => 1, #
	dir => "" #parent directory of .bed/.wig files
);
my $usage = "perl bed_wig.pl [-b] [-w] job_control_file [> output_file_root]\n"
            ."-b ignore bed files\n"
            ."-w ignore wig files\n"
            ."job file parameter defaults:\n"
            .'$alignment=alignment_ID'."\n"
            .'$flowcell=FLOWCELL_ID'."\n"
            .'$lanes=1|2|3|4|5|6|7|8'."\n"
            .'$chromosomes=1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22|X|Y|M'."\n"
            .'$mismatches_regex=chr([0-9XYMT]{1,2})_(\d+)_(\d+)_mismatches_([ACGT]{4})\.bed'."\n"
            .'$seq_density_regex=chr([0-9XYMT]{1,2})_(\d+)_(\d+)_seq_density_([ACGT]{4})\.wig'."\n"
            .'$lanes_regex=lane_(1|2|3|4|5|6|7|8)_\d+$'."\n"
            .'$dir=<parent directory of files to process> default is same directory as job file'
            .'$go starts processing of files in curernt parent directory)'."\n";
#options
my %opts; #sequence, genomic intervals, log
getopts("bw", \%opts);
if ($opts{b} && $opts{w}) {die "Neither bed files nor wig files\n"};
unless (@ARGV) {
	print $usage;
	exit;
}
my $input_file = "";
my $input_file1 = $ARGV[0];
my ($vol,$dir,$file) = File::Spec->splitpath($input_file1);
my $default_dir = $vol.$dir;
my $out_file_root;
my $bed_file;
my %bed;
my $wig_file;
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
$bed_file = $out_file_root.".bed";
$wig_file = $out_file_root.".wig";
while($ARGV[$#ARGV] eq ">" || $ARGV[$#ARGV] eq $gt){
    pop(@ARGV);
}
unless ($opts{b}){
    open BED, ">", $bed_file or die "can't open $bed_file: $!";
    print BED "alignment\tflowcell\tlane\tbarcode\tchrom\tchromStart\tchromEnd\tstrand\trefNCBI\tobserved\tclass\treads\tvar_reads\tf_var_reads\tr_var_reads\tscore\n";
}
unless ($opts{w}){
    open WIG, ">", $wig_file or die "can't open $wig_file: $!";
	my $header = $params{include_wig_alignment}?"alignment\t":"";
	$header .= $params{include_wig_flowcell}?"flowcell\t":"";
	$header .= $params{include_wig_lane}?"lane\t":"";
	$header .= $params{include_wig_barcode}?"barcode\t":"";
	$header .= $params{include_wig_chrom}?"chrom\t":"";
	$header .= $params{include_wig_chromStart}?"chromStart\t":"";
	$header .= $params{include_wig_chromEnd}?"chromEnd\t":"";
	$header .= $params{include_wig_reads}?"reads\t":"";
	print WIG "$header\n" if $header;
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
        if(/^\s*\$go/){
            unless(-d $params{dir}){
                $params{dir} = $default_dir;
            }
            process_dir($params{dir});
        }
    }
}
print "done\n";
exit;

sub process_dir{
    my $chromosomes = $params{chromosomes};
    my $mismatches_regex = $params{mismatches_regex};
    my $seq_density_regex = $params{seq_density_regex};
    my $lanes = $params{lanes};
    my $lanes_regex = $params{lanes_regex};
    my @x = @_;
    use File::Find;
    find sub{
        my $found = $_;
        my $d = $File::Find::dir;
        if($d =~ /$lanes_regex/i){
            my $lane = $1;
            if($found =~ /$mismatches_regex/i && !$opts{b}){
                my $roi_chr = $1;
                my $roi_start = $2;
                my $roi_end = $3;
                my $barcode;
                $barcode = $4 or $barcode = "";
                my $bed_regex = '(chr[0-9XYMT]{1,2}(?:_random)?)\s(\d+)\s(\d+)\s(INS:|DEL:)?([ACTG-]+)->([ACTG-]+)\((\d+):(\d+):\S+%\[F:(\d+)\S*\|R:(\d+)\S*\]\)\s+(\d+)\s+([+-])';
                open(B,"< $_") or die "Can't open $_ $!\n";
                print "processing mismatch file $d/$_\n";
                while(<B>){
                    if(/$bed_regex/i){
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
                        my $bed_key = "$params{alignment}|$params{flowcell}|$lane|$barcode|$chrom|$chromStart";
                        printf BED ("%s\t%s\t%s\t%s\t%s\t%u\t%u\t%s\t%s\t%s\t%s\t%u\t%u\t%u\t%u\t%u\n")
                        ,$params{alignment},$params{flowcell},$lane,$barcode
                        ,$chrom,$chromStart,$chromEnd,$strand,$refNCBI,$observed,$class,$reads,$var_reads,$f_var_reads,$r_var_reads,$score;

                    }
                }
            }
            if($found =~ /$seq_density_regex/i && !$opts{w}){
                my $roi_chr = $1;
                my $roi_start = $2;
                my $roi_end = $3;
                my $barcode;
                $barcode = $4 or $barcode = "";
                my $wig_regex = '(chr[0-9XYMT]{1,2}(?:_random)?)\s(\d+)\s(\d+)\s(\d+)';
                open(W,"< $_") or die "Can't open $_ $!\n";
                print "processing density file $d/$_\n";
                while(<W>){
                    if(/$wig_regex/i){
                        my $chrom = $1;
                        my $chromStart = $2; # 0 based position
                        my $chromEnd = $3; #1 based position
                        my $reads = $4;
						my $record = $params{include_wig_alignment}?"$params{alignment}\t":"";
						$record .= $params{include_wig_flowcell}?"$params{flowcell}\t":"";
						$record .= $params{include_wig_lane}?"$lane\t":"";
						$record .= $params{include_wig_barcode}?"$barcode\t":"";
						$record .= $params{include_wig_chrom}?"$chrom\t":"";
						$record .= $params{include_wig_chromStart}?"$chromStart\t":"";
						$record .= $params{include_wig_chromEnd}?"$chromEnd\t":"";
						$record .= $params{include_wig_reads}?"$reads":"";
						print WIG "$record\n" if $record;
                    }
                }
            }
        }
    }, shift;
}

#change parameter ($key value)
sub parameter {
	if(/\s*\$([^\s=>]+){1}[\s=>]+(.+){1}$/){
		my $key = lc($1);
		my $value = $2;
        #strip leading and trailing spaces
        $value =~ s/^\s*//;
        $value =~ s/\s*$//;
        #remove leading and trailing " or  '
        $value =~ s/^['"]?//;
        $value =~ s/['"]?$//;
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
