#!/usr/bin/perl -w
use strict;
use warnings;
use POSIX qw(ceil floor);
use File::Spec;

my $usage = "perl count_tag_reads.pl job_file [output_file]\n"
	. "job_file :\n"
	. "[\$base_dir = base directory that contains flowcell dirs]\n"
	. "\tdefault: ".'/data/storage-1-00/solexa'."\n"
	. "[\$seq_dir = directory within each flowcell directory that contains sequence files]\n"
	. "\tdefault: ".'Data/C1-36_Firecrest1.8.28_25-01-2008_solexa/Bustard1.8.28_25-01-2008_solexa'."\n"
	. "[\$seq_file_left = regular expression for part of file name before lane number]\n"
	. "\tdefault: ".'s_'."\n"
	. "[\$seq_file_right = regular expression for part of file name after lane number]\n"
	. "\tdefault: ".'_(\d+)_seq.txt'."\n"
	. "[\$seq_line = regular expression for line in sequence file]\n"
	. "\tdefault: ".'^(\d+)\s+(\d+)\s+-?(\d+)\s+-?(\d+)\s+([.ACGTacgt]+).*$'."\n"
	. "flowcell_directory lane_number [[lane_number] ...]\n"
	. "[flowcell_directory lane_number [[lane_number] ...]]\n"
	. "[...]\n";
my %params = (
	base_dir => '/data/storage-1-00/solexa', #dir that contains flowcell dirs
    seq_dir => 'Data/C1-36_Firecrest1.8.28_25-01-2008_solexa/Bustard1.8.28_25-01-2008_solexa', #dir that contains sequence files
	seq_file_left => 's_',
	seq_file_right => '_(\d+)_seq.txt',
	seq_line => '^(\d+)\s+(\d+)\s+-?(\d+)\s+-?(\d+)\s+([.ACGTacgt]+).*$'
);
@ARGV or die $usage;
my $outfile = "";
if (@ARGV>1){
	$outfile = pop @ARGV;
	open OUT,">",$outfile or die "Can't open $outfile: $!";
}
@ARGV or die "No job file\n$usage";
my $files = 0;
my $lines = 0;
my %counts;
#read job file
while(<>){
    chomp;
    if($_){
        #comment
        if (/^\s*#/) { }
        #parameter
        elsif (/^\s*\$/) {
            parameter();
        }
        #flowcell and lanes
        elsif(/^\s*(\S+|(?:".+"))\s+((?:\d+[^0-9#]+)*\d*)/g){
            my $cell = $1;
            my @lanes = split /\D+/,$2;
            my $dir = File::Spec->catdir($params{base_dir},$cell,$params{seq_dir});
			-d $dir or die "$dir not found\n";
			my $seqfileregex = $params{seq_file_left}.join("|",@lanes).$params{seq_file_right};
            print "reading cell $1 lanes ".join(", ",@lanes)."\n";
			use File::Find;
			find sub{
				print STDERR $seqfileregex."\n";
				print STDERR $_."\n";
				if(/$seqfileregex/i){
					my $lane = $1;
					my $tile = $2;
					my $seqfile = $_;
					open(SEQFILE,"< $seqfile") or die "Can't open $seqfile $!\n";
					$files++;
					while(<SEQFILE>){
						chomp;
						if(/$params{seq_line}/i){
							$lines++;
							my $input_line = $_;
							my ($lane,$tile,$tileX,$tileY,$seq) = ($1,$2,$3,$4,$5);
							my $tag = substr($seq,0,4);
							$counts{$cell."_".$lane."_".$tag}++;
						}
						else{print STDERR "unrecognized sequence in $seqfile\n\t$_\n";}
					}
				}
			}, $dir;
        }
    }
}
print STDERR "$files files, $lines lines\n";
print STDERR map { "$_\t$counts{$_}\n" } sort {$a cmp $b} keys %counts;
print OUT map { "$_\t$counts{$_}\n" } sort {$a cmp $b} keys %counts if $outfile;
exit;

#change parameter ($key value)
sub parameter {
	if(/\s*\$([^\s=>#]+){1}[\s=>#]+([^#]+)#*$/){
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
		if (exists($params{$key})) {
			if ($params{$key} ne $value){
				print STDERR sprintf("%s changed from %s to %s\n",$key,$params{$key},$value) unless $params{$key} eq $value;
				$params{$key} = $value;
			}
		}
		else {
			print STDERR sprintf("%s isn't a setting", $key);
		}
	}
} # !parameter
