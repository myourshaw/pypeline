#!/usr/bin/perl -w
use strict;
use warnings;
use File::Glob ':glob';
use File::Spec;

#EDIT US
my %ids = (
	AATA => 13,
	ACGA => 14,
	ACTC => 15,
	CAGA => 16,
);
my $map = "genome_accurate";
my $flowcell = "090213_HWI-EAS335_309HH";
my $lane = 7;
my $dir = 'C:\!bak\solexa\090213_HWI-EAS335_309HH\Reports\genome_sequencing_custom_bfast_accurate_genome_alignment_hg18\\';
my $wig = 1;
my $bed = 1;
#EDIT THEM

if($wig){
	open WIG, ">", "$dir$map.$flowcell.$lane.wig" or die;
	print WIG "id\tmap\tflowcell\tlane\tbarcode\tchrom\tchromStart\tchromEnd\treads\n";
	my $wig_regex = '(chr[0-9XYMT]{1,2}(?:_random)?)\s(\d+)\s(\d+)\s(\d+)';
	foreach my $file(glob($dir."seq_density_*.wig")){
		if($file=~/.+seq_density_([ACGT]{4})T\.wig/){
			my $barcode=$1;
			my $id = $ids{$barcode};
			print "$file\n";
			open IN, "<$file" or die;
			while(<IN>) {
				chomp;
				if(/$wig_regex/){
					my $chrom = $1;
					my $chromStart = $2; # 0 based position
					my $chromEnd = $3; #1 based position
					my $reads = $4;
					printf WIG ("%s\t%s\t%s\t%s\t%s\t%s\t%u\t%u\t%u\n"
						,$id,$map,$flowcell,$lane,$barcode,$chrom,$chromStart,$chromEnd,$reads);
				}
			}
		}
	}
}
if($bed){
	open BED, ">", "$dir$map.$flowcell.$lane.bed" or die;
	print BED "id\tmap\tflowcell\tlane\tbarcode\tchrom\tchromStart\tchromEnd\tstrand\trefNCBI\tobserved\tclass\treads\tvar_reads\tf_var_reads\tr_var_reads\tscore\titemRGB\n";
	my $bed_regex = '(chr[0-9XYMT]{1,2}(?:_random)?)\s(\d+)\s(\d+)\s(INS:|DEL:)?([ACTG-]+)->([ACTG-]+)\((\d+):(\d+):\S+%\[F:(\d+)\S*\|R:(\d+)\S*\]\)\s+(\d+)\s+([+-])\s+(\d+)\s+(\d+)\s+(\d+,\d+,\d+)';
	foreach my $file(glob($dir."mismatches_*.bed")){
		if($file=~/.+mismatches_([ACGT]{4})T\.bed/){
			my $barcode=$1;
			my $id = $ids{$barcode};
			print "$file\n";
			open IN, "<$file" or die;
			while(<IN>) {
				chomp;
				if(/$bed_regex/){
					my $chrom = $1;
					my $chromStart = $2;
					my $chromEnd = $3;
					my $class = defined($4)?$4 eq "INS:"?"insertion":"deletion":"single";
					my $refNCBI = $5;
					my $observed = $6;
					my $reads = $7;
					my $var_reads = $8;
					my $f_var_reads = $9;
					my $r_var_reads = $10;
					my $score = $11;
					my $strand = $12;
					my $thickStart = $13;
					my $thickEnd = $14;
					my $itemRGB = $15;
					printf BED ("%s\t%s\t%s\t%s\t%s\t%s\t%u\t%u\t%s\t%s\t%s\t%s\t%u\t%u\t%u\t%u\t%u\t%s\n")
						,$id,$map,$flowcell,$lane,$barcode
						,$chrom,$chromStart,$chromEnd,$strand,$refNCBI,$observed,$class,$reads,$var_reads,$f_var_reads,$r_var_reads,$score,$itemRGB;
				}
			}
		}
	}
}
