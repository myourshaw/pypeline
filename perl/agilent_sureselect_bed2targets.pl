#!/usr/bin/perl -w
use strict;
use warnings;

#/Users/myourshaw/lab/analysis/sequencing/exome/agilent/SureSelectHumanAllExonKit_S0274956_Broad_Build37_20100910_annotation/027495_D_BED_20100910.bed

print "usage: perl agilent_sureselect_bed2targets.pl <list of sureselect library bed files" unless @ARGV;
my @bed_files = @ARGV;

my %ucsc2n = {
chr1 => 1,
chr1_gl000191_random => ,
chr1_gl000192_random => ,
chr10 => 10,
chr11 => 11,
chr12 => 12,
chr13 => 13,
chr14 => 14,
chr15 => 15,
chr16 => 16,
chr17 => 17,
chr17_gl000205_random => ,
chr18 => 18,
chr19 => 19,
chr19_gl000209_random => ,
chr2 => 2,
chr20 => 20,
chr21 => 21,
chr22 => 22,
chr3 => 3,
chr4 => 4,
chr5 => 5,
chr6 => 6,
chr7 => 7,
chr8 => 8,
chr8_gl000196_random => ,
chr9 => 9,
chr9_gl000201_random => ,
chrUn_gl000211 => ,
chrUn_gl000212 => ,
chrUn_gl000213 => ,
chrUn_gl000214 => ,
chrUn_gl000215 => ,
chrUn_gl000217 => ,
chrUn_gl000219 => ,
chrUn_gl000220 => ,
chrUn_gl000221 => ,
chrUn_gl000222 => ,
chrUn_gl000223 => ,
chrUn_gl000227 => ,
chrUn_gl000228 => ,
chrUn_gl000229 => ,
chrUn_gl000236 => ,
chrUn_gl000238 => ,
chrUn_gl000241 => ,
chrUn_gl000244 => ,
chrUn_gl000246 => ,
chrUn_gl000249 => ,
chrX => 23,
chrY => 24,
}

foreach my $bed_file(@bed_files){
	my $range_ucsc_file = "$bed_file.ucsc.ranges.targets";
	my $range_ensembl_file = "$bed_file.ensembl.ranges.targets";
	my $positions_ucsc_file = "$bed_file.ucsc.positions.targets";
	my $positions_ensembl_file = "$bed_file.ensembl.positions.targets";
	my %targets;
	open BED, "<", $bed_file or die "can't open $bed_file $!";
	while (<BED>){
		if(/chr\S+\s\d+\s\d+/){
			my ($chrom, $chromStart, $chromEnd, $name, $score, $strand) = split;
			if (!defined $chrom | !defined $chromStart | !defined $chromEnd){
			print;
				
			}
			$targets{chrom2n($chrom)."\t$chromStart\t$chromEnd"}++;
			#last;
		}
	}
	my @targets = map {[split(/\t/,$_)]} keys %targets;
	@targets = sort {@$a[0] <=> @$b[0] || @$a[1] <=> @$b[1] || @$a[2] <=> @$b[2]} @targets;
	open RANGEUCSC,">",$range_ucsc_file;
	open RANGEENS,">",$range_ensembl_file;
	open POSUCSC,">",$positions_ucsc_file;
	open POSENS,">",$positions_ensembl_file;
	my $that = "";
	foreach my $target(@targets){
		my ($chr,$chromStart,$chromEnd) = @$target;
		$chromStart++;
		print RANGEUCSC n2ucsc($chr)."\t$chromStart\t$chromEnd\n";
		print RANGEENS n2ensembl($chr)."\t$chromStart\t$chromEnd\n";
		for my $pos ($chromStart .. $chromEnd){
			if("$chr:$pos" ne $that){
				print POSUCSC n2ucsc($chr)."\t$pos\n";
				print POSENS n2ensembl($chr)."\t$pos\n";
			}
			$that = "$chr:$pos";
		}
	}
}

sub chrom2n{
	my $chr = uc(shift);
	$chr =~ s/^chr//i;
	return $chr eq "X" ? 23 : $chr eq "Y" ? 24 : $chr eq "M" ? 25 : $chr eq "MT" ? 25 : $chr;
}
sub n2ucsc{
	my $n = shift;
	return $n  == 23 ? "chrX" : $n  == 24 ? "chrY" : $n == 25 ? "chrM" : "chr$n";	
}
sub n2ensembl{
	my $n = shift;
	return $n  == 23 ? "X" : $n  == 24 ? "Y" : $n == 25 ? "MT" : $n;	
}
