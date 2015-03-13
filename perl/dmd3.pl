#!/usr/bin/perl -w
use strict;
use warnings;
use feature "switch";

#sudo  ./gfServer -stepSize=5 start localhost 17779 hg18.2bit &
#./gfPcr -maxSize=10000 -out=bed localhost 17779 /data/genomes /home/myourshaw/Documents/Lab/DMD/final/DMD_primers_gfPcr.txt /home/myourshaw/Documents/Lab/DMD/final/DMD_primers_gfPcr.bed
#./gfPcr -maxSize=10000 -out=fa localhost 17779 /data/genomes /home/myourshaw/Documents/Lab/DMD/final/DMD_primers_gfPcr.txt /home/myourshaw/Documents/Lab/DMD/final/DMD_primers_gfPcr.fa

use POSIX qw(ceil floor);
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
#BioPerl
use Bio::Seq;
use Bio::SeqIO;
use Bio::Tools::Primer3;
#EnsEMBL Perl API
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::Utils::Sequence qw (sequence_with_ambiguity);

my $usage = "perl dmd3.pl input_file\nCreate gfPcr or isPcr batch file or raindance primer list from primer3 output\n";
unless (@ARGV || $ARGV[0] ne 'help') {
	print $usage;
	exit;
}
my $omit_dupes=1;
my $primer="";# = 1;#="" for primer3's first choice
my $out_type = "raindance";
my $pcr_type="gfPcr"; #"gfPcr" | isPcr
my $no_out_file=0;
my $maxProductSize = 10000; #used by gfPcr
my $Amplicon_Name='';
my $Fwd_Primer_sequence='';
my $Rev_Primer_sequence='';
my $Amplicon_sequence='';
my $seq='';
my @primer_left=(0,0);
my @primer_right=(0,0);
my $primer_product_size=0;
my($num,$name,$chr,$start,$strand_symbol,$end,$roiStart,$roiEnd);
my %forward_primers;
my %reverse_primers;
my %both_primers;
my %primer_names;
my %missing_primers;
unless($no_out_file){
	my $out = $ARGV[0]."_$out_type".".txt";
	my $pcr_out = $ARGV[0]."_$pcr_type".".txt";
	open OUT, ">", $out;
	open PCR, ">", $pcr_out;
}
print "Amplicon_Name\tFwd_Primer_sequence\tRev_Primer_sequence\tAmplicon_sequence\n";
print OUT "Amplicon_Name\tFwd_Primer_sequence\tRev_Primer_sequence\tAmplicon_sequence\n" unless $no_out_file;
while(<>){
	chomp;
	if(/^=$/){
		$Amplicon_Name = sprintf("%s_%s_%s:%u%s%u[%u%s%u]",$num,$name,$chr,$start+$primer_left[0]-1,$strand_symbol,$start+$primer_left[0]-1+$primer_product_size-1,$roiStart,$strand_symbol,$roiEnd);
		$forward_primers{"$Fwd_Primer_sequence"}++;
		$reverse_primers{"$Rev_Primer_sequence"}++;
		$both_primers{"$Fwd_Primer_sequence|$Fwd_Primer_sequence"}++;
		$primer_names{"$Fwd_Primer_sequence|$Fwd_Primer_sequence"}.="$Amplicon_Name|";
		unless ($Fwd_Primer_sequence && $Rev_Primer_sequence){
			$missing_primers{$Amplicon_Name}++;
		}
		unless($omit_dupes && $both_primers{"$Fwd_Primer_sequence|$Fwd_Primer_sequence"}>1){
			given ($out_type){
				my $record;
				when ("raindance"){
					$Amplicon_sequence=substr($seq,$primer_left[0]-1,$primer_product_size);
					$record = sprintf("%s\t%s\t%s\t%s\n",$Amplicon_Name,$Fwd_Primer_sequence,$Rev_Primer_sequence,$Amplicon_sequence);
					print $record;
					print OUT $record unless $no_out_file;
				}
				default{
					$Amplicon_sequence=substr($seq,$primer_left[0]-1,$primer_product_size);
					$record = sprintf("%s\t%s\t%s\t%s\n",$Amplicon_Name,$Fwd_Primer_sequence,$Rev_Primer_sequence,$Amplicon_sequence);
					print $record;
					print OUT $record unless $no_out_file;
				}
			}
			given ($pcr_type){
				my $record;
				when ("gfPcr"){
					$record = sprintf("%s\t%s\t%s\t%u\n",$Amplicon_Name,$Fwd_Primer_sequence,$Rev_Primer_sequence,$maxProductSize);
					print PCR $record unless $no_out_file;
				}
				when ("isPcr"){
					$record = sprintf("%s\t%s\t%s\n",$Amplicon_Name,$Fwd_Primer_sequence,$Rev_Primer_sequence);
					print PCR $record unless $no_out_file;
				}
			}
		}
		$Amplicon_Name='';
		$Fwd_Primer_sequence='';
		$Rev_Primer_sequence='';
		$Amplicon_sequence='';
		$seq='';
		@primer_left=(0,0);
		@primer_right=(0,0);
		$primer_product_size=0;
	}
	elsif(/^PRIMER_SEQUENCE_ID=(\S+)/i){
		my $id=$1;
		$id =~ /([^_]+)_([^_]+)_([0-9MXY]{1,2}):(\d+)([+-])(\d+)\[(\d+)[+-](\d+)/i;
		($num,$name,$chr,$start,$strand_symbol,$end,$roiStart,$roiEnd)=($1,$2,$3,$4,$5,$6,$7,$8);
	}
	elsif(/^SEQUENCE=([A-Z]+)/i){
		$seq=$1;
	}
	elsif(/^PRIMER_LEFT_(\d?)_?SEQUENCE=([ACGT]+)/i){
		if($1 eq $primer){
		   $Fwd_Primer_sequence=$2;
		}
	}
	elsif(/^PRIMER_RIGHT_(\d?)_?SEQUENCE=([ACGT]+)/i){
		if($1 eq $primer){
			$Rev_Primer_sequence=$2;
		}
	}
	elsif(/^PRIMER_(\d?)_?LEFT=(\d+),(\d+)/i){
		if($1 eq $primer){
			@primer_left=($2,$3);
		}
	}
	elsif(/^PRIMER_(\d?)_?RIGHT=(\d+),(\d+)/i){
		if($1 eq $primer){
			@primer_right=($2,$3);
		}
	}
	elsif(/^PRIMER_(\d?)_?PRODUCT_SIZE=(\d+)/i){
		if($1 eq $primer){
			$primer_product_size=$2;
		}
	}
}
if(%missing_primers){
	printf "\n%u MISSING PRIMERS\n",scalar keys %missing_primers ;
	print map {"$_\n"} sort keys %missing_primers;
}
if(%forward_primers){
	printf "\n%u DUPLICATE FORWARD PRIMERS\n",scalar keys %forward_primers;
	print map {$_&&$forward_primers{$_}>1?"$_\t$forward_primers{$_}\n":""} sort keys %forward_primers;
}
if(%reverse_primers){
	printf "\n%u DUPLICATE REVERSE PRIMERS\n",scalar keys %reverse_primers;
	print map {$_&&$reverse_primers{$_}>1?"$_\t$reverse_primers{$_}\n":""} sort keys %reverse_primers;
}
if(%both_primers){
	printf "\n%u DUPLICATE PRIMER PAIRS\n",scalar keys %both_primers;
	print map {$_&&$both_primers{$_}>1?"$primer_names{$_}\t$_\t$both_primers{$_}\n":""} sort keys %both_primers;
}
