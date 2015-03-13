#!/usr/bin/perl -w
use strict;
use warnings;

use Bio::Seq;
use Bio::SeqIO;
use Getopt::Std;

printf "%d %s\n", my $i++, $_ for @INC;
#get parameters from command line
my %opts;
getopt('fdope', \%opts);

#print usage if no file parameter
if(!$opts{f}&&!$opts{d})
{
	print "Create PieceMaker input file from fasta file of gene sequences.\n";
	print "Description format: >ACSL4:chromosome:ncbi36:x:108770622:108863877:-1\n";
	print "Gene format: upstreamEXONintron...EXONdownstream\n";
	print "Usage: perl f2p -d|-f <file> [-o <file>] [-g <parent>] [-p n[,m]] [-i x[,y]]\n";
	print "\t-f <fasta file>\n";
	print "\t-d <directory of fasta files>\n";
	print "\t-d or -f required. -f trumps -d\n";
	print "\t-o <PieceMaker file> (default is fasta file or directory with '.piecemaker' added\n";
	print "\t-p <n,m> padding of n bases downstream from start and m bases upstream from end (default 400,500)\n";
	print "\t-e <x,y> exons only, including x bases upstream of exons 2-ultimate, y bases downstream of exons 1-penultimate (default 25,10)\n";
}
else
{
	#set values from user-specificed parameters
	#my $input_file = $opts{f};
	#my $directory = $opts{d};
	my @opt_p = $opts{p} ? split(/,/,$opts{p}) : (0,0);
	my $pad_up = $opt_p[0] || 400;
	my $pad_down = $opt_p[1] || 500;
	my @opt_e = $opts{e} ? split(/,/,$opts{e}) : (0,0);
	my $exon_up = $opt_e[0] || 25;
	my $exon_down = $opt_e[1] || 10;
	
	#process a file or all files in a directory
	my($file,$output_file);
	if($opts{f})
	{
		$output_file = $opts{o} || $opts{f} . ".piecemaker";
		open (PIECEMAKERFILE,">",$output_file)
			or die "Couldn't open '$output_file' for writing.\n";
		CreatePieceMakerFile($opts{f});
		close(PIECEMAKERFILE)
			or die "Couldn't close '$output_file', stopped";
	}
	else
	{
		my @dir = split(/\\/,$opts{d});
		my $filename = $opts{o} || $dir[$#dir] . ".piecemaker";
		$output_file = "$opts{d}\\$filename";
		open (PIECEMAKERFILE,">",$output_file)
			or die "Couldn't open '$output_file' for writing.\n";
		opendir(DIR, $opts{d}) or die "can't open $opts{d}: $!";
		while (defined($file = readdir(DIR)))
		{
			next if $file =~ /^\.?$/;
			next if $file =~ /^\.\.?$/;
			CreatePieceMakerFile("$opts{d}\\$file");
		}
		close(PIECEMAKERFILE)
			or die "Couldn't close '$output_file', stopped";
	}
	sub CreatePieceMakerFile
	{
		my $input_file = $_[0];
		#get sequences from fasta file
		my $format = "fasta";
		my ($sequence, @sequences);
		my $sequence_file = Bio::SeqIO->new(-file => $input_file, -format => $format );
		my $i = 0;
		while ($sequence = $sequence_file->next_seq)
		{
			$sequences[$i++] = $sequence;
		}
		
		#create PieceMaker input file	
		if($i > 0)
		{
			$i--;
			print PIECEMAKERFILE "#from $format file: $input_file\n"
				or die "Couldn't write to '$output_file'.\n";
			print PIECEMAKERFILE "#Target_name\tParent_gene\tStart\tROI_begin\tROI_end\tSequence\n"
				or die "Couldn't write to '$output_file', stopped";
			print "#from $format file: $input_file\n";
			print "#Target_name\tParent_gene\tStart\tROI_begin\tROI_end\n";
			for(my $j=0; $j<=$i; $j++)
			{
				$sequence = $sequences[$j];
		    	my $length= $sequence->length;
		    	my @id = split(/:/,$sequence->id);
		    	my $gene = $id[0];
		    	my $chromosome = uc($id[3]);
				my $range_begin = $id[4];
		    	my $range_end = $id[5];
		    	my $strand = $id[6] =~ "-" ? "-" : "+";
		    	#$sequence->id=~m/(.*):chromosome:(\d*)\:(\d*)\-(\d*) +5'pad=(\d*) +3'pad=(\d*)/i;
		    	#my $chromosome=$1;
				#my $range_begin=$2;
		    	#my $range_end=$3;
		    	#my $pad_5=$4;
		    	#my $pad_3=$5;
		    	my ($roi_begin,$roi_end);
		    	#parse sequence
				if(!$opts{e})
			    {
			    	#include introns
			    	if ($pad_up + $pad_down >= $length){die $gene." is too short for padding of ".$pad_up.",".$pad_down.", stopped";}		    		
			        $roi_begin = $pad_up < 0 ? 1 : $pad_up+1;
			        $roi_end = $length - ($pad_down < 0 ? 0 : $pad_down);
			    	print PIECEMAKERFILE $gene."\t$gene\t$range_begin\t$roi_begin\t$roi_end\t".uc($sequence->seq)."\n"
			    		or die "Couldn't write to '$output_file', stopped";	    	
			    	print $gene."\t$gene\t$range_begin\t$roi_begin\t$roi_end\n";
		        }
			    else
			    {
			    	#exons only - exons must be in upper case
			    	#$sequence->seq =~ m/^(?<up>[a-z]+)?(?:(?<exon>[A-Z]+)(?<intron>[a-z]+(?=[A-Z]))?)*(?<downstream>[a-z]+)?$/;
		        }
			}		
		}
		print "Created '$output_file' file with " . ++$i . " sequences.\n";
	}
}