#!/usr/bin/perl -w
#"C:\Users\Michael\Documents\Lab\DMD\NM_004006_ENSG00000198947.txt"
use strict;
use warnings;
use Bio::Seq;
use Bio::SeqIO;

use constant ID => "ENSG00000198947";
my $chr;
my $start;
my $end;
my $strand;
my $seq = "";
my @seqlines;
my $fragment_number = 0;
my $head_in = '>chromosome:NCBI\d+:(d{1,2}|X|Y|M):(\d+):(\d+):(-?1)';
my $display_idf = "%s_%u";
my $descf = "%s:%u-%u:%s";
my $seqio_obj = Bio::SeqIO->new(-file => '>sequence.fasta', -format => 'fasta' );
while(<>){
	chomp;
	if(/$head_in/){
		($chr,$start,$end,$strand) = ($1,$2,$3,$4==1?"+":"-");
		if($seq){
			out();
			$seq = "";
			$fragment_number = 0;
		}
	}
	elsif(/^[A-Za-z]+$/){
		$seq.=$_;
	}
}
if($seq){
	out();
}
sub out{
if($strand eq "-"){
	#$seq = reverse($seq);
}
$seq =~ s/([A-Z])([a-z])/$1\n$2/g;
$seq =~ s/([a-z])([A-Z])/$1\n$2/g;
@seqlines = split "\n",$seq;
foreach my $s(@seqlines){
	if($strand eq "-"){
		$start = $end-length($s)+1;
	}
	else{
		$end = $start+length($s)-1;
	}
	my $display_id = sprintf($display_idf,ID,$fragment_number++);
	my $desc = sprintf($descf,$chr,$start,$end,$strand);
	my $seq_obj = Bio::Seq->new(-display_id => $display_id,-desc => $desc,-seq => $s,-alphabet => 'dna' );
	$seqio_obj->write_seq($seq_obj);
	if($strand eq "-"){
		$end = $start-1;
	}
	else{
		$start = $end+1;
	}

}

}
