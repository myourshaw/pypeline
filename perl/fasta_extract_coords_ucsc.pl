#!/usr/bin/perl -w
use strict;
use warnings;

use Bio::Seq;
use Bio::SeqIO;

my $file = shift;
print "id\tregion\tchrom\tstart\tend\tstrand\tlength\n";
my $seqio_obj = Bio::SeqIO->new(-file => $file, -format => "fasta" );
my $seq_obj;
while($seq_obj = $seqio_obj->next_seq){
	$seq_obj->primary_id =~ /hg\d+_\S+_\S+_(\d+)/;
	my $no = int($1/2+1);
	$no .= $1%2?"i":"e";
	$seq_obj->desc =~/range=(\S+):(\d+)-(\d+).*strand=(\S)/;
	printf "%s\t%s\t%s\t%u\t%u\t%s\t%u\n", $seq_obj->primary_id,$no,$1,$2,$3,$4,$3-$2+1;
}
