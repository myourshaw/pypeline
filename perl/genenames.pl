#!/usr/bin/perl -w
use strict;
use warnings;
use Bio::SeqIO;
my $file = shift;
my %genes;
my $seqin = Bio::SeqIO->new(-file => $file,
                                       -format => 'fasta');
while(my $seq = $seqin->next_seq){
    $seq->desc =~ /^\|.*?\|.*?\|.*?\|([^|]+)/;
    if(!exists($genes{$1})){

        $genes{$1} = $seq->description;
        print "$1\n";
    }
}
