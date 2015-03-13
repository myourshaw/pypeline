#!/usr/bin/perl -w
#extracts SNP data from Entrez SNP flat file
#http://www.ncbi.nlm.nih.gov/books/bv.fcgi?rid=helpsnpfaq.section.Search.Locating_SNPs_in_a_C
use strict;
use warnings;
my $id;
my $type;
my $alleles;
my $chr;
my $pos;
print "snpId\tsnpType\talleles\tchr\tposition\n";
while(<>){
    if(/(\S+) \| Homo sapiens \| 9606 \| (\S+) \| /){
        $id = $1;
        $type = $2;
    }
    if(/^SNP \| alleles='(\S+)' \| /){
        $alleles = $1;
    }
    if(/^CTG \| assembly=reference \| chr=([0-9XYMT]{1,2}) \| chr-pos=(\d+) \| /){
        $chr = $1;
        $pos = $2;
        if(defined($id)){
            print "$id\t$type\t$alleles\t$chr\t$pos\n";
            undef($id);
        }
     }
}
