#!/usr/bin/perl -w
use strict;
use warnings;

use Bio::EnsEMBL::Registry;

my $reg = 'Bio::EnsEMBL::Registry';

$reg->load_registry_from_db(
    -host => 'cortex.local',
    -port => 3306,
    -user => 'ensembl',
    -pass => 'ensembl',
);

my $enst = 'ENST00000222270'; 
my $ga    = $reg->get_adaptor( 'Homo_sapiens', 'Core', 'Gene' );
my $ta    = $reg->get_adaptor( 'Homo_sapiens', 'Core', 'Transcript' );
my $va = $reg->get_adaptor('homo_sapiens', 'variation', 'variation');
my $vfa = $reg->get_adaptor("human","variation","variationfeature");
my $sa = $reg->get_adaptor("human","core","slice");

my $slice = $sa->fetch_by_region('chromosome', '19', 36211374, 36211374);

my $hgnc;

foreach my $gene ( @{ $slice->get_all_Genes } ) {
    my $ensg = $gene->stable_id();#ENSG00000105663
    my @dbentries = @{$gene->get_all_DBEntries()};
    foreach my $dbentry(@dbentries){
        print "$dbentry->{dbname}\t$dbentry->{primary_id}\t$dbentry->{display_id}\n";
    }
    my @entries = grep {$_->database eq 'HGNC'} @{$gene->get_all_DBEntries()};
    if(scalar @entries) {
        $hgnc = $entries[0]->display_id;
    }
    $hgnc = undef if defined($hgnc) && $hgnc eq '-';
 }




  foreach my $vf (@{$vfa->fetch_all_by_Slice($slice)}) {
    print $vf->start(), '-', $vf->end(), ' ', $vf->allele_string(), "\n";
  }



my $t = $ta->fetch_by_stable_id($enst);
my $genes = $t->get_overlapping_Genes();
while ( my $gene = shift @{$genes} ) {
    my $ensg = $gene->stable_id();#ENSG00000105663
    my $hgnc = $t->{_gene_hgnc};
    if(!defined($hgnc && $gene)) {
        my @entries = grep {$_->database eq 'HGNC'} @{$gene->get_all_DBEntries()};
        if(scalar @entries) {
            $hgnc = $entries[0]->display_id;
        }
    }
    $hgnc = undef if defined($hgnc) && $hgnc eq '-';
}



