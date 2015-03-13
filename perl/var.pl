#!/usr/bin/perl -w
use strict;

use warnings;

use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
);

#>ACSL4:chromosome:ncbi36:x:108770622:108863877:-1
my $chr = "X";
my $start = 31042729-1000;
my $end = 37522716+1000;
my $strand = -1;

my $slice_adaptor = $registry->get_adaptor( 'Human', 'Core', 'slice');
my $var_adaptor = $registry->get_DBAdaptor( 'Human', 'variation' );

my $slice = $slice_adaptor->fetch_by_region('chromosome',"X",$start,$end,$strand); #get chromosome slice

use Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor;

my $variation_feature_adaptor = $var_adaptor->get_VariationFeatureAdaptor();

 my $variation_features = $variation_feature_adaptor->fetch_all_by_Slice($slice);

foreach my $vf (@{$variation_features}){
     #print "Variation: ", $vf->variation_name, " with alleles ", $vf->allele_string, " in chromosome ", $slice->seq_region_name, " and position ", $vf->start,"-",$vf->end,"\n";
     if($strand==1){
          #??
          printf( "%s\t%s\t%s\t%u\t%u\n",$vf->variation_name,$vf->allele_string,$slice->seq_region_name,$end-$vf->start+1,$end-$vf->end) #chromStart,chromEnd coord between nucleotides (1st = 0,1);
     }
     else{
          printf( "%s\t%s\t%s\t%u\t%u\n",$vf->variation_name,$vf->allele_string,$slice->seq_region_name,$end-$vf->end,$end-$vf->start+1) #chromStart,chromEnd coord between nucleotides (1st = 0,1);
     }
}

1;
