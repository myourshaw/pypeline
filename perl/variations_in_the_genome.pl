#!/usr/bin/perl -w
use strict;
use warnings;

use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

# connect to Variation database
my $dbVariation = Bio::EnsEMBL::Variation::DBSQL::DBAdaptor->new
  (-host   => 'ensembldb.ensembl.org',
   -dbname => 'danio_rerio_variation_47_7a',
   -species => 'zebrafish',
   -group   => 'variation',
   -user   => 'anonymous');

# connect to Core database
my $dbCore = Bio::EnsEMBL::DBSQL::DBAdaptor->new
  (-host   => 'ensembldb.ensembl.org',
   -dbname => 'danio_rerio_core_47_7a',
   -species => 'zebrafish',
   -group   => 'core',
   -user   => 'anonymous');


my $slice_adaptor = $dbCore->get_SliceAdaptor(); #get the database adaptor for Slice objects
my $slice = $slice_adaptor->fetch_by_region('chromosome',25); #get chromosome 25 in zebrafish
my $vf_adaptor = $dbVariation->get_VariationFeatureAdaptor(); #get adaptor to VariationFeature object
my $vfs = $vf_adaptor->fetch_all_by_Slice($slice); #return ALL variations defined in $slice
foreach my $vf (@{$vfs}){
    print "Variation: ", $vf->variation_name, " with alleles ", $vf->allele_string, " in chromosome ", $slice->seq_region_name, " and position ", $vf->start,"-",$vf->end,"\n";
}
exit 0;


