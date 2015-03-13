#!/usr/bin/perl -w
use strict;
use warnings;

use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

# connect to Variation database
my $dbVariation = Bio::EnsEMBL::Variation::DBSQL::DBAdaptor->new
  (-host   => 'ensembldb.ensembl.org',
   -dbname => 'gallus_gallus_variation_47_2e',
   -species => 'chicken',
   -group   => 'variation',
   -user   => 'anonymous');

# connect to Core database
my $dbCore = Bio::EnsEMBL::DBSQL::DBAdaptor->new
  (-host   => 'ensembldb.ensembl.org',
   -dbname => 'gallus_gallus_core_47_2e',
   -species => 'chicken',
   -group   => 'core',
   -user   => 'anonymous');

my $stable_id = 'ENSGALT00000007843'; #this is the stable_id of a chicken transcript
my $transcript_adaptor = $dbCore->get_TranscriptAdaptor(); #get the adaptor to get the Transcript from the database
my $transcript = $transcript_adaptor->fetch_by_stable_id($stable_id); #get the Transcript object
my $trv_adaptor = $dbVariation->get_TranscriptVariationAdaptor; #get the adaptor to get TranscriptVariation objects

my $trvs = $trv_adaptor->fetch_all_by_Transcripts([$transcript]); #get ALL effects of Variations in the Transcript

foreach my $tv (@{$trvs}){
    print "SNP :",$tv->variation_feature->variation_name, " has a consequence/s ", join(",",@{$tv->consequence_type}), " in transcript ", $stable_id, "\n";
    #print the name of the variation and the effect (consequence_type) of the variation in the Transcript
}
exit 0;


