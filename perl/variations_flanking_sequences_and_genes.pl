#!/usr/bin/perl -w
use strict;
use warnings;

use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

# connect to Variation database
my $dbVar = Bio::EnsEMBL::Variation::DBSQL::DBAdaptor->new
  (-host   => 'ensembldb.ensembl.org',
   -dbname => 'homo_sapiens_variation_47_36i',
   -species => 'human',
   -group   => 'variation',
   -user   => 'anonymous');

# connect to Core database
my $dbCore = Bio::EnsEMBL::DBSQL::DBAdaptor->new
  (-host   => 'ensembldb.ensembl.org',
   -dbname => 'homo_sapiens_core_47_36i',
   -species => 'human',
   -group   => 'core',
   -user   => 'anonymous');

my $va_adaptor = $dbVar->get_VariationAdaptor; #get the different adaptors for the different objects needed
my $vf_adaptor = $dbVar->get_VariationFeatureAdaptor;
my $gene_adaptor = $dbCore->get_GeneAdaptor;  
my @rsIds = qw(rs1367827 rs1367830);
foreach (@rsIds){
# get Variation object
  my $var = $va_adaptor->fetch_by_name($_); #get the Variation from the database using the name
  &get_VariationFeatures($var);
}


sub get_VariationFeatures{
	my $var = shift;
        # get all VariationFeature objects: might be more than 1 !!!
	foreach my $vf (@{$vf_adaptor->fetch_all_by_Variation($var)}){
	    print $vf->variation_name(),","; # print rsID
            print $vf->allele_string(),","; # print alleles
            print join(",",@{$vf->get_consequence_type()}),","; # print consequenceType
	    print substr($var->five_prime_flanking_seq,-10) , "[",$vf->allele_string,"]"; #print the allele string
	    print substr($var->three_prime_flanking_seq,0,10), ","; # print RefSeq
	    print $vf->seq_region_name, ":", $vf->start,"-",$vf->end; # print position in Ref in format Chr:start-end
            &get_TranscriptVariations($vf); # get Transcript information
	}
}

sub get_TranscriptVariations{
	my $vf = shift;
        # get all TranscriptVariation objects: might be more than 1 !!!
	my $transcript_variations = $vf->get_all_TranscriptVariations; #get ALL the effects of the variation in different Transcripts
	if (defined $transcript_variations){
	    foreach my $tv (@{$transcript_variations}){
		print ",", $tv->pep_allele_string if (defined $tv->pep_allele_string);# the AA change, but only if it is in a coding region
		my $gene = $gene_adaptor->fetch_by_transcript_id($tv->transcript->dbID);
		print ",",$gene->stable_id if (defined $gene->external_name); # and the external gene name	   

	    }
	}
	print "\n";
}
exit 0;


