#!/usr/bin/perl -w
use strict;

use warnings;

use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
);

#>ACSL4:chromosome:ncbi36:x:108770622:108863877:-1
my $chr = "X";
my $start = 108770622;
my $end = 108863877;
my $strand = -1;


my $slice_adaptor = $registry->get_DBAdaptor( 'Human', 'Core');
my $var_adaptor = $registry->get_DBAdaptor( 'Human', 'variation' );

my $slice2 = sequence_with_ambiguity($slice_adaptor,$var_adaptor,$chr,$start,$end,$strand);
my $sequence = $slice2->seq();
print length $sequence."\n";
print $sequence;

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

#my $va_adaptor = $dbVar->get_VariationAdaptor; #get the different adaptors for the different objects needed
#my $vf_adaptor = $dbVar->get_VariationFeatureAdaptor;
#my $gene_adaptor = $dbCore->get_GeneAdaptor;  

#sequence_with_ambiguity

#  Arg[1]      : Bio::EnsEMBL::DBSQL::DBAdaptor $dbCore
#  Arg[2]      : Bio::EnsEMBL::Variation::DBSQL::DBAdaptor $dbVar
#  Arg[3]      : string $chr (chromosome)
#  Arg[4]      : int $start
#  Arg[5]      : int $end
#  Arg[6]      : int $strand

use Bio::EnsEMBL::Variation::Utils::Sequence qw(sequence_with_ambiguity);

my $slice = sequence_with_ambiguity($slice_adaptor,$var_adaptor,$chr,$start,$end,$strand);
my $seq = $slice->seq();
print $seq;
print "\n";
print length $seq;
print "\n";
print "the sequence with ambiguity code for your region is: ",$slice->seq();
print $slice->seq()->length;


