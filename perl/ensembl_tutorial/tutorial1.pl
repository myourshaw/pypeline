#!/usr/bin/perl -w
use strict;

use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'cortex.local',
    -user => 'ensembl',
    -pass => 'ensembl',
    -port => 3306,
    #-host => 'myourshaw-dev.genome.ucla.edu',
    #-user => 'ensembl',
    #-pass => 'ensembl'
);
#
#my @db_adaptors = @{ $registry->get_all_DBAdaptors() };
#
#foreach my $db_adaptor (@db_adaptors) {
#    my $db_connection = $db_adaptor->dbc();
#
#    printf(
#        "species/group\t%s/%s\ndatabase\t%s\nhost:port\t%s:%s\n\n",
#        $db_adaptor->species(),   $db_adaptor->group(),
#        $db_connection->dbname(), $db_connection->host(),
#        $db_connection->port()
#    );
#}

#my $genome_db_adaptor = Bio::EnsEMBL::Registry->get_adaptor(
#    'Multi', 'compara', 'GenomeDB');

my $stable_id = 'ENST00000528762';

my $transcript_adaptor =
  $registry->get_adaptor( 'Human', 'Core', 'Transcript' );
my $transcript = $transcript_adaptor->fetch_by_stable_id($stable_id);

my $translation = $transcript->translation();

print $transcript->translation()->stable_id(), "\n";
my $protein_seq = $translation->seq();
print $protein_seq,         "\n";


my $pfeatures = $translation->get_all_ProteinFeatures();
my $seg_features    = $translation->get_all_ProteinFeatures('Seg');
my $domain_features = $translation->get_all_DomainFeatures();
while ( my $pfeature = shift @{$pfeatures} ) {
    my $logic_name = $pfeature->analysis()->logic_name();

    printf(
        "%d-%d %s %s %s\n",
        $pfeature->start(), $pfeature->end(), $logic_name,
        $pfeature->interpro_ac(),
        $pfeature->idesc()
    );
}

print 'done';
