#!/usr/bin/perl -w
use strict;

use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
);

#Retrieve and print DBEntries for a gene, its transcripts and its translations: 

# Define a helper subroutine to print DBEntries
sub print_DBEntries
{
    my $db_entries = shift;

    foreach my $dbe ( @{$db_entries} ) {
        printf "\tXREF %s (%s)\n", $dbe->display_id(), $dbe->dbname();
    }
}

my $gene_adaptor = $registry->get_adaptor( 'Human', 'Core', 'Gene' );

# Get the 'COG6' gene from human
my $gene = $gene_adaptor->fetch_by_display_label('COG6');

print "GENE ", $gene->stable_id(), "\n";
print_DBEntries( $gene->get_all_DBEntries() );

foreach my $transcript ( @{ $gene->get_all_Transcripts() } ) {
    print "TRANSCRIPT ", $transcript->stable_id(), "\n";
    print_DBEntries( $transcript->get_all_DBEntries() );

    # Watch out: pseudogenes have no translation
    if ( defined $transcript->translation() ) {
        my $translation = $transcript->translation();

        print "TRANSLATION ", $translation->stable_id(), "\n";
        print_DBEntries( $translation->get_all_DBEntries() );
    }
}

#Shorten the above using
print_DBEntries( $gene->get_all_DBLinks() );
