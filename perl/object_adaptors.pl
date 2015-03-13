#!/usr/bin/perl -w
use strict;

use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
);

my $gene_adaptor  = $registry->get_adaptor( 'Human', 'Core', 'Gene' );
my $slice_adaptor = $registry->get_adaptor( 'Human', 'Core', 'Slice' );

