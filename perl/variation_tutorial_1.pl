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

my $slice_adaptor = $registry->get_adaptor( 'zebrafish', 'core');
my $slice = $slice_adaptor->fetch_by_region('chromosome',25); #get chromosome 25 in zebrafish
print "!";
