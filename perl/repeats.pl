#!/usr/bin/perl -w
use strict;

use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
);

my $slice_adaptor = $registry->get_adaptor( 'Human', 'Core', 'Slice' );

# Obtain a slice covering the entire chromosome X
my $slice = $slice_adaptor->fetch_by_region( 'chromosome', 'X' );

my @repeats = @{ $slice->get_all_RepeatFeatures() };

foreach my $repeat (@repeats) {
    printf( "%s %d-%d\n",
        $repeat->display_id(), $repeat->start(), $repeat->end() );
}

#Hard-masking replaces sequence in repeat regions with Ns
#Soft-masking replaces sequence in repeat regions with lower-case sequence
my $unmasked_seq   = $slice->seq();
my $hardmasked_seq = $slice->get_repeatmasked_seq();
my $softmasked_seq = $slice->get_repeatmasked_seq( undef, 1 );

# Soft-mask sequence using TRF results only
my $tandem_masked_seq = $slice->get_repeatmasked_seq( ['TRF'], 1 );
