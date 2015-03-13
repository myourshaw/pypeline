#!/usr/bin/perl -w
use strict;

use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
);

my $slice_adaptor = $registry->get_adaptor( 'Human', 'Core', 'Slice' );
my $tr_adaptor    = $registry->get_adaptor( 'Human', 'Core', 'Transcript' );
my $daf_adaptor   = $registry->get_adaptor( 'Human', 'Core', 'DnaAlignFeature' );

# Get a slice of chromosome 20, 10MB-11MB
my $slice = $slice_adaptor->fetch_by_region( 'chromosome', '20', 10e6, 11e6 );

# Fetch all of the transcripts overlapping chromosome 20, 10MB-11MB
my @transcripts = @{ $tr_adaptor->fetch_all_by_Slice($slice) };
while ( my $tr = shift @transcripts ) {
    my $dbID      = $tr->dbID();
    my $start     = $tr->start();
    my $end       = $tr->end();
    my $strand    = $tr->strand();
    my $stable_id = $tr->stable_id();

    print "Transcript $stable_id [$dbID] $start-$end ($strand)\n";
}

# Fetch all of the DNA-DNA alignments overlapping chromosome 20, 10MB-11MB
my @dafs = @{ $daf_adaptor->fetch_all_by_Slice($slice) };
while ( my $daf = shift @dafs ) {
    my $dbID     = $daf->dbID();
    my $start    = $daf->start();
    my $end      = $daf->end();
    my $strand   = $daf->strand();
    my $hseqname = $daf->hseqname();

    print "DNA Alignment $hseqname [$dbID] $start-$end ($strand)\n";
}

# Fetch a transcript by its internal identifier
my $transcript = $tr_adaptor->fetch_by_dbID(100);

# Fetch a dnaAlignFeature by its internal identifiers
my $dna_align_feat = $daf_adaptor->fetch_by_dbID(100);

