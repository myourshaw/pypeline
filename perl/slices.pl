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

# Obtain a slice covering the entire clone AL359765.6
$slice = $slice_adaptor->fetch_by_region( 'clone', 'AL359765.6' );

# Obtain a slice covering an entire NT contig
$slice = $slice_adaptor->fetch_by_region( 'supercontig', 'NT_011333' );

# Obtain a slice covering the region from 1MB to 2MB (inclusively) of
# chromosome 20
$slice = $slice_adaptor->fetch_by_region( 'chromosome', '20', 1e6, 2e6 );

#Obtain a slice that contains the sequence of the gene specified by its stable Ensembl ID
#It also returns 5000bp of flanking sequence at both the 5' and 3' ends
$slice = $slice_adaptor->fetch_by_gene_stable_id( 'ENSG00000099889', 5e3 );

# Retrieve slices of every chromosome in the database
my @slices = @{ $slice_adaptor->fetch_all('chromosome') };

# Retrieve slices of every BAC clone in the database
@slices = @{ $slice_adaptor->fetch_all('clone') };

#Break up larger slices into smaller component slices. 
use Bio::EnsEMBL::Utils::Slice qw(split_Slices);

my $slices = $slice_adaptor->fetch_all('chromosome');

# Base pair overlap between returned slices
my $overlap = 0;

# Maximum size of returned slices
my $max_size = 100_000;

# Break chromosomal slices into smaller 100k component slices
$slices = split_Slices( $slices, $max_size, $overlap );

#Obtain sequence from a slice 
my $sequence = $slice->seq();
print $sequence, "\n";

$sequence = $slice->subseq( 100, 200 );

#Query the slice for information about itself: 
# The method coord_system() returns a Bio::EnsEMBL::CoordSystem object
my $coord_sys  = $slice->coord_system()->name();
my $seq_region = $slice->seq_region_name();
my $start      = $slice->start();
my $end        = $slice->end();
my $strand     = $slice->strand();

print "Slice: $coord_sys $seq_region $start-$end ($strand)\n";

#Obtain a list of genes which overlap a slice: 
my $gene_adaptor  = $registry->get_adaptor( 'Human', 'Core', 'Gene' );

my @genes = @{ $gene_adaptor->fetch_all_by_Slice($slice) };

# Another way of doing the same thing:
@genes = @{ $slice->get_all_Genes() };
