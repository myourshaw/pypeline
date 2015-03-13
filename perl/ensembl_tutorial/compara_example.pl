#!/usr/bin/perl -w
use strict;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::SimpleAlign;
use Bio::AlignIO;
use Bio::LocatableSeq;
use Getopt::Long;

my $usage = qq{
perl DumpMultiAlign.pl
  Getting help:
    [--help]

  For the query slice:
    [--species species]
        Query species. Default is "human"
    [--coord_system coordinates_name]
        Query coordinate system. Default is "chromosome"
    --seq_region region_name
        Query region name, i.e. the chromosome name
    --seq_region_start start
    --seq_region_end end

  For the alignments:
    [--method_link_type method_link_name]
        The type of alignment. Default is "GERP_CONSERVATION_SCORE"
    [--species_set species_set]
        The species set used to get those alignments. Default is
        "mammal". The names should correspond to the name of the
        core database in the registry_configuration_file or any of its
        aliases

  Ouput:
    [--output_format clustalw|fasta|...]
        The type of output you want. "clustalw" is the default.
    [--output_file filename]
        The name of the output file. By default the output is the
        standard output
};


my $method_link_type = 'GERP_CONSERVATION_SCORE';
my $species_set  = 'mammals';
my $species = "human";
my $coord_system = "chromosome";
my $seq_region = "14";
my $seq_region_start = 75000000;
my $seq_region_end = 75100000;
my $output_file = undef;
my $output_format = "clustalw";
my $help;

GetOptions(
    "help" => \$help,
    "species=s" => \$species,
    "coord_system=s" => \$coord_system,
    "seq_region=s" => \$seq_region,
    "seq_region_start=i" => \$seq_region_start,
    "seq_region_end=i" => \$seq_region_end,
    "method_link_type=s" => \$method_link_type,
    "species_set=s" => \$species_set,
    "output_format=s" => \$output_format,
    "output_file=s" => \$output_file);

# Print Help and exit
if ($help) {
    print $usage;
    exit(0);
}

if ($output_file) {
    open(STDOUT, ">$output_file") or die("Cannot open $output_file");
}

Bio::EnsEMBL::Registry->load_registry_from_db(
    -host => 'cortex.local', -user => 'ensembl', -pass => 'ensembl', -port => 3306);


my $reg = 'Bio::EnsEMBL::Registry';
my $mlss_adap = $reg->get_adaptor('Multi', 'compara', 'MethodLinkSpeciesSet')
    or die "Failed to connect to compara database\n";

my $mlss = $mlss_adap->fetch_by_method_link_type_species_set_name($method_link_type, $species_set)
    or die "Failed to fetch MLSS for ".$method_link_type." and ".$species_set."\n";

my $cs_adap = $reg->get_adaptor('Multi', 'compara', 'ConservationScore')
    or die "Failed to fetch conservation adaptor\n";

## Getting all the Bio::EnsEMBL::Compara::GenomeDB objects
#my $genome_dbs;
#my $genome_db_adaptor = Bio::EnsEMBL::Registry->get_adaptor(
#    'Multi', 'compara', 'GenomeDB');
#
#throw("Cannot connect to Compara") if (!$genome_db_adaptor);
#
#foreach my $this_species (split(":", $set_of_species)) {
#    my $this_meta_container_adaptor = Bio::EnsEMBL::Registry->get_adaptor(
#        $this_species, 'core', 'MetaContainer');
#
#    throw("Registry configuration file has no data for connecting to <$this_species&")
#        if (!$this_meta_container_adaptor);
#
#    my $this_production_name = $this_meta_container_adaptor->get_production_name;
#
#    # Fetch Bio::EnsEMBL::Compara::GenomeDB object
#    my $genome_db = $genome_db_adaptor->fetch_by_name_assembly($this_production_name);
#
#    # Add Bio::EnsEMBL::Compara::GenomeDB object to the list
#    push(@$genome_dbs, $genome_db);
#}
#
## Getting Bio::EnsEMBL::Compara::MethodLinkSpeciesSet object
#my $method_link_species_set_adaptor = Bio::EnsEMBL::Registry->get_adaptor(
#    'Multi', 'compara', 'MethodLinkSpeciesSet');
#
#my $method_link_species_set =
#    $method_link_species_set_adaptor->fetch_by_method_link_type_GenomeDBs(
#      $alignment_type, 
#      $genome_dbs);
#
#throw("The database do not contain any $alignment_type data for $set_of_species!")
#    if (!$method_link_species_set);

# Fetching the query Slice:
my $slice_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species, 'core', 'Slice');

throw("Registry configuration file has no data for connecting to <$species>")
    if (!$slice_adaptor);

my $slice = $slice_adaptor->fetch_by_region('toplevel', $seq_region, $seq_region_start, $seq_region_end);

throw("No Slice can be created with coordinates $seq_region:$seq_region_start-".
    "$seq_region_end") if (!$slice);

## Fetching all the GenomicAlignBlock corresponding to this Slice:
#my $genomic_align_block_adaptor = Bio::EnsEMBL::Registry->get_adaptor(
#    'Multi', 'compara', 'GenomicAlignBlock');

my $conservation_scores = $cs_adap->fetch_all_by_MethodLinkSpeciesSet_Slice($mlss, $slice);

my $scores = $cs_adap->fetch_all_by_MethodLinkSpeciesSet_Slice(
    $mlss,                      # our MLSS for the conservation metric and the set of species
    $slice,                             # our slice
    ($slice->end - $slice->start + 1),  # the number of scores we want back (one for each base)
);


#my $all_aligns;
#
## Get a Bio::SimpleAlign object from every GenomicAlignBlock
#foreach my $this_genomic_align_block (@$genomic_align_blocks) {
#    my $simple_align = $this_genomic_align_block->get_SimpleAlign;
#    push(@$all_aligns, $simple_align);
#}
#
## print all the genomic alignments using a Bio::AlignIO object
#my $alignIO = Bio::AlignIO->newFh(
#    -interleaved => 0,
#    -fh => \*STDOUT,
#    -format => $output_format,
#    -idlength => 10
#);
#  
#foreach my $this_align (@$all_aligns) {
#    print $alignIO $this_align;
#}

exit;


