#!/usr/bin/perl -w
use strict;
use warnings;

#BioPerl
use Bio::Seq;
use Bio::SeqIO;
#EnsEMBL Perl API
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::Utils::Sequence qw (sequence_with_ambiguity);
#connect to EnsEMBL database
my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db(
	-host => 'ensembldb.ensembl.org',
	-user => 'anonymous');
my $exon_adaptor = $registry->get_adaptor('human', 'core', 'exon');
my $gene_adaptor = $registry->get_adaptor('human', 'core', 'gene');
my $repeat_adaptor = $registry->get_adaptor('human', 'core', 'repeatfeature');
my $slice_adaptor = $registry->get_adaptor('human', 'core', 'slice');
my $transcript_adaptor = $registry->get_adaptor('human', 'core', 'transcript');
my $id = 'MYO1D';
my @genes = @{ $gene_adaptor->fetch_all_by_external_name($id) };
foreach my $gene(@genes){
    my @transcripts = @{ $gene->get_all_Transcripts };
    my @exons = @{ $gene->get_all_Exons };
    print @transcripts,@genes;
 }
