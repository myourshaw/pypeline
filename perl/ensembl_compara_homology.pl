#!/usr/bin/perl -w

use strict;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::AlignIO;
use Getopt::Long;

Bio::EnsEMBL::Registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org', -user => 'anonymous');

# first you have to get a GeneMember object. In case of homology is a gene, in 
# case of family it can be a gene or a protein

# get the MemberAdaptor
my $genemember_adaptor = Bio::EnsEMBL::Registry->get_adaptor('Multi','compara','GeneMember');

# fetch a Member
my $gene_member = $genemember_adaptor->fetch_by_stable_id('ENSG00000004059');

# print out some information about the Member
print $gene_member->source_name, ": ", $gene_member->dnafrag->name, " ( ", $gene_member->dnafrag_start, " - ", $gene_member->dnafrag_end, " ): ", $gene_member->description, "\n";

my $taxon = $gene_member->taxon;
print "common_name ", $taxon->common_name,"\n";
print "genus ", $taxon->genus,"\n";
print "species ", $taxon->species,"\n";
print "binomial ", $taxon->binomial,"\n";
print "classification ", $taxon->classification,"\n";

# then you get the homologies where the member is involved

my $homology_adaptor = Bio::EnsEMBL::Registry->get_adaptor('Multi', 'compara', 'Homology');
my $homologies = $homology_adaptor->fetch_all_by_Member($gene_member);

# That will return a reference to an array with all homologies (orthologues in
# other species and paralogues in the same one)
# Then for each homology, you can get all the Members implicated

foreach my $homology (@{$homologies}) {
  # You will find different kind of description
  # see ensembl-compara/docs/docs/schema_doc.html for more details

  print $homology->description," ", $homology->taxonomy_level,"\n";

  # And if they are defined dN and dS related values

  print " dn ", $homology->dn // '',"\n";
  print " ds ", $homology->ds // '',"\n";
  print " dnds_ratio ", $homology->dnds_ratio // '',"\n";
}

my $homology = $homologies->[0]; # take one of the homologies and look into it

foreach my $member (@{$homology->get_all_Members}) {

  # each AlignedMember contains both the information on the SeqMember and in
  # relation to the homology

  print join(" ", map { $member->$_ } qw(stable_id taxon_id)), "\n";
  print join(" ", map { $member->$_ } qw(perc_id perc_pos perc_cov)), "\n";

}

use Bio::AlignIO;

my $simple_align = $homology->get_SimpleAlign();
my $alignIO = Bio::AlignIO->newFh(
    -interleaved => 0,
    -fh => \*STDOUT,
    -format => "clustalw",
    -idlength => 20);

print $alignIO $simple_align;


my $family_adaptor = Bio::EnsEMBL::Registry->get_adaptor('Multi','compara','Family');
my $families = $family_adaptor->fetch_all_by_GeneMember($gene_member);

foreach my $family (@{$families}) {
    print join(" ", map { $family->$_ }  qw(description description_score))."\n";

    foreach my $member (@{$family->get_all_Members}) {
        print $member->stable_id," ",$member->taxon_id,"\n";
    }

    my $simple_align = $family->get_SimpleAlign();
    my $alignIO = Bio::AlignIO->newFh(
        -interleaved => 0,
        -fh          => \*STDOUT,
        -format      => "phylip",
        -idlength    => 20);

    print $alignIO $simple_align;

    $simple_align = $family->get_SimpleAlign(-seq_type => 'cds');
    $alignIO = Bio::AlignIO->newFh(
        -interleaved => 0,
        -fh          => \*STDOUT,
        -format      => "phylip",
        -idlength    => 20);

    print $alignIO $simple_align;
}
