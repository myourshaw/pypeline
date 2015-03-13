=head1 LICENSE

 Copyright (c) 2011-2012 Michael Yourshaw.  All rights reserved.                                                                      

=head1 CONTACT                                                                                                       

 Michael Yourshaw <myourshaw@ucla.edu>
    
=cut

=head1 NAME

 GeneIDs

=head1 SYNOPSIS

 mv GeneIDs.pm ~/.vep/Plugins
 perl variant_effect_predictor.pl -i variations.vcf --plugin GeneIDs

=head1 DESCRIPTION

 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that
 adds the following new columns:
 ENSG, Gene_Description, RefSeq_summary, Entrez_Gene_Name, UniProtKB_AC, UniProt_ID, Gene_Ontology.
 
 Requires vw plugin

=cut

package GeneIDs;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);
use Bio::EnsEMBL::Variation::Utils::VEP qw(@OUTPUT_COLS);

use vw;
use VAX qw(get_tva_info get_unique);

sub new {
    my $class = shift;
    my $self = $class->SUPER::new(@_);
    
    my $goa = $self->{config}->{reg}->get_adaptor('Multi', 'Ontology', 'GOTerm');
    $self->{goa} = $goa;

    return $self;
}

sub version {
    return '2.3';
}

sub feature_types {
    return ['Bio::EnsEMBL::Transcript', 'Bio::EnsEMBL::Funcgen::RegulatoryFeature', 'Bio::EnsEMBL::Funcgen::MotifFeature'];
}

sub get_header_info {
    my @new_output_cols = qw(
        ENSG
        Gene_Description
        RefSeq_summary
        Entrez_Gene_Name
        UniProt_ID
        Gene_Ontology
    );
    @OUTPUT_COLS = (@OUTPUT_COLS, @new_output_cols);

    return {
        ENSG => "Gene stable ID (Ensembl)",
        Gene_Description => "A short description of the gene (Ensembl)",
        RefSeq_summary => "Gene summary (RefSeq)",
        Entrez_Gene_Name => "The Entrez Gene name of the relevant gene (Ensembl)",
        UniProt_ID => "The UniProt ID of the relevant protein (Ensembl)",
        Gene_Ontology => "GO slim IDs and terms associated with the relevant transcript (Ensembl)",
    };
}

sub run {
    my ($self, $tva, $line_hash) = @_;

    my %tva_info = %{get_tva_info(@_)};

    my $gene = $tva_info{gene};;
    if (defined $gene){
        $line_hash->{ENSG} = $tva_info{ensg};

        if (defined $gene->{description}) {
            $line_hash->{Gene_Description} = $gene->{description}
        }
        
        my $xrefs = $gene->get_all_DBLinks(); 
        foreach my $xref (@{$xrefs}) {
            if ($xref->dbname() eq 'EntrezGene') {
                $line_hash->{Entrez_Gene_Name} = $xref->display_id();
#                $line_hash->{Entrez_Gene_ID} = $xref->primary_id();
                last;
            }
        }
        
        $xrefs = defined($tva_info{transcript}) ?
            $tva_info{transcript}->get_all_DBLinks() :
            $gene->get_all_DBLinks(); 
        my @uniprotkb_ac;
        my @uniprot_id;
        my @omim_id;
        my @hgmd_id;
        my @go_names = ();
        foreach my $xref (@{$xrefs}) {
            if ($xref->dbname() eq 'Uniprot/SWISSPROT') {
                push(
                    @uniprotkb_ac,
                    $xref->primary_id()
                );
                push(
                    @uniprot_id,
                    $xref->display_id()
                );
            }
            if ($xref->dbname() =~/^MIM/) {
                push(
                    @omim_id,
                    $xref->display_id()
                );
            }
            #if ($xref->dbname() =~ /HGMD/) {
            #    push(
            #        @hgmd_id,
            #        $xref->display_id()
            #    );
            #}
            elsif ($xref->dbname() eq 'goslim_goa') {
                my $go_id = $xref->display_id();
                my $go_term;
                my $go_name;
                my $go_definition;
                if (defined($self->{goa})) {
                    $go_term = $self->{goa}->fetch_by_accession($go_id);
                    $go_name = $go_term->name();
                    $go_definition = $go_term->definition();
                }
                if (defined($go_name)) {
                    push(@go_names, "$go_name [$go_id]");
                }
                else {
                    push(@go_names, "$go_id");
                }
            }
        }
        @uniprotkb_ac = @{get_unique(\@uniprotkb_ac)};
        $line_hash->{UniProtKB_AC} = $uniprotkb_ac[0];
        @uniprot_id = @{get_unique(\@uniprot_id)};
        $line_hash->{UniProt_ID} = join(',',@uniprot_id);
        #@omim_id = @{get_unique(\@omim_id)};
        #$line_hash->{OMIM_ID} = join(',',@omim_id);
        #@hgmd_id = @{get_unique(\@hgmd_id)};
        #$line_hash->{HGMD_ID} = join(',',@hgmd_id);
        my $gene_ontology = join('|', @{get_unique(\@go_names)});
        $line_hash->{Gene_Ontology} = $gene_ontology;

        my @data;
        my $query = "CALL $vw::vw_database.ensg2refseq_summary('$tva_info{ensg}')";
        my $qh = $vw::vw_conn->prepare($query);
        $qh->execute() or die "Unable to execute $query: $DBI::errstr\n";
        while (my @row = $qh->fetchrow_array()){
            push @data, $row[0];
        }
        $line_hash->{RefSeq_summary} = @data ? join('|', @data) : '';
    }
    return {};
}


1;

