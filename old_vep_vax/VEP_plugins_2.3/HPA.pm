=head1 LICENSE

 Copyright (c) 2011-2011 Michael Yourshaw.  All rights reserved.                                                                      

=head1 CONTACT                                                                                                       

 Michael Yourshaw <myourshaw@ucla.edu>
    
=cut

=head1 NAME

 HPA

=head1 SYNOPSIS

 mv HPA.pm ~/.vep/Plugins
 perl variant_effect_predictor.pl -i variations.vcf --plugin HPA,cortex.local,3306,vw,vw,mysql,vw

=head1 DESCRIPTION

 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that
 adds the following new columns:
 HPA_<tissue>_<cell_type>, ..., HPA_subcellular_location.
 All fields are uri escaped for ';='.
 
 Requires vw plugin

=cut

package HPA;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);
use Bio::EnsEMBL::Variation::Utils::VEP qw(@OUTPUT_COLS);

use vw;
use VAX qw(get_tva_info);

sub new {
    my $class = shift;
    my $self = $class->SUPER::new(@_);

    return $self;
}

sub version {
    return '2.3';
}

sub feature_types {
    return ['Bio::EnsEMBL::Transcript', 'Bio::EnsEMBL::Funcgen::RegulatoryFeature', 'Bio::EnsEMBL::Funcgen::MotifFeature'];
}

sub get_header_info {
    my ($self) = @_;

    my %tissue_cell_type_columns;
    $tissue_cell_type_columns{HPA_subcellular_location} = "Subcellular localisation of proteins based on immunofluorescently stained cells (Human Protein Atlas)";
    my @new_output_cols = ();
    my $query = "CALL $vw::vw_database.hpa_tissue_cell_type_list()";
    my $qh = $vw::vw_conn->prepare($query);
    $qh->execute() or die "Unable to execute $query: $DBI::errstr\n";
    while (my @row = $qh->fetchrow_array()){
        push @new_output_cols, 'HPA_'.unspace($row[0]).'_'.unspace($row[1]);
        $tissue_cell_type_columns{'HPA_'.unspace($row[0]).'_'.unspace($row[1])} = "Expression profile for protein in human $row[0] $row[1] based on immunohistochemisty using tissue micro arrays (Human Protein Atlas)";
    }
    @new_output_cols = sort(@new_output_cols);
    push @new_output_cols, 'HPA_subcellular_location';
    @OUTPUT_COLS = (@OUTPUT_COLS, @new_output_cols);

    return \%tissue_cell_type_columns;
}

sub run {
    my ($self, $tva, $line_hash) = @_;

    my %tva_info = %{get_tva_info(@_)};

    if (defined $tva_info{gene}){
        my @data;
        
        my $query = "CALL $vw::vw_database.ensg2hpa_subcellular_location('$tva_info{ensg}')";
        my $qh = $vw::vw_conn->prepare($query);
        $qh->execute() or die "Unable to execute $query: $DBI::errstr\n";
        while (my @row = $qh->fetchrow_array()){
            push @data, $row[0];
        }
        $line_hash->{HPA_subcellular_location} = @data ? join('|', @data) : '';

        $query = "CALL $vw::vw_database.ensg2hpa_tissue('$tva_info{ensg}')";
        $qh = $vw::vw_conn->prepare($query);
        $qh->execute() or die "Unable to execute $query: $DBI::errstr\n";
        while (my @row = $qh->fetchrow_array()){
            $line_hash->{'HPA_'.unspace($row[0]).'_'.unspace($row[1])} = $row[2];
        }
    }
    return {};
}

sub unspace{
    my $s = shift;
    $s =~ s/, |[\s"]/_/g;
    return $s;
}


1;

