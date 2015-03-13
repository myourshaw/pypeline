=head1 LICENSE

 Copyright (c) 2011-2012 Michael Yourshaw.  All rights reserved.                                                                      

=head1 CONTACT                                                                                                       

 Michael Yourshaw <myourshaw@ucla.edu>
    
=cut

=head1 NAME

 ExtraCols

=head1 SYNOPSIS

 mv ExtraCols.pm ~/.vep/Plugins
 perl variant_effect_predictor.pl -i variations.vcf --plugin ExtraCols

=head1 DESCRIPTION

 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that
 adds the VEP's Extra data to additional output columns and
 expands the SIFT, PolyPhen, and Condel entries to prediction and score:
 SIFT_prediction, SIFT_score, PolyPhen_prediction, PolyPhen_score, Condel_prediction, Condel_score.

=cut

package ExtraCols;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);
use Bio::EnsEMBL::Variation::Utils::VEP qw(@OUTPUT_COLS);

sub new {
    my $class = shift;
    my $self = $class->SUPER::new(@_);
    
    return $self;
}

sub version {
    return '2.6';
}

sub feature_types {
    return ['Transcript', 'RegulatoryFeature', 'MotifFeature', 'Intergenic', 'Gene', 'Exon'];
}

sub variant_feature_types {
    return ['VariationFeature', 'StructuralVariationFeature'];
}

sub get_header_info {
    my @new_output_cols = qw(
        CANONICAL
        CCDS
        HGNC
        ENSP
        HGVSc
        HGVSp
        EXON
        INTRON
        DOMAINS
        SIFT
        PolyPhen
        Condel
        MOTIF_NAME
        MOTIF_POS
        HIGH_INF_POS
        MOTIF_SCORE_CHANGE
        SIFT_prediction
        SIFT_score
        PolyPhen_prediction
        PolyPhen_score
        Condel_prediction
        Condel_score
        Conservation
    );
    @OUTPUT_COLS = (@OUTPUT_COLS, @new_output_cols);

    return {
        SIFT_prediction => "SIFT prediction (Ensembl)",
        SIFT_score => "SIFT score (Ensembl)",
        PolyPhen_prediction => "PolyPhen prediction (Ensembl)",
        PolyPhen_score => "PolyPhen score (Ensembl)",
        Condel_prediction => "Condel SIFT/PolyPhen consensus prediction (Ensembl)",
        Condel_score => "Condel SIFT/PolyPhen consensus score (Ensembl)",
    };
}

sub run {
    my ($self, $bvfoa, $line_hash) = @_;
    my $bvf = $bvfoa->base_variation_feature;
   foreach my $key (keys %{$line_hash->{Extra} || {}}){
        my $value = $line_hash->{Extra}->{$key};
        if($key eq 'SIFT' && $value =~ /([^(]+)\(([0-9.-]+)/){
            $line_hash->{SIFT_prediction} = $1;
            $line_hash->{SIFT_score} = $2;
        }
        elsif($key eq 'PolyPhen' && $value =~ /([^(]+)\(([0-9.-]+)/){
            $line_hash->{PolyPhen_prediction} = $1;
            $line_hash->{PolyPhen_score} = $2;
        }
        elsif($key eq 'Condel' && $value =~ /([^(]+)\(([0-9.-]+)/){
            $line_hash->{Condel_prediction} = $1;
            $line_hash->{Condel_score} = $2;
        }
        $line_hash->{$key} = $value;
    }
    return {};
}


1;

