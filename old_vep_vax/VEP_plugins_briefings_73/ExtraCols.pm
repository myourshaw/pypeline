=head1 LICENSE

 Copyright (c) 2011-2013 Michael Yourshaw.  All rights reserved.                                                                      

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
    return '73';
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
        SYMBOL
        ENSP
        HGVSc
        HGVSp
        EXON
        INTRON
        DOMAINS
        SIFT
        SIFT_prediction
        SIFT_score
        PolyPhen
        PolyPhen_prediction
        PolyPhen_score
        Condel
        Condel_prediction
        Condel_score
        CAROL
        CAROL_prediction
        CAROL_score
        FATHMM
        FATHMM_prediction
        FATHMM_score
        Conservation
        MOTIF_NAME
        MOTIF_POS
        HIGH_INF_POS
        MOTIF_SCORE_CHANGE
        GMAF
        AFR_MAF
        AMR_MAF
        AMR_MAF
        EUR_MAF
        AA_MAF
        EA_MAF
        CLIN_SIG
        BIOTYPE
        PUBMED
        ALLELE_NUM
        GO_VEP
        LRT_score
        LRT_pred
        MutationTaster_score
        MutationTaster_pred
        MutationAssessor_score
        MutationAssessor_pred
        FATHMM_pred
        GERP++_RS
        phyloP
        1000Gp1_AC
        1000Gp1_AF
        1000Gp1_AFR_AC
        1000Gp1_AFR_AF
        1000Gp1_EUR_AC
        1000Gp1_EUR_AF
        1000Gp1_AMR_AC
        1000Gp1_AMR_AF
        1000Gp1_ASN_AC
        1000Gp1_ASN_AF
        ESP6500_AA_AF
        ESP6500_EA_AF
    );
    @OUTPUT_COLS = (@OUTPUT_COLS, @new_output_cols);

    return {
        SIFT_prediction => "SIFT prediction (Ensembl)",
        SIFT_score => "SIFT score (Ensembl)",
        PolyPhen_prediction => "PolyPhen prediction (Ensembl)",
        PolyPhen_score => "PolyPhen score (Ensembl)",
        Condel_prediction => "Condel SIFT/PolyPhen consensus prediction (Ensembl)",
        Condel_score => "Condel SIFT/PolyPhen consensus score (Ensembl)",
        CAROL_prediction => "CAROL prediction (Ensembl)",
        CAROL_score => "CAROL score (Ensembl)",
        FATHMM_prediction => "FATHMM prediction (Ensembl)",
        FATHMM_score => "FATHMM score (Ensembl)",
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
        elsif($key eq 'CAROL' && $value =~ /([^(]+)\(([0-9.-]+)/){
            $line_hash->{CAROL_prediction} = $1;
            $line_hash->{CAROL_score} = $2;
        }
        elsif($key eq 'FATHMM' && $value =~ /([^(]+)\(([0-9.-]+)/){
            $line_hash->{FATHMM_prediction} = $1;
            $line_hash->{FATHMM_score} = $2;
        }
        if ($key eq 'GO') {
            $line_hash->{GO_VEP} = $value;
        }
        else{
            $line_hash->{$key} = $value;
        }
    }
    return {};
}


1;

