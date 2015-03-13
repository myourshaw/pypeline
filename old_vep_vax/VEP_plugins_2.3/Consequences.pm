=head1 LICENSE

 Copyright (c) 2011-2012 Michael Yourshaw.  All rights reserved.                                                                      

=head1 CONTACT                                                                                                       

 Michael Yourshaw <myourshaw@ucla.edu>
    
=cut

=head1 NAME

 Consequences

=head1 SYNOPSIS

 mv Consequences.pm ~/.vep/Plugins
 perl variant_effect_predictor.pl -i variations.vcf --plugin Consequences

=head1 DESCRIPTION

 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that
 adds the following new columns:
 Consequence_severest, Consequence_rank, Consequences_all. 

=cut

package Consequences;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);
use Bio::EnsEMBL::Variation::Utils::VEP qw(@OUTPUT_COLS);

use VAX qw(get_consequence_info);

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
    my @new_output_cols = qw(
        Consequence_severest
        Consequence_rank
        Consequences_all
    );
    @OUTPUT_COLS = (@OUTPUT_COLS, @new_output_cols);

    return {
        Consequence_severest => "Severest consequence to transcript of variant",
        Consequence_rank => "Rank of severest consequence (lower is more severe)",
        Consequences_all => "All consequences of variant to all transcripts",
    };
}

sub run {
    my ($self, $tva, $line_hash) = @_;
    
    my %consequences_info = %{get_consequence_info(@_)};

    $line_hash->{Consequence_severest} = $consequences_info{consequence_severest};
    $line_hash->{Consequence_rank} = $consequences_info{consequence_severest_rank};
    $line_hash->{Consequences_all} = join ",", @{$consequences_info{consequences_sorted}};
    return {};
}

1;
