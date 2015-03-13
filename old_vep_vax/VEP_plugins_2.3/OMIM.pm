=head1 LICENSE

 Copyright (c) 2011-2012 Michael Yourshaw.  All rights reserved.                                                                      

=head1 CONTACT                                                                                                       

 Michael Yourshaw <myourshaw@ucla.edu>
    
=cut

=head1 NAME

 OMIM

=head1 SYNOPSIS

 mv OMIM.pm ~/.vep/Plugins
 perl variant_effect_predictor.pl -i variations.vcf --plugin OMIM

=head1 DESCRIPTION

 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that
 adds the following new column:
 OMIM_Disorder.
 All fields are uri escaped for ';='.
 
 Requires vw plugin
 
=cut

package OMIM;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);
use Bio::EnsEMBL::Variation::Utils::VEP qw(@OUTPUT_COLS);

use vw;
use VAX qw(get_tva_info get_unique);

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
        OMIM_Disorder
    );
    @OUTPUT_COLS = (@OUTPUT_COLS, @new_output_cols);

    return {
        OMIM_Disorder => "OMIM disorder(s) for gene (| separated) (OMIM)",
    };
}

sub run {
    my ($self, $tva, $line_hash) = @_;

    my %tva_info = %{get_tva_info(@_)};

    my @DISEASES_PHENOTYPES;
    
    if (defined $tva_info{gene}){
        my @data;
        
        my $query = "CALL $vw::vw_database.ensg2omim('$tva_info{ensg}')";
        my $qh = $vw::vw_conn->prepare($query);
        $qh->execute() or die "Unable to execute $query: $DBI::errstr\n";
        while (my @row = $qh->fetchrow_array()){
            push @data, $row[0];
            
            push(@DISEASES_PHENOTYPES, split(/\|/,$row[0]));
        }
        $line_hash->{OMIM_Disorder} = @data ? join('|', @data) : '';
    }

    my @diseases_phenotypes = defined($line_hash->{DISEASES_PHENOTYPES}) ? split(/\|/,$line_hash->{DISEASES_PHENOTYPES}) : ();
    my @union_diseases_phenotypes = (@diseases_phenotypes, @DISEASES_PHENOTYPES);
    my @unique_union_diseases_phenotypes = @{get_unique(\@union_diseases_phenotypes)};
    if(@unique_union_diseases_phenotypes){
        @unique_union_diseases_phenotypes = sort(@unique_union_diseases_phenotypes);
        my $diseases_phenotypes_str = join('|',@unique_union_diseases_phenotypes);
        $line_hash->{DISEASES_PHENOTYPES} = $diseases_phenotypes_str;
    }
    
    return {};
}


1;

