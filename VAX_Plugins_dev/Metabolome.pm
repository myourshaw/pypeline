=head1 LICENSE

 Copyright (c) 2013 Aliz R. Rao.  All rights reserved.                                                                      

=head1 CONTACT                                                                                                       

 Michael Yourshaw <myourshaw@ucla.edu>
    
=cut

=head1 NAME

 Metabolome

=head1 SYNOPSIS

 mv Metabolome.pm ~/.vep/Plugins
 perl variant_effect_predictor.pl -i variations.vcf --plugin vw[,host,port,user,password,mysql,database] --plugin Metabolome

=head1 DESCRIPTION

 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that
 adds the following new columns:
 metabolic_function, metabolites

 Requires that the vw plugin be in the Plugins directory and database installed with the VAX installer.
 
 Requires that the VAX.pm module be in the Plugins directory

 References:
    (1) Wishart DS, Tzur D, Knox C, et al., HMDB: the Human Metabolome Database. Nucleic Acids Res. 2007 Jan;35(Database issue):D521-6. 17202168 
    (2) Wishart DS, Knox C, Guo AC, et al., HMDB: a knowledgebase for the human metabolome. Nucleic Acids Res. 2009 37(Database issue):D603-610. 18953024 
    (3) Wishart DS, Jewison T, Guo AC, Wilson M, Knox C, et al., HMDB 3.0 — The Human Metabolome Database in 2013. Nucleic Acids Res. 2013. Jan 1;41(D1):D801-7. 23161693

=head1 PARAMETERS

=cut

package Metabolome;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);
use Bio::EnsEMBL::Variation::Utils::VEP qw(@OUTPUT_COLS);

use vw;
use VAX qw(get_bvfoa_info);

sub new {
    my $class = shift;
    my $self = $class->SUPER::new(@_);
    ($self->{hmdb_version}) = @{$self->{params}} ? @{$self->{params}} : '';

    return $self;
}

sub version {
    return '75';
}

sub feature_types {
    return ['Transcript', 'Gene', 'Exon',];
}

sub get_header_info {
    my $self = shift;
    my $hmdb_version = $self->{hmdb_version};
    
    my @output_cols = qw(
        metabolic_function
        metabolites
    );
    @OUTPUT_COLS = (@OUTPUT_COLS, @output_cols);
    
    return {
            metabolic_function => "Metabolic function of the protein (HMDB $hmdb_version)",
            metabolites => "Metabolites of the protein (HMDB $hmdb_version)"
        };
}

sub run {
    my ($self, $bvfoa, $line_hash) = @_;
    my $bvf = $bvfoa->base_variation_feature;    
    my %bvfoa_info = %{get_bvfoa_info(@_)};
    my $input_line = $bvfoa_info{_line};
    
    if (defined $bvfoa_info{hgnc}){
        my @data;

        my $query = "CALL $vw::vw_database.hgnc2metabolome('$bvfoa_info{hgnc}')";
        my $qh = $vw::vw_conn->prepare($query);
        $qh->execute() or die "Unable to execute $query: $DBI::errstr\n";
        my $hash_ref = $qh->fetchrow_hashref();
        $line_hash->{metabolic_function} = defined($hash_ref->{metabolic_function}) ? $hash_ref->{metabolic_function} : '';
        $line_hash->{metabolites} = defined($hash_ref->{metabolites}) ? $hash_ref->{metabolites} : '';
    }
    return {};
}


1;

