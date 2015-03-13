=head1 LICENSE

 Copyright (c) 2011-2012 Michael Yourshaw.  All rights reserved.                                                                      

=head1 CONTACT                                                                                                       

 Michael Yourshaw <myourshaw@ucla.edu>
    
=cut

=head1 NAME

 VCFCols

=head1 SYNOPSIS

 mv VCFCols.pm ~/.vep/Plugins
 perl variant_effect_predictor.pl -i variations.vcf --plugin VCFCols

=head1 DESCRIPTION

 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that
 adds all input VCF columns as the first columns of the output:
 CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, [FORMAT, GT[, GT ...]].
 
=cut

package VCFCols;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);
use Bio::EnsEMBL::Variation::Utils::VEP qw(@OUTPUT_COLS);
use VAX qw(get_tva_info);

sub new {
    my $class = shift;
    my $self = $class->SUPER::new(@_);
    
    my $input_file = $self->{config}->{input_file};
    open VCF, '<', $input_file or die "Can't open $input_file ($!)";
    while(<VCF>){
        chomp;
        if(/^\#CHROM\s+POS\s+ID\sREF\sALT/){
            my @vcf_cols = split("\t");
            $vcf_cols[0] = substr($vcf_cols[0],1);
            $self->{_vcf_cols} = \@vcf_cols;
            @OUTPUT_COLS = (@vcf_cols, @OUTPUT_COLS);
            last;
        }
    }
    close VCF;

    return $self;
}

sub version {
    return '2.3';
}

sub feature_types {
    return ['Bio::EnsEMBL::Transcript', 'Bio::EnsEMBL::Funcgen::RegulatoryFeature', 'Bio::EnsEMBL::Funcgen::MotifFeature'];
}

sub get_header_info {
    my $self = shift;
    
    my %vcf_cols;
    if (defined($self->{_vcf_cols})){
        foreach my $vcf_col(@{$self->{_vcf_cols}}){
            if(uc($vcf_col) eq 'CHROM'){
                $vcf_cols{CHROM} = "Chromosome (VCF)";
            }
            elsif(uc($vcf_col) eq 'POS'){
                $vcf_cols{'POS'} = "Position (VCF)";
            }
            elsif(uc($vcf_col) eq 'ID'){
                $vcf_cols{'ID'} = "ID (VCF)";
            }
            elsif(uc($vcf_col) eq 'REF'){
                $vcf_cols{'REF'} = "Reference allele (VCF)";
            }
            elsif(uc($vcf_col) eq 'ALT'){
                $vcf_cols{'ALT'} = "Alternate allele (VCF)";
            }
            elsif(uc($vcf_col) eq 'QUAL'){
                $vcf_cols{'QUAL'} = "Quality score (VCF)";
            }
            elsif(uc($vcf_col) eq 'FILTER'){
                $vcf_cols{'FILTER'} = "Filter (VCF)";
            }
            elsif(uc($vcf_col) eq 'INFO'){
                $vcf_cols{'INFO'} = "Info (VCF)";
            }
            elsif(uc($vcf_col) eq 'FORMAT'){
                $vcf_cols{'FORMAT'} = "Format of genotype columns (VCF)";
            }
            else{
                $vcf_cols{$vcf_col} = "$vcf_col genotype (VCF)";
            }
        }
        
    }
    return \%vcf_cols;
}

sub run {
    my ($self, $tva, $line_hash) = @_;

    my %tva_info = %{get_tva_info(@_)};
    my $input_line = $tva_info{_line};

    my @vcf_data = split("\t", $input_line);
    map {$line_hash->{$self->{_vcf_cols}[$_]} = $vcf_data[$_]} 0..$#{$self->{_vcf_cols}};
    return {};
}


1;

