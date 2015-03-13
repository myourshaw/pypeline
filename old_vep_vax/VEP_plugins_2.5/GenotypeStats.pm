=head1 LICENSE

 Copyright (c) 2011-2012 Michael Yourshaw.  All rights reserved.                                                                      

=head1 CONTACT                                                                                                       

 Michael Yourshaw <myourshaw@ucla.edu>
    
=cut

=head1 NAME

 GenotypeStats

=head1 SYNOPSIS

 mv GenotypeStats.pm ~/.vep/Plugins
 perl variant_effect_predictor.pl -i variations.vcf --plugin GenotypeStats

=head1 DESCRIPTION

 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that
 adds the following new columns:
 GT_ALLELE_COUNT,GT_REF_ALLELE_COUNT,GT_ALT_ALLELE_COUNT,GT_REF_ALLELE_FREQUENCY,GT_ALT_ALLELE_FREQUENCY,GT_SAMPLE_COUNT,GT_HET_SAMPLE_COUNT,GT_HOM_REF_SAMPLE_COUNT,GT_HOM_ALT_SAMPLE_COUNT.
 
=cut

package GenotypeStats;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);
use Bio::EnsEMBL::Variation::Utils::VEP qw(@OUTPUT_COLS);

use VAX qw(get_bvfoa_info);

sub new {
    my $class = shift;
    my $self = $class->SUPER::new(@_);
    
    return $self;
}

sub version {
    return '2.5';
}

sub feature_types {
    return ['Transcript', 'RegulatoryFeature', 'MotifFeature', 'Intergenic', 'Gene', 'Exon'];
}

sub variant_feature_types {
    return ['VariationFeature', 'StructuralVariationFeature'];
}

sub get_header_info {
    my @new_output_cols = qw(
        GT_ALLELE_COUNT
        GT_REF_ALLELE_COUNT
        GT_ALT_ALLELE_COUNT
        GT_REF_ALLELE_FREQUENCY
        GT_ALT_ALLELE_FREQUENCY
        GT_SAMPLE_COUNT
        GT_HET_SAMPLE_COUNT
        GT_HOM_REF_SAMPLE_COUNT
        GT_HOM_ALT_SAMPLE_COUNT
    );
    @OUTPUT_COLS = (@OUTPUT_COLS, @new_output_cols);

    return {
        GT_ALLELE_COUNT => "number of alleles in called genotypes (VCF file)",
        GT_REF_ALLELE_COUNT => "number of reference alleles in called genotypes (VCF file)",
        GT_ALT_ALLELE_COUNT => "number of alternate alleles in called genotypes (VCF file)",
        GT_REF_ALLELE_FREQUENCY => "reference allele frequency in called genotypes (VCF file)",
        GT_ALT_ALLELE_FREQUENCY => "alternate allele frequency in called genotypes (VCF file)",
        GT_SAMPLE_COUNT => "number of samples  in called genotypes (VCF file)",
        GT_HET_SAMPLE_COUNT => "number of heterozygous samples  in called genotypes (VCF file)",
        GT_HOM_REF_SAMPLE_COUNT => "number of homozygous reference samples  in called genotypes  (VCF file)",
        GT_HOM_ALT_SAMPLE_COUNT => "number of homozygous alternate samples  in called genotypes (VCF file)",
    };
}

sub run {
    my ($self, $bvfoa, $line_hash) = @_;
    my $bvf = $bvfoa->base_variation_feature;    
    my %bvfoa_info = %{get_bvfoa_info(@_)};

    my $input_line = $bvfoa_info{_line};

    my @vcf_data = split("\t", $input_line);
    
    my ($CHROM,$POS,$ID,$REF,$ALT,$QUAL,$FILTER,$INFO,$FORMAT) = 0..8;
    
    if($#vcf_data>8){

        my @gts = @vcf_data[9..$#vcf_data];
        my @alt_list = split(/,/,$vcf_data[$ALT]);
        my @these_alleles = ($vcf_data[$REF], @alt_list);

        my $sample_count = 0;
        my $allele_count = 0;
        my $ref_allele_count = 0;
        my $alt_allele_count = 0;
        my $het_sample_count = 0;
        my $hom_ref_sample_count = 0;
        my $hom_alt_sample_count = 0;
        my $ref_allele_frequency = 0.0;
        my $alt_allele_frequency = 0.0;

        foreach my $this_gt(@gts){
            my @split_gt = split(/:/,$this_gt);
            my $gt = $split_gt[0];
            my ($allele1x, $phase, $allele2x) = ('.','/','.');
            if($gt =~ /([0-9.]+)([\/\|])([0-9.]+)/){
                ($allele1x, $phase, $allele2x) = ($1,$2,$3);
            }
            if($allele1x ne '.' && $allele2x ne '.'){
                my ($allele1, $allele2) = ($these_alleles[$allele1x], $these_alleles[$allele2x]);
                $sample_count+=1;
                $allele_count+=2;
                if ($allele1 eq $vcf_data[$REF] and $allele2 eq $vcf_data[$REF]){
                    $hom_ref_sample_count+=1;
                    $ref_allele_count+=2;
                }
                elsif ($allele1 ne $vcf_data[$REF] and $allele2 ne $vcf_data[$REF]){
                    $hom_alt_sample_count+=1;
                    $alt_allele_count+=2;
                }
                elsif (($allele1 ne $vcf_data[$REF] and $allele2 eq $vcf_data[$REF]) or ($allele2 ne $vcf_data[$REF] and $allele1 eq $vcf_data[$REF])){
                    $het_sample_count+=1;
                    $alt_allele_count+=1;
                    $ref_allele_count+=1;
                }
            }
        }
        if($allele_count != 0){
            $ref_allele_frequency = $ref_allele_count/$allele_count;
            $alt_allele_frequency = $alt_allele_count/$allele_count;
        }
        $line_hash->{GT_ALLELE_COUNT} = $allele_count;
        $line_hash->{GT_REF_ALLELE_COUNT} = $ref_allele_count;
        $line_hash->{GT_ALT_ALLELE_COUNT} = $alt_allele_count;
        $line_hash->{GT_REF_ALLELE_FREQUENCY} = $ref_allele_frequency;
        $line_hash->{GT_ALT_ALLELE_FREQUENCY} = $alt_allele_frequency;
        $line_hash->{GT_SAMPLE_COUNT} = $sample_count;
        $line_hash->{GT_HET_SAMPLE_COUNT} = $het_sample_count;
        $line_hash->{GT_HOM_REF_SAMPLE_COUNT} = $hom_ref_sample_count;
        $line_hash->{GT_HOM_ALT_SAMPLE_COUNT} = $hom_alt_sample_count;
    }

    return {};
}


1;

