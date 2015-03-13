=head1 LICENSE

 Copyright (c) 2011-2014 Michael Yourshaw.  All rights reserved.                                                                      

=head1 CONTACT                                                                                                       

 Michael Yourshaw <myourshaw@ucla.edu>
    
=cut

=head1 NAME

 GenotypeFile

=head1 SYNOPSIS

 mv GenotypeFile.pm ~/.vep/Plugins
 perl variant_effect_predictor.pl -i variations.vcf --plugin GenotypeFile

=head1 DESCRIPTION

 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that
 creates an additional output file <output_file>.gt.txt with one row per sample genotype
 and the following new columns:
 SAMPLE, GT, ALLELE1, PHASE, ALLELE2, ZYGOSITY.
 ZYGOSITY = '' for ./.
            0 for homozygous REF
            1 for heterozygous REF/ALT
            2 for homozygous ALT
            >2 for either allele not REF or ALT such as het two different non-refs

=cut

package GenotypeFile;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);
use Bio::EnsEMBL::Variation::Utils::VEP qw(@OUTPUT_COLS);

use VAX qw(get_bvfoa_info);

sub new {
    my $class = shift;
    my $self = $class->SUPER::new(@_);
    
    my $input_file = $self->{config}->{input_file};
    open VCF, '<', $input_file or die "Can't open $input_file ($!)";
    while(<VCF>){
        if(/^\#CHROM\s+POS\s+ID\sREF\sALT/){
            chomp;
            my @vcf_cols = split("\t");
            $vcf_cols[0] = substr($vcf_cols[0],1);
            $self->{_vcf_cols} = \@vcf_cols;
            if ($#vcf_cols>8){
                my @gt_cols = @vcf_cols[9..$#vcf_cols];
                $self->{_vcf_gt_cols} = \@gt_cols;
            }
            last;
        }
    }
    close VCF;

    if($self->{config}->{output_file} =~ /stdout/i) {
        $self->{genotype_output_file} = '~/vax.gt.txt';
    }
    else {
        $self->{genotype_output_file} = $self->{config}->{output_file}.'.gt.txt';
    }
    
    open(GTOUTFH, ">", $self->{genotype_output_file}) or die("ERROR: Could not write to genotype output file ", $self->{genotype_output_file}, "\n");
    $self->{genotype_output_file_handle} = *GTOUTFH;
   
    $self->{genotype_header_written} = 0;
    
    return $self;
}

sub version {
    return '75';
}

sub feature_types {
    return ['Transcript', 'RegulatoryFeature', 'MotifFeature', 'Intergenic'];
}

sub get_header_info {
    return {};
}

sub run {
    my ($self, $bvfoa, $line_hash) = @_;
    my $bvf = $bvfoa->base_variation_feature;    
    my %line = %{$line_hash};
    $line{Extra} = join ';', map { $_.'='.$line{Extra}->{$_} } keys %{ $line{Extra} || {} };
    
    my %bvfoa_info = %{get_bvfoa_info(@_)};

    my $input_line = $bvfoa_info{_line};

    my @vcf_data = split("\t", $input_line);
    map {$line{$self->{_vcf_cols}[$_]} = $vcf_data[$_]} 0..$#{$self->{_vcf_cols}};
    
    if(defined($self->{_vcf_gt_cols})
       and defined($line{REF})
       and defined($line{ALT})
       ){
        my @output_columns = (@OUTPUT_COLS,'SAMPLE','GT','ALLELE1','PHASE','ALLELE2','ZYGOSITY');
        
        my $fh = $self->{genotype_output_file_handle};
        unless($self->{genotype_header_written}){
            if(substr($output_columns[0],0,1) ne '#'){
                $output_columns[0] = '#'.$output_columns[0];
            }
            print $fh join("\t", @output_columns)."\n" ;
            $self->{genotype_header_written} = 1;
        }

        if ($#vcf_data>8){
            my @samples = @{$self->{_vcf_cols}}[9..$#{$self->{_vcf_cols}}];
            my @gts = @vcf_data[9..$#vcf_data];
            my @alt_list = split(/,/,$line{ALT});
            my @these_alleles = ($line{REF}, @alt_list);
            for my $i(0..$#samples){
                $line{SAMPLE} = $samples[$i];
                my $this_gt = $gts[$i];
                my @split_gt = split(/:/,$this_gt);
                $line{GT} = $split_gt[0];
                my ($allele1x, $phase, $allele2x) = ('.','/','.');
                if($line{GT} =~ /([0-9.]+)([\/\|])([0-9.]+)/){
                    ($allele1x, $line{PHASE}, $allele2x) = ($1,$2,$3);
                }
                ($line{ALLELE1}, $line{ALLELE2}) = ($allele1x ne '.' ? $these_alleles[$allele1x] : '.', $allele2x ne '.' ? $these_alleles[$allele2x] : '.');
                $line{ZYGOSITY} = $line{ALLELE1} eq '.' || $line{ALLELE2} eq '.' ? ''
                    : $line{ALLELE1} eq $line{REF} && $line{ALLELE2} eq $line{REF} ? 0
                    : $line{ALLELE1} eq $line{ALT} && $line{ALLELE2} eq $line{ALT} ? 2
                    : ($line{ALLELE1} eq $line{REF} && $line{ALLELE2} eq $line{ALT}) || ($line{ALLELE1} eq $line{ALT} && $line{ALLELE2} eq $line{REF}) ? 1
                    : $allele1x + $allele2x;
                my $output = join("\t", map {defined($line{$_}) ? $line{$_} : ''} @output_columns);
                print $fh $output."\n";
            }
        }
    }

    return {};
}


1;

