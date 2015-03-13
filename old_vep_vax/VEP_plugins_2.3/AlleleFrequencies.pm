=head1 LICENSE

 Copyright (c) 2011-2012 Michael Yourshaw.  All rights reserved.                                                                      

=head1 CONTACT                                                                                                       

 Michael Yourshaw <myourshaw@ucla.edu>
    
=cut

=head1 NAME

 AlleleFrequencies

=head1 SYNOPSIS

 mv AlleleFrequencies.pm ~/.vep/Plugins
 perl variant_effect_predictor.pl -i variations.vcf --plugin AlleleFrequencies,cortex.local,3306,vw,vw,mysql,vw

=head1 DESCRIPTION

 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that
 adds the following new columns:
 CLN,G5,G5A,SCS,in_dbsnp135,in_1kg,in_nhlbi,in_niehs,allele_count,ref_allele_count,alt_allele_count,ref_allele_frequency,alt_allele_frequency,sample_count,het_sample_count,hom_ref_sample_count,hom_alt_sample_count.
 
=head1 PARAMETERS

=cut

package AlleleFrequencies;

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
    my @new_output_cols = qw(
        CLN
        G5
        G5A
        SCS
        in_dbsnp135
        in_1kg
        in_nhlbi
        in_niehs
        allele_count
        ref_allele_count
        alt_allele_count
        ref_allele_frequency
        alt_allele_frequency
        sample_count
        het_sample_count
        hom_ref_sample_count
        hom_alt_sample_count
    );
    @OUTPUT_COLS = (@OUTPUT_COLS, @new_output_cols);

    return {
        CLN => "Variant is Clinical(LSDB,OMIM,TPA,Diagnostic) (dbSNP135)",
        G5 => ">5% minor allele frequency in 1+ populations (dbSNP135)",
        G5A => ">5% minor allele frequency in each and all populations (dbSNP135)",
        SCS => "Variant Clinical Significance, 0 - unknown, 1 - untested, 2 - non-pathogenic, 3 - probable-non-pathogenic, 4 - probable-pathogenic, 5 - pathogenic, 6 - drug-response, 7 - histocompatibility, 255 - other (dbSNP135)",
        in_dbsnp135 => "variant is in dbSNP 135 (dbSNP135)",
        in_1kg => "variant is in 1000 Genomes Phase 1 Integrated Variant Call Set updated 2011-12-09 (1000Genomes)",
        in_nhlbi => "variant is in Exome Variant Server,  NHLBI Exome Sequencing Project (ESP), Seattle, WA (URL: http://evs.gs.washington.edu/EVS/) ESP5400 December 10, 2011 (NHLBI)",
        in_niehs => "variant is in NIEHS Environmental Genome Project, Seattle, WA (URL: http://evs.gs.washington.edu/niehsExome/) version 0.0.6 September 30, 2011 (NIEHS)",
        allele_count => "number of alleles in called genotypes (dbSNP135,1000Genomes,NHLBI,NIEHS)",
        ref_allele_count => "number of reference alleles in called genotypes (dbSNP135,1000Genomes,NHLBI,NIEHS)",
        alt_allele_count => "number of alternate alleles in called genotypes (dbSNP135,1000Genomes,NHLBI,NIEHS)",
        ref_allele_frequency => "reference allele frequency (dbSNP135,1000Genomes,NHLBI,NIEHS)",
        alt_allele_frequency => "alternate allele frequency (dbSNP135,1000Genomes,NHLBI,NIEHS)",
        sample_count => "number of samples in combined datasets (dbSNP135,1000Genomes,NHLBI,NIEHS)",
        het_sample_count => "number of heterozygous samples in combined datasets (dbSNP135,1000Genomes,NHLBI,NIEHS)",
        hom_ref_sample_count => "number of homozygous reference samples in combined datasets  (dbSNP135,1000Genomes,NHLBI,NIEHS)",
        hom_alt_sample_count => "number of homozygous alternate samples in combined datasets (dbSNP135,1000Genomes,NHLBI,NIEHS)",
    };
}

sub run {
    my ($self, $tva, $line_hash) = @_;
    
    my %tva_info = %{get_tva_info(@_)};
    my $input_line = $tva_info{_line};

    my ($chrom,$pos,$id,$ref,$alt) = split(/\t/,$input_line);

    if (defined $chrom && defined $pos && defined $ref && defined $alt){
        my $query = "CALL $vw::vw_database.coord2af('$chrom',$pos,'$ref','$alt')";
        my $qh = $vw::vw_conn->prepare($query);
        $qh->execute() or die "Unable to execute $query: $DBI::errstr\n";
        while (my @row = $qh->fetchrow_array()){
            $line_hash->{CLN} = defined($row[0]) ? $row[0] : '';
            $line_hash->{G5} = defined($row[1]) ? $row[1] : '';
            $line_hash->{G5A} = defined($row[2]) ? $row[2] : '';
            $line_hash->{SCS} = defined($row[3]) ? $row[3] : '';
            $line_hash->{in_dbsnp135} = defined($row[4]) ? $row[4] : '';
            $line_hash->{in_1kg} = defined($row[5]) ? $row[5] : '';
            $line_hash->{in_nhlbi} = defined($row[6]) ? $row[6] : '';
            $line_hash->{in_niehs} = defined($row[7]) ? $row[7] : '';
            $line_hash->{allele_count} = defined($row[8]) ? $row[8] : '';
            $line_hash->{ref_allele_count} = defined($row[9]) ? $row[9] : '';
            $line_hash->{alt_allele_count} = defined($row[10]) ? $row[10] : '';
            $line_hash->{ref_allele_frequency} = defined($row[11]) ? $row[11] : '';
            $line_hash->{alt_allele_frequency} = defined($row[12]) ? $row[12] : '';
            $line_hash->{sample_count} = defined($row[13]) ? $row[13] : '';
            $line_hash->{het_sample_count} = defined($row[14]) ? $row[14] : '';
            $line_hash->{hom_ref_sample_count} = defined($row[15]) ? $row[15] : '';
            $line_hash->{hom_alt_sample_count} = defined($row[16]) ? $row[16] : '';
            last;
        }
    }
    return {};
}


1;

