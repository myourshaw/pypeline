=head1 LICENSE

 Copyright (c) 2011-2014 Michael Yourshaw.  All rights reserved.                                                                      

=head1 CONTACT                                                                                                       

 Michael Yourshaw <myourshaw@ucla.edu>
    
=cut

=head1 NAME

 EVS

=head1 SYNOPSIS

 mv EVS.pm ~/.vep/Plugins
 perl variant_effect_predictor.pl -i variations.vcf --plugin vw[,host,port,user,password,mysql,database]  --plugin EVS,NHLBI-EVS

=head1 DESCRIPTION

 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that
 adds the following new columns, derived from the NHLBI Exome Sequencing Project (ESP) Exome Variant Server:
  EA_AC_allele1,
  EA_AC_allele2,
  AA_AC_allele1,
  AA_AC_allele2,
  TAC_allele1,
  TAC_allele2,
  EA_GTC_this,
  AA_GTC_this,
  GTC_this,
  total_EA_GTC_this,
  total_AA_GTC_this,
  total_GTC_this,
  EA_allele1_freq,
  EA_allele2_freq,
  AA_allele1_freq,
  AA_allele2_freq,
  EA_AA_allele1_freq,
  EA_AA_allele2_freq,
  EA_GT_freq,
  AA_GT_freq,
  EA_AA_GT_freq
 
 Requires that the VAX.pm module be in the Plugins directory
 
 References:
   (1) Exome Variant Server, NHLBI GO Exome Sequencing Project (ESP), Seattle, WA (URL: http://evs.gs.washington.edu/EVS/)

=head1 PARAMETERS

    nhlbi_version (default: NHLBI-EVS)

=cut

package EVS;
use strict;
use warnings;

#require Exporter;
#our @ISA = qw(Exporter);
#our @EXPORT_OK = qw();

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);
use Bio::EnsEMBL::Variation::Utils::VEP qw(@OUTPUT_COLS);

use vw;
use VAX qw(get_bvfoa_info);

sub new {
    my $class = shift;
    my $self = $class->SUPER::new(@_);
    ($self->{nhlbi_version}) = @{$self->{params}} ? @{$self->{params}} : 'NHLBI-EVS';
    return $self;
}

sub version {
    return '75';
}

sub feature_types {
    return ['Transcript', 'RegulatoryFeature', 'MotifFeature', 'Intergenic', 'Gene', 'Exon'];
}

sub variant_feature_types {
    return ['VariationFeature', 'StructuralVariationFeature'];
}

sub get_header_info {
    my $self = shift;
    my $nhlbi_version = $self->{nhlbi_version};
    
    my %columns = (
        oALT_NHLBI => "Original NHLBI-EVS ALT field, possibly including multiple alleles ($nhlbi_version)",
        EA_AC_allele1 => "European American Allele Count for allele2 ($nhlbi_version)",
        EA_AC_allele2 => "European American Allele Count for allele2 ($nhlbi_version)",
        AA_AC_allele1 => "African American Allele Count for allele1 ($nhlbi_version)",
        AA_AC_allele2 => "African American Allele Count for allele2 ($nhlbi_version)",
        TAC_allele1 => "Total Allele Count of European American and African American for allele1 ($nhlbi_version)",
        TAC_allele2 => "Total Allele Count of European American and African American for allele2 ($nhlbi_version)",
        #MAF_EA => "EA minor Allele Frequency in percent  ",
        #MAF_AA => "AA minor Allele Frequency in percent ($nhlbi_version)",
        #MAF_All => "All minor Allele Frequency in percent ($nhlbi_version)",
        #EA_GTC_this => "European American Sample Count for this genotype ($nhlbi_version)",
        #AA_GTC_this => "African American Samp Count for this genotype ($nhlbi_version)",
        #GTC_this => "Total Sample Count of European American and African American for this genotype ($nhlbi_version)",
        #AA => "chimpAllele ($nhlbi_version)",
        AA_AC => "African American Allele Count in the order of AltAlleles,RefAllele. For INDELs, A1, A2, or An refers to the N-th alternate allele while R refers to the reference allele. ($nhlbi_version)",
        #AA_AGE => "Esitmated Variant Age in kilo years for the African American Population ($nhlbi_version)",
        AA_GTC => "African American Genotype Counts in the order of listed GTS ($nhlbi_version)",
        #CA => "clinicalAssociation ($nhlbi_version)",
        #CDS_SIZES => "Coding DNA Sizes ($nhlbi_version)",
        #CG => "consScoreGERP ($nhlbi_version)",
        #CP => "scorePhastCons ($nhlbi_version)",
        #DBSNP => "dbSNP version which established the rs_id ($nhlbi_version)",
        #DP => "Average Sample Read Depth ($nhlbi_version)",
        EA_AC => "European American Allele Count in the order of AltAlleles,RefAllele. For INDELs, A1, A2, or An refers to the N-th alternate allele while R refers to the reference allele. ($nhlbi_version)",
        #EA_AGE => "Esitmated Variant Age in kilo years for the European American Population ($nhlbi_version)",
        EA_GTC => "European American Genotype Counts in the order of listed GTS ($nhlbi_version)",
        #EXOME_CHIP => "Whether a SNP is on the Illumina HumanExome Chip ($nhlbi_version)",
        #FG => "functionGVS ($nhlbi_version)",
        #GL => "geneList ($nhlbi_version)",
        #GS => "granthamScore ($nhlbi_version)",
        GTC => "Total Genotype Counts in the order of listed GTS ($nhlbi_version)",
        GTS => "Observed Genotypes. For INDELs, A1, A2, or An refers to the N-th alternate allele while R refers to the reference allele. ($nhlbi_version)",
        #GWAS_PUBMED => "PubMed records for GWAS hits ($nhlbi_version)",
        #HGVS_CDNA_VAR => "HGVS Coding DNA Variant ($nhlbi_version)",
        #HGVS_PROTEIN_VAR => "HGVS Protein Variant ($nhlbi_version)",
        #MAF => "Minor Allele Frequency in percent in the order of EA,AA,All ($nhlbi_version)",
        #PH => "polyPhen2 result including prediction class and score ($nhlbi_version)",
        TAC => "Total Allele Count in the order of AltAlleles,RefAllele For INDELs, A1, A2, or An refers to the N-th alternate allele while R refers to the reference allele. ($nhlbi_version)",
        EA_GTC_total => "Total Sample Count of European American for all genotypes ($nhlbi_version)",
        AA_GTC_total => "Total Sample Count of African American for all genotypes ($nhlbi_version)",
        GTC_total => "Total Sample Count of European American and African American for all genotypes ($nhlbi_version)",
        EA_allele1_freq => "European American Allele Frequency for allele1 ($nhlbi_version)",
        EA_allele2_freq => "European American Allele Frequency for for allele2 ($nhlbi_version)",
        AA_allele1_freq => "African American Allele Frequency for for allele1 ($nhlbi_version)",
        AA_allele2_freq => "African American Allele Frequency for for allele2 ($nhlbi_version)",
        EA_AA_allele1_freq => "European American and African American Allele Frequency for for allele1 ($nhlbi_version)",
        EA_AA_allele2_freq => "European American and African American Allele Frequency for for allele2 ($nhlbi_version)",
        EA_GT_freq => "European American Genotype Frequency for this genotype ($nhlbi_version)",
        AA_GT_freq => "African American Genotype Frequency for this genotype ($nhlbi_version)",
        EA_AA_GT_freq => "European American and African American Genotype Frequency for this genotype ($nhlbi_version)",
    );

    my @output_cols = sort keys %columns;
    $self->{output_cols} = \@output_cols;
    @OUTPUT_COLS = (@OUTPUT_COLS, @output_cols);

    return \%columns;

}

sub run {
    my ($self, $bvfoa, $line_hash) = @_;
    my $bvf = $bvfoa->base_variation_feature;    
    my %bvfoa_info = %{get_bvfoa_info(@_)};
    my $input_line = $bvfoa_info{_line};

    my @output_cols = @{$self->{output_cols}};

    my ($chrom,$pos,$id,$ref,$alt) = split(/\t/,$input_line);

    if (defined $chrom && defined $pos && defined $ref && defined $alt){
        my $query = "CALL $vw::vw_database.get_nhlbi_evs('$chrom',$pos,'$ref','$alt')";
        my $qh = $vw::vw_conn->prepare($query);
        $qh->execute() or die "Unable to execute $query: $DBI::errstr\n";
        my $hash_ref = $qh->fetchrow_hashref();
        foreach my $c(@output_cols){
            $line_hash->{$c} = defined($hash_ref->{$c}) ? $hash_ref->{$c} : '';
            
        }
    }    
    return {};
}


1;
