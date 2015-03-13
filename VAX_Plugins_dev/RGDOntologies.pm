=head1 LICENSE

 Copyright (c) 2011-2014 Michael Yourshaw.  All rights reserved.                                                                      

=head1 CONTACT                                                                                                       

 Michael Yourshaw <myourshaw@ucla.edu>
    
=cut

=head1 NAME

 dbSNP

=head1 SYNOPSIS

 mv RGDOntologies.pm ~/.vep/Plugins
 perl variant_effect_predictor.pl -i variations.vcf --plugin vw[,host,port,user,password,mysql,database]  --plugin RGDOntologies,rgd_version

=head1 DESCRIPTION

 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that
 adds the following new columns, derived from ontologies gosted by the RAD Genome Database (RGD):
  biological_process,
  cellular_component,
  ChEBI_ontology,
  human_phenotype,
  molecular_function,
  mammalian_phenotype,
  neuro_behavioral_ontology,
  pathway_ontology,
  RGD_disease_ontology,
 
 Requires that the VAX.pm module be in the Plugins directory

=head1 PARAMETERS

    rgd_version (default: RGD)

=cut

package RGDOntologies;
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
    ($self->{rgd_version}) = @{$self->{params}} ? @{$self->{params}} : 'RGD';
    return $self;
}

sub version {
    return '75';
}

sub feature_types {
    return ['Transcript', 'Gene', 'Exon'];
}

sub variant_feature_types {
    return ['VariationFeature', 'StructuralVariationFeature'];
}

sub get_header_info {
    my $self = shift;
    my $rgd_version = $self->{rgd_version};
    my %columns = (
        biological_process => "Biological Process ontology. Gene Ontology ($rgd_version)",
        cellular_component => "Cellular Component ontology. Gene Ontology ($rgd_version)",
        ChEBI_ontology => "Chemical Entities of Biological Interest (ChEBI) ontology. Degtyarenko et al. (2008) ChEBI: a database and ontology for chemical entities of biological interest. Nucleic Acids Res. 36, D344–D350 ($rgd_version)",
        human_phenotype => "Human Phenotype Ontology (HPO). http://www.human-phenotype-ontology.org/ ($rgd_version)",
        molecular_function => "Molecular Function ontology. Gene Ontology ($rgd_version)",
        mammalian_phenotype => "Mammalian Phenotype ontology. Source:MGI ($rgd_version)",
        neuro_behavioral_ontology => "Neuro Behavior Ontology (NBO). http://bioportal.bioontology.org/ontologies/1621 ($rgd_version)",
        pathway_ontology => "Pathway Ontology (PW), currently being developed at the Rat Genome Database. Pathway Ontology (PW), is currently being developed at the Rat Genome Database ($rgd_version)",
        RGD_disease_ontology => "RGD's Disease Ontology (RDO) is the MEDIC vocabulary developed and maintained by the Comparative Toxicogenomics Database (CTD). http://ctdbase.org/reports/CTD_diseases.obo.gz ($rgd_version)",
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

    if (defined $bvfoa_info{hgnc}){
        my $query = "CALL $vw::vw_database.hgnc2rgd_ontologies('$bvfoa_info{hgnc}')";
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
