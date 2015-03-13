#!/usr/bin/perl -w
use strict;
use warnings;
use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'cortex.local',
    -user => 'ensembl',
    pass => 'ensembl',
    port => 3306,
    
);

# Fetch a variation object
my $var_adaptor = $registry->get_adaptor('human', 'variation', 'variation');
#my $var = $var_adaptor->fetch_by_name('rs1421085');
my $gene = 'PCSK1';
my @variants = ('rs7482144','rs11692396 ','rs4687319','rs7194883', 'rs121434357', 'rs61752123');

print "phenotype_description\tvariation_names\tsource_name\tp_value\tassociated_variant_risk_allele\n";

# Fetch all the variation annotations associated with the variation
my $va_adaptor = $registry->get_adaptor('homo_sapiens', 'variation', 'variationannotation');
foreach my $v(@variants){
    my $var = $var_adaptor->fetch_by_name($v);
    #foreach my $va (@{$va_adaptor->fetch_all_by_associated_gene($gene)}) {
    foreach my $va (@{$va_adaptor->fetch_all_by_Variation($var)}) {
        print $va->phenotype_description . "\t" . $va->variation_names . "\t" . $va->source_name;

        if (defined($va->p_value)){
            print "\t" . $va->p_value;
        }
        
        if (defined($va->associated_variant_risk_allele)){
            my @risk_allele_array = split(/\-/, $va->associated_variant_risk_allele);
            if (scalar(@risk_allele_array)>1){
                print "\t".$risk_allele_array[1];
            }
        }
        
        print ".\n";
    }
}

#Variation rs266695 is associated with the phenotype 'Early onset extreme obesity' in the source Open Access GWAS Database
#Use of uninitialized value $risk_allele in print at /home/myourshaw/lab/pypeline/vax23/variant_annotation_test.pl line 34.
# with a p-value of 0.0268.
#Variation rs113318 is associated with the phenotype 'Type II Diabetes Mellitus-Fasting insulin' in the source Open Access GWAS Database with a p-value of 0.00025.
#Variation rs1995186 is associated with the phenotype 'Type II Diabetes Mellitus-Fasting insulin' in the source Open Access GWAS Database with a p-value of 0.00054.
#Variation CM971141 is associated with the phenotype 'Annotated by HGMD' in the source HGMD-PUBLIC.
#Variation CM074404 is associated with the phenotype 'Annotated by HGMD' in the source HGMD-PUBLIC.
#Variation CM034110 is associated with the phenotype 'Annotated by HGMD' in the source HGMD-PUBLIC.
#Variation CD034167 is associated with the phenotype 'Annotated by HGMD' in the source HGMD-PUBLIC.
#Variation COSM82279 is associated with the phenotype 'COSMIC:tumour_site:ovary' in the source COSMIC.
#Variation COSM78277 is associated with the phenotype 'COSMIC:tumour_site:ovary' in the source COSMIC.
#Variation rs6235 is associated with the phenotype 'Proinsulin levels' in the source NHGRI_GWAS_catalog with a p-value of 1e-26.
#The risk allele is G.
#Variation rs6234 is associated with the phenotype 'BODY MASS INDEX QUANTITATIVE TRAIT LOCUS 12' in the source OMIM.
#Variation rs6235 is associated with the phenotype 'BODY MASS INDEX QUANTITATIVE TRAIT LOCUS 12' in the source OMIM.
#Variation rs6232 is associated with the phenotype 'OBESITY, SUSCEPTIBILITY TO' in the source OMIM.
#The risk allele is .
#Variation rs6234 is associated with the phenotype 'PROPROTEIN CONVERTASE, SUBTILISIN/KEXIN-TYPE, 1' in the source OMIM.
#Variation rs6235 is associated with the phenotype 'PROPROTEIN CONVERTASE, SUBTILISIN/KEXIN-TYPE, 1' in the source OMIM.
