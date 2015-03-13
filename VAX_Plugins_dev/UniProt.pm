=head1 LICENSE

 Copyright (c) 2011-2014 Michael Yourshaw.  All rights reserved.                                                                      

=head1 CONTACT                                                                                                       

 Michael Yourshaw <myourshaw@ucla.edu>
    
=cut

=head1 NAME

 UniProt

=head1 SYNOPSIS

 mv UniProt.pm ~/.vep/Plugins
 perl variant_effect_predictor.pl -i variations.vcf --plugin vw[,host,port,user,password,mysql,database] --plugin UniProt[,version]

=head1 DESCRIPTION

 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that
 adds the following new columns:
 VARIANT, MUTAGEN, SITES, OTHER_OVERLAPPING_FEATURES, ALLERGEN, ALTERNATIVE_PRODUCTS,
 CATALYTIC_ACTIVITY, CAUTION, COFACTOR, DE, DEVELOPMENTAL_STAGE, DISEASE, DOMAIN,
 ENZYME_REGULATION, FUNCTION, GeneNames, INDUCTION, INTERACTION, KEGG,
 KEYWORDS, MIM_gene, MIM_phenotype, MISCELLANEOUS, PATHWAY, Pathway_Interaction, PE, POLYMORPHISM,
 PTM, Reactome, RecName, RefSeq_NM, RefSeq_NP, RNA_EDITING, SEQUENCE_CAUTION, SIMILARITY,
 SUBCELLULAR_LOCATION, SUBUNIT, TISSUE_SPECIFICITY, UCSC, WEB_RESOURCE.
 
 Requires that the vw plugin be in the Plugins directory and database installed with the VAX installer.
 
 Requires that the VAX.pm module be in the Plugins directory
 
 References:
    (1) Apweiler R, Martin MJ, O'Donovan C et al.
        Update on activities at the Universal Protein Resource (UniProt $uniprot_version) in 2013
        Nucleic Acids Research 2013;41:D43-D47

=cut

package UniProt;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);
use Bio::EnsEMBL::Variation::Utils::VEP qw(@OUTPUT_COLS);

use vw;
use VAX qw(get_bvfoa_info get_unique);

sub new {
    my $class = shift;
    my $self = $class->SUPER::new(@_);
    my $config = $self->{config};
    ($self->{uniprot_version}) = @{$self->{params}} ? @{$self->{params}} : '';
    
    return $self;
}

sub version {
    return '75';
}

sub feature_types {
    return ['Bio::EnsEMBL::Transcript'];
}

sub get_header_info {
    my $self = shift;
    my $uniprot_version = $self->{uniprot_version};

    my @output_cols = qw(
        VARIANT
        MUTAGEN
        SITES
        OTHER_OVERLAPPING_FEATURES    
        ALLERGEN
        ALTERNATIVE_PRODUCTS
        CATALYTIC_ACTIVITY
        CAUTION
        COFACTOR
        DE
        DEVELOPMENTAL_STAGE
        DISEASE
        DOMAIN
        ENZYME_REGULATION
        FUNCTION
        GeneNames
        INDUCTION
        INTERACTION
        KEGG
        KW
        MIM_gene
        MIM_phenotype
        MISCELLANEOUS
        PATHWAY
        Pathway_Interaction
        PE
        POLYMORPHISM
        PTM
        Reactome
        RecName
        RefSeq_NM
        RefSeq_NP
        RNA_EDITING
        SEQUENCE_CAUTION
        SIMILARITY
        SUBCELLULAR_LOCATION
        SUBUNIT
        TISSUE_SPECIFICITY
        UCSC
        WEB_RESOURCE
    );
    @OUTPUT_COLS = (@OUTPUT_COLS, @output_cols);

    return {
        VARIANT => "Description of a natural variant of the protein (UniProt $uniprot_version)",
        MUTAGEN => "Site which has been experimentally altered by mutagenesis. (UniProt $uniprot_version)",
        SITES => "ACT_SITE - Amino acid(s) involved in the activity of an enzyme; BINDING - Binding site for any chemical group (co-enzyme, prosthetic group, etc.); CA_BIND - Extent of a calcium-binding region; DISULFID - Disulfide bond; DNA_BIND - Extent of a DNA-binding region; METAL - Binding site for a metal ion; NP_BIND - Extent of a nucleotide phosphate-binding region; SITE - Any interesting single amino-acid site on the sequence, that is not defined by another feature key. It can also apply to an amino acid bond which is represented by the positions of the two flanking amino acids; ZN_FING - Extent of a zinc finger region (UniProt $uniprot_version)",
        OTHER_OVERLAPPING_FEATURES => "Other protein features not in VARIANT, MUTAGEN, SITES, which overlap with the position of the relevant amino acid. (UniProt $uniprot_version)",
        ACs => "ACcession number(s) associated with an entry (UniProt $uniprot_version)",
        ALLERGEN => "Information relevant to allergenic proteins (UniProt $uniprot_version)",
        ALTERNATIVE_PRODUCTS => "Description of the existence of related protein sequence(s) produced by alternative splicing of the same gene, alternative promoter usage, ribosomal frameshifting or by the use of alternative initiation codons",
        BIOPHYSICOCHEMICAL_PROPERTIES => "Description of the information relevant to biophysical and physicochemical data and information on pH dependence, temperature dependence, kinetic parameters, redox potentials, and maximal absorption (UniProt $uniprot_version)",
        BIOTECHNOLOGY => "Description of the use of a specific protein in a biotechnological process (UniProt $uniprot_version)",
        CATALYTIC_ACTIVITY => "Description of the reaction(s) catalyzed by an enzyme (UniProt $uniprot_version)",
        CAUTION => "Warning about possible errors and/or grounds for confusion (UniProt $uniprot_version)",
        COFACTOR => "Description of non-protein substance required by an enzyme to be active (UniProt $uniprot_version)",
        DE => "General descriptive information about the sequence stored. This information is generally sufficient to identify the protein precisely. (UniProt $uniprot_version)",
        DEVELOPMENTAL_STAGE => "Description of the developmentally-specific expression of mRNA or protein (UniProt $uniprot_version)",
        DISEASE => "Description of the disease(s) associated with a deficiency of a protein (UniProt $uniprot_version)",
        DISRUPTION_PHENOTYPE => "Description of the effects caused by the disruption of the gene coding for the protein. Note that we only describe effects caused the complete absence of a gene and thus a protein in vivo (null mutants caused by random or target deletions, insertions of a transposable element etc.) To avoid description of phenotypes due to partial or dominant negative mutants, missense mutations are not described in this comment, but in FT MUTAGEN instead. (UniProt $uniprot_version)",
        DOMAIN => "Description of the domain(s) present in a protein (UniProt $uniprot_version)",
        DRs => "Database cross-Reference pointers to information in external data resources (UniProt $uniprot_version)",
        #ENSG => "Database of automatically annotated sequences of large genomes gene identifier (Ensembl database) (UniProt $uniprot_version)",
        #ENSP => "Database of automatically annotated sequences of large genomes protein identifier (Ensembl database) (UniProt $uniprot_version)",
        #ENST => "Database of automatically annotated sequences of large genomes transcript identifier (Ensembl database) (UniProt $uniprot_version)",
        ENZYME_REGULATION => "Description of an enzyme regulatory mechanism (UniProt $uniprot_version)",
        FUNCTION => "General description of the function(s) of a protein (UniProt $uniprot_version)",
        #Gene => "The official gene name (UniProt $uniprot_version)",
        GeneNames => "(a.k.a gene symbols). The name(s) used to represent a gene (UniProt $uniprot_version)",
        #GO => "Gene Ontology (GO) database accession number/primary key (UniProt $uniprot_version)",
        #GO_term => "Gene Ontology (GO) database; this field is a 1-letter abbreviation for one of the 3 ontology aspects, separated from the GO term by a column. If the term is longer than 46 characters, the first 43 characters are indicated followed by 3 dots ('...'). The abbreviations for the 3 distinct aspects of the ontology are P (biological Process), F (molecular Function), and C (cellular Component) (UniProt $uniprot_version)",
        #HGNC => "Human gene nomenclature database (HGNC); the gene designation. If the gene designation is not available, a dash ('-') is used (UniProt $uniprot_version)",
        INDUCTION => "Description of the effects of environmental factors on the gene expression (UniProt $uniprot_version)",
        INTERACTION => "Interaction with other protein(s) (UniProt $uniprot_version)",
        KEGG => "Kyoto encyclopedia of genes and genomes database accession number/primary key (UniProt $uniprot_version)",
        KW => "List of controlled vocabulary which summarises the content of an entry (UniProt $uniprot_version)",
        MASS_SPECTROMETRY => "Reports the exact molecular weight of a protein or part of a protein as determined by mass spectrometric methods (UniProt $uniprot_version)",
        MIM_gene => "Mendelian Inheritance in Man Database (MIM) gene database accession number/primary key (UniProt $uniprot_version)",
        MIM_phenotype => "Mendelian Inheritance in Man Database (MIM) phenotype database accession number/primary key(UniProt $uniprot_version)",
        MISCELLANEOUS => "Any relevant information that doesn't fit in any other defined sections (UniProt $uniprot_version)",
        PATHWAY => "Description of associated metabolic pathways (UniProt $uniprot_version)",
        Pathway_Interaction => "NCI-Nature Pathway Interaction Database 'full pathway name' (UniProt $uniprot_version)",
        PE => "The evidence of the existence of a protein translated from this transcript (1-Evidence at protein level, 2-Evidence at transcript level, 3-Inferred from homology, 4-Predicted, 5-Uncertain) (UniProt $uniprot_version)",
        PHARMACEUTICAL => "Description of the use of a protein as a pharmaceutical drug (UniProt $uniprot_version)",
        POLYMORPHISM => "Description of polymorphism(s) (UniProt $uniprot_version)",
        PTM => "Description of post-translational modifications (UniProt $uniprot_version)",
        Reactome => "Curated resource of core pathways and reactions in human biology (Reactome) name of the pathway (UniProt $uniprot_version)",
        RecName => "The name recommended by the UniProt consortium (UniProt $uniprot_version)",
        RefSeq_NM => "NCBI reference sequences nucleotide sequence identifier (UniProt $uniprot_version)",
        RefSeq_NP => "NCBI reference sequences database accession number/primary key (UniProt $uniprot_version)",
        RNA_EDITING => "Description of amino acid change(s) due to RNA editing (UniProt $uniprot_version)",
        SEQUENCE_CAUTION => "Description of protein sequence reports that differ from the sequence that is shown in UniProtKB due to conflicts that are not described in FT CONFLICT lines, such as frameshifts, erroneous gene model predictions, etc. (UniProt $uniprot_version)",
        SIMILARITY => "Description of the sequence similaritie(s) with other proteins (UniProt $uniprot_version)",
        SQ => "Amino acid sequence (UniProt $uniprot_version)",
        SUBCELLULAR_LOCATION => "Description of the subcellular location of the mature protein (UniProt $uniprot_version)",
        SUBUNIT => "Description of the quaternary structure of a protein (UniProt $uniprot_version)",
        TISSUE_SPECIFICITY => "Description of the tissue-specific expression of mRNA or protein (UniProt $uniprot_version)",
        UCSC => "UCSC genome browser database accession number/primary key (UniProt $uniprot_version)",
        WEB_RESOURCE => "Links to related web resource(s) or database(s) (UniProt $uniprot_version)",
    };
}

sub run {
    my ($self, $bvfoa, $line_hash) = @_;
    my $bvf = $bvfoa->base_variation_feature;    
    my %bvfoa_info = %{get_bvfoa_info(@_)};

    my @DISEASES_PHENOTYPES;
        
    #get transcript annotation
    if(defined($bvfoa_info{enst})){
        my $uniprot = get_uniprot_from_vw($bvfoa_info{enst});
        map {$line_hash->{$_} = $uniprot->{$_} if $_ ne 'ID' && $_ ne 'ENST' && $_ ne 'ENSG' && $_ ne 'ENSP'} keys %{$uniprot};
        
        if(defined($uniprot->{DISEASE})){
            push(@DISEASES_PHENOTYPES, split(/\|/,$uniprot->{DISEASE}));
        }
    }

    #get overlapping features on protein
    if(defined($bvfoa_info{translation})
       && defined($bvfoa_info{enst})
       && defined($bvfoa_info{altered_aa_start})
       && defined($bvfoa_info{altered_aa_end})
       ){
        my @uniprot_feature_mutagen;
        my @uniprot_feature_variant;
        my @uniprot_feature_sites;
        my @uniprot_feature_other;
        my @uniprot_site_features = qw(
            ACT_SITE
            BINDING
            CA_BIND
            DISULFID
            DNA_BIND
            METAL
            NP_BIND
            SITE
            ZN_FING
        );
        my $uniprot_feature = get_uniprot_feature_from_vw($bvfoa_info{enst}, $bvfoa_info{altered_aa_start}, $bvfoa_info{altered_aa_end});
        #some keys get special handling
        foreach my $key (keys %{$uniprot_feature}){
            if ($key eq 'MUTAGEN'){
                push(@uniprot_feature_mutagen, join(':', @{$uniprot_feature->{$key}}));
            }
            if ($key eq 'VARIANT'){
                push(@uniprot_feature_variant, join(':', @{$uniprot_feature->{$key}}));
            }
            elsif(grep {$key eq $_} @uniprot_site_features){
                push(@uniprot_feature_sites, "$key:".join(':', @{$uniprot_feature->{$key}}));
            }
            else{
                push(@uniprot_feature_other, "$key:".join(':', @{$uniprot_feature->{$key}}));
            }
        }
        if ((defined(\@uniprot_feature_mutagen))
            && (scalar(@uniprot_feature_mutagen) > 0))
        {
            $line_hash->{MUTAGEN} = join( ';', @uniprot_feature_mutagen);
        }
        if ((defined(\@uniprot_feature_variant))
            && (scalar(@uniprot_feature_variant) > 0))
        {
            $line_hash->{VARIANT} = join( ';', @uniprot_feature_variant);
        }
        if ((defined(\@uniprot_feature_sites))
            && (scalar(@uniprot_feature_sites) > 0))
        {
            $line_hash->{SITES} = join( ';', @uniprot_feature_sites);
        }
        if ((defined(\@uniprot_feature_other))
            && (scalar(@uniprot_feature_other) > 0))
        {
            $line_hash->{OTHER_OVERLAPPING_FEATURES} = join( ';', @uniprot_feature_other);
        }
     } #if(defined($translation))

    #add phenotypes to the DISEASES_PHENOTYPES portmanteau column
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

sub get_uniprot_from_vw{
    
    my ($enst) = @_;
    
    my %data;
    
    my $query = "CALL $vw::vw_database.enst2uniprot('$enst')";
    my $qh = $vw::vw_conn->prepare($query);
    $qh->execute() or die "Unable to execute $query: $DBI::errstr\n";
    while (my @row = $qh->fetchrow_array()){
        #$data{topic} = value
        $data{$row[0]} = $row[1];
    }
   
   return \%data;
}

sub get_uniprot_feature_from_vw{
    
    my ($enst, $aaStart, $aaEnd) = @_;
    
    my %data;
    
    my $query = "CALL $vw::vw_database.enst2uniprot_feature('$enst', $aaStart, $aaEnd)";
    my $qh = $vw::vw_conn->prepare($query);
    $qh->execute() or die "Unable to execute $query: $DBI::errstr\n";
    while (my @row = $qh->fetchrow_array()){
        #$data{feature} = (aaStart,aaEnd,description)
        my @feature = @row[1..3];
        $data{$row[0]} = \@feature;
    }
   
   return \%data;
}


1;

