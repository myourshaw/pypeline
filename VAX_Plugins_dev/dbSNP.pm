=head1 LICENSE

 Copyright (c) 2011-2014 Michael Yourshaw.  All rights reserved.                                                                      

=head1 CONTACT                                                                                                       

 Michael Yourshaw <myourshaw@ucla.edu>
    
=cut

=head1 NAME

 dbSNP

=head1 SYNOPSIS

 mv dbSNP.pm ~/.vep/Plugins
 perl variant_effect_predictor.pl -i variations.vcf --plugin vw[,host,port,user,password,mysql,database]  --plugin dbSNP[,version]

=head1 DESCRIPTION

 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that
 adds the following new columns, derived from dbSNP:
  dbSNP138,
  SAO,
  SSR,
  PM,
  MUT,
  G5A,
  G5,
  CDA,
  CAF,
  COMMON,
  CLNSIG,
  CLNDBN,
  CNKMI,
 
 Requires that the VAX.pm module be in the Plugins directory

=head1 PARAMETERS

    dbsnp_version (default: dbSNP)

=cut

package dbSNP;
use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

#require Exporter;
#our @ISA = qw(Exporter);
#our @EXPORT_OK = qw();

use Bio::EnsEMBL::Variation::Utils::VEP qw(@OUTPUT_COLS);

use vw;
use VAX qw(get_bvfoa_info);

sub new {
    my $class = shift;
    my $self = $class->SUPER::new(@_);
    ($self->{dbsnp_version}) = @{$self->{params}} ? @{$self->{params}} : '';

    return $self;
}

sub version {
    return '75';
}

sub feature_types {
    return ['Transcript', 'RegulatoryFeature', 'MotifFeature', 'Gene', 'Exon', 'Intergenic'];
}

sub variant_feature_types {
    return ['VariationFeature', 'StructuralVariationFeature'];
}

sub get_header_info {
    my $self = shift;
    my $dbsnp_version = $self->{dbsnp_version};
    
    my %columns = (
        RS => "dbSNP ID (i.e. rs number) (dbSNP $dbsnp_version)",
        #RSPOS => "Chr position reported in dbSNP (dbSNP $dbsnp_version)",
        #RV => "RS orientation is reversed (dbSNP $dbsnp_version)",
        #VP => "Variation Property.  Documentation is at ftp://ftp.ncbi.nlm.nih.gov/snp/specs/dbSNP_BitField_latest.pdf (dbSNP $dbsnp_version)",
        #GENEINFO => "Pairs each of gene symbol:gene id.  The gene symbol and id are delimited by a colon (:) and each pair is delimited by a vertical bar (|) (dbSNP $dbsnp_version)",
        #dbSNPBuildID => "First dbSNP Build for RS (dbSNP $dbsnp_version)",
        SAO => "Variant Allele Origin: 0 - unspecified, 1 - Germline, 2 - Somatic, 3 - Both (dbSNP $dbsnp_version)",
        SSR => "Variant Suspect Reason Codes (may be more than one value added together) 0 - unspecified, 1 - Paralog, 2 - byEST, 4 - oldAlign, 8 - Para_EST, 16 - 1kg_failed, 1024 - other (dbSNP $dbsnp_version)",
        #WGT => "Weight, 00 - unmapped, 1 - weight 1, 2 - weight 2, 3 - weight 3 or more (dbSNP $dbsnp_version)",
        #VC => "Variation Class (dbSNP $dbsnp_version)",
        PM => "Variant is Precious(Clinical,Pubmed Cited) (dbSNP $dbsnp_version)",
        #TPA => "Provisional Third Party Annotation(TPA) (currently rs from PHARMGKB who will give phenotype data) (dbSNP $dbsnp_version)",
        #PMC => "Links exist to PubMed Central article (dbSNP $dbsnp_version)",
        #S3D => "Has 3D structure - SNP3D table (dbSNP $dbsnp_version)",
        #SLO => "Has SubmitterLinkOut - From SNP-SubSNP-Batch.link_out (dbSNP $dbsnp_version)",
        #NSF => "Has non-synonymous frameshift A coding region variation where one allele in the set changes all downstream amino acids. FxnClass = 44 (dbSNP $dbsnp_version)",
        #NSM => "Has non-synonymous missense A coding region variation where one allele in the set changes protein peptide. FxnClass = 42 (dbSNP $dbsnp_version)",
        #NSN => "Has non-synonymous nonsense A coding region variation where one allele in the set changes to STOP codon (TER). FxnClass = 41 (dbSNP $dbsnp_version)",
        #REF => "Has reference A coding region variation where one allele in the set is identical to the reference sequence. FxnCode = 8 (dbSNP $dbsnp_version)",
        #SYN => "Has synonymous A coding region variation where one allele in the set does not change the encoded amino acid. FxnCode = 3 (dbSNP $dbsnp_version)",
        #U3 => "In 3' UTR Location is in an untranslated region (UTR). FxnCode = 53 (dbSNP $dbsnp_version)",
        #U5 => "In 5' UTR Location is in an untranslated region (UTR). FxnCode = 55 (dbSNP $dbsnp_version)",
        #ASS => "In acceptor splice site FxnCode = 73 (dbSNP $dbsnp_version)",
        #DSS => "In donor splice-site FxnCode = 75 (dbSNP $dbsnp_version)",
        #INT => "In Intron FxnCode = 6 (dbSNP $dbsnp_version)",
        #R3 => "In 3' gene region FxnCode = 13 (dbSNP $dbsnp_version)",
        #R5 => "In 5' gene region FxnCode = 15 (dbSNP $dbsnp_version)",
        #OTH => "Has other variant with exactly the same set of mapped positions on NCBI refernce assembly. (dbSNP $dbsnp_version)",
        #CFL => "Has Assembly conflict. This is for weight 1 and 2 variant that maps to different chromosomes on different assemblies. (dbSNP $dbsnp_version)",
        #ASP => "Is Assembly specific. This is set if the variant only maps to one assembly (dbSNP $dbsnp_version)",
        MUT => "Is mutation (journal citation, explicit fact): a low frequency variation that is cited in journal and other reputable sources (dbSNP $dbsnp_version)",
        VLD => "Is Validated.  This bit is set if the variant has 2+ minor allele count based on frequency or genotype data. (dbSNP $dbsnp_version)",
        G5A => "5% minor allele frequency in each and all populations (dbSNP $dbsnp_version)",
        G5 => "5% minor allele frequency in 1+ populations (dbSNP $dbsnp_version)",
        #HD => "Marker is on high density genotyping kit (50K density or greater).  The variant may have phenotype associations present in dbGaP. (dbSNP $dbsnp_version)",
        #GNO => "Genotypes available. The variant has individual genotype (in SubInd table). (dbSNP $dbsnp_version)",
        #KGValidated => "1000 Genome validated (dbSNP $dbsnp_version)",
        #KGPhase1 => "1000 Genome phase 1 (incl. June Interim phase 1) (dbSNP $dbsnp_version)",
        #KGPilot123 => "1000 Genome discovery all pilots 2010(1,2,3) (dbSNP $dbsnp_version)",
        #KGPROD => "Has 1000 Genome submission (dbSNP $dbsnp_version)",
        #OTHERKG => "non-1000 Genome submission (dbSNP $dbsnp_version)",
        #PH3 => "HAP_MAP Phase 3 genotyped: filtered, non-redundant (dbSNP $dbsnp_version)",
        CDA => "Variation is interrogated in a clinical diagnostic assay (dbSNP $dbsnp_version)",
        #LSD => "Submitted from a locus-specific database (dbSNP $dbsnp_version)",
        #MTP => "Microattribution/third-party annotation(TPA:GWAS,PAGE) (dbSNP $dbsnp_version)",
        #OM => "Has OMIM/OMIA (dbSNP $dbsnp_version)",
        #NOC => "Contig allele not present in variant allele list. The reference sequence allele at the mapped position is not present in the variant allele list, adjusted for orientation. (dbSNP $dbsnp_version)",
        #WTD => "Is Withdrawn by submitter If one member ss is withdrawn by submitter, then this bit is set.  If all member ss' are withdrawn, then the rs is deleted to SNPHistory (dbSNP $dbsnp_version)",
        #NOV => "Rs cluster has non-overlapping allele sets. True when rs set has more than 2 alleles from different submissions and these sets share no alleles in common. (dbSNP $dbsnp_version)",
        CAF => "An ordered, comma delimited list of allele frequencies based on 1000Genomes, starting with the reference allele followed by alternate alleles as ordered in the ALT column. Where a 1000Genomes alternate allele is not in the dbSNPs alternate allele set, the allele is added to the ALT column.  The minor allele is the second largest value in the list, and was previuosly reported in VCF as the GMAF.  This is the GMAF reported on the RefSNP and EntrezSNP pages and VariationReporter (dbSNP $dbsnp_version)",
        COMMON => "RS is a common SNP.  A common SNP is one that has at least one 1000Genomes population with a minor allele of frequency = 1% and for which 2 or more founders contribute to that minor allele frequency. (dbSNP $dbsnp_version)",
        #CLNHGVS => "Variant names from HGVS.    The order of these variants corresponds to the order of the info in the other clinical  INFO tags. (dbSNP $dbsnp_version)",
        #CLNALLE => "Variant alleles from REF or ALT columns.  0 is REF, 1 is the first ALT allele, etc.  This is used to match alleles with other corresponding clinical (CLN) INFO tags.  A value of -1 indicates that no allele was found to match a corresponding HGVS allele name. (dbSNP $dbsnp_version)",
        #CLNSRC => "Variant Clinical Chanels (dbSNP $dbsnp_version)",
        CLNORIGIN => "Allele Origin. One or more of the following values may be added: 0 - unknown; 1 - germline; 2 - somatic; 4 - inherited; 8 - paternal; 16 - maternal; 32 - de-novo; 64 - biparental; 128 - uniparental; 256 - not-tested; 512 - tested-inconclusive; 1073741824 - other (dbSNP $dbsnp_version)",
        #CLNSRCID => "Variant Clinical Channel IDs (dbSNP $dbsnp_version)",
        CLNSIG => "Variant Clinical Significance, 0 - Uncertain significance, 1 - not provided, 2 - Benign, 3 - Likely benign, 4 - Likely pathogenic, 5 - Pathogenic, 6 - drug response, 7 - histocompatibility, 255 - other (dbSNP $dbsnp_version)",
        #CLNDSDB => "Variant disease database name (dbSNP $dbsnp_version)",
        #CLNDSDBID => "Variant disease database ID (dbSNP $dbsnp_version)",
        CLNDBN => "Variant disease name (dbSNP $dbsnp_version)",
        #CLNACC => "Variant Accession and Versions (dbSNP $dbsnp_version)",
        CLINVAR => "Vatiant is in ClinVar (dbSNP $dbsnp_version)",
        CNKMI => "Common variant with no known medical impact (dbSNP $dbsnp_version)",
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
        my $query = "CALL $vw::vw_database.get_dbsnp('$chrom',$pos,'$ref','$alt')";
        my $qh = $vw::vw_conn->prepare($query);
        $qh->execute() or die "Unable to execute $query: $DBI::errstr\n";
        my $hash_ref = $qh->fetchrow_hashref();
        foreach my $c(@output_cols){
            $line_hash->{$c} = defined($hash_ref->{$c}) ? $hash_ref->{$c} : '';
        }

        $query = "CALL $vw::vw_database.get_dbsnp_clinvar('$chrom',$pos,'$ref','$alt')";
        $qh = $vw::vw_conn->prepare($query);
        $qh->execute() or die "Unable to execute $query: $DBI::errstr\n";
        my $numrows = $qh->rows;
            $line_hash->{'CLINVAR'} = $numrows > 0 ? 1 : 0;

        $query = "CALL $vw::vw_database.get_dbsnp_common_no_known_medical_impact('$chrom',$pos,'$ref','$alt')";
        $qh = $vw::vw_conn->prepare($query);
        $qh->execute() or die "Unable to execute $query: $DBI::errstr\n";
        #$numrows = $qh->numrows();
        $numrows = $qh->rows;
            $line_hash->{'CNKMI'} = $numrows > 0 ? 1 : 0;
    }
       
    return {};
}


1;
