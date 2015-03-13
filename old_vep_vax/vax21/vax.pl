#!/usr/bin/perl

#--vcf_out --buffer_size 5000 -input_file /Volumes/storage/vw/snp1_test.vcf --format vcf -output_file /Volumes/storage/vw/snp1_test.vcf.vep --force_overwrite --terms SO --sift=b --polyphen=b --condel=b --regulatory --hgnc --hgvs --protein --gene --check_ref --host myourshaw-dev.genome.ucla.edu --port=3306 --user ensembl --password ensembl --no_progress
#--vcf_out --buffer_size 5000 -input_file /Volumes/storage/vw/indel1_test.vcf --format vcf -output_file /Volumes/storage/vw/indel1_test.vcf.vep --force_overwrite --terms SO --sift=b --polyphen=b --condel=b --regulatory --hgnc --hgvs --protein --gene --check_ref --host myourshaw-dev.genome.ucla.edu --port=3306 --user ensembl --password ensembl --no_progress
#--vcf_out --buffer_size 5000 -input_file /Volumes/storage/vw/snp_test.vcf --format vcf -output_file /Volumes/storage/vw/snp_test.vcf.vep --force_overwrite --terms SO --sift=b --polyphen=b --condel=b --regulatory --hgnc --hgvs --protein --gene --check_ref --host myourshaw-dev.genome.ucla.edu --port=3306 --user ensembl --password ensembl --no_progress
#--vcf_out --buffer_size 5000 -input_file /Volumes/storage/vw/indel_test.vcf --format vcf -output_file /Volumes/storage/vw/indel_test.vcf.vep --force_overwrite --terms SO --sift=b --polyphen=b --condel=b --regulatory --hgnc --hgvs --protein --gene --check_ref --host myourshaw-dev.genome.ucla.edu --port=3306 --user ensembl --password ensembl --no_progress

#/Volumes/storage/vw/gmd_splice_site_test.vcf
#--vcf_out --buffer_size 1 -input_file /Volumes/storage/vw/gmd_splice_site_test.vcf --format vcf -output_file /Volumes/storage/vw/gmd_splice_site_test.vcf.vep --force_overwrite --terms SO --sift=b --polyphen=b --condel=b --regulatory --hgnc --hgvs --protein --gene --check_ref --host myourshaw-dev.genome.ucla.edu --port=3306 --user ensembl --password ensembl --no_progress

#--vcf_out --buffer_size 5000 -input_file /Volumes/storage/vw/gmd_distinct_variant.vcf --format vcf -output_file /Volumes/storage/vw/gmd_distinct_variant.vcf.vep --force_overwrite --terms SO --sift=b --polyphen=b --condel=b --regulatory --hgnc --hgvs --protein --gene --check_ref --host myourshaw-dev.genome.ucla.edu --port=3306 --user ensembl --password ensembl --no_progress

#cluster
#--vcf_out --buffer_size 1000 -input_file /scratch0/tmp/myourshaw/gmd/vcfs/tmp_GL0Q4I/gmd32.analysis_ready.vcf.part.00000000 --format vcf -output_file /scratch0/tmp/myourshaw/gmd/vcfs/tmp_GL0Q4I/gmd32.analysis_ready.vcf.part.00000000.vax --force_overwrite --terms SO --sift=b --polyphen=b --condel=b --regulatory --hgnc --hgvs --protein --gene --check_ref --host cortex.local --port=3306 --user ensembl --password ensembl --species human --no_progress --vw_database vw --vw_host myourshaw-dev.genome.ucla.edu --vw_password vw --vw_platform mysql --vw_port 3306 --vw_user vw


=head1 LICENSE

  Copyright (c) 1999-2011 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=head1 NAME

Variant Effect Predictor - a script to predict the consequences of genomic variants

http://www.ensembl.org/info/docs/variation/vep/vep_script.html

Version 2.1

by Will McLaren (wm2@ebi.ac.uk)
=cut

use strict;
use Getopt::Long;
use FileHandle;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor;
use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(MAX_DISTANCE_FROM_TRANSCRIPT);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::Variation::Utils::Sequence qw(unambiguity_code);
use Storable qw(nstore_fd fd_retrieve);

# we need to manually include all these modules for caching to work
use Bio::EnsEMBL::CoordSystem;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::DBSQL::GeneAdaptor;
use Bio::EnsEMBL::DBSQL::SliceAdaptor;
use Bio::EnsEMBL::DBSQL::TranslationAdaptor;
use Bio::EnsEMBL::DBSQL::TranscriptAdaptor;
use Bio::EnsEMBL::DBSQL::MetaContainer;
use Bio::EnsEMBL::DBSQL::CoordSystemAdaptor;

#MY
use Bio::Matrix::IO;
use List::Util qw(first);
use LWP::Simple;
use HTML::TokeParser;
use SWISS::Entry;	    	
use SWISS::CCs;
use SWISS::DEs;
use File::Basename;
use Tie::IxHash;
use Memoize;
use Memoize::ExpireLRU;
use DBI;
use DBD::mysql;
#MY

#MY
#number of recent results to cache using memoize.
#cached results prevent unnecessary data retrieval
#or calculations for adjacent input SNPs
#improves throughput ~4X
my $cachesize = 50;

#cache vw searches
tie my %get_uniprot_from_vw_cache => 'Memoize::ExpireLRU',
    CACHESIZE         => $cachesize;
memoize 'get_uniprot_from_vw',
    SCALAR_CACHE => [ HASH => \%get_uniprot_from_vw_cache ];

sub get_uniprot_from_vw_normalizer {
    my %hash = (@_);
    join( ',', $hash{enst} );
}

tie my %get_uniprot_feature_from_vw_cache => 'Memoize::ExpireLRU',
    CACHESIZE         => $cachesize;
memoize 'get_uniprot_feature_from_vw',
    SCALAR_CACHE => [ HASH => \%get_uniprot_feature_from_vw_cache ];

sub get_uniprot_feature_from_vw_normalizer {
    my %hash = (@_);
    join( ',', $hash{enst}, $hash{aaStart}, $hash{aaEnd} );
}

tie my %get_kegg_from_vw_cache => 'Memoize::ExpireLRU',
    CACHESIZE         => $cachesize;
memoize 'get_kegg_from_vw',
    SCALAR_CACHE => [ HASH => \%get_kegg_from_vw_cache ];

sub get_kegg_from_vw_normalizer {
    my %hash = (@_);
    join( ',', $hash{ensg} );
}

tie my %get_omim_from_vw_cache => 'Memoize::ExpireLRU',
    CACHESIZE         => $cachesize;
memoize 'get_omim_from_vw',
    SCALAR_CACHE => [ HASH => \%get_omim_from_vw_cache ];

sub get_omim_from_vw_normalizer {
    my %hash = (@_);
    join( ',', $hash{ensg} );
}

#cache uniprot searches
tie my %uniprot_cache => 'Memoize::ExpireLRU',
    CACHESIZE         => $cachesize;
memoize 'get_uniprot_record',
    SCALAR_CACHE => [ HASH => \%uniprot_cache ],
    NORMALIZER   => 'get_uniprot_record_normalizer';

sub get_uniprot_record_normalizer {
    my %hash = (@_);
    join( ',', $hash{id} );
}

#cache ncbi searches
tie my %ncbi_cache => 'Memoize::ExpireLRU',
    CACHESIZE      => $cachesize;
memoize 'get_ncbi_record',
    SCALAR_CACHE => [ HASH => \%ncbi_cache ],
    NORMALIZER   => 'get_ncbi_record_normalizer';

sub get_ncbi_record_normalizer {
    my %hash = (@_);
    join( ',', $hash{query} );
}

#cache omes calculation
tie my %omes_cache => 'Memoize::ExpireLRU',
    CACHESIZE      => $cachesize;
memoize 'calculate_omes',
    SCALAR_CACHE => [ HASH => \%omes_cache ],
    NORMALIZER   => 'calculate_omes_normalizer';

sub calculate_omes_normalizer {
    my %hash = (@_);
    join( ',',
        $hash{position_of_interest},
        $hash{protein_id}, $hash{comparison_type} );
}

#cache uniprot annotations
tie my %uniprot_annotations_cache => 'Memoize::ExpireLRU',
    CACHESIZE         => $cachesize;
memoize 'get_uniprot_annotations',
    SCALAR_CACHE => [ HASH => \%uniprot_annotations_cache ];

#cache ortholog stuff
tie my %get_model_phenotypes_from_model_genomic_coords_cache =>'Memoize::ExpireLRU',
    CACHESIZE         => $cachesize;
memoize 'get_model_phenotypes_from_model_genomic_coords',
    SCALAR_CACHE => [ HASH => \%get_model_phenotypes_from_model_genomic_coords_cache ];
tie my %get_model_genomic_coords_of_model_amino_acid_cache =>'Memoize::ExpireLRU',
    CACHESIZE         => $cachesize;
memoize 'get_model_genomic_coords_of_model_amino_acid',
    SCALAR_CACHE => [ HASH => \%get_model_genomic_coords_of_model_amino_acid_cache ];
tie my %determine_context_alignments_cache =>'Memoize::ExpireLRU',
    CACHESIZE         => $cachesize;
memoize 'determine_context_alignments',
    SCALAR_CACHE => [ HASH => \%determine_context_alignments_cache ];
tie my %determine_aligned_homolog_residues_cache =>'Memoize::ExpireLRU',
    CACHESIZE         => $cachesize;
memoize 'determine_aligned_homolog_residues',
    SCALAR_CACHE => [ HASH => \%determine_aligned_homolog_residues_cache ];
#MY

# debug
#use Time::HiRes qw(tv_interval gettimeofday);

# output columns
my @OUTPUT_COLS = qw(
    Uploaded_variation
    Location
    Allele
    Gene
    Feature
    Feature_type
    Consequence
    Consequence_rank
    Consequences_all
    cDNA_position
    CDS_position
    Protein_position
    Amino_acids
    Codons
    Existing_variation
    HGNC
    ENSP
    HGVSc
    HGVSp
    SIFT_prediction
    SIFT_score
    PolyPhen_prediction
    PolyPhen_score
    Condel_prediction
    Condel_score
    MATRIX
    HIGH_INF_POS
    Alignment_Score_Change
    C_blosum
    Context_Conservation
    Amino_Acids_In_Orthologues
    Orthologue_Species
    Reference_Splice_Site
    Variant_Splice_Site
    Protein_Length
    Protein_Length_Decrease
    Protein_Sequence_Lost
    Protein_Length_Increase
    Protein_Sequence_Gained
    Gene_Description
    Entrez_Gene_Name
    Entrez_Gene_ID
    Overlapping_Protein_Domains
    Gene_Ontology
    Phenotypes_Position
    Phenotypes_Gene
    HGMD
    OMIM_Disorder
    KEGG_Pathway
    UniProt_ID
    UniProtKB_AC
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
    GO
    GO_term
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
    Overlapping_Protein_Features_Orthologues
    Protein_OMES_Orthologues
    Protein_OMES_Family
    Canonical_Transcript
);

# global vars
my $VERSION = '2.1';

# set output autoflush for progress bars
$| = 1;

# configure from command line opts
my $config = &configure(scalar @ARGV);

#MY
my @vcf_header;
#MY

#ZILA
# dynamically get species common and latin names
my ($speciesListCommon_ptr,$speciesListLatin_ptr) = &getEnsemblSpecies();

# generate distance array from phylogenetic tree
my $distHash_ptr =
    &generatePhyloTreeDistArray($config->{species}, $speciesListCommon_ptr, $speciesListLatin_ptr);
my %distHash = %{$distHash_ptr};
#ZILA

# run the main sub routine
&main($config);

# this is the main sub-routine - it needs the configured $config hash
sub main {
    my $config = shift;
    
    debug("Starting...") unless defined $config->{quiet};
    
    my ($include_regions, $transcript_cache);
    
    # scan file if requested
    $include_regions = &scan_file($config) if defined($config->{scan});
    
    # build transcript cache upfront if requested
    $transcript_cache = &cache_transcripts($config, $include_regions) if defined($config->{upfront});
    
    # create a hash to hold slices so we don't get the same one twice
    my %slice_cache = ();
    
    # load slices from the transcript cache if we have it
    # saves us fetching them again
    %slice_cache = %{&build_slice_cache($config, $transcript_cache)} if defined($transcript_cache);
    
    my @new_vfs;
    my %vf_hash;
#MY
    my @new_vcfs;
    my %vcf_hash;
#MY
    
    my $line_number = 0;
    my ($vf_count, $total_vf_count);
    my $in_file_handle = $config->{in_file_handle};
    
    # read the file
    while(<$in_file_handle>) {
        chomp;
        
        $line_number++;
      
#MY
        #VCF header
        if ($config->{vcf_out} && /^\#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO/){
            @vcf_header = (split /\t/, $_);
            $vcf_header[0] =~ s/^\#//;
            @OUTPUT_COLS = (@vcf_header, @OUTPUT_COLS);
            $config->{out_file_handle} = &get_out_file_handle($config);
        }
#MY
        # header line?
        next if /^\#/;
      
#MY
        #VCF
        my %vcf_record;
        if ($config->{vcf_out} && (($config->{format} =~ /vcf/i) || scalar(@vcf_header)>0)){
            my @vcf_fields = (split /\s+/, $_);
            map {$vcf_record{$vcf_header[$_]} = $vcf_fields[$_]} 0..$#vcf_fields;
        }
#MY     
        # some lines (pileup) may actually parse out into more than one variant
        foreach my $sub_line(@{&parse_line($config, $_)}) {
            # get the sub-line into named variables
            my ($chr, $start, $end, $allele_string, $strand, $var_name) = @{$sub_line};
            
            next if defined($config->{chr}) && !$config->{chr}->{$chr};
            
            # non-variant line from VCF
            next if $chr eq 'non-variant';
            
            # fix inputs
            $chr =~ s/chr//ig unless $chr =~ /^chromosome$/i;
            $chr = 'MT' if $chr eq 'M';
            $strand = ($strand =~ /\-/ ? "-1" : "1");
            $allele_string =~ tr/acgt/ACGT/;
            
            # sanity checks
            unless($start =~ /^\d+$/ && $end =~ /^\d+$/) {
              warn("WARNING: Start $start or end $end coordinate invalid on line $line_number\n") unless defined $config->{quiet};
              next;
            }
            
            unless($allele_string =~ /([ACGT-]+\/*)+/) {
              warn("WARNING: Invalid allele string $allele_string on line $line_number\n") unless defined $config->{quiet};
              next;
            }
            
            # now get the slice
            my $slice;
            
            # don't get slices if we're using cache
            # we can steal them from transcript objects later
            if((!defined($config->{cache}) && !defined($config->{whole_genome})) || defined($config->{check_ref})) {
                
                # check if we have fetched this slice already
                if(defined $slice_cache{$chr}) {
                    $slice = $slice_cache{$chr};
                }
                
                # if not create a new one
                else {
                    
                    $slice = &get_slice($config, $chr);
                    
                    # if failed, warn and skip this line
                    if(!defined($slice)) {
                        warn("WARNING: Could not fetch slice named $chr on line $line_number\n") unless defined $config->{quiet};
                        next;
                    }    
                    
                    # store the hash
                    $slice_cache{$chr} = $slice;
                }
            } #if((!defined($config->{cache}) && !defined($config->{whole_genome})) || defined($config->{check_ref}))
            
            # check reference allele if requested
            if(defined $config->{check_ref}) {
                my $ref_allele = (split /\//, $allele_string)[0];
                
                my $ok = 0;
                my $slice_ref_allele;
                
                # insertion, therefore no ref allele to check
                if($ref_allele eq '-') {
                    $ok = 1;
                }
                else {
                    my $slice_ref = $slice->sub_Slice($start, $end, $strand);
                    
                    if(!defined($slice_ref)) {
                        warn "WARNING: Could not fetch sub-slice from $start\-$end\($strand\) on line $line_number" unless defined $config->{quiet};
                    }
                    
                    else {
                        $slice_ref_allele = $slice_ref->seq;
                        $ok = ($slice_ref_allele eq $ref_allele ? 1 : 0);
                    }
                }
                
                if(!$ok) {
                    warn
                        "WARNING: Specified reference allele $ref_allele ",
                        "does not match Ensembl reference allele",
                        ($slice_ref_allele ? " $slice_ref_allele" : ""),
                        " on line $line_number" unless defined $config->{quiet};
                    next;
                }
            } # if(defined $config->{check_ref})
           
            # create a new VariationFeature object
            my $new_vf = Bio::EnsEMBL::Variation::VariationFeature->new(
              -start => $start,
              -end => $end,
              -slice => $slice,           # the variation must be attached to a slice
              -allele_string => $allele_string,
              -strand => $strand,
              -map_weight => 1,
              -adaptor => $config->{vfa},           # we must attach a variation feature adaptor
              -variation_name => (defined $var_name ? $var_name : $chr.'_'.$start.'_'.$allele_string),
            );
            
            if(defined $config->{whole_genome}) {
                push @{$vf_hash{$chr}{int($start / $config->{chunk_size})}{$start}}, $new_vf;
                $vf_count++;
                $total_vf_count++;
#MY
                push @{$vcf_hash{$chr}{int($start / $config->{chunk_size})}{$start}}, \%vcf_record;
#MY
                
                if($vf_count == $config->{buffer_size}) {
                    debug("Read $vf_count variants into buffer") unless defined($config->{quiet});
                    
                    $include_regions ||= &regions_from_hash($config, \%vf_hash);
                    
                    &check_existing_hash($config, \%vf_hash) if defined($config->{check_existing});
#MY
                    #&whole_genome_fetch($config, \%vf_hash, $transcript_cache, $include_regions);
                    &whole_genome_fetch($config, \%vf_hash, $transcript_cache, $include_regions, \%vcf_hash);
#MY
                    
                    debug("Processed $total_vf_count total variants") unless defined($config->{quiet});
                    
                    undef $include_regions unless defined($config->{scan});
                    %vf_hash = ();
                    $vf_count = 0;
#MY
                    %vcf_hash = ();
#MY
                }
            } #if(defined $config->{whole_genome})
            else {
#MY
                #&print_consequences($config, [$new_vf]);
                &print_consequences($config, [$new_vf], [\%vcf_record]);
#MY
                $vf_count++;
                $total_vf_count++;
                debug("Processed $vf_count variants") if $vf_count =~ /0$/ && defined($config->{verbose});
            } #else
        } #foreach my $sub_line(@{&parse_line($config, $_)})
    } #while(<$in_file_handle>)
    
    # if in whole-genome mode, finish off the rest of the buffer
    if(defined $config->{whole_genome} && %vf_hash) {
        debug("Read $vf_count variants into buffer") unless defined($config->{quiet});
        $include_regions ||= &regions_from_hash($config, \%vf_hash);
        &check_existing_hash($config, \%vf_hash) if defined($config->{check_existing});
#MY
        #&whole_genome_fetch($config, \%vf_hash, $transcript_cache, $include_regions);
        &whole_genome_fetch($config, \%vf_hash, $transcript_cache, $include_regions, \%vcf_hash);
#MY
    } #if(defined $config->{whole_genome} && %vf_hash)
    
#MY
    $config->{vw_conn}->disconnect();
#MY
    
    debug("Executed ", defined($Bio::EnsEMBL::DBSQL::StatementHandle::count_queries) ? $Bio::EnsEMBL::DBSQL::StatementHandle::count_queries : 'unknown number of', " SQL statements") if defined($config->{count_queries}) && !defined($config->{quiet});
    
    debug("Finished!") unless defined $config->{quiet};
} #main

# takes a listref of variation features and prints out consequence information
sub print_consequences {
    my $config = shift;
    my $vfs = shift;
#MY
    my $vcfs = shift;
#MY
   
    my $out_file_handle = $config->{out_file_handle};
    
    # method name for consequence terms
    my $term_method = $config->{terms}.'_term';
    
    my ($vf_count, $vf_counter);
    $vf_count = scalar @$vfs;

#MY
    #foreach my $new_vf(@$vfs) {
    for(my $i=0;$i<scalar(@$vfs);$i++){
        my $new_vf = @$vfs[$i];
#MY
        &progress($config, $vf_counter++, $vf_count) unless $vf_count == 1;
        
        # find any co-located existing VFs
        my $existing_vf = $new_vf->{existing};
        $existing_vf ||= &find_existing($config, $new_vf) if defined $config->{check_existing};
        
        # initiate line hash for this variation
        my $line = {
            Uploaded_variation  => $new_vf->variation_name,
            Location            => $new_vf->seq_region_name.':'.&format_coords($new_vf->start, $new_vf->end),
            Existing_variation  => $existing_vf,
            Extra               => {},
        };
        
#MY
        my $vcf = @$vcfs[$i];
        if($config->{vcf_out} && scalar(@vcf_header)>0){
            map {$line->{$_} = $vcf->{$_}} keys %{$vcf};
        }
#MY
        # force empty hash into object's transcript_variations if undefined from whole_genome_fetch
        # this will stop the API trying to go off and fill it again
        $new_vf->{transcript_variations} ||= {} if defined $config->{whole_genome};
#MY
        &my_vf_annotations($config, $line, $new_vf);
#MY
        
        # regulatory stuff
        if(!defined $config->{coding_only} && defined $config->{regulatory}) {
            
            for my $rfv (@{ $new_vf->get_all_RegulatoryFeatureVariations }) {
                
                my $rf = $rfv->regulatory_feature;
                
                $line->{Feature_type}   = 'RegulatoryFeature';
                $line->{Feature}        = $rf->stable_id;
                
                # this currently always returns 'RegulatoryFeature', so we ignore it for now
                #$line->{Extra}->{REG_FEAT_TYPE} = $rf->feature_type->name;
                
#MY
                &my_rf_annotations($config, $line, $rf);
#MY
#MY
                my_rfv_annotations($config, $line, $rfv);
#MY
                for my $rfva (@{ $rfv->get_all_alternate_RegulatoryFeatureVariationAlleles }) {
#MY
                    $line->{Allele}         = $rfva->variation_feature_seq;
                    my @consequences = map {$_->$term_method || $_->SO_term}
                        @{$rfva->get_all_OverlapConsequences};
                    my @ranks = map {$_->rank} @{$rfva->get_all_OverlapConsequences};
                    my %consequences_ranks;
                    map {$consequences_ranks{$consequences[$_]} = $ranks[$_]} 0..$#ranks;
                    my @consequences_sorted = sort {$consequences_ranks{$a} <=> $consequences_ranks{$b}} keys %consequences_ranks;
                    my $consequence_severest = $consequences_sorted[0];
                    my $consequence_severest_rank = $consequences_ranks{$consequence_severest};
                    $line->{Consequence} = $consequence_severest;
                    $line->{Consequence_rank} = $consequence_severest_rank;
                    $line->{Consequences_all} = join ",", @consequences_sorted;
                    
                    #$line->{Consequence}    = join ',', 
                    #    map { $_->$term_method || $_->display_term } 
                    #        @{ $rfva->get_all_OverlapConsequences };
                            
                    &my_rfva_annotations($config, $line, $rfva);
#MY
                    print_line($line);
                }
            }
            
            for my $mfv (@{ $new_vf->get_all_MotifFeatureVariations }) {
                
                my $mf = $mfv->motif_feature;
                
                $line->{Feature_type}   = 'MotifFeature';
                $line->{Feature}        = $mf->binding_matrix->name;
               
#MY
                &my_mf_annotations($config, $line, $mf);
#MY
#MY
                &my_mfv_annotations($config, $line, $mfv);
#MY
                for my $mfva (@{ $mfv->get_all_alternate_MotifFeatureVariationAlleles }) {
                   
#MY
                    #$line->{Extra}->{MATRIX} = $mf->binding_matrix->description.'_'.$mf->display_label,
                    $line->{MATRIX} = $mf->binding_matrix->description.'_'.$mf->display_label,
                    #$line->{Extra}->{MATRIX} =~ s/\s+/\_/g;
                    $line->{MATRIX} =~ s/\s+/\_/g;
#MY

                    my $high_inf_pos = $mfva->in_informative_position;

                    if (defined $high_inf_pos) {
#MY
                        #$line->{Extra}->{HIGH_INF_POS}  = ($high_inf_pos ? 'Y' : 'N');
                        $line->{HIGH_INF_POS}  = ($high_inf_pos ? 'Y' : 'N');
#MY
                    }
                    
                    $line->{Allele}         = $mfva->variation_feature_seq;
#MY
                    my @consequences = map {$_->$term_method || $_->SO_term}
                        @{$mfva->get_all_OverlapConsequences};
                    my @ranks = map {$_->rank} @{$mfva->get_all_OverlapConsequences};
                    my %consequences_ranks;
                    map {$consequences_ranks{$consequences[$_]} = $ranks[$_]} 0..$#ranks;
                    my @consequences_sorted = sort {$consequences_ranks{$a} <=> $consequences_ranks{$b}} keys %consequences_ranks;
                    my $consequence_severest = $consequences_sorted[0];
                    my $consequence_severest_rank = $consequences_ranks{$consequence_severest};
                    $line->{Consequence} = $consequence_severest;
                    $line->{Consequence_rank} = $consequence_severest_rank;
                    $line->{Consequences_all} = join ",", @consequences_sorted;

                    #$line->{Consequence}    = join ',', 
                    #    map { $_->$term_method || $_->display_term } 
                    #        @{ $mfva->get_all_OverlapConsequences };
                            
                    &my_mfva_annotations($config, $line, $mfva);
#MY
                    print_line($line);
                }
            }
        }
        
        
        # get TVs
        my $tvs = $new_vf->get_all_TranscriptVariations;
        
        # no TVs (intergenic) or only most severe
        if(!@$tvs || defined($config->{most_severe}) || defined($config->{summary})) {
            if(defined($config->{summary})) {
                $line->{Consequence} = join ",", @{$new_vf->consequence_type($config->{terms}) || $new_vf->consequence_type};
            }
            else {
                $line->{Consequence} = $new_vf->display_consequence($config->{terms}) || $new_vf->display_consequence;
            }
            
            &print_line($line);
        }
        
        else {
            foreach my $tv(@$tvs) {
                
                next if(defined $config->{coding_only} && !($tv->affects_transcript));
                
                my $t = $tv->transcript;
                
                $line->{Feature_type}       = 'Transcript';
                $line->{Feature}            = $t->stable_id if defined $t;
                $line->{cDNA_position}      = &format_coords($tv->cdna_start, $tv->cdna_end);
                $line->{CDS_position}       = &format_coords($tv->cds_start, $tv->cds_end);
                $line->{Protein_position}   = &format_coords($tv->translation_start, $tv->translation_end);
                
                # get gene
                my $gene;
                
                if(defined($config->{gene})) {
                    $line->{Gene} = $tv->transcript->{_gene_stable_id};
                    
                    if(!defined($line->{Gene})) {
                        $gene = $config->{ga}->fetch_by_transcript_stable_id($t->stable_id);
                        $line->{Gene}= $gene->stable_id;
                    }
                }
                
#MY
                &my_tv_annotations($config, $line, $tv);
#MY
                foreach my $tva(@{$tv->get_all_alternate_TranscriptVariationAlleles}) {
                    
                    # basic stuff
                    $line->{Allele}         = $tva->variation_feature_seq;
                    $line->{Amino_acids}    = $tva->pep_allele_string;
                    $line->{Codons}         = $tva->display_codon_allele_string;
#MY
                    #$line->{Consequence}    = join ",", map {$_->$term_method || $_->display_term} @{$tva->get_all_OverlapConsequences};
#MY                    
                    # HGNC
                    if(defined $config->{hgnc}) {
                        my $hgnc;
                        $hgnc = $tv->transcript->{_gene_hgnc};
                        
                        if(!defined($hgnc)) {
                            if(!defined($gene)) {
                                $gene = $config->{ga}->fetch_by_transcript_stable_id($tv->transcript->stable_id);
                            }
                            
                            my @entries = grep {$_->database eq 'HGNC'} @{$gene->get_all_DBEntries()};
                            if(scalar @entries) {
                                $hgnc = $entries[0]->display_id;
                            }
                        }
                        
                        $hgnc = undef if $hgnc eq '-';
#MY                        
                        #$line->{Extra}->{HGNC} = $hgnc if defined($hgnc);
                        $line->{HGNC} = $hgnc if defined($hgnc);
#MY
                    }
                    
                    # protein ID
                    if(defined $config->{protein} && $t->translation) {
#MY
                        #$line->{Extra}->{ENSP} = $t->translation->stable_id;
                        $line->{ENSP} = $t->translation->stable_id;
#MY
                    }
                    
                    # HGVS
                    if(defined $config->{hgvs}) {
#MY
                        #$line->{Extra}->{HGVSc} = $tva->hgvs_coding if defined($tva->hgvs_coding);
                        $line->{HGVSc} = $tva->hgvs_coding if defined($tva->hgvs_coding);
                        #$line->{Extra}->{HGVSp} = $tva->hgvs_protein if defined($tva->hgvs_protein);
                        $line->{HGVSp} = $tva->hgvs_protein if defined($tva->hgvs_protein);
#MY
                    }
                    
                    foreach my $tool (qw(SIFT PolyPhen Condel)) {
                        my $lc_tool = lc($tool);
                        
                        if (my $opt = $config->{$lc_tool}) {
                            my $want_pred   = $opt =~ /^p/i;
                            my $want_score  = $opt =~ /^s/i;
                            my $want_both   = $opt =~ /^b/i;
                            
                            if ($want_both) {
                                $want_pred  = 1;
                                $want_score = 1;
                            }
                            
                            next unless $want_pred || $want_score;
                            
                            my $pred_meth   = $lc_tool.'_prediction';
                            my $score_meth  = $lc_tool.'_score';
#MY
                            my $pred_col   = $tool.'_prediction';
                            my $score_col  = $tool.'_score';
#MY
                            
                            my $pred = $tva->$pred_meth;
                            
                            if($pred) {
                                
                                if ($want_pred) {
                                    $pred =~ s/\s+/\_/;
#MY
                                    #$line->{Extra}->{$tool} = $pred;
                                    $line->{$pred_col} = $pred;
#MY
                                }
                                    
                                if ($want_score) {
                                    my $score = $tva->$score_meth;
                                    
                                    if(defined $score) {
                                        if($want_pred) {
#MY
                                            #$line->{Extra}->{$tool} .= "($score)";
                                            $line->{$score_col} = $score;
#MY
                                        }
                                        else {
#MY
                                            #$line->{Extra}->{$tool} = $score;
                                            $line->{$score_col} = $score;
#MY
                                        }
                                    }
                                }
                            }
                        }
                    } #foreach my $tool (qw(SIFT PolyPhen Condel))
#MY
                    &my_tva_annotations($config, $line, $tv, $tva);
#MY
                    &print_line($line);
                } #foreach my $tva(@{$tv->get_all_alternate_TranscriptVariationAlleles})
            } #foreach my $tv(@$tvs)
        } #
    } #for(my $i=0;$i<scalar(@$vfs);$i++) was foreach my $new_vf(@$vfs)
    
    &end_progress($config) unless $vf_count == 1;
} #sub print_consequences

# prints a line from the hash
sub print_line {
    my $line = shift;

#MY
    #$line->{Extra} = join ';', map { $_.'='.$line->{Extra}->{$_} } keys %{ $line->{Extra} || {} };
    #my $output = join "\t", map { $line->{$_} || '-' } @OUTPUT_COLS;
    my $output = join "\t", map { $line->{$_} || '.' } @OUTPUT_COLS;
#MY

    my $fh = $config->{out_file_handle};

    print $fh "$output\n";

    # clear out the Extra column for the next line
    $line->{Extra} = {};
} #print_line

# sets up configuration hash that is used throughout the script
sub configure {
    my $args = shift;
    
    my $config = {};
    
    GetOptions(
        $config,
        'help',                    # displays help message
        
        # input options,
        'config=s',                # config file name
        'input_file=s',            # input file name
        'format=s',                # input file format
        
        # DB options
        'species=s',               # species e.g. human, homo_sapiens
        'registry=s',              # registry file
        'host=s',                  # database host
        'port=s',                  # database port
        'user=s',                  # database user name
        'password=s',              # database password
        'db_version=i',            # Ensembl database version to use e.g. 62
        'genomes',                 # automatically sets DB params for e!Genomes
        #'no_disconnect',           # disables disconnect_when_inactive
        
        # runtime options
        'most_severe',             # only return most severe consequence
        'summary',                 # only return one line per variation with all consquence types
        'buffer_size=i',           # number of variations to read in before analysis
        'chunk_size=s',            # size in bases of "chunks" used in internal hash structure
        'check_ref',               # check supplied reference allele against DB
        'check_existing',          # find existing co-located variations
        'check_alleles',           # only attribute co-located if alleles are the same
        'failed=i',                # include failed variations when finding existing
        'no_whole_genome',         # disables now default whole-genome mode
        'whole_genome',            # proxy for whole genome mode - now just warns user
        'gp',                      # read coords from GP part of INFO column in VCF (probably only relevant to 1KG)
        'chr=s',                   # analyse only these chromosomes, e.g. 1-5,10,MT
        
        # verbosity options
        'verbose',                 # print out a bit more info while running
        'quiet',                   # print nothing to STDOUT (unless using -o stdout)
        'no_progress',             # don't display progress bars
        
        # output options
        'output_file=s',           # output file name
        'force_overwrite',         # force overwrite of output file if already exists
        'terms=s',                 # consequence terms to use e.g. NCBI, SO
        'coding_only',             # only return results for consequences in coding regions
        'protein',                 # add e! protein ID to extra column
        'hgnc',                    # add HGNC gene ID to extra column
        'hgvs',                    # add HGVS names to extra column
        'sift=s',                  # SIFT predictions
        'polyphen=s',              # PolyPhen predictions
        'condel=s',                # Condel predictions
        'gene',                    # force gene column to be populated (disabled by default, enabled when using cache)
        'regulatory',              # enable regulatory stuff
        
        # cache stuff
        'cache',                   # use cache
        'write_cache',             # enables writing to the cache
        'build=s',                 # builds cache from DB from scratch; arg is either all (all top-level seqs) or a list of chrs
        'scan',                    # scan the whole input file at the beginning to get regions
        'upfront',                 # fetch transcripts and prefetch upfront before analysis starts (requires scan)
        'prefetch',                # prefetch exons, translation, introns, codon table etc for each transcript
        'strip',                   # strips adaptors etc from objects before caching them
        'rebuild=s',               # rebuilds cache by reading in existing then redumping - probably don't need to use this any more
        'dir=s',                   # dir where cache is found (defaults to $HOME/.vep/)
        'cache_region_size=i',     # size of region in bases for each cache file
        'no_slice_cache',          # tell API not to cache features on slice
        'standalone',              # standalone mode uses minimal set of modules installed in same dir, no DB connection
        'skip_db_check',           # don't compare DB parameters with cached
        'compress=s',              # by default we use zcat to decompress; user may want to specify gzcat or "gzip -dc"
        
        # debug
        'cluck',                   # these two need some mods to Bio::EnsEMBL::DBSQL::StatementHandle to work. Clucks callback trace and SQL
        'count_queries',           # counts SQL queries executed

#MY
        'vcf_out', #VCF output format, with annotations in INFO column
        #NGS-SNP options
        'flanking=i', #amount of flanking genomic sequence to write on each side of SNPs when -f option is used (Optional; default is 100 bases)
        'flanking_output=s', #write genomic flanking sequence for each SNP to this file (Optional)
        'comparison_species=s', #names of species to use when assessing sequence conservation, default all
        'compara_db=s', #the name of the compara database to be used, default 'Multi'
        'model=s', #the model species to use when filling in the 'Model_Annotations' column in the output,default value is 'homo_sapiens'
        'scoring_matrix=s', #the scoring matrix file to use for amino acid comparisons, default blosum62.mat
        'omes_script=s', #the location of the OMES_score.pl script
        'ncbi_script=s', #the location of the ncbi_search.pl script
        'uniprot_script=s', #the location of the ebi_fetch.pl script
        'run_omes', #perform the calculations used to fill in the 'Protein_OMES_Orthologues' and 'Protein_OMES_Family'values (SLOW), default off
        'no_ncbi', #do not perform annotation steps that involve retrieving data from NCBI
        'no_uniprot', #do not perform annotation steps that involve retrieving data from UniProt
        'flanking_for_context_conservation=i', #amount of flanking protein sequence on either side of SNP-affected residue to use for determining conservation of region containing coding SNP, default 10
        'degenerate_flanking', #indicate known SNP sites within flanking sequence as lowercase IUPAC DNA bases (Optional; default is to return unmodified reference sequence)
        'is_fake_vfa', #
        'omes_family', #whether or not omes should be calculated using protein family
        'omes_orthologues', #whether or not omes should be calculated using orthologues
        'omes_minimum_seqs=i', #the minimum number of sequences needed to calculate OMES score, after filtering by identity, default 10
        'omes_maximum_identity=i', #the maximum pairwise percent identity allowed between two sequences used in OMES calculation, default 90
        'omes_max_seqs=i', #the maximum number of sequences to be submitted to OMES calculation. Set to undef to use all, default 40
        #VW options
        'vw_host=s',
        'vw_port=s',
        'vw_user=s',
        'vw_password=s',
        'vw_platform=s',
        'vw_database=s',
#MY
    );
    
    # print usage message if requested or no args supplied
    if(defined($config->{help}) || !$args) {
        &usage;
        exit(0);
    }
    
    # config file?
    if(defined $config->{config}) {
        
        open CONFIG, $config->{config} or die "ERROR: Could not open config file \"".$config->{config}."\"\n";
        
        while(<CONFIG>) {
            next if /^\#/;
            my ($key, $value) = split /\s+|\=/;
            $key =~ s/^\-//g;
            $config->{$key} = $value unless defined $config->{$key};
        }
        
        close CONFIG;
    }

    # can't be both quiet and verbose
    die "ERROR: Can't be both quiet and verbose!" if defined($config->{quiet}) && defined($config->{verbose});
    
    # check file format
    if(defined $config->{format}) {
        die "ERROR: Unrecognised input format specified \"".$config->{format}."\"\n" unless $config->{format} =~ /pileup|vcf|guess/i;
    }
    
#MY
    $config->{vw_host} ||= 'myourshaw-dev.genome.ucla.edu';
    $config->{vw_port} ||= 3306;
    $config->{vw_user} ||= 'vw';
    $config->{vw_password} ||= 'vw';
    $config->{vw_platform} ||= 'mysql';
    $config->{vw_database} ||= 'vw';

    $config->{flanking} ||= 100;
    $config->{flanking_output} ||= undef;
    $config->{compara_db} ||= 'Multi';
    $config->{model} ||= 'homo_sapiens';
    $config->{scoring_matrix} ||= File::Spec->catfile(dirname(__FILE__), 'blosum62.mat');
    $config->{omes_script} ||= File::Spec->catfile(dirname(__FILE__), 'OMES_score.pl');
    $config->{ncbi_script} ||= File::Spec->catfile(dirname(__FILE__), 'ncbi_search.pl');
    $config->{uniprot_script} ||= File::Spec->catfile(dirname(__FILE__), 'ebi_fetch.pl');
    $config->{flanking_for_context_conservation} ||= 10;
    $config->{omes_family} ||= 1;
    $config->{omes_orthologues} ||= 1;
    $config->{omes_minimum_seqs} ||= 10;
    $config->{omes_maximum_identity} ||= 90;
    $config->{omes_max_seqs} ||= 40;
    $config->{degenerate_flanking} ||= undef;
    $config->{is_fake_vfa} ||= 0;
    #get scoring matrix object
    my $parser = Bio::Matrix::IO->new(
        -format => 'scoring',
        -file   => $config->{scoring_matrix}
    );
    $config->{matrix} = $parser->next_matrix;

    #determine maximum possible alignment score change value, for normalization
    #defaults to 15
    $config->{max_alignment_score_change}
        = get_max_alignment_score_change($config->{matrix});
    my %species_short = (
        Ailuropoda_melanoleuca => 'Am',
        Anolis_carolinensis => 'Ac',
        Bos_taurus => 'Bt',
        Caenorhabditis_elegans => 'Ce',
        Callithrix_jacchus => 'Cj',
        Canis_familiaris => 'Cf',
        Cavia_porcellus => 'Cp',
        Choloepus_hoffmanni => 'Ch',
        Ciona_intestinalis => 'Ci',
        Ciona_savignyi => 'Cs',
        Danio_rerio => 'Dr',
        Dasypus_novemcinctus => 'Dn',
        Dipodomys_ordii => 'Do',
        Drosophila_melanogaster => 'Dm',
        Echinops_telfairi => 'Et',
        Equus_caballus => 'Ec',
        Erinaceus_europaeus => 'Ee',
        Felis_catus => 'Fc',
        Gallus_gallus => 'Gga',
        Gasterosteus_aculeatus => 'Ga',
        Gorilla_gorilla => 'Gg',
        Homo_sapiens => 'Hs',
        Loxodonta_africana => 'La',
        Macaca_mulatta => 'Mm',
        Macropus_eugenii => 'Me',
        Meleagris_gallopavo => 'Mg',
        Microcebus_murinus => 'Mmur',
        Monodelphis_domestica => 'Md',
        Mus_musculus => 'Mmus',
        Myotis_lucifugus => 'Ml',
        Nomascus_leucogenys => 'Nl',
        Ochotona_princeps => 'Op',
        Ornithorhynchus_anatinus => 'Oa',
        Oryctolagus_cuniculus => 'Oc',
        Oryzias_latipes => 'Ol',
        Otolemur_garnettii => 'Og',
        Pan_troglodytes => 'Pt',
        Pongo_abelii => 'Pa',
        Procavia_capensis => 'Pc',
        Pteropus_vampyrus => 'Pv',
        Rattus_norvegicus => 'Rn',
        Saccharomyces_cerevisiae => 'Sc',
        Sorex_araneus => 'Sa',
        Spermophilus_tridecemlineatus => 'St',
        Sus_scrofa => 'Ss',
        Taeniopygia_guttata => 'Tg',
        Takifugu_rubripes => 'Tr',
        Tarsius_syrichta => 'Ts',
        Tetraodon_nigroviridis => 'Tn',
        Tupaia_belangeri => 'Tb',
        Tursiops_truncatus => 'Tt',
        Vicugna_pacos => 'Vp',
        Xenopus_tropicalis => 'Xt',    );
    $config->{species_short} = \%species_short;
    my @exclude_uniprot_cc = (
        'Copyright',
        'ALTERNATIVE PRODUCTS',
        'BIOPHYSICOCHEMICAL PROPERTIES',
        'BIOTECHNOLOGY',
        'MASS SPECTROMETRY',
        'PHARMACEUTICAL',
        'SEQUENCE CAUTION',
        'TOXIC DOSE'
        );
    $config->{exclude_uniprot_cc} = \@exclude_uniprot_cc;
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
    $config->{uniprot_site_features} = \@uniprot_site_features;

#MY
    # connection settings for Ensembl Genomes
    if($config->{genomes}) {
        $config->{host} ||= 'mysql.ebi.ac.uk';
        $config->{port} ||= 4157;
    }
    
    # connection settings for main Ensembl
    else {
        $config->{species} ||= "homo_sapiens";
        $config->{host}    ||= 'ensembldb.ensembl.org';
        $config->{port}    ||= 5306;
    }
    
    # output term
    if(defined $config->{terms}) {
        die "ERROR: Unrecognised consequence term type specified \"".$config->{terms}."\" - must be one of ensembl, so, ncbi\n" unless $config->{terms} =~ /ensembl|display|so|ncbi/i;
        if($config->{terms} =~ /ensembl|display/i) {
            $config->{terms} = 'display';
        }
        else {
            $config->{terms} = uc($config->{terms});
        }
    }
    
    # check nsSNP tools
    foreach my $tool(grep {defined $config->{lc($_)}} qw(SIFT PolyPhen Condel)) {
        die "ERROR: Unrecognised option for $tool \"", $config->{lc($tool)}, "\" - must be one of p (prediction), s (score) or b (both)\n" unless $config->{lc($tool)} =~ /^(s|p|b)/;
        
        die "ERROR: $tool not available for this species\n" unless $config->{species} =~ /human|homo/i;
        
        die "ERROR: $tool not available in standalone mode\n" if defined($config->{standalone});
    }
    
    # force quiet if outputting to STDOUT
    if(defined($config->{output_file}) && $config->{output_file} =~ /stdout/i) {
        delete $config->{verbose} if defined($config->{verbose});
        $config->{quiet} = 1;
    }
    
    # summarise options if verbose
    if(defined $config->{verbose}) {
        my $header =<<INTRO;
#----------------------------------#
# ENSEMBL VARIANT EFFECT PREDICTOR #
#----------------------------------#

version $VERSION

By Will McLaren (wm2\@ebi.ac.uk)

Configuration options:

INTRO
        print $header;
        
        my $max_length = (sort {$a <=> $b} map {length($_)} keys %$config)[-1];
        
        foreach my $key(sort keys %$config) {
            print $key.(' ' x (($max_length - length($key)) + 4)).$config->{$key}."\n";
        }
        
        print "\n".("-" x 20)."\n\n";
    }
    
    # set defaults
    $config->{user}              ||= 'anonymous';
    $config->{buffer_size}       ||= 5000;
    $config->{chunk_size}        ||= '50kb';
    $config->{output_file}       ||= "variant_effect_output.txt";
    $config->{tmpdir}            ||= '/tmp';
    $config->{format}            ||= 'guess';
    $config->{terms}             ||= 'display';
    $config->{gene}              ||= 1 unless defined($config->{whole_genome});
    $config->{cache_region_size} ||= 1000000;
    $config->{dir}          ||= join '/', ($ENV{'HOME'}, '.vep');
    $config->{compress}          ||= 'zcat';
    
    # warn users still using whole_genome flag
    if(defined($config->{whole_genome})) {
        debug("INFO: Whole-genome mode is now the default run-mode for the script. To disable it, use --no_whole_genome") unless defined($config->{quiet});
    }
    
    $config->{whole_genome}      = 1 unless defined $config->{no_whole_genome};
    $config->{include_failed}    = 1 unless defined $config->{include_failed};
    $config->{chunk_size}        =~ s/mb?/000000/i;
    $config->{chunk_size}        =~ s/kb?/000/i;
    $config->{cache_region_size} =~ s/mb?/000000/i;
    $config->{cache_region_size} =~ s/kb?/000/i;
    
    # cluck and display executed SQL?
    $Bio::EnsEMBL::DBSQL::StatementHandle::cluck = 1 if defined($config->{cluck});
    
    # standalone needs cache, can't use HGVS
    if(defined($config->{standalone})) {
        $config->{cache} = 1;
        
        die("ERROR: Cannot generate HGVS coordinates in standalone mode") if defined($config->{hgvs});
        
        die("ERROR: Cannot analyse regulatory features in standalone mode") if defined($config->{regulatory});
    }
    
    # no_slice_cache, prefetch and whole_genome have to be on to use cache or upfront
    if(defined($config->{cache}) || defined($config->{upfront})) {
        $config->{prefetch} = 1;
        $config->{no_slice_cache} = 1;
        $config->{whole_genome} = 1;
        $config->{strip} = 1;
        
        # scan should also be on for upfront
        $config->{scan} = 1 if defined($config->{upfront});
    }
    
    $config->{build} = $config->{rebuild} if defined($config->{rebuild});
    
    # force options for full build
    if(defined($config->{build})) {
        $config->{prefetch} = 1;
        $config->{gene} = 1;
        $config->{hgnc} = 1;
        $config->{no_slice_cache} = 1;
        $config->{cache} = 1;
        $config->{strip} = 1;
        $config->{write_cache} = 1;
    }
    
    # connect to databases
    $config->{reg} = &connect_to_dbs($config);
    $config->{vw_conn} = &connect_to_vw($config);
    # complete dir with species name and db_version
    $config->{dir} .= '/'.(
        join '/', (
            defined($config->{standalone}) ? $config->{species} : ($config->{reg}->get_alias($config->{species}) || $config->{species}),
            $config->{db_version} || $config->{reg}->software_version
        )
    );
    
    # warn user cache directory doesn't exist
    if(!-e $config->{dir}) {
        
        # if using write_cache
        if(defined($config->{write_cache})) {
            debug("INFO: Cache directory ", $config->{dir}, " not found - it will be created") unless defined($config->{quiet});
        }
        
        # want to read cache, not found
        elsif(defined($config->{cache})) {
            die("ERROR: Cache directory ", $config->{dir}, " not found");
        }
    }
    
    # suppress warnings that the FeatureAdpators spit if using no_slice_cache
    Bio::EnsEMBL::Utils::Exception::verbose(1999) if defined($config->{no_slice_cache});
    
    # get adaptors
    if(defined($config->{cache})) {
        
        # try and load adaptors from cache
        if(!&load_dumped_adaptor_cache($config)) {
            &get_adaptors($config);
            &dump_adaptor_cache($config) if defined($config->{write_cache});
        }
        
        # check cached adaptors match DB params
        else {
            my $dbc = $config->{sa}->{dbc};
        
            my $ok = 1;
            
            if($dbc->{_host} ne $config->{host}) {
                
                # ens-livemirror, useastdb and ensembldb should all have identical DBs
                unless(
                    (
                        $dbc->{_host} eq 'ens-livemirror'
                        || $dbc->{_host} eq 'ensembldb.ensembl.org'
                        || $dbc->{_host} eq 'useastdb.ensembl.org'
                    ) && (
                        $config->{host} eq 'ens-livemirror'
                        || $config->{host} eq 'ensembldb.ensembl.org'
                        || $config->{host} eq 'useastdb.ensembl.org'
                    )
                ) {
                    $ok = 0;
                }
                
                # but we still need to reconnect
                debug("INFO: Defined host ", $config->{host}, " is different from cached ", $dbc->{_host}, " - reconnecting to host") unless defined($config->{quiet});
                
                &get_adaptors($config);
            }
            
            if(!$ok) {
                if(defined($config->{skip_db_check})) {
                    debug("INFO: Defined host ", $config->{host}, " is different from cached ", $dbc->{_host}) unless defined($config->{quiet});
                }
                else {
                    die "ERROR: Defined host ", $config->{host}, " is different from cached ", $dbc->{_host}, ". If you are sure this is OK, rerun with -skip_db_check flag set";
                }
            }
        }
    }
    else {
        &get_adaptors($config);
        &dump_adaptor_cache($config) if defined($config->{write_cache})
    }
    
    # get terminal width for progress bars
    unless(defined($config->{quiet})) {
        my $width;
        
        # module may not be installed
        eval {
            use Term::ReadKey;
        };
        
        if(!$@) {
            my ($w, $h);
            
            # module may be installed, but e.g.
            eval {
                ($w, $h) = GetTerminalSize();
            };
            
            $width = $w if defined $w;
        }
        
        $width ||= 60;
        $width -= 12;
        $config->{terminal_width} = $width;
    }
    
    # jump out to build cache if requested
    if(defined($config->{build})) {
        
        # build the cache
        debug("Building cache for ".$config->{species}) unless defined($config->{quiet});
        &build_full_cache($config, $config->{rebuild});
        
        # exit script
        debug("Finished building cache") unless defined($config->{quiet});
        exit(0);
    }
    
    # warn user DB will be used for SIFT/PolyPhen/Condel
    if(defined($config->{cache}) && (defined($config->{sift}) || defined($config->{polyphen}) || defined($config->{condel}) || defined($config->{hgvs}) || defined($config->{regulatory}))) {
        debug("INFO: Database will be accessed for SIFT/PolyPhen/Condel, HGVS and regulatory features") unless defined($config->{quiet});
    }
    
    # get list of chrs if supplied
    if(defined($config->{chr})) {
        my %chrs;
        
        foreach my $val(split /\,/, $config->{chr}) {
            my @nnn = split /\-/, $val;
            
            foreach my $chr($nnn[0]..$nnn[-1]) {
                $chrs{$chr} = 1;
            }
        }
        
        $config->{chr} = \%chrs;
    }
    
    # get input file handle
    $config->{in_file_handle} = &get_in_file_handle($config);
    
    # configure output file
#MY
    if (! $config->{vcf_out}){
        $config->{out_file_handle} = &get_out_file_handle($config);
    }
#MY
    
    return $config;
}

# connects to DBs; in standalone mode this just loads registry module
sub connect_to_dbs {
    my $config = shift;
    
    # get registry
    my $reg = 'Bio::EnsEMBL::Registry';
    
    unless(defined($config->{standalone})) {
        # load DB options from registry file if given
        if(defined($config->{registry})) {
            debug("Loading DB config from registry file ", $config->{registry}) unless defined($config->{quiet});
            $reg->load_all(
                $config->{registry},
                $config->{verbose},
                undef,
                $config->{no_slice_cache}
            );
        }
        
        # otherwise manually connect to DB server
        else {
            $reg->load_registry_from_db(
                -host       => $config->{host},
                -user       => $config->{user},
                -pass       => $config->{password},
                -port       => $config->{port},
                -db_version => $config->{db_version},
                -species    => $config->{species} =~ /^[a-z]+\_[a-z]+/i ? $config->{species} : undef,
                -verbose    => $config->{verbose},
                -no_cache   => $config->{no_slice_cache},
            );
        }
        
        eval { $reg->set_reconnect_when_lost() };
        
        if(defined($config->{verbose})) {
            # get a meta container adaptors to check version
            my $core_mca = $reg->get_adaptor($config->{species}, 'core', 'metacontainer');
            my $var_mca = $reg->get_adaptor($config->{species}, 'variation', 'metacontainer');
            
            if($core_mca && $var_mca) {
                debug(
                    "Connected to core version ", $core_mca->get_schema_version, " database ",
                    "and variation version ", $var_mca->get_schema_version, " database"
                );
            }
        }
    }
    
    return $reg;
}

# get adaptors from DB
sub get_adaptors {
    my $config = shift;
    
    die "ERROR: No registry" unless defined $config->{reg};
    
    $config->{vfa} = $config->{reg}->get_adaptor($config->{species}, 'variation', 'variationfeature');
    $config->{tva} = $config->{reg}->get_adaptor($config->{species}, 'variation', 'transcriptvariation');
    
    # get fake ones for species with no var DB
    if(!defined($config->{vfa})) {
        $config->{vfa} = Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor->new_fake($config->{species});
    }
    else {
        $config->{vfa}->db->include_failed_variations($config->{include_failed}) if defined($config->{vfa}->db) && $config->{vfa}->db->can('include_failed_variations');
    }
    
    $config->{sa}  = $config->{reg}->get_adaptor($config->{species}, 'core', 'slice');
    $config->{ga}  = $config->{reg}->get_adaptor($config->{species}, 'core', 'gene');
    $config->{ta}  = $config->{reg}->get_adaptor($config->{species}, 'core', 'transcript');
    $config->{mca} = $config->{reg}->get_adaptor($config->{species}, 'core', 'metacontainer');
    $config->{csa} = $config->{reg}->get_adaptor($config->{species}, 'core', 'coordsystem');
#MY
    #adaptors used by NGS-SNP
    $config->{ma} = $config->{reg}->get_adaptor('Multi', 'compara', 'Member');
    $config->{ha} = $config->{reg}->get_adaptor('Multi', 'compara', 'Homology');
    $config->{fa} = $config->{reg}->get_adaptor('Multi', 'compara', 'Family');
    $config->{goa} = $config->{reg}->get_adaptor('Multi', 'Ontology', 'GOTerm');
    $config->{translation_adaptor} = $config->{reg}->get_adaptor($config->{species}, 'core', 'translation');
    #these adaptors are used for the Model_Annotations field.
    #if they cannot be created then this field is not filled in.
    if ( defined( $config->{model} ) ) {
        $config->{model_translation_adaptor} = $config->{reg}->get_adaptor($config->{model}, 'core', 'translation');
        $config->{model_transcript_adaptor} = $config->{reg}->get_adaptor($config->{model}, 'core', 'transcript');
        $config->{model_gene_adaptor} = $config->{reg}->get_adaptor($config->{model}, 'core', 'gene');
        $config->{model_variation_adaptor} = $config->{reg}->get_adaptor($config->{model}, 'variation', 'variation');
        $config->{model_variationfeature_adaptor} = $config->{reg}->get_adaptor($config->{model}, 'variation', 'variationfeature');
        $config->{model_slice_adaptor} = $config->{reg}->get_adaptor($config->{model}, 'core', 'slice');
        $config->{model_variationannotation_adaptor} = $config->{reg}->get_adaptor($config->{model}, 'variation', 'variationannotation');

        #if (   ( !defined( $config->{model_translation_adaptor} ) )
        #    || ( !defined( $config->{model_transcript_adaptor} ) )
        #    || ( !defined( $config->{model_gene_adaptor} ) )
        #    || ( !defined( $config->{model_variation_adaptor} ) )
        #    || ( !defined( $config->{model_variationfeature_adaptor} ) )
        #    || ( !defined( $config->{model_slice_adaptor} ) )
        #    || ( !defined( $config->{model_variationannotation_adaptor} ) ) )
        #{
        #    message( $config->{verbose}, $config->{log_file},
        #        "Unable to get adaptors for the species '$config->{model}' specified using the '-model' option.\n"
        #    );
        #}
    }
#MY
    # cache schema version
    $config->{mca}->get_schema_version if defined $config->{mca};
    
    # check we got slice adaptor - can't continue without a core DB
    die("ERROR: Could not connect to core database\n") unless defined $config->{sa};
}

# gets file handle for input
sub get_in_file_handle {
    my $config = shift;

    # define the filehandle to read input from
    my $in_file_handle = new FileHandle;
    
    if(defined($config->{input_file})) {
        
        # check defined input file exists
        die("ERROR: Could not find input file ", $config->{input_file}, "\n") unless -e $config->{input_file};
        
        if($config->{input_file} =~ /\.gz$/){
            $in_file_handle->open($config->{compress}." ". $config->{input_file} . " | " ) or die("ERROR: Could not read from input file ", $config->{input_file}, "\n");
        }
        else {
            $in_file_handle->open( $config->{input_file} ) or die("ERROR: Could not read from input file ", $config->{in_file}, "\n");
        }
    }
    
    # no file specified - try to read data off command line
    else {
        $in_file_handle = 'STDIN';
        debug("Reading input from STDIN (or maybe you forgot to specify an input file?)...") unless defined $config->{quiet};
    }
    
    return $in_file_handle;
}

# gets file handle for output and adds header
sub get_out_file_handle {
    my $config = shift;
    
    # define filehandle to write to
    my $out_file_handle = new FileHandle;
    
    # check if file exists
    if(-e $config->{output_file} && !defined($config->{force_overwrite})) {
        die("ERROR: Output file ", $config->{output_file}, " already exists. Specify a different output file with --output_file or overwrite existing file with --force_overwrite\n");
    }
    
    if($config->{output_file} =~ /stdout/i) {
        $out_file_handle = *STDOUT;
    }
    else {
        $out_file_handle->open(">".$config->{output_file}) or die("ERROR: Could not write to output file ", $config->{output_file}, "\n");
    }
    
    # make header
    my $time = &get_time;
    my $db_string = $config->{mca}->dbc->dbname." on ".$config->{mca}->dbc->host if defined $config->{mca};
    $db_string .= "\n## Using cache in ".$config->{dir} if defined($config->{cache});
    my $version_string =
        "Using API version ".$config->{reg}->software_version.
        ", DB version ".(defined $config->{mca} && $config->{mca}->get_schema_version ? $config->{mca}->get_schema_version : '?');

my $input_file_modified_time = localtime((stat $config->{in_file_handle})[9]);

    my $header =<<HEAD;
##ENSEMBL VARIANT EFFECT PREDICTOR v$VERSION
##NGS-SNP
##FILE: annotate_SNPs.pl
##AUTH: Paul Stothard (stothard\@ualberta.ca)
##DATE: May 16, 2011
##Input file : $config->{input_file}
##Input file time : $input_file_modified_time
##Output file : $config->{output_file}
##Output produced at : $time
##Connected to : $db_string
##$version_string
##VEP column keys:
##Uploaded_variation : VCF ID or as chromosome_start_alleles
##Location : in standard coordinate format (chr:start or chr:start-end)
##Allele : the variant allele used to calculate the consequence
##Gene : Ensembl stable ID of affected gene
##Feature : Ensembl stable ID of feature
##Feature_type : type of feature. Currently one of Transcript, RegulatoryFeature, MotifFeature.
##Consequence : the most severe SO consequence type of this variation
##Consequence_rank : Rank order of consequence (<=9 may be most interesting)
##Consequences_all : All the consequences of the variation
##cDNA_position : Relative position in cDNA - base pair position in cDNA sequence
##CDS_position : Relative position in CDS - base pair position in coding sequence
##Protein_position : Relative position in protein - amino acid position in protein
##Amino_acids : Amino acid change - only given if the variation affects the protein-coding sequence
##Codons : the alternate codons with the variant base highlighted as bold (HTML) or upper case (text)
##Existing_variation : identifier of existing variation
##HGNC : HGNC gene identifier
##ENSP : Ensembl protein identifer
##HGVSc : HGVS coding sequence name
##HGVSp : HGVS protein sequence name
##SIFT_prediction : SIFT prediction
##SIFT_score : SIFT score
##PolyPhen_prediction : PolyPhen prediction
##PolyPhen_score : PolyPhen score
##Condel_prediction : Condel SIFT/PolyPhen consensus prediction
##Condel_score : Condel SIFT/PolyPhen consensus score
##MATRIX : The source and identifier of a transcription factor binding profile aligned at this position
##HIGH_INF_POS : A flag indicating if the variant falls in a high information position of a transcription factor binding profile
##Alignment_Score_Change : the alignment score for the variant amino acid vs. the orthologous amino acids minus the alignment score for the reference amino acid vs. the orthologous amino acids. When there are multiple variant amino acids, the most extreme difference is given. A positive value indicates that the variant amino acid better resembles the orthologues than does the reference amino acid, whereas a negative value indicates that the reference amino acid better resembles the orthologues than does the variant amino acid. The value is scaled to between -1 and 1
##C_blosum : a measure of the conservation of the reference amino acid with the aligned amino acids in orthologous sequences. The alignment score for the reference amino acid vs. the orthologous amino acids is divided by the alignment score that would be obtained if all the orthologous residues matched the reference. Higher values tend to be associated with changes to the amino acid having a greater functional consequence. The formula used is equivalent to the C_blosum formula given in Kowarsch A et al. (2010 PLoS Comput Biol 6(9): e1000923), except that in NGS-SNP scoring matrices other than BLOSUM62 can be used
##Context_Conservation : the average percent identity obtained when the region of the reference protein containing the SNP-affected residue is aligned with the orthologous region from other species. The size of the region examined is determined by the -cf option, which specifies how much sequence on either side of the SNP to examine. For example, if '-cf 10' is used, the size of the region is 10 + 1 + 10 = 21
##Amino_Acids_In_Orthologues : the amino acids aligned with the reference amino acid in orthologous sequences
##Orthologue_Species : the species from which sequences were obtained to generate the 'Amino_Acid_In_Orthologues', 'Alignment_Score_Change', 'C_blosum', and 'Context_Conservation' values. The order of the species matches the order used to generate the 'Amino_Acid_In_Orthologues' value. Appreviations - Ailuropoda_melanoleuca=Am,Anolis_carolinensis=Ac,Bos_taurus=Bt,Caenorhabditis_elegans=Ce,Callithrix_jacchus=Cj,Canis_familiaris=Cf,Cavia_porcellus=Cp,Choloepus_hoffmanni=Ch,Ciona_intestinalis=Ci,Ciona_savignyi=Cs,Danio_rerio=Dr,Dasypus_novemcinctus=Dn,Dipodomys_ordii=Do,Drosophila_melanogaster=Dm,Echinops_telfairi=Et,Equus_caballus=Ec,Erinaceus_europaeus=Ee,Felis_catus=Fc,Gallus_gallus=Gga,Gasterosteus_aculeatus=Ga,Gorilla_gorilla=Gg,Homo_sapiens=Hs,Loxodonta_africana=La,Macaca_mulatta=Mm,Macropus_eugenii=Me,Meleagris_gallopavo=Mg,Microcebus_murinus=Mmur,Monodelphis_domestica=Md,Mus_musculus=Mmus,Myotis_lucifugus=Ml,Nomascus_leucogenys=Nl,Ochotona_princeps=Op,Ornithorhynchus_anatinus=Oa,Oryctolagus_cuniculus=Oc,Oryzias_latipes=Ol,Otolemur_garnettii=Og,Pan_troglodytes=Pt,Pongo_abelii=Pa,Procavia_capensis=Pc,Pteropus_vampyrus=Pv,Rattus_norvegicus=Rn,Saccharomyces_cerevisiae=Sc,Sorex_araneus=Sa,Spermophilus_tridecemlineatus=St,Sus_scrofa=Ss,Taeniopygia_guttata=Tg,Takifugu_rubripes=Tr,Tarsius_syrichta=Ts,Tetraodon_nigroviridis=Tn,Tupaia_belangeri=Tb,Tursiops_truncatus=Tt,Vicugna_pacos=Vp,Xenopus_tropicalis=Xt
##Reference_Splice_Site : the sequence of the splice site that is altered by the SNP. The splice site bases (i.e. the first two and last two bases in the intron) are given as they appear on the sense strand of the reference sequence. This value is reported when the functional class is 'ESSENTIAL_SPLICE_SITE'
##Variant_Splice_Site : the sequence of the splice site that is altered by the SNP. The splice site bases (i.e. the first two and last two bases in the intron) are given as they appear on the sense strand of the variant sequence. This value is reported when the functional class is 'ESSENTIAL_SPLICE_SITE'
##Protein_Length : Number of amino acids in translation of the transcript of the reference sequence
##Protein_Length_Decrease : gives the length in amino acids of the protein segment that is lost because of an allele that introduces a stop codon (i.e. functional class is 'STOP_GAINED'). The value given in parentheses is the length of the lost protein segment expressed as a percentage of the length of the reference protein
##Protein_Sequence_Lost : the peptide sequence removed from the reference sequence because of an allele that introduces a stop codon (i.e. functional class is 'STOP_GAINED')
##Protein_Length_Increase : gives the length in amino acids of the protein segment that is gained because of an allele that removes a stop codon (i.e. functional class is 'STOP_LOST'). The value given in parentheses is the length of the gained protein segment expressed as a percentage of the length of the reference protein
##Protein_Sequence_Gained : the peptide sequence added to the reference sequence because of an allele that removes a stop codon (i.e. functional class is 'STOP_LOST')
##Gene_Description : a short description of the relevant gene
##Entrez_Gene_Name : the Entrez Gene name of the relevant gene
##Entrez_Gene_ID : the Entrez Gene ID of the relevant gene
##UniProt_ID : the UniProt ID of the relevant protein
##Overlapping_Protein_Domains : protein domains that overlap with the position of the affected amino acid
##Gene_Ontology : GO slim IDs and terms associated with the relevant transcript
##Phenotypes_Position : phenotypic information associated with known variation in the model species at the SNP position. If the input SNPs are from the model species then their locations are used to identify known variations at the same locations, and any phenotypic information linked to these variations is reported. If the input SNPs are not from the model species and the input SNP alters a protein, then protein alignment is used to find the orthologous genomic region from the model. Model species variations that alter this region are obtained, and any phenotypic information linked to these variations is reported (Ensembl)
##Phenotypes_Gene : phenotypic information associated with the model species gene (NCBI)
##HGMD : HGMD ID
##OMIM_Disorder : (OMIM)
##KEGG_Pathway : (KEGG)
##UniProt_ID : The first item on the ID line is the entry name of the sequence. This name is a useful means of identifying a sequence, but it is not a stable identifier as is the accession number (UniProt)
##UniProtKB_AC : The 'primary accession number' (UniProt)
##VARIANT : Description of a natural variant of the protein (UniProt)
##MUTAGEN : Site which has been experimentally altered by mutagenesis. (UniProt)
##SITES : ACT_SITE - Amino acid(s) involved in the activity of an enzyme; BINDING - Binding site for any chemical group (co-enzyme, prosthetic group, etc.); CA_BIND - Extent of a calcium-binding region; DISULFID - Disulfide bond; DNA_BIND - Extent of a DNA-binding region; METAL - Binding site for a metal ion; NP_BIND - Extent of a nucleotide phosphate-binding region; SITE - Any interesting single amino-acid site on the sequence, that is not defined by another feature key. It can also apply to an amino acid bond which is represented by the positions of the two flanking amino acids; ZN_FING - Extent of a zinc finger region (UniProt)
##OTHER_OVERLAPPING_FEATURES : Other protein features not in VARIANT, MUTAGEN, SITES, which overlap with the position of the relevant amino acid. (UniProt)
##ACs : ACcession number(s) associated with an entry (UniProt)
##ALLERGEN : Information relevant to allergenic proteins (UniProt)
##ALTERNATIVE_PRODUCTS : Description of the existence of related protein sequence(s) produced by alternative splicing of the same gene, alternative promoter usage, ribosomal frameshifting or by the use of alternative initiation codons
##BIOPHYSICOCHEMICAL_PROPERTIES : Description of the information relevant to biophysical and physicochemical data and information on pH dependence, temperature dependence, kinetic parameters, redox potentials, and maximal absorption (UniProt)
##BIOTECHNOLOGY : Description of the use of a specific protein in a biotechnological process (UniProt)
##CATALYTIC_ACTIVITY : Description of the reaction(s) catalyzed by an enzyme (UniProt)
##CAUTION : Warning about possible errors and/or grounds for confusion (UniProt)
##COFACTOR : Description of non-protein substance required by an enzyme to be active (UniProt)
##DE : General descriptive information about the sequence stored. This information is generally sufficient to identify the protein precisely. (UniProt)
##DEVELOPMENTAL_STAGE : Description of the developmentally-specific expression of mRNA or protein (UniProt)
##DISEASE : Description of the disease(s) associated with a deficiency of a protein (UniProt)
##DISRUPTION_PHENOTYPE : Description of the effects caused by the disruption of the gene coding for the protein. Note that we only describe effects caused the complete absence of a gene and thus a protein in vivo (null mutants caused by random or target deletions, insertions of a transposable element etc.) To avoid description of phenotypes due to partial or dominant negative mutants, missense mutations are not described in this comment, but in FT MUTAGEN instead. (UniProt)
##DOMAIN : Description of the domain(s) present in a protein (UniProt)
##DRs : Database cross-Reference pointers to information in external data resources (UniProt)
##ENSG : Database of automatically annotated sequences of large genomes gene identifier (Ensembl database) (UniProt)
##ENSP : Database of automatically annotated sequences of large genomes protein identifier (Ensembl database) (UniProt)
##ENST : Database of automatically annotated sequences of large genomes transcript identifier (Ensembl database) (UniProt)
##ENZYME_REGULATION : Description of an enzyme regulatory mechanism (UniProt)
##FUNCTION : General description of the function(s) of a protein (UniProt)
##Gene : The official gene name (UniProt)
##GeneNames : (a.k.a gene symbols). The name(s) used to represent a gene (UniProt)
##GO : Gene Ontology (GO) database accession number/primary key (UniProt)
##GO_term : Gene Ontology (GO) database; this field is a 1-letter abbreviation for one of the 3 ontology aspects, separated from the GO term by a column. If the term is longer than 46 characters, the first 43 characters are indicated followed by 3 dots ('...'). The abbreviations for the 3 distinct aspects of the ontology are P (biological Process), F (molecular Function), and C (cellular Component) (UniProt)
##HGNC : Human gene nomenclature database (HGNC); the gene designation. If the gene designation is not available, a dash ('-') is used (UniProt)
##INDUCTION : Description of the effects of environmental factors on the gene expression (UniProt)
##INTERACTION : Interaction with other protein(s) (UniProt)
##KEGG : Kyoto encyclopedia of genes and genomes database accession number/primary key (UniProt)
##KEYWORDS : List of controlled vocabulary which summarises the content of an entry (UniProt)
##MASS_SPECTROMETRY : Reports the exact molecular weight of a protein or part of a protein as determined by mass spectrometric methods (UniProt)
##MIM_gene : Mendelian Inheritance in Man Database (MIM) gene database accession number/primary key (UniProt)
##MIM_phenotype : Mendelian Inheritance in Man Database (MIM) phenotype database accession number/primary key(UniProt)
##MISCELLANEOUS : Any relevant information that doesn't fit in any other defined sections (UniProt)
##PATHWAY : Description of associated metabolic pathways (UniProt)
##Pathway_Interaction : NCI-Nature Pathway Interaction Database 'full pathway name' (UniProt)
##PE : The evidences of the existence of a protein (1-Evidence at protein level, 2-Evidence at transcript level, 3-Inferred from homology, 4-Predicted, 5-Uncertain) (UniProt)
##PHARMACEUTICAL : Description of the use of a protein as a pharmaceutical drug (UniProt)
##POLYMORPHISM : Description of polymorphism(s) (UniProt)
##PTM : Description of post-translational modifications (UniProt)
##Reactome : Curated resource of core pathways and reactions in human biology (Reactome) name of the pathway (UniProt)
##RecName : The name recommended by the UniProt consortium (UniProt)
##RefSeq_NM : NCBI reference sequences nucleotide sequence identifier (UniProt)
##RefSeq_NP : NCBI reference sequences database accession number/primary key (UniProt)
##RNA_EDITING : Description of amino acid change(s) due to RNA editing (UniProt)
##SEQUENCE_CAUTION : Description of protein sequence reports that differ from the sequence that is shown in UniProtKB due to conflicts that are not described in FT CONFLICT lines, such as frameshifts, erroneous gene model predictions, etc. (UniProt)
##SIMILARITY : Description of the sequence similaritie(s) with other proteins (UniProt)
##SQ : Amino acid sequence (UniProt)
##SUBCELLULAR_LOCATION : Description of the subcellular location of the mature protein (UniProt)
##SUBUNIT : Description of the quaternary structure of a protein (UniProt)
##TISSUE_SPECIFICITY : Description of the tissue-specific expression of mRNA or protein (UniProt)
##UCSC : UCSC genome browser database accession number/primary key (UniProt)
##WEB_RESOURCE : Links to related web resource(s) or database(s) (UniProt)##Overlapping_Protein_Features_Orthologues : protein features from the model species gene that overlap with the position of the relevant amino acid (UniProt)
##Protein_OMES_Orthologues : reports the residues in the reference sequence that are most mutationally correlated with the altered residue, as determined using the Observed Minus Expected Squared (OMES) method described in Fodor and Aldrich (2004 Proteins 56(2): 211-221). Orthologous protein sequences are obtained from Ensembl and aligned with the reference using Muscle. All column pairs in the resulting multiple alignment are then compared to identify correlated residues. The output for each correlated residue is in the form 'score(X)(Y)' where 'score' is the raw OMES score, 'X' is the position of the SNP-altered residue in the reference sequence, and 'Y' is the position of the correlated residue in the reference sequence. 'Y' may sometimes be given as '-', which indicates that the correlated position detected in the multiple alignment is absent from the reference sequence. The value for this key is only calculated if the -run_omes option is specified
##Protein_OMES_Family : determined in the same manner as the 'Protein_OMES_Orthologues' value, except that non-orthologous sequences belonging to the same protein family as the reference sequence are also used in the calculation. The value for this key is only calculated if the -run_omes option is specified
##Canonical_Transcript : 1 if the transcript is the Ensembl canonical transcript - the longest CDS, if the gene has translated transcripts, or the longest cDNA
HEAD

    # add headers
    print $out_file_handle $header;
    
    # add column headers
    print $out_file_handle '#', (join "\t", @OUTPUT_COLS);
    print $out_file_handle "\n";
    
    return $out_file_handle;
}

# parses a line of input
sub parse_line {
    my $config = shift;
    my $line   = shift;
    
    my @data = (split /\s+/, $_);
    
    # pileup: chr1 60 T A
    if(
       ($config->{format} =~ /pileup/i) ||
       (
            $data[0] =~ /(chr)?\w+/ &&
            $data[1] =~ /\d+/ &&
            $data[2] =~ /^[ACGTN-]+$/ &&
            $data[3] =~ /^[ACGTNRYSWKM*+\/-]+$/
        )
    ) {
        my @return = ();
        
        if($data[2] ne "*"){
            my $var;
            
            if($data[3] =~ /^[A|C|G|T]$/) {
                $var = $data[3];
            }
            else {
                ($var = unambiguity_code($data[3])) =~ s/$data[2]//ig;
            }
            if(length($var)==1){
                push @return, [$data[0], $data[1], $data[1], $data[2]."/".$var, 1, undef];
            }
            else{
                for my $nt(split //,$var){
                    push @return, [$data[0], $data[1], $data[1], $data[2]."/".$nt, 1, undef];
                }
            }
        }
        else{ #indel
            my @genotype=split /\//,$data[3];
            foreach my $allele(@genotype){
                if(substr($allele,0,1) eq "+") { #ins
                    push @return, [$data[0], $data[1]+1, $data[1], "-/".substr($allele,1), 1, undef];
                }
                elsif(substr($allele,0,1) eq "-"){ #del
                    push @return, [$data[0], $data[1], $data[1]+length($data[3])-4, substr($allele,1)."/-", 1, undef];
                }
                elsif($allele ne "*"){
                    warn("WARNING: invalid pileup indel genotype: $line\n") unless defined $config->{quiet};
                    push @return, ['non-variant'];
                }
            }
        }
        return \@return;
    }
    
    # VCF: 20      14370   rs6054257 G     A      29    0       NS=58;DP=258;AF=0.786;DB;H2          GT:GQ:DP:HQ
    elsif(
        ($config->{format} =~ /vcf/i) ||
        (
            $data[0] =~ /(chr)?\w+/ &&
            $data[1] =~ /\d+/ &&
            $data[3] =~ /^[ACGTN-]+$/ &&
            $data[4] =~ /^([\.ACGTN-]+\,?)+$/
        )
    ) {
        
        # non-variant line in VCF, return dummy line
        if($data[4] eq '.') {
            return [['non-variant']];
        }
        
        # get relevant data
        my ($chr, $start, $end, $ref, $alt) = ($data[0], $data[1], $data[1], $data[3], $data[4]);
        
        if(defined $config->{gp}) {
            $chr = undef;
            $start = undef;
            
            foreach my $pair(split /\;/, $data[7]) {
                my ($key, $value) = split /\=/, $pair;
                if($key eq 'GP') {
                    ($chr,$start) = split /\:/, $value;
                    $end = $start;
                }
            }
            
            unless(defined($chr) and defined($start)) {
                warn "No GP flag found in INFO column" unless defined $config->{quiet};
                return [['non-variant']];
            }
        }
        
        # adjust end coord
        $end += (length($ref) - 1);
        
        # find out if any of the alt alleles make this an insertion or a deletion
        my ($is_indel, $is_sub, $ins_count, $total_count);
        foreach my $alt_allele(split /\,/, $alt) {
            $is_indel = 1 if $alt_allele =~ /D|I/;
            $is_indel = 1 if length($alt_allele) != length($ref);
            $is_sub = 1 if length($alt_allele) == length($ref);
            $ins_count++ if length($alt_allele) > length($ref);
            $total_count++;
        }
        
        # multiple alt alleles?
        if($alt =~ /\,/) {
            if($is_indel) {
                
                my @alts;
                
                if($alt =~ /D|I/) {
                    foreach my $alt_allele(split /\,/, $alt) {
                        # deletion (VCF <4)
                        if($alt_allele =~ /D/) {
                            push @alts, '-';
                        }
                        
                        elsif($alt_allele =~ /I/) {
                            $alt_allele =~ s/^I//g;
                            push @alts, $alt_allele;
                        }
                    }
                }
                
                else {
                    $ref = substr($ref, 1);
                    $ref = '-' if $ref eq '';
                    $start++;
                    
                    foreach my $alt_allele(split /\,/, $alt) {
                        $alt_allele = substr($alt_allele, 1);
                        $alt_allele = '-' if $alt_allele eq '';
                        push @alts, $alt_allele;
                    }
                }
                
                $alt = join "/", @alts;
            }
            
            else {
                # for substitutions we just need to replace ',' with '/' in $alt
                $alt =~ s/\,/\//;
            }
        }
        
        else {
            if($is_indel) {
                # deletion (VCF <4)
                if($alt =~ /D/) {
                    my $num_deleted = $alt;
                    $num_deleted =~ s/\D+//g;
                    $end += $num_deleted - 1;
                    $alt = "-";
                    $ref .= ("N" x ($num_deleted - 1)) unless length($ref) > 1;
                }
                
                # insertion (VCF <4)
                elsif($alt =~ /I/) {
                    $ref = '-';
                    $alt =~ s/^I//g;
                    $start++;
                }
                
                # insertion or deletion (VCF 4+)
                else {
                    # chop off first base
                    $ref = substr($ref, 1);
                    $alt = substr($alt, 1);
                    
                    $start++;
                    
                    if($ref eq '') {
                        # make ref '-' if no ref allele left
                        $ref = '-';
                    }
                    
                    # make alt '-' if no alt allele left
                    $alt = '-' if $alt eq '';
                }
            }
        }
        
        return [[$chr, $start, $end, $ref."/".$alt, 1, ($data[2] eq '.' ? undef : $data[2])]];
        
    }
    
    # our format
    else {
        # we allow commas as delimiter so re-split
        @data = (split /\s+|\,/, $_);
        return [\@data];
    }
}

# takes a hash of VFs and fetches consequences by pre-fetching overlapping transcripts
# from database and/or cache
sub whole_genome_fetch {
    my $config = shift;
    my $vf_hash = shift;
    my $transcript_cache = shift;
    my $include_regions = shift;
#MY
    my $vcf_hash = shift;
    my @finished_vcfs;
#MY
    
    my $up_down_size = MAX_DISTANCE_FROM_TRANSCRIPT;
    
    my (%vf_done, @finished_vfs, %seen_trs);
    
    # convert regions to cached sizes
    my $converted_regions = &convert_regions($config, $include_regions) if defined($config->{cache});
    
    foreach my $chr(sort {$a <=> $b} keys %$vf_hash) {
        if(defined($config->{standalone}) && !-e $config->{dir}.'/'.$chr) {
            debug("No cache found for chromsome $chr") unless defined($config->{quiet});
            next;
        }
        
        my $slice_cache;
        
        debug("Analyzing chromosome $chr") unless defined($config->{quiet});
        
        my $use_regions = defined($config->{cache}) ? $converted_regions : $include_regions;
        my ($count_from_db, $count_from_cache, $count_duplicates) = (0, 0, 0);
        
        if(!defined($transcript_cache->{$chr})) {
            
            # no regions defined (this probably shouldn't happen)
            if(!defined($use_regions->{$chr})) {
                
                # spoof regions covering whole chromosome
                my $start = 1;
                my $end = $config->{cache_region_size};
                my $slice = &get_slice($config, $chr);
                
                if(defined($slice)) {
                    while($start < $slice->end) {
                        push @{$use_regions->{$chr}}, $start.'-'.$end;
                        $start += $config->{cache_region_size};
                        $end += $config->{cache_region_size};
                    }
                }
            } #if(!defined($use_regions->{$chr}))
            
            # check we have defined regions
            if(defined($use_regions->{$chr})) {
                my $region_count = scalar @{$use_regions->{$chr}};
                my $counter;
                
                debug("Reading transcript data from cache and/or database") unless defined($config->{quiet});
                
                foreach my $region(sort {(split /\-/, $a)[0] <=> (split /\-/, $b)[1]} @{$use_regions->{$chr}}) {
                    &progress($config, $counter++, $region_count);
                    
                    # skip regions beyond the end of the chr
                    next if defined($slice_cache->{$chr}) && (split /\-/, $region)[0] > $slice_cache->{$chr}->length;
                    
                    # force quiet so other methods don't mess up the progress bar
                    my $quiet = $config->{quiet};
                    $config->{quiet} = 1;
                    
                    # try and load cache from disk if using cache
                    my $tmp_cache;
                    if(defined($config->{cache})) {
                        $tmp_cache = &load_dumped_transcript_cache($config, $chr, $region);
                        $count_from_cache += scalar @{$tmp_cache->{$chr}} if defined($tmp_cache->{$chr});
                    }
                    
                    # no cache found on disk or not using cache
                    if(!defined($tmp_cache->{$chr})) {
                        
                        if(defined($config->{standalone})) {
                            debug("WARNING: Could not find cache for $chr\:$region") unless defined($config->{quiet});
                            next;
                        }
                        
                        # spoof temporary region hash
                        my $tmp_hash;
                        push @{$tmp_hash->{$chr}}, $region;
                        
                        $tmp_cache = &cache_transcripts($config, $tmp_hash);
                        
                        # make it an empty arrayref that gets cached
                        # so we don't get confused and reload next time round
                        $tmp_cache->{$chr} ||= [];
                        
                        $count_from_db += scalar @{$tmp_cache->{$chr}};
                        
                        # dump to disk if writing to cache
                        &dump_transcript_cache($config, $tmp_cache, $chr, $region) if defined($config->{write_cache});
                    }
                    
                    # add loaded transcripts to main cache
                    if(defined($tmp_cache->{$chr})) {
                        while(my $tr = shift @{$tmp_cache->{$chr}}) {
                            
                            # track already added transcripts by dbID
                            my $dbID = $tr->dbID;
                            if($seen_trs{$dbID}) {
                                $count_duplicates++;
                                next;
                            }
                            $seen_trs{$dbID} = 1;
                            
                            push @{$transcript_cache->{$chr}}, $tr;
                        }
                    }
                    
                    undef $tmp_cache;
                    
                    # restore quiet status
                    $config->{quiet} = $quiet;
                    
                    # build slice cache
                    $slice_cache = &build_slice_cache($config, $transcript_cache) unless defined($slice_cache->{$chr});
                } #foreach my $region(sort {(split /\-/, $a)[0] <=> (split /\-/, $b)[1]} @{$use_regions->{$chr}})
                
                &end_progress($config);
            } #if(defined($use_regions->{$chr}))
        } #if(!defined($transcript_cache->{$chr}))
        
        # skip chr if no cache
        next unless defined($transcript_cache->{$chr});
        
        # copy slice from transcript to slice cache
        $slice_cache = &build_slice_cache($config, $transcript_cache) unless defined($slice_cache->{$chr});
        
        my $tr_count = scalar @{$transcript_cache->{$chr}};
        
        debug("Retrieved $tr_count transcripts ($count_from_cache cached, $count_from_db DB, $count_duplicates duplicates)") unless defined($config->{quiet});
        debug("Analyzing variants") unless defined($config->{quiet});
        
        my $tr_counter;
        
        while($tr_counter < $tr_count) {
            
            &progress($config, $tr_counter, $tr_count);
            
            my $tr = $transcript_cache->{$chr}->[$tr_counter++];
            
            # do each overlapping VF
            my $s = $tr->start - $up_down_size;
            my $e = $tr->end + $up_down_size;
            
            # get the chunks this transcript overlaps
            my %chunks;
            $chunks{$_} = 1 for (int($s/$config->{chunk_size})..int($e/$config->{chunk_size}));
            map {delete $chunks{$_} unless defined($vf_hash->{$chr}{$_})} keys %chunks;
            
            foreach my $chunk(keys %chunks) {
                foreach my $pos(grep {$_ >= $s && $_ <= $e} keys %{$vf_hash->{$chr}{$chunk}}) {
                    foreach my $vf(@{$vf_hash->{$chr}{$chunk}{$pos}}) {
                        
                        # pinch slice from slice cache if we don't already have it
                        $_->{slice} ||= $slice_cache->{$chr} for @{$vf_hash->{$chr}{$chunk}{$pos}};
                        
                        my $tv = Bio::EnsEMBL::Variation::TranscriptVariation->new(
                            -transcript => $tr,
                            -variation_feature => $vf,
                            -adaptor => $config->{tva},
                            -no_ref_check => 1
                        );
                        
                        # prefetching stuff here prevents doing loads at the
                        # end and makes progress reporting more useful
                        $tv->_prefetch_for_vep;
                        
                        $vf->add_TranscriptVariation($tv);
                    }
                }
            }
        } #while($tr_counter < $tr_count)
        
        # sort results into @finished_vfs array
        foreach my $chunk(sort {$a <=> $b} keys %{$vf_hash->{$chr}}) {
            foreach my $pos(sort {$a <=> $b} keys %{$vf_hash->{$chr}{$chunk}}) {
                
                # pinch slice from slice cache if we don't already have it
                $_->{slice} ||= $slice_cache->{$chr} for @{$vf_hash->{$chr}{$chunk}{$pos}};
                
                # add to final array
                push @finished_vfs, @{$vf_hash->{$chr}{$chunk}{$pos}};
#MY 
                push @finished_vcfs, @{$vcf_hash->{$chr}{$chunk}{$pos}};
#MY 
            }
        }
        
        &end_progress($config);
        
        debug("Calculating and writing output") unless defined($config->{quiet});
#MY
        #&print_consequences($config, \@finished_vfs);
        &print_consequences($config, \@finished_vfs, \@finished_vcfs);
        undef @finished_vfs;
        undef @finished_vcfs;
#MY
        
        # clean hash
        delete $vf_hash->{$chr};
#MY
        delete $vcf_hash->{$chr};
#MY
        
        delete $transcript_cache->{$chr} if defined($config->{cache});
    } #foreach my $chr(sort {$a <=> $b} keys %$vf_hash)
} #sub whole_genome_fetch

# gets existing VFs for a vf_hash
sub check_existing_hash {
    my $config = shift;
    my $vf_hash = shift;
    my $variation_cache;
    
    debug("Checking for existing variations") unless defined($config->{quiet});
    
    my ($chunk_count, $counter);
    $chunk_count += scalar keys %{$vf_hash->{$_}} for keys %{$vf_hash};
    
    foreach my $chr(keys %{$vf_hash}) {
        
        my %loaded_regions;
        
        foreach my $chunk(keys %{$vf_hash->{$chr}}) {
            &progress($config, $counter++, $chunk_count);
            
            # get the VFs for this chunk
            my ($start, $end);
            
            # work out start and end using chunk_size
            $start = $config->{chunk_size} * $chunk;
            $end = $config->{chunk_size} * ($chunk + 1);
            
            # using cache?
            if(defined($config->{cache})) {
                my $tmp_regions;
                push @{$tmp_regions->{$chr}}, $start.'-'.$end;
                
                my $converted_regions = &convert_regions($config, $tmp_regions);
                
                foreach my $region(@{$converted_regions->{$chr}}) {
                
                    unless($loaded_regions{$region}) {
                        my $tmp_cache = &load_dumped_variation_cache($config, $chr, $region);
                        
                        # load from DB if not found in cache
                        if(!defined($tmp_cache->{$chr})) {
                            if(defined($config->{standalone})) {
                                debug("WARNING: Could not find variation cache for $chr\:$region") unless defined($config->{quiet});
                                next;
                            }
                            
                            $tmp_cache->{$chr} = &get_variations_in_region($config, $chr, $region);
                            &dump_variation_cache($config, $tmp_cache, $chr, $region) if defined($config->{write_cache});
                        }
                        
                        # merge tmp_cache with the main cache
                        foreach my $key(keys %{$tmp_cache->{$chr}}) {
                            $variation_cache->{$chr}->{$key} = $tmp_cache->{$chr}->{$key};
                            delete $tmp_cache->{$chr}->{$key};
                        }
                        
                        # clear memory
                        undef $tmp_cache;
                        
                        # record this region as fetched
                        $loaded_regions{$region} = 1;
                    }
                }
            }
            
            # no cache, get all variations in region from DB
            else {
                $variation_cache->{$chr} = &get_variations_in_region($config, $chr, $start.'-'.$end);
            }
            
            # now compare retrieved vars with vf_hash
            foreach my $pos(keys %{$vf_hash->{$chr}->{$chunk}}) {
                foreach my $var(@{$vf_hash->{$chr}->{$chunk}->{$pos}}) {
                    my @found;
                    
                    if(defined($variation_cache->{$chr})) {
                        if(my $existing_vars = $variation_cache->{$chr}->{$pos}) {
                            foreach my $existing_var(@$existing_vars) {
                                push @found, $existing_var->[0] unless &is_var_novel($config, $existing_var, $var);
                            }
                        }
                    }
                    
                    $var->{existing} = join ",", @found;
                    $var->{existing} ||= '-';
                }
            }
        }
        
        delete $variation_cache->{$chr};
    }
    
    &end_progress($config);
}

# gets a slice from the slice adaptor
sub get_slice {
    my $config = shift;
    my $chr = shift;
    
    return undef unless defined($config->{sa}) && defined($chr);
    
    my $slice;
    
    # first try to get a chromosome
    eval { $slice = $config->{sa}->fetch_by_region('chromosome', $chr); };
    
    # if failed, try to get any seq region
    if(!defined($slice)) {
        $slice = $config->{sa}->fetch_by_region(undef, $chr);
    }
    
    return $slice;
}




# METHODS THAT DEAL WITH "REGIONS"
##################################

# scans file to get all slice bits we need
sub scan_file() {
    my $config = shift;
    
    my $in_file_handle = $config->{in_file_handle};
    
    my %include_regions;
    
    debug("Scanning input file") unless defined($config->{quiet});
    
    while(<$in_file_handle>) {
        chomp;
        
        # header line?
        next if /^\#/;
        
        # some lines (pileup) may actually parse out into more than one variant)
        foreach my $sub_line(@{&parse_line($config, $_)}) {
        
            # get the sub-line into named variables
            my ($chr, $start, $end, $allele_string, $strand, $var_name) = @{$sub_line};
            $chr =~ s/chr//ig unless $chr =~ /^chromosome$/i;
            $chr = 'MT' if $chr eq 'M';
            
            next if defined($config->{chr}) && !$config->{chr}->{$chr};
            
            $include_regions{$chr} ||= [];
            
            &add_region($start, $end, $include_regions{$chr});
        }
    }
    
    # close filehandle and recycle
    close $in_file_handle;
    $config->{in_file_handle} = &get_in_file_handle($config);
    
    # merge regions
    &merge_regions(\%include_regions);
    
    return \%include_regions;
}

# gets regions from VF hash
sub regions_from_hash {
    my $config = shift;
    my $vf_hash = shift;
    
    my %include_regions;
    
    # if using cache we just want the regions of cache_region_size
    # since that's what we'll get from the cache (or DB if no cache found)
    if(defined($config->{cache})) {
        
        my $region_size = $config->{cache_region_size};
        
        foreach my $chr(keys %$vf_hash) {
            $include_regions{$chr} = [];
            my %temp_regions;
            
            foreach my $chunk(keys %{$vf_hash->{$chr}}) {
                foreach my $pos(keys %{$vf_hash->{$chr}{$chunk}}) {
                    my ($s, $e) = ($pos - MAX_DISTANCE_FROM_TRANSCRIPT, $pos + MAX_DISTANCE_FROM_TRANSCRIPT);
                    
                    my $low = int ($s / $region_size);
                    my $high = int ($e / $region_size) + 1;
                    
                    for my $i($low..($high - 1)) {
                        $temp_regions{(($i * $region_size) + 1).'-'.(($i + 1) * $region_size)} = 1;
                    }
                }
            }
            
            @{$include_regions{$chr}} = keys %temp_regions;
        }
    }
    
    # if no cache we don't want to fetch more than is necessary, so find the
    # minimum covered region of the variations in the hash
    else {
        foreach my $chr(keys %$vf_hash) {
            $include_regions{$chr} = [];
            
            foreach my $chunk(keys %{$vf_hash->{$chr}}) {
                foreach my $pos(keys %{$vf_hash->{$chr}{$chunk}}) {
                    &add_region($_->start, $_->end, $include_regions{$chr}) for @{$vf_hash->{$chr}{$chunk}{$pos}};
                }
            }
        }
        
        # merge regions
        &merge_regions(\%include_regions);
    }
    
    return \%include_regions;
}

# adds a region to region list, expanding existing one if overlaps
sub add_region {
    my $start = shift;
    my $end = shift;
    my $region_list = shift;
    
    # fix end for insertions
    $end = $start if $end < $start;
    
    my $added = 0;
    my $i = 0;
    
    while ($i < scalar @$region_list) {
        my ($region_start, $region_end) = split /\-/, $region_list->[$i];
        
        if($start <= $region_end && $end >= $region_start) {
            my $new_region_start = ($start < $end ? $start : $end) - MAX_DISTANCE_FROM_TRANSCRIPT;
            my $new_region_end = ($start > $end ? $start : $end) + MAX_DISTANCE_FROM_TRANSCRIPT;
            
            $region_start = $new_region_start if $new_region_start < $region_start;
            $region_end = $new_region_end if $new_region_end > $region_end;
            
            $region_list->[$i] = $region_start.'-'.$region_end;
            $added = 1;
        }
        
        $i++;
    }
    
    unless($added) {
        push @{$region_list}, ($start - MAX_DISTANCE_FROM_TRANSCRIPT).'-'.($end + MAX_DISTANCE_FROM_TRANSCRIPT);
    }
}

# merges overlapping regions from scans
sub merge_regions {
    my $include_regions = shift;
    
    # now merge overlapping regions
    foreach my $chr(keys %$include_regions) {
        my $max_index = $#{$include_regions->{$chr}};
        my (@new_regions, %skip);
        
        for my $i(0..$max_index) {
            next if $skip{$i};
            my ($s, $e) = split /\-/, $include_regions->{$chr}[$i];
            
            for my $j(($i+1)..$max_index) {
                next if $skip{$j};
                my ($ns, $ne) = split /\-/, $include_regions->{$chr}[$j];
                
                if($s <= $ne && $e >= $ns) {
                    $s = $ns if $ns < $s;
                    $e = $ne if $ne > $e;
                    
                    $skip{$j} = 1;
                }
            }
            
            push @new_regions, $s.'-'.$e;
        }
        
        # replace original
        $include_regions->{$chr} = \@new_regions;
        
        $config->{region_count} += scalar @new_regions;
    }
    
    return $include_regions;
}

# converts regions as determined by scan_file to regions loadable from cache
sub convert_regions {
    my $config = shift;
    my $regions = shift;
    
    return undef unless defined $regions;
    
    my $region_size = $config->{cache_region_size};
    
    my %new_regions;
    
    foreach my $chr(keys %$regions) {
        my %temp_regions;
        
        foreach my $region(@{$regions->{$chr}}) {
            my ($s, $e) = split /\-/, $region;
            
            my $low = int ($s / $region_size);
            my $high = int ($e / $region_size) + 1;
            
            for my $i($low..($high - 1)) {
                $temp_regions{(($i * $region_size) + 1).'-'.(($i + 1) * $region_size)} = 1;
            }
        }
        
        @{$new_regions{$chr}} = keys %temp_regions;
    }
    
    return \%new_regions;
}





# CACHE METHODS
###############

# get transcripts for slices
sub cache_transcripts {
    my $config = shift;
    my $include_regions = shift;
    
    my $transcript_cache;
    my $i;
    
    debug("Caching transcripts") unless defined($config->{quiet});
    
    foreach my $chr(keys %$include_regions) {
        
        my $slice = &get_slice($config, $chr);
        
        next unless defined $slice;
        
        # prefetch some things
        $slice->is_circular;
        
        # trim bumf off the slice
        delete $slice->{coord_system}->{adaptor} if defined($config->{write_cache});
        
        # no regions?
        if(!scalar @{$include_regions->{$chr}}) {
            my $start = 1;
            my $end = $config->{cache_region_size};
            
            while($start < $slice->end) {
                push @{$include_regions->{$chr}}, $start.'-'.$end;
                $start += $config->{cache_region_size};
                $end += $config->{cache_region_size};
            }
        }
        
        my $region_count;
        
        if(scalar keys %$include_regions == 1) {
            my ($chr) = keys %$include_regions;
            $region_count = scalar @{$include_regions->{$chr}};
            debug("Caching transcripts for chromosome $chr") unless defined($config->{quiet});
        }
        
        foreach my $region(@{$include_regions->{$chr}}) {
            &progress($config, $i++, $region_count || $config->{region_count});
            
            my ($s, $e) = split /\-/, $region;
            
            # sanity check start and end
            $s = 1 if $s < 1;
            $e = $slice->end if $e > $slice->end;
            
            # get sub-slice
            my $sub_slice = $slice->sub_Slice($s, $e);
            
            # add transcripts to the cache, via a transfer to the chrom's slice
            if(defined($sub_slice)) {
                foreach my $gene(@{$sub_slice->get_all_Genes(undef, undef, 1)}) {
                    my $gene_stable_id = $gene->stable_id;
                    
                    foreach my $tr(map {$_->transfer($slice)} @{$gene->get_all_Transcripts}) {
                        $tr->{_gene_stable_id} = $gene_stable_id;
                        
                        if(defined($config->{prefetch})) {
                            $tr->{_gene} = $gene;
                            &prefetch_transcript_data($config, $tr);
                            delete $tr->{_gene};
                        }
                        
                        # strip some unnecessary data from the transcript object
                        &clean_transcript($tr) if defined($config->{write_cache});
                        
                        push @{$transcript_cache->{$chr}}, $tr;
                    }
                }
            }
        }
    }
    
    &end_progress($config);
    
    return $transcript_cache;
}

# gets rid of extra bits of info attached to the transcript that we don't need
sub clean_transcript {
    my $tr = shift;
    
    foreach my $key(qw(display_xref external_db external_display_name external_name external_status created_date status description edits_enabled modified_date)) {
        delete $tr->{$key} if defined($tr->{$key});
    }
    
    # clean all attributes but miRNA
    if(defined($tr->{attributes})) {
        my @new_atts;
        foreach my $att(@{$tr->{attributes}}) {
            push @new_atts, $att if $att->{code} eq 'miRNA';
        }
        $tr->{attributes} = \@new_atts;
    }
    
    $tr->{analysis} = {};
}

# build slice cache from transcript cache
sub build_slice_cache {
    my $config = shift;
    my $transcript_cache = shift;
    
    my %slice_cache;
    
    foreach my $chr(keys %$transcript_cache) {
        $slice_cache{$chr} = $transcript_cache->{$chr}[0]->slice;
        
        # reattach adaptor to the coord system
        $slice_cache{$chr}->{coord_system}->{adaptor} ||= $config->{csa};
    }
    
    return \%slice_cache;
}

# pre-fetches per-transcript data
sub prefetch_transcript_data {
    my $config = shift;
    my $tran = shift;
    
    # introns, translateable_seq, mapper
    $tran->{_variation_effect_feature_cache}->{introns} ||= $tran->get_all_Introns;
    $tran->{_variation_effect_feature_cache}->{translateable_seq} ||= $tran->translateable_seq;
    $tran->{_variation_effect_feature_cache}->{mapper} ||= $tran->get_TranscriptMapper;
    
    # peptide
    unless ($tran->{_variation_effect_feature_cache}->{peptide}) {
        my $translation = $tran->translate;
        $tran->{_variation_effect_feature_cache}->{peptide} = $translation ? $translation->seq : undef;
    }
    
    # codon table
    unless ($tran->{_variation_effect_feature_cache}->{codon_table}) {
        # for mithocondrial dna we need to to use a different codon table
        my $attrib = $tran->slice->get_all_Attributes('codon_table')->[0];
        
        $tran->{_variation_effect_feature_cache}->{codon_table} = $attrib ? $attrib->value : 1;
    }
    
    # gene HGNC
    if(defined $config->{hgnc}) {
        # get from gene cache if found already
        if(defined($tran->{_gene}->{_hgnc})) {
            $tran->{_gene_hgnc} = $tran->{_gene}->{_hgnc};
        }
        else {
            my @entries = grep {$_->database eq 'HGNC'} @{$tran->{_gene}->get_all_DBEntries()};
            if(scalar @entries) {
                $tran->{_gene_hgnc} = $entries[0]->display_id;
            }
            
            $tran->{_gene_hgnc} ||= '-';
            
            # cache it on the gene object too
            $tran->{_gene}->{_hgnc} = $tran->{_gene_hgnc};
        }
    }
    
    return $tran;
}

# dumps out transcript cache to file
sub dump_transcript_cache {
    my $config = shift;
    my $transcript_cache = shift;
    my $chr = shift;
    my $region = shift;
    
    debug("Dumping cached transcript data") unless defined($config->{quiet});
    
    # clean the slice adaptor before storing
    &clean_slice_adaptor($config);
    
    &strip_transcript_cache($config, $transcript_cache);
    
    $config->{reg}->disconnect_all;
    
    my $dir = $config->{dir}.'/'.$chr;
    my $dump_file = $dir.'/'.($region || "dump").'.gz';
    
    # make directory if it doesn't exist
    if(!(-e $dir)) {
        system("mkdir -p ".$dir);
    }
    
    debug("Writing to $dump_file") unless defined($config->{quiet});
    
    # storable
    open my $fh, "| gzip -c > ".$dump_file or die "ERROR: Could not write to dump file $dump_file";
    nstore_fd($transcript_cache, $fh);
    close $fh;
}

# loads in dumped transcript cache to memory
sub load_dumped_transcript_cache {
    my $config = shift;
    my $chr = shift;
    my $region = shift;
    
    my $dir = $config->{dir}.'/'.$chr;
    my $dump_file = $dir.'/'.($region || "dump").'.gz';
    
    return undef unless -e $dump_file;
    
    debug("Reading cached transcript data for chromosome $chr".(defined $region ? "\:$region" : "")." from dumped file") unless defined($config->{quiet});
    
    open my $fh, $config->{compress}." ".$dump_file." |" or return undef;
    my $transcript_cache = fd_retrieve($fh);
    close $fh;
    
    return $transcript_cache;
}

# strips cache
sub strip_transcript_cache {
    my $config = shift;
    my $cache = shift;
    
    foreach my $chr(keys %$cache) {
        foreach my $tr(@{$cache->{$chr}}) {
            foreach my $exon(@{$tr->{_trans_exon_array}}) {
                delete $exon->{adaptor};
                delete $exon->{slice}->{adaptor};
            }
            
            delete $tr->{adaptor};
            delete $tr->{slice}->{adaptor};
        }
    }
}

# cleans slice adaptor before storing in cache
sub clean_slice_adaptor{
    my $config = shift;
    
    # clean some stuff off the slice adaptor
    $config->{sa}->{asm_exc_cache} = {};
    $config->{sa}->{sr_name_cache} = {};
    $config->{sa}->{sr_id_cache} = {};
    delete $config->{sa}->{db}->{seq_region_cache};
    delete $config->{sa}->{db}->{name_cache};
}


# dump adaptors to cache
sub dump_adaptor_cache {
    my $config = shift;
    
    $config->{reg}->disconnect_all;
    
    my $dir = $config->{dir};
    my $dump_file = $dir.'/adaptors.gz';
    
    # make directory if it doesn't exist
    if(!(-e $dir)) {
        system("mkdir -p ".$dir);
	}
	
    open my $fh, "| gzip -c > ".$dump_file or die "ERROR: Could not write to dump file $dump_file";
    nstore_fd($config, $fh);
    close $fh;
}

# load dumped adaptors
sub load_dumped_adaptor_cache {
    my $config = shift;
    
    my $dir = $config->{dir};
    my $dump_file = $dir.'/adaptors.gz';
    
    return undef unless -e $dump_file;
    
    debug("Reading cached adaptor data") unless defined($config->{quiet});
    
    open my $fh, $config->{compress}." ".$dump_file." |" or return undef;
    my $cached_config = fd_retrieve($fh);
    close $fh;
    
    $config->{$_} = $cached_config->{$_} for qw(sa ga ta vfa tva mca csa);
    
    return 1;
}

# dumps cached variations to disk
sub dump_variation_cache {
    my $config = shift;
    my $v_cache = shift;
    my $chr = shift;
    my $region = shift;
    
    my $dir = $config->{dir}.'/'.$chr;
    my $dump_file = $dir.'/'.($region || "dump").'_var.gz';
    
    # make directory if it doesn't exist
    if(!(-e $dir)) {
        system("mkdir -p ".$dir);
    }
    
    open DUMP, "| gzip -c > ".$dump_file or die "ERROR: Could not write to adaptor dump file $dump_file";
    
    foreach my $pos(keys %{$v_cache->{$chr}}) {
        foreach my $v(@{$v_cache->{$chr}->{$pos}}) {
            my ($name, $source, $start, $end, $as, $strand) = @$v;
            
            print DUMP join " ", (
                $name,
                $source == 1 ? '' : $source,
                $start,
                $end == $start ? '' : $end,
                $as,
                $strand == 1 ? '' : $strand,
            );
            print DUMP "\n";
        }
    }
    
    close DUMP;    
}

# loads dumped variation cache
sub load_dumped_variation_cache {
    my $config = shift;
    my $chr = shift;
    my $region = shift;
    
    my $dir = $config->{dir}.'/'.$chr;
    my $dump_file = $dir.'/'.($region || "dump").'_var.gz';
    
    return undef unless -e $dump_file;
    
    open DUMP, $config->{compress}." ".$dump_file." |" or return undef;
    
    my $v_cache;
    
    while(<DUMP>) {
        chomp;
        my ($name, $source, $start, $end, $as, $strand) = split / /, $_;
        $source ||= 1;
        $end ||= $start;
        $strand ||= 1;
        
        my @v = ($name, $source, $start, $end, $as, $strand);
        push @{$v_cache->{$chr}->{$start}}, \@v;
    }
    
    close DUMP;
    
    return $v_cache;
}

# builds a full cache for this species
sub build_full_cache {
    my $config = shift;
    my $rebuild = shift;
    
    my @slices;
    
    if($config->{build} =~ /all/i) {
        @slices = @{$config->{sa}->fetch_all('toplevel')};
    }
    else {
        foreach my $val(split /\,/, $config->{build}) {
            my @nnn = split /\-/, $val;
            
            foreach my $chr($nnn[0]..$nnn[-1]) {
                my $slice = &get_slice($config, $chr);
                push @slices, $slice if defined($slice);
            }
        }
    }
    
    foreach my $slice(@slices) {
        my $chr = $slice->seq_region_name;
        
        my $regions;
        
        # for progress
        my $region_count = int($slice->end / $config->{cache_region_size}) + 1;
        my $counter = 0;
        
        # initial region
        my ($start, $end) = (1, $config->{cache_region_size});
        
        debug((defined($config->{rebuild}) ? "Rebuild" : "Creat")."ing cache for chromosome $chr") unless defined($config->{quiet});
        
        while($start < $slice->end) {
            
            &progress($config, $counter++, $region_count);
            
            # store quiet status
            my $quiet = $config->{quiet};
            $config->{quiet} = 1;
            
            # store transcripts
            $regions->{$chr} = [$start.'-'.$end];
            my $tmp_cache = ($rebuild ? &load_dumped_transcript_cache($config, $chr, $start.'-'.$end) : &cache_transcripts($config, $regions));
            $tmp_cache->{$chr} ||= [];
            
            &dump_transcript_cache($config, $tmp_cache, $chr, $start.'-'.$end);
            undef $tmp_cache;
            
            # store variations
            my $variation_cache;
            $variation_cache->{$chr} = &get_variations_in_region($config, $chr, $start.'-'.$end);
            $variation_cache->{$chr} ||= {};
            
            &dump_variation_cache($config, $variation_cache, $chr, $start.'-'.$end);
            undef $variation_cache;
            
            # restore quiet status
            $config->{quiet} = $quiet;
            
            # increment by cache_region_size to get next region
            $start += $config->{cache_region_size};
            $end += $config->{cache_region_size};
        }
        
        &end_progress($config);
        
        undef $regions;
    }
}

# format coords for printing
sub format_coords {
    my ($start, $end) = @_;
    
    if(!defined($start)) {
        return '-';
    }
    elsif(!defined($end)) {
        return $start;
    }
    elsif($start == $end) {
        return $start;
    }
    elsif($start > $end) {
        return $end.'-'.$start;
    }
    else {
        return $start.'-'.$end;
    }
}




# METHODS TO FIND CO-LOCATED / EXISTING VARIATIONS
##################################################

# finds an existing VF in the db
sub find_existing {
    my $config = shift;
    my $new_vf = shift;
    
    if(defined($new_vf->adaptor->db)) {
        
        my $sth = $new_vf->adaptor->db->dbc->prepare(qq{
            SELECT variation_name, source_id, seq_region_start, seq_region_end, allele_string, seq_region_strand
            FROM variation_feature
            WHERE seq_region_id = ?
            AND seq_region_start = ?
            AND seq_region_end = ?
            ORDER BY source_id ASC
        });
        
        $sth->execute($new_vf->slice->get_seq_region_id, $new_vf->start, $new_vf->end);
        
        my @v;
        for my $i(0..5) {
            $v[$i] = undef;
        }
        
        $sth->bind_columns(\$v[0], \$v[1], \$v[2], \$v[3], \$v[4], \$v[5]);
        
        my @found;
        
        while($sth->fetch) {
            push @found, $v[0] unless &is_var_novel($config, \@v, $new_vf);
        }
        
        $sth->finish();
        
        return (scalar @found ? join ",", @found : undef);
    }
    
    return undef;
}

# compare a new vf to one from the cache / DB
sub is_var_novel {
    my $config = shift;
    my $existing_var = shift;
    my $new_var = shift;
    
    my $is_novel = 1;
    
    $is_novel = 0 if $existing_var->[2] == $new_var->start && $existing_var->[3] == $new_var->end;
    
    if(defined($config->{check_alleles})) {
        my %existing_alleles;
        
        $existing_alleles{$_} = 1 for split /\//, $existing_var->[4];
        
        my $seen_new = 0;
        foreach my $a(split /\//, $new_var->allele_string) {
            reverse_comp(\$a) if $new_var->seq_region_strand ne $existing_var->[5];
            $seen_new = 1 unless defined $existing_alleles{$a};
        }
        
        $is_novel = 1 if $seen_new;
    }
    
    return $is_novel;
}

# gets all variations in a region
sub get_variations_in_region {
    my $config = shift;
    my $chr = shift;
    my $region = shift;
    
    my ($start, $end) = split /\-/, $region;
    
    my %variations;
    
    if(defined($config->{vfa}->db)) {
        my $sth = $config->{vfa}->db->dbc->prepare(qq{
            SELECT vf.variation_name, vf.source_id, vf.seq_region_start, vf.seq_region_end, vf.allele_string, vf.seq_region_strand
            FROM variation_feature vf, seq_region s
            WHERE s.seq_region_id = vf.seq_region_id
            AND s.name = ?
            AND vf.seq_region_start >= ?
            AND vf.seq_region_start <= ?
        });
        
        $sth->execute($chr, $start, $end);
        
        my @v;
        for my $i(0..5) {
            $v[$i] = undef;
        }
        
        $sth->bind_columns(\$v[0], \$v[1], \$v[2], \$v[3], \$v[4], \$v[5]);
        
        while($sth->fetch) {
            my @v_copy = @v;
            push @{$variations{$v[2]}}, \@v_copy;
        }
        
        $sth->finish();
    }
    
    return \%variations;
}




# DEBUG AND STATUS METHODS
##########################

# gets time
sub get_time() {
    my @time = localtime(time());

    # increment the month (Jan = 0)
    $time[4]++;

    # add leading zeroes as required
    for my $i(0..4) {
        $time[$i] = "0".$time[$i] if $time[$i] < 10;
    }

    # put the components together in a string
    my $time =
         ($time[5] + 1900)."-".
         $time[4]."-".
         $time[3]." ".
        $time[2].":".
        $time[1].":".
        $time[0];

    return $time;
}

# prints debug output with time
sub debug {
    my $text = (@_ ? (join "", @_) : "No message");
    my $time = get_time;
    
    print $time." - ".$text.($text =~ /\n$/ ? "" : "\n");
}

# update or initiate progress bar
sub progress {
    my ($config, $i, $total) = @_;
    
    return if defined($config->{quiet}) || defined($config->{no_progress});
    
    my $width = $config->{terminal_width};
    my $percent = int(($i/$total) * 100);
    my $numblobs = (($i/$total) * $width) - 2;
    
    # this ensures we're not writing to the terminal too much
    return if(defined($config->{prev_prog})) && $numblobs.'-'.$percent eq $config->{prev_prog};
    $config->{prev_prog} = $numblobs.'-'.$percent;
    
    printf("\r% -${width}s% 1s% 10s", '['.('=' x $numblobs).($numblobs == $width - 2 ? '=' : '>'), ']', "[ " . $percent . "% ]");
}

# end progress bar
sub end_progress {
    my $config = shift;
    return if defined($config->{quiet}) || defined($config->{no_progress});
    &progress($config, 1,1);
    print "\n";
    delete $config->{prev_prog};
}

# outputs usage message
sub usage {
    my $usage =<<END;
#----------------------------------#
# ENSEMBL VARIANT EFFECT PREDICTOR #
#----------------------------------#

version $VERSION

By Will McLaren (wm2\@ebi.ac.uk)

http://www.ensembl.org/info/docs/variation/vep/vep_script.html

Usage:
perl variant_effect_predictor.pl [arguments]

Options
=======

--help                 Display this message and quit
--verbose              Display verbose output as the script runs [default: off]
--quiet                Suppress status and warning messages [default: off]
--no_progress          Suppress progress bars [default: off]

--config               Load configuration from file. Any command line options
                       specified overwrite those in the file [default: off]

-i | --input_file      Input file - if not specified, reads from STDIN. Files
                       may be gzip compressed.
--format               Alternative input file format - one of "pileup", "vcf"
-o | --output_file     Output file. Write to STDOUT by specifying -o STDOUT - this
                       will force --quiet [default: "variant_effect_output.txt"]
--force_overwrite      Force overwriting of output file [default: quit if file
                       exists]

-t | --terms           Type of consequence terms to output - one of "ensembl", "SO",
                       "NCBI" [default: ensembl]
 
--sift=[p|s|b]         Add SIFT [p]rediction, [s]core or [b]oth [default: off]
--polyphen=[p|s|b]     Add PolyPhen [p]rediction, [s]core or [b]oth [default: off]
--condel=[p|s|b]       Add Condel SIFT/PolyPhen consensus [p]rediction, [s]core or
                       [b]oth [default: off]

NB: SIFT, PolyPhen and Condel predictions are currently available for human only

--regulatory           Look for overlaps with regulatory regions. The script can
                       also call if a variant falls in a high information position
                       within a transcription factor binding site. Output lines have
                       a Feature type of RegulatoryFeature or MotifFeature. Requires
                       database connection. [default: off]
                       
NB: Regulatory consequences are currently available for human and mouse only

--hgnc                 If specified, HGNC gene identifiers are output alongside the
                       Ensembl Gene identifier [default: off]
--hgvs                 Output HGVS identifiers (coding and protein). Requires database
                       connection [default: off]
--protein              Output Ensembl protein identifer [default: off]
--gene                 Force output of Ensembl gene identifer - disabled by default
                       unless using --cache or --no_whole_genome [default: off]

--coding_only          Only return consequences that fall in the coding region of
                       transcripts [default: off]
--most_severe          Ouptut only the most severe consequence per variation.
                       Transcript-specific columns will be left blank. [default: off]
--summary              Ouptut only a comma-separated list of all consequences per
                       variation. Transcript-specific columns will be left blank.
                       [default: off]

--check_ref            If specified, checks supplied reference allele against stored
                       entry in Ensembl Core database [default: off]
--check_existing       If specified, checks for existing co-located variations in the
                       Ensembl Variation database [default: off]
--check_alleles        If specified, the alleles of existing co-located variations
                       are compared to the input; an existing variation will only
                       be reported if no novel allele is in the input (strand is
                       accounted for) [default: off]

--chr [list]           Select a subset of chromosomes to analyse from your file. Any
                       data not on this chromosome in the input will be skipped. The
                       list can be comma separated, with "-" characters representing
                       an interval [default: off]
--gp                   If specified, tries to read GRCh37 position from GP field in the
                       INFO column of a VCF file. Only applies when VCF is the input
                       format and human is the species [default: off]

--species              Species to use [default: "human"]
--host                 Manually define database host [default: "ensembldb.ensembl.org"]
-u | --user            Database username [default: "anonymous"]
--port                 Database port [default: 5306]
--password             Database password [default: no password]
--genomes              Sets DB connection params for Ensembl Genomes [default: off]
--registry             Registry file to use defines DB connections [default: off]
                       Defining a registry file overrides above connection settings.
--db_version=[number]  Force script to load DBs from a specific Ensembl version. Not
                       advised due to likely incompatibilities between API and DB

--no_whole_genome      Run in old-style, non-whole genome mode [default: off]
--buffer_size          Sets the number of variants sent in each batch [default: 5000]
                       Increasing buffer size can retrieve results more quickly
                       but requires more memory. Only applies to whole genome mode.
                       
--cache                Enables read-only use of cache [default: off]
--dir [directory]      Specify the base cache directory to use [default: "\$HOME/.vep/"]
--write_cache          Enable writing to cache [default: off]
--build [all|list]     Build a complete cache for the selected species. Build for all
                       chromosomes with --build all, or a list of chromosomes (see
                       --chr). DO NOT USE WHEN CONNECTED TO PUBLIC DB SERVERS AS THIS
                       VIOLATES OUR FAIR USAGE POLICY [default: off]
                       
--compress             Specify utility to decompress cache files - may be "gzcat" or
                       "gzip -dc Only use if default does not work [default: zcat]
                       
--skip_db_check        ADVANCED! Force the script to use a cache built from a different
                       database than specified with --host. Only use this if you are
                       sure the hosts are compatible (e.g. ensembldb.ensembl.org and
                       useastdb.ensembl.org) [default: off]
--cache_region_size    ADVANCED! The size in base-pairs of the region covered by one
                       file in the cache. [default: 1MB]
END

    print $usage;
}

#MY
sub get_kegg_from_vw{
   
    my ($ensg) = @_;
    
    my @data;
    
    my $query = "CALL $config->{vw_database}.ensg2kegg('$ensg')";
    my $qh = $config->{vw_conn}->prepare($query);
    $qh->execute() or die "Unable to execute $query: $DBI::errstr\n";
    while (my @row = $qh->fetchrow_array()){
        if (substr($row[0],0,3) eq 'hsa'){
            print;
        }
        push @data, $row[0];
    }
   
   return \@data;
}

sub get_omim_from_vw{

    my ($ensg) = @_;
    
    my @data;
    
    my $query = "CALL $config->{vw_database}.ensg2omim('$ensg')";
    my $qh = $config->{vw_conn}->prepare($query);
    $qh->execute() or die "Unable to execute $query: $DBI::errstr\n";
    while (my @row = $qh->fetchrow_array()){
        push @data, $row[0];
    }
   
   return \@data;
}

sub get_uniprot_from_vw{
    
    my ($enst) = @_;
    
    my %data;
    
    my $query = "CALL $config->{vw_database}.enst2uniprot('$enst')";
    my $qh = $config->{vw_conn}->prepare($query);
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
    
    my $query = "CALL $config->{vw_database}.enst2uniprot_feature('$enst', $aaStart, $aaEnd)";
    my $qh = $config->{vw_conn}->prepare($query);
    $qh->execute() or die "Unable to execute $query: $DBI::errstr\n";
    while (my @row = $qh->fetchrow_array()){
        #$data{feature} = (aaStart,aaEnd,description)
        my @feature = @row[1..3];
        $data{$row[0]} = \@feature;
    }
   
   return \%data;
}

#connect to vw database
sub connect_to_vw(){
    my $config = shift;
    my $dsn = "dbi:$config->{vw_platform}:$config->{vw_database}:$config->{vw_host}:$config->{port}";
    my $vw_conn = DBI->connect($dsn, $config->{vw_user}, $config->{vw_password})
        or die "Unable to connect: $DBI::errstr\n";
    return $vw_conn;
}

#VariationFeature
sub my_vf_annotations(){
    my $config = shift;
    my $line = shift;
    my $new_vf = shift;
}

#RegulatoryFeature
sub my_rf_annotations(){
    my $config = shift;
    my $line = shift;
    my $rf = shift;
}

#RegulatoryFeatureVariation
sub my_rfv_annotations(){
    my $config = shift;
    my $line = shift;
    my $rfv = shift;
}

#RegulatoryFeatureVariationAllele
sub my_rfva_annotations(){
    my $config = shift;
    my $line = shift;
    my $rfva = shift;
}

#MotifFeature
sub my_mf_annotations(){
    my $config = shift;
    my $line = shift;
    my $mf = shift;
}

#MotifFeatureVariation
sub my_mfv_annotations(){
    my $config = shift;
    my $line = shift;
    my $mfv = shift;
}

#MotifFeatureVariationAllele
sub my_mfva_annotations(){
    my $config = shift;
    my $line = shift;
    my $mfva = shift;
}

#TranscriptVariation
sub my_tv_annotations(){
    my $config = shift;
    my $line = shift;
    my $tv = shift;
}

#TranscriptVariationAllele
sub my_tva_annotations(){
    my $config = shift;
    my $line = shift;
    debug("my_tva_annotations $line->{Location}") unless defined($config->{quiet});
    my $tv = shift;
    my $tva = shift;
    my $vf = $tv->{variation_feature};
    my $chrom = $vf->seq_region_name;
    my $chrom_start = $vf->start;
    my $chrom_end = $vf->end;
    my $chrom_strand = $vf->strand;
    my %genomic_coords = (
        chr    => $chrom,
        start  => $chrom_start,
        end    => $chrom_end,
        strand => $chrom_strand
    );
    my $reference_allele = $tv->{reference_allele}->variation_feature_seq;
    my $non_reference_allele = $tva->variation_feature_seq;
    my $transcript = $tv->transcript;
    my $transcript_id = $transcript->{stable_id};
    my $gene = $config->{ga}->fetch_by_transcript_stable_id($transcript_id);
    my $DAS = $gene->get_all_DBEntries;
    my $is_canonical = $transcript_id eq $gene->canonical_transcript()->stable_id;
    $line->{Canonical_Transcript} = $is_canonical ? 1 : 0;
    my $gene_id = $gene->{stable_id};
    my $gene_description = $gene->{description};
    $line->{Gene_Description} = $gene_description;
    $line->{KEGG_Pathway} = join('|', @{get_kegg_from_vw($gene_id)});
    $line->{OMIM_Disorder} = join('|', @{get_omim_from_vw($gene_id)});
    my $uniprot = get_uniprot_from_vw($transcript_id);
    map {$line->{$_} = $uniprot->{$_} if $_ ne 'ID'} keys %{$uniprot};
    my $translation = $transcript->translation;
    my $protein_id;
    my $protein_sequence;
    my $protein_length;
    my $altered_aa_start;
    my $altered_aa_end;
    my $altered_base_start;
    my $altered_base_end;
    my $pep_allele_string;
    my $amino_acid_reference;
    my $amino_acid_variant;
    my $changes_protein = 0;
    my $entrez_gene_name;
    my $entrez_gene_id;
    my @uniprotkb_ac;
    my @uniprot_id;
    my @omim_id;
    my @hgmd_id;
    my @go_names = ();
    my $xrefs = $gene->get_all_DBLinks(); 
    foreach my $xref (@{$xrefs}) {
        if ($xref->dbname() eq 'EntrezGene') {
            $entrez_gene_name = $xref->display_id();
            $entrez_gene_id = $xref->primary_id();
            last;
        }
    }
    $line->{Entrez_Gene_Name} = $entrez_gene_name;
    $line->{Entrez_Gene_ID} = $entrez_gene_id;
    $xrefs = $transcript->get_all_DBLinks(); 
    foreach my $xref (@{$xrefs}) {
        if ($xref->dbname() eq 'Uniprot/SWISSPROT') {
            push(
                @uniprotkb_ac,
                $xref->primary_id()
            );
            push(
                @uniprot_id,
                $xref->display_id()
            );
        }
        if ($xref->dbname() =~/^MIM/) {
            push(
                @omim_id,
                $xref->display_id()
            );
        }
        if ($xref->dbname() =~ /HGMD/) {
            push(
                @hgmd_id,
                $xref->display_id()
            );
        }
        elsif ($xref->dbname() eq 'goslim_goa') {
            my $go_id = $xref->display_id();
            my $go_term;
            my $go_name;
            my $go_definition;
            if (defined($config->{goa})) {
                $go_term = $config->{goa}->fetch_by_accession($go_id);
                $go_name = $go_term->name();
                $go_definition = $go_term->definition();
            }
            if (defined($go_name)) {
                push(@go_names, "[$go_id]:$go_name");
            }
            else {
                push(@go_names, "$go_id");
            }
        }
    }
    @uniprotkb_ac = @{get_unique(\@uniprotkb_ac)};
    $line->{UniProtKB_AC} = $uniprotkb_ac[0];
    @uniprot_id = @{get_unique(\@uniprot_id)};
    $line->{UniProt_ID} = join(',',@uniprot_id);
    @omim_id = @{get_unique(\@omim_id)};
    #$line->{OMIM} = join(',',@omim_id);
    @hgmd_id = @{get_unique(\@hgmd_id)};
    $line->{HGMD} = join(',',@hgmd_id);
    my $gene_ontology = join(';', @{get_unique(\@go_names)});
    $line->{Gene_Ontology} = $gene_ontology;
    my $Allele= $tva->{variation_feature_seq};
    my $Amino_acids = $tva->pep_allele_string;
    my $Codons = $tva->display_codon_allele_string;
    my $term_method = $config->{terms}.'_term';
    my @c_so = map {$_->SO_term} @{$tva->get_all_OverlapConsequences};
    my @consequences = map {$_->$term_method || $_->SO_term} @{$tva->get_all_OverlapConsequences};
    my @ranks = map {$_->rank} @{$tva->get_all_OverlapConsequences};
    my %consequences_ranks_so;
    my %consequences_ranks;
    map {$consequences_ranks_so{$c_so[$_]} = $ranks[$_]} 0..$#ranks;
    map {$consequences_ranks{$consequences[$_]} = $ranks[$_]} 0..$#ranks;
    my @consequences_sorted = sort {$consequences_ranks{$a} <=> $consequences_ranks{$b}} keys %consequences_ranks;
    my $consequence_severest = $consequences_sorted[0];
    my $consequence_severest_rank = $consequences_ranks{$consequence_severest};
    $line->{Consequence} = $consequence_severest;
    $line->{Consequence_rank} = $consequence_severest_rank;
    $line->{Consequences_all} = join ",", @consequences_sorted;
    if(defined($translation)){
        $protein_id = $translation->{stable_id};
        $protein_sequence = $translation->{seq};
        $protein_length = length($protein_sequence);
        $line->{Protein_Length} = $protein_length;
        $altered_aa_start = $tv->{translation_start};
        $altered_aa_end = $tv->{translation_end};
        $altered_base_start = $tv->{cdna_start};
        $altered_base_end = $tv->{cdna_end};
        $pep_allele_string = $tva->pep_allele_string;
        if (defined($pep_allele_string)) {
            #reference is first
            if ($pep_allele_string =~ m/([^\/]+)\/([^\/]+)/) {
                $amino_acid_reference = $1;
                $amino_acid_variant = $2;
                if ($amino_acid_reference ne $amino_acid_variant)
                {
                    $changes_protein = 1;
                }
            }
        }
        
        #get overlapping features on protein from UniProt
        my @uniprot_feature_mutagen;
        my @uniprot_feature_variant;
        my @uniprot_feature_sites;
        my @uniprot_feature_other;
        if (defined($altered_aa_start) && defined($altered_aa_end)){
            my $uniprot_feature = get_uniprot_feature_from_vw($transcript_id, $altered_aa_start, $altered_aa_end);
            foreach my $key (keys %{$uniprot_feature}){
                if ($key eq 'MUTAGEN'){
                    push(@uniprot_feature_mutagen, join(':', @{$uniprot_feature->{$key}}));
                }
                if ($key eq 'VARIANT'){
                    push(@uniprot_feature_variant, join(':', @{$uniprot_feature->{$key}}));
                }
                elsif(grep {$key eq $_} @{$config->{uniprot_site_features}}){
                    push(@uniprot_feature_sites, "$key:".join(':', @{$uniprot_feature->{$key}}));
                }
                else{
                    push(@uniprot_feature_other, "$key:".join(':', @{$uniprot_feature->{$key}}));
                }
            }
            if ((defined(\@uniprot_feature_mutagen))
                && (scalar(@uniprot_feature_mutagen) > 0))
            {
                $line->{MUTAGEN} = join( ';', @uniprot_feature_mutagen);
            }
            if ((defined(\@uniprot_feature_variant))
                && (scalar(@uniprot_feature_variant) > 0))
            {
                $line->{VARIANT} = join( ';', @uniprot_feature_variant);
            }
            if ((defined(\@uniprot_feature_sites))
                && (scalar(@uniprot_feature_sites) > 0))
            {
                $line->{SITES} = join( ';', @uniprot_feature_sites);
            }
            if ((defined(\@uniprot_feature_other))
                && (scalar(@uniprot_feature_other) > 0))
            {
                $line->{OTHER_OVERLAPPING_FEATURES} = join( ';', @uniprot_feature_other);
            }
        }

        #get a list of aligned orthologous resides from model species
        #will be undef if species and model are same
        #results are used in get_overlapping_protein_features_from_UniProt_orthologues
        #and get_model_phenotypes
        my $aligned_model_orthologue_proteins =
            get_model_orthologous_residues(
                $gene_id,
                $protein_id,
                $altered_aa_start,
                $altered_aa_end
            );

        my ($aligned_homolog_residues, $homolog_species) =
            determine_aligned_homolog_residues(
                $changes_protein,
                $gene_id,
                $protein_id,
                $altered_aa_start,
                $altered_aa_end,
                $amino_acid_reference,
                \%consequences_ranks_so, #$consequence
                $config->{comparison_species} #$comparison_species
                );
        my $alignment_score_change =
            determine_alignment_score_change(
                $config->{matrix},
                $config->{max_alignment_score_change},
                $aligned_homolog_residues,
                $amino_acid_reference,
                $amino_acid_variant
                );
        $line->{Alignment_Score_Change} = $alignment_score_change;

        my $reference_amino_acid_conservation =
            determine_reference_amino_acid_conservation(
                $aligned_homolog_residues,
                $amino_acid_reference
                );
        my $reference_amino_acid_score =
            determine_reference_amino_acid_score(
                $config->{matrix},
                $aligned_homolog_residues,
                $amino_acid_reference
                );
        $line->{C_blosum} = $reference_amino_acid_score;

        my ($context_alignments, $homolog_species_context) =
            determine_context_alignments(
                $changes_protein,
                $altered_aa_start,
                $protein_sequence,
                $gene_id,
                $protein_id
                );
        my $context_conservation = #context_average_percent_identity
            determine_context_average_percent_identity(
                $context_alignments
                );
        $line->{Context_Conservation} = $context_conservation;

        my $overlapping_protein_domains =
            determine_overlapping_protein_domains(
                $translation,
                $altered_aa_start,
                $altered_aa_end
                );
        #protein domains (from Ensembl) overlapping with affected amino acid
        if ((defined( $overlapping_protein_domains))
            && (scalar(@{$overlapping_protein_domains}) > 0))
        {
            s/\;//g for (@{$overlapping_protein_domains});
            $line->{Overlapping_Protein_Domains} =
                join( ';', @{$overlapping_protein_domains});
        }

        # order the orthologues according to genetic distance from my species
        if ((defined($aligned_homolog_residues))
            && (scalar( @{$aligned_homolog_residues}) > 0)
            && (defined($homolog_species))
            && (scalar(@{$homolog_species}) > 0))
        {
            ($aligned_homolog_residues, $homolog_species) =
              &orderOrthologues($aligned_homolog_residues, $homolog_species, $distHash_ptr);
        }
        #Information about protein site conservation
        my $first_long;
        if ((defined($aligned_homolog_residues))
            && (scalar(@{$aligned_homolog_residues}) > 0))
        {
            #separate AAs by semicolons if any are more than one AA
            $first_long = first {length($_)>1} @{$aligned_homolog_residues};
            $line->{Amino_Acids_In_Orthologues} =
                join(defined($first_long) ? ',' : '', @{$aligned_homolog_residues});
         }
        if ((defined($homolog_species))
            && (scalar(@{$homolog_species}) > 0))
        {
            my @species_short;
            for my $sp (@{$homolog_species}){
                push @species_short, exists($config->{species_short}->{$sp}) ?
                    $config->{species_short}->{$sp} : $sp;
            }
            my @foo = map {exists($config->{species_short}->{$_}) ?
                $config->{species_short}->{$_} : $_} @{$homolog_species};
            $line->{Orthologue_Species} =
                join(',', @{$homolog_species});
            $line->{Orthologue_Species} =
                join(defined($first_long) ? ',' : '', map {exists($config->
                    {species_short}->{$_}) ? $config->{species_short}->{$_} : $_}
                     @{$homolog_species});
        }

        my $model_phenotype_annotations; #array ptr
        if (defined($config->{model})) {
            $model_phenotype_annotations =
                get_model_phenotypes(
                    \%genomic_coords,
                    $protein_id,
                    $altered_aa_start,
                    $altered_aa_end,
                    $aligned_model_orthologue_proteins
                    );
        }
        #Phenotypic information associated with known variation at site (or orthologous site) in the model species
        #using information from Ensembl
        if ((defined($model_phenotype_annotations))
            && (scalar(@{$model_phenotype_annotations}) > 0)
            )
        {
            s/[\|\;]//g for (@{$model_phenotype_annotations});
            my $phenotypes_position = join( '|', @{$model_phenotype_annotations});
            $line->{Phenotypes_Position} = $phenotypes_position;
        }

        #determine how STOP_GAINED SNP changes protein length
        #this may need to be altered for indels
        my ($protein_length_decrease, $protein_sequence_lost) =
            determine_effect_of_stop_gained_on_protein(
                \%consequences_ranks_so,
                $protein_sequence,
                $altered_aa_start
                );
        #Protein_Length_Decrease
        if (defined($protein_length_decrease)) {
            $line->{Protein_Length_Decrease} = $protein_length_decrease;
        }
        #Protein_Sequence_Lost
        if (defined($protein_sequence_lost)) {
            $line->{Protein_Sequence_Lost} = $protein_sequence_lost;
        }

        #determine how STOP_LOST SNP changes protein length
        #this may need to be altered for indels
        my ($protein_sequence_gained, $protein_length_increase) =
            determine_effect_of_stop_lost_on_protein(
                \%consequences_ranks_so,
                $protein_sequence,
                $altered_base_start,
                $transcript,
                $amino_acid_variant #$transcript_snp_reads
                );
        #Protein_Length_Increase
        if (defined($protein_length_increase)) {
            $line->{Protein_Length_Increase} = $protein_length_increase;
        }
        #Protein_Sequence_Gained
        if (defined($protein_sequence_gained)) {
            $line->{Protein_Sequence_Gained} = $protein_sequence_gained;
        }

        #determine how ESSENTIAL_SPLICE_SITE SNP changes
        #splice site
        my ($reference_splice_site, $variant_splice_site) =
            determine_effect_of_essential_splice_site_on_splice_site(
                \%consequences_ranks_so,
                $transcript,
                \%genomic_coords,
                $reference_allele, #$Chromosome_Reference
                $non_reference_allele
                );
        #Reference_Splice_Site
        if (defined($reference_splice_site)) {
            $line->{Reference_Splice_Site} = $reference_splice_site;
        }
        #Variant_Splice_Site
        if (defined($variant_splice_site)) {
            $line->{Variant_Splice_Site} = $variant_splice_site;
        }

        my ($correlated_family, $correlated_orthologues);
        if ($config->{run_omes}) {
            ($correlated_family, $correlated_orthologues) =
                determine_if_correlated_amino_acid(
                    $protein_id,
                    $changes_protein,
                    $altered_aa_start,
                    $amino_acid_reference,
                    $gene_id
                    );
        }
        #Protein_OMES_Orthologues
        if (defined($correlated_orthologues)) {
            my $protein_omes_orthologues = $correlated_orthologues;
            $line->{Protein_OMES_Orthologues} = $protein_omes_orthologues;
        }
        #Protein_OMES_Family
        if (defined($correlated_family)) {
            my $protein_omes_family = $correlated_family;
            $line->{Protein_OMES_Family} = $protein_omes_family;
        }
     } #if(defined($translation))
} #my_tva_annotations

sub replace_str{
    #uses zero-based indexing
    my ($string, $replacement, $start, $end) = @_;
    my $a = substr($string,0,$start);
    my $z = substr($string,$end+1);
    my $replaced = $a . $replacement . $z;
    return $a . $replacement . $z;
}
# Perl trim function to remove whitespace from the start and end of the string
sub trim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}
# Left trim function to remove leading whitespace
sub ltrim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	return $string;
}
# Right trim function to remove trailing whitespace
sub rtrim($)
{
	my $string = shift;
	$string =~ s/\s+$//;
	return $string;
}
#MY

#NGS-SNP
sub get_unique_seqs {
    my $seq  = shift;
    my %seen = ();
    my @uniq = grep { !$seen{ $_->id }++ } @{$seq};
    return \@uniq;
}

sub get_unique {
    my $list = shift;
    my %seen = ();
    my @uniq = grep { !$seen{$_}++ } @{$list};
    return \@uniq;
}

#Determine effect of 'ESSENTIAL_SPLICE_SITE' consequences on splice site
sub determine_effect_of_essential_splice_site_on_splice_site {
    
    my ($consequences_ranks_so, $transcript, $genomic_coords, $reference_allele, $non_reference_allele) = @_;
#TODO: indel
    
    #return values, as DD~AA
    my ($reference_splice_site, $variant_splice_site);

    #if ($consequence_ensembl eq 'ESSENTIAL_SPLICE_SITE' ) {
    if (exists $consequences_ranks_so->{splice_acceptor_variant}
        || exists $consequences_ranks_so->{splice_donor_variant}){

        my $ref = $genomic_coords->{strand}==1 ? $reference_allele : complement($reference_allele);
        my $alt = $genomic_coords->{strand}==1 ? $non_reference_allele : complement($non_reference_allele);
        $ref = $ref eq '-' ? '' : $ref;
        $alt = $alt eq '-' ? '' : $alt;

        my $introns = $transcript->get_all_Introns();

        foreach my $intron (@{$introns}) {

            $intron = $intron->transform('toplevel');

            my $intron_strand = $intron->strand();
            my $intron_start  = $intron->start();
            my $intron_end    = $intron->end();
            my $intron_length = $intron->length();
            my $intron_seq    = $intron->seq();
            
            if (($genomic_coords->{start} >= $intron_start && $genomic_coords->{start} <= $intron_end)
                || ($genomic_coords->{end} <= $intron_end && $genomic_coords->{end} <= $intron_end)){
                my $ref_seq = $intron_strand == 1 ? $intron_seq : complement($intron_seq);
                my $rel_start = $genomic_coords->{start}-$intron_start;
                my $rel_end = $genomic_coords->{end}-$intron_start;
                my $alt_seq;
                if ($ref eq '' && $genomic_coords->{start} == $genomic_coords->{end}+1){
                    #ins
                    if ($rel_start >= 0 && $rel_start <= $intron_length){
                        #insertion point adjacent to or in intron
                        $alt_seq = replace_str($ref_seq, $alt, $rel_start, $rel_end);
                    }
                    else{
                        #insertion point is outside intron
                        last;
                    }
                }
                elsif ($alt eq ''){
                    #del
                    if ($rel_start < 0){
                        $rel_start = 0;
                    }
                    if ($rel_end > $intron_length-1){
                        $rel_end = $intron_length-1;
                    }
                    $alt_seq = replace_str($ref_seq, $alt, $rel_start, $rel_end);
                }
                else{
                    $alt_seq = replace_str($ref_seq, $alt, $rel_start, $rel_end);
                }
                if($intron_strand == -1){
                    $alt_seq = complement($alt_seq);
                }
                my $donor = $intron_strand == 1 ?
                    substr($ref_seq,0,2) :
                    reverse substr($ref_seq,-2,2);
                my $acceptor = $intron_strand == 1 ?
                    substr($ref_seq,-2,2) :
                    reverse substr($ref_seq,0,2);
                $reference_splice_site = $donor.'~'.$acceptor;
                my $alt_donor = $intron_strand == 1 ?
                    substr($alt_seq,0,2) :
                    reverse substr($alt_seq,-2,2);
                my $alt_acceptor = $intron_strand == 1 ?
                    substr($alt_seq,-2,2) :
                    reverse substr($alt_seq,0,2);
                $variant_splice_site = $alt_donor.'~'.$alt_acceptor;
                last;
            }
            #
            ##The intron can be on the forward (1) or reverse (-1) strand.
            ##12..34, i.e. GT..AG if intron strand = 1
            ##43..21, i.e. CT..AC if intron strand = -1
            #my $variant_position_in_splice_site = undef;
            #
            #if ($intron_strand == 1) {
            #    if ($genomic_coords->{start} == $intron_start) {
            #        $variant_position_in_splice_site = 1;
            #    }
            #    elsif ($genomic_coords->{start} == ($intron_start + 1))
            #    {
            #        $variant_position_in_splice_site = 2;
            #    }
            #    elsif ($genomic_coords->{start} == ($intron_end - 1)) {
            #        $variant_position_in_splice_site = 3;
            #    }
            #    elsif ($genomic_coords->{start} == $intron_end) {
            #        $variant_position_in_splice_site = 4;
            #    }
            #    else {
            #        next;
            #    }
            #}
            #elsif ($intron_strand == -1) {
            #    if ($genomic_coords->{start} == $intron_start) {
            #        $variant_position_in_splice_site = 4;
            #    }
            #    elsif ($genomic_coords->{start} == ($intron_start + 1))
            #    {
            #        $variant_position_in_splice_site = 3;
            #    }
            #    elsif ($genomic_coords->{start} == ($intron_end - 1)) {
            #        $variant_position_in_splice_site = 2;
            #    }
            #    elsif ($genomic_coords->{start} == $intron_end) {
            #        $variant_position_in_splice_site = 1;
            #    }
            #    else {
            #        next;
            #    }
            #}
            #
            ##The intron sequence is in sense orientation already
            #my $first = get_subsequence($intron_seq, 1, 2);
            #my $second = get_subsequence($intron_seq, ($intron_length - 1),$intron_length);
            #
            #if ((!defined($first)) || (!defined($second))) {
            #    return;
            #}
            #
            #$reference_splice_site = $first . $second;
            #
            ##move $non_reference_allele to forward strand
            ##$non_reference_allele is on the strand given by $Chromosome_Strand
            #my $variant_base_in_splice_site;
            #if ($genomic_coords->{strand} == 1) {
            #    $variant_base_in_splice_site = $non_reference_allele;
            #}
            #elsif ($genomic_coords->{strand} == -1) {
            #    $variant_base_in_splice_site = reverse_complement($non_reference_allele);
            #}
            #else {
            #    die("Unexpected 'Chromosome_Strand' value encountered");
            #}
            #
            ##then move from forward strand to sense strand
            #if ($intron_strand == -1) {
            #    $variant_base_in_splice_site
            #        = reverse_complement($variant_base_in_splice_site);
            #}
            #
            #if (length($reference_allele) == 1 && $reference_allele ne '-'
            #    && length($non_reference_allele) == 1 && $non_reference_allele ne '-'){
            #    #SNV
            #    #check to confirm that SNP reference base matches reference base in splice site
            #    my $reference_base_from_splice_site
            #        = get_base($variant_position_in_splice_site,
            #        $reference_splice_site );
            #
            #    my $reference_base_from_snp_forward_strand;
            #    if ($genomic_coords->{strand} == 1) {
            #        $reference_base_from_snp_forward_strand
            #            = $reference_allele;
            #    }
            #    elsif ($genomic_coords->{strand} == -1) {
            #        $reference_base_from_snp_forward_strand
            #            = reverse_complement( $reference_allele);
            #    }
            #    else {
            #        die("Unexpected 'Chromosome_Strand' value encountered");
            #    }
            #
            #    my $reference_base_from_splice_site_forward_strand;
            #    if ($intron_strand == -1) {
            #        $reference_base_from_splice_site_forward_strand
            #            = reverse_complement($reference_base_from_splice_site);
            #    }
            #    else {
            #        $reference_base_from_splice_site_forward_strand
            #            = $reference_base_from_splice_site;
            #    }
            #    
            #    if (uc($reference_base_from_splice_site_forward_strand) ne
            #        uc($reference_base_from_snp_forward_strand))
            #    {
            #        die("Unexpected bases encountered");
            #    }
            #
            #    $variant_splice_site
            #        = replace_base($variant_position_in_splice_site,
            #        $reference_splice_site, $variant_base_in_splice_site);
            #}
            #else{
            #    #indel
            #    my $seq = $intron_seq;
            #}
             
        } #foreach my $intron (@{$introns})
    } #if (exists $consequences_ranks_so->{splice_acceptor_variant} || exists $consequences_ranks_so->{splice_donor_variant})

    return ($reference_splice_site, $variant_splice_site);
}
sub get_subsequence {
    my $sequence = shift;
    my $start    = shift;
    my $end      = shift;

    if ( ( $start < 1 ) || ( $start > $end ) || ( $end > length($sequence) ) )
    {

        #this happens rarely
        return undef;
    }

    return substr( $sequence, $start - 1, ( $end - $start + 1 ) );
}

#Determine effect of 'STOP_LOST' consequences on protein
sub determine_effect_of_stop_lost_on_protein {

    my ($consequences_ranks_so, $protein_sequence, $altered_base_start, $transcript, $transcript_snp_reads) = @_;

    #return values
    my ($protein_sequence_gained, $protein_length_increase);
    
    #NOTE: $codon_table should be set to match codon table used to
    #generate reference protein.
    my $codon_table = 1;

    #if ($consequence_ensembl eq 'STOP_LOST') {
    if (exists $consequences_ranks_so->{stop_lost}){
        if (defined($protein_sequence)) {
            my $reference_protein_length
                = length($protein_sequence);
            if ($protein_sequence =~ m/\*$/) {
                $reference_protein_length--;
            }

            my $transcript_sequence = $transcript->seq()->seq();
            if ((defined($transcript_sequence))
                && (defined($altered_base_start)))
            {

                my $translation_start = $transcript->cdna_coding_start;
#TODO: indel?
                my $variant_transcript_sequence = replace_base(
                    $altered_base_start,
                    $transcript_sequence,
                    $transcript_snp_reads
                );

                #get the variant sequence starting with the start codon
                my $variant_transcript_starting_with_start_codon
                    = substr($variant_transcript_sequence, $translation_start - 1, 1)
                    . substr($variant_transcript_sequence, $translation_start,
                    length($variant_transcript_sequence) - $translation_start);

                #create a sequence object to represent variant sequence and translate
                my $variant_transcript_starting_with_start_codon_seq_obj
                    = Bio::Seq->new(
                    -seq => $variant_transcript_starting_with_start_codon,
                    -alphabet => 'dna'
                    );
                my $variant_transcript_starting_with_start_codon_prot_obj
                    = $variant_transcript_starting_with_start_codon_seq_obj
                    ->translate(undef, undef, undef, $codon_table);
                my $variant_transcript_starting_with_start_codon_protein_sequence
                    = $variant_transcript_starting_with_start_codon_prot_obj
                    ->seq();

                #get the reference sequence starting with the start codon
                my $reference_transcript_starting_with_start_codon
                    = substr($transcript_sequence, $translation_start - 1, 1)
                    . substr( $transcript_sequence, $translation_start,
                    length($transcript_sequence) - $translation_start);

                #create a sequence object to represent reference sequence and translate
                my $reference_transcript_starting_with_start_codon_seq_obj
                    = Bio::Seq->new(
                    -seq => $reference_transcript_starting_with_start_codon,
                    -alphabet => 'dna'
                    );
                my $reference_transcript_starting_with_start_codon_prot_obj
                    = $reference_transcript_starting_with_start_codon_seq_obj
                    ->translate(undef, undef, undef, $codon_table);
                my $reference_transcript_starting_with_start_codon_protein_sequence
                    = $reference_transcript_starting_with_start_codon_prot_obj
                    ->seq();

                #need to examine translations to identify added region
                if ($reference_transcript_starting_with_start_codon_protein_sequence
                    =~ m/^([^\*]+)/)
                {
                    my $reference_cds_translation = $1;

                    if (!(length($reference_cds_translation)
                            == $reference_protein_length))
                    {
                        #something is wrong with translations--don't finish calculation
                        return;
                    }

                    if ($variant_transcript_starting_with_start_codon_protein_sequence
                        =~ m/\Q$reference_cds_translation\E([^\*]*\*?)/)
                    {
                        $protein_sequence_gained = $1;
                        $protein_length_increase = length($protein_sequence_gained);
                        if ($protein_sequence_gained =~ m/\*$/ )
                        {
                            $protein_length_increase--;
                        }

                        my $percentage_length_change = sprintf("%.0f",
                            ($protein_length_increase / $reference_protein_length) * 100);
                        $protein_length_increase =
                            $protein_length_increase . "($percentage_length_change)";
                    }
                }
            }
        }
    }
    return ($protein_sequence_gained, $protein_length_increase);
}

#Determine effect of 'STOP_GAINED' consequence on protein
sub determine_effect_of_stop_gained_on_protein {

    my ($consequences_ranks_so, $protein_sequence, $altered_aa_start) = @_;
    
    #return values
    my ($protein_length_decrease, $protein_sequence_lost);

    #if ($consequence_ensembl eq 'STOP_GAINED') {
    if (exists $consequences_ranks_so->{stop_gained}){
        if (defined($protein_sequence)) {
            my $reference_protein_length
                = length($protein_sequence);
            if ( $protein_sequence =~ m/\*$/) {
                $reference_protein_length--;
            }
            if (defined($altered_aa_start)) {
                my $length_lost = $reference_protein_length
                    - $altered_aa_start + 1;
                my $percentage_length_change = sprintf( "%.0f",
                    ($length_lost / $reference_protein_length) * 100);
                $protein_length_decrease
                    = $length_lost . "($percentage_length_change)";

                $protein_sequence_lost = substr($protein_sequence,
                    $altered_aa_start - 1, 1)
                    . substr($protein_sequence, $altered_aa_start,
                    length($protein_sequence) - $altered_aa_start);
            }
        }
    }
    return ($protein_length_decrease, $protein_sequence_lost);
}

#Get overlapping protein features orthologous proteins, from UniProt
#undef if model and species are same
sub get_overlapping_protein_features_from_UniProt_orthologues {
    my $aligned_model_orthologue_proteins = shift;

    my @uniprot_overlapping_protein_features_orthologues;

    foreach my $orthologous_protein(@{$aligned_model_orthologue_proteins})
    {

        if ((!defined( $orthologous_protein->{uniprot_id}))
            || (!defined( $orthologous_protein->{sequence}))
            || (!defined( $orthologous_protein->{start}))
            || (!defined( $orthologous_protein->{end})))
        {
            next;
        }

        
        foreach my $id (@{$orthologous_protein->{uniprot_id}}) {

            my $record = get_uniprot_record(
                script => $config->{uniprot_script},
                id     => $id,
                db     => 'uniprotkb',
                format => 'SWISS',
                style  => 'raw'
            );

            my $overlapping_features
                = get_overlapping_features_from_uniprot_record(
                record            => $record,
                expected_sequence => $orthologous_protein->{sequence},
                region_start      => $orthologous_protein->{start},
                region_end        => $orthologous_protein->{end},
                id                => $id
                );

            if ( defined($overlapping_features) ) {
                push(@uniprot_overlapping_protein_features_orthologues,
                    @{$overlapping_features}
                );
            }
        }
    }
    @uniprot_overlapping_protein_features_orthologues
        = @{get_unique(\@uniprot_overlapping_protein_features_orthologues)};
    return \@uniprot_overlapping_protein_features_orthologues;
}

#Get overlapping protein features from UniProt
sub get_overlapping_protein_features_from_UniProt {
    my ($uniprot_id, $protein_sequence, $altered_aa_start, $altered_aa_end) = @_;

    if ((!defined($uniprot_id))
        || (!defined($protein_sequence))
        || (!defined($altered_aa_start) )
        || (!defined($altered_aa_end)))
    {
        return;
    }

    my @uniprot_overlapping_protein_features;

    foreach my $id (@{$uniprot_id}) {

        my $record = get_uniprot_record(
            script => $config->{uniprot_script},
            id     => $id,
            db     => 'uniprotkb',
            format => 'SWISS',
            style  => 'raw'
        );

        my $overlapping_features
            = get_overlapping_features_from_uniprot_record(
            record            => $record,
            expected_sequence => $protein_sequence,
            region_start      => $altered_aa_start,
            region_end        => $altered_aa_end,
            id                => $id
            );
        if (defined($overlapping_features)) {
            push(
                @uniprot_overlapping_protein_features,
                @{$overlapping_features}
            );
        }
    }
    @uniprot_overlapping_protein_features
        = @{get_unique(\@uniprot_overlapping_protein_features)};
    return \@uniprot_overlapping_protein_features;
}

sub get_overlapping_features_from_uniprot_record {
    my %args = (@_);

    if (!defined( $args{record})) {
        return undef;
    }
#
    local $/ = "\n//\n";
    my $fullParse=1;
    # Read the entry
    my  $entry = SWISS::Entry->fromText($args{record}, $fullParse);
    my $seq = $entry->SQ;
    if (!defined($seq)) {

        #sometimes the accession obtained from Ensembl
        #doesn't work for retrieval, e.g. ANK36_HUMAN
        return undef;
    }
    if (uc($seq) ne uc($args{expected_sequence})) {

        #sometimes the sequences don't match
        #do not transfer features if this is the case
        return undef;
    }
    my @FTs = $entry->FTs->elements();
    my @features = ();
    for my $ft(@FTs){
        my %feature = (
            key => @{$ft}[0],
            start =>  @{$ft}[1],
            end => @{$ft}[2],
            description => @{$ft}[3]
        );
        if (overlaps(
                s1 => $args{region_start},
                e1 => $args{region_end},
                s2 => $feature{start},
                e2 => $feature{end})
            ){
            #push(@features, "$key:$start:$end:$description");
            push(@features, \%feature);
        }
    }

    if (scalar(@features) > 0) {
        return \@features;
    }
    
    return undef;
}

#Get UniProt record from EBI
sub get_uniprot_record {
    my %args = (@_);

    #create temp file for output
    my $tmp_output = new File::Temp();
    my $tmp_output_filename = $tmp_output->filename;

    my $command
        = 'perl '
        . $args{script}
        . " -i '$args{id}'"
        . " -o $tmp_output_filename"
        . " -d '$args{db}'"
        . " -f '$args{format}'"
        . " -s '$args{style}'";

    my $result = system($command);
    if ( $result != 0 ) {
        die("The following command failed: '$command'\n");
    }

    close($tmp_output) or die("Cannot close file : $!");

    local $/ = undef;
    open( my $FILE, '<', $tmp_output_filename )
        or die("Cannot open file '$tmp_output_filename': $!");

    my $output = <$FILE>;

    close($FILE) or die("Cannot close file : $!");
    return $output;
}

#Get information about model orthologue from NCBI Gene database.
#Rationale is that model genes have much more detailed annotation.
sub get_model_orthologue_info_from_NCBI {
    my ($entrez_gene_id, $gene_id) = @_;

    #both of the following can be used to retrieve Entrez Gene records
    my $model_entrez_gene_id;
    my $model_refseq_peptide_accession;

    if ( lc($config->{species}) eq lc($config->{model})) {
        if (defined( $entrez_gene_id)) {
            $model_entrez_gene_id = $entrez_gene_id;
        }
        else {
            return;
        }
    }
    else {
        if (!defined($gene_id)) {
            return;
        }

        #$options->{ma} is a Bio::EnsEMBL::Compara::DBSQL::MemberAdaptor
        #$member is a Bio::EnsEMBL::Compara::Member
        my $member = $config->{ma}->fetch_by_source_stable_id('ENSEMBLGENE', $gene_id);

        #Rarely the $member object may be undef
        if (!defined($member)) {
            return;
        }

        #$options->{ha} is a Bio::EnsEMBL::Compara::DBSQL::HomologyAdaptor
        #$homologies is a list of Bio::EnsEMBL::Compara::Homology objects

        my $homologies = [];
        push(
            @{$homologies},
            @{$config->{ha}->fetch_all_by_Member_paired_species($member,
                    $config->{model}, ['ENSEMBL_ORTHOLOGUES'])}
        );

        foreach my $homology (@{$homologies}) {

     #$homologues is an array ref of (2) Bio::EnsEMBL::Compara::Member objects

            my $homologues     = undef;
            my $homologue_gene = undef;

            eval {
                $homologues     = $homology->gene_list();
                $homologue_gene = $$homologues[1]->get_Gene;
            };
            if ($@) {
                next;
            }

            my $xrefs = $homologue_gene->get_all_DBLinks();

            #try to obtain Entrez Gene ID first
            foreach my $xref (@{$xrefs}) {
                if ($xref->dbname() eq 'EntrezGene') {
                    $model_entrez_gene_id = $xref->primary_id();
                }
            }
            if (defined($model_entrez_gene_id)) {
                last;
            }

            #try to obtain RefSeq protein ID second
            foreach my $xref (@{$xrefs}) {
                if ($xref->dbname() eq 'RefSeq_peptide') {
                    $model_refseq_peptide_accession = $xref->primary_id();
                }
            }
            if (defined($model_refseq_peptide_accession)) {
                last;
            }

        }

    }

    my $xml;
    if (defined($model_entrez_gene_id)) {
        $xml = get_ncbi_record(
            script => $config->{ncbi_script},
            db     => 'gene',
            query  => $model_entrez_gene_id . '[UID]'
        );
    }
    elsif (defined($model_refseq_peptide_accession)) {
        $xml = get_ncbi_record(
            script => $config->{ncbi_script},
            db     => 'gene',
            query  => $model_refseq_peptide_accession . '[Protein Accession]'
        );
    }

    my $parsed = get_info_from_ncbi_gene_record($xml);

    my $entrez_gene_model_kegg_pathways = $parsed->{kegg};
    my $entrez_gene_model_phenotypes = $parsed->{phenotypes};
    my $entrez_gene_model_interactions_count = $parsed->{interaction_count};
    return ($entrez_gene_model_kegg_pathways, $entrez_gene_model_phenotypes, $entrez_gene_model_interactions_count);
}

#this parsing should be improved
sub get_info_from_ncbi_gene_record {
    my $xml = shift;

    my %results
        = ( kegg => undef, phenotypes => undef, interaction_count => undef );

    if ( !defined($xml) ) {
        return \%results;
    }

    #KEGG pathways
    #<Gene-commentary_text>KEGG pathway: Apoptosis</Gene-commentary_text>
    my $kegg_re
        = "\Q<Gene-commentary_text>KEGG pathway: \E(.+)\Q</Gene-commentary_text>\E";

    $results{kegg} = get_unique_re_matches( $kegg_re, $xml );

#phenotypes
#<Gene-commentary_type value="phenotype">19</Gene-commentary_type>
#<Gene-commentary_heading>Congenital bilateral absence of vas deferens</Gene-commentary_heading>
    my $phenotypes_re
        = "\Q<Gene-commentary_type value=\"phenotype\">19</Gene-commentary_type>\E"
        . '[\s\n]*'
        . "\Q<Gene-commentary_heading>\E(.+)\Q</Gene-commentary_heading>\E";

    $results{phenotypes} = get_unique_re_matches( $phenotypes_re, $xml );

    #interactions
    #<Gene-commentary_heading>Interactions</Gene-commentary_heading>
    #to
    #<Gene-commentary_heading>Pathways</Gene-commentary_heading>
    my $interactions_section;
    if ( $xml
        =~ m/<Gene\-commentary_heading>Interactions<\/Gene\-commentary_heading>([\s\S]+)<Gene\-commentary_heading>Pathways<\/Gene\-commentary_heading>/
        )
    {
        $interactions_section = $1;
    }
    if ( defined($interactions_section) ) {
        my $interactions_re
            = "\Q<Gene-commentary_type value=\"generif\">18</Gene-commentary_type>\E"
            . '[\s\n]*'
            . "\Q<Gene-commentary_text>\E(.+)\Q</Gene-commentary_text>\E";

        my $interactions = get_unique_re_matches( $interactions_re,
            $interactions_section );
        if ( defined($interactions) ) {
            $results{interaction_count} = scalar( @{$interactions} );
        }

    }
    return \%results;
}

sub get_unique_re_matches {
    my $re      = shift;
    my $text    = shift;
    my @matches = ();
    while ( $text =~ m/$re/g ) {
        push( @matches, $1 );
    }

    if ( ( scalar(@matches) ) == 0 ) {
        return undef;
    }
    return get_unique( \@matches );
}

sub get_ncbi_record {
    my %args = (@_);

    my $tmp_output          = new File::Temp();
    my $tmp_output_filename = $tmp_output->filename;

    my $command
        = 'perl '
        . $args{script}
        . " -q '$args{query}'"
        . " -o $tmp_output_filename"
        . " -d '$args{db}'"
        . " -r xml" . " -m 1";

    my $result = system($command);
    if ( $result != 0 ) {
        die("The following command failed: '$command'\n");
    }

    close($tmp_output) or die("Cannot close file : $!");

    local $/ = undef;
    open( my $FILE, '<', $tmp_output_filename )
        or die("Cannot open file '$tmp_output_filename': $!");

    my $output = <$FILE>;

    close($FILE) or die("Cannot close file : $!");
    return $output;
}

#whether the affected amino acid position is correlated with other amino acid
#positions, according to the OMES procedure described in
#Fodor and Aldrich PROTEINS: Structure, Function, and Bioinformatics 56:211–221
#and
#Andreas Kowarsch et al, PLOS Computational Biology, 6:e1000923
sub determine_if_correlated_amino_acid {
    
    my ($protein_id, $changes_protein, $altered_aa_start, $amino_acid_reference, $gene_id) = @_;
    
    #return values
    my ($correlated_family, $correlated_orthologues);

    if ((!defined($protein_id))
        || (!defined($changes_protein))
        || (!defined( $config->{omes_script})))
    {
        return;
    }

    #A multiple protein alignment is needed for the calculation used to determine
    #correlated mutations

    if ($config->{omes_family}) {
        #######################
        #Option 1: Family object
        #One option is to obtain a Family object.
        #Families are clusters of proteins including all the EnsEMBL proteins plus all the metazoan SwissProt and SP-Trembl entries
        #The drawback of this option is that it won't necessarily include only orthologues. The paper by Kowarsch et al only
        #used orthologues.
        #The inclusion of all family members could prevent a correlated site specific to the orthologues from being detected.
        #For example, paralogous members could have a slightly different catalytic site that leads to different correlated residues.
        my @seqs_for_omes = ();

        my $results = get_aligned_proteins_for_omes_from_family(
            $protein_id,
            $altered_aa_start,
            $amino_acid_reference,
            $gene_id);

        if ((defined($results))
            && (defined($results->{aligned_seqs})))
        {
            my $is_first = 1;
            foreach my $aligned_seq (@{$results->{aligned_seqs}}) {
                if ($is_first) {
                    if ($aligned_seq->id ne $results->{altered_protein_id})
                    {
                        #can occur if altered protein not in group returned
                        #using Ensembl Gene ID
                        return;
                    }
                }

                #remove gaps since OMES script will redo the alignments
                my $sequence = $aligned_seq->seq();
                $sequence =~ s/\-//g;

                push(@seqs_for_omes,
                    '>' . $aligned_seq->id . "\n" . $sequence . "\n");
                $is_first = 0;
            }

            if (scalar(@seqs_for_omes) > 0) {

                if ((defined($config->{omes_max_seqs}))
                    && (scalar(@seqs_for_omes) > $config->{omes_max_seqs}))
                {
                    @seqs_for_omes =
                        @seqs_for_omes[0 .. ($config->{omes_max_seqs} - 1)];
                }

                my $omes_result = calculate_omes(
                    input_seqs           => \@seqs_for_omes,
                    min_seqs             => $config->{omes_minimum_seqs},
                    max_ident            => $config->{omes_maximum_identity},
                    position_of_interest => $altered_aa_start,
                    omes_script          => $config->{omes_script},
                    protein_id           => $protein_id,
                    comparison_type      => 'family'
                );

                $correlated_family = $omes_result;
            }
        }
    }

    #######################
    if ($config->{omes_orthologues}) {

        #######################
        #Option 2: Use a list of Bio::EnsEMBL::Compara::Homology objects.
        #Homology objects store orthologous and paralogous relationships between Members
        #This approach can be used to obtain orthologues and their sequences.
        #Drawback is that the sequences will need to be aligned
        #######################
        my @seqs_for_omes = ();

        my $results = get_protein_orthologues_for_omes_from_homology(
            $protein_id,
            $altered_aa_start,
            $amino_acid_reference,
            $gene_id);

        if ((defined($results)) && (defined($results->{seqs}))) {
            my $is_first = 1;
            foreach my $seq (@{$results->{seqs}}) {
                if ($is_first) {
                    if ($seq->id ne $results->{altered_protein_id}) {
                        #can occur if altered protein not in group returned
                        #using Ensembl Gene ID
                        return;
                    }
                }

                #remove gaps since they were not added by multiple alignment
                my $sequence = $seq->seq();
                $sequence =~ s/\-//g;

                push(@seqs_for_omes,
                    '>' . $seq->id . "\n" . $sequence . "\n");
                $is_first = 0;
            }

            if (scalar(@seqs_for_omes) > 0) {

                if ((defined($config->{omes_max_seqs}))
                    && (scalar(@seqs_for_omes) > $config->{omes_max_seqs})
                    )
                {
                    @seqs_for_omes =
                        @seqs_for_omes[0 .. ($config->{omes_max_seqs} - 1)];
                }

                my $omes_result = calculate_omes(
                    input_seqs           => \@seqs_for_omes,
                    min_seqs             => $config->{omes_minimum_seqs},
                    max_ident            => $config->{omes_maximum_identity},
                    position_of_interest => $altered_aa_start,
                    omes_script          => $config->{omes_script},
                    protein_id           => $protein_id,
                    comparison_type      => 'orthologues'
                );

                $correlated_orthologues = $omes_result;
            }
        }
    }
    return ($correlated_family, $correlated_orthologues);
}

#Runs OMES_score script and parses output
sub calculate_omes {
    my %args = (@_);

    my $reference_seq = $args{input_seqs}->[0];
    $reference_seq =~ s/>[^\n]+\n//;
    $reference_seq =~ s/\s//g;

    my $reference_length      = length($reference_seq);
    my $top_positions_to_keep = 5;

    #write sequences to temporary file
    my $tmp_input          = new File::Temp();
    my $tmp_input_filename = $tmp_input->filename;

    foreach my $sequence ( @{ $args{input_seqs} } ) {
        print $tmp_input $sequence;
    }
    close($tmp_input) or die("Cannot close file : $!");

    #create temp file for OMES_score.pl score output
    my $tmp_output          = new File::Temp();
    my $tmp_output_filename = $tmp_output->filename;

    #create temp file for OMES_score.pl alignment output
    my $tmp_output_alignment          = new File::Temp();
    my $tmp_output_alignment_filename = $tmp_output_alignment->filename;

    #run OMES_score.pl
    my $omes_command
        = 'perl '
        . $args{omes_script}
        . " -i $tmp_input_filename"
        . " -o $tmp_output_filename"
        . " -a $tmp_output_alignment_filename"
        . " -n $top_positions_to_keep"
        . " -p $args{max_ident}"
        . " -m $args{min_seqs}"
        . " -residue $args{position_of_interest}";

    my $result = system($omes_command);
    if ( $result != 0 ) {
        die("The following command failed: '$omes_command'\n");
    }

    close($tmp_output)           or die("Cannot close file : $!");
    close($tmp_output_alignment) or die("Cannot close file : $!");

    open( my $OMES_FILE, '<', $tmp_output_filename )
        or die("Cannot open file '$tmp_output_filename': $!");

    #sample output is:
    #5.305,32(12),37(13)
    #where first number is OMES score, and positions in parentheses are residues
    #in reference sequence.
    my @omes_scores = ();
    while ( my $line = <$OMES_FILE> ) {

      #a message starting with '#' indicates insufficient data for calculation
      #or some other issue preventing data return
        if ( $line =~ m/^\s*\#/ ) {
            return undef;
        }
        my @values = split( /,/, $line );

        #should always give three values
        if ( scalar(@values) == 3 ) {
            my $score   = $values[0];
            my $column1 = $values[1];
            my $column2 = $values[2];
            my $residue1;
            my $residue2;

            #column values can contain '-' if a residue in
            #the first sequence (reference sequence) is not
            #in the column
            if ( $column1 =~ m/\((.+?)\)/ ) {
                $residue1 = $1;
            }
            if ( $column2 =~ m/\((.+?)\)/ ) {
                $residue2 = $1;
            }

            if ((      ( defined($residue1) )
                    && ( $residue1 eq $args{position_of_interest} )
                )
                || (   ( defined($residue2) )
                    && ( $residue2 eq $args{position_of_interest} ) )
                )
            {
                push( @omes_scores, $score . "($residue1)($residue2)" );
            }
        }
    }

    close($OMES_FILE) or die("Cannot close file : $!");

    if ( scalar(@omes_scores) > 0 ) {
        return join( ',', @omes_scores );
    }

    return undef;
}

#Returns list of Bio::Seq objects describing protein sequence of orthologues to reference protein
sub get_protein_orthologues_for_omes_from_homology {
    
    my ($protein_id, $altered_aa_start, $amino_acid_reference, $gene_id) = @_;
    
    my %results = (
        seqs                            => [],
        altered_protein_id              => undef,
        altered_aa_position_in_sequence => undef,
        altered_aa                      => undef
    );

    $results{altered_protein_id} = $protein_id;
    $results{altered_aa_position_in_sequence} = $altered_aa_start;
    $results{altered_aa} = $amino_acid_reference;

    my $member = $config->{ma}
        ->fetch_by_source_stable_id('ENSEMBLGENE', $gene_id);

    #Rarely the $member object may be undef
    if (!defined($member)) {
        return;
    }

    my $homologies
        = $config->{ha}->fetch_all_by_Member_method_link_type($member,
        'ENSEMBL_ORTHOLOGUES');

    foreach my $homology (@{$homologies}) {

     #$homologues is an array ref of (2) Bio::EnsEMBL::Compara::Member objects
     #$align is Bio::SimpleAlign object

        my $align = undef;

        eval { $align = $homology->get_SimpleAlign(); };
        if ($@) {
            return;
        }

        if (!(defined($align))) {
            return;
        }

        my @aligned_seqs = $align->each_seq();

        #expect two sequences to be added, one of which is reference
        push(@{$results{seqs}}, @aligned_seqs);
    }

    #$results{seqs} will contain many copies of reference--remove duplicates based on id
    $results{seqs} = get_unique_seqs($results{seqs});

    #sort so that reference sequence is first
    $results{seqs} = sort_seqs_so_reference_first($protein_id, $results{seqs} );

    return \%results;
}

#Returns list of Bio::Seq objects describing protein sequences from same sequence family
#aligned to reference protein sequence
sub get_aligned_proteins_for_omes_from_family {
    
    my ($protein_id, $altered_aa_start, $amino_acid_reference, $gene_id) = @_;
    
    my %results = (
        aligned_seqs                     => undef,
        altered_protein_id               => undef,
        altered_aa_position_in_sequence  => undef,
        altered_aa_position_in_alignment => undef,
        altered_aa                       => undef
    );

    $results{altered_protein_id} = $protein_id;
    $results{altered_aa_position_in_sequence} = $altered_aa_start;
    $results{altered_aa} = $amino_acid_reference;

    my $member = $config->{ma}
        ->fetch_by_source_stable_id('ENSEMBLGENE', $gene_id);

    #Rarely the $member object may be undef
    if ( !defined($member) ) {
        return;
    }

    my $families = $config->{fa}->fetch_all_by_Member($member);

    #there can be multiple families
    #will use first family only
    my $align;
    if ( scalar(@{$families}) > 0) {
        my $family = $families->[0];
        $align = $family->get_SimpleAlign();
    }

    if (defined($align)) {

    #determine position of affected residue in alignment
    #get the column containing the residue encoded by the SNP-containing codon
        my $col = undef;

     #column_from_residue_number can throw an exception if the altered residue
     #in the reference is the one that encodes the stop codon
        eval {
            $col
                = $align->column_from_residue_number(
                $results{altered_protein_id},
                $results{altered_aa_position_in_sequence});};
        if ($@) {
            return;
        }

        $results{altered_aa_position_in_alignment} = $col;

        my @aligned_seqs = $align->each_seq();

        #sort so that reference sequence is first
        @aligned_seqs = @{sort_seqs_so_reference_first($protein_id,\@aligned_seqs)};

        $results{aligned_seqs} = \@aligned_seqs;
    }

    return \%results;

}

sub sort_seqs_so_reference_first {
    my $reference_id  = shift;
    my $array_of_seqs = shift;

    my @sorted = map { $_->[0] }

        sort { $b->[1] <=> $a->[1] }
        map {
        [   $_,
            get_sort_seqs_so_reference_first_value( $_->id, $reference_id )
        ]
        } @{$array_of_seqs};

    return \@sorted;

}

sub get_sort_seqs_so_reference_first_value {
    my $sequence_id  = shift;
    my $reference_id = shift;

    if ( $sequence_id eq $reference_id ) {
        return 1;
    }
    return -1;
}

#adds information about known model phenotypes to consequence
#undef if model and species are the same
sub get_model_phenotypes {
    my ($genomic_coords, $protein_id, $altered_aa_start, $altered_aa_end, $aligned_model_orthologue_proteins) = @_;

    my @model_phenotype_annotations;
    my @aligned_phenotypes = ();

    if ((!defined($config->{model_translation_adaptor}))
        || (!defined( $config->{model_transcript_adaptor}))
        || (!defined( $config->{model_variation_adaptor}))
        || (!defined( $config->{model_variationfeature_adaptor}))
        || (!defined( $config->{model_slice_adaptor}))
        || (!defined( $config->{model_variationannotation_adaptor})))
    {
        return;
    }

    #if input SNPs are from model then look for known variants
    #and phenotypes at that position
    if ( lc( $config->{species} ) eq lc( $config->{model} ) ) {
        push(
            @aligned_phenotypes,
            @{get_model_phenotypes_from_model_genomic_coords($genomic_coords)}
        );

    }

    #otherwise transform non-model coordinate to model
    else {
        if ((!defined($protein_id))
            || (!defined($altered_aa_start))
            || (!defined($altered_aa_end)))
        {
            return;
        }

        foreach my $model_orthologous_residue (
            @{$aligned_model_orthologue_proteins})
        {

            my $genomic_coords
                = get_model_genomic_coords_of_model_amino_acid($model_orthologous_residue);

            if (defined($genomic_coords)) {
                push(
                    @aligned_phenotypes,
                    @{get_model_phenotypes_from_model_genomic_coords($genomic_coords)}
                );
            }
        }
    }

    my $unique_phenotypes = get_unique( \@aligned_phenotypes );
    push(
        @model_phenotype_annotations,
        @{$unique_phenotypes}
    );
    return \@model_phenotype_annotations;
}

#returns array of hashes describing model orthologous residues
#apparently this doesn't do anything if species and model are the same
sub get_model_orthologous_residues {
    my ($gene_id, $protein_id, $altered_aa_start, $altered_aa_end) = @_;
    if ((!defined($gene_id))
        || (!defined($protein_id))
        || (!defined($altered_aa_start))
        || (!defined($altered_aa_end))
        || (!defined( $config->{model_translation_adaptor})))
    {
        return;
    }
 
    #return value
    my @model_orthologous_residues = ();

    #$config->{ma} is a Bio::EnsEMBL::Compara::DBSQL::MemberAdaptor
    #$member is a Bio::EnsEMBL::Compara::Member
    my $member = $config->{ma}
        ->fetch_by_source_stable_id('ENSEMBLGENE', $gene_id);

    #Rarely the $member object may be undef
    if (!defined($member)) {
        return;
    }

    #$config->{ha} is a Bio::EnsEMBL::Compara::DBSQL::HomologyAdaptor
    #$homologies is a list of Bio::EnsEMBL::Compara::Homology objects

    my $homologies = [];

    #fetch the homology relationships where the given member is implicated
    #in pair with another member from the paired species. Member species and
    #paired species should be different
    push(
        @{$homologies},
        @{$config->{ha}->fetch_all_by_Member_paired_species($member,
                $config->{model}, ['ENSEMBL_ORTHOLOGUES'])}
    );

    foreach my $homology (@{$homologies}) {

     #$homologues is an array ref of (2) Bio::EnsEMBL::Compara::Member objects
     #$aln is Bio::SimpleAlign object

        #added eval block on 2009-10-22 because
        #$aln = $homology->get_SimpleAlign();
        #throws exception for ENSBTAT00000060548
        my $homologues = undef;
        my $taxon1     = undef;
        my $taxon2     = undef;
        my $aln        = undef;

        eval {
            $homologues = $homology->gene_list();
            $taxon1     = $$homologues[0]->taxon;
            $taxon2     = $$homologues[1]->taxon;
            $aln        = $homology->get_SimpleAlign();
        };
        if ($@) {
            next;
        }

        if ( !( defined($aln) ) ) {
            next;
        }

        #confirm that reference protein is in the alignment
        if (!(scalar($aln->each_seq_with_id($protein_id)) == 1))
        {
            next;
        }

#TODO:verify indel
        #get the column containing the residue encoded by the SNP-containing codon
        my $col = undef;
        my $col_end = undef;

        #This function gives the position in the alignment
        #(i.e. column number) of the given residue number in the
        #sequence with the given name.
        #column_from_residue_number can throw an exception if the altered residue
        #in the reference is the one that encodes the stop codon
        eval {
            $col = $aln->column_from_residue_number(
                $protein_id,
                $altered_aa_start
            );
        };
        if ($@) {
            next;
        }
        eval {
            $col_end = $aln->column_from_residue_number(
                $protein_id,
                $altered_aa_end
            );
        };
        if ($@) {
            next;
        }

        #Bio::SimpleAlign->slice
        #Creates a slice from the alignment inclusive of start and
        #end columns.  Sequences with no residues in the slice are
        #excluded from the new alignment and a warning is printed.
        #Slice beyond the length of the sequence does not do
        #padding.
        #get an alignment containing only the desired column
        #my $sub_align = $aln->slice( $col, $col );
        my $sub_align = undef;
        #prevent annoying warnings; seems to require Bioperl-live, not Bioperl-1.2.3
        if (!defined($config->{verbose})){
            $SIG{'__WARN__'} = sub {return 1};
        }
        eval { $sub_align = $aln->slice($col, $col_end); };
        if ($@) {
            next;
        }

        #confirm that there are two sequences in the subalignment
        #(sequences are excluded if the alignment slice contains no residues from the sequence)
        if (!( scalar( $sub_align->each_seq()) == 2)) {
            next;
        }

        my $seq1 = $sub_align->get_seq_by_pos(1);
        my $seq2 = $sub_align->get_seq_by_pos(2);

        #in all the examples I've examined, $seq1 is the reference sequence
        if ($seq1->display_id eq $protein_id) {

            #want to get Entrez Gene ID from gene object
            my $model_gene = $$homologues[1]->get_Gene;
            my $entrez_gene_id;
            my $refseq_protein_id;
            foreach my $xref (@{$model_gene->get_all_DBLinks()}) {
                if ($xref->dbname() eq 'EntrezGene') {
                    $entrez_gene_id = $xref->primary_id();
                }
                elsif ($xref->dbname() eq 'RefSeq_peptide') {
                    $refseq_protein_id = $xref->display_id();
                }
            }

            #obtain the the Ensembl ID of the model protein
            my $model_translation_ID = $seq2->display_id;

            #obtain the position in the model protein of the amino acid that aligns with the SNP-altered residue
            my $model_aligned_aa_start = $seq2->start;
            my $model_aligned_aa_end   = $seq2->end;

            my %orthologue = (
                id             => $model_translation_ID,
                start          => $model_aligned_aa_start,
                end            => $model_aligned_aa_end,
                sequence       => undef,
                entrez_gene_id => $entrez_gene_id,
                refseq_peptide => undef,   #will only be used to try to get NCBI Gene record
                uniprot_id => [],
                go         => undef
            );

            #attempt to fill in usefull IDs and other info
            my $model_translation_object
                = $config->{model_translation_adaptor}
                ->fetch_by_stable_id($model_translation_ID);
            $orthologue{sequence} = $model_translation_object->seq();

            my $xrefs = $model_translation_object->get_all_DBLinks();

            my @go_names = ();
            foreach my $xref (@{$xrefs}) {
                if ($xref->dbname() eq 'Uniprot/SWISSPROT') {
                    push(@{$orthologue{uniprot_id}}, $xref->display_id());
                }
                elsif ($xref->dbname() eq 'RefSeq_peptide') {
                    $orthologue{refseq_peptide} = $xref->display_id();
                }
                elsif ($xref->dbname() eq 'goslim_goa') {
                    my $go_id = $xref->display_id();
                    my $go_term;
                    my $go_name;
                    my $go_definition;
                    if (defined($config->{goa})) {
                        $go_term = $config->{goa}->fetch_by_accession($go_id);
                        $go_name       = $go_term->name();
                        $go_definition = $go_term->definition();
                    }

                    if (defined($go_name)) {
                        push(@go_names, "[$go_id]:$go_name");
                    }
                    else {
                        push(@go_names, "$go_id");
                    }
                }
            }
            $orthologue{go} = get_unique( \@go_names );

            if (!defined($orthologue{refseq_peptide})) {
                $orthologue{refseq_peptide} = $refseq_protein_id;
            }

            push(@model_orthologous_residues, \%orthologue);

        }
        else {
            die("Unexpected sequence ordering in alignment.");
        }
    }
    return \@model_orthologous_residues;
}

#returns array of phenotypes associated with model genomic coordinates
sub get_model_phenotypes_from_model_genomic_coords {
    my $genomic_coords = shift;

    my @phenotypes = ();

    my $slice = $config->{model_slice_adaptor}->fetch_by_region(
        undef,
        $genomic_coords->{chr},
        $genomic_coords->{start},
        $genomic_coords->{end},
        $genomic_coords->{strand}
    );

    my $vfs = $config->{model_variationfeature_adaptor}
        ->fetch_all_by_Slice($slice);

    foreach my $vf ( @{$vfs} ) {

        my $variation_name = $vf->variation_name();
        my $variation      = $config->{model_variation_adaptor}
            ->fetch_by_name($variation_name);

        my $vas = $config->{model_variationannotation_adaptor}
            ->fetch_all_by_Variation($variation);

        foreach my $va ( @{$vas} ) {

            my @phenotype_descriptors = ();
            if ( defined( $va->source_name() ) ) {
                push( @phenotype_descriptors,
                    'Source: ' . $va->source_name() );
            }
            if ( defined( $va->phenotype_name() ) ) {
                push( @phenotype_descriptors,
                    'Phenotype_name: ' . $va->phenotype_name() );
            }
            if ( defined( $va->phenotype_description() ) ) {
                push( @phenotype_descriptors,
                    'Description: ' . $va->phenotype_description() );
            }
            if ( defined( $va->variation_names() ) ) {
                push( @phenotype_descriptors,
                    'Variation_names: ' . $va->variation_names() );
            }
            if ( scalar(@phenotype_descriptors) > 0 ) {
                push( @phenotypes, join( ' ', @phenotype_descriptors ) );
            }
        }
    }
    return \@{get_unique(\@phenotypes)};
}

#returns hash containing chromosome, start, end, and strand of genomic sequence encoding
#amino acid described by $model_translation_id, $amino_acid_start, $amino_acid_end
sub get_model_genomic_coords_of_model_amino_acid {
    my $model_orthologous_residue = shift;

    my $model_translation_id = $model_orthologous_residue->{id};
    my $amino_acid_start     = $model_orthologous_residue->{start};
    my $amino_acid_end       = $model_orthologous_residue->{end};

    #obtain the model protein translation object
    my $model_translation = $config->{model_translation_adaptor}
        ->fetch_by_stable_id($model_translation_id);
    my $model_transcript = $config->{model_transcript_adaptor}
        ->fetch_by_translation_stable_id($model_translation_id);

    #    print "\$amino_acid_start is $amino_acid_start\n";
    #    print "\$amino_acid_end is $amino_acid_end\n";

    #determine the position on the transcript
    my $trmapper = $model_transcript->get_TranscriptMapper();

    my $model_slice = $config->{model_slice_adaptor}
        ->fetch_by_transcript_stable_id($model_transcript->stable_id);

    my $model_chromosome_name = $model_slice->seq_region_name();
    my $model_slice_start     = $model_slice->start();
    my $model_slice_end       = $model_slice->end();
    my $model_slice_strand    = $model_slice->strand();

#a list of Bio::EnsEMBL::Mapper::Gap and Bio::EnsEMBL::Mapper::Coordinate objects
    my @model_genomic_coords
        = $trmapper->pep2genomic($amino_acid_start, $amino_acid_end);

    #not sure how to handle more complicated coordinates, so skip
    if (scalar(@model_genomic_coords) > 1) {
        return;
    }

    my $model_aligned_aa_start_genomic = $model_genomic_coords[0]->start;
    my $model_aligned_aa_end_genomic   = $model_genomic_coords[0]->end;
    my $model_aligned_strand_genomic   = $model_genomic_coords[0]->strand;

#check that the slice start and end is consistent with the position obtained using pep2genomic
    if ((!( $model_slice_start <= $model_aligned_aa_start_genomic))
        || (!($model_slice_end >= $model_aligned_aa_end_genomic)))
    {
        return;
    }

    my %model_genomic_cords = (
        chr    => $model_chromosome_name,
        start  => $model_aligned_aa_start_genomic,
        end    => $model_aligned_aa_end_genomic,
        strand => $model_aligned_strand_genomic
    );
    return \%model_genomic_cords;

}
sub determine_overlapping_protein_domains {
    my ($translation, $altered_aa_start, $altered_aa_end) = @_;
    if ((!defined($translation))
        || (!defined($altered_aa_start))
        || (!defined($altered_aa_end)))
    {
        return;
    }

    my @overlapping_domains = ();

    my $pfeatures = $translation->get_all_DomainFeatures();
    while (my $pfeature = shift @{$pfeatures}) {
        my $overlaps    = 0;
        my $start       = $pfeature->start();
        my $end         = $pfeature->end();
        my $description = $pfeature->analysis()->logic_name();

        #interpro_ac() and idesc() seem to give whitespace sometimes
        if ((defined($pfeature->interpro_ac()))
            && ($pfeature->interpro_ac() =~ m/\S/))
        {
            $description = $description . ' ' . $pfeature->interpro_ac();
        }
        if ((defined($pfeature->idesc()))
            && ($pfeature->idesc() =~ m/\S/))
        {
            $description = $description . ' ' . $pfeature->idesc();
        }

        if (overlaps(
                s1 => $start,
                e1 => $end,
                s2 => $altered_aa_start,
                e2 => $altered_aa_end)
            )
        {

#   print "Protein feature '$description' starts at $start and ends at $end and overlaps with " .$altered_aa_start} . " to " .$altered_aa_end} . "\n";
            push( @overlapping_domains, $description );
        }
        else {

#   print "Protein feature '$description' starts at $start and ends at $end and does not overlap with " .$altered_aa_start} . " to " .$altered_aa_end} . "\n";
        }
    }
    

    return (scalar(@overlapping_domains) > 0) ?
        get_unique( \@overlapping_domains ) :
        undef;
}

sub overlaps {
    my %args = (@_);

    if (   ( $args{s1} >= $args{s2} )
        && ( $args{e1} <= $args{e2} ) )
    {
        return 1;
    }
    elsif (( $args{s1} <= $args{s2} )
        && ( $args{e1} >= $args{e2} ) )
    {
        return 1;
    }
    elsif (
        ( $args{s1} <= $args{s2} )
        && (   ( $args{e1} <= $args{e2} )
            && ( $args{e1} >= $args{s2} ) )
        )
    {
        return 1;
    }
    elsif (
        ( $args{e1} >= $args{e2} )
        && (   ( $args{s1} >= $args{s2} )
            && ( $args{s1} <= $args{e2} ) )
        )
    {
        return 1;
    }
    return 0;
}

sub determine_context_average_percent_identity {
    my ($context_alignments) = @_;

    if ((!defined($context_alignments))
        || (scalar(@{$context_alignments}) == 0))
    {
        return;
    }
    
    #return value
    my $context_conservation;

    my $identity_sum   = 0;
    my $identity_count = 0;
    foreach my $subsequence_alignment (@{$context_alignments})
    {
        my $percent_identity = get_percent_identity_between_two_seqs(
            $subsequence_alignment->{reference},
            $subsequence_alignment->{homolog}
        );
        if (defined($percent_identity)) {
            $identity_sum = $identity_sum + $percent_identity;
            $identity_count++;
        }
    }

    if ( $identity_count > 0 ) {
        $context_conservation
            = sprintf( "%.1f", ( $identity_sum / $identity_count ) );
    }
    return $context_conservation;
}

sub determine_context_alignments {
    my ($changes_protein, $altered_aa_start, $protein_sequence, $gene_id, $protein_id) = @_;

    if (! $changes_protein) {
        return;
    }
    
    #return values
    my (@context_alignments, @homolog_species_context);
    
    my $context_start = $altered_aa_start
        - $config->{flanking_for_context_conservation};
    my $context_end = $altered_aa_start
        + $config->{flanking_for_context_conservation};

    if ( $context_start < 1 ) {
        $context_start = 1;
    }
    if ($context_end > length($protein_sequence)) {
        $context_end = length($protein_sequence);
    }

    #$config->{ma} is a Bio::EnsEMBL::Compara::DBSQL::MemberAdaptor
    #$member is a Bio::EnsEMBL::Compara::Member
    my $member = $config->{ma}
        ->fetch_by_source_stable_id('ENSEMBLGENE', $gene_id);

    #Rarely the $member object may be undef
    if (!defined($member)) {
        return;
    }

    #$config->{ha} is a Bio::EnsEMBL::Compara::DBSQL::HomologyAdaptor
    #$homologies is a list of Bio::EnsEMBL::Compara::Homology objects

    my $homologies = [];

    if (defined($config->{comparison_species})) {

        #obtain orthologues from species of interest
        foreach my $species (keys(%{$config->{comparison_species}})) {
            push(
                @{$homologies},
                @{$config->{ha}->fetch_all_by_Member_paired_species($member,
                    $species, ['ENSEMBL_ORTHOLOGUES'])}
            );
        }
    }
    else {

        #obtain orthologues from all species
        push(
            @{$homologies},
            @{$config->{ha}->fetch_all_by_Member_method_link_type($member,
                    'ENSEMBL_ORTHOLOGUES')}
        );
    }

    foreach my $homology (@{$homologies}) {

        #$homologues is an array ref of (2) Bio::EnsEMBL::Compara::Member objects
        #$aln is Bio::SimpleAlign object

        #added eval block on 2009-10-22 because
        #$aln = $homology->get_SimpleAlign();
        #throws exception for ENSBTAT00000060548
        my $homologues = undef;
        my $taxon1     = undef;
        my $taxon2     = undef;
        my $aln        = undef;

        eval {
            $homologues = $homology->gene_list();
            $taxon1     = $$homologues[0]->taxon;
            $taxon2     = $$homologues[1]->taxon;
            $aln        = $homology->get_SimpleAlign();
        };
        if ($@) {
            next;
        }

        if (!(defined($aln))) {
            next;
        }

        #confirm that reference protein is in the alignment
        if (!(scalar($aln->each_seq_with_id($protein_id)) == 1))
        {
            next;
        }

        #get an alignment containing the SNP-affected column and flanking sequence
        my $sub_align_context = undef;

        eval {
            my $col_context_start
                = $aln->column_from_residue_number(
                $protein_id,
                $context_start);

            my $col_context_end = $aln->column_from_residue_number(
                $protein_id, $context_end);

            $sub_align_context
                = $aln->slice($col_context_start, $col_context_end);
        };
        if ($@) {
            next;
        }

        #confirm that there are two sequences in the subalignment
        #(sequences are excluded if the alignment slice contains no residues from the sequence)
        my $context_string = undef;
        if ((defined($sub_align_context))
            && (scalar($sub_align_context->each_seq()) == 2))
        {
            my $seq1_context = $sub_align_context->get_seq_by_pos(1);
            my $seq2_context = $sub_align_context->get_seq_by_pos(2);

            my %seq_hash = (reference => undef, homolog => undef);

           #in all the examples I've examined, $seq1 is the reference sequence
            if ($seq1_context->display_id eq $protein_id) {

                $seq_hash{reference} = $seq1_context->seq();
                $seq_hash{homolog}   = $seq2_context->seq();

                push( @context_alignments, \%seq_hash );
                push(
                    @homolog_species_context,
                    space_to_underscore( get_taxon_name($taxon2) )
                );
            }
            else {
                die("Unexpected sequence ordering in alignment.");
            }
        }
    }
    return (\@context_alignments, \@homolog_species_context);
}

sub determine_alignment_score_change {
    my ($matrix, $max_alignment_score_change, $aligned_homolog_residues, $amino_acid_reference, $amino_acid_variant) = @_;

    if ((! defined($aligned_homolog_residues))
        || (scalar(@{$aligned_homolog_residues}) == 0)) {
        return;
    }

    my $reference_score = get_average_similarity_score(
        $matrix,
        $amino_acid_reference,
        $aligned_homolog_residues
    );
    my $variant_score = get_average_similarity_score(
        $matrix,
        $amino_acid_variant,
        $aligned_homolog_residues
    );

    my $alignment_score_change = sprintf("%.3f",
        ($variant_score - $reference_score) / $max_alignment_score_change);
    return $alignment_score_change;
}

sub determine_reference_amino_acid_conservation {
    my ($aligned_homolog_residues, $amino_acid_reference) = @_;

    if ((! defined($aligned_homolog_residues))
        || (scalar(@{$aligned_homolog_residues}) == 0)) {
        return;
    }

    my $reference_amino_acid_conservation = get_percent_identity_column(
        $amino_acid_reference,
        $aligned_homolog_residues
    );

    if (defined($reference_amino_acid_conservation)) {
        my $reference_amino_acid_conservation
            = sprintf("%.1f", $reference_amino_acid_conservation);
    }
    return $reference_amino_acid_conservation;
}

sub determine_reference_amino_acid_score {
    my ($matrix, $aligned_homolog_residues, $amino_acid_reference) = @_;
 
    if ((! defined($aligned_homolog_residues))
        || (scalar(@{$aligned_homolog_residues}) == 0)) {
        return;
    }

    my $reference_amino_acid_score = get_score_column(
        $matrix,
        $amino_acid_reference,
        $aligned_homolog_residues
    );

    if (defined($reference_amino_acid_score)) {
        $reference_amino_acid_score
            = sprintf("%.1f", $reference_amino_acid_score);
    }
    return $reference_amino_acid_score;
}

sub determine_aligned_homolog_residues {
    my ($changes_protein, $gene_id, $protein_id, $altered_aa_start, $altered_aa_end, $amino_acid_reference, $consequences_ranks_so, $comparison_species) = @_;
    if ( (!$changes_protein)
        || (!defined($gene_id))
        || (!defined($protein_id))
        || (!defined($altered_aa_start))
        || (!defined($altered_aa_end))
        || (!defined($amino_acid_reference)))
    {
        return;
    }
    #return values
    my @aligned_homolog_residues;
    my @homolog_species;
    
    if (!($changes_protein)) {
        return;
    }

    #$config->{ma} is a Bio::EnsEMBL::Compara::DBSQL::MemberAdaptor
    #$member is a Bio::EnsEMBL::Compara::Member
    my $member = $config->{ma}
        ->fetch_by_source_stable_id('ENSEMBLGENE', $gene_id);

    #Rarely the $member object may be undef
    if (!defined($member)) {
        return;
    }

    #$config->{ha} is a Bio::EnsEMBL::Compara::DBSQL::HomologyAdaptor
    #$homologies is a list of Bio::EnsEMBL::Compara::Homology objects

    my $homologies = [];
    if (defined($config->{comparison_species})) {

        #obtain orthologues from species of interest
        foreach my $species (keys( %{$config->{comparison_species}})) {
            push(
                @{$homologies},
                @{$config->{ha}
                    ->fetch_all_by_Member_paired_species($member,
                    $species, ['ENSEMBL_ORTHOLOGUES'])
                }
            );
        }
    }
    else {

        #obtain orthologues from all species
        push(
            @{$homologies},
            @{$config->{ha}->fetch_all_by_Member_method_link_type($member,
                    'ENSEMBL_ORTHOLOGUES')}
        );
    }

    foreach my $homology (@{$homologies}) {

        #$homologues is an array ref of (2) Bio::EnsEMBL::Compara::Member objects
        #$aln is Bio::SimpleAlign object

        #added eval block on 2009-10-22 because
        #$aln = $homology->get_SimpleAlign();
        #throws exception for ENSBTAT00000060548
        my $homologues = undef;
        my $taxon1     = undef;
        my $taxon2     = undef;
        my $aln        = undef; #Bio::SimpleAlign object

        eval {
            $homologues = $homology->gene_list();
            $taxon1     = $$homologues[0]->taxon;
            $taxon2     = $$homologues[1]->taxon;
            $aln        = $homology->get_SimpleAlign();
        };
        if ($@) {
            next;
        }

        if (!(defined($aln))) {
            next;
        }

        #confirm that reference protein is in the alignment
        #Gets an array of Seq objects from the
        #alignment, the contents being those sequences
        #with the given name (there may be more than one)
        if (!(scalar($aln->each_seq_with_id($protein_id)) == 1))
        {
            next;
        }

#TODO:verify indel
        #get the column containing the residue encoded by the SNP-containing codon
        my $col = undef;
        my $col_end = undef;

        #This function gives the position in the alignment
        #(i.e. column number) of the given residue number in the
        #sequence with the given name.
        #column_from_residue_number can throw an exception if the altered residue
        #in the reference is the one that encodes the stop codon
        eval {
            $col = $aln->column_from_residue_number(
                $protein_id,
                $altered_aa_start
            );
        };
        if ($@) {
            next;
        }
        eval {
            $col_end = $aln->column_from_residue_number(
                $protein_id,
                $altered_aa_end
            );
        };
        if ($@) {
            next;
        }

        #Bio::SimpleAlign->slice
        #Creates a slice from the alignment inclusive of start and
        #end columns.  Sequences with no residues in the slice are
        #excluded from the new alignment and a warning is printed.
        #Slice beyond the length of the sequence does not do
        #padding.
        #get an alignment containing only the desired column
        #my $sub_align = $aln->slice( $col, $col );
        my $sub_align = undef;
        #prevent annoying warnings; seems to require Bioperl-live, not Bioperl-1.2.3
        if (!defined($config->{verbose})){
            $SIG{'__WARN__'} = sub {return 1};
        }
        eval { $sub_align = $aln->slice($col, $col_end); };
        if ($@) {
            next;
        }

        #confirm that there are two sequences in the subalignment
        #(sequences are excluded if the alignment slice contains no residues from the sequence)
        if (!(scalar($sub_align->each_seq()) == 2)) {
            next;
        }

        my $seq1 = $sub_align->get_seq_by_pos(1);
        my $seq2 = $sub_align->get_seq_by_pos(2);

        #in all the examples I've examined, $seq1 is the reference sequence
        if ($seq1->display_id eq $protein_id) {
            if (uc($seq1->seq()) ne
                uc($amino_acid_reference))
            {
#TO: Do we need to deal with this side effect somewhere?
                #if the residue in the alignment is selenocysteine (U), change the reference to U
                #as U in reference seems to be given as * in string returned by $con->pep_allele_string()
                if ((uc($seq1->seq()) eq 'U')
                    && ($amino_acid_reference eq '*'))
                {
                    $amino_acid_reference = 'U';

                    #assume base change in U-coding residue is NON_SYNONYMOUS_CODING
                    #$con->consequence_type seems to be incorrect for sites encoding 'U'
                    
                    #if ($consequence eq 'STOP_LOST') {
                    if (exists $consequences_ranks_so->{stop_lost}){
                        #$consequence = 'NON_SYNONYMOUS_CODING';
                        delete $consequences_ranks_so->{stop_lost};
#TODO: correct rank
                        $consequences_ranks_so->{non_synonymous_codon} = 0;
                    }
                }
                else {
#TODO: something other than ignore things like this
#VCF: 17	39340795	splice_acceptor_variant	ACGGCAGCAGCTGGACATACCACAGCTGGGGTGGCAGGTGGTCTGACAGCAGAGTGGG	A	.	.	.
#protein_id: ENSP00000381489
#seq1->display_id: ENSP00000381489
#amino_acid_reference: RPLCCQTTCHPSCGMSSCCR
#seq1->seq: RPLCCQTTCHP---SCGMSSCCR
                    #die("The alignment '"
                    #        . $seq1->seq()
                    #        . "' and reference '$amino_acid_reference' amino acids do not match for sequence $protein_id."
                    #);
                    next;
                }
            }
            push(
                @aligned_homolog_residues,
                $seq2->seq()
            );
            push(
                @homolog_species,
                space_to_underscore(get_taxon_name($taxon2))
            );
        }
        else {
            die("Unexpected sequence ordering in alignment.");
        }
    }
    return (\@aligned_homolog_residues, \@homolog_species)
}

#return flanking from forward strand of genome
sub get_flanking_sequence {
    my $options = shift;
    my $snp     = shift;

    my $start = $snp->{Chromosome_Position} - $config->{flanking};
    if ( $start < 1 ) {
        $start = 1;
    }
    my $end = $snp->{Chromosome_Position} + $config->{flanking};

    my $left_slice
        = $config->{sa}->fetch_by_region( undef, $snp->{Chromosome}, $start,
        $snp->{Chromosome_Position} - 1 );
    my $snp_site_slice = $config->{sa}->fetch_by_region(
        undef, $snp->{Chromosome},
        $snp->{Chromosome_Position},
        $snp->{Chromosome_Position}
    );
    my $right_slice
        = $config->{sa}->fetch_by_region( undef, $snp->{Chromosome},
        $snp->{Chromosome_Position} + 1, $end );

    #do not indicate known SNPs in flanking if user does not want this info
    #or if database of known SNPs is not available
    if (   ( !defined( $config->{degenerate_flanking} ) )
        || ( $options->{is_fake_vfa} ) )
    {
        $snp->{flanking}
            = uc( $left_slice->seq() ) . '['
            . $snp->{allele_string} . ']'
            . uc( $right_slice->seq() );
        return;
    }

    #indicate known SNPs using lowercase IUPAC characters

    my $combined = uc( $left_slice->seq() )

        #        . lc( $snp_site_slice->seq() )
        . lc( $snp->{Chromosome_Reference} ) . uc( $right_slice->seq() );

    my $bases_to_snp = length( $left_slice->seq() ) + 1;

    my $flank_start
        = $snp->{Chromosome_Position} - length( $left_slice->seq() );
    my $flank_end
        = $snp->{Chromosome_Position} + length( $right_slice->seq() );

    my $slice
        = $options->{sa}
        ->fetch_by_region( undef, $snp->{Chromosome}, $flank_start,
        $flank_end );

    #return all variations defined in $slice
    my $vfs = $options->{vfa}->fetch_all_by_Slice($slice);
    foreach my $vf ( @{$vfs} ) {

        #without transform() start will be '1'
        $vf = $vf->transform('toplevel');
        my @alleles              = split( /\//, $vf->allele_string );
        my $vf_strand            = $vf->strand;
        my $vf_start_on_flanking = $vf->start - $flank_start + 1;
        my $vf_end_on_flanking   = $vf->end - $flank_start + 1;

        #don't handle indels
        if ( $vf_start_on_flanking != $vf_end_on_flanking ) {
            next;
        }

        my $iupac_snp = convert_base_list_to_iupac( \@alleles );

        #genomic sequence will represent forward strand.
        #need to convert SNP alleles to forward strand if necessary.
        if ( $vf_strand == -1 ) {
            $iupac_snp = complement($iupac_snp);
        }

        #incorporate genomic reference base into iupac code
        $iupac_snp
            = convert_base_list_to_iupac(
            [ $iupac_snp, substr( $combined, $vf_start_on_flanking - 1, 1 ) ]
            );

        $combined = replace_base( $vf_start_on_flanking, $combined,
            lc($iupac_snp) );

    }

    my $iupac_snp = convert_base_list_to_iupac(
        [   $snp->{Chromosome_Reference},
            $snp->{Chromosome_Reads},
            uc( get_base( length( $left_slice->seq() ) + 1, $combined ) )
        ]
    );
    my $square_bracket_snp = convert_iupac_to_square_brackets($iupac_snp);

    $snp->{flanking}
        = replace_base( $bases_to_snp, $combined, $square_bracket_snp );

}

sub get_base {

    #first base is position '1'
    my $position = shift;
    my $sequence = shift;
    return substr( $sequence, $position - 1, 1 );
}

sub replace_base {

    #first base is position '1'
    my $position    = shift;
    my $sequence    = shift;
    my $replacement = shift;

    my $upstream = substr( $sequence, 0, $position - 1 );
    my $snp      = substr( $sequence, $position - 1, 1 );
    my $downstream = substr( $sequence, $position, length($sequence) - $position );

    my $new_sequence = $upstream . $replacement . $downstream;
    return $new_sequence;
}

sub space_to_underscore {
    my $text = shift;
    $text =~ s/\s{1,}/_/g;
    return $text;
}

#determines maximum possible alignment score change value
sub get_max_alignment_score_change {
    my $matrix   = shift;
    my @residues = (
        'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K',
        'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X', '*'
    );
    my $max_score = undef;
    foreach my $row_residue (@residues) {
        my $row_max = undef;
        my $row_min = undef;

        foreach my $column_residue (@residues) {
            my $score = $matrix->get_entry( $row_residue, $column_residue );
            if ( defined($score) ) {
                if ( ( !defined($row_max) ) || ( $score > $row_max ) ) {
                    $row_max = $score;
                }
                if ( ( !defined($row_min) ) || ( $score < $row_min ) ) {
                    $row_min = $score;
                }
            }
        }

        if (   ( !defined($max_score) )
            || ( ( $row_max - $row_min ) > $max_score ) )
        {
            $max_score = $row_max - $row_min;
        }
    }
    return $max_score;
}

#compares $residue to each amino acid in $array_of_aligned
#and determines average similarity score
sub get_average_similarity_score {
    my $matrix           = shift;
    my $residue          = shift;
    my $array_of_aligned = shift;

    my $score_sum   = 0;
    my $score_count = 0;

    foreach my $aligned ( @{$array_of_aligned} ) {
        my $score = $matrix->get_entry( uc($residue), uc($aligned) );
        if ( defined($score) ) {
            $score_sum = $score_sum + $score;
            $score_count++;
        }
    }

    my $score = undef;
    if ( $score_count == 0 ) {
        $score = undef;
    }
    else {
        $score = $score_sum / $score_count;
        $score = sprintf( "%.2f", $score );
    }
    return $score;
}

#compares $residue to each amino acid in $array_of_aligned
#and determines percentage of residues in $array_of_aligned
#that are identical.
#$array_of_aligned is expected to only contain residues (no gaps)
#nonetheless gaps are skipped by this function.
#
#G vs GA will give 50
#G vs -G will give 100
#
#Is equivalent to "Cident" formula in "Correlated Mutations: A Hallmark of Phenotypic Amino Acid Substitutions"
#by Kowarsch et al., 2010
sub get_percent_identity_column {
    my $residue          = shift;
    my $array_of_aligned = shift;

    my $match_count   = 0;
    my $checked_count = 0;

    foreach my $aligned ( @{$array_of_aligned} ) {
        if ( ( is_gap($residue) ) || ( is_gap($aligned) ) ) {
            next;
        }
        $checked_count++;
        if ( uc($residue) eq uc($aligned) ) {
            $match_count++;
        }
    }
    if ( $checked_count != 0 ) {
        my $percent_identity
            = sprintf( "%.1f", ( $match_count / $checked_count ) * 100 );

        #   print "residue: $residue\n";
        #   print "aligned: " . join(',', @{$array_of_aligned}) . "\n";
        #   print "percent_identity: $percent_identity\n";

        return $percent_identity;
    }
    return undef;

}

#compares $residue to each amino acid in $array_of_aligned
#and determines sum of scores using scoring matrix.
#this value is then divided by the score obtained if all
#residues in $array_of_aligned were to match $residue
#$array_of_aligned is expected to only contain residues (no gaps)
#nonetheless gaps are skipped by this function.
#
#G vs GA will give score(G,G) + score(G,A) / 2 * score(G,G)
#G vs -G will give score(G,G) / score(G,G)
#
#Is equivalent to "Cblosum" formula in "Correlated Mutations: A Hallmark of Phenotypic Amino Acid Substitutions"
#by Kowarsch et al., 2010
sub get_score_column {
    my $matrix           = shift;
    my $residue          = shift;
    my $array_of_aligned = shift;

    my $score_sum     = 0;
    my $max_score_sum = 0;

    my $max_score = $matrix->get_entry( uc($residue), uc($residue) );
    if ( !defined($max_score) ) {
        return undef;
    }

    foreach my $aligned ( @{$array_of_aligned} ) {
        if ( ( is_gap($residue) ) || ( is_gap($aligned) ) ) {
            next;
        }

        my $score = $matrix->get_entry( uc($residue), uc($aligned) );
        if ( defined($score) ) {
            $score_sum = $score_sum + $score;
        }

        $max_score_sum = $max_score_sum + $max_score;

    }
    if ( $max_score_sum != 0 ) {
        my $score = sprintf( "%.1f", $score_sum / $max_score_sum );

        #   print "residue: $residue\n";
        #   print "aligned: " . join(',', @{$array_of_aligned}) . "\n";
        #   print "column score: $score\n";

        return $score;
    }
    return undef;

}

#determines percent identity between two sequences
#by looking for character matches.
#
#Aligning gaps are ignored
#
#--G
#--G will give 100 match
#
#-GA
#-GG will give 50 match
sub get_percent_identity_between_two_seqs {

    my $seq1 = shift;
    my $seq2 = shift;

    my @seq1 = split( //, $seq1 );
    my @seq2 = split( //, $seq2 );

    my $match_count   = 0;
    my $checked_count = 0;
    for ( my $i = 0; $i < scalar(@seq1); $i++ ) {
        if ( ( is_gap( $seq1[$i] ) ) || ( is_gap( $seq2[$i] ) ) ) {
            next;
        }

        $checked_count++;
        if ( uc( $seq1[$i] ) eq uc( $seq2[$i] ) ) {
            $match_count++;
        }
    }
    if ( $checked_count != 0 ) {
        my $percent_identity
            = sprintf( "%.1f", ( $match_count / $checked_count ) * 100 );

        #   print "seq1: $seq1\n";
        #   print "seq2: $seq2\n";
        #   print "percent_identity: $percent_identity\n";

        return $percent_identity;
    }
    return undef;
}

sub get_degenerate_hash_ref {
    my %degenerate_bases = (
        A => ['A'],
        C => ['C'],
        G => ['G'],
        T => ['T'],
        R => [ 'A', 'G' ],
        Y => [ 'C', 'T' ],
        S => [ 'C', 'G' ],
        W => [ 'A', 'T' ],
        K => [ 'G', 'T' ],
        M => [ 'A', 'C' ],
        B => [ 'C', 'G', 'T' ],
        D => [ 'A', 'G', 'T' ],
        H => [ 'A', 'C', 'T' ],
        V => [ 'A', 'C', 'G' ],
        N => [ 'A', 'C', 'G', 'T' ]
    );

    return \%degenerate_bases;
}

sub get_bases_to_degenerate_hash_ref {
    my %bases_to_degenerate = (
        A    => 'A',
        C    => 'C',
        G    => 'G',
        T    => 'T',
        AG   => 'R',
        CT   => 'Y',
        CG   => 'S',
        AT   => 'W',
        GT   => 'K',
        AC   => 'M',
        CGT  => 'B',
        AGT  => 'D',
        ACT  => 'H',
        ACG  => 'V',
        ACGT => 'N'
    );

    return \%bases_to_degenerate;
}

sub is_gap {
    my $residue = shift;
    if ( $residue =~ m/[\.\-]/ ) {
        return 1;
    }
    return 0;
}

sub reverse_complement {
    my $seq = shift;
    $seq = reverse complement($seq);
    return $seq;
}

sub complement {
    my $seq = shift;
    $seq =~ tr/gatcryswkmbdhvnGATCRYSWKMBDHVN/ctagyrswmkvhdbnCTAGYRSWMKVHDBN/;
    return $seq;
}

sub get_taxon_name {
    #Bio::EnsEMBL::Compara::NCBITaxon
    my $taxon = shift;
    if (defined( $taxon->binomial)) {
        return $taxon->binomial;
    }
    if ( defined($taxon->name)) {
        return $taxon->name;
    }
    return '?';
}

#NGS-SNP

#ZILA
sub generatePhyloTreeDistArray {
	my $species = shift;
	my $speciesListCommon_ptr = shift;
  my @speciesListCommon = @{$speciesListCommon_ptr};
	my $speciesListLatin_ptr = shift;
  my @speciesListLatin = @{$speciesListLatin_ptr};
  #return will be ptr to hash of species:distance
  my %distanceHash;
  
	# get index of our species, default to Homo_sapiens
 	my $speciesIndex = first { lc($speciesListLatin[$_]) eq lc($species) } 0..$#speciesListLatin;
  if (! defined($speciesIndex)){
    $species = 'Homo_sapiens';
    $speciesIndex = first { lc($speciesListLatin[$_]) eq lc($species) } 0..$#speciesListLatin;
  }
	# check defined species exists
	die("ERROR: Could not find species \"", $species, "\" in Ensembl species list\n") unless defined $speciesIndex;

	# default phylogenetic tree; branch lengths for Ensembl species phylogenetic tree obtained from http://tinyurl.com/ensembltree on 2011-08-20
  my $default_phyloTree = lc("(((((((((((((((((((((((Homo_sapiens:0.0067,Pan_troglodytes:0.006667):0.00225,Gorilla_gorilla:0.008825):0.00968,Pongo_abelii:0.018318):0.00717,Nomascus_leucogenys:0.025488):0.00717,(Macaca_mulatta:0.007853,?Papio_hamadryas:0.007637):0.029618):0.021965,Callithrix_jacchus:0.066131):0.05759,Tarsius_syrichta:0.137823):0.011062,(Microcebus_murinus:0.092749,Otolemur_garnettii:0.129725):0.035463):0.015494,Tupaia_belangeri:0.186203):0.004937,(((((Mus_musculus:0.084509,Rattus_norvegicus:0.091589):0.197773,Dipodomys_ordii:0.211609):0.022992,Cavia_porcellus:0.225629):0.01015,Spermophilus_tridecemlineatus:0.148468):0.025746,(Oryctolagus_cuniculus:0.114227,Ochotona_princeps:0.201069):0.101463):0.015313):0.020593,((((Vicugna_pacos:0.107275,(Tursiops_truncatus:0.064688,(Bos_taurus:0.061796,?Ovis_aries:0.061796):0.061796):0.025153):0.0201675,Sus_scrofa:0.079):0.0201675,((Equus_caballus:0.109397,(Felis_catus:0.098612,(Ailuropoda_melanoleuca:0.051229,Canis_familiaris:0.051229):0.051229):0.049845):0.006219,(Myotis_lucifugus:0.14254,Pteropus_vampyrus:0.113399):0.033706):0.004508):0.011671,(Erinaceus_europaeus:0.221785,Sorex_araneus:0.269562):0.056393):0.021227):0.023664,(((Loxodonta_africana:0.082242,Procavia_capensis:0.155358):0.02699,Echinops_telfairi:0.245936):0.049697,(Dasypus_novemcinctus:0.116664,Choloepus_hoffmanni:0.096357):0.053145):0.006717):0.234728,(Monodelphis_domestica:0.125686,Macropus_eugenii:0.122008):0.2151):0.071664,Ornithorhynchus_anatinus:0.456592):0.109504,((((Gallus_gallus:0.041384,Meleagris_gallopavo:0.041384):0.041384,Anas_platyrhynchos:0.082768):0.082768,Taeniopygia_guttata:0.171542):0.199223,Anolis_carolinensis:0.489241):0.105143):0.172371,Xenopus_tropicalis:0.855573):0.311354,(((Tetraodon_nigroviridis:0.224159,Takifugu_rubripes:0.203847):0.195181,(Gasterosteus_aculeatus:0.316413,Oryzias_latipes:0.48197):0.05915):0.32564,Danio_rerio:0.730752):0.147949):0.526688,?Petromyzon_marinus:0.526688),(Ciona_savignyi:0.8,Ciona_intestinalis:0.8)Cionidae:0.6)Chordata:0.2,(?Apis_mellifera:0.9,(((?Aedes_aegypti:0.25,?Culex_quinquefasciatus:0.25):0.25,?Anopheles_gambiae:0.5)Culicinae:0.2,Drosophila_melanogaster:0.8)Diptera:0.1)Endopterygota:0.7)Coelomata:0.1,Caenorhabditis_elegans:1.7)Bilateria:0.3,Saccharomyces_cerevisiae:1.9)Fungi_Metazoa_group:0.3)");
	# dynamically get tree with branch lengths
	my $url = 'http://tinyurl.com/ensembltree';
	my $content = get $url;
  my $phyloTree = defined $content ? lc($content) : $default_phyloTree;

	my @distanceToSpecies = ((-1) x scalar(@speciesListLatin));
  for (my $i = 0; $i < scalar(@speciesListLatin); $i++){
    $distanceHash{$speciesListLatin[$i]} = -1;
  }
	my $specieslc;
	my $currspecieslc;

	# trim branches of tree that aren't in our species list
	my @trimmedtree = split(/([\(\):,])/, $phyloTree);
	for(my $i = 0; $i < scalar(@trimmedtree); $i++) {
		my $item = $trimmedtree[$i];
		$item =~ s/\?//;
		my $itemIndex = first { lc($speciesListLatin[$_]) eq $item } 0..$#speciesListLatin;
		if (($item =~ m/[a-z_]+/)  && !(defined $itemIndex)) {
			#print "Item $item not found. Removing from tree.\n";
			my @replacement = ("","","");
			splice(@trimmedtree, $i, 3, @replacement);			
		}
	}
	$phyloTree = join("", @trimmedtree);
	#print $phyloTree."\n";

	# calculate distances between our species and every other species
	for(my $i = 0; $i < scalar(@speciesListLatin); $i++) {
		
		# set distance to our own species to zero
		if ( $i == $speciesIndex) {
			$distanceToSpecies[$i] = 0;
      $distanceHash{$speciesListLatin[$i]} = 0;
		}
		else {
			# does our species and the current species exist in our phylogenetic tree?
			$specieslc = lc($species);
			$currspecieslc = lc($speciesListLatin[$i]);
			if (($phyloTree =~ m/$specieslc/) && ($phyloTree =~ m/$currspecieslc/)) {
				#print "Found ".$specieslc." and ".$currspecieslc."\n";

				# make a copy of our tree
				# e.g. ((((T:0.7,F:0.3):0.1,(G:0.4,O:0.5):0.6):0.2;D:0.8):0.9;A:1.1)
				my $tree = $phyloTree;

				# trim all leaves of the tree that aren't our species or the current species
				# e.g. if our species are F and D --> ((((,F:0.3):0.1,(,):0.6):0.2;D:0.8):0.9,)
				for(my $j = 0; $j < scalar(@speciesListLatin); $j++) {
					my $jspecies = lc($speciesListLatin[$j]);
					if (($jspecies ne $specieslc) && ($jspecies ne $currspecieslc)) {
						$tree =~ s/$jspecies:\d.[\d]+//g;
					}
				}

				# trim all branches of the tree that aren't on the path between the two species of interest
				# e.g. if our species are F and D --> (((,F:0.3):0.1,):0.2;D:0.8)
				while ($tree =~ m/\(,\):\d.[\d]+/) {
					$tree =~ s/\(,\):\d.[\d]+//g;
				}
				while ($tree =~ m/\(,\)/) {
					$tree =~ s/\(,\)//g;
				}
				# find set of parentheses that encloses our species and discard everything outside
				# code adapted from http://www.perlmonks.org/?node_id=660316
				my @queue = ( $tree );			
				my $regex = qr/
					(			# start of bracket 1
					\(			# match an opening parentheses
						(?:               
						[^\(\)]++	# one or more non parentheses, non backtracking
							|                  
						(?1)		# recurse to parentheses 1
						)*                 
					\)			# match a closing parentheses
					)			# end of parentheses 1
					/x;

				$" = "\n\t";

				my @potentials;
				while( @queue )
				{
					my $string = shift @queue;
					my @groups = $string =~ m/$regex/g;
					#print "Found:\n\t@groups\n\n" if @groups;
					if (($string =~ m/$specieslc/) && ($string =~ m/$currspecieslc/)) {
						push @potentials, $string;
					}
					unshift @queue, map { s/^\(//; s/\)$//; $_ } @groups;
				}
				#print "Potentials: ".join("\n", @potentials)."\n";
				my $shortest = $potentials[0];
				foreach my $potential (@potentials) {
					$shortest = $potential if length($potential) < length($shortest);
				    }
				#print "Shortest: (".$shortest.")\n";

				# add together all remaining numbers in our string
				# e.g. 0.3 + 0.1 + 0.2 + 0.8 --> 1.4
				my @splittree = split(/([\(\):,])/, $shortest);
				my $sum = 0;
				foreach my $item (@splittree) {
					if ($item =~ /\d.[\d]+/) {	# if it's a number
						$sum += $item;	# branch lengths to scale
						#$sum += 0.1;	# fixed branch lengths
					}
				}
				#print "Distance between ".$speciesListCommon->[$speciesIndex]." and ".$speciesListCommon->[$i]." is: ".$sum."\n";
				$distanceToSpecies[$i] = $sum;
        $distanceHash{$speciesListLatin[$i]} = $sum;
			}
		}
	}
  for (my $i = 0; $i < scalar(@speciesListLatin); $i++){
    $distanceHash{$speciesListLatin[$i]} = $distanceToSpecies[$i];
  }
	return \%distanceHash;
}

#rearrange orthologues by evolutionary distance of species
sub orderOrthologues {
    my $AAsInOrthologues_ptr = shift;
    my @AAsInOrthologues = @{$AAsInOrthologues_ptr};
    my $orthologues_ptr = shift;
    my @orthologues = @{$orthologues_ptr};
    #species:distance
    my $distHash_ptr = shift;
    my %distHash = %{$distHash_ptr};
    
    #uniqueify species/AAs
    my %distinctOrthologues;
    for (my $i=0;$i<scalar(@orthologues);$i++){
      $distinctOrthologues{$orthologues[$i]}->{$AAsInOrthologues[$i]}++;
    }
    
    #sort species by evolutionary distance
    my %positives = map {$_=>$distHash{$_}} grep {$distHash{$_}>=0} keys %distHash;
    my @sortedSpecies = sort {$positives{$a} <=> $positives{$b}} keys %positives;
    #append any species with no or negative distance
    while(my($species,$aas) = each(%distinctOrthologues)){
      push @sortedSpecies, $species unless exists($positives{$species});
    }
    
    #push species and AAs to new arrays ordered by evolutionary distance
    my @sortedAAsInOrthologues = ();
    my @sortedOrthologues = ();
    for my $thisSpecies(@sortedSpecies){
        while (my($aa, $count) = each(%{$distinctOrthologues{$thisSpecies}})){
            push @sortedOrthologues, $thisSpecies;
            push @sortedAAsInOrthologues, $aa;
        }
    }
    return \(@sortedAAsInOrthologues, @sortedOrthologues);
}

sub getEnsemblSpecies {
    # assign arrays of species common and latin names
    #defaults downloaded from http://www.ensembl.org/info/about/species.html on 2011-08-20
    my @speciesListCommon = qw(Alpaca Gorilla Pig Anole_Lizard Guinea_Pig Pika Armadillo Hedgehog Platypus Horse Rabbit Bushbaby Human Rat C.elegans Hyrax S.cerevisiae C.intestinalis Kangaroo_rat C.savignyi Shrew Cat Lesser_hedgehog_tenrec Sloth Chicken Macaque Squirrel Chimpanzee Marmoset Stickleback Cow Medaka Tarsier Dog Megabat Tetraodon Dolphin Microbat Tree_Shrew Mouse Turkey Elephant Mouse_Lemur Wallaby Fruitfly Opossum X.tropicalis Fugu Orangutan Zebra_Finch Gibbon Panda Zebrafish);
    my @speciesListLatin = qw(Vicugna_pacos Gorilla_gorilla Sus_scrofa Anolis_carolinensis Cavia_porcellus Ochotona_princeps Dasypus_novemcinctus Erinaceus_europaeus Ornithorhynchus_anatinus Equus_caballus Oryctolagus_cuniculus Otolemur_garnettii Homo_sapiens Rattus_norvegicus Caenorhabditis_elegans Procavia_capensis Saccharomyces_cerevisiae Ciona_intestinalis Dipodomys_ordii Ciona_savignyi Sorex_araneus Felis_catus Echinops_telfairi Choloepus_hoffmanni Gallus_gallus Macaca_mulatta Spermophilus_tridecemlineatus Pan_troglodytes Callithrix_jacchus Gasterosteus_aculeatus Bos_taurus Oryzias_latipes Tarsius_syrichta Canis_familiaris Pteropus_vampyrus Tetraodon_nigroviridis Tursiops_truncatus Myotis_lucifugus Tupaia_belangeri Mus_musculus Meleagris_gallopavo Loxodonta_africana Microcebus_murinus Macropus_eugenii Drosophila_melanogaster Monodelphis_domestica Xenopus_tropicalis Takifugu_rubripes Pongo_abelii Taeniopygia_guttata Nomascus_leucogenys Ailuropoda_melanoleuca Danio_rerio);
    my $html = get("http://uswest.ensembl.org/info/about/species.html");
    #my $html = get("foo");
    if (defined $html){
      $html =~ m{Ensembl [Ss]pecies[</a-z0-9>\s]*<table(.*?)</table>}s;
      if (defined $1){
          @speciesListCommon = ();
          @speciesListLatin = ();
          my $speciestable = $1;
          while ($speciestable =~ m{<a href=.*?>(.*?)</a>[)]?<br />(.*?)</td>}g) {
            if (!($1 =~ m{preview })) { 	# skip species that are only preview
              my $commonname = $1;
              $commonname =~ tr[ ][_];
              my $latinname;
              my $temp = $2;
              if ($temp =~ m{<i>(.*?)</i>}) {
                $latinname = $1;
                $latinname =~ tr[ ][_];
                if ($commonname eq $latinname) {
                  $commonname =~ m{([a-zA-Z]*?)_([a-z]*)};
                  $commonname = substr($1, 0, 1).".".$2;
                }
              }
              else {
                $latinname = $commonname;
                $commonname =~ m{([a-zA-Z]*?)_([a-z]*)};
                $commonname = substr($1, 0, 1).".".$2;
              }
              push(@speciesListCommon, $commonname);
              push(@speciesListLatin, $latinname);
            }		
          }
      }
    }
    #print '"'.join('", "',@speciesListCommon).'"'."\n";
    #print '"'.join('", "',@speciesListLatin).'"'."\n";
    return \(@speciesListCommon,@speciesListLatin);
}
#ZILA

#TPAIGE
sub get_uniprot_annotations {
    my $upid = shift;
    
    
    my $record = get_uniprot_record(
        script => $config->{uniprot_script},
        id     => $upid,
        db     => 'uniprotkb',
        format => 'SWISS',
        style  => 'raw'
    );
    
    local $/ = "\n//\n";
    my $fullParse=1;
    
    my %topics;

    my  $entry = SWISS::Entry->fromText($record, $fullParse);
    
    #description
    foreach my $de ($entry->DEs->elements) {
        push @{$topics{DE}}, $de->text;
    }
    #protein existence
    if (defined($entry->PE->text) && length(trim($entry->PE->text))>0){
        push @{$topics{PE}}, trim($entry->PE->text);
    }

    foreach my $kw ($entry->KWs->elements) {
        push @{$topics{KW}}, $kw->text;
    }
    for my $CC ($entry->CCs->elements()) {
        if (grep {$CC->topic eq $_} @{$config->{exclude_uniprot_cc}}){
            next;
        }
        if ($CC->topic eq 'RNA EDITING') {
            push @{$topics{$CC->topic}}, $CC->note;
        }
        else{
            push @{$topics{$CC->topic}}, $CC->comment;
        }
    }

  return \%topics;
}
#TPAIGE
