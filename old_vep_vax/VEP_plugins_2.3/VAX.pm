=head1 LICENSE

 Copyright (c) 2011-2012 Michael Yourshaw.  All rights reserved.

 For license, please contact

   myourshaw@ucla.edu

=head1 CONTACT

 Please email comments or questions to myourshaw@ucla.edu
 
=cut

=head1 NAME

VAX - Methods used by the Ensembl Variant Effect Predictor VAX plugins

=head1 SYNOPSIS

  use VAX qw([get_unique] [...]);

=head1 METHODS

=cut


package VAX;

use strict;
use warnings;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(
    &get_tva_info
    &get_consequence_info
    &get_unique
    &replace_str
    &trim
    &ltrim
    &rtrim
    &strip
    &lstrip
    &rstrip
    &overlaps
    &get_base
    &replace_base
    &space_to_underscore
    &complement
    &reverse_complement
    &get_taxon_name
    &is_gap
);

sub get_tva_info{
    my ($self, $tva, $line_hash) = @_;
    my $config = $self->{config};
    
    my %tva_info;
    
    if ($tva->isa('Bio::EnsEMBL::Variation::TranscriptVariationAllele')){
        $tva_info{variation} = $tva->transcript_variation;
        $tva_info{feature} = $tva->variation_feature;
        $tva_info{transcript} = $tva->transcript;
        $tva_info{enst} = $tva_info{transcript}->stable_id;
        $tva_info{gene} = $config->{ga}->fetch_by_transcript_stable_id($tva_info{enst});
        $tva_info{ensg} = $tva_info{transcript}->{_gene_stable_id};
        $tva_info{hgnc} = $tva_info{transcript}->{_gene_hgnc};
        $tva_info{translation} = $tva_info{transcript}->translation;
        if(defined($tva_info{translation})){
            $tva_info{ensp} = $tva_info{translation}->{stable_id};
            $tva_info{protein_sequence} = $tva_info{translation}->{seq};
            $tva_info{protein_length} = length($tva_info{protein_sequence});
            $tva_info{altered_aa_start} = $tva_info{variation}->{translation_start};
            $tva_info{altered_aa_end} = $tva_info{variation}->{translation_end};
            $tva_info{altered_base_start} = $tva_info{variation}->{cdna_start};
            $tva_info{altered_base_end} = $tva_info{variation}->{cdna_end};
            $tva_info{pep_allele_string} = $tva->pep_allele_string;
            if (defined($tva_info{pep_allele_string})) {
                #reference is first
                if ($tva_info{pep_allele_string} =~ m/([^\/]+)\/([^\/]+)/) {
                    $tva_info{amino_acid_reference} = $1;
                    $tva_info{amino_acid_variant} = $2;
                    $tva_info{changes_protein} = ($tva_info{amino_acid_reference} ne $tva_info{amino_acid_variant}) ? 1 : 0;
                }
            }
        }
    }
    elsif($tva->isa('Bio::EnsEMBL::Variation::RegulatoryFeatureVariationAllele')){
        $tva_info{variation} = $tva->regulatory_feature_variation;
        $tva_info{feature} = $tva->regulatory_feature;
        $tva_info{nearest_genes} = $config->{ga}->fetch_nearest_Gene_by_Feature($tva_info{feature});
        $tva_info{gene} = @{$tva_info{nearest_genes}}[0];
    }
    elsif($tva->isa('Bio::EnsEMBL::Variation::MotifFeatureVariationAllele')){
        $tva_info{variation} = $tva->motif_feature_variation;
        $tva_info{feature} = $tva->motif_feature;
        $tva_info{nearest_genes} = $config->{ga}->fetch_nearest_Gene_by_Feature($tva_info{feature});
        $tva_info{gene} = @{$tva_info{nearest_genes}}[0];
    }
    else{
        warn "unrecognized tva type $tva->type";
    }
    $tva_info{base_variation_feature} = $tva_info{variation}->{base_variation_feature};
    $tva_info{_line} = $tva_info{base_variation_feature}->{_line};
    if(!defined($tva_info{ensg})) {
        $tva_info{ensg} = $tva_info{gene} ? $tva_info{gene}->stable_id : undef;
    }
    if(!defined($tva_info{hgnc} && $tva_info{gene})) {
        my @entries = grep {$_->database eq 'HGNC'} @{$tva_info{gene}->get_all_DBEntries()};
        if(scalar @entries) {
            $tva_info{hgnc} = $entries[0]->display_id;
        }
    }
    $tva_info{hgnc} = undef if defined($tva_info{hgnc}) && ($tva_info{hgnc} eq '' || $tva_info{hgnc} eq '-');
    
    $tva_info{chrom} = $tva_info{feature}->seq_region_name;
    $tva_info{chrom_start} = $tva_info{feature}->start;
    $tva_info{chrom_end} = $tva_info{feature}->end;
    $tva_info{chrom_strand} = $tva_info{feature}->strand;
    if(defined($tva_info{chrom})
    && defined($tva_info{chrom_start})
    && defined($tva_info{chrom_end})
    && defined($tva_info{chrom_strand})){
        my %genomic_coords = (
            chr    => $tva_info{chrom},
            start  => $tva_info{chrom_start},
            end    => $tva_info{chrom_end},
            strand => $tva_info{chrom_strand}
        );
        $tva_info{genomic_coords} = \%genomic_coords;
    }
    $tva_info{reference_allele} = $tva_info{variation}->{reference_allele}->variation_feature_seq;
    $tva_info{non_reference_allele} = $tva->variation_feature_seq;
   
    return \%tva_info;
} #get_tva_info

sub get_consequence_info{
    my ($self, $tva, $line_hash) = @_;
    my $config = $self->{config};
    
    my %consequences_info;
    
    my $term_method = $config->{terms}.'_term';
    my @c_so = map {$_->SO_term} @{$tva->get_all_OverlapConsequences};
    my @consequences = map {$_->$term_method || $_->SO_term} @{$tva->get_all_OverlapConsequences};
    my @consequences_so = map {$_->SO_term || $_->SO_term} @{$tva->get_all_OverlapConsequences};
    my @consequences_ensembl = map {$_->display_term || $_->SO_term} @{$tva->get_all_OverlapConsequences};
    my @consequences_ncbi = map {$_->NCBI_term || $_->SO_term} @{$tva->get_all_OverlapConsequences};
    my @ranks = map {$_->rank} @{$tva->get_all_OverlapConsequences};
    my %consequences_ranks;
    my %consequences_ranks_so;
    my %consequences_ranks_ensembl;
    my %consequences_ranks_ncbi;
    map {$consequences_ranks{$consequences[$_]} = $ranks[$_]} 0..$#ranks;
    map {$consequences_ranks_so{$consequences_so[$_]} = $ranks[$_]} 0..$#ranks;
    map {$consequences_ranks_ensembl{$consequences_ensembl[$_]} = $ranks[$_]} 0..$#ranks;
    map {$consequences_ranks_ncbi{$consequences_ncbi[$_]} = $ranks[$_]} 0..$#ranks;
    my @consequences_sorted = sort {$consequences_ranks{$a} <=> $consequences_ranks{$b}} keys %consequences_ranks;
    my @consequences_so_sorted = sort {$consequences_ranks_so{$a} <=> $consequences_ranks_so{$b}} keys %consequences_ranks_so;
    my @consequences_ensembl_sorted = sort {$consequences_ranks_ensembl{$a} <=> $consequences_ranks_ensembl{$b}} keys %consequences_ranks_ensembl;
    my @consequences_ncbi_sorted = sort {$consequences_ranks_ncbi{$a} <=> $consequences_ranks_ncbi{$b}} keys %consequences_ranks_ncbi;
    
    $consequences_info{consequences_sorted} = \@consequences_sorted;
    $consequences_info{consequence_severest} = $consequences_sorted[0];
    $consequences_info{consequence_severest_rank} = $consequences_ranks{$consequences_info{consequence_severest}};
    $consequences_info{consequences_ranks} = \%consequences_ranks;
    
    $consequences_info{consequences_so_sorted} = \@consequences_so_sorted;
    $consequences_info{consequence_so_severest} = $consequences_so_sorted[0];
    $consequences_info{consequence_so_severest_rank} = $consequences_ranks_so{$consequences_info{consequence_so_severest}};
    $consequences_info{consequences_ranks_so} = \%consequences_ranks_so;
    
    $consequences_info{consequences_ensembl_sorted} = \@consequences_ensembl_sorted;
    $consequences_info{consequence_ensembl_severest} = $consequences_ensembl_sorted[0];
    $consequences_info{consequence_ensembl_severest_rank} = $consequences_ranks_ensembl{$consequences_info{consequence_ensembl_severest}};
    $consequences_info{consequences_ranks_ensembl} = \%consequences_ranks_ensembl;
    
    $consequences_info{consequences_ncbi_sorted} = \@consequences_ncbi_sorted;
    $consequences_info{consequence_ncbi_severest} = $consequences_ncbi_sorted[0];
    $consequences_info{consequence_ncbi_severest_rank} = $consequences_ranks_so{$consequences_info{consequence_ncbi_severest}};
    $consequences_info{consequences_ranks_ncbi} = \%consequences_ranks_ncbi;
    
    return \%consequences_info;
}

sub get_unique{
    my $list = shift;
    my %seen = ();
    my @uniq = grep { !$seen{$_}++ } @{$list};
    return \@uniq;
}

sub replace_str{
    #uses zero-based indexing
    my ($string, $replacement, $start, $end) = @_;
    my $a = substr($string,0,$start);
    my $z = substr($string,$end+1);
    my $replaced = $a . $replacement . $z;
    return $a . $replacement . $z;
}
# Perl trim function to remove whitespace from the start and end of the string
sub trim($){
	my $string = shift;
	$string =~ s/^\s+|\s+$//g;
	return $string;
}
# Left trim function to remove leading whitespace
sub ltrim($){
	my $string = shift;
	$string =~ s/^\s+//g;
	return $string;
}
# Right trim function to remove trailing whitespace
sub rtrim($){
	my $string = shift;
	$string =~ s/\s+$//g;
	return $string;
}

# Perl strip function to remove regex from the start and end of the string
sub strip($$){
	my ($string,$strip) = @_;
	$string =~ s/^$strip+|$strip+$//g;
	return $string;
}
# Left strip function to remove leading regex
sub lstrip($$){
	my ($string,$strip) = @_;
	$string =~ s/^$strip+//g;
	return $string;
}
# Right strip function to remove trailing regex
sub rstrip($$){
	my ($string,$strip) = @_;
	$string =~ s/$strip+$//g;
	return $string;
}

sub overlaps{
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

sub get_base{

    #first base is position '1'
    my $position = shift;
    my $sequence = shift;
    return substr( $sequence, $position - 1, 1 );
}

sub replace_base{

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

sub space_to_underscore{
    my $text = shift;
    $text =~ s/\s{1,}/_/g;
    return $text;
}

sub complement {
    my $seq = shift;
    $seq =~ tr/gatcryswkmbdhvnGATCRYSWKMBDHVN/ctagyrswmkvhdbnCTAGYRSWMKVHDBN/;
    return $seq;
}
sub reverse_complement {
    my $seq = shift;
    $seq = reverse complement($seq);
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

sub is_gap {
    my $residue = shift;
    if ( $residue =~ m/[\.\-]/ ) {
        return 1;
    }
    return 0;
}


1;
