=head1 LICENSE

 Copyright (c) 2011-2014 Michael Yourshaw.  All rights reserved.                                                                      

=head1 CONTACT                                                                                                       

 Michael Yourshaw <myourshaw@ucla.edu>
    
=cut

=head1 NAME

 Conservation

=head1 SYNOPSIS

 mv OMES.pm ~/.vep/Plugins
 perl variant_effect_predictor.pl -i variations.vcf --plugin OMES

=head1 DESCRIPTION

 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that
 adds the following new columns:
 Protein_OMES_Orthologues, Protein_OMES_Family.
 
 WARNING: This plugin is untested.
 
 =head1 PARAMETERS
    parameters are in the form key=value, comma-separated
    'compara_db=s', #the name of the compara database to be used, default 'Multi'
    'omes_script=s', #the location of the OMES_score.pl script
    'omes_family', #whether or not omes should be calculated using protein family
    'omes_orthologues', #whether or not omes should be calculated using orthologues
    'omes_minimum_seqs=i', #the minimum number of sequences needed to calculate OMES score, after filtering by identity, default 10
    'omes_maximum_identity=i', #the maximum pairwise percent identity allowed between two sequences used in OMES calculation, default 90
    'omes_max_seqs=i', #the maximum number of sequences to be submitted to OMES calculation. Set to undef to use all, default 40

=cut

package OMES;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);
use Bio::EnsEMBL::Variation::Utils::VEP qw(@OUTPUT_COLS);
use File::Basename;

use VAX qw(get_bvfoa_info);

sub new {
    my $class = shift;
    my $self = $class->SUPER::new(@_);
    my $config = $self->{config};
    my %params;
    foreach my $param(@{$self->{params}}){
        my ($key,$value) = split(/=/,$param,1);
        $params{$key} = $value;
    }
    $params{compara_db} ||= 'Multi';
    $params{omes_script} ||= File::Spec->catfile(dirname(__FILE__), 'OMES_score.pl');
    $params{omes_family} ||= 1;
    $params{omes_orthologues} ||= 1;
    $params{omes_minimum_seqs} ||= 10;
    $params{omes_maximum_identity} ||= 90;
    $params{omes_max_seqs} ||= 40;
    $self->{params} = \%params;

    #adaptors used by NGS-SNP
    $self->{ma} = $config->{reg}->get_adaptor($params{compara_db}, 'compara', 'Member');
    $self->{ha} = $config->{reg}->get_adaptor($params{compara_db}, 'compara', 'Homology');
    $self->{fa} = $config->{reg}->get_adaptor($params{compara_db}, 'compara', 'Family');
    
    return $self;
}

sub version {
    return '75';
}

sub feature_types {
    return ['Bio::EnsEMBL::Transcript'];
}

sub get_header_info {
    my @new_output_cols = qw(
        Protein_OMES_Orthologues
        Protein_OMES_Family
    );
    @OUTPUT_COLS = (@OUTPUT_COLS, @new_output_cols);

    return {
        Protein_OMES_Orthologues => "reports the residues in the reference sequence that are most mutationally correlated with the altered residue, as determined using the Observed Minus Expected Squared (OMES) method described in Fodor and Aldrich (2004 Proteins 56(2): 211-221). Orthologous protein sequences are obtained from Ensembl and aligned with the reference using Muscle. All column pairs in the resulting multiple alignment are then compared to identify correlated residues. The output for each correlated residue is in the form 'score(X)(Y)' where 'score' is the raw OMES score, 'X' is the position of the SNP-altered residue in the reference sequence, and 'Y' is the position of the correlated residue in the reference sequence. 'Y' may sometimes be given as '-', which indicates that the correlated position detected in the multiple alignment is absent from the reference sequence. The value for this key is only calculated if the -run_omes option is specified",
        Protein_OMES_Family => "determined in the same manner as the 'Protein_OMES_Orthologues' value, except that non-orthologous sequences belonging to the same protein family as the reference sequence are also used in the calculation. The value for this key is only calculated if the -run_omes option is specified",
    };
}

sub run {
    my ($self, $bvfoa, $line_hash) = @_;
    my $bvf = $bvfoa->base_variation_feature;    my $config = $self->{config};
    
    my %bvfoa_info = %{get_bvfoa_info(@_)};

    if(defined($bvfoa_info{translation})
       && defined($bvfoa_info{ensg})
       && defined($bvfoa_info{ensp})
       && defined($bvfoa_info{amino_acid_reference})
       && defined($bvfoa_info{altered_aa_start})
       ){
        
        #OMES
        my ($correlated_family, $correlated_orthologues);
         ($correlated_family, $correlated_orthologues) =
            determine_if_correlated_amino_acid(
                $self,
                $bvfoa_info{ensp},
                $bvfoa_info{changes_protein},
                $bvfoa_info{altered_aa_start},
                $bvfoa_info{amino_acid_reference},
                $bvfoa_info{ensg}
                );
        #Protein_OMES_Orthologues
        if (defined($correlated_orthologues)) {
            my $protein_omes_orthologues = $correlated_orthologues;
            $line_hash->{Protein_OMES_Orthologues} = $protein_omes_orthologues;
        }
        #Protein_OMES_Family
        if (defined($correlated_family)) {
            my $protein_omes_family = $correlated_family;
            $line_hash->{Protein_OMES_Family} = $protein_omes_family;
        }
     } #if(defined($translation))

    return {};
}

#whether the affected amino acid position is correlated with other amino acid
#positions, according to the OMES procedure described in
#Fodor and Aldrich PROTEINS: Structure, Function, and Bioinformatics 56:211–221
#and
#Andreas Kowarsch et al, PLOS Computational Biology, 6:e1000923
sub determine_if_correlated_amino_acid {
    
    my ($self, $protein_id, $changes_protein, $altered_aa_start, $amino_acid_reference, $gene_id) = @_;
    
    #return values
    my ($correlated_family, $correlated_orthologues);

    if ((!defined($protein_id))
        || (!defined($changes_protein))
        || (!defined( $self->{params}->{omes_script})))
    {
        return;
    }

    #A multiple protein alignment is needed for the calculation used to determine
    #correlated mutations

    if ($self->{params}->{omes_family}) {
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
            $self,
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

                if ((defined($self->{params}->{omes_max_seqs}))
                    && (scalar(@seqs_for_omes) > $self->{params}->{omes_max_seqs}))
                {
                    @seqs_for_omes =
                        @seqs_for_omes[0 .. ($self->{params}->{omes_max_seqs} - 1)];
                }

                my $omes_result = calculate_omes(
                    input_seqs           => \@seqs_for_omes,
                    min_seqs             => $self->{params}->{omes_minimum_seqs},
                    max_ident            => $self->{params}->{omes_maximum_identity},
                    position_of_interest => $altered_aa_start,
                    omes_script          => $self->{params}->{omes_script},
                    protein_id           => $protein_id,
                    comparison_type      => 'family'
                );

                $correlated_family = $omes_result;
            }
        }
    }

    #######################
    if ($self->{params}->{omes_orthologues}) {

        #######################
        #Option 2: Use a list of Bio::EnsEMBL::Compara::Homology objects.
        #Homology objects store orthologous and paralogous relationships between Members
        #This approach can be used to obtain orthologues and their sequences.
        #Drawback is that the sequences will need to be aligned
        #######################
        my @seqs_for_omes = ();

        my $results = get_protein_orthologues_for_omes_from_homology(
            $self,
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

                if ((defined($self->{params}->{omes_max_seqs}))
                    && (scalar(@seqs_for_omes) > $self->{params}->{omes_max_seqs})
                    )
                {
                    @seqs_for_omes =
                        @seqs_for_omes[0 .. ($self->{params}->{omes_max_seqs} - 1)];
                }

                my $omes_result = calculate_omes(
                    input_seqs           => \@seqs_for_omes,
                    min_seqs             => $self->{params}->{omes_minimum_seqs},
                    max_ident            => $self->{params}->{omes_maximum_identity},
                    position_of_interest => $altered_aa_start,
                    omes_script          => $self->{params}->{omes_script},
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
    
    my ($self, $protein_id, $altered_aa_start, $amino_acid_reference, $gene_id) = @_;
    
    my %results = (
        seqs                            => [],
        altered_protein_id              => undef,
        altered_aa_position_in_sequence => undef,
        altered_aa                      => undef
    );

    $results{altered_protein_id} = $protein_id;
    $results{altered_aa_position_in_sequence} = $altered_aa_start;
    $results{altered_aa} = $amino_acid_reference;

    my $member = $self->{ma}
        ->fetch_by_source_stable_id('ENSEMBLGENE', $gene_id);

    #Rarely the $member object may be undef
    if (!defined($member)) {
        return;
    }

    my $homologies
        = $self->{ha}->fetch_all_by_Member_method_link_type($member,
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
    
    my ($self, $protein_id, $altered_aa_start, $amino_acid_reference, $gene_id) = @_;
    
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

    my $member = $self->{ma}
        ->fetch_by_source_stable_id('ENSEMBLGENE', $gene_id);

    #Rarely the $member object may be undef
    if ( !defined($member) ) {
        return;
    }

    my $families = $self->{fa}->fetch_all_by_Member($member);

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


1;

