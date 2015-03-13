#FILE: OMES_score.pl
#AUTH: Xiaoping Liao (xliao2@cs.ualberta.ca)
#AUTH: Paul Stothard (stothard@ualberta.ca)
#DATE: October 3, 2010
#VERS: 1.0

use warnings;
use strict;
use File::Temp;
use File::Basename;
use Getopt::Long;
use Data::Dumper;
use DB_File;

#sort keys by third value in key name
#for example 3.009 in 1_20_3.009.
#then break ties by first value
#then break ties by second value
$DB_BTREE->{'compare'} = sub {
    my ( $key1, $key2 ) = @_;
    my @key1 = split( /_/, $key1 );
    my @key2 = split( /_/, $key2 );
    $key2[2] <=> $key1[2] || $key1[0] <=> $key2[0] || $key1[1] <=> $key2[1];
};

my %options = (
    input                       => undef,
    output                      => undef,
    multiple_alignment_file     => undef,
    number_of_results_to_report => undef,
    percent_identity_threshold  => 90,
    minimum_sequence_number     => 30,
    stretcher_command           => 'stretcher',
    t_coffee_command            => 't_coffee',
    muscle_command              => 'muscle',
    use_t_coffee                => undef,
    residue                     => undef,
    verbose                     => undef,
    help                        => undef
);

GetOptions(
    'i=s'          => \$options{input},
    'o=s'          => \$options{output},
    'a=s'          => \$options{multiple_alignment_file},
    'n=i'          => \$options{number_of_results_to_report},
    'p=i'          => \$options{percent_identity_threshold},
    'm=i'          => \$options{minimum_sequence_number},
    'stretcher=s'  => \$options{stretcher_command},
    't_coffee=s'   => \$options{t_coffee_command},
    'muscle=s'     => \$options{muscle_command},
    'use_t_coffee' => \$options{use_t_coffee},
    'residue=i'    => \$options{residue},
    'v'            => \$options{verbose},
    'h|help'       => \$options{help}
);

if ( defined( $options{help} ) ) {
    print_usage();
    exit(0);
}

# check the options
check_option(
    option_name     => 'i',
    option_value    => $options{input},
    must_be_defined => 1
);
check_option(
    option_name     => 'o',
    option_value    => $options{output},
    must_be_defined => 1
);
check_option(
    option_name     => 'a',
    option_value    => $options{multiple_alignment_file},
    must_be_defined => 0
);
check_option(
    option_name  => 'n',
    option_value => $options{number_of_results_to_report},
    type_int     => 1,
    min_value    => 1
);
check_option(
    option_name  => 'p',
    option_value => $options{percent_identity_threshold},
    type_int     => 1,
    min_value    => 1
);
check_option(
    option_name  => 'm',
    option_value => $options{minimum_sequence_number},
    type_int     => 1,
    min_value    => 1
);
check_option(
    option_name     => 'stretcher',
    option_value    => $options{stretcher_command},
    must_be_defined => 1
);
check_option(
    option_name     => 't_coffee',
    option_value    => $options{t_coffee_command},
    must_be_defined => 1
);
check_option(
    option_name     => 'muscle',
    option_value    => $options{muscle_command},
    must_be_defined => 1
);
check_option(
    option_name     => 'residue',
    option_value    => $options{residue},
    must_be_defined => 0,
    type_int        => 1,
    min_value       => 1
);

# read the input fasta format data into parallel arrays
# @file_data stores the sequences data
# @description_data stores the description of each sequence
my @file_data        = extract_sequence_from_fasta_data( $options{input} );
my @description_data = extract_description_from_fasta_data( $options{input} );

# call stretcher program to calculate the pairwise global alignment percent identity
my @pairwise_alignment_identity = ();

message( $options{verbose}, "Performing pairwise alignments.\n" );

for my $i ( 0 .. ( @file_data - 1 ) ) {
    for my $j ( ( $i + 1 ) .. ( @file_data - 1 ) ) {

        message( $options{verbose},
            "Aligning '$description_data[$i]' and '$description_data[$j]'.\n"
        );

        my $tmp_seqa = new File::Temp()->filename;
        my $tmp_seqb = new File::Temp()->filename;
        write_to_temp_seq( $file_data[$i], $tmp_seqa );
        write_to_temp_seq( $file_data[$j], $tmp_seqb );
        my $tmp_output = new File::Temp()->filename;

        my $cmd
            = "$options{stretcher_command} $tmp_seqa $tmp_seqb $tmp_output -auto y -sprotein1 -sprotein2";
        my $result = system($cmd);
        if ( $result != 0 ) {
            die("The following command failed: '$cmd'\n");
        }

        open( my $TEMP_FILE, '<', $tmp_output )
            or die("Cannot open file '$tmp_output': $!");
        while (<$TEMP_FILE>) {
            chomp $_;

            my $identity = $1 if /Identity:\s+(\d+\/\d+).*/;
            $pairwise_alignment_identity[$i][$j] = $identity if $identity;
        }
        close($TEMP_FILE);
    }
}

# make it a symmetric matrix
for my $i ( 0 .. ( @file_data - 1 ) ) {
    $pairwise_alignment_identity[$i][$i] = 0;
    for my $j ( ( $i + 1 ) .. ( @file_data - 1 ) ) {
        $pairwise_alignment_identity[$j][$i]
            = $pairwise_alignment_identity[$i][$j];
    }
}

# format the data
for my $i ( 0 .. ( @file_data - 1 ) ) {
    for my $j ( 0 .. ( @file_data - 1 ) ) {
        $pairwise_alignment_identity[$i][$j]
            = sprintf( "%4.3f", eval( $pairwise_alignment_identity[$i][$j] ) )
            + 0;
    }
}

message( $options{verbose},
    "Removing sequences that are more than $options{percent_identity_threshold} percent identical.\n"
);

# find the index of sequence that should be deleted
my $index = find_index( @pairwise_alignment_identity,
    $options{percent_identity_threshold} );
my @delete_seq_label = ();
while ( $index != -1 ) {
    push( @delete_seq_label, $index );

    for my $j ( 0 .. ( @file_data - 1 ) ) {
        $pairwise_alignment_identity[$index][$j] = 0;
    }
    for my $i ( 0 .. ( @file_data - 1 ) ) {
        $pairwise_alignment_identity[$i][$index] = 0;
    }
    $index = find_index( @pairwise_alignment_identity,
        $options{percent_identity_threshold} );
}

my @seq_label            = ( 0 .. ( @file_data - 1 ) );
my $seq_label_ref        = \@seq_label;
my $delete_seq_label_ref = \@delete_seq_label;

message( $options{verbose},
    scalar(@delete_seq_label) . " sequences were removed.\n" );

# @label stores the index of sequences to be kept
my @label = set_difference( $seq_label_ref, $delete_seq_label_ref );

# after filtering, temp file $final_fasta stores the final seq data for t_coffee program

my $final_fasta = new File::Temp()->filename;
open( my $TEMP_FASTA, '>', $final_fasta )
    or die("Cannot open file '$final_fasta': $!");
for my $i ( 0 .. ( @label - 1 ) ) {
    print $TEMP_FASTA '>' . $description_data[ $label[$i] ] . "\n";
    print $TEMP_FASTA $file_data[ $label[$i] ] . "\n";
}
close($TEMP_FASTA);

# if the number of seq in the $final_fasta is smaller than the default value,
# writes '#minimum_sequence_number not reached' to output file and exit
open( my $OUTPUT, ">" . "$options{output}" )
    || die "cannot open $options{output} :$!";
if ( @label < $options{minimum_sequence_number} ) {
    print $OUTPUT '#minimum_sequence_number not reached';
    message( $options{verbose},
        "Fewer than $options{minimum_sequence_number} sequences remain--exiting.\n"
    );
    exit(0);
}

message( $options{verbose},
    "Generating multiple alignment for remaining sequences.\n" );

# call t_coffee to get the alignment file
my $tmp_output   = new File::Temp()->filename;
my $tmp_treefile = new File::Temp()->filename;
my $cmd;
if ( $options{use_t_coffee} ) {
    $cmd
        = "$options{t_coffee_command} $final_fasta -outfile=$tmp_output -output=fasta -outorder=input -newtree=$tmp_treefile > /dev/null 2>&1";
}
else {

#Note that the Muscle program used to support a '-stable' option, which caused
#sequence order to be preserved. This option is no longer supported, and in
#the versions that do support it the ordering is said to be unreliable.

  #Note that Muscle does not recognize 'U', even though it claims to in its
  #documentation. Unrecognized characters are replaced with 'X' in the output.

#for accuracy
#$cmd = "$options{muscle_command} -seqtype protein -in $final_fasta -out $tmp_output -quiet";

    #for speed
    $cmd
        = "$options{muscle_command} -seqtype protein -in $final_fasta -out $tmp_output -quiet -maxiters 1 -diags -sv -distance1 kbit20_3";
}
my $result = system($cmd);
if ( $result != 0 ) {
    die("The following command failed: '$cmd'\n");
}

# if the Muscle program is used, need to sort the sequences so that the correct sequence comes first.
if ( !$options{use_t_coffee} ) {
    my %hash_seq            = ();
    my @file_data_alignment = extract_sequence_from_fasta_data($tmp_output);
    my @description_data_alignment
        = extract_description_from_fasta_data($tmp_output);
    for my $i ( 0 .. ( @file_data_alignment - 1 ) ) {
        $hash_seq{ $description_data_alignment[$i] }
            = $file_data_alignment[$i];
    }
    open( my $TEMP_FASTA, '>', $tmp_output )
        or die("Cannot open file '$tmp_output': $!");
    for my $i ( 0 .. ( @label - 1 ) ) {

        my $title    = $description_data[ $label[$i] ];
        my $sequence = $hash_seq{ $description_data[ $label[$i] ] };

        $title    =~ s/\s+$//;
        $sequence =~ s/\s+$//;

        print $TEMP_FASTA '>' . $title . "\n";
        print $TEMP_FASTA $sequence . "\n";
    }
    close($TEMP_FASTA);
}

# if options multiple_alignment_file is specified, output the alignment file.
if ( $options{multiple_alignment_file} ) {

    message( $options{verbose},
        "Multiple alignment will be written to $options{multiple_alignment_file}.\n"
    );

    open( my $TEMP, $tmp_output )
        or die("Cannot open file '$tmp_output': $!");
    open( my $ALIGNMENT_FILE, '>', $options{multiple_alignment_file} )
        or die("Cannot open file $options{multiple_alignment_file}: $!");
    while (<$TEMP>) { print $ALIGNMENT_FILE $_; }
    close($TEMP);
    close($ALIGNMENT_FILE);
}

#@final_file_data stores the final sequence data
my @final_file_data = extract_sequence_from_fasta_data("$tmp_output");

# initial variables
my $number_protein = 0;
my $length_protein = 0;

# read the first line to store some additional information for output
my $first_line = $final_file_data[0];

$length_protein = length( $final_file_data[0] );
$number_protein = scalar(@final_file_data);
my @element = split( //, $first_line );
my $count   = 0;
my %info    = ();
for my $i ( 0 .. ( @element - 1 ) ) {
    $info{ $i + 1 } = '-' if $element[$i] eq '-';
    if ( $element[$i] ne '-' ) {
        $count++;
        $info{ $i + 1 } = $count;
    }

}

#Caculate the OMES score, and use a hash table score_hash to store the final score
my $tmp_hash_file = new File::Temp()->filename;
tie(my %score_hash,
    'DB_File', $tmp_hash_file, O_RDWR | O_CREAT,
    0640, $DB_BTREE
) || die $!;

if ( defined( $options{residue} ) ) {

    #determine position of residue in aligned sequence
    my $residue_location_in_alignment = undef;
    foreach my $key ( keys(%info) ) {
        if ( $info{$key} eq $options{residue} ) {
            $residue_location_in_alignment = $key;
            last;
        }
    }

    if ( !defined($residue_location_in_alignment) ) {

        #This happens when the altered residue was a stop codon in reference
        #and thus is not in the alignment
        print $OUTPUT '#desired residue not present in alignment';
        message( $options{verbose},
            "Residue $options{residue} is not in the alignment--exiting.\n" );
        untie %score_hash;
        exit(0);
    }

    message( $options{verbose},
        "Calculating OMES score for residue $options{residue} from multiple alignment of length $length_protein and consisting of $number_protein sequences.\n"
    );

    my $i = $residue_location_in_alignment;

    message( $options{verbose}, "processing column $i.\n" );

    for my $j ( 1 .. ( $i - 1 ), ( $i + 1 ) .. $length_protein ) {

        my $total_score = omes( $i, $j, \@final_file_data );
        $score_hash{ $i . '_' . $j . '_' . $total_score } = 1;

    }

}
else {
    message( $options{verbose},
        "Calculating OMES scores from multiple alignment of length $length_protein and consisting of $number_protein sequences.\n"
    );

    for my $i ( 1 .. $length_protein ) {
        message( $options{verbose}, "processing column $i.\n" );

        for my $j ( $i + 1 .. $length_protein ) {

            my $total_score = omes( $i, $j, \@final_file_data );
            $score_hash{ $i . '_' . $j . '_' . $total_score } = 1;

        }
    }

}

#print results
if ( $options{number_of_results_to_report} ) {
    message( $options{verbose},
        "Choosing the top $options{number_of_results_to_report} scores to report.\n"
    );
    my $written_count = 0;
    my $previous_value;
    while ( my ( $key, $value ) = each %score_hash ) {

        my @split_value = split( /_/, $key );

        #write out all results tied for last place
        if ( $written_count >= $options{number_of_results_to_report} ) {
            unless ( ( defined($previous_value) )
                && ( $previous_value == $split_value[2] ) )
            {
                last;
            }
        }

        print $OUTPUT $split_value[2] . ","
            . $split_value[0]
            . "($info{$split_value[0]})" . ','
            . $split_value[1]
            . "($info{$split_value[1]})" . "\n";

        $written_count++;
        $previous_value = $split_value[2];
    }

}
else {
    while ( my ( $key, $value ) = each %score_hash ) {
        my @split_value = split( /_/, $key );

        print $OUTPUT $split_value[2] . ","
            . $split_value[0]
            . "($info{$split_value[0]})" . ','
            . $split_value[1]
            . "($info{$split_value[1]})" . "\n";
    }
}

untie %score_hash;

close($OUTPUT);

message( $options{verbose}, "Results written to $options{output}.\n" );

sub check_option {
    my %args = (@_);

    if ( $args{must_be_defined} ) {
        if ( !defined( $args{option_value} ) ) {
            print "Option '$args{option_name}' must be defined.\n";
            print_usage();
            exit(1);
        }
    }

    if ( !defined( $args{option_value} ) ) {
        return;
    }

    if ( $args{type_int} ) {
        if ( $args{option_value} =~ m/[^\d\+\-]/ ) {
            die("Option '$args{option_name}' must be an integer.");
        }
    }

    if ( $args{type_real} ) {
        if ( $args{option_value} =~ m/[^\d\+\-\.]/ ) {
            die("Option '$args{option_name}' must be a real number.");
        }
    }

    if ( defined( $args{max_value} ) ) {
        if ( $args{option_value} > $args{max_value} ) {
            die("Option '$args{option_name}' must be less than or equal to $args{max_value}."
            );
        }
    }

    if ( defined( $args{min_value} ) ) {
        if ( $args{option_value} < $args{min_value} ) {
            die("Option '$args{option_name}' must be greater than or equal to $args{min_value}."
            );
        }
    }

    if ( defined( $args{allowed_values} ) ) {
        my $is_allowed = 0;
        foreach my $allowed_value ( @{ $args{allowed_values} } ) {
            if ( lc( $args{option_value} ) eq lc($allowed_value) ) {
                $is_allowed = 1;
                last;
            }
        }
        if ( !$is_allowed ) {
            die("Option '$args{option_name}' must match one of the allowed values: "
                    . join( ",", @{ $args{allowed_values} } )
                    . "." );
        }
    }
}

sub check_symbol {
    my $symbol = shift;
    my $flag   = 1;
    $flag = 0 if ( $symbol eq "-" );
    $flag = 0 if ( $symbol eq "." );
    $flag = 0 if ( $symbol eq "X" );
    $flag = 0 if ( $symbol eq "x" );
    return $flag;

}

sub extract_description_from_fasta_data {
    my ($filename) = @_;
    my @all_description = ();
    open( my $GET_FILE_DATA, $filename )
        or die "Cannot open file '$filename':$!\n\n";
    while (<$GET_FILE_DATA>) {
        chomp $_;
        if ( $_ =~ m/^>/ ) {
            my $description = $_;
            $description =~ s/^>//;
            $description =~ s/\s+$//;
            push( @all_description, $description );
        }
    }
    return @all_description;
}

sub extract_sequence_from_fasta_data {
    my ($filename) = @_;

    # Initialize variables
    open( my $GET_FILE_DATA, $filename )
        or die "Cannot open file '$filename':$!\n\n";

    # Declare and initialize variables
    my $description = <$GET_FILE_DATA>;
    my @file_data   = ();
in: while (1) {
        my $sequence = '';
    out: while ( my $line = <$GET_FILE_DATA> ) {
            if ( eof($GET_FILE_DATA) ) {
                push @file_data, $sequence .= $line;
                last in;
            }
            chomp $line;

            # discard blank line
            if ( $line =~ /^\s*$/ ) {
                next;

                # discard comment line
            }
            elsif ( $line =~ /^\s*#/ ) {
                next;
            }
            $sequence .= $line if not $line =~ /^>/;

            # discard fasta header line

            if ( $line =~ /^>/ ) {
                push @file_data, $sequence;
                last out;
            }

        }
    }
    close($GET_FILE_DATA);
    return @file_data;
}

sub find_index {
    my $threshold      = pop;
    my @identity_array = (@_);
    my %index_hash     = ();
    $index_hash{0} = 0;
    for my $i ( 1 .. ( @identity_array - 1 ) ) {
        my $count = 0;
        for my $j ( 0 .. ( @identity_array - 1 ) ) {
            $count++ if $identity_array[$i][$j] > $threshold / 100;
        }
        $index_hash{$i} = $count;
    }
    my @value
        = ( sort { $index_hash{$b} <=> $index_hash{$a} } keys %index_hash );
    return $value[0] if $index_hash{ $value[0] } != 0;
    return -1 if $index_hash{ $value[0] } == 0;
}

sub set_difference {
    my ( $list_ref, $stop_ref ) = @_;
    my (%temp);

    @temp{@$list_ref} = ();
    foreach (@$stop_ref) {
        delete $temp{$_};
    }
    sort { $a <=> $b } keys %temp;
}

sub write_to_temp_seq {
    my $seq  = shift;
    my $file = shift;
    open( my $TEMP_FILE, '>', $file )
        or die("Cannot open file '$file': $!");
    print $TEMP_FILE $seq;
}

sub message {
    my $verbose = shift;
    my $message = shift;
    if ($verbose) {
        print $message;
    }
}

sub omes {
    my $i               = shift;
    my $j               = shift;
    my $final_file_data = shift;
    my $hash            = {};
    my $hash1           = {};
    foreach ( @{$final_file_data} ) {
        chomp;
        if (    check_symbol( substr( $_, $i - 1, 1 ) )
            and check_symbol( substr( $_, $j - 1, 1 ) ) )
        {
            $hash->{ substr( $_, $i - 1, 1 ) . substr( $_, $j - 1, 1 ) }++;
            $hash->{'validnumber'}++;
            $hash1->{$i}->{ substr $_, $i - 1, 1 }++;
            $hash1->{$j}->{ substr $_, $j - 1, 1 }++;
        }
    }
    my $count       = 0;
    my $total_score = 0;
    for my $key_i ( keys %{ $hash1->{$i} } ) {
        for my $key_j ( keys %{ $hash1->{$j} } ) {
            $count++;
            my $N_valid = $hash->{'validnumber'};
            my $Nex = $hash1->{$i}{$key_i} * $hash1->{$j}{$key_j} / $N_valid;
            my $Nobs;
            $Nobs = $hash->{"$key_i$key_j"}
                if $hash->{"$key_i$key_j"};
            $Nobs = 0 if not $hash->{"$key_i$key_j"};
            my $ometempscore = ( $Nobs - $Nex ) * ( $Nobs - $Nex ) / $N_valid;
            $total_score = $total_score + $ometempscore;
        }
    }
    $total_score = sprintf( "%4.3f", $total_score );
    return $total_score;
}

sub print_usage {
    print <<BLOCK;
USAGE:
perl OMES_score.pl [-arguments]

 -i [FILE]     : file containing protein sequences in fasta format (Required).
 -o [FILE]     : the output file to create (Required).
 -a [FILE]     : write the multiple alignment obtained from T_Coffee to this 
                 file (Optional).
 -n [INTEGER]  : write this number of top OMES scores to the output
                 file. The number of scores written out can exceed
                 this value if there are multiple column pairs sharing
                 the bottom OMES score (Optional; default is to report
                 all results).
 -p [INTEGER]  : remove sequences from the calculation such that no two 
                 sequences exhibit a pairwise percent identity greater than 
                 this value (Optional; default is 90).
 -m [INTEGER]  : do not provide scores if the number of sequences suitable for 
                 calculating the OMES value is less than this number 
                 (Optional; default is 30).
 -stretcher [STRING] : command used to run EMBOSS stretcher (Optional; default
                 is 'stretcher').
 -t_coffee [STRING]  : command used to run T_Coffee (Optional; default is 
                 't_coffee').
 -muscle   [STRING]  : command used to run Muscle (Optional; default
                 is to use 'muscle').
 -use_t_coffee : use T-Coffee instead of Muscle to build the multiple
                 alignment (Optional; default is to use Muscle). 
 -residue [INTEGER]  : only calculate scores for this position in the first
                 sequence (Optional; default is to calculate value for all
                 positions). 
 -v            : provide progress messages (Optional).
 -help         : provide list of arguments and exit (Optional).
 
perl OMES_score.pl -i seqs.fasta -o scores.csv -a alignment.fasta -n 100 -m 3 -p 80 -t
BLOCK
}

