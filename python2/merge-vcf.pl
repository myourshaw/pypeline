#!/usr/bin/perl -w

use strict;
use warnings;
use Carp;
use Vcf;

my $opts = parse_params();
merge_vcf_files($opts);

exit;

#--------------------------------

sub error
{
    my (@msg) = @_;
    if ( scalar @msg ) { croak join('',@msg); }
    die
        "About: Merge the bgzipped and tabix indexed VCF files. (E.g. bgzip file.vcf; tabix -p vcf file.vcf.gz)\n",
        "Usage: merge-vcf [OPTIONS] file1.vcf file2.vcf.gz ... > out.vcf\n",
        "Options:\n",
        "   -c, --chromosomes <list|file>           Do only the given chromosomes (comma-separated list or one chromosome per line in a file).\n",
        "   -d, --remove-duplicates                 If there should be two consecutive rows with the same chr:pos, print only the first one.\n",
        "   -H, --vcf-header <file>                 Use the VCF header\n",
        "   -h, -?, --help                          This help message.\n",
        "\n";
}


sub parse_params
{
    my $opts = {};
    while (my $arg=shift(@ARGV))
    {
        if ( $arg eq '-d' || $arg eq '--remove-duplicates' ) { $$opts{rm_dups}=1; next; }
        if ( $arg eq '-H' || $arg eq '--vcf-header' ) { $$opts{vcf_header}=shift(@ARGV); next; }
        if ( $arg eq '-c' || $arg eq '--chromosomes' ) { $$opts{chromosomes}=shift(@ARGV); next; }
        if ( $arg eq '-?' || $arg eq '-h' || $arg eq '--help' ) { error(); }
        if ( -e $arg ) { push @{$$opts{files}},$arg; next; }
        error("Unknown parameter or non-existent file \"$arg\". Run -? for help.\n");
    }
    if ( !exists($$opts{files}) ) { error() }
    return $opts;
}


# Returns the common prefix of the files.
sub common_prefix
{
    my ($files) = @_;
    my @paths;
    my $len = -1;
    for my $file (@$files)
    {
        my @path = split(m{/+},$file);
        if ( $len<0 || $len>scalar @path ) { $len=scalar @path; }
        push @paths, \@path;
    }
    my @common;
    for (my $i=0; $i<$len; $i++)
    {
        my $identical=1;
        for (my $ifile=1; $ifile<scalar @paths; $ifile++)
        {
            if ( $paths[$ifile]->[$i] ne $paths[0]->[$i] ) { $identical=0; last; }
        }
        if ( !$identical ) { last; }
        push @common, $paths[0]->[$i];
    }
    return join('/+',@common);
}


sub read_chrom_list
{
    my ($opts) = @_;

    my @chroms = ();
    if ( exists($$opts{chromosomes}) ) 
    { 
        if ( -e $$opts{chromosomes} )
        {
            open(my $chrms,'<',$$opts{chromosomes}) or error("$$opts{chromosomes}: $!");
            while (my $line=<$chrms>)
            {
                chomp($line);
                push @chroms, $line;
            }
            close($chrms);
        }
        else
        {
            @chroms = split(/,/,$$opts{chromosomes}); 
        }
    }
    return (@chroms);
}


sub init_cols
{
    my ($opts,$vcf_out) = @_;

    my $prefix;
    my @chroms = read_chrom_list($opts);
    my @vcfs;
    my @cols;
    my %has_chrom;
    my %col_names;
    my $icol = 9;

    if ( !$$opts{has_col_names} ) { $prefix = common_prefix($$opts{files}); }

    # Go through all files and read header, obtain list of chromosomes. The file names will be used for columns, unless
    #   they were read from the header.
    for my $file (@{$$opts{files}})
    {
        my $vcf = Vcf->new(file=>$file);
        $vcf->parse_header();
        $vcf->close();
        push @vcfs, $vcf;

        # Update the list of known chromosomes
        if ( !exists($$opts{chromosomes}) )
        {
            my $chrms = $vcf->get_chromosomes();
            for my $chr (@$chrms)
            {
                if ( exists($has_chrom{$chr}) ) { next; }
                $has_chrom{$chr} = 1;
                push @chroms, $chr;
            }
        }

        my $col_prefix = '';
        if ( !$$opts{has_col_names} )
        {
            # Make the column names nice - strip common prefix and the suffix .vcf.gz
            $col_prefix = $file; 
            $col_prefix =~ s{^/*$prefix/*}{};
            $col_prefix =~ s/\.gz$//i;
            $col_prefix =~ s/\.vcf$//i;
            $col_prefix .= '_';
        }

        # Create good names for the columns in the merged vcf file
        my @vcf_cols = @{$$vcf{columns}};
        $$vcf{__col_names} = [];
        for my $col (@vcf_cols[9..$#vcf_cols])
        {
            my $col_name = $col;
            if ( $$opts{has_col_names} ) 
            {
                if ( $icol >= @{$$vcf_out{columns}} ) { error("Fewer columns in the header than in the VCF files total.\n"); }
                $col_name = $$vcf_out{columns}[$icol];
                $icol++;

                if ( exists($col_names{$col_name}) ) { error("The column names not unique in the header: $col_name\n"); }
            }
            else
            {
                if ( exists($col_names{$col_name}) ) { $col_name = $col_prefix.$col; }
                if ( exists($col_names{$col_name}) ) { warn("FIXME: the column name [$col_name] not unique.\n"); }
            }
            warn("Using column name '$col_name' for $file:$col\n");
            $col_names{$col_name} = 1;

            push @cols, $col_name;
            push @{$$vcf{__col_names}}, $col_name;
        }
    }

    if ( $$opts{has_col_names} && $icol!=@{$$vcf_out{columns}} ) { error("More columns in the header than in the VCF files total.\n"); }

    $$opts{vcfs} = \@vcfs;
    $$opts{cols} = \@cols;
    $$opts{chroms} = \@chroms;
}


sub merge_vcf_files
{
    my ($opts) = @_;

    # Create output VCF
    my $vcf_out;
    if ( $$opts{vcf_header} )
    {
        $vcf_out = Vcf->new(file=>$$opts{vcf_header});
        $vcf_out->parse_header();
        if ( $$vcf_out{columns} && @{$$vcf_out{columns}} ) { $$opts{has_col_names}=1; }
    }
    else
    {
        $vcf_out = Vcf->new();
    }

    my $ignore_dups = $$opts{rm_dups} ? 1 : 0;

    init_cols($opts,$vcf_out);
    my @chroms = @{$$opts{chroms}};
    my @cols   = @{$$opts{cols}};
    my @vcfs   = @{$$opts{vcfs}};

    # Get the header of the output VCF ready
    $vcf_out->add_columns(@cols);
    if ( !$$vcf_out{has_header} )
    {
        # To get the missig fields filled by the default values
        for my $hline (@{$vcfs[0]->{header_lines}})
        {
            $vcf_out->add_header_line($hline);
        }
    }
    $vcf_out->recalc_ac_an(2);
    $vcf_out->add_header_line({key=>'INFO',ID=>'AC',Number=>-1,Type=>'Integer',Description=>'Allele count in genotypes'});
    $vcf_out->add_header_line({key=>'INFO',ID=>'AN',Number=>1,Type=>'Integer',Description=>'Total number of alleles in called genotypes'});
    print $vcf_out->format_header();

    # Go through all VCF files simultaneously and output each line, one chromosome at a time.
    for my $chr (@chroms)
    {
        # Open files
        for my $vcf (@vcfs) 
        { 
            delete($$vcf{last_line});
            $vcf->open_tabix($chr); 
            advance_position($vcf);
        }

        while (1)
        {
            my $pos = get_min_position(\@vcfs);
            if ( !defined $pos ) { last; }

            my %out;
            $out{CHROM}  = $chr;
            $out{POS}    = $pos;
            $out{ID}     = '.';
            $out{ALT}    = [];
            $out{QUAL}   = $$vcf_out{defaults}{QUAL};    # It is not clear how the QUAL should be merged, output empty for now
            $out{FILTER} = ['.'];
            $out{FORMAT} = [];
            my %format;

            my %ref_alt_map = ();
            # Find out the REFs and ALTs: in VCFv4.0, the REFs can differ and ALTs must be converted
            for my $vcf (@vcfs)
            {
                my $line = $$vcf{last_line};
                if ( !$line or $pos ne $$line{POS} ) { next; }
                my $ref = $$line{REF};
                for my $alt (@{$$line{ALT}}) { $ref_alt_map{$ref}{$alt}=1; }
            }
            # Do the REF,ALT conversion only when necessary
            my $new_ref; 
            if ( scalar keys %ref_alt_map > 1 ) { $new_ref = $vcf_out->fill_ref_alt_mapping(\%ref_alt_map); }
            for my $vcf (@vcfs)
            {
                my $line = $$vcf{last_line};

                # If this file does not have a record for this position, than for all its columns output undef gtype
                if ( !$line or $pos ne $$line{POS} )
                {
                    for (my $i=0; $i<@{$$vcf{__col_names}}; $i++)
                    {
                        my $name = $$vcf{__col_names}->[$i];
                        $out{gtypes}{$name}{GT} = $$vcf_out{defaults}{GT};
                    }
                    next;
                }

                if ( $$line{ID} ne '.' && $out{ID} eq '.' ) { $out{ID}=$$line{ID}; }

                # Remember the FORMAT fields
                for my $field (@{$$line{FORMAT}}) { $format{$field} = 1; }

                my $ref = $$line{REF};

                # Now fill in the genotype information for each column
                for (my $i=0; $i<@{$$vcf{__col_names}}; $i++)
                {
                    my $gt = $$vcf{columns}->[$i+9];

                    # This is to convert 0/1 to G/C
                    my ($a,$sep,$b) = $vcf->parse_alleles($line,$gt);
                    if ( defined $new_ref )
                    {
                        $a = $ref_alt_map{$ref}{$a};
                        $b = $ref_alt_map{$ref}{$b};
                    }
                    if ( !defined $a ) { $a=''; }
                    if ( !defined $b ) { $b=''; }

                    my $name = $$vcf{__col_names}->[$i];
                    $out{gtypes}{$name} = $$line{gtypes}{$$vcf{columns}->[$i+9]};
                    $out{gtypes}{$name}{GT}  = $a.$sep.$b;
                }
                $out{REF} = defined $new_ref ? $new_ref : $ref;
                advance_position($vcf,$ignore_dups);
            }

            # The GT field must come as first
            delete($format{GT});
            $out{FORMAT} = ['GT'];
            for my $key (keys %format) { push @{$out{FORMAT}},$key; }

            $vcf_out->format_genotype_strings(\%out);
            print $vcf_out->format_line(\%out);
        }
    }
}


sub advance_position
{
    my ($vcf,$ignore_dups) = @_;

    if ( exists($$vcf{last_line}) && !$$vcf{last_line} ) { return; }

    my $line;
    while (!$line)
    {
        $line = $vcf->next_data_hash();
        if ( !$line )
        {
            $$vcf{last_line} = $line;
            return;
        }
        if ( !$$vcf{last_line} ) { last; }

        if ( $$vcf{last_line}{POS} eq $$line{POS} ) 
        {
            print STDERR "The position appeared twice: $$vcf{file} .. $$line{CHROM}:$$line{POS}\n";

            # This is the only reason for the while loop: if ignoring dups, get the next line
            if ( $ignore_dups ) { undef($line); }
        }
        elsif ( $$vcf{last_line}{POS} > $$line{POS})
        { 
            error("Wrong order: $$vcf{file} .. $$line{CHROM}:$$line{POS} comes after $$vcf{last_line}{CHROM}:$$vcf{last_line}{POS}\n");
        }
    }

    $$vcf{last_line} = $line;

    return;
}


sub get_min_position
{
    my ($vcfs) = @_;
    my ($pos,$ref);
    for my $vcf (@$vcfs)
    {
        my $line = $$vcf{last_line};
        if ( !$line ) { next; }

        # Designate this position as the minimum of all the files if:
        # .. is this the first file?
        if ( !defined $pos )
        {
            $pos = $$line{POS};
            $ref = $$line{REF};

            next;
        }
        
        # .. has this file lower position?
        if ( $pos>$$line{POS} )
        {
            $pos = $$line{POS};
            $ref = $$line{REF};

            next;
        }
    }
    return $pos;
}
