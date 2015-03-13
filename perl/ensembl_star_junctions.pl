#!/usr/bin/perl

#-output /scratch1/tmp/myourshaw/ensembl/star/ensembl_74_star_junctions.txt

#create intron file for STAR RNASeq aligner
#Chr \tab\ Start \tab\ End \tab\ Strand(+or-)
#where Start and End are first and last bases of the introns (1-based chromosome coordinates).
#usage to generate genome index used in pass 1
#STAR --runMode genomeGenerate --genomeDir $genomeDir --genomeFastaFiles b37.fa --runThreadN <n> \
#--sjdbFileChrStartEnd /path/to/ensembl_75_star_junctions.txt --sjdbOverhang 75;

use Getopt::Long;
use FileHandle;
use File::Basename qw( fileparse );
use File::Path qw( make_path );
use File::Spec;
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

my @OUTPUT_COLS = qw(
    Chr
    Start
    End
    Strand
);
my $null = '';

# configure from command line opts
my $config = &configure(scalar @ARGV);

&main($config);
#print;

# this is the main sub-routine - it needs the configured $config hash
sub main {
    my $config = shift;
    
    my ( $outfile, $directories ) = fileparse $config->{output};
    if ( !$outfile ) {
        $outfile = 'ensembl_star_junctions.txt';
        $config->{output} = File::Spec->catfile( $config->{output}, $outfile );
    }
    
    if ( !-d $directories ) {
        make_path $directories or die "Failed to create path: $directories";
    }

    print "create STAR junction file\n";
    print "getting all introns ...\n";
    my @all_introns;
    open OUT,">",$config->{output} or die "can't open $config->{output}";
    #TODO: does the file require/allow a header?
    #print OUT '#' . (join "\t", @OUTPUT_COLS) . "\n";

    my @genes = @{$config->{ga}->fetch_all()};
    my @stable_ids = map {$_->{stable_id}} @genes;
    undef @genes;
    #@genes = @genes[-50..-1];
    my $gene_count = scalar(@stable_ids);
    my $this_gene_number = 0;
    print "$gene_count Ensembl genes\n";
    foreach my $ensg(@stable_ids){
        $this_gene_number++;
        #print "$this_gene_number\t$ensg\n";
        my $gene = $config->{ga}->fetch_by_stable_id($ensg);
		my @transcripts = @{$gene->get_all_Transcripts};
		if (@transcripts) {
            foreach my $transcript(@transcripts){
                my @introns = @{ $transcript->get_all_Introns() };
                if (@introns) {
                    foreach my $intron(@introns){
                        my ($chrom, $intronStart, $intronEnd, $strand);
                        $chrom = $intron->{slice}->{seq_region_name};
                        $start = $intron->start();
                        $end = $intron->end();
                        $strand = $intron->strand();
                        push @all_introns, {zero_chr=>chr2zerochr($chrom), start=>$start, end=>$end, strand=>$strand};
                    }
                }
            }
        }
    }
	if($#all_introns){
		print "sorting $#all_introns introns and eliminating duplicates ...\n";
		my @sorted_introns = sort {$a->{zero_chr} cmp $b->{zero_chr} || $a->{start} <=> $b->{start} || $a->{end} <=> $b->{end} || $b->{strand} <=> $a->{strand}} @all_introns;
        my $unique_intron_count = 0;
		my $that_chr;
		my $this_chr;
		my $that_start;
		my $this_start;
		my $that_end;
		my $this_end;
		my $that_strand;
		my $this_strand;
		foreach my $intron(@sorted_introns){
			$this_chr = $intron->{zero_chr};
			$this_start = $intron->{start};
			$this_end = $intron->{end};
			$this_strand = $intron->{strand};
			if (!defined($that_chr) || $this_chr ne $that_chr || $this_start != $that_start || $this_end != $that_end || $this_strand != $that_strand){
                my $chrom = zerochr2chr($this_chr);
                my $star_strand = $this_strand == 1 ? '+' : '-';
                print OUT "$chrom\t$this_start\t$this_end\t$star_strand\n";
				$unique_intron_count++;
				$that_chr = $this_chr;
				$that_start = $this_start;
				$that_end = $this_end;
                $that_strand = $this_strand;
			}
		}
        print "finished, output $unique_intron_count unique intron regions\n";
    }
} #main

#convert chromosome to number
sub chr2zerochr {
    my $chr = uc(shift);
    return $chr =~ /^\d$/ ? '0'.$chr : $chr;
}
#convert number to chromosome
sub zerochr2chr {
    my $chr = shift;
    return $chr  =~ /^0\d$/ ? substr($chr,1) : $chr;
}


# sets up configuration hash that is used throughout the script
sub configure {
    my $args = shift;
    
    my $config = {};
    
    GetOptions(
        $config,
        'host=s',
        'port=i',
        'user=s',
        'pass:s',
        'output=s',
        'species=s',
    );
    $config->{host} ||= "cortex.local";
    $config->{port} ||= 3306;
    $config->{user} ||= "ensembl";
    $config->{pass} ||= "ensembl";
    #$config->{output} ||= "ensembl_star_junctions.txt";
    $config->{species} ||= 'homo_sapiens';
    
    # connect to databases
    $config->{reg} = &connect_to_dbs($config);

    &get_adaptors($config);
   
    return $config;
}

# connects to DBs
sub connect_to_dbs {
    my $config = shift;
    # get registry
    my $reg = 'Bio::EnsEMBL::Registry';
    $reg->load_registry_from_db(
        -host => $config->{host},
        -user => $config->{user},
        -pass => $config->{pass},
        -port => $config->{port}
        );
    eval { $reg->set_reconnect_when_lost() };
    
    return $reg;        
}

# get adaptors from DB
sub get_adaptors {
    my $config = shift;
    
    die "ERROR: No registry" unless defined $config->{reg};
    
    $config->{sa}  = $config->{reg}->get_adaptor($config->{species}, 'core', 'slice');
    $config->{ga}  = $config->{reg}->get_adaptor($config->{species}, 'core', 'gene');
    $config->{ta}  = $config->{reg}->get_adaptor($config->{species}, 'core', 'transcript');
    $config->{vfa} = $config->{reg}->get_adaptor($config->{species}, 'variation', 'variationfeature');
    $config->{tva} = $config->{reg}->get_adaptor($config->{species}, 'variation', 'transcriptvariation');
    # check we got slice adaptor - can't continue without a core DB
    die("ERROR: Could not connect to core database\n") unless defined $config->{ga};
}
