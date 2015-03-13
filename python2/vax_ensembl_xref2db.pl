#!/usr/bin/perl

#-host cortex.local -port 3306 -user ensembl -pass ensembl -output -host cortex -port 3306 -user ensembl -pass ensembl -output /scratch0/tmp/myourshaw/ensembl/xref_ensembl-68.txt
#-output /scratch1/tmp/myourshaw/ensembl/xref/xref_ensembl_73.txt
#-species mus_musculus -output /share/apps/myourshaw/vax/75/ensembl_xref/ensembl_xref_75_mus_musculus_with_dups.txt
#-species mus_musculus -output /scratch1/vax/ontology/go_gene/ensembl_xref_75_mus_musculus_tmp.txt


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

my @OUTPUT_COLS = qw(
    ENSG
    ENST
    ENSP
    dbname
    primary_id
    display_id
    chrom
    geneStart
    geneEnd
    strand
    txStart
    txEnd
    cdsStart
    cdsEnd
    biotype
    is_canonical
);
my $null = '';

# configure from command line opts
my $config = &configure(scalar @ARGV);

#&main($config);
#print;

# this is the main sub-routine - it needs the configured $config hash
#sub main {
    #my $config = shift;
    
    open OUT,">",$config->{output} or die "can't open $config->{output}";
    print OUT '#' . (join "\t", @OUTPUT_COLS) . "\n";

    my @genes = @{$config->{ga}->fetch_all()};
    my @stable_ids = map {$_->{stable_id}} @genes;
    #my @stable_ids = qw(ENSMUSG00000038007 ENSMUSG00000091609);
    undef @genes;
    #@genes = @genes[-50..-1];
    my $gene_count = scalar(@stable_ids);
    my $this_gene_number = 0;
    print "##$gene_count Ensembl genes\n";
    print "#gene_number\tENSG\n";
    foreach my $ensg(@stable_ids){
        $this_gene_number++;
        print "$this_gene_number\t$ensg\t";
        my ($chrom, $geneStart, $geneEnd, $strand, $biotype);
        my $gene = $config->{ga}->fetch_by_stable_id($ensg);
        #print "\n<".$gene->{is_current}.">\t<".$gene->{is_known}.">\t<".$gene->{is_reference}.">\t".$gene->{is_stored}."\n";
        #$ensg = $gene->{stable_id};
        $chrom = $gene->{slice}->{seq_region_name};
        $geneStart = $gene->{start};
        $geneEnd = $gene->{end};
        $strand = $gene->{strand};
        $biotype = $gene->{biotype};
        my $xrefs = $gene->get_all_DBLinks();
        foreach my $xref(@{$xrefs}){
            #print ".";
            my %line;
            $line{ENSG} = $ensg;
            $line{chrom} = $chrom;
            $line{geneStart} = $geneStart;
            $line{geneEnd} = $geneEnd;
            $line{strand} = $strand;
            $line{dbname} = $xref->{dbname};
            $line{primary_id} = $xref->primary_id();
            $line{display_id} = $xref->display_id();
            $line{biotype} = $biotype;
            $line{is_canonical} = '';
            my $output = join "\t", map { $line{$_} || $null } @OUTPUT_COLS;
            print OUT "$output\n";
        }
        print "\n";
        my ($enst, $ensp, $txStart, $txEnd, $cdsStart, $cdsEnd, $tx_biotype, $tx_is_canonical);
		my @transcripts = @{$gene->get_all_Transcripts};
		if (@transcripts) {
            foreach my $transcript(@transcripts){
                #print "\t$transcript->{stable_id}\t";
                $enst = $transcript->{stable_id};
                $txStart = $transcript->{start};
                $txEnd = $transcript->{end};
                $cdsStart = $transcript->coding_region_start;
                $cdsEnd = $transcript->coding_region_end;
                $tx_biotype = $transcript->{biotype};
                $tx_is_canonical = $transcript->{is_canonical};
                my $translation = $transcript->translation;
                if(defined($translation)){
                    $ensp = $translation->{stable_id};
                }
                my $xrefs = $transcript->get_all_DBLinks();
                foreach my $xref(@{$xrefs}){
                    #print ".";
                    my %line;
                    $line{ENSG} = $ensg;
                    $line{ENST} = $enst;
                    $line{ENSP} = $ensp;
                    $line{dbname} = $xref->{dbname};
                    $line{primary_id} = $xref->primary_id();
                    $line{display_id} = $xref->display_id();
                    $line{chrom} = $chrom;
                    $line{geneStart} = $geneStart;
                    $line{geneEnd} = $geneEnd;
                    $line{strand} = $strand;
                    $line{txStart} = $txStart;
                    $line{txEnd} = $txEnd;
                    $line{cdsStart} = $cdsStart;
                    $line{cdsEnd} = $cdsEnd;
                    $line{biotype} = $tx_biotype;
                    $line{is_canonical} = defined($tx_is_canonical) ? $tx_is_canonical ? 1 : 0 : '';
                    my $output = join "\t", map { $line{$_} || $null } @OUTPUT_COLS;
                    print OUT "$output\n";
                }
                #print "\n";
            }
        }
    }
    print "##closing output file\n";
    close OUT or die "can't close $config->{output}";
    open(FH,">",$config->{done_file}) or die "Can't create $config->{done_file}: $!";
    close(FH);
    print "##ensembl_db2xref done\n";
#    return;
#}

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
    $config->{species} ||= 'homo_sapiens';
    
    # connect to databases
    $config->{reg} = &connect_to_dbs($config);

    &get_adaptors($config);
   
    #info for file headers
    my $start_time = &get_time;

	my $core_mca = $config->{reg}->get_adaptor($config->{species}, 'core', 'metacontainer');
    my $ensembl_api_version = $config->{reg}->software_version;
    my $ensembl_db_version = (defined $core_mca && $core_mca->get_schema_version ? $core_mca->get_schema_version : '?');
	$config->{config_options_header} .= "## Connected to ".$core_mca->dbc->dbname." on ".$core_mca->dbc->host."\n" if defined $core_mca;
    
    $config->{output} ||= "ensembl_xref".$ensembl_db_version."_with_dups.txt";
    my ($v,$d,$f) = File::Spec->splitpath($config->{output});
    $config->{done_file} = File::Spec->join($v, $d, '.'.$f.'.done');
    if (-e $config->{done_file}) {
        print sprintf('Ensembl xref2db already done. To reschedule "rm %s"', $config->{done_file});
        exit;
    }
    
    
	my $prolog =<<"PROLOG_END";
##Creating Ensembl xref text file $config->{output}.
##Output produced at $start_time.
##Using API version $ensembl_api_version, DB version $ensembl_db_version.
##This file will contain some duplicate entries.
##to remove them in a MySQL database, run:
##CREATE TABLE `db`.`ensembl_xref` AS
##SELECT DISTINCT ENSG,ENST,ENSP,dbname,primary_id,display_id,chrom,geneStart,geneEnd,strand,txStart,txEnd,cdsStart,cdsEnd,biotype,is_canonical
##FROM `db`.`ensembl_xref_with_dups`;
##UPDATE `db`.`ensembl_xref` SET `ENST` = '' WHERE `ENST` IS NULL;
##UPDATE `db`.`ensembl_xref` SET `ENSP` = '' WHERE `ENSP` IS NULL;
##ALTER TABLE `db`.`ensembl_xref` ADD PRIMARY KEY (ENSG, ENST, ENSP, dbname, primary_id, display_id);
##DROP TABLE `db`.`ensembl_xref_with_dups`;
PROLOG_END

print $prolog;

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
