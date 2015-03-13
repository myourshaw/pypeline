#!/usr/bin/perl -w
use strict;

use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
);

sub feature2string
{
    my $feature = shift;

    my $stable_id  = $feature->stable_id();
    my $seq_region = $feature->slice->seq_region_name();
    my $start      = $feature->start();
    my $end        = $feature->end();
    my $strand     = $feature->strand();

    return sprintf( "%s: %s:%d-%d (%+d)",
        $stable_id, $seq_region, $start, $end, $strand );
}

my $slice_adaptor = $registry->get_adaptor( 'Human', 'Core', 'Slice' );
my $slice = $slice_adaptor->fetch_by_region( 'chromosome', 'X', 1e6, 10e6 );

my @genes = @{ $slice->get_all_Genes() };
while ( my $gene = shift @genes ) {
    my $gstring = feature2string($gene);
    print "$gstring\n";

    my @transcripts = @{ $gene->get_all_Transcripts() };
    while ( my $transcript = shift @transcripts ) {
        my $tstring = feature2string($transcript);
        print "\t$tstring\n";

        my @exons = @{ $transcript->get_all_Exons() };
        while ( my $exon = shift @exons ) {
            my $estring = feature2string($exon);
            print "\t\t$estring\n";
        }
        
        #Introns are not defined explicitly in the database but can be obtained by the Transcript method get_all_Introns()
         my @introns = @{ $transcript->get_all_Introns() };
         while ( my $intron = shift @introns ) {
            my $istring = feature2string($intron);
            print "\t\t$istring\n";
        }
   }
}


@genes = @{ $slice->get_all_Genes() };
while ( my $gene = shift @genes ) {
    my $gstring = feature2string($gene);
    print "$gstring\n";

    my @transcripts = @{ $gene->get_all_Transcripts() };
    while ( my $transcript = shift @transcripts ) {
        # The spliced_seq() method returns the concatenation of the exon
        # sequences.  This is the cDNA of the transcript
        print "cDNA: ", $transcript->spliced_seq(), "\n";
        
        # The translateable_seq() method returns only the CDS of the transcript
        print "CDS: ", $transcript->translateable_seq(), "\n";
        
        # UTR sequences are obtained via the five_prime_utr() and
        # three_prime_utr() methods
        my $fiv_utr = $transcript->five_prime_utr();
        my $thr_utr = $transcript->three_prime_utr();
        
        print "5' UTR: ", ( defined $fiv_utr ? $fiv_utr->seq() : "None" ), "\n";
        print "3' UTR: ", ( defined $thr_utr ? $thr_utr->seq() : "None" ), "\n";
        
        # The peptide sequence is obtained from the translate() method.  If the
        # transcript is non-coding, undef is returned.
        my $peptide = $transcript->translate();
        
        print "Translation: ", ( defined $peptide ? $peptide->seq() : "None" ), "\n";
    }
}

