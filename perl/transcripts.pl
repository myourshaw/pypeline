#!/usr/bin/perl -w
use strict;
use warnings;
use Bio::Seq;

my $gene_name = "DCX";
use Bio::EnsEMBL::Registry;
my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db(-host => 'ensembldb.ensembl.org',-user => 'anonymous');
my $gene_adaptor = $registry->get_adaptor( 'human', 'core', 'gene');
my $exon_adaptor = $registry->get_adaptor( 'human', 'core', 'exon');
my $slice_adaptor = $registry->get_adaptor( 'human', 'core', 'slice');

my @genes;
my $gene;
my @exons;
#my $exon;
my @transcripts;
#my $transcript;
my $slice;
my $five_utr;
my $three_utr;
my $five_utr_start;
my $five_utr_end;
my $three_utr_start;
my $three_utr_end;
my $strand;

@genes = @{$gene_adaptor->fetch_all_by_external_name($gene_name)};
$gene = $genes[0];
$strand = $gene->strand;
printf "gene: %s strand: %s\n",$gene->external_name, $strand==1?"+":"-";
@transcripts = @{$gene->get_all_Transcripts};
for (my $t=0;$t<@transcripts;$t++)
{
    my $transcript = $transcripts[$t];
    print "transcript ", $t+1, ": ",$transcript->start,"-",$transcript->end,"\n";
    $five_utr = $transcript->five_prime_utr();
    $three_utr = $transcript->three_prime_utr();
    if(defined $five_utr){
        $five_utr_start = $strand==1?$transcript->start:$transcript->end-$five_utr->length+1;
        $five_utr_end = $strand==1?$transcript->start+$five_utr->length-1:$transcript->end;
        if($strand==1){
            print "\t5' UTR: ", $five_utr_start, "-", $five_utr_end, "\n";
        }
    }
    if(defined $three_utr){
        $three_utr_start = $strand==1?$transcript-$three_utr->length+1->end:$transcript->start;
        $three_utr_end = $strand==1?$transcript->end:$transcript->start+$three_utr->length-1;
        if($strand!=1){
            print "\t3' UTR: ", $three_utr_start, "-", $three_utr_end, "\n";
        }
    }
    @exons = @{$transcript->get_all_Exons};
    for (my $e=$strand==1?0:(@exons-1);$strand==1?$e<@exons:$e>=0;$strand==1?$e++:$e--){
        my $exon = $exons[$e];
        my $exon_number = $strand==1?$e+1:@exons-$e;
        if(defined($exons[$e]->coding_region_start($transcript))){
            if($exon->coding_region_start($transcript) == $exon->start && $exon->coding_region_end($transcript) == $exon->end){
                print "\texon ", $exon_number, ": ", $exon->start,"-",$exon->end, " coding\n";
            }
            else{
                print "\texon ", $exon_number, ": ", $exon->start,"-",$exon->end, "\n";
                if($exon->coding_region_start($transcript) > $exon->start){
                    my $utr_start = $exon->start;
                    my $utr_end = $exon->coding_region_start($transcript)-1;
                    print "\t\t", 4+$exon->strand, "' UTR: ", $utr_start, "-", $utr_end, "\n";
                }
                print "\t\tcoding: ", $exon->coding_region_start($transcript), " ", $exon->coding_region_end($transcript), "\n";
                if($exon->coding_region_end($transcript) < $exon->end){
                    my $utr_start = $exon->coding_region_end($transcript)+1;
                    my $utr_end = $exon->end;
                    print "\t\t", 4-$exon->strand, "' UTR: ", $utr_start, "-", $utr_end, "\n";
                }
            }
        }
        else{
            print "\texon ", $exon_number, ": ", $exon->start,"-",$exon->end, " non-coding\n";
        }
    }
    if(defined $five_utr){
        my $five_utr_start = $strand==1?$transcript->start:$transcript->end-$five_utr->length+1;
        my $five_utr_end = $strand==1?$transcript->start+$five_utr->length-1:$transcript->end;
        if($strand!=1){
            print "\t5' UTR: ", $five_utr_start, "-", $five_utr_end, "\n";
        }
    }
    if(defined $three_utr){
        if($strand==1){
            print "\t3' UTR: ", $three_utr_start, "-", $three_utr_end, "\n";
        }
    }
}


