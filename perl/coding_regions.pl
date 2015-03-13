#!/usr/bin/perl -w
use strict;
use warnings;
use Bio::Seq;

my $jump_to_input_file = '';
my $output_file = '';
my $format = 'fasta';
my $include_gene_5_prime_utr = 0;
my $include_gene_3_prime_utr = 0;
my $include_exon_5_prime_utr = 0;
my $include_exon_3_prime_utr = 0;
my $include_exon_cds = 1;
my $merge_exon_cds = 1;
my $include_introns = 0;
my $roi_extra_upstream = 0;
my $roi_extra_downstream = 0;
my $roi_extra_upstream_exon = 0;
my $roi_extra_downstream_exon = 0;
my $roi_extra_upstream_exon_1 = 0;
my $roi_extra_downstream_exon_1 = 0;
my $roi_extra_upstream_intron = 0;
my $include_transcribed_strand = 1;
my $include_reverse_complement_strand = 1;
my $padding_upstream = 0;
my $padding_downstream= 0;
my $show_variations = 1;
my $min_length = 0;
my $max_length = 0;
my $warnings_only = 0;

#read from a list of filenames on the command line
#read from SDTIN if no filenames were given
#"-" indicates STDIN
#"someprogram |" indicates the output of another program
while(<>){
    if (/\$/){
    print @_;
    }
    else{
        print;
    }
}
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
my $utr_start;
my $utr_end;
my $strand;
my %coordinate;
my @starts;
my @ends;
my @ranges;

@genes = @{$gene_adaptor->fetch_all_by_external_name($gene_name)};
$gene = $genes[0];
$strand = $gene->strand;
printf "gene: %s strand: %s\n",$gene->external_name, $strand==1?"+":"-";
@transcripts = @{$gene->get_all_Transcripts};
for (my $t=0;$t<@transcripts;$t++){
    my $transcript = $transcripts[$t];
    print "transcript ", $t+1, ": ",$transcript->start,"-",$transcript->end,"\n";
    if(defined($include_gene_5_prime_utr)){
        $five_utr = $transcript->five_prime_utr();
        if(defined($five_utr)){
            $utr_start = $strand==1?$transcript->start:$transcript->end-$five_utr->length+1;
            $utr_end = $strand==1?$transcript->start+$five_utr->length-1:$transcript->end;
            print "\t5' UTR: ", $utr_start, "-", $utr_end, "\n";
            push @ranges, new Bio::Range(-start=>$utr_start, -end=>$utr_end, -strand=>$strand);
        }
    }
    if(defined($include_exon_5_prime_utr) || defined($include_exon_3_prime_utr) || defined($include_exon_cds)){
        @exons = @{$transcript->get_all_Exons};
        for (my $e=$strand==1?0:(@exons-1);$strand==1?$e<@exons:$e>=0;$strand==1?$e++:$e--){
            my $exon = $exons[$e];
            my $exon_number = $strand==1?$e+1:@exons-$e;
            print "\texon ", $exon_number, ": ", $exon->start,"-",$exon->end, "\n";
            if(defined($exons[$e]->coding_region_start($transcript))){
                if(defined($include_exon_5_prime_utr)){
                    if($strand==1){
                        if($exon->coding_region_start($transcript) > $exon->start){
                            $utr_start = $exon->start;
                            $utr_end = $exon->coding_region_start($transcript)-1;
                        }
                    }
                    else{
                        if($exon->coding_region_end($transcript) < $exon->end){
                            $utr_start = $exon->coding_region_end($transcript)+1;
                            $utr_end = $exon->end;
                        }
                    }
                    print "\t\t", "5' UTR: ", $utr_start, "-", $utr_end, "\n";
                    push @ranges, new Bio::Range(-start=>$utr_start, -end=>$utr_end, -strand=>$strand);
                }
                if(defined($include_exon_cds)){
                    print "\t\tcoding: ", $exon->coding_region_start($transcript), " ", $exon->coding_region_end($transcript), "\n";
                    push @ranges, new Bio::Range(-start=>$exon->coding_region_start($transcript), -end=>$exon->coding_region_end($transcript), -strand=>$strand);
                }
                if(defined($include_exon_3_prime_utr)){
                    if($strand==1){
                        if($exon->coding_region_end($transcript) < $exon->end){
                            $utr_start = $exon->coding_region_end($transcript)+1;
                            $utr_end = $exon->end;
                        }
                    }
                    else{
                        if($exon->coding_region_start($transcript) > $exon->start){
                            $utr_start = $exon->start;
                            $utr_end = $exon->coding_region_start($transcript)-1;
                        }
                    }
                    print "\t\t", "3' UTR: ", $utr_start, "-", $utr_end, "\n";
                    push @ranges, new Bio::Range(-start=>$utr_start, -end=>$utr_end, -strand=>$strand);
                }
            }
            elsif (defined($include_exon_5_prime_utr) || defined($include_exon_3_prime_utr)){
                print "\texon ", $exon_number, ": ", $exon->start,"-",$exon->end, " non-coding\n";
                push @ranges, new Bio::Range(-start=>$exon->start, -end=>$exon->end, -strand=>$strand);
            }
        }
    }
    if(defined($include_gene_3_prime_utr)){
        $three_utr = $transcript->three_prime_utr();
        if(defined $three_utr){
            $utr_start = $strand==1?$transcript-$three_utr->length+1->end:$transcript->start;
            $utr_end = $strand==1?$transcript->end:$transcript->start+$three_utr->length-1;
            print "\t3' UTR: ", $utr_start, "-", $utr_end, "\n";
            push @ranges, new Bio::Range(-start=>$utr_start, -end=>$utr_end, -strand=>$strand);
        }
    }
}

for(my $i=$#ranges;$i>0;$i--){
    if($ranges[$i]->start >= $ranges[$i-1]->start && $ranges[$i]->start <= $ranges[$i-1]->end){
        if($ranges[$i]->end > $ranges[$i-1]->end){
            $ranges[$i-1]->end($ranges[$i]->end);
        }
        splice(@ranges,$i,1);
    }
}

