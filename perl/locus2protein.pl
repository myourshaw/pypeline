#!/usr/bin/perl -w
use strict;

use warnings;
use Getopt::Long;
use FileHandle;
use File::Spec;

my %opts;
my $output;
GetOptions ('output=s' => \$output);

unless (@ARGV) {
	print("requires list of genomic coordinates, e.g. 2:9630329 2:9630569-9630570\n");
	exit;
}
#my @coords = qw(2:9630329 2:9630569-9630570);
my @coords = @ARGV;


use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'cortex.local',
    -user => 'ensembl',
    -pass => 'ensembl'
);

my $sa = $registry->get_adaptor( 'homo_sapiens', 'core', 'slice' );
my $ga = $registry->get_adaptor('homo_sapiens', 'core', 'gene');

print "#coord\tchrom_start\tchrom_end\tref\tensg\tenst\tensp\thgnc\tcanonical\tprotein_length\taa_start\taa_end\tprotein_seq\n";

for my $coord(@coords){
    my @c = split /:|-/, $coord;
    my $chrom = $c[0];
    my $start = $c[1];
    my $end = defined($c[2]) ? $c[2] : $c[1];
    my $slice = $sa->fetch_by_region('chromosome', $chrom, $start, $end);
    my $ref = $slice->seq();
    my @genes = @{$slice->get_all_Genes()};
    while (my $gene = shift @genes) {
        my $ensg = $gene->stable_id();
        $gene->external_db('HGNC');
        my $hgnc = $gene->external_name() || '';
        my $transcripts = $gene->get_all_Transcripts();
        my $enst;
        my $ensp;
        while (my $tx = shift @{$transcripts}){
            #if (defined($tx) && $tx->biotype eq 'protein_coding'){
            if (defined($tx)){
                $enst = $tx->stable_id();
                my $tr = $tx->translation();
                if(defined($tr)){
                    $ensp = $tr->stable_id;
                    my $peptide = $tr->seq();
                    if(defined($peptide)){
                        my $transcript = $tx->transform('chromosome'); 
                        my $tm = Bio::EnsEMBL::TranscriptMapper->new($transcript);
                        my $strand = $tx->strand;
                        my @pep_coords = $tm->genomic2pep($start, $end, $tx->strand);
                        foreach my $pep_coord(@pep_coords){
                            if ($pep_coord->start != $start && $pep_coord->end != $end){
                                my $peptide_start = $pep_coord->start;
                                my $peptide_end = $pep_coord->end;
                                my $marked_peptide;
                                if($peptide_start <= length($peptide) && $peptide_end < length($peptide)){
                                    my $pep_left = substr($peptide,0,$peptide_start-1);
                                    my $pep_mid = substr($peptide,$peptide_start-1,$peptide_end-$peptide_start+1);
                                    my $pep_right = substr($peptide,$peptide_end);
                                    $marked_peptide = $pep_left."[".$pep_mid."]".$pep_right;
                                }
                                else{
                                    $marked_peptide = $peptide.'[]';
                                }
                                printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%s\n", $coord, $chrom, $start, $end, $ref, $ensg, $enst, $ensp, $hgnc, $tx->is_canonical(), length($peptide), $peptide_start, $peptide_end, $marked_peptide);
                            }
                        }
                    }
                }
            }
        }
    }
}
