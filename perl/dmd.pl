use strict;
use warnings;
use Getopt::Std;
use File::Spec;
use Date::Calc qw(Today_and_Now Delta_YMDHMS);
use POSIX qw(ceil floor);
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
#BioPerl
use Bio::Seq;
use Bio::SeqIO;
#EnsEMBL Perl API
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::Utils::Sequence qw (sequence_with_ambiguity);

use constant PRIMER3 => "0";
use constant GENE =>  "DMD";
use constant ENSMBL_GENE_ID => "ENSG00000198947";
use constant VEGA_GENE_ID => "OTTHUMG00000021336";

#connect to EnsEMBL database
my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db(
	-host => 'ensembldb.ensembl.org',
	-user => 'anonymous');
my $exon_adaptor = $registry->get_adaptor('human', 'core', 'exon');
my $gene_adaptor = $registry->get_adaptor('human', 'core', 'gene');
my $vega_adaptor = $registry->get_adaptor('human', 'vega', 'gene');
my $repeat_adaptor = $registry->get_adaptor('human', 'core', 'repeatfeature');
my $slice_adaptor = $registry->get_adaptor('human', 'core', 'slice');
my $transcript_adaptor = $registry->get_adaptor('human', 'core', 'transcript');
my @db_adaptors = @{ $registry->get_all_DBAdaptors() };

my $id = ENSMBL_GENE_ID;
my $gene = $gene_adaptor->fetch_by_stable_id($id);
if ($gene) {
	my %exon;
	my @gene_exons = @{$gene->get_all_Exons};
	if($gene->strand==1){
		@gene_exons = sort {$a->{start} <=> $b->{start} || $a->{end} <=> $b->{end}} @gene_exons;
	}
	else{
		@gene_exons = sort {$b->{start} <=> $a->{start} || $b->{end} <=> $a->{end}} @gene_exons;
	}
	for(my$k=0;$k<@gene_exons;$k++){
		$exon{$gene_exons[$k]->stable_id}=$k+1;
	}
    my @transcripts = @{$gene->get_all_Transcripts};
	if($gene->strand==1){
		@transcripts = sort {$a->{start} <=> $b->{start} || $a->{end} <=> $b->{end}} @transcripts;
	}
	else{
		@transcripts = sort {$b->{start} <=> $a->{start} || $b->{end} <=> $a->{end}} @transcripts;
	}
	if(PRIMER3){
		for(my $j=0;$j<@gene_exons;$j++){
			#my $ex = $gene_exons[$j];
			my $slice = $slice_adaptor->fetch_by_region('chromosome', $gene_exons[$j]->slice->seq_region_name,$gene_exons[$j]->start,$gene_exons[$j]->end,$gene_exons[$j]->strand);
			printf "PRIMER_SEQUENCE_ID=%s(%u)\n",$gene_exons[$j]->stable_id,abs($slice->end-$slice->start)+1;
			printf "SEQUENCE=%s\n=\n",$slice->seq;
		}
	}
	else{
		print "name\tid\tchr\tstart\tend\tlength\texons\tstrand\tbiotype\n";
		printf "%s\t%s\t%s\t%u\t%u\t%u\t%u\t%s\t%s\n",$gene->external_name,$gene->stable_id,$gene->slice->seq_region_name,$gene->start,$gene->end,$gene->end-$gene->start+1,$#gene_exons+1,$gene->strand==1?"+":"-",$gene->biotype;
		for(my $j=0;$j<@gene_exons;$j++){
			my $ex = $gene_exons[$j];
			printf "exon_%u\t%s\t%s\t%u\t%u\t%u\t\t%s\t%s\n",$exon{$ex->stable_id},$ex->stable_id,$ex->slice->seq_region_name,$ex->start,$ex->end,$ex->end-$ex->start+1,$ex->strand==1?"+":"-","";#,$ex->biotype;
		}
		for(my $i=0;$i<@transcripts;$i++){
			my $tr = $transcripts[$i];
			my @exons = @{$tr->get_all_Exons};
			if($gene->strand==1){
				@exons = sort {$a->{start} <=> $b->{start} || $a->{end} <=> $b->{end}} @exons;
			}
			else{
				@exons = sort {$b->{start} <=> $a->{start} || $b->{end} <=> $a->{end}} @exons;
			}
			printf "transcript_%u\t%s\t%s\t%u\t%u\t%u\t%u\t%s\t%s\n",$i+1,$tr->stable_id,$tr->slice->seq_region_name,$tr->start,$tr->end,$tr->end-$tr->start+1,$#exons+1,$tr->strand==1?"+":"-",$tr->biotype;
			for(my $j=0;$j<@exons;$j++){
				my $ex = $exons[$j];
				printf "exon_%u\t%s\t%s\t%u\t%u\t%u\t\t%s\t%s\n",$exon{$ex->stable_id},$ex->stable_id,$ex->slice->seq_region_name,$ex->start,$ex->end,$ex->end-$ex->start+1,$ex->strand==1?"+":"-","";#,$ex->biotype;
			}
		}
	}
}
else {
	printf STDERR "%s isn't a gene\n", $id;
}
