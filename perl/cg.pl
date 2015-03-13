#!/usr/bin/perl -w
use strict;
use warnings;

#BioPerl
use Bio::Seq;
use Bio::SeqIO;
#EnsEMBL Perl API
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::Utils::Sequence qw (sequence_with_ambiguity);

#MYO1D amplicons
#chr17:27843618-27846171 chr17:27950242-27956650 chr17:27989732-27990114 chr17:28004817-28010621 chr17:28062962-28072425 chr17:28089194-28089653 chr17:28096139-28096546 chr17:28099897-28100264 chr17:28106456-28106919 chr17:28111295-28118983 chr17:28122096-28129875 chr17:28131584-28132048

use constant USAGE => "cg computes CG content over a range of sequence(s) and a specified window.\n"
	 . "usage:\n"
	 . "perl cg.pl window <range list>\n"
	 . "<range list> = chr17:27843618-27846171 chr17:27950242-27956650 ...\n"
;
unless (@ARGV || ($ARGV[0] ne 'help' && $ARGV[0] ne '?')) {
	print STDERR USAGE;
	exit;
}

#connect to EnsEMBL database
my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db(
	-host => 'ensembldb.ensembl.org',
	-user => 'anonymous');
#my $exon_adaptor = $registry->get_adaptor('human', 'core', 'exon');
#my $gene_adaptor = $registry->get_adaptor('human', 'core', 'gene');
#my $vega_adaptor = $registry->get_adaptor('human', 'vega', 'gene');
#my $repeat_adaptor = $registry->get_adaptor('human', 'core', 'repeatfeature');
my $slice_adaptor = $registry->get_adaptor('human', 'core', 'slice');
#my $transcript_adaptor = $registry->get_adaptor('human', 'core', 'transcript');
#my @db_adaptors = @{ $registry->get_all_DBAdaptors() };
my $window = shift;
print "chr\tpos\tbase\tCG($window)\n";
while(my $range = shift){
	if ($range =~ /^\s*(?:chr)?[\s]*([0-9XYMT]{1,2})[\s:]+([,\d]+)\D+([,\d]+)(?:\s+([01+-]{1}))?/i){
		if ($1 && $2 && $3) {
			my $chr = $1;
			my $start = $2;
			my $end = $3;
			my $strand = defined($4)
			 && $4 == "+" ? 1 : defined($4)
			 && $4 == "-" ? -1 : $4;
			$start =~ s/,//g;
			$end =~ s/,//g;
			$start = $start <= $end ? $start : $end;
			$end = $end >= $start ? $end : $start;
			my $sliceStart = $start-$window/2;
			my $sliceEnd = $end+$window/2;
			my $slice = $slice_adaptor->fetch_by_region('chromosome', $chr, $sliceStart>0?$sliceStart:1, $sliceEnd, $strand);
			if($slice){
				my ($lpad,$rpad)=('','');
				#pad in case out of bounds
				if($sliceStart<$slice->start){$lpad=reverse(substr($slice->seq,0,$slice->start-$sliceStart))}
				if($sliceEnd>$slice->end){$rpad=reverse(substr($slice->seq,0,$slice->end-$sliceEnd))}
				my $seq = $lpad.$slice->seq.$rpad;
				#print "$range window: $window\n";
				#print "chr\tpos\tbase\tCG\n";
				for (my $i=$window/2; $i <= length($seq)-$window/2-1; $i++){
					my $subseq = substr($seq, $i-$window/2, $window);
					my $cg = $subseq =~ tr/cgCG/cgCG/;
					my $content = $cg / $window;
					printf("%s\t%u\t%s\t%f\n",$chr,$start+$i-$window/2,substr($seq,$i+1*$window%2,1),$content);
				}
			}
		}
	}
}
