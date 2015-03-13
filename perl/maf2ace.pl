#!/usr/bin/perl -w
use strict;
use warnings;
use File::Spec;
use File::Basename;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use IO::Compress::Gzip qw(gzip $GzipError);
require File::Temp;
use File::Temp ();
use File::Temp qw/ :seekable /;

use Bio::SeqIO;
use Bio::AlignIO;

#test data
#"C:\Users\Michael\solexa\solexa_datasets\Reports\Mike_Yourshaw_targeted_sequencing\09-11-08\080108_HWI-EAS172_200G8_s_2_all_seq.fasta.bfast\bfast\targeted_chr17_27843500-28132200\bfast.report.080108_HWI-EAS172_200G8_s_2_all_seq.fasta.maf.gz"
#"C:\Users\Michael\solexa\solexa_datasets\Reports\Mike_Yourshaw_targeted_sequencing\09-11-08\080108_HWI-EAS172_200G8_s_2_all_seq.fasta.bfast\bfast\targeted_chr17_27843500-28132200\bfast.report.080108_HWI-EAS172_200G8_s_2_all_seq.fasta.maf"

my $file = shift;
my ($maf_file,$dir,$ext) = fileparse($file,'\.maf');


if($ext eq ".gz"){
	my $fh = new File::Temp() or croak();
	my $fname = $fh->filename;
	gunzip $file => $fh or die "gunzip failed for $file: $GunzipError\n";
	$fh->seek(0,SEEK_SET);
}

my $maf = Bio::AlignIO->new(-file => $file);
my $ace_file = File::Spec->catfile($dir,$maf_file.'.ace');
my $ace = Bio::AlignIO->new('-file' => ">$ace_file", '-format' => 'mase');
while (my $aln = $maf->next_aln()){
    $ace->write_aln($aln);

}

exit;
my $seqio_object = Bio::SeqIO->new(-file => $file);
my $seq_object = $seqio_object->next_seq;
