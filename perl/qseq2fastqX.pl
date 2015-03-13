#!perl -w
use strict;
use File::Glob;




my @run_ids = ("090813_HWUSI-EAS172_42MEGb");
my @lanes = (1,2);
my @reads = (1,2);
my @read_lengths = (76,76);
my $barcode_length = 5;
my @barcode_trim = (0,1); #characters to drop from barcode (left,right)
my @valid_barcodes = ();

my $qseq_record_regex = '(?<machine>\S+)\s+(?<run>\S+)\s+(?<lane>\d+)\s+(?<tile>\d+)\s+(?<x>\d+)\s+(?<y>\d+)\s+(?<index>\d+)\s+(?<read>[12]{1})\s+(?<sequence>[ACGT.]+)\s+(?<quality>\S+)\s+(?<filter>[01]{1})';
foreach my $run_id(@run_ids){
    foreach my $lane(@lanes){
        foreach my $read(@reads){
            my $qseq_glob = "~/Lab/solexa_datasets/$run_id/s_$lane"."_$read"."_????_qseq.txt";
            #my $qseq_dir = "/home/solexa/solexa_datasets/$run_id/Data/Intensities/BaseCalls/s_?_?_????_qseq.txt";
            my $fastq_dir = "~/solexa_datasets/$run_id/fastq";
            foreach my $path(glob($qseq_glob)){
                open IN,"<",$path;
                my $record_count = 0;
                while(<IN>){
                    $record_count++;
                    if(/$qseq_record_regex/){
                        my %qqq = %+;
                        $+{lane} eq $lane or die "record: $record_count, file lane: $lane != record lane: $+{lane}";
                        $+{read} eq $read or die "record: $record_count, file read: $read != record read: $+{read}";
                        my $barcode = substr($+{sequence},$barcode_trim[0],$barcode_length-$barcode_trim[1]);
                        my $q = substr($+{quality},$barcode_length-1);
                        my $q1 = $+{quality};
                        my $barcode_qualty = substr($+{quality},$barcode_trim[0],$barcode_length-$barcode_trim[1]);
                        (my $sequence = substr($+{sequence},$barcode_length)) =~ s/\./N/g;
                        my $s = $+{sequence};
                        my $quality = substr($+{quality},$barcode_length);
                       printf OUT "@%s:%s:%s:%s:%s:%s:%s:%s:%s:%s\n%s:",($run_id,$lane,$read,$barcode,$barcode_qualty,$+{tile},$+{x},$+{y},$+{index},$+{filter},$sequence);
                    }
                }
            }
        }
    }
}
my ($machine,$run,$lane,$tile,$x,$y,$index,$read,$sequence,$quality,$filter) = ();

#print length("CAGATCTGAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATATCGTATTTCGTCATCTTCTTGAAACACA");
#HWUSI-EASXXX	1	1	75	23	914	0	1	CAGATCAGTTGTTGCAGTGGAATAGTCAGGTTAAATTTAATGTGACCGTTTATCGCAATCTGCCGAACACTCGCGA	ababbbbabbbbbab_ba`abbbbbbba`aabb``abbbbbbbba_^`^`bb`X_]^abb`U__]ND^]\]`Z]]T	1
my $qseq_file_regex = 's_[12]_[12]_\d{4}_qseq\.txt';
my @qseq_files;
foreach my $path(@qseq_files){
}
print "done\n";
while (<>) {
	chomp;
	my @items = split /\t/;
	print "@","$items[0]:$items[2]:$items[3]:$items[4]:$items[5]#$items[6]/$items[7]\n";
	print "$items[8]\n";
	print "+\n";
	print "$items[9]\n";
}
