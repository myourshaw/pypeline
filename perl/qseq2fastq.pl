#!perl -w
use strict;


# Author: Michael Yourshaw

use strict;
use warnings;
use Getopt::Std;

select STDERR; $| = 1;  # make unbuffered

my %opts;
my $version = '0.1.1';
my $usage = qq{
Usage: qseq2fastq.pl <input control-file>

	This script will convert Illumina output files (*qseq files)
	to the BFAST fastq multi-end format.  For single-end reads 
	that do not have more than one end (neither mate-end nor
	paired end), the output format is strictly fastq format.
	For multi-end data (mate-end, paired end ect.) each end is
	outputted consecutively in 5'->3' ordering in fastq blocks.
	Therefore, the output reamins strictly conforming to the 
	fastq format, but keeps intact the relationships between
	sequences that all originate from the same peice of DNA.
	We assume for paired end data that all data is able to be
	paired.

	All .qseq files will be inferred from the QSEQ_DIR, LANES, and  READS
        parameters in the control file. qseq file names must be in the form
            s_<lane>_<read>_<tile>_qseq.txt
        
        Control file parameters (case-insensitive, terminated by # or end of line,
                                lists are comma delimited)
            QSEQ_DIR #directory containing qseq files, default current working directory
            SAMPLE_ID #library sample id
            RUN_ID #default is machine_time from qseq record
            LANES #list of lanes to be converted and merged
            READS=1,2 #paired end/mate pair read ids, default 1,2
            READ_LENGTHS=76,76 #number of bases in read, same order as READS
            BARCODE_LENGTH=5 #length of barcode at 5'end, default 0 (no barcode)
            BARCODE_TRIM=0,1 #number of bases to trim from barcode (5',3')
            VALID_BARCODES=AATA,ACGA,AAGC,TGAG #list of valid trimmed barcodes
            GO=QSEQ2FASTQ #run program with current parameters, default QSEQ2FASTQ
};
my %params = (
    QSEQ_DIR => '.',
    SAMPLE_ID => '',
    RUN_ID => '',
    LANES => '', 
    READS => '1,2',
    READ_LENGTHS => '',
    BARCODE_LENGTH => 0,
    BARCODE_TRIM =>'0,1',
    VALID_BARCODES => '',
    GO => 'QSEQ2FASTQ'
);
getopts ('', \%opts);
die($usage) if (@ARGV < 1);

while(<>){
    chomp;
    if(!/^\s*#/ && /\s*([^\s=#]+){1}[\s=#]+([^#]+)#*$/){
        my $key = uc($1);
        my $value = $2;
        $value =~ s/^['"]{1}|['"]{1}$//g;
        if($key eq "GO" && $value eq "QSEQ2FASTQ"){
            qseq2fastq();
        }
        elsif (exists($params{$key})) {
            print STDERR "$key=$value\n";
            $params{$key} = $value;
        }
        else {
            print STDERR "WARNING: invalid parameter $key=$value ignored\n";
        }
    }
}

sub qseq2fastq{
    if(exists $params{QSEQ_DIR} && -d $params{QSEQ_DIR}){
	my @lane_file_list = ();
        my $l = 0;
        foreach my $lane (split /,/, $params{LANES}){
	    my @read_file_list = ();
	    my $r = 0;
            foreach my $read (split /,/, $params{READS}){
		my $file_regex = sprintf("s_%s_%s_.*_qseq\.txt",$lane,$read);
		my @files = ();
		get_files($params{QSEQ_DIR},$file_regex,\@files);
		if(($r == 0 && scalar @files > 0) || scalar @files == @{$read_file_list[$r-1]}){
		    $read_file_list[$r++] = \@files;
		 }
		 else{
		     die "No files or wrong number of files for lane $lane read $read in $params{QSEQ_DIR}\n";
		 }                
            }
            $lane_file_list[$l++] = \@read_file_list;
        }
	for(my $l=0;$l<@lane_file_list;$l++){ #lane
	    my @in_filehandles = ();
	    my @out_filehandles = ();
	    #for(my $t=0;$t<@{$lane_file_list[0]};$t++){ #tile
		for(my $r=0;$r<@{$lane_file_list[$l]};$r++){ #read
		    open my $in_fh,"<",${$lane_file_list[$l]}[$r];
		    push @in_filehandles,<$in_fh>;
		    open my $out_fh,"<","${$lane_file_list[$l]}[$r].fastq";
		    push @out_filehandles,<$out_fh>;
		}
	    #}
	    #for(my $t=0;$t<@{$lane_file_list[0]};$t++){ #tile
		for(my $f=0;$f<@in_filehandles;$f++){
		    while(<$in_filehandles[$f]>){
			print;
		    }
		}
	    #}
	}
	foreach my $lane_files(@lane_file_list){
	    foreach my $read_file(@{$lane_files}){
	    }
	    print;
	}
        return 1;
    }
    else{
        print STDERR "$params{QSEQ_DIR} doesn't exist\n";
        return 0;
    }
}

sub get_files {
	my $dir = shift;
	my $file_regex = shift;
	my $files = shift;
	local *DIR;
	opendir(DIR, "$dir/") or die("Error.  Can't open $dir\n");
	@$files = grep !/^\.\.?$/, readdir DIR;
	@$files = grep /^$file_regex$/, @$files;
	for(my $i=0;$i<scalar(@$files);$i++) {
		@$files[$i] = $dir."/".@$files[$i];
	}
	close(DIR);
}
