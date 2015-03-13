#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FileHandle;

my $opt = {};
GetOptions(
    $opt,
    'help',                    # displays help message
    'species=s',               # species e.g. human, homo_sapiens
    'host=s',                  # database host
    'port=s',                  # database port
    'user=s',                  # database user name
    'password=s',              # database password
    'mode|m=s',                  # all or canonical
    'output_file|o=s',         # output file name
    'amino_acids=s',             #comma-separated list of pos=aa specifications, e.g., 349=Y,778=E,674=R
) or die "ERROR: Failed to parse command-line flags\n";
# set defaults
$opt->{species} ||= "homo_sapiens";
$opt->{host} ||= 'cortex.local';
$opt->{port} ||= 3306;
$opt->{user} ||= 'ensembl';
$opt->{password} ||= 'ensembl';
$opt->{mode} ||= 'canonical';
$opt->{output_file} ||= "STDOUT";

    my $usage =<<END;
Usage:
perl protein_lengths.pl [options]

Options
=======

--help                 Display this message and quit
-o | --output_file     Output file [default: "STDOUT"]
-m | --mode            Mode ("all" for all transcripts,
                        "canonical" for canonical transcript and longest protein
                        if canonical is not longest) [default: "canonical"] 
--species [species]    Species to use [default: "homo_sapiens"]
--host                 Manually define database host [default: "cortex.local"]
--port                 Database port [default: "3306"]
--user                 Database username [default: "ensembl"]
--password             Database password [default: "ensembl"]
--amino_acids          comma-separated list of pos=aa specifications, e.g., 349=Y,778=E,674=R
END

if ($opt->{help}){
    print $usage;
    exit(0);
}


my $out_file_handle = new FileHandle;
if($opt->{output_file} =~ /stdout/i) {
    $out_file_handle = *STDOUT;
}
else {
    $out_file_handle->open(">".$opt->{output_file}) or die("ERROR: Could not write to output file ", $opt->{output_file}, "\n");
}

my %aa = ();
my @aas = split /,/,$opt->{amino_acids};
for my $a(@aas){
    my @pos_aa = split /=/, $a;
    $aa{$pos_aa[0]} = $pos_aa[1]; 
}

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

my $start_time = get_time();
print $out_file_handle "## Protein lengths run started $start_time\n";

use Bio::EnsEMBL::Registry;

my $reg = 'Bio::EnsEMBL::Registry';

$reg->load_registry_from_db(
    -host => $opt->{host},
    -port => $opt->{port},
    -user => $opt->{user},
    -pass => $opt->{password}
);

my $core_mca = $reg->get_adaptor('homo_sapiens', 'core', 'metacontainer');
my $schema_version = $core_mca->get_schema_version;
print $out_file_handle "## Ensembl database schema version: $schema_version\n";

my $ga = $reg->get_adaptor('homo_sapiens', 'core', 'gene');

my @genes = @{$ga->fetch_all()};

print $out_file_handle "#ensg\tenst\tensp\thgnc\tprotein_length\tcanonical\n";

while (my $gene = shift @genes) {
    my $ensg = $gene->stable_id();
    $gene->external_db('HGNC');
    my $hgnc = $gene->external_name() || '';
    if ($opt->{mode} eq 'all'){
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
                    my $protein = $tr->seq();
                    if(defined($protein)){
                        my $found = 1;
                        while ((my ($pos, $a)) = each (%aa)){
                            if((length($protein) < $pos) || (substr($protein, $pos-1, 1) ne $a)){
                                $found = 0;
                                last
                            }
                        }
                        if ($found == 1){
                            printf $out_file_handle "%s\t%s\t%s\t%s\t%d\t%d\n", $ensg, $enst, $ensp, $hgnc, length($protein),$tx->is_canonical();
                        }
                        else{
                            
                        }
                    }
                }
            }
        }
    }
    else{
        #get canonical transcript
        my $canonicaltx = $gene->canonical_transcript();
        my $canonicalprotein;
        if (defined($canonicaltx) && $canonicaltx->biotype eq 'protein_coding'){
            my $enst = $canonicaltx->stable_id();
            my $tr = $canonicaltx->translation();
            if(defined($tr)){
                my $ensp = $tr->stable_id;
                $canonicalprotein = $tr->seq();
                if(defined($canonicalprotein)){
                    my $found = 1;
                    #if($ensp eq 'ENSP00000354876'){
                    #    print;
                    #}
                    foreach my $pos (sort keys %aa){
                        if((length($canonicalprotein) < $pos) || (substr($canonicalprotein, $pos-1, 1) ne $aa{$pos})){
                            $found = 0;
                            last;
                        }
                    }
                    if ($found == 1){
                        printf $out_file_handle "%s\t%s\t%s\t%s\t%d\t1\n", $ensg, $enst, $ensp, $hgnc, length($canonicalprotein);
                        print "$canonicalprotein\n";
                        print length($canonicalprotein)."\n";
                        print substr($canonicalprotein, 349, 1)."\n";
                        print substr($canonicalprotein, 674, 1)."\n";
                        print substr($canonicalprotein, 778, 1)."\n";
                    }
                }
            }
        }
    }
}

