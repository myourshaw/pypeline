#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FileHandle;
use Bio::EnsEMBL::Registry;

#--host cortex.local --user ensembl --password ensembl

        #Ensembl version
        # On 06/10/2014 19:14, Michael Yourshaw wrote:
        # Is there a programatic way to get:
        #  1. The version number of the current Ensembl?
        #  2. The genome build for a species in the current version?
        #  3. The genome build for a species in any specified version number?
        #  4. The most recent version for a given species/genome build?
        # Solutions with either mysql or the perl API will do.
        #
        # Hi Michael,
        #
        # Most of this information is available via the meta table.
        #
        # The API will allow you to fetch this information with the following queries:
        # my $meta_container = $registry->get_adaptor('human', 'core', 'MetaContainer');
        # my $api_version = $meta_container->get_schema_version();
        # my $genebuild_version = $meta_container->get_genebuild();
        # my $coord_system_adaptor = $registry->get_adaptor('human', 'core', 'CoordSystem');
        # ###my $genome_version = $coord_system_adaptor->fetch_top_level();
        # my $genome_version = $coord_system_adaptor->fetch_by_name('chromosome')->version;
        #
        # It is however not possible to retrieve information about previous builds from the current database.
        # To do this, you would need to check the genebuild_version, which indicates when the genebuild was last updated.
        # Based on the data, find the archive prior to that date.
        #
        # It is also possible to retrieve this information from the REST API.
        #
        # http://rest.ensembl.org/info/software?content-type=application/json
        # will return the current Ensembl version
        #
        # http://rest.ensembl.org/info/assembly/homo_sapiens?content-type=application/json
        # will return basic information about the current genebuild and assembly for a given species.
        
my $opt = {};
GetOptions(
    $opt,
    'help',                    # displays help message
    'host=s',                  # database host
    'port=i',                  # database port
    'user=s',                  # database user name
    'species=s',               # species e.g. human, homo_sapiens
    'db_version=i',            # database version e.g. 77
    'password=s',              # database password
    'output_file|o=s',         # output file name
    'debug!'                    #debug mode, print alternate db info
) or die "ERROR: Failed to parse command-line flags\n";
# set defaults
$opt->{host} ||= 'ensembldb.ensembl.org'; #cortex.local
#$opt->{host} ||= ' useastdb.ensembl.org'; #current and previous versions only
$opt->{port} ||= 3306;
$opt->{user} ||= 'anonymous'; #ensembl
$opt->{password} ||= ''; #ensembl
$opt->{species} ||= 'homo_sapiens';
$opt->{db_version} ||= 0;
$opt->{output_file} ||= 'STDOUT';
$opt->{debug} ||= 0;

my $usage =<<END;
Returns current Ensembl version and genome version (tab-separated) for installed Ensembl Perl API, by default in STDOUT
Usage:
perl ensembl_version.pl [options]

Options
=======

--help                 Display this message and quit
-o | --output_file     Output file [default: "STDOUT"]
--host                 Manually define database host [default: "ensembldb.ensembl.org"]
--port                 Database port [default: "3337"]
--user                 Database username [default: "anonymous"]
--password             Database password [default: ""]
--species [species]    Species to use [default: "homo sapiens"]
--db_version [db_version]    Database version to use [default: ""]
--debug                Print additional debug info
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

my $reg = 'Bio::EnsEMBL::Registry';
if ($opt->{db_version}) {
    $reg->load_registry_from_db(
        -host => $opt->{host},
        -port => $opt->{port},
        -user => $opt->{user},
        -pass => $opt->{password},
        -species => $opt->{species},
        -db_version => $opt->{db_version},
    );
}
else{
    $reg->load_registry_from_db(
        -host => $opt->{host},
        -port => $opt->{port},
        -user => $opt->{user},
        -pass => $opt->{password},
        -species => $opt->{species},
    );
}


my $core_mca = $reg->get_adaptor($opt->{species}, 'core', 'metacontainer');
my $api_version = $core_mca->get_schema_version;
my $coord_system_adaptor = $reg->get_adaptor($opt->{species}, 'core', 'CoordSystem');
my $genome_version = $coord_system_adaptor->fetch_by_name('chromosome')->version;
print $out_file_handle "$api_version\t$genome_version";

if ($opt->{debug}) {

    my %db_params =
    (
        cortex_human => {(host => 'cortex.local', port => 3306, user => 'ensembl', pass => 'ensembl', species => 'homo_sapiens',)},
        cortex_mouse => {(host => 'cortex.local', port => 3306, user => 'ensembl', pass => 'ensembl', species => 'mus_musculus',)},
        ensembl_3306_human => {(host => 'ensembldb.ensembl.org', port => 3306, user => 'anonymous', pass => '', species => 'homo_sapiens',)},
        ensembl_3306_mouse => {(host => 'ensembldb.ensembl.org', port => 3306, user => 'anonymous', pass => '', species => 'mus_musculus',)},
    );
    if ($opt->{db_version}) {
        $db_params{requested} = {(host => $opt->{host}, port => $opt->{port}, user => $opt->{user}, pass => $opt->{password}, species => $opt->{species}, db_version => $opt->{db_version},)};
    }
    else{
        $db_params{requested} = {(host => $opt->{host}, port => $opt->{port}, user => $opt->{user}, pass => $opt->{password}, species => $opt->{species}, db_version => $opt->{db_version},)};
    }
    
    print "\nhost\tport\tuser\tpassword\tspecies\tapi_version\tgenebuild_version\tgenome_version\n";
    
    foreach my $db(keys %db_params) {
        my %p = %{$db_params{$db}};
        my $registry = 'Bio::EnsEMBL::Registry';
        $registry->load_registry_from_db(
           -host => $p{host},
           -port => $p{port},
           -user => $p{user},
           -pass => $p{pass},
           -species => $p{species},
           );
        my $meta_container = $registry->get_adaptor($p{species}, 'core', 'MetaContainer');
        my $api_version = $meta_container->get_schema_version();
        my $genebuild_version = $meta_container->get_genebuild();
        my $coord_system_adaptor = $registry->get_adaptor($p{species}, 'core', 'CoordSystem');
        my $genome_version = $coord_system_adaptor->fetch_by_name('chromosome')->version;
        
        print "$p{host}\t$p{port}\t$p{user}\t$p{pass}\t$p{species}\t$api_version\t$genebuild_version\t$genome_version\n";
    }
}

1;
