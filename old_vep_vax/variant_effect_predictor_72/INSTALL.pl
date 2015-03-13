#!/usr/bin/perl

use Getopt::Long;
use LWP::Simple qw($ua getstore get);
use File::Listing qw(parse_dir);
use File::Path;
use File::Copy;
use Archive::Tar;
use Net::FTP;
use Cwd;

$| = 1;
our $VERSION = 72;

# CONFIGURE
###########

my ($DEST_DIR, $ENS_CVS_ROOT, $API_VERSION, $BIOPERL_URL, $CACHE_URL, $FTP_USER, $help);

GetOptions(
  'DESTDIR|d=s'  => \$DEST_DIR,
  'VERSION|v=i'  => \$API_VERSION,
  'BIOPERL|b=s'  => \$BIOPERL_URL,
  'CACHEURL|u=s' => \$CACHE_URL,
  'CACHEDIR|c=s' => \$CACHE_DIR,
  'HELP|h'       => \$help
);

if(defined($help)) {
  usage();
  exit(0);
}

my $default_dir_used;

# check if $DEST_DIR is default
if(defined($DEST_DIR)) {
  print "Using non-default installation directory $DEST_DIR - you will probably need to add $DEST_DIR to your PERL5LIB\n";
  $default_dir_used = 0;
}
else {
  $DEST_DIR ||= '.';
  $default_dir_used = 1;
}

my $lib_dir = $DEST_DIR;

$DEST_DIR       .= '/Bio';
$ENS_CVS_ROOT ||= 'http://cvs.sanger.ac.uk/cgi-bin/viewvc.cgi/';
$BIOPERL_URL  ||= 'http://bioperl.org/DIST/BioPerl-1.6.1.tar.gz';
$API_VERSION  ||= 72;
$CACHE_URL    ||= "ftp://ftp.ensembl.org/pub/release-$API_VERSION/variation/VEP";
$CACHE_DIR    ||= $ENV{HOME} ? $ENV{HOME}.'/.vep' : 'cache';
$FTP_USER     ||= 'anonymous';

our $prev_progress = 0;

print "\nHello! This installer is configured to install v$API_VERSION of the Ensembl API for use by the VEP.\nIt will not affect any existing installations of the Ensembl API that you may have.\n\nIt will also download and install cache files from Ensembl's FTP server.\n\n";


# CHECK EXISTING
################

print "Checking for installed versions of the Ensembl API...";

# test if the user has the API installed
my $has_api = {
  'ensembl' => 0,
  'ensembl-variation' => 0,
  'ensembl-functgenomics' => 0,
};

eval q{
  use Bio::EnsEMBL::Registry;
};

my $installed_version;

unless($@) {
  $has_api->{ensembl} = 1;
  
  $installed_version = Bio::EnsEMBL::Registry->software_version;
}

eval q{
  use Bio::EnsEMBL::Variation::Utils::VEP;
};

$has_api->{'ensembl-variation'} = 1 unless $@;

eval q{
  use Bio::EnsEMBL::Funcgen::RegulatoryFeature;
};

$has_api->{'ensembl-functgenomics'} = 1 unless $@;


print "done\n";

my $total = 0;
$total += $_ for values %$has_api;

my $message;

if($total == 3) {
  
  if(defined($installed_version)) {
    if($installed_version == $API_VERSION) {
      $message = "It looks like you already have v$API_VERSION of the API installed.\nYou shouldn't need to install the API";
    }
    
    elsif($installed_version > $API_VERSION) {
      $message = "It looks like this installer is for an older distribution of the API than you already have";
    }
    
    else {
      $message = "It looks like you have an older version (v$installed_version) of the API installed.\nThis installer will install a limited set of the API v$API_VERSION for use by the VEP only";
    }
  }
  
  else {
    $message = "It looks like you have an unidentified version of the API installed.\nThis installer will install a limited set of the API v$API_VERSION for use by the VEP only"
  }
}

elsif($total > 0) {
  $message = "It looks like you already have the following API modules installed:\n\n".(join "\n", grep {$has_api->{$_}} keys %$has_api)."\n\nThe VEP requires the ensembl, ensembl-variation and optionally ensembl-functgenomics modules";
}

if(defined($message)) {
  print "$message\n\nSkip to the next step (n) to install cache files\n\nDo you want to continue installing the API (y/n)? ";
  
  my $ok = <>;
  
  if($ok !~ /^y/i) {
    print " - skipping API installation\n";
    goto CACHE;
  }
}



# SETUP
#######

print "\nSetting up directories\n";

# check if install dir exists
if(-e $DEST_DIR) {
  print "Destination directory $DEST_DIR already exists.\nDo you want to overwrite it (if updating VEP this is probably OK) (y/n)? ";
  
  my $ok = <>;
  
  if($ok !~ /^y/i) {
    print "Exiting\n";
    exit(0);
  }
  
  else {
    unless($default_dir_used) {
      print "WARNING: You are using a non-default install directory.\nPressing \"y\" again will remove $DEST_DIR and its contents!!!\nAre you really, really sure (y/n)? ";
      $ok = <>;
      
      if($ok !~ /^y/i) {
        print "Exiting\n";
        exit(0);
      }
    }
    
    # try to delete the existing dir
    rmtree($DEST_DIR) or die "ERROR: Could not delete directory $DEST_DIR\n";
  }
}

mkdir($DEST_DIR) or die "ERROR: Could not make directory $DEST_DIR\n";
mkdir($DEST_DIR.'/tmp') or die "ERROR: Could not make directory $DEST_DIR/tmp\n";

# set up a user agent's proxy
$ua->env_proxy;

# enable progress
eval q{
  $ua->show_progress(1);
};



# API
#####

print "\nDownloading required files\n";

# set up the URLs
my $ensembl_url_tail = '.tar.gz?root=ensembl&view=tar&only_with_tag=branch-ensembl-';

foreach my $module(qw(ensembl ensembl-variation ensembl-functgenomics)) {
  my $url = $ENS_CVS_ROOT.$module.$ensembl_url_tail.$API_VERSION;
  
  print " - fetching $module\n";
  
  my $target_file = $DEST_DIR.'/tmp/'.$module.'.tar.gz';
  
  if(!-e $target_file) {
    unless(getstore($url, $target_file) == 200) {
      die "ERROR: Failed to fetch $module from $url - perhaps you have a proxy/firewall? Set the http_proxy ENV variable if you do\nError code: $response\n";
    }
  }
  
  print " - unpacking $target_file\n";
  unpack_tar("$DEST_DIR/tmp/$module.tar.gz", "$DEST_DIR/tmp/");
  
  print " - moving files\n";
  
  if($module eq 'ensembl') {
    move("$DEST_DIR/tmp/$module/modules/Bio/EnsEMBL", "$DEST_DIR/EnsEMBL") or die "ERROR: Could not move directory\n".$!;
  }
  elsif($module eq 'ensembl-variation') {
    move("$DEST_DIR/tmp/$module/modules/Bio/EnsEMBL/Variation", "$DEST_DIR/EnsEMBL/Variation") or die "ERROR: Could not move directory\n".$!;
  }
  elsif($module eq 'ensembl-functgenomics') {
    move("$DEST_DIR/tmp/$module/modules/Bio/EnsEMBL/Funcgen", "$DEST_DIR/EnsEMBL/Funcgen") or die "ERROR: Could not move directory\n".$!;
  }
  
  rmtree("$DEST_DIR/tmp/$module") or die "ERROR: Failed to remove directory $DEST_DIR/tmp/$module\n";
}



# BIOPERL
#########

# now get BioPerl
print " - fetching BioPerl\n";

$bioperl_file = (split /\//, $BIOPERL_URL)[-1];

my $target_file = $DEST_DIR.'/tmp/'.$bioperl_file;

unless(getstore($BIOPERL_URL, $target_file) == 200) {
  die "ERROR: Failed to fetch BioPerl from $BIOPERL_URL - perhaps you have a proxy/firewall?\nError code: $response\n";
}

print " - unpacking $target_file\n";
unpack_tar("$DEST_DIR/tmp/$bioperl_file", "$DEST_DIR/tmp/");

print " - moving files\n";

$bioperl_file =~ /(bioperl.+?)\.tar\.gz/i;
my $bioperl_dir = $1;

opendir BIO, "$DEST_DIR/tmp/$bioperl_dir/Bio/";
move("$DEST_DIR/tmp/$bioperl_dir/Bio/$_", "$DEST_DIR/$_") for readdir BIO;
closedir BIO;

rmtree("$DEST_DIR/tmp") or die "ERROR: Failed to remove directory $DEST_DIR/tmp\n";


# TEST
######

print "\nTesting VEP script\n";

my $test_vep = `perl variant_effect_predictor.pl --help 2>&1`;

$test_vep =~ /ENSEMBL VARIANT EFFECT PREDICTOR/ or die "ERROR: Testing VEP script failed with the following error\n$test_vep\n";

print " - OK!\n";



# CACHE FILES
#############

CACHE:

print "\nThe VEP can either connect to remote or local databases, or use local cache files.\n";
print "Cache files will be stored in $CACHE_DIR\n";
print "Do you want to install any cache files (y/n)? ";

my $ok = <>;

if($ok !~ /^y/i) {
  print "Exiting\n";
  exit(0);
}

# check cache dir exists
if(!(-e $CACHE_DIR)) {
  print "Cache directory $CACHE_DIR does not exists - do you want to create it (y/n)? ";

  my $ok = <>;

  if($ok !~ /^y/i) {
    print "Exiting\n";
    exit(0);
  }
  
  mkdir($CACHE_DIR) or die "ERROR: Could not create directory $CACHE_DIR\n";
}

mkdir($CACHE_DIR.'/tmp') unless -e $CACHE_DIR.'/tmp';

# get list of species
print "\nDownloading list of available cache files\n";

my $num = 1;
my $species_list;
my @files;

if($CACHE_URL =~ /^ftp/i) {
  $CACHE_URL =~ m/(ftp:\/\/)?(.+?)\/(.+)/;
  my $ftp = Net::FTP->new($2) or die "ERROR: Could not connect to FTP host $2\n$@\n";
  $ftp->login($FTP_USER) or die "ERROR: Could not login as $FTP_USER\n$@\n";
  
  foreach my $sub(split /\//, $3) {
    $ftp->cwd($sub) or die "ERROR: Could not change directory to $sub\n$@\n";
  }
  
  push @files, grep {$_ =~ /tar.gz/} $ftp->ls;
}
else {
  push @files, map {$_->[0]} grep {$_->[0] =~ /tar.gz/} @{parse_dir(get($CACHE_URL))};
}

# if we don't have a species list, we'll have to guess
if(!scalar(@files)) {
  print "Could not get current species list - using predefined list instead\n";
  
  @files = (
    "bos_taurus_vep_$API_VERSION.tar.gz",
    "danio_rerio_vep_$API_VERSION.tar.gz",
    "homo_sapiens_vep_$API_VERSION.tar.gz",
    "homo_sapiens_vep_$API_VERSION\_sift_polyphen.tar.gz",
    "mus_musculus_vep_$API_VERSION.tar.gz",
    "rattus_norvegicus_vep_$API_VERSION.tar.gz",
  );
}

foreach my $file(@files) {
  $species_list .= $num++." : ".$file."\n";
}

print "The following species/files are available; which do you want (can specify multiple separated by spaces): \n$species_list\n? ";

foreach my $file(split /\s+/, <>) {
  my $file_path = $files[$file - 1];
  
  my ($species, $file_name);
  
  if($file_path =~ /\//) {
    ($species, $file_name) = (split /\//, $file_path);
  }
  else {
    $file_name = $file_path;
    $file_name =~ m/^(\w+?)\_vep/;
    $species = $1;
  }
  
  # check if user already has this species and version
  if(-e "$CACHE_DIR/$species/$API_VERSION") {
    print "\nWARNING: It looks like you already have the cache for $species (v$API_VERSION) installed.\nIf you continue the existing cache will be overwritten.\nAre you sure you want to continue (y/n)? ";
    
    my $ok = <>;
    
    if($ok !~ /^y/i) {
      print " - skipping $species\n";
      next;
    }
    
    rmtree("$CACHE_DIR/$species/$API_VERSION") or die "ERROR: Could not delete directory $CACHE_DIR/$species/$API_VERSION\n";
  }
  
  my $target_file = "$CACHE_DIR/tmp/$file_name";
  
  print " - downloading $CACHE_URL/$file_path\n";
  
  unless(getstore("$CACHE_URL/$file_path", $target_file) == 200) {
    die "ERROR: Failed to fetch cache file $file_name from $CACHE_URL/$file_path - perhaps you have a proxy/firewall? Set the http_proxy ENV variable if you do\nError code: $response\n";
  }  
  
  print " - unpacking $file_name\n";
  
  
  unpack_tar($target_file, $CACHE_DIR.'/tmp/');
  
  # does species dir exist?
  if(!-e "$CACHE_DIR/$species") {
    mkdir("$CACHE_DIR/$species") or die "ERROR: Could not create directory $CACHE_DIR/$species\n";
  }
  
  # move files
  opendir CACHEDIR, "$CACHE_DIR/tmp/$species/";
  move("$CACHE_DIR/tmp/$species/$_", "$CACHE_DIR/$species/$_") for readdir CACHEDIR;
  closedir CACHEDIR;
}

# cleanup
rmtree("$CACHE_DIR/tmp") or die "ERROR: Could not delete directory $CACHE_DIR/tmp\n";

print "\nSuccess\n";


# SUBS
######

# unpack a tarball
sub unpack_tar {
  my ($arch_file, $dir) = @_;
  
  my $prev_dir = getcwd;
  chdir($dir);
  
  $arch_file =~ s/.+\///g;
  my $tar = Archive::Tar->new;
  $tar->read($arch_file);
  $tar->extract;
  
  chdir($prev_dir);
}

# update or initiate progress bar
sub progress {
  my ($i, $total) = @_;
  
  my $width = 60;
  my $percent = int(($i/$total) * 100);
  my $numblobs = (($i/$total) * $width) - 2;

  return unless $numblobs != $prev_progress;
  $prev_progress = $numblobs;
  
  printf("\r% -${width}s% 1s% 10s", '['.('=' x $numblobs).($numblobs == $width - 2 ? '=' : '>'), ']', "[ " . $percent . "% ]");
}

sub usage {
    my $usage =<<END;
#---------------#
# VEP INSTALLER #
#---------------#

version $VERSION

By Will McLaren (wm2\@ebi.ac.uk)

http://www.ensembl.org/info/docs/variation/vep/vep_script.html#installer

Usage:
perl INSTALL.pl [arguments]

Options
=======

-h | --help        Display this message and quit
-d | --DESTDIR     Set destination directory for API install (default = './')
-v | --VERSION     Set API version to install (default = $VERSION)
-c | --CACHEDIR    Set destination directory for cache files (default = '$HOME/.vep/')
END

    print $usage;
}
