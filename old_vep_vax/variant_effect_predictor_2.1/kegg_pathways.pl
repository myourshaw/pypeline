#!/usr/bin/perl -w
use strict;

use SOAP::Lite;
use DBI;
use DBD::MySQL;
use Getopt::Long;
Getopt::Long::Configure('no_ignore_case');

#-o /Users/myourshaw/lab/pypeline/variant_effect_predictor_2.1/kegg_pathways.txt -u myourshaw -p m.cha3ly

sub SOAP::Serializer::as_ArrayOfstring{
  my ($self, $value, $name, $type, $attr) = @_;
  return [$name, {'xsi:type' => 'array', %$attr}, $value];
}

sub SOAP::Serializer::as_ArrayOfint{
  my ($self, $value, $name, $type, $attr) = @_;
  return [$name, {'xsi:type' => 'array', %$attr}, $value];
}

sub unquote($){
    my $str = shift;
    $str =~ s/'/\\'/g;
    return $str;
}

my %options;
GetOptions(
    'o|output=s' => \$options{output},
    'host' => \$options{host},
    'port' => \$options{port},
    'u|user=s' => \$options{user},
    'p|password=s' => \$options{password},
    'platform=s' => \$options{platform},
    'database=s' => \$options{database},
    'table=s' => \$options{table},
);
$options{host} ||= "myourshaw-dev.genome.ucla.edu";
$options{port} ||= "3306";
$options{platform} ||= "mysql";
$options{database} ||= "vw";
$options{table} ||= "kegg_gene_pathway";

my $create_tables = 0;
if (defined($options{host}) && defined($options{user}) && defined($options{password})){
    $create_tables = 1;
}

my $dsn = "dbi:$options{platform}:$options{database}:$options{host}:$options{port}";
my $insert = "INSERT INTO $options{table} (gene_id, path_id, pathway) VALUES ('%s', '%s', '%s')";
my $truncate = "TRUNCATE TABLE $options{table}";
my $conn;
if ($create_tables){
    $conn = DBI->connect($dsn, $options{user}, $options{password})
        or die "Unable to connect: $DBI::errstr\n";
    my $qh = $conn->prepare($truncate);
    $qh->execute() or die "Unable to execute $truncate: $DBI::errstr\n";
}

my $write_files = 0;
if (defined($options{output})){
    $write_files = 1;
}
if ($write_files){
    open OUT, ">", $options{output}
        or die "can't open $options{output}\n";
    print OUT "#gene_id\tpath_id\tpathway\n";
}

my $wsdl = 'http://soap.genome.jp/KEGG.wsdl';
my $paths = SOAP::Lite
             -> service($wsdl)
             -> list_pathways("hsa");
foreach my $path (@{$paths}) {
    my $genes = SOAP::Lite
             -> service($wsdl)
             -> get_genes_by_pathway($path->{entry_id});
    foreach my $gene_id(@{$genes}){
        my $path_id = $path->{entry_id};
        (my $pathway = $path->{definition}) =~ s/ - Homo sapiens \(human\)$//;
        if ($create_tables){
            my $query = sprintf($insert, ($gene_id, $path_id, unquote($pathway)));
            my $qh = $conn->prepare($query);
            $qh->execute() or die "Unable to execute $query: $DBI::errstr\n";
        }
        if ($write_files){
            print OUT "$gene_id\t$path_id\t$pathway\n";
        }
    }
}

