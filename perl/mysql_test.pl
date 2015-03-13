#! /usr/bin/perl

use strict;
use warnings;

use DBI;

my $dbh=DBI->connect(
  "DBI:mysql:database=vax_cluster_74;host=cortex.local",
  "vax",
  "vax"
) || die "Error connecting to database: $!\n";

#my $sth = $dbh->prepare();

my $sth = $dbh->prepare("CALL `vax_cluster_74`.`get_nhlbi_evs`('9', 37783990, 'T', 'G')");

$sth->execute();

if (my $hash_ref = $sth->fetchrow_hashref()) {
    my $foo = $hash_ref->{EA_AA_GT_freq};
    print $foo;
}
$sth->finish;

while (my $hash_ref = $sth->fetchrow_hashref()) {
    my $delim = "";
    foreach my $key (keys (%{$hash_ref}))
    {
        $hash_ref->{$key} = "" if !defined ($hash_ref->{$key}); # NULL value?
        print $delim, $hash_ref->{$key};
        $delim = ",";
    }
    print "\n";
}


while (my $row = $sth->fetchrow_arrayref) {
  print $row->[1] . "\n";
}

exit;
