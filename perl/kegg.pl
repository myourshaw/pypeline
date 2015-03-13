#!/usr/bin/perl -w
use strict;

#http://www.genome.jp/kegg/docs/keggapi_manual.html
#http://www.genome.jp/kegg/soap/doc/keggapi_manual.html

use SOAP::Lite;

my $wsdl = 'http://soap.genome.jp/KEGG.wsdl';

my $serv = SOAP::Lite->service($wsdl);

my $offset = 1;
my $limit = 5;

my $top5 = $serv->get_best_neighbors_by_gene('eco:b0002', $offset, $limit);

foreach my $hit (@{$top5}) {
  print "$hit->{genes_id1}\t$hit->{genes_id2}\t$hit->{sw_score}\n";
}
