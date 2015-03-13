#!/usr/bin/perl -w
use strict;


use MRS::Client;

# Create a new MRS Client object and point to the right server
my $client = MRS::Client->new ( search_url => 'http://localhost:18081/', blast_url => '
http://localhost:18082/', clustal_url => 'http://localhost:18083/');

#Receive an approximate count of the results by calling the count method
print $client->db ('embl')->find ('(kw:transporter OR de:transporter) AND oc:"Viridiplantae" AND de:complete')->count;

#form the query and store the result.
my $query = $client->db ('embl')->find ('(kw:transporter OR de:transporter) AND oc:"Viridiplantae" AND de:complete');

#Now obtain the results by examining the data structure and call the next method
while (my $record = $query->next) {
print $record . "\n";
}