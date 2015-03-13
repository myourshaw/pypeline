#!/usr/bin/perl -w
use strict;
use warnings;
#EnsEMBL Perl API
use Bio::EnsEMBL::Registry;
my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db(
	-host => 'ensembldb.ensembl.org',
	-user => 'anonymous');
my $slice_adaptor = $registry->get_adaptor('human', 'core', 'slice');
for(my $i=1;$i<26;$i++){
    my $slice = $slice_adaptor->fetch_by_region('chromosome', 22);
    my @genes = @{$slice->get_all_Genes};
    foreach my $gene (@genes){
        printf "%s\t%s\n",
        $gene->external_name,
        $gene->stable_id,
        $gene->start,
        $gene->end;
    }
}
