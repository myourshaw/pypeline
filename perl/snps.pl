#!/usr/bin/perl -w
use strict;
use warnings;
#EnsEMBL Perl API
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::Utils::Sequence qw (sequence_with_ambiguity);
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db(
	-host => 'ensembldb.ensembl.org',
	-user => 'anonymous');

my $exon_adaptor = $registry->get_adaptor('human', 'core', 'exon');
my $gene_adaptor = $registry->get_adaptor('human', 'core', 'gene');
my $repeat_adaptor = $registry->get_adaptor('human', 'core', 'repeatfeature');
my $slice_adaptor = $registry->get_adaptor('human', 'core', 'slice');
my $transcript_adaptor = $registry->get_adaptor('human', 'core', 'transcript');
my $variation_feature_adaptor = $registry->get_adaptor('human', 'variation','variationfeature');

my $chr;
my $start;
my $end;
my $ref;
my $var;
my %mismatches;
while(<>){
    chomp;
    if (/^\s*(?:chr)?[\s]*([0-9XYMT]{1,2})[\s:]+([,\d]+)\D+([,\d]+)\s+(\S+)\s+(\S+)/i){
        if(defined($1) && defined($2) && defined($3) && defined($4) && defined($5)){
            $chr = uc($chr);
            if ($chr eq "M"){
                $chr = "MT";
            }
            $chr = $1;
            $start = $2;
            $end = $3;
            $ref = $4;
            $var = $5;
            $start =~ s/,//g;
            $end =~ s/,//g;
            $start = $start <= $end ? $start : $end;
            $end = $end >= $start ? $end : $start;
            $mismatches{$_} = {chr=>$chr,start=>$start,end=>$end,ref=>$ref,var=>$var};
        }
    }
}
print "done\n";
foreach my $mismatch (sort keys %mismatches){
    my $slice = $slice_adaptor->fetch_by_region('chromosome', $mismatches{$mismatch}->{chr}, $mismatches{$mismatch}->{start}, $mismatches{$mismatch}->{end});

}

sub chr2numchr{
    my $chr = uc(shift);
    return $chr eq "X" ? 23 : $chr eq "Y" ? 24 : ($chr eq "M" || $chr eq "MT") ? 25 : $chr;
}

sub numchr2chr{
    my $chr = shift;
    return $chr == 23 ? "X" : $chr == 24 ? "Y" : $chr == 25 ? "M" : $chr;
}
