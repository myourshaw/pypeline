#!/usr/bin/perl -w
use strict;
use warnings;
my %individuals;
my %results;
# "C:\Users\Michael\Documents\Lab\Autism\MYO1D sequencing\solexa\haplotype_in.txt"
my @in = split(/_/,$ARGV[0]);
my $outfile = $in[0].".data";
open OUT,">",$outfile;
while(<>)
{
	if(/(\S+)\s+(\S+)\s+(\S+)\s+(\S?)/){
		my $snp = $1;
		my $individual = $2;
		my $value = $3;
		my $hetero = $4;
		$individuals{$individual}++;
		push @{$results{$snp}{$individual}}, ($value,$hetero);
	}
}
print OUT "pedigree\tindividual\tfather\tmother\tsex\taffectation";
print OUT map {"\t$_"} sort {$a cmp $b} keys %results;
print OUT "\n";
foreach my $individual(sort keys %individuals){
	if($individual ne "refNCBI"){
		my $pedigree = substr($individual,0,6);
		my $line = "$pedigree\t$individual\t0\t0\t1\t2";
		foreach my $snp(sort keys %results){
			my $ref = $results{$snp}{"refNCBI"}[0];
			my $allele1 = $results{$snp}{$individual}[0]?$results{$snp}{$individual}[0]:$ref;
			my $allele2 = $results{$snp}{$individual}[1]?$results{$snp}{$individual}[1]eq"*"?$ref:$results{$snp}{$individual}[1]:$ref;
			#$allele1=~tr/[A,C,G,T]/[1,2,3,4]/;
			#$allele2=~tr/[A,C,G,T]/[1,2,3,4]/;
			#$allele1=substr($allele1,0,1);
			#$allele2=substr($allele2,0,1);
			#if($allele1!~/[1234]/){
				#$allele1=0
			#};
			#if($allele2!~/[1234]/){
				#$allele2=0
			#};
			$line .= "\t$allele1 $allele2";
		}
		print OUT "$line\n";
	}
}
