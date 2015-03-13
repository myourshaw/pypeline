#
#	LLldb.pm
#Tue Jan 15 12:13:18 CET 2002
#

package LLldb;
require 5.000;
@ISA = qw(LLocator);

$ldbLocation = '/ldb/ldb';

#
#	<§> dependencies
#

use TempFileNames;
use LLocator;
use Set;

#
#	<§> data structures
#

#
#	<§> standard implementation
#

sub new { my ($class, $path, $aliases) = @_;
	my $self = bless($class->SUPER::new($path, $aliases), $class);
	# if no superclass is in place use the following instead
	#my $self = bless({}, $class);
	# own initialization follows
	return $self;
}

#
#	<§> method implementation
#

sub digestChromosome { my ($self, $chromosome, $chromosomeName) = @_;
	my $markerList = [];

	$chromosome =~ s{^([^\s]+)\s+(\d+\.\d+).{17}([^\s]+)}
		{ push(@{$markerList}, {
			name => grep(/\Q$1/, ('ptr', 'qtr'))? "$chromosomeName$1": $1,
			phyPos => $2, band => "$chromosomeName$3"
		}) }omge;
	my $length = $markerList->[$#$markerList]->{phyPos};

	foreach $marker (@{$markerList})
	{
		my $posName = sprintf("%.6f", $marker->{pos} = $marker->{phyPos} / $length);
		$marker->{pos} = $posName, $marker->{p} = "${posName}c$chromosomeName";
		$marker->{c} = $chromosomeName;
	}
	return $markerList;
}

# <!> names are converted through the alias table, but are ultimately stored as uppercase strings to canonicalize names; mb: <i> make optional

sub grepLociFromChromosome { my ($self, $chromosomeName, $loci) = @_;
	Log("Scrutinizing chromosome $chromosomeName", 1);
	my $fullPath = firstDef($self->path(), '/mnt/cdrom/ldb')."/chrom$chromosomeName/gmap";
	if (! -e $fullPath) {
		Log("Specification file: $fullPath does not exist.");
		return undef;
	}
	my $aliases = $self->aliases();
	my $chromosome = readFile($fullPath);
	my $lociList = $self->digestChromosome($chromosome, $chromosomeName);
	my $foundLoci = [grep {
		my $tname = defined($aliases->{$_->{name}})? $aliases->{$_->{name}}: $_->{name};
		my $flag = defined($loci->{$tname});
		($flag? ($loci->{$tname} = 1, $_->{name} = uc($tname)) : 0), $flag
	} @{$lociList}];

	return { c => $lociList, l => $foundLoci };
}


1;
