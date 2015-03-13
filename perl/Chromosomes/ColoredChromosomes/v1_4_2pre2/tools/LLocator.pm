#
#	LLocator.pm
#Tue Jan 15 12:13:25 CET 2002
#

package LLocator;
require 5.000;
#@ISA = qw(SuperClass);

#
#	<§> dependencies
#

#use SuperClass;

#
#	<§> data structures
#

#
#	<§> standard implementation
#

sub new { my ($class, $path, $aliases) = @_;
	#my $self = bless($class->SUPER::new(@args), $class);
	# if no superclass is in place use the following instead
	my $self = bless({}, $class);
	# own initialization follows

	$self->{path} = $path;
	$self->{aliases} = $aliases;

	return $self;
}

#
#	<§> method implementation
#

sub path { my ($self) = @_; return $self->{path}; }
sub aliases { my ($self) = @_; return $self->{aliases}; }

# <!> names are converted through the alias table, but are ultimately stored as uppercase strings to canonicalize names; mb: <i> make optional

sub grepLociFromChromosome { my ($self, $loci, $chromosomeName) = @_;
	die "The method grepLociFromChromosome is called for with an abstract class.\n";
	# this method should retrieve information for a chromosome and return
	# the loci from $loci

	# should return a dictionary:
	# { c => $lociList, l => $foundLoci }
	# when $lociList is an array containing all loci from the respective source
	# on the chromosome and $foundLoci contains the loci from $loci found on this
	# chromosome
	# The format of each array entry is a dictionary as follows:
	#
	# name : the locus name; ${chromosomeName}ptr/qtr are needed for the gaps
	# phyPos : a relative positino ranging from ptr to qtr (0 to 1)
	# band : the name of the band on the chromosome

}


sub calculateGaps { my ($self, $loci) = @_;
	my $allGaps = [];

	foreach $chromosomeName (1..22, 'X', 'Y') # chrom 'Y' is not supported by ldb
	{
		my $g = $self->grepLociFromChromosome($chromosomeName, $loci);
		my $l = [sort { $a->{phyPos} <=> $b->{phyPos} } @{$g->{l}}];
		unshift(@{$l}, grep { $_->{name} eq "${chromosomeName}ptr" } @{$g->{c}});
		push(@{$l}, grep { $_->{name} eq "${chromosomeName}qtr" } @{$g->{c}});
		my $i;

		for ($i = 1; $i < @{$l}; $i++)
		{
#			printf("Gap: %.4f %s\n", $l->[$i]->{phyPos} - $l->[$i - 1]->{phyPos},
#				"[$l->[$i - 1]->{name} - $l->[$i]->{name}]");
			push(@{$allGaps}, { gap => $l->[$i]->{phyPos} - $l->[$i - 1]->{phyPos},
				start => $l->[$i - 1]->{name},
				stop => $l->[$i]->{name}
			});
		}
	}
	return $allGaps;
}

sub searchLoci { my ($self, $loci) = @_;
	my $positionList = [];

	$loci = dictWithKeys($loci);

	foreach $chromosomeName (1..22, 'X', 'Y') # chrom 'Y' is not supported by ldb
	{
		my $g = $self->grepLociFromChromosome($chromosomeName, $loci);

		push(@{$positionList}, grep {
			$flag = defined($loci->{$_->{name}});
			($flag? $loci->{$_->{name}} = 1: 0);
			$flag
		} @{$g->{c}});
	}
	return $positionList;
}

1;
