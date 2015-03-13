#
#	CCAnnotation.pm
#Wed Nov  7 10:27:11 CET 2001
#

package CCAnnotation;
require 5.000;
#@ISA = qw(SuperClass);

#
#	<§> dependencies
#

#use SuperClass;
use Set;

#
#	<§> data structures
#

# <!> defined also in coloredChromosomes.pl -> unify in another structure
$chromNameRE = '(?:\\d*)[A-Z]?';

#
#	<§> standard implementation
#

#	the following keys are usually defined:
#	global: the whole config-dict as read from the main config file
#	placement: the sub-dict in the global config which defines the current chrom. layout
#	map: positions and sizes where chromosomes have been placed
#	labels: labels/makers to be drawn inside and outside the chromosome
#		these are read globally, since this may be a lot of data (and should
#		not be read more than once)
#	Colors: a global color managing object

sub newWithStreamChromosomeConfig { my ($class, $ps, $c, $gc) = @_;
	#my $self = bless($class->SUPER::new(@args), $class);
	# if no superclass is in place use the following instead
	my $self = bless({}, $class);
	# own initialization follows

	$self->{ps} = $ps;
	$self->{c} = $c;
	foreach $k (keys %{$gc}) {
		$self->{$k} = $gc->{$k};
	}
	return $self;
}

#
#	<§> method implementation
#

sub positionFromSpec { my ($self, $spec) = @_;
	# several notations are supported

	# accept positions as 'FcX' here
	# where X is the chromosome name, F is a float and [pq] indicates the arm
	my ($p) = ($spec =~ m{^([01]|\d*\.\d+)c($chromNameRE)}ogs);
	if (defined($p)) {
		return $p;
	}

	my ($pos1, $arm1, $arm2, $band, $pos2) = ($spec =~
		m{^(?:(?:(\d*\.\d+)([pq]))
		  |(?:[\dA-Z]*([pq:])?((?:[\dA-Z]+)(?:\.\d+)?)(?:\:(\d*(?:\.\d+)?))?)$)
		}ox);
	return ($arm1, $pos1) if (defined($arm1));
	$band = $self->{c}{chromosome}{bandingDict}{$arm2.$band};
	if (!defined($band))
	{	$spec =~ m{(.*):.*}og;
		print STDERR "Unknown band:$1. Known bands are: ",
			join(' ', sort(keys(%{$self->{c}{chromosome}{bandingDict}}))), "\n";
		return undef;
	}
	return ($arm2, $band->{start} + $pos2 * ($band->{stop} - $band->{start}));
}

sub absolutePositionFromSpec { my ($self, $spec) = @_;
	my ($p, $q) = ($self->{position}{p}, $self->{position}{q});
	my $size = $self->{position}{size};
	# this spec is of the form POScCHR where
	# POS is a fraction and CHR is the chrom name
	if ($spec =~ m{^(\d*\.\d+)c(\d+|X|Y)}ogs)
	{
		return $1 * $size;
	}
	my ($arm, $pos) = $self->positionFromSpec($spec, $chromosome);
	# <N> either no p/q separtion of chromosomes exists
	# or bands are specified uniformly (banding vs pBanding/qBanding)
	return ($arm eq ':' || $size > 0)
	? (1 - $pos) * $size
	: ($arm eq 'q'? $q - $pos * $q: $q + $pos *$p);
}

# <i><!> find a clean solution to merge functionality with absolutePositionFromSpec
sub absoluteSizeOfBand { my ($self, $bandName) = @_;
	my ($p, $q) = ($self->{position}{p}, $self->{position}{q});
	my $size = $self->{position}{size};
	my $band = $self->{c}{chromosome}{bandingDict}{$bandName};
	# size take precendce over $p/$q
	my $referenceSize = $size > 0? $size
	: (substr($bandName, 0, 1) eq 'q'? $q: $p);
#	my $referenceSize = ($p > 0 || $q > 0)
#	? (substr($bandName, 0, 1) eq 'q'? $q: $p): $size;
	return $referenceSize * ($band->{stop} - $band->{start});
}

sub midBandPosition { my ($self, $spec, $position) = @_;
	my $name = $self->{c}{chromosome}{name};
	return $self->absolutePositionFromSpec($name.$spec
		.':'.firstDef($position, '0.5'));
}

sub draw { my ($self) = @_;
	die "The method 'draw' is called for the abstract class ".
		"CCAnnotation which should not happen.";
}

1;

=head1 NAME

CCAnnotation.pm - A class from which all annotations used with coloredChromosomes.pl are to be derived

=head1 SYNOPSIS

 Subclass by:
 @ISA = qw(CCAnnotation);
 use CCAnnotation;

=head1 DESCRIPTION

This class is an abstract superclass from which all classes implementing annotations for I<coloredChromosomes.pl> are derived. Each subclass has to implement the method I<draw> which is to draw a single chromosome. Each instance of CCAnnotation (or a subclass) can assume that an instance variable called 'ps' is initialized. This is a Postscript object which supports many postscript operators (see B<Postscript>). When I<draw> is called the instance may safely assume the shape of the current chromosome to be the chromosomal shape, which is clipped to if the annotation is internal.

The following methods can be used to determine positions along the chromosome:

=over 4

=item *

positionFromSpec($spec) - return the a relative positions (0 to 1) of the position I<$spec> along the chromosome

=item *

absolutePositionFromSpec - return the vertical position in points of the specified position

=item *

absoluteSizeOfBand($bandName) - return the size of a band from the current chromosome in points

=item *

midBandPosition($bandName) - return the vertical middle position of the band from the current chromosome named I<$bandName> in points.

=back

The format of specifications for chromosomal locations is detailed in the man page of B<coloredChromosomes.pl>. The implemented annotation classes can be referred on as examples on how to explore the use of these methods.

=head1 SEE ALSO

coloredChromosomes.pl(1) Postscript(3)
