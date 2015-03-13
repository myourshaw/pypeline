#
#	CCAitags.pm
#Wed Nov 14 09:51:40 CET 2001
#

package CCAitags;
require 5.000;
@ISA = qw(CCAnnotation);

#
#	<§> dependencies
#

use CCAnnotation;
use CCcolors;
use Set;

#
#	<§> data structures
#

#
#	<§> standard implementation
#

sub newWithStreamChromosomeConfig { my ($class, @args) = @_;
	my $self = bless($class->SUPER::newWithStreamChromosomeConfig(@args), $class);
	# if no superclass is in place use the following instead
	# my $self = bless({}, $class);

	# own initialization follows
	$self->{mylabels} = $self->{labels}{byChromosome}{$self->{position}{chromosomeName}};
	$self->{colorManager} = CCcolors->new(firstDef($self->{annotation}{colors},
		$self->{labels}{colors}));

	return $self;
}

#
#	<§> method implementation
#

sub draw { my ($self) = @_;
	my ($x, $y, $w, $size) = ($self->{c}{x}, $self->{c}{yBottom},
		$self->{c}{w}, $self->{c}{size});
	my $hw = $w / 2.0;
	my $drawingSize = $self->{gmap}{maxh} * $self->{gmap}{scalingFactor}
		* $self->{annotation}{defaultColorHeight};

	$self->{ps}->gsave();

	# <i> we sort labels here to draw small values last (idea: p-values)
	# extension: ascending order vs. descending order vs. compare function

	# sort by applying the function $sf and draw in that order
	my $sfDef = $self->{annotation}{sortBy};
	my $sf = defined($sfDef)
		? eval 'sub { my $v = $_[0]; '.$sfDef.' }'
		: sub { - $_[0] };
	foreach $label (sort { $sf->($a->{v}) <=> $sf->($b->{v}) } @{$self->{mylabels}})
	{	my $labelValue = $label->{v};
		# is it a color band or something else?
		if (!defined($labelValue)) {
			$labelValue = $self->{colorManager}->defaultValue();
			next if (!defined($labelValue));
		}
		my $psColor = $self->{colorManager}->colorForValue($labelValue);
		my $pos = $self->absolutePositionFromSpec($label->{p});
		my $yp = $y + $pos - $drawingSize / 2.0;

#printf STDERR "ColorRect:%.2f %.2f %.2f %.2f [r%.2fg%.2fb%.2f]\n", $x - $hw, $yp, $w, $drawingSize, $psColor->{r}, $psColor->{g}, $psColor->{b};
		$self->{ps}->colorRect($psColor, $x - $hw, $yp, $w, $drawingSize);
	}
	$self->{ps}->grestore();
}

1;

#	$placementConfig->{labels} = defined($labelPlacement)?
#		$config->{labelPlacements}{$labelPlacement} :
#		(!defined($placementConfig->{labels})?
#			$config->{labelPlacements}{standard}:
#			$config->{labelPlacements}{$placementConfig->{labels}});
