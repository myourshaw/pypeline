#
#	Statistics::Histogram.pm
#Mon Mar 12 13:50:00 MET 2001
#

package Statistics::Histogram;
require 5.000;

#
#	<§> dependencies
#

#
#	<§> methods
#

sub newWithData { my ($class, $data) = @_;
	my $self = {};
	# own initialization follows
	bless($self, $class);
	$self->{data} = $data;
	return $self;
}

sub histgromWithSlitLength { my ($self, $slit, $slitStart) = @_;
	my $slits = {};
	($self->{slitStart}, $self->{slit}) = ($slitStart, $slit);
	foreach $el (@{$self->{data}})
	{
		$slits->{ int(($el - $slitStart) / $slit) }++;
	}
	return $self->{histogram} = $slits;
}

sub printHistogram { my ($self) = @_;
	foreach $slot (sort { $a <=> $b } keys %{$self->{histogram}})
	{
		printf("%.4f\t%d\n",
			($slot * $self->{slit}) + $self->{slitStart},
			$self->{histogram}->{$slot});
	}
}

1;
