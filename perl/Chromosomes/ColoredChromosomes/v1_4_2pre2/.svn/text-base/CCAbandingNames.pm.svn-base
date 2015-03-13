#
#	CCAbandingNames.pm
#Mon Nov 12 09:11:43 CET 2001
#

package CCAbandingNames;
require 5.000;
@ISA = qw(CCAnnotation);

#
#	<§> dependencies
#

use CCAnnotation;

#
#	<§> data structures
#

#
#	<§> standard implementation
#

sub newWithStreamChromosomeConfig { my ($class, @args) = @_;
	my $self = bless($class->SUPER::newWithStreamChromosomeConfig(@args), $class);
	# if no superclass is in place use the following instead
	#my $self = bless({}, $class);
	# own initialization follows

	$self->{bconfig} = defined($self->{annotation}{config})
		? $self->{global}{bandNamePlacements}{$self->{annotation}{config}}
		: $self->{annotation};

	return $self;
}

#
#	<§> method implementation
#

sub	psDrawBandNames { my ($stream, $chromosome, $chromPos, $config) = @_;
}

sub draw { my ($self) = @_;
	my ($x, $y, $w, $size) = ($self->{c}{x}, $self->{c}{yBottom}, $self->{c}{w}, $self->{c}{size});

#	return if (!grep { $self->{c}{tag} eq $_ } @{$self->{bconfig}{tags}});

	$self->{ps}->gsave();
	$self->{ps}->setlinewidth(0);

	my $font = $self->{bconfig}->{font};
	$self->{ps}->setfont($font, { name => 'Helvetica', size => 6 });
	$self->{ps}->setrgbcolor($self->{bconfig}{color});

	my $tx = $x - ( $w/2 + $self->{bconfig}{distance} );
	my ($lastInset, $lineDistance) = (0, $font->{size}/3);

	# set default label inset size
	$self->{bconfig}{labelInset} = $font->{size} * 3
		if (!defined($self->{bconfig}{labelInset}));

	# sort bandingnames according to position to allow proper placement
	my $bd = $self->{c}{chromosome}{bandingDict};

	my @sortedBandNames = sort { $self->midBandPosition($a) <=>
		$self->midBandPosition($b)} keys %{$bd};

	my $i;
	foreach $bandName (@sortedBandNames)
	{
		my $size = $self->absoluteSizeOfBand($bandName);
		my $ty = $y + $self->midBandPosition($bandName);
		my $didRescale = 0;
		my $ttx = $tx;		# do we have to shift?

		# this is a hack to clean up band names
		# <i><N> should be replaced by an architctual change
		$bandName =~ s/^://os;
		if ($size < $font->{size})
		{
			if (uc($self->{bconfig}{rescaleFont}) eq 'YES')
			{
				$didRescale = 1;
				$self->{ps}->gsave(), $self->{ps}->setfont(
					{ name => $font->{name}, size => $size });
			} else {	# inset labels
				$lastInset = $self->{bconfig}{labelInset} - $lastInset;
				$ttx -= $lastInset;
				$self->{ps}->strokeline($ttx + $lineDistance, $ty, $tx, $ty)
					if ($lastInset);
			}
		} else {
			$lastInset = 0;	# reset the inset: normal bands are not insetted
		}
		$self->{ps}->showRyC($bandName, $ttx, $ty);
		$self->{ps}->grestore() if ($didRescale);
	}
	$self->{ps}->grestore();
}

1;
