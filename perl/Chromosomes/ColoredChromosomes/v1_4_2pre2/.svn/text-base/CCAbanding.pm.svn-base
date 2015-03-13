#
#	CCAbanding.pm
#Wed Nov  7 13:52:13 CET 2001
#

package CCAbanding;
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
	return $self;
}

#
#	<§> method implementation
#

sub drawCentromere { my ($self, $w, $h) = @_;
	my $stream = $self->{ps};
	my $hw = $w / 2.0;
	my $l = $w / 10.0;
	$stream->gsave();

	$stream->newpath();								# define centromere region
	$stream->setgray(0.6666);
	$stream->rect(- $hw, 0, $w, $h);
	$stream->gsave();
		$stream->fill();
	$stream->grestore();
	$stream->clip();
	$stream->newpath();								# draw inside telomere region
	$stream->setgray(0.2222);
	$stream->setlinewidth($l);
	my $i;
	for ($i = 0; $i < 10; $i++)
	{	my $x_ = - 1.5 *$w + $i * 2 * $l;
		$stream->line($x_, - $h, $x_ + 4 * $h, 2 * $h);
	}
	$stream->stroke();

	$stream->grestore();
}

sub drawBanding { my ($self, $w, $l, $banding) = @_;
	my $placement = $self->{placement};
	my $defColors = $placement->{bandingColors};
	my $hw = $w / 2.0;
	@{$banding} = sort { $a->{start} <=> $b->{start} } @{$banding};
	my $lastColor = 'black';
	my $chrome = lc($self->{annotation}{drawBandingWithVisualEffect}) eq 'chrome';
	$self->{Colors}->addColors($defColors);

	foreach $band (@{$banding})
	{
		if (defined($band->{color})) {
			if ($band->{color} eq 'stem') {
				next();
				}
			elsif ($band->{color} eq 'centromere') {
				$self->drawCentromere($w, $l * ($band->{stop} - $band->{start}));
				next();
				}
			}

		# oscillate automatically if nothing different is defined
		$color = defined($band->{color})? $band->{color}:
			($lastColor eq 'black'? 'white': 'black');
		my $rgbColor = $self->{Colors}->rgbColor($color);
		if ($chrome && $color eq 'black')
		{
			$self->{ps}->hshiningRect($self->{Colors}->rgbColor('chrome'),
				$self->{Colors}->rgbColor('white'), 0.4,
				- $hw, $l * $band->{start}, $w, $l * ($band->{stop} - $band->{start}), .2);
		} else {
			$self->{ps}->colorRect($rgbColor, - $hw, $l * $band->{start},
				$w, $l * ($band->{stop} - $band->{start}));
		}
		$lastColor = $color;
	}
}

sub draw { my ($self) = @_;
	my ($x, $y, $p, $q, $w, $size) = ($self->{position}{x}, $self->{position}{y},
		$self->{position}{p}, $self->{position}{q}, $self->{position}{w},
		$self->{position}{size});
	my $bdict = $self->{c}{chromosome};	# bandingdict
	my $hw = $w / 2.0;
	$self->{ps}->gsave();

	$self->{ps}->clip();								# clip to chromosome path

	$self->{ps}->gsave();
		if (!$self->{position}{p}) {	# no centromere
			$self->{ps}->translate($x, $self->{position}{yBottom} + $size);
			$self->{ps}->scale(1, -1);
			$self->drawBanding($w, $size, $bdict->{banding});
		} else {
			$self->{ps}->translate($x, $y);
			$self->drawBanding($w, $p, $bdict->{pBanding});
			$self->{ps}->scale(1, -1);
			$self->drawBanding($w, $q, $bdict->{qBanding});
		}
	$self->{ps}->grestore();
	
	$self->{ps}->grestore();
}

1;
