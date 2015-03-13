#
#	CCAsideline.pm
#Thu Jan 17 11:48:18 CET 2002
#

package CCAsideline;
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
	$self->{colorManager} = CCcolors->new($self->{labels}{colors});

	return $self;
}

#
#	<§> method implementation
#

sub tidentity { my ($self, $v, $h) = @_;
	$v = $v > $h->{max}? $h->{max}: ($v < $h->{min}? $h->{min}: $v);
	return ($v - $h->{min}) / ($h->{max} - $h->{min});
}
sub tlogarithmic { my ($self, $v, $h) = @_;
	return 1 if ($v <= 0);	# we return the normalized max
	# log to the base of 10
	my $t = - log($v) /  2.30258509299405;
	$t = $t > $h->{max}? $h->{max}: ($t < $h->{min}? $h->{min}: $t);
	return ($t - $h->{min}) / ($h->{max} - $h->{min});
}
sub toneMinus { my ($self, $v, $h) = @_;
	my $t = 1 - $v;
	$t = $t > $h->{max}? $h->{max}: ($t < $h->{min}? $h->{min}: $t);
	return ($t - $h->{min}) / ($h->{max} - $h->{min});
}


# we do a - log($v) and normalize for {annotation}{logTransform}{min}
# t: short for transform
sub t { my ($self, $v) = @_;
	if (defined($self->{annotation}{logTransform})) {
		return 1 if ($v <= 0);	# we return the normalized max
		# log to the base of 10
		my $t = log($v) /  2.30258509299405 / $self->{annotation}{logTransform}{min};
		return $t > 1? 1: $t;
	}
	return 1 - $v;
}

sub num {((@_ ? $_[0] : $_) =~ /^(\d+\.?\d*)/ => 0)[0]}

sub draw { my ($self) = @_;
	my ($x, $y, $w, $size) = ($self->{c}{x}, $self->{c}{yBottom},
		$self->{c}{w}, $self->{c}{size});
#	my ($x, $y, $w) = ($self->{position}{x}, $self->{position}{y}, $self->{position}{w});
	my $hw = $w / 2.0;
	my $kernelWidth = firstDef($self->{annotation}{valueBin}, 0.005);
	# the calculation used to determine $kernelStep is a graphical compromise
	# absolute steps should be equal on each page
	# since chromosomes might be strechted differently on different pages,
	# the absolute (graphical) step size is taken from the page harbouring the
	# (physically) largest chromosome which may result in finer physical mapping
	# on other pages. $kernelStep is a fraction of the chromosome (a relative unit).
	my $kernelStep = $kernelWidth * ($self->{gmap}{maxh}
		/ $self->{c}{chromosome}{size});
	my $valueHeight = $size * $kernelStep;
	my $xs = $x + $hw + $self->{annotation}{inset};	# where to start drawing x-wise
	my $dw = $self->{annotation}{width};		# width of area to draw in
	my $t = 't'.firstDef($self->{annotation}{transformation}, 'oneMinus');
	my $th = firstDef($self->{annotation}{transformationHints}, { min => 0, max => 1});

	$self->{ps}->gsave();

	if (defined($self->{annotation}{drawCutoffLineAt})) {
		my $cutoffLineWidth = firstDef($self->{annotation}{cutoffLineWidth}, .25);
		my @cos = ref($self->{annotation}{drawCutoffLineAt})
		 ? @{$self->{annotation}{drawCutoffLineAt}}
		 : ($self->{annotation}{drawCutoffLineAt});

		$self->{ps}->setlinewidth($cutoffLineWidth);
		$self->{ps}->setrgbcolor($self->{colorManager}->rgbColor(
			$self->{annotation}{cutoffLineColor}));
		foreach $coValue (@cos) {
			my $cox = $xs + &{$t}($self, $coValue, $th) * $dw;
			$self->{ps}->line($cox, $y, $cox, $y + $size);
			$self->{ps}->stroke();
		}
	}

	my $sidelinecolor = $self->{annotation}{color};
	my $color = $self->{colorManager}->rgbColor(firstDef($sidelinecolor, 'black'));
	my $sidelineWidth = firstDef($self->{annotation}{sidelineWidth}, .7);
	$self->{ps}->setrgbcolor($color);
	$self->{ps}->setlinewidth($sidelineWidth);

	# we sort by positions
	my @labelsS = defined($self->{mylabels})
		? sort { num($a->{p}) <=> num($b->{p}) } (@{$self->{mylabels}})
		: ();
	my $argmaxDef = $self->{annotation}{argmaxInBin};
	my $argmaxFct = defined($argmaxDef)
		? eval 'sub { my $v = $_[0]; '.$argmaxDef.' }'
		: sub { - $_[0] };
	my ($i, $j, $min, $preMin);
	for ($i = 0; $i * $kernelStep < 1; $i++, $preMin = $min) {
		# we restart over each time, since we splice away the beginning afterwards
		for ($j = 0; $j < int(@labelsS)
		     && num($labelsS[$j]{p}) <= ($i + 1) * $kernelStep; $j++)
			{ }
		my ($maxValue, $argmax);
		foreach $l (splice(@labelsS, 0, $j)) {
			if (!defined($maxValue) || $argmaxFct->($l->{v}) > $maxValue) {
				$maxValue = $argmaxFct->($l->{v});
				$argmax = $l->{v};
			}
		}
		$argmax = $self->{annotation}{defaultValue} if (!defined($argmax));
		my $xn = $xs + &{$t}($self, $argmax, $th) * $dw;
		if ($i) {
			$self->{ps}->lineto($xn, $y + $i * $valueHeight);
		} else {
			$self->{ps}->moveto($xn, $y + $i * $valueHeight);
		}
		$self->{ps}->lineto($xn, $y + ($i + 1)* $valueHeight)
	}
	$self->{ps}->stroke();

	$self->{ps}->grestore();
}

1;
