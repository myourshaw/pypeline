#
#	CCcolors.pm
#Creation Date
#

package CCcolors;
require 5.000;
#@ISA = qw(SuperClass);

#
#	<§> dependencies
#

use Set;
#use SuperClass;

#
#	<§> data structures
#


$standardColors = {
	black => { r => 0, g => 0, b => 0 },
	white => { r => 1, g => 1, b=> 1 },
	grey => { r => 0.6666, g => 0.6666, b=> 0.6666 },
	chrome => { r => 0.078, g => 0.027, b => 0.329 },
	acen => { r => 1, g => 0.5, b => 0.5 },
	gneg => { r => 1, g => 1, b=> 1 },
	gpos25 => { r => 0.75, g => 0.75, b => 0.75 },
	gpos50 => { r => 0.50, g => 0.50, b => 0.50 },
	gpos75 => { r => 0.25, g => 0.25, b => 0.25 },
	gpos100 => { r => 0, g => 0, b => 0 },
	gvar => { r => 0, g => 0, b => 0 },
	stalk => { r => 0.50, g => 0.50, b => 0.50 }
};

#
#	<§> standard implementation
#

sub new { my ($class, $colorConfig) = @_;
	#my $self = bless($class->SUPER::new(@args), $class);
	# if no superclass is in place use the following instead
	my $self = bless({}, $class);

	# own initialization follows
	$self->{colors} = { %{$standardColors} };
	$self->{config} = $colorConfig;
	return $self;
}

#
#	<§> method implementation
#

sub addColors { my ($self, $newColors) = @_;
	mergeDict2dict($newColors, $self->{colors});
}
sub rgbColor { my ($self, $def) = @_;
	return $def if (ref($def) eq 'HASH');
	defined($self->{colors}{$def})
	|| die "Color $def is unknown. Known colors are:".
	   join(' ', keys(%{$self->{colors}})). "\n";
	return $self->{colors}{$def};
}

sub colorForValue { my ($self, $v) = @_;
	my ($colorConfig, $color, $i) = ($self->{config});
	# <!> we assume $colorConfig tb sorted
	# search for the interpolation port; <i> binary search
	for ($i = 1; $i < @{$colorConfig->{interpolationList}}; $i++)
	{	last if ($v <= $colorConfig->{interpolationList}[$i]{value});
	}
	return $colorConfig->{specialColor} if ($i >= @{$colorConfig->{interpolationList}});

	my $c1 = $colorConfig->{interpolationList}[$i - 1]{color};
	my $c2 = $colorConfig->{interpolationList}[$i]{color};
	my $cr = {};
	$i = '';

	foreach $component ('r', 'g', 'b')
	{
		$cr->{$component} = (1 - $v) * $c1->{$component} + $v * $c2->{$component};
	}
	return $cr;
#<t> return Postscript::blendedColorFromColorsAndFraction($c1, $c2, $v);
}

sub defaultValue { my ($self) = @_;
	return $self->{config}{defaultValue};
}

1;
