#
#	CCAplain.pm
#Wed Nov 14 14:13:25 CET 2001
#

package CCAplain;
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

sub draw { my ($self) = @_;
	my ($x, $y, $p, $q, $w) = ($self->{position}{x}, $self->{position}{y},
		$self->{position}{p}, $self->{position}{q}, $self->{position}{w});
	my $bdict = $self->{c}{chromosome};	# bandingdict
	my $hw = $w / 2.0;

	$self->{ps}->setrgbcolor($self->{Colors}->rgbColor($self->{annotation}{color}));
	$self->{ps}->fill();

}

1;
