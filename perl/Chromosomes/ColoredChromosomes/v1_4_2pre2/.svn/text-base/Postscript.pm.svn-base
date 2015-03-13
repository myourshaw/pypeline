#
#	Postscript.pm
#Wed Jan  3 11:14:46 MET 2001
#

package Postscript;
require 5.000;

#
#	<§> dependencies
#

use TempFileNames;
use Set;

#
#	<§> data structures
#

#	<!> unsatisfactory
$templatePath = '/LocalDeveloper/Libraries/privatePerl/Postscript';

#%!PS-Adobe-2.0 F-2.0

$stdHeader = <<STD_HEADER;
%!PS-Adobe-2.0
%%Title: TITLE
%%Creator: AUTHOR 
%%CreationDate: DATE
%%DocumentFonts: FONTS
%%Pages: PAGECOUNT
%%BoundingBox: BBOX
%%Orientation: Portrait
MEDIA
%%EndComments

STD_HEADER

$stdFooter = <<STD_FOOTER;
%%Trailer
%%EOF
STD_FOOTER

%sizes = (
	A4 => { width => 596, height => 842 },
	A5 => { width => 842/2, height => 596 },
	A6 => { width => 596/2, height => 842/2 },
);

%stdcolors = (
	black => { r => 0, g => 0, b => 0 },
	white => { r => 1, g => 1, b => 1 }
);

sub mediaStringFor { my ($media) = @_;
	return sprintf("%%%%DocumentMedia: %s %d %d 0 () ()", $media,
		$sizes{$media}{width},
		$sizes{$media}{height}
	);
}

#
#	<§> standard implementation
#

sub newWithProperties { my ($class, $properties) = @_;
	# my $self = $class->SUPER::new($properties);
	# own initialization follows
	my $self = bless({}, $class);
	$self->{persistent} = uc($properties->{persistent}) eq 'YES';
	
	if (defined($properties->{path})) {
		$self->{file} = $self->{path} = $properties->{path};
		open($self->{file}, '>'.$self->{path});
	} else {
		$self->{file} = $properties->{file};
	}

	if ($self->{persistent}) {
		$self->{tmp} = tempFileName('/tmp/postscript_pm', '.eps');
		open($self->{tmp}, '>'.$self->{tmp});
		
	}
	$self->{properties} = $properties->{properties};
	# will be outomatically recognized: $self->{properties}{PAGECOUNT} = 1;
	return $self;
}

#
#	<§> method implementation
#

sub addPSpackage { my ($self, $name) = @_;
	if (!defined($self->{templates}->{$name}))
	{
		$self->{templates}->{$name} = 1;
		my $resourcePath = resourcePath("Postscript/$name.ps");
#print STDERR $resourcePath;
		$self->addProlog(readFile($resourcePath));
	}
}

# this is written after the construction of the file
sub writeHeader { my ($self) = @_;
	return if (defined($self->{headerWritten}));
	$self->{headerWritten} = 1;

	my $defaultDict = { TITLE => 'Automatically generated postscript',
		AUTHOR => 'The automagical postscript author',
		DATE => ''.localtime(time),
		BBOX => '0 0 596 842',	# DIN A4
		FONTS => join(' ', keys(%{$self->{fonts}})),
		PAGECOUNT => $self->{properties}{PAGECOUNT} - 1,
		MEDIA => ''
		};
	my $realDict = mergeDict2dict($self->{properties}, $defaultDict);
	print {$self->{file}} mergeDictToString($realDict, $stdHeader);
}

sub addProlog { my ($self, $prolog) = @_;
	$self->{prolog} .= $prolog;
}

sub write { my ($self, $text) = @_;
	# automatically start a new page <i><!> internal interface for PAGECOUNT
	$self->newpage() if (!defined($self->{properties}{PAGECOUNT}));
	if ($self->{persistent}) {
		print {$self->{tmp}} $text;
	} else {
		$self->{text} .= $text;
	}
}

sub newpage { my ($self) = @_;
	my $newPageCount = ++$self->{properties}{PAGECOUNT};
	$self->write("showpage\n") if ($self->{properties}{PAGECOUNT} > 1);
	$self->write("%%Page: $newPageCount $newPageCount\n");
}

sub printFooter { my ($self) = @_;
	my $realFooter = mergeDictToString(
		{ FONTS => join(' ', keys(%{$self->{fonts}})),
		  PAGECOUNT => $self->{properties}{PAGECOUNT} - 1 }, $stdFooter);
	$self->write($realFooter);
}

sub close { my ($self) = @_;
	$self->write("showpage\n");
	$self->printFooter();

	# do the actual file construction
	$self->writeHeader();
	print {$self->{file}} "%%BeginProlog\n";
	print {$self->{file}} $self->{prolog};
	print {$self->{file}} "%%EndProlog\n";

	if ($self->{persistent}) {
		close($self->{tmp});
		defined($self->{path})
		 || die "Can't use a filehandle in persistent mode\n";
		close($self->{file});
		system("cat $self->{tmp} >> $self->{file}");
	} else {
		print {$self->{file}} $self->{text} if (defined($self->{file}));
	}

	if (defined($self->{path})) {
		close($self->{file});
	} # otherwise the caller has to close file (he gave us a handle)
}

# $c describes the graphics context to be established prior to inclusion
#	txp: translate prior to all other transformations
#	typ: "

sub epsInclude { my ($self, $path, $c) = @_;
	$self->addPSpackage('epsInclude');

	$self->gsave();
	$self->write("BeginEPSF\n");
	# this path is not too complex
	# a clip path can be set up before
	# <i> callback function with information on the epsfile itself
	if (defined($c->{clipRect})) {
		$self->rect($c->{clipRect}{x}, $c->{clipRect}{y},
			$c->{clipRect}{width}, $c->{clipRect}{height});
		$self->clip();
	}
	$self->translate($c->{txp}, $c->{typ}) if (defined($c->{txp}));
	$self->rotate($c->{angle}) if (defined($c->{angle}));
	$self->scale($c->{sx}, $c->{sy}) if (defined($c->{sx}));
	$self->translate($c->{txa}, $c->{tya}) if (defined($c->{txa}));
	$self->write("%%BeginDocument: $path\n");
	$self->write(readFile($path));
	$self->write("%%EndDocument\n");
	$self->write("EndEPSF\n");
	$self->grestore();
}

sub arc { my ($self, $x, $y, $r, $strt, $stop) = @_;
	$strt = 0 if (!defined($strt));
	$stop = 360 if (!defined($stop));
	$self->write(sprintf("%.4f %.4f %.4f %.4f %.4f arc\n", $x, $y, $r, $strt, $stop));
}
sub arcfill { my ($self, $x, $y, $r, $strt, $stop) = @_;
	$self->arc($x, $y, $r, $strt, $stop);
	$self->fill();
}

sub arcn { my ($self, $x, $y, $r, $strt, $stop) = @_;
	$strt = 0 if (!defined($strt));
	$stop = 360 if (!defined($stop));
	$self->write(sprintf("%.4f %.4f %.4f %.4f %.4f arcn\n", $x, $y, $r, $strt, $stop));
}

sub moveto { my ($self, $x, $y) = @_;
	$self->write(sprintf("%.4f %.4f moveto\n", $x, $y, $r));
}

sub rlineto { my ($self, $x, $y) = @_;
	$self->write(sprintf("%.4f %.4f rlineto\n", $x, $y, $r));
}
sub lineto { my ($self, $x, $y) = @_;
	$self->write(sprintf("%.4f %.4f lineto\n", $x, $y, $r));
}
sub line { my ($self, $x, $y, $x1, $y1) = @_;
	$self->moveto($x, $y);
	$self->lineto($x1, $y1);
}
sub strokeline { my ($self, $x, $y, $x1, $y1) = @_;
	$self->line($x, $y, $x1, $y1);
	$self->stroke();
}

sub setgray { my ($self, $gray) = @_;
	$self->write(sprintf("%.4f setgray\n", $gray));
}

sub setlinewidth { my ($self, $width) = @_;
	$self->write(sprintf("%.4f setlinewidth\n", $width));
}

sub rect { my ($self, $x, $y, $w, $h) = @_;
	$self->moveto($x, $y);
	$self->lineto($x + $w, $y);
	$self->lineto($x + $w, $y + $h);
	$self->lineto($x, $y + $h);
	$self->closepath();
}
sub rectstroke { my ($self, $x, $y, $w, $h) = @_;
	$self->rect($x, $y, $w, $h);
	$self->stroke();
}
sub rectfill { my ($self, $x, $y, $w, $h) = @_;
	$self->rect($x, $y, $w, $h);
	$self->fill();
}

sub stroke { my ($self) = @_;
	$self->write("stroke\n");
}
sub fill { my ($self) = @_;
	$self->write("fill\n");
}

sub gsave { my ($self) = @_;
	$self->write("gsave\n");
}
sub grestore { my ($self) = @_;
	$self->write("grestore\n");
}
sub clip { my ($self) = @_;
	$self->write("clip\n");
}
sub newpath { my ($self) = @_;
	$self->write("newpath\n");
}
sub closepath { my ($self) = @_;
	$self->write("closepath\n");
}

# <i> to be localized in Locale.pm
%psencoding = ( "ä" => 228, 'ö' => 246, 'ü' => 252, 'ß' => 223,
	'Ä' => 196, 'Ö' => 214, 'Ü' => 220 );
sub reencode { my ($self, $text) = @_;
	$text =~ s{.}{ord($&) >= 128? sprintf("\\%o", $psencoding{$&}): $&}seg;
	return $text;
}

sub show { my ($self, $text, $flipped) = @_;
	if ($flipped) {
		$self->gsave();
		$self->scale(1, -1);
	}
	$self->write(sprintf("(%s) show\n", $self->reencode($text)));
	$self->grestore() if ($flipped);
}

sub showxy { my ($self, $text, $x, $y, $flipped) = @_;
	$self->moveto($x, $y);
	$self->show($text, $flipped);
}

sub centeredShow { my ($self, $text, $x, $y, $flipped) = @_;
	$self->addPSpackage('centeredText');
	$self->write(sprintf("(%s) dup stringwidth pop 2 div %.4f exch sub ".
		"%.4f 2 index charheight exch pop 2 div sub moveto\n", $text, $x, $y));
	if ($flipped) {
		$self->gsave();
		$self->scale(1, -1);
	}
	$self->write("show\n");
	$self->grestore() if ($flipped);
}
sub ycenteredShow { my ($self, $text, $x, $y) = @_;
	$self->addPSpackage('centeredText');
	$self->write(sprintf("(%s) %.4f ".
		"%.4f 2 index charheight exch pop 2 div sub moveto show\n", $text, $x, $y));
}

# show right aligned
sub showr { my ($self, $text, $x, $y, $flipped) = @_;
	$text = $self->reencode($text);
	$self->write(sprintf("(%s) stringwidth pop %.4f exch sub ".
		"%.4f moveto\n", $text, $x, $y));
	$self->show($text, $flipped);
}

#	show rightaligned, centered on y-coordinate
sub showRyC { my ($self, $text, $x, $y) = @_;
	$self->addPSpackage('centeredText');
	$self->write(sprintf("(%s) dup stringwidth pop %.4f exch sub ".
		"%.4f 2 index charheight exch pop 2 div sub moveto show\n", $text, $x, $y));
}

sub showRotated { my ($self, $text, $x, $y, $angle) = @_;
	$self->gsave();
	$self->moveto($x, $y);
	$self->rotate($angle);
	$self->show($text);
	$self->grestore();
}

sub setfont { my ($self, $fontname, $size, $noEnscriptEncoding) = @_;
	if (ref($fontname) eq 'HASH' || ref($size) eq 'HASH')
	{	my ($font, $default) = ($fontname, $size);
		($fontname, $size) = defined($font)?
			($font->{name}, $font->{size}): ($default->{name}, $default->{size});
		$self->setrgbcolor($font->{color}) if (defined($font->{color}));
	}
	# do we have to change the encoding? (to allow for umlauts e.g.)
	if ( !$noEnscriptEncoding && !$self->{fontProperties}{$fontname}{didChangeEncoding}) {
		$self->{fontProperties}{$fontname}{didChangeEncoding} = 1;
		$self->addPSpackage('enscriptEncoding');
		$self->write("/$fontname /${fontname}_ne MF\n");
	}

	$self->{fonts}->{$fontname} = 0;
	$fontname .= '_ne' if (!$noEnscriptEncoding);
	$self->write(sprintf("/$fontname findfont %.4f scalefont setfont\n", $size));

}

sub setrgbcolor { my ($self, $c) = @_;
	$c = $stdcolors{$c} if (defined($stdcolors{$c}));
	$self->write(sprintf("%.4f %.4f %.4f setrgbcolor\n", $c->{r}, $c->{g}, $c->{b}));
}

sub rotate { my ($self, $a) = @_;
	$self->write(sprintf("%.4f rotate\n", $a));
}
sub scale { my ($self, $sx, $sy) = @_;
	$self->write(sprintf("%.4f %.4f scale\n", $sx, $sy));
}
sub translate { my ($self, $x, $y) = @_;
	$self->write(sprintf("%.4f %.4f translate\n", $x, $y));
}

sub colorRect { my ($self, $c, $x, $y, $w, $h) = @_;
	$self->newpath();
	$self->setrgbcolor($c);
	$self->rect($x, $y, $w, $h);
	$self->fill();
}

sub setBlendedColorFromColorsAndFraction { my ($self, $c1, $c2, $f) = @_;
	$self->setrgbcolor(blendedColorFromColorsAndFraction($c1, $c2, $f));
}

sub blendedColorFromColorsAndFraction { my ($c1, $c2, $f) = @_;
	my $cr = {};
	foreach $component ('r', 'g', 'b')
	{
		$cr->{$component} = (1 - $f) * $c1->{$component} + $f * $c2->{$component};
	}

	return $cr;
}
sub hblendedRect { my ($self, $c1, $c2, $x, $y, $w, $h, $stepSize) = @_;
	my ($i, $count);
	$stepSize = .5 if (!defined($stepSize));
	# adjust $stepSize to exactly fit to the borders of the rect
	$count = int($w / $stepSize + .5), $stepSize = $w / $count;

	for ($i = 0; $i < $count; $i++)
	{
		$self->newpath();
		$self->setrgbcolor(blendedColorFromColorsAndFraction($c1, $c2, $i/($count - 1)));
		$self->rect($x + $i * $stepSize, $y, $stepSize, $h);
		$self->fill();
	}
}
sub setInterpolatedColor { my ($self, $ilist, $v) = @_;
	my ($color, $i, $j);
	# <!> we assume $colorConfig tb sorted
	# <!> the case of a single color in the list is not handled
	# search for the interpolation port; <i> binary search
	for ($i = 1; $i < @{$ilist}; $i++)
	{	last if ($v <= $ilist->[$i]{value});
	}
	$j = $i - 1;
	
	if ($i >= @{$ilist}) { # <N> level off at max
		$self->setrgbcolor($ilist->[$j]{color});
	} else {
		$self->setBlendedColorFromColorsAndFraction($ilist->[$j]{color}, $ilist->[$i]{color},
		 ($v - $ilist->[$j]{value}) / ($ilist->[$i]{value} - $ilist->[$j]{value}));
	}
}

# colors are blended from $c1 to $c2 to $c1
# $f is the fraction of $h where $c1 is to appear
sub hshiningRect { my ($self, $c1, $c2, $f, $x, $y, $w, $h, $stepSize) = @_;
	$self->hblendedRect($c1, $c2, $x, $y , $w * $f, $h, $stepSize);
	$self->hblendedRect($c2, $c1, $x + $w * $f, $y , $w - $w * $f, $h, $stepSize);
}

sub curve { my ($self, $x1, $y1, $x2, $y2, $x3, $y3, $x4, $y4) = @_;
	$self->moveto($x1, $y1);
	$self->write(sprintf("%.4f %.4f %.4f %.4f %.4f %.4f curveto\n", $x2, $y2, $x3, $y3, $x4, $y4));
}

# draw a rect with rounded corners where $r indicates the radius of corners

sub bubble { my ($self, $x, $y, $w, $h, $r) = @_;
	$self->moveto($x, $y + $r);
	# lower left corner
	$self->arc($x + $r, $y + $r, $r, 180, 270);
	# lower border
	$self->lineto($x + $w - $r, $y);
	# lower right corner
	$self->arc($x + $w - $r, $y + $r, $r, 270, 360);
	# right border
	$self->lineto($x + $w, $y + $h - $r);
	# upper right corner
	$self->arc($x + $w - $r, $y + $h - $r, $r, 0, 90);
	# upper border
	$self->lineto( $x + $r, $y + $h);
	# upper left corner
	$self->arc($x + $r, $y + $h - $r, $r, 90, 180);
	# left border
	$self->closepath();
}

#
# This method arranges pages on bigger pages such that a book
# can be build by cutting and folding the bigger pages
# e.g. print four A6 pages on an A4 page each, such that no rearrangment
# is necessary
#
# $c is the config dict:
#	log2up:	logarithm to base 2 of how many pages are to be tiled on a basic page
#	format:	size of a basic page
#

use POSIX;

sub printBook { my ($self, $countOfPages, $callback, $context, $c) = @_;
	if (!defined($sizes{uc($c->{format})})) {
		die "Undefined format: $c->{format}. Defined formats are "
			.join(' ', keys(%sizes))."\n";
	}
	my ($w, $h) = ($sizes{uc($c->{format})}->{width}, $sizes{uc($c->{format})}->{height});
	# portrait or landscape?
	my ($width, $height) = (($c->{log2up} % 1)? ($h, $w): ($w, $h));
	$self->{properties}{BBOX} = sprintf("0 0 %d %d", $width, $height);
	$self->{properties}{MEDIA} = mediaStringFor(uc($c->{format}));

	my $up = 1 << $c->{log2up};
	Log("Print book configured with $up pages per basic page.");
	my ($rows, $cols) = (POSIX::floor(sqrt($up)), POSIX::ceil(sqrt($up)));
	# size of the subpages tiled on a basic page
	my ($swidth, $sheight) = ($width / $cols, $height / $rows);

	Log("Printing: $rows rows, $cols cols.");
	my $ppbp = $rows * $cols;	# pages per basic page
	my ($i, $bpages, $jr, $jc);
	# $countOfPages must be rounded appropriately
	# $cop is the rounded version of $countOfPages we round to 2 * $ppbp dt the double
	# sided concept of this book printing <i> tb configurable
	my $pos = $countOfPages % (2 * $ppbp);	# page overshoot
	my $cop = $countOfPages + ($pos? (2 * $ppbp - $pos): 0);
	# my $cop = $countOfPages + (($countOfPages % $ppbp)? ($ppbp - ($countOfPages % $ppbp)): 0);
	Log("Pages to print: $cop");
	my $coph = $countOfPages / 2;

	for ($i = 0, $bpages = $cop / $ppbp; $i < $bpages; $i++) {
		$self->newpage() if ($i);

		$self->gsave();
		$self->setlinewidth(0);
		$self->setgray(.33333);
		for ($jr = 1; $jr < $rows; $jr++) { # iterate rows
			for ($jc = 1; $jc < $cols; $jc++) {	# iterate columns
				$self->strokeline($jc / $cols * $width, 0, $jc / $cols * $width, $height);
				$self->strokeline(0, $jr / $rows * $height, $width, $jr / $rows * $height);
			}
		}
		$self->grestore();

		for ($jr = 0; $jr < $rows; $jr++) { # iterate rows
			for ($jc = 0; $jc < $cols; $jc++) {	# iterate columns
				# which page in the sequence of 1 .. $countOfPages - 1
				# has to be printed at the current location $i, $jc, $jr?
				# $i is the major index which advances at steps as follows
				# We determine the minimum page number on the basic page (mpnpbp) 
				# to be printed:

				# e.g. the first sample page for 8 up
				# |-----------------------------------------|
				# | cop - 1	| 0		| cop - 3	| 2			|
				# |-----------------------------------------|
				# | cop - 5 | 4		| cop - 7	| 6			|
				# |-----------------------------------------|
				# we first print even pages as measured from mpnpbp
				# and afterwards odd pages which are to be printed on the back of
				# the even pages

				# mpnpbp(i) = ($i % ($cop/2)) * $ppbp + int(2 * $i / $cop)
				# we define the basic index (bi) to reach into the page set according to
				# jc and jr: bi = (mpnpbp + jc - 1 + jr * cols)
				# the actual page number can then be computed as follows:
				# pn(i, jc, jr) = (jc % 2)? bi : (cop - 1 - bi)

				# you, the source reader are herewith challenged to improve on this
				# code and its explanation

				my $odd = int(2 * $i / $bpages);	# indicate whether we are in the second
													# batch of uneven pages
				my $mpnpbp = ($i % ($bpages/2)) * $ppbp;
				my $bi = ($mpnpbp + $jc + $jr * $cols);
				my $pn = ($jc % 2)^($odd)? $bi + ($odd? 1: -1): ($cop - $pos - 1 - $bi);

				next if ( ((($jc % 2)^$odd && $pn >= $coph))
					|| (!($jc % 2)^$odd && $pn < $coph) );	# this is an overshoot page

				Log("pn: $pn mpnpbp: $mpnpbp bi: $bi odd: $odd  bp: $i jr: $jr jc: $jc", 3);

				# draw the page
				$self->gsave();
				# move to subpage origin
#				$self->translate($jc * $swidth, $jr * $sheight);
				my ($mleft, $mright) = ($c->{margins}{left}, $c->{margins}{right});
				if (defined($c->{margins}{inner})) {
					$mleft =   ($pn%2)? $c->{margins}{inner}: $c->{margins}{outer};
					$mright = !($pn%2)? $c->{margins}{inner}: $c->{margins}{outer};
				}
				my $layout = { width => $swidth - $mright - $mleft,
					height => $sheight - $c->{margins}{bottom} - $c->{margins}{top},
					margins => $c->{margins}
				};

				if (uc($c->{flip}) eq 'YES') {
					$self->scale(1, -1);
					$self->translate($jc * $swidth + $mright,
						- $height + $jr * $sheight + $c->{margins}{top});
				} else {
					$self->translate($jc * $swidth + $mright,
						$jr * $sheight + $c->{margins}{bottom});
				}
				$callback->($self, $pn, $layout, $context);
				$self->grestore();
			}
		}
	}
}

1;

=head1 NAME

Postscript.pm - A class to build of postscript files using a procedural perl interface

=head1 SYNOPSIS

 use Postscript;

 $ps = Postscript->newWithProperties($properties);

 $ps->lineto(100, 140); # for example
 ...
 $ps->close();


=head1 DESCRIPTION

Instantiate a Postscript object. I<$properties> is a hash which which for which the key I<file> indicates a path where to write the final Postscript file. Then issue arbitrary Postscript commands. Finally invoke the close method which writes a valid Postscript file including header and footer. Some more complex methods (like centered text) automatically add Postscript "functions" to the file which are included at most once. The following list documents the relevant implemented methods which for most the time correspond with a Postscript operator. Therefore more competent explenation of these methods can be saught for in various Postscript documentations.

=head2 Postscript code producing methods

=over 4

=item *

arc($x, $y, $r, $strt, $stop) - issue the postscript operator 'arc' with the arguments on stack

=item *

arcfill($x, $y, $r, $strt, $stop) - first issue the arc() method, then the fill() method

=item *

arcn($x, $y, $r, $strt, $stop) - issue the postscript operator 'arcn' (which draws the other way round than 'arc')

=item *

moveto($x, $y) - issue the postscript operator 'moveto'

=item *

rlineto($x, $y) - issue the postscript operator 'rlineto'

=item *

lineto($x, $y) - issue the postscript operator 'lineto'

=item *

line($x, $y, $x1, $y1) - issue the postscript operator 'line'

=item *

strokeline($x, $y, $x1, $y1) - issue the postscript operator 'line', then 'stroke'

=item *

setgray($gray) - issue the postscript operator 'setgray'

=item *

setlinewidth($width) - issue the postscript operator 'setlinewidth'

=item *

rect($x, $y, $w, $h) - issue the postscript operators 'moveto', 'lineto', 'closepath' as to produce a rectangle with lower left corner ($x, $y) and width I<$w>, height I<$h>.

=item *

rectstroke($x, $y, $w, $h) - issue rect(), then stroke()

=item *

rectfill($x, $y, $w, $h) - issue rect(), then fill()

=item *

stroke() - issue the postscript operator 'stroke'

=item *

fill() - issue the postscript operator 'fill'

=item *

gsave() - issue the postscript operator 'gsave'

=item *

grestore() - issue the postscript operator 'grestore'

=item *

clip() - issue the postscript operator 'clip'

=item *

newpath() - issue the postscript operator 'newpath'

=item *

closepath() - issue the postscript operator 'closepath'

=item *

show($text) - issue the postscript operator 'show'

=item *

centeredShow($text, $x, $y) - show the text I<$text> centered by first calculating text width and then moving accordingly

=item *

ycenteredShow($text, $x, $y) - a corresponding method to ceneteredShow which, this time centers vertically by calculating the height of the text

=item *

showRyC($text, $x, $y) - show the text right aligned and centered vertically

=item *

showRotated($text, $x, $y, $angle) - show the text rotated by the angle I<$angle>. The graphic context is saved before and restored after the operation

=item *

setfont($fontname, $size) - issue the postscript operators 'findfont' and 'setfont'.

=item *

setrgbcolor($c) - issue the postscript operator 'setrgbcolor' when I<$c> is a hash containing the keys I<r>, I<g> and I<b>

=item *

rotate($a) - issue the postscript operator 'rotate'

=item *

scale($sx, $sy) - issue the postscript operator 'scale'

=item *

translate($x, $y) - issue the postscript operator 'translate'

=item *

colorRect($c, $x, $y, $w, $h) - issue the methods B<newPath()>, B<setrgbcolor()>, B<rect()> and B<fill()>

=item *

setBlendedColorFromColorsAndFraction($c1, $c2, $f) - issue B<setrgbcolor()> with a color lying I<$f>-way between I<$c1> and I<$c2>.

=item *

hblendedRect($c1, $c2, $x, $y, $w, $h, $stepSize) - draw a rectangle starting with color I<$c1> lefthand and ending in color I<$c2> righthand. I<$stepSize> indicates how many points each interpolated rectangle is wide.

=item *

hshiningRect($c1, $c2, $f, $x, $y, $w, $h, $stepSize) - by the use of the method B<hblendedRect()> B<hshiningRect> first blends from color I<$c1> to color I<$c2> and then back to color I<$c1>. Where within the specified rect the color I<$c2> is drawn is determind by I<$f> giving a fraction relative to rectangle width

=item *

curve($x1, $y1, $x2, $y2, $x3, $y3, $x4, $y4) - issue the postscript operators 'moveto' (with arguments I<$x1> and I<$x2>) and 'curveto' (with the rest of the arguments)

=item *

bubble($x, $y, $w, $h, $r) - draw a rect with rounded corners the radius of which is I<$r>

=back




