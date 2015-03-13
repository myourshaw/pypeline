#!/usr/local/bin/perl

use Getopt::Long;
use PropertyList;
use TempFileNames;
use Postscript;
use Set;
use CCcolors;
use CCAitags;

# what are chromosome names allowed to look like?
$chromNameRE = '(?:\\d*)[A-Z]?';

sub psChromosomePath { my ($stream, $chromPos, $chromosome) = @_;
	my ($x, $y, $p, $q, $w, $size) = ($chromPos->{x}, $chromPos->{y},
		$chromPos->{p}, $chromPos->{q}, $chromPos->{w}, $chromPos->{size});
	my $hw = $w / 2.0;
	my $r = $w / 5.0;
	my $r1 = $r;
	$stream->newpath();								# define telomere region

	# <!><i> this requires the banding dict to be in p/q subgroups
	if (uc($chromosome->{acrocentric}) eq 'YES')
	{
		# the chromosome satellite is always the outermost band
		# there is only one single band beyond the stem
		my @pBands = sort { $b->{start} <=> $a->{start} }
			@{$chromosome->{pBanding}};
		my $satellite = $pBands[0];
		my $satelliteHeight = ($satellite->{stop} - $satellite->{start}) * $p;
		$r1 = $satelliteHeight / 3.0 if (3 * $r > $satelliteHeight);
		$stream->bubble($x - $hw, $y + $p * $satellite->{start},
			$w, $satelliteHeight, $r1);

		my $stem = $pBands[1]; #$chromosome->{bandingDict}->{pstem};
		my $stemHeight = ($stem->{stop} - $stem->{start}) * $p;
		my $stemStart = $stem->{start} * $p;
		$stream->gsave();
			$stream->setlinewidth(0);
			$stream->strokeline($x - $hw * 0.5,
				$y + $stemStart + $stemHeight * 0.3333,
				$x + $hw * 0.5, $y + $stemStart + $stemHeight * 0.3333);
			$stream->strokeline($x - $hw * 0.5,
				$y + $stemStart + $stemHeight * 0.6666,
				$x + $hw * 0.5, $y + $stemStart + $stemHeight * 0.6666);
		$stream->grestore();
		$p -= $satelliteHeight + $stemHeight;
	}
	# the chromosome does not have arms
	if (!defined($chromosome->{p})) {
		$stream->bubble($x - $hw, $y - $size, $w, $size, $r1);
	} else {
		$stream->bubble($x - $hw, $y, $w, $p, $r1);
		$stream->bubble($x - $hw, $y - $q, $w, $q, $r);
	}
}

sub	drawAnnotations { my ($ps, $gc, $annotations, $doClip) = @_;
	my $tag = $gc->{map}{positions}{$gc->{name}}{tag};
	foreach $annotation (@{$annotations}) {
		next if (!(grep { $tag && $tag eq $_ } @{$annotation->{tags}})
			&& (defined($tag) || int(@{$annotation->{tags}})));
		$ps->gsave();
			$ps->clip() if ($doClip);			# clip to chromosome path
			my $moduleName = "CCA$annotation->{name}";
			my $annotator = $moduleName->newWithStreamChromosomeConfig($ps,
				$gc->{map}{positions}{$gc->{name}},
				{%{$gc}, annotation => $annotation->{config}} );
			$annotator->draw();
		$ps->grestore();
	}
}

sub drawChromosome { my ($stream, $chromPos, $chromosome, $placementConfig,
	$labels, $labelConfig, $placementMap, $gc) = @_;

	$stream->gsave();
	psChromosomePath($stream, $chromPos, $chromosome);

	drawAnnotations($stream, $gc, $gc->{placement}{internalAnnotations}, 1);

	# set outline color if defined in the configuration
	$stream->setrgbcolor($Colors->rgbColor($placementConfig->{outlineColor}))
		if (defined($placementConfig->{outlineColor}));
	# set outline width if defined in the configuration
	$stream->setlinewidth($placementConfig->{outlineWidth})
		if (defined($placementConfig->{outlineWidth}));
	$stream->stroke();							# draw the chromosome outline

	# draw annotations outside or over the chromosome shape
	drawAnnotations($stream, $gc, $gc->{placement}{annotations}, 0);
	$stream->grestore();
}

#	creates a dict of chromosome names
#	also creates a band name dict inside each chromosome

sub chromosomeDictFromConfig { my ($config) = @_;
	my $chromosomes = {};
	foreach $chromosome (@{$config->{chromosomes}})
	{
		$chromosomes->{$chromosome->{name}} = $chromosome;

		# <!> this is a persistent change in the data structures of $config
		next if (defined($chromosome->{bandingDict}));

		# the banding dict contains only sizes scaled to the whole chromosome
		if (defined($chromosome->{pBanding})) {
			my ($p, $q) = ($chromosome->{p}, $chromosome->{q});
			my ($rp, $rq) = ($p / ($p + $q), $q / ($p + $q));
			foreach $pband (@{$chromosome->{pBanding}})
			{
				$chromosome->{bandingDict}{'p'.$pband->{name}} =
					{ %{$pband}, start => $rp - $pband->{start} * $rp,
								 stop => $rp - $pband->{stop} * $rp};
			}
			foreach $qband (@{$chromosome->{qBanding}})
			{
				$chromosome->{bandingDict}{'q'.$qband->{name}} =
					{ %{$qband}, start => $rp + $qband->{start} * $rq,
								 stop => $rp + $qband->{stop} * $rq};
			}
		}
		# chromosomes are specified in a uniform manner
		# p/q specification is not necessary
		# bands without arm specification are prefixed with ':'
		foreach $band (@{$chromosome->{banding}})
		{
			my $prefix = ($band->{name} =~ m{^[pq]}o)? '': ':';
			$chromosome->{bandingDict}{$prefix.$band->{name}} = $band;
		}

		my $doCreateCentromere = uc($chromosome->{doCreateCentromere}) eq 'YES';
		# if all band are presented in a single list and no p/q length
		# is specified we calculate it now this feature is turn on by
		# the 'doCreateCentromere' eq 'YES' entry in the chromosome dict
		if ($doCreateCentromere) {
			my $bd = $chromosome->{bandingDict};
			my @pArmSizes = sort { $a <=> $b }
				map { ($bd->{$_}{start}, $bd->{$_}{stop}) }
				grep { /^p/o } keys(%{$bd});
			my $rpsize = $pArmSizes[$#pArmSizes] - $pArmSizes[0];
			# create a pBanding array with names stripped from their arm name
			$chromosome->{pBanding} = [map { {name => substr($_, 1),
				start => ($bd->{$_}{start} - $pArmSizes[0]) / $rpsize,
				stop => ($bd->{$_}{stop} - $pArmSizes[0]) / $rpsize} }
				grep { /^p/o } keys(%{$bd})];
			my @qArmSizes = sort { $a <=> $b }
				map { ($bd->{$_}{start}, $bd->{$_}{stop}) }
				grep { /^q/o } keys(%{$bd});
			my $rqsize = $qArmSizes[$#qArmSizes] - $qArmSizes[0];
			$chromosome->{qBanding} = [map { {name => substr($_, 1),
				start => ($bd->{$_}{start} - $qArmSizes[0]) / $rqsize,
				stop => ($bd->{$_}{stop} - $qArmSizes[0]) / $rqsize} }
				grep { /^q/o } keys(%{$bd})];
			$chromosome->{p} = $chromosome->{size} * $rpsize;
			$chromosome->{q} = $chromosome->{size} * $rqsize;
			# undef size to trigger centromere drawing <%><!>
			#$chromosome->{size} = undef();
			Log("p/q size chrom. [$chromosome->{name}]:\n\t"
			   ."p:$chromosome->{p} q:$chromosome->{q}", 4);
		}
	}
	return $chromosomes;
}

sub iterateChromosomes { my ($placement, $iSpec) = @_;
#$chromSubSub, $chromSub,
#	$prepareGroup, $groupSub, $prepareLane, $laneSub,
#	$prepareSubgroup, $concludeSubgroup) = @_;
	my $d = sub {}; # dummy subroutine
	my $defaultImp = { prepareLane => $d, prepareGroup => $d, prepareSubgroup => $d,
		chromSubSub => $d, concludeSubgroup => $d, chromSub => $d, groupSub => $d,
		laneSub => $d, preparePage => $d, concludePage => $d
	};
	mergeDict2dict($iSpec, $defaultImp);
	$iSpec = $defaultImp;

	if (ref($iSpec->{lanes}) eq 'ARRAY') {
		# <p> page1 is created automagically
		# <!> we modify $iSpec and should copy before
		$iSpec->{pages} = { page1 => $iSpec->{lanes} };
	} else {
		$iSpec->{pages} = $iSpec->{lanes};
	}

	foreach $page (sort keys %{$iSpec->{pages}}) {
		$iSpec->{preparePage}->($page, $iSpec);
		foreach $lane (@{$iSpec->{pages}{$page}})
		{	$iSpec->{prepareLane}->($lane, $iSpec);
			foreach $group (@{$lane})
			{	$iSpec->{prepareGroup}->($group, $iSpec);
				foreach $chromosome (@{$group})
				{
					if (ref($chromosome) eq 'ARRAY')
					{
						$iSpec->{prepareSubgroup}->($chromosome, $iSpec);
						# we have reached a subgroup
						foreach $chrom (@{$chromosome})
						{
							$iSpec->{chromSubSub}->($chrom, $iSpec);
						}
						$iSpec->{concludeSubgroup}->($chromosome, $iSpec);
					} else {
						$iSpec->{chromSub}->($chromosome, $iSpec);
					}
				}
				$iSpec->{groupSub}->($group, $iSpec);
			}
			$iSpec->{laneSub}->($lane, $iSpec);
		}
		$iSpec->{concludePage}->($page, $iSpec);
	}
}

sub maxR { my ($hold, $new) = @_;
	$$hold = $new if (!defined($$hold) || $new > $$hold); }
sub minR { my ($hold, $new) = @_;
	$$hold = $new if (!defined($$hold) || $new < $$hold); }

sub determineChromosomeSizes { my ($placement, $lanes) = @_;
	my $pmap = {};	# the placement map
	my $i;
	my $chromosomes = chromosomeDictFromConfig($config);

	my $chromSubSub = sub {	my ($chromosome) = @_;
		my ($chromName, $tag) = ($chromosome =~ m{^($chromNameRE)(?:_(.*))?}o);
		my ($p, $q, $size) = ($chromosomes->{$chromName}{p},
			$chromosomes->{$chromName}{q}, $chromosomes->{$chromName}{size});
		# overall max
		maxR(\$pmap->{maxp}, $p);
		maxR(\$pmap->{maxq}, $q);
		maxR(\$pmap->{maxh}, firstDef($size, $p + $q));
		$chromosomes->{$chromName}{placeToBottom} = 1 if (!defined($q));
		# garuantee the existance of {size}
		$chromosomes->{$chromName}{size} = $p + $q
			if (!$chromosomes->{$chromName}{size});

		# store the page on which the largest chromosome resides to allow
		# extraction of the corresponding scaling factor
		$pmap->{maxhPage} = $page if ($pmap->{maxh} == firstDef($size, $p + $q));
		# per page max
		maxR(\$pmap->{$page}{maxp}, $p);
		maxR(\$pmap->{$page}{maxq}, $q);
		maxR(\$pmap->{$page}{maxh}, firstDef($size, $p + $q));
		# per line max
		maxR(\$pmap->{$page}{map}[$i]{maxp}, $p);
		maxR(\$pmap->{$page}{map}[$i]{maxq}, $q);
	};

#	my $lanePrepare = sub { ($maxp, $maxq) = (0, 0); };
	my $laneSub = sub { my ($lane) = @_;
		my ($lmxp, $lmxq) = ($pmap->{$page}{map}[$i]{maxp},
			$pmap->{$page}{map}[$i]{maxq});
		$pmap->{$page}{allMaxp} += $lmxp;
		$pmap->{$page}{allMaxq} += $lmxq;
		$pmap->{$page}{allMaxh} += $pmap->{$page}{maxh};

		$pmap->{$page}{map}[$i]{maxh} = Set::max(
			$pmap->{$page}{maxh}, # chromosomes are placed to the bottom of the lane
			$lmxp + $lmxq);

		$i++;
	};
	my $concludePage = sub { my ($pageName, $spec) = @_;
		Log("H:$pmap->{$page}{allMaxh} P+Q:".
			$pmap->{$page}{allMaxp} + $pmap->{$page}{allMaxq}, 3);
		$pmap->{$page}{maxExtend} = Set::max(
			$pmap->{$page}{allMaxh},
			$pmap->{$page}{allMaxp} + $pmap->{$page}{allMaxq});

		maxR(\$pmap->{maxExtend}, $pmap->{$page}{maxExtend});
		$pmap->{$page}{scalingFactor} = ( $placement->{paperSize}{height}
			- $placement->{borderTop}
			- $placement->{borderButtom}
			- @{$spec->{pages}{$pageName}} * 3 * $placement->{namesize})
			/ $pmap->{$page}{maxExtend};
#		$pmap->{$page}{maxh} = $pmap->{$page}{maxp} + $pmap->{$page}{maxq};
		# also compute a global scaling factor applicable to all chromosomes
		minR(\$pmap->{scalingFactor}, $pmap->{$page}{scalingFactor});
	};

	iterateChromosomes($placement,
		{
			chromSubSub => $chromSubSub, chromSub => $chromSubSub,
			laneSub => $laneSub,
			preparePage => sub { $i = 0; },
			concludePage => $concludePage,
			lanes => $lanes	# what to iterate
		}
	);
	return $pmap;
}

sub scalingFactor { my ($pconfig, $pmap, $page) = @_;
	return (uc($pconfig->{globalScalingFactor}) eq 'YES')?
		$pmap->{scalingFactor}: $pmap->{$page}{scalingFactor};
}

sub determineAbsoluteChromosomePositions { my ($pmap, $placement, $lanes, $config) = @_;
	my $chromosomes = chromosomeDictFromConfig($config);
	# $yb: bottom of chromosomes, $yt: centromere position
	my ($i, $x, $y, $yb, $yt, $subgroupX);

	my $preparePage = sub {
		$y = $placement->{paperSize}{height} - $placement->{borderTop};
		$i = 0;
	};
	my $lanePrepare = sub { my ($lane) = @_;
		my $sf = scalingFactor($placement, $pmap, $page);
		$x = $placement->{borderLeft};
		$yt = $y - $pmap->{$page}{map}[$i]{maxp} * $sf;	# centromere position
		$yb = $y - $pmap->{$page}{map}[$i]{maxh} * $sf;
		$y = $yb - 3 * $placement->{namesize};
	};

	# called from closure <!>
	my $genericChromosomePlacement = sub { my ($taggedName) = @_;
		my ($chromName, $tag) = ($taggedName =~ m{^($chromNameRE)(?:_(.*))?}o);
		my $c = $chromosomes->{$chromName};
		my $sf = scalingFactor($placement, $pmap, $page);
		my $ypos = $yt;	# by default placed to the centromere
		my $yBottom = $yb;

		if ($c->{placeToBottom}) {
			$ypos = $yb + firstDef($c->{size}, $c->{q}) * $sf;	# bottom line
		} else {
			$yBottom += ($yt - $yb) - $c->{q} * $sf;
		}
		$pmap->{$page}{positions}{$taggedName} = {
			chromosome => $c, tag => $tag, chromosomeName => $chromName, name => '',
			x => $x, y => $ypos, # centromere position
			yBottom => $yBottom,		# bottom of chromosome
			p => $c->{p} * $sf,
			q => $c->{q} * $sf,
			size => $c->{size} * $sf,
			w => $placement->{chromosomeWidth}
		};
		Log("Abs. Chromsize[$chromName]: ".$c->{size} * $sf, 5);
	};

	my $chromSub = sub { my ($taggedName) = @_;
		$genericChromosomePlacement->($taggedName);# genericChromosomePlacement();
		push(@{$pmap->{$page}{labelPositions}},
			{ text => $taggedName, x => $x, y => $y + 2 * $placement->{namesize}});
		$x += $placement->{spacing};
	};
	my $chromSubSub = sub {	my ($taggedName) = @_;
		$genericChromosomePlacement->($taggedName); #genericChromosomePlacement();
		$x += $placement->{subgroupSpacing};
	};
	my $startSubgroup = sub { $subgroupX = $x; };
	my $concludeSubgroup = sub { my ($group) = @_;
		# presume this exists
		my ($text) = ($group->[0] =~ m{^($chromNameRE)(?:_(.*))?}o);
		$x -= $placement->{subgroupSpacing};
		push(@{$pmap->{$page}{labelPositions}},
			{ text => $text, x => ($x + $subgroupX) / 2,
			  y => $y + 2 * $placement->{namesize}});
		$x += $placement->{spacing};
	};
	my $groupSub = sub { $x += $placement->{intergroupSpacing}; };

	iterateChromosomes($placement,
		{	prepareLane => $lanePrepare, prepareSubgroup => $startSubgroup,
			chromSubSub => $chromSubSub, concludeSubgroup => $concludeSubgroup,
			chromSub => $chromSub, groupSub => $groupSub, laneSub => sub { $i++; },
			preparePage => $preparePage,

			lanes => $lanes
		}
	);
}

sub placeChromosomes { my ($placement, $config, $lanes) = @_;
	my $chromosomes = chromosomeDictFromConfig($config);
	my $pmap = determineChromosomeSizes($placement, $lanes);
	determineAbsoluteChromosomePositions($pmap, $placement, $lanes, $config);
	return $pmap;
}

sub drawChromosomes { my ($stream, $gc) = @_;
	my ($config, $placementConfig, $placementMap, $colors) =
		($gc->{global}, $gc->{placement}, $gc->{map}, $gc->{colors});
	my $chromosomes = chromosomeDictFromConfig($config);

	# save streamstate
	$stream->gsave();
	# introduce a font for drawing labels
	$stream->setfont('Helvetica', $placementConfig->{namesize});

	# draw the background
	if (defined($placementConfig->{backgroundColor}))
	{
		$stream->colorRect($Colors->rgbColor($placementConfig->{backgroundColor}),
			0, 0,
			$placementConfig->{paperSize}{width},
			$placementConfig->{paperSize}{height});
	}

	# load all required modules
	foreach $annotation (
		@{$gc->{placement}{internalAnnotations}},
		@{$gc->{placement}{annotations}}) {
		Log("Annotating with module: $annotation->{name}", 2);
		my $moduleName = "CCA$annotation->{name}";
		eval "use $moduleName";
	}

	# draw the chromosomes
	foreach $placedChromosomeName (sort keys %{$placementMap->{positions}})
	{
		my $c = $placementMap->{positions}{$placedChromosomeName};
		Log(sprintf("Chromosome: $placedChromosomeName (%.2f, %.2f)",
			$c->{x}, $c->{y}), 1);
		drawChromosome($stream, $c, $c->{chromosome}, $placementConfig,
			$colors->{byChromosome}{$placedChromosomeName}, $colors,
			$placementMap,
				{ %{$gc}, name => $placedChromosomeName, position => $c });
	}

	my $textColor;
	if (defined($placementConfig->{backgroundColor})
	&& $placementConfig->{backgroundColor} eq 'black') {
		$textColor = $Colors->rgbColor('white');
	} else {
		$textColor = $Colors->rgbColor('black');
	}
	$stream->setrgbcolor($textColor);
	foreach $l (@{$placementMap->{labelPositions}})
	{
		$stream->centeredShow($l->{text}, $l->{x}, $l->{y});
	}

	$stream->grestore();
}

sub layoutAndDrawChromosomes { my ($stream, $gc) = @_;
	$gc->{gmap} = placeChromosomes($gc->{placement}, $gc->{global},
		$gc->{placement}{lanes});
	my @pages;

	if (ref($gc->{placement}{lanes}) eq 'ARRAY') {
		@pages = ( 'page1' ); # page1 is created automagically <p>
	} elsif (ref($gc->{placement}{lanes}) eq 'HASH') {
		@pages = sort keys %{$gc->{placement}{lanes}};
	} else {
		die "Lanes must be either specified by array or hash.\n";
	}

	my $i = 0;
	foreach $page (@pages) {
		$stream->newpage() if ($i++);
		# this option allows to draw a bounding rectangle to force a specific bbox
		if (defined($gc->{boundingBox})) {
			$stream->gsave();
			$stream->setrgbcolor(0, 0, 0);
			my ($x, $y, $w, $h) = ($gc->{boundingBox}
				=~ m{(\d+)\s+(\d+)\s+(\d+)\s+(\d+)}os);
			Log("Bounding Box imposed: $x $y $w $h", 2);
			$stream->rectstroke($x, $y, $w, $h);
			$stream->grestore();
		}
		$gc->{map} = $gc->{gmap}{$page};
		drawChromosomes($stream, $gc);
	}
}

sub sortLabels { my ($labels, $labelConfig) = @_;
	my $c;
	my $tag = $labelConfig->{tag};

	foreach $label (@{$labels->{labels}})
	{
		# this has been warned about before <!>
		next if (!defined($label->{p}));
		# accept positions either as 'X[pq]' or as 'FcX'
		# where X is the chromosome name, F is a float and [pq] indicates the arm
		if (!defined($c = ($label->{p} =~ m{^($chromNameRE)[pq]}ogs)[0]))
		{	$c = ($label->{p} =~ m{^(?:[01]|\d*\.\d+)c($chromNameRE)}ogs)[0];
		}
		if (defined($c)) {
			$c .= "_$tag" if (defined($tag) && $tag ne '');
			push(@{$labels->{byChromosome}{$c}}, $label);
		} else {
			Log("Warning: unknown chromosome spec for label: $label->{p}", 4);
		}
	}
}

# all labels are mapped and sorted

sub readLabels { my ($colorfile, $labelMapPath, $labelConfig) = @_;
	my ($colors, $labelMap);
	$colors = propertyFromString(readFile($colorfile))
		if (defined($colorfile) && -e $colorfile);
	$labelMap = propertyFromString(readFile($labelMapPath))
		if (defined($labelMapPath) && -e $labelMapPath);

	# map labels if a map is given
	if (defined($labelMap)) {
		foreach $l (@{$colors->{labels}}) {
			Log("Couldn't map $l->{l}", 1) if (!defined($labelMap->{$l->{l}}));
			$l->{p} = $labelMap->{$l->{l}}{p} if (defined($labelMap->{$l->{l}}));
		}
	}

	# sort labels in chromosome folders
	sortLabels($colors, $labelConfig);

	return $colors;
}

@defaultPaths = ('configs');

if ( defined($ENV{HOME}) ) {
	push @defaultPaths, "$ENV{HOME}/Dokumente/Genetics/configs";
	push @defaultPaths, "$ENV{HOME}/Documents/Genetics/configs";
}

if ( defined($ENV{CCCONFIGPATH}) ) {
	push @defaultPaths, "$ENV{CCCONFIGPATH}";
}

#main $#ARGV @ARGV %ENV

	$result = GetOptions('help' => \$help, 'h' => \$help, 'chromosomeSpec=s' => \$cSpec,
		'placement=s' => \$placementName, 'o=s' => \$output, 'labelFile=s' => \$colorfile,
		'labelMap=s' => \$labelMap,
		'labelPlacement=s' => \$labelPlacement, 'boundingBox=s' => \$bbox );

	# note: $labelPlacement is currently undefined, and has no role in the script beyond
	#       this point.  It is intended to be developed into an equivalent to the 'labels'
	#	key in the placements specifications.
	$labelPlacement = '';
	# note: the line above is solely to avoid spurious warnings about $labelPlacement not
	#	being used.  It can be safely removed at anytime.

	if ($help || !$result)
	{	printf("USAGE: %s [--chromosomeSpec path] [--placement name] [--labelFile path] ".
			"[--labelMap path] [--labelPlacement name]\n",
			($0 =~ m{/?([^/]*)$}o));
		exit(!$result);
	}
	initLog(6);

	
	($configFile = readFileFirstLocation(firstDef($cSpec, "./humanChromosomes.cfg"),
		@defaultPaths)) || die "No configfile found.\n";
	$config = propertyFromString($configFile);
	$placementName = 'standard' if (!defined($placementName));
	$output = "/tmp/chromosomes.ps" if (!defined($output));

	Log("Output file:$output", 1);
	Log("Specification file:$cSpec", 2);
	Log("Chromosome placement specification:$placementName", 2);

	open(PS, ">$output");
	$psStream = Postscript->newWithProperties( { file => \*PS,
		properties => { AUTHOR => 'coloredChromosomes.pl '.
			'(c) 2001-2004, Stefan Boehringer '.
			'<stefan.boehringer@uni-essen.de>' } } );

	$placementConfig = $config->{placements}{$placementName}
	|| die "Placement scheme $placementName not known.\n".
		"Known placements are: ".
		join(' ', sort keys(%{$config->{placements}})). "\n";

	$Colors = CCcolors->new();	# create a global object
	layoutAndDrawChromosomes($psStream, {
			placement => $placementConfig,
			global => $config,
			labels => readLabels($colorfile, $labelMap, $placementConfig->{labels}),
			Colors => $Colors,
			boundingBox => $bbox
		}
	);
	$psStream->close();
	close(PS);
	#system("open $output");
exit(0);

# Convert:
#	pod2man coloredChromosomes.pl > man1/coloredChromosomes.pl.1
# View:
#	man -M . coloredChromosomes.pl
# Combined:
#  pod2man --center "Genetics Lib" coloredChromosomes.pl > man1/coloredChromosomes.pl.1 ; man -M . coloredChromosomes.pl

=head1 NAME

coloredChromosomes.pl - A program to automize complex ideogram drawing

=head1 SYNOPSIS

coloredChromosomes.pl [B<--help>] [B<--chromosomeSpec>=I<path>] [B<--placement>=I<placement name>] [B<--labelFile>=I<labelFile>] [B<--labelMap>=I<path>] [B<--labelPlacement>=I<placement name>] [B<--o>=I<output path>]

=head1 DESCRIPTION

B<coloredChromosomes.pl> is a program designed to draw chromosomal ideograms. It reads a specification file as given by B<--chromosomeSpec> which describes a chromosomal layout and banding pattern for some organism. The option B<--placement> chooses a placement pattern for the chromosomes on the paper out of a set of patterns present in the configuration file. The contents of the configuration file is described below. B<--labelFile> denotes the path of an additional file which holds annotation information to draw within or alongside the chromomal shapes. For one such I<labelFile> the option B<--labelPlacement> specifies a section in the main configuration file which gives further options to control annotation. If annotations are complex B<--labelMap> enables a further level of indirection by representing chromosomal locations by name which can be referred upon in the I<labelMap>.

=head2 The main specification file (MSF)

In this section the format of the specification file is explained. It is recommended to modify the supplied specification file rather than building one from scratch. Also you may find it helpful to change various different parameters to understand their effect, if in doubt. The main specification file is in PropertyList format (see L<PropertyList>). A global, unnamed dictionary contains the list I<chromosomes> and the dictionaries I<placements>, I<labelPlacements>, I<bandNamePlacements>.

=head2 MSF: I<chromosomes>

The entry I<chromosomes> in the B<MSF> is a list, each entry of which describes a single chromosome. Each entry is a dictionary containing the keys I<name>, I<p>, I<q>, I<pBanding> and I<qBanding>. The value of I<name> is a string which is used to refer to this chromosome. Also this string is printed under chromosome in the final program output. I<p> is a number which specifies the length of the p-arm of this chromosomes in arbitrary units. I<q> is the length of the q-arm. The units are arbitrary, but have to be consistent over all chromosomes such that relative sizes are properly drawn. I<pBanding> and I<qBanding> refer to lists each, which describe the banding patterns of the respective chromosomal arms. Each list entry is a dictionary containing the keys I<name>, I<start>, I<stop> and I<color>. The following example describes a single chromosome:

 {
 	name = "1";
 	p = "375.3212";
 	q = "389.9321";
 	pBanding = (
 		{ name = "11"; color = centromere;
 			stop = "0.0122"; start = "0.0000"; },
 		{ name = "12"; stop = "0.0365";
 			start = "0.0122"; color = black; }
 	);
 	qBanding = (
 		{ name = "11"; color = centromere;
 			stop = "0.0211"; start = "0.0000"; }
 	);
 }

Each entry in either I<pBanding> or I<qBanding> is a dictionary and describes a single band. The keys are I<name>, I<color>, I<start> and I<stop>. Again, I<name> is used both for reference of the band and as a name listed in the figure produced. I<color> may be either a name of a color which is then defined elsewhere or a dictionary containing the keys I<r>, I<g> and I<b> specifying color in RGB-space (values range from 0 to 1). There are predefined colors with names I<black>, I<white>, I<grey> and I<centromere>. The color I<centromere> is special in that it doesn't refer to a plain color, but to a pattern, which can be used to hightlight the centromeric region.
Banding data for human chromosomes is supplied with the standard configuration file.

=head2 MSF: I<placements>

Another section in the B<MSF> describes the actual placement of chromosomes on a page. In the main dictionary the entry I<placements> refers to a dictionary which contains the names of arbitrary many placements (as keys) which refer to dictionaries containing the actual placement specification. This way one can collect many useful placement schemes in one configuration file and refer to them by the I<--placement> option. If no placement is specified by the I<--placement>-option the placement with name I<standard> is chosen.

=head2 MSF: I<placements:lanes>

The layout for chromosomes is specified, by grouping chromosomes together and designating a lane to diplay groups of chromosomes. The key used to specify the layout is I<lanes>.
Example:

	lanes = (
		( (1, 2, 3), (4, 5) ),
		( (6, 7, 8, 9, 10, 11, 12) ),
		( (13, 14, 15), (16, 17, 18) ),
		( (19, 20, 21, 22), (X), (Y) )
	);

I<lanes> is a list, for which each entry is a list itself, specifying a lane of chromosome. The above example will display chromosomes 1,2,3,4 and 5 in the first lane, chromosomes 6 to 12 in lane 2 and so on. The fact that each lane is a list of lists means that within each lane chromosomes are grouped. In the example chromosomes 1,2 and 3 are grouped together in lane 1 much as chromosomes 4 and 5 are (a further complication of subgrouping is explained below: Advanced topics).

=head2 MSF: I<placements:distances>

The actual sizes of chromosomes is determined by their relative sizes together with absolute sizes specifying borders and other elements. The following table lists keys of a placement which specify absolute sizes. All measurements are in points (1/72th of an inch).

   Table 1: Absolute sizes to determine chromosome placements

   Key                Default  Explanation of value
                      value
 ---------------------------------------------------------------------
   borderLeft         70       Left border of the printed page
   borderTop          70       Top border of the printed page
   borderButtom       70       Bottom border of the printed page
   chromosomeWidth    10       Width of a chromosome
   spacing            70       Space between chromosomes
   intergroupSpacing  40       Additional space between groups of
                               chromosomes
   paperSize          dict     A dictionary containing the keys
                               width and height of the whole page
                               e.g. paperSize = { width = 595;
                               height = 842; }; for a DIN A4 page
   outlineWidth       0        The width of the chromosome outline.
                               For postscript output a size of 0 means
                               that the line is drawn as small as the
                               output device permits ("1 pixel").
   namesize           14       The font size of chromosome names

Figure 1 illustrates most sizes graphically. Provided that only 3 chromosomes are displayed altogether within a sinlge lane and that chromosome 1 and 2 are grouped the following identities hold: [a]: borderTop, [b]: borderLeft, [c]: spacing, [d]: intergroupSpacing, [e]: borderBottom, [f]: chromosomeWidth. Note, that the right border is determined by the sizes [b], [c], [d] and [f] and therefore doesn't need to be specified.


              ^
              |                                    |[f]|
             [a]                                   \   /
              |                                     | |
             \/
             .-.           .-.                      .-.
<-   [b]  -> | | <- [c] -> | | <- [c] ->  <- [d] -> | | 
             | |           | |                      | |
             | |           | |                      | |
             `-'           `-'                      `-'
             .-.           .-.                      .-.
             | |           | |                      | |
             | |           | |                      | |
             | |           | |                      | |
             `-'           `-'                      `-'

              1             2                        3
              ^
              |
             [e]
              |
             \/

  Figure 1: Graphical illustration of size specifications

=head2 MSF: I<placements:colors>

Additionally to sizes, also colors may be specified. I<backgroundColor> and I<outlineColor> specify the relevant colors. Chromosome names are drawn with the color I<outlineColor>.


=head2 Annotations

In the context of this software package all graphical elements besides chromosomal shapes and names are called annotations. Annotations are furtherly subdivided into internal and external annotations. The difference is that internal annotations are confined to the chromosomal shape, whereas external annotations are free to draw anywhere on the page. This behaviour is garuanteed by the Postscript clip operator. Any annotation is performed by a specific Perl class. For comprehensiveness the readily supplied annotation classes are discussed here, rather than in the respective class documentation.
For each placement configuration there are two lists with keys I<internalAnnotations> and I<annotations> to denote which kinds of annotation should be applied. Each list contains dictionaries which contain the keys I<name> and I<config>. I<name> specifies the type of annotation and I<config> provides additional information forwarded to the annotation module. Possible annotations are enumerated below.

=head2 Annotations: banding

The I<banding> annotation is an internal annotation which draws a banding pattern within the chromosome shape as provided by the chromosome definitions (see B<MSF:> I<chromosomes>). The I<config> dictionary may contain the key I<drawBandingWithVisualEffect> when the value may be I<chrome>. This option imitates some chrome effect.

=head2 Annotations: bandingNames

This type of annotation draws band names to the left hand side of chromosomes. The config dictionary must contain the key I<config> where the value is a name which refers toanother dictionary. A key with this name must exist in the dictionary I<bandNamePlacements> which in turn is located in the main dictionary. This way a specification of how to place banding names can be shared amongst several chromosomal placement schemes.
The annotation dictionary therefore could read:

 {
 	name = bandingNames;
 	config = { config = standard; };
 }

The config dictionary referred to by name must contain the following keys: I<font>, I<distatnce>, I<color>, I<tags>. Font refers to a dictionary containing font name and size (e.g. 'font = { name = Helvetica; size = 4; }'). I<distance> denotes the space between chromosome and the band name (distance [a] in Fig. 2). I<color> specifies which color is used to draw the lines and names. Again it may be a name or a dictionary specifying RGB values of the color. Finally I<tags> specifies to which chromosomes band names should be attached (this is only relevant when using subgroups; s. below). For most cases the value 'tags = ("");' should suffice. The size [c] (Fig.2) is the length of the line connecting band name and middle of the band if an overlap has to be avoided. This distance can be specified by the optional I<labelInset> key and is chosen as 3 times the font size as a default value.


                          .-.
                          | |
            band1 <-[a]-> | |
   band2 ---------        | | 
         |<-[c]->|        `-'
                          .-.
                          | |
                          | |
                          | |
                          `-'

 Figure 2: Distances involved in placing band names.

=head2 Annotations: plain

A very simple kind of internal annotation is the I<plain> annotation. In the config dict it expects a key I<color> which denotes a color to fill the interior of the chromosome (for the key tags, see "Advanced Topics" below). For example:

 {
 	name = plain;
 	config = { color = grey; };
 	tags = ( "" );
 }

=head2 Annotations: itags

To highlight distinct chromosomal positions the internal annotation I<itag> may be used. The config dictionary contains the key I<defaultColorHeight> which denotes the vertical size of rectangular marks printed inside the chromosomal shape. This size is relative to the largest chromosomal height on the page. Since chromosomes are scaled according to the available space the relative size of marks will keep constant no matter to what size chromosomes are streched to.

 {	name = itags;
	config = {
		/* height in percent of maximal chromosome length */
		defaultColorHeight = 0.0005;
 	};
 	tags = ( "" );
 }

The information where to place marks and what color to give them is stored in a different file. This file is specified via the B<--labelFile> option. This file is in PropertyList format and has the following format. A main dictionary contains a list referred to by I<labels>, and the dictionaries call I<colors> and I<labelProperties>. The following example illustrates the file format:

 {
 	labels = (
 		{ v = "1.000000"; tag = "H";
 		  l = "D10S582"; p = "10p21.3:0.5"; },
 		{ v = "1.000000";
 		  l = "D17S953"; p = "10p21.3:0.6"; }
 	);
 	colors = {
 		interpolationList = (
 			{ value = "0"; color =
 			  { g = "1"; b = "0"; r = "1"; }; },
 			{ value = "0.05"; color =
 			  { g = "1"; b = "0"; r = "0"; }; },
 			{ value = "0.06"; color =
 			  { g = "0"; b = "0"; r = "1"; }; },
 			{ value = "1"; color =
 			  { g = "0"; b = "1"; r = "0"; }; }
 		);
 		defaultValue = "0";
 		specialColor = { g = "0"; b = "1"; r = "0"; };
 	};
 	labelProperties = { selectSmallerEqual = ".05"; };
 }

Each label is a dictionary with the keys I<v>, I<p>, I<l> and I<tag>. I<v> is a value between 0 and 1 and represents a quantitative information about a locus the position of which is specified by I<p>. I<p> may be either a chromosomal position in standard nomenclature or a numerical value (for a detailed explanation see B<Representing chromosomal locations>). I<l> and I<tag> are textual annotations not used for internal tagging. The I<color> section of the label file describes how values are mapped to colors. One entry in this dictionary is an I<interpolationList>. Each entry in that list is a dictionary with a key I<value> and a key I<color>. If the value of a label matches one of the list the respective color is chosen. If it is not in the list the color is linerly interpolated between the entries with values matching most closely. The I<defaultValue> designates the value in case of missing value for an label entry. I<specialColor> is chosen if the label value is bigger than 1. I<labelProperties> are only used for the I<etag> annotation.

=head2 Annotations: etags

The external annotation I<etags> use the same information as I<itags> to draw textual information at specified positions besides chromosomes. A specification for I<etags> annotations may look as follows:

 {
 	name = etags;
 	config = { config = standard; };
 	tags = ( "" );
 }

The value of I<config> within the I<config> dictionary specifies an entry within the global dictionary with key I<labelPlacements> (compare description of B<bandingNames> annotations). An exmaple is shown:

 standard = {
 	height = 6;
 	chromosomeDistance = 15;
 	curveInset = 3;
 	font = { name = Helvetica; size = 4; }
 	doPrintLabelText = YES;
 };

Figure 3 illustrates the distances involved in drawing the labels. [a] corresponds to I<curveInset>, [b] to height and [c] to I<chromosomeDistance>. The I<height> parameter allows to define a rectangular region centered around the lable to be free of other label texts. If two rectangular regions associated with different labels overlap the position of the labels is changed to eliminate the overlap. This is done by a locally optimizing algorithm which minimizes the amount of movement. A bezier curve is drawn which connects chromosomal position and label text. You can suppress label printing (just showing the bezier curves) by setting I<doPrintLabelText> to I<NO>.

  .-.                                    ____
  | |                                   /  ^
  | |                                __/   |
  | | <-[a]-> --------- <-[a]-> label__   [b]
  `-' <-         [c]         ->        \   |
  .-.                                   \_\/__
  | |                                    
  | |
  | |
  `-'

 Figure 3: Distances involved in placing labels.

=head2 Annotations: sideline

This annotation is to represent some numerical value ranging across chromosomal locations.

=head2 Representing chromosomal locations

The annotations I<itags> and I<etags> need to specify chromosomal locations. There are several options to give these positions. In the label file each label entry may contains the key I<p> which can have either of two formats. First it is possible to use "ChromosomeArmBandname:fraction" where I<Chromosome> must be contained in the I<chromosomes> dictionary, I<Arm> must be either 'p' or 'q', I<Bandname> must be contained in the band name dictionary of the chromosome. Finally I<fraction> indicates where within a band the position is to be localized. I<fraction> ranges from 0 to 1 with 0 indicating the centromeric start of the band and 1 denoting the most telomeric position within a band.
Alternatively the positions I<p> may be specified by a single fraction ranging over the whole chromosome. The format is as follows I<fraction> immidiately followd by 'c' follow by a chromosome name (e.g. '0.2172c11'). A fraction of 1 indicates the most telomeric position on the 'p' arm and a fraction of 0 points to the telomeric position of the 'q' arm.
A further level of indirection is possible. The option B<--mapFile> specifies a file containing positions only. A global dictionary contains labels as keys which point to a directory which must contain the key I<p> which specifies a position by either of the formats explicated above. An example is shown:

 {
 	D11S4111 = { p = "0.811944c11"; band = "11q23.1"; };
 	D11S4112 = { p = "0.999278c11"; band = "11q25"; };
 	D12S1070 = { p = "0.725909c12"; band = "12q23.1"; };
 }

If no position is given in the label file as supplied to the I<itags> and I<etags> annotations the program tries to look up the label name in the map file dictionary. This way you may generate a label file from some experimental data and keep that file separate from the positioning information.

=head2 Advanced topics: subgrouping

With the capabilities presented so far only simple ideograms are possible. For example it is difficult to produce diploid ideograms. To accomplish arrangements like this it is possible to establish a variant of the chromosomal placements (B<MSF: I<placements>>). The following examples illustrates the concept of subgrouping:

 lanes = (
 	( ((1, 1_A), (2, 2_A), (3, 3_A) ), ( (4, 4_A), (5, 5_A)) ),
 	( ((6, 6_A), (7, 7_A), (8, 8_A), (9, 9_A), (10, 10_A),
 	   (11, 11_A), (12, 12_A)) ),
 	( ((13, 13_A), (14, 14_A), (15, 15_A)), ((16, 16_A),
 	   (17, 17_A), (18, 18_A))),
 	( ((19, 19_A), (20, 20_A), (21, 21_A), (22, 22_A)),
 	  ((X, X_A)), ((Y, Y_A)) )
 );

In contrast with simple grouping, subgrouping uses lists within groups of chromomes. A new distance in the I<placement> dictionary called I<intergroupSpacing> specifies the distance between chromosomes in subgroups. Also you have the possibiliy to attach an alphanumerical tag with an underscore to a chromosome name. This is to later refer to specific chromosomes should several chromosomes with the same name should exist or a specific set of chromosomes need to be distinguished from the others (this feature is also possible when no subgrouping is used). In the list of annotations each dictionary contains an entry called I<tags> which is a list of strings. These strings designate tags of chromosomes to which the annotation in question should be applied. The string "" is used to specify plain chromsome names. In the above example the string "A" would apply to all second chromosomes within each subgroup.

=head1 INTERNALS

=head2 How the space is allotted

Horizontal space is dispersed by the absolute sizes decribed in the paragraph B<MSF: placements:distances>. This gives you fine control over the horizontal arrangement of chromosomes.
The vertical space, however, is allotted relative to chromosomal sizes in order to use all available space on the page. First all absolute vertical sizes are collected and subtracted from the available vertical space. Next the maximal chromosomal sizes of each lane are calculated. The size of the lanes is then calculated by being proportional to the relative sizes and filling the available space.

=head2 How to extend this software package

If you want to write new annotion classes you have to subclass the class CCAnnotation. The available methods are documented there. You should also have a look at he other annotation classes as example reference. Especially the I<plain> annotation is as simple as an annotation can be.

=head1 FILES

If not specified by the B<--chromosomeSpec> option, the specification file is sought for in the following positions:

=over 2

=item * 

F<./humanChromosomes.cfg>

=item * 

F<${HOME}/Dokumente/Genetics/configs/humanChromosomes.cfg>

=item * 

F<${HOME}/Documents/Genetics/configs/humanChromosomes.cfg>

=item * 

F<${CCCONFIGPATH}/humanChromosomes.cfg>

=back


