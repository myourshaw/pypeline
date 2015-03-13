#!/usr/local/bin/perl
#	CCtk.pl
#Mon Mar 24 15:12:50 CET 2003

use Getopt::Long;
use PropertyList;
use TempFileNames;
use TkForms;
use Set;
use Data::Dumper;

# config manipulation
sub setFieldForName { my ($c, $fieldValueName, $fieldValue, $fieldName) = @_;
	foreach $f (@{$c->{fields}}) {
		if ($f->{name} eq $fieldName) {
			$f->{$fieldValueName} = $fieldValue;
			last;
		}
	}
}

sub saveValues { my ($dict) = @_;
	my $stringRep = stringFromProperty($dict);
	my $savePath = firstDef($defaultsPath, map { -e $_? "$_/CCdefaults.cfg": undef() } ("$ENV{HOME}/Library/Configs", '.'));
	Log("Saving config file to: $savePath");
	Log("Defaults:\n$stringRep", 3);
	writeFile($savePath, $stringRep);
}

# action callback
sub doAction { my ($d, $mainWindow, $n) = @_;
	my $w = $d->{mainWindow};
	my $vs = $w->collectValues();
	SWITCH: {
		if ($d->{action} eq 'save') { saveValues($vs); last SWITCH; }
		if ($d->{action} eq 'cancel') { $w->destroy(); last SWITCH; }
		if ($d->{action} eq 'run') {
			my $flags = join(' ', ($vs->{mapfile} ne ''? "--labelMap $vs->{mapfile}": '')
				,($vs->{labelfile} ne ''? "--labelFile $vs->{labelfile}": ''));
			my $cmd = "perl -w ./coloredChromosomes.pl --chromosomeSpec \"$vs->{config}\" $flags "
				."--placement $vs->{placement} --o \"$vs->{output}\"";
			Log("$cmd", 1);
			system($cmd);

			my $split = splitPathDict($vs->{output});
			Log("ps --> pdf"), system("cd $split->{dir} ; ps2pdf $split->{file}")
				if ($vs->{pdf});

			my $showfile = $vs->{pdf}? "$split->{dir}/$split->{base}.pdf": $vs->{output};
			my $viewer = firstDef($vs->{viewer}, 'gv');
			# this is a bad hack to differetiate UNIX/Windows <!><%>
			my $detach = (-e '/usr')? ' &': '';
			Log("Showing..."), system("\"$viewer\" \"$showfile\"".$detach) if ($vs->{show});

			#$w->destroy();
			last SWITCH;
		}
	}
}

# this function is triggered if a field is modified and requests a field action
sub fieldAction { my ($d, $mainWindow, $n) = @_;
	my $w = $d->{mainWindow};

	Log("Field did change :".$n);
	SWITCH: {
		if ($n eq 'config') {
			my $c = propertyFromString(readFile($w->valueForField($n)));
			$w->setPropertiesForField('placement', [keys %{$c->{placements}}]);
		last SWITCH; }
	}
	return 1;
}


#main $#ARGV @ARGV %ENV
	$result = GetOptions('help' => \$help, 'h' => \$help, 'defaults=s' => \$defaultsPath);
	initLog(5);

	if ($help || !$result)
	{	printf("USAGE: %s\n", ($0 =~ m{/?([^/]*)$}o));
		exit(!$result);
	}

	my $plistFile = readFileFirstLocation("CCtk.cfg",
		(".", "$ENV{HOME}/Library/Configs", "/Library/Configs"));
	die "Config file CCtk.cfg not found" if (!defined($plistFile));
	$c = propertyFromString($plistFile);
	$c->{actionHandler} = \&doAction;
	$c->{fieldHandler} = \&fieldAction;

	my $defaultsFile = readFileFirstLocation("CCdefaults.cfg",
		('.',  "$ENV{HOME}/Library/Configs", "/Library/Configs"));
	my $defaults = defined($defaultsFile)? propertyFromString($defaultsFile): {};
Log($defaultsFile);
	my $ccconfig = propertyFromString(readFileFirstLocation(firstDef(
		$defaults->{chromosomesConfigPath}, 'humanChromosomes.cfg'),
		('.', './configs', './Configs',
			"$ENV{HOME}/Library/Configs", "$ENV{HOME}/Library/Configs/Genetics",
			'/Library/Configs', '/Library/Configs/Genetics')));

	setFieldForName($c, 'options', [keys %{$ccconfig->{placements}}], 'placement');

	$f = TkForms->newWithConfig($c);
	$f->setValues($defaults);
	$f->run();

exit(0);
