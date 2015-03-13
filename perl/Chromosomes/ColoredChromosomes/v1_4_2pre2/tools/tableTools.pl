#!/usr/local/bin/perl

use Getopt::Long;
use Set;
use PropertyList;
use TempFileNames;

sub plistColumn { my ($set, $columnName) = @_;
	my $indent = 2;
	print "{\n\t$columnName = (\n"; # header
	print "\t" x $indent;
	my $l; foreach $row (@{$set->{list}})
	{
		print ", " if ($l++);
		print "\"", $row->{$columnName}, "\"";
	}
	print "\n\t);\n}\n";			# footer
}

sub	requote { my ($set, $quoteChar) = @_;
	foreach $row (@{$set->{list}})
	{
		foreach $column (@{$set->{factors}})
		{
			my ($text) = ($row->{$column} =~ m{^[\'\"]?(.*?)[\'\"]?$}o);
			if ($text =~ m{\s}o)
			{	$row->{$column} = $quoteChar. $text. $quoteChar;
			} else {
				$row->{$column} = $text;
			}
		}
	}
	foreach $factor (@{$set->{factors}})
	{
		push(@{$set->{printFactors}},  ($factor =~ m{\s}o)?
			$quoteChar. $factor. $quoteChar: $factor);
	}
}

sub constraint { my ($set, $select) = @_;
	my $i, $r;
	$select =~ s/C_([A-Za-z_0-9]+)/\$r->{$1}/g;
	print STDERR "select expression: $select\n";

	for ($i = 0; $i < @{$set->{list}}; )
	{	my $r = $set->{list}->[$i];
		if (!eval $select)
		{	splice(@{$set->{list}}, $i, 1);
		} else { $i++; }
	}
}

sub sortSet { my ($set, $columns) = @_;
	my $cmp;
	my @columnList = split(/\s+/, $columns);
	print STDERR "Sort by: $column\n";
	$set->{list} = [sort { foreach $c (@columnList)
		{ return $cmp if ($cmp = $a->{$c} cmp $b->{$c}); } }
		@{$set->{list}}];
}

#main $#ARGV @ARGV %ENV
	$result = GetOptions('stripWhiteSpace' => \$strip,
		'plist' => \$plist,			# write
		'readPlist:s' => \$readPlist,
		'plistColumn=s' => \$column, 'h' => \$help, 'help' => \$help,
		'quoteSpace=s' => \$requote, 'select=s' => \$select, 'sortBy=s' => \$sortBy,
		'project=s' => \$pcolumns, 'supplementary=s' => \$supplementPath,
		'postprocess=s' => \$postprocess
	);

	if ($help || !$result || @ARGV < 1)
	{	printf("USAGE: %s [--plistColumn col] [--stripWhiteSpace] [--plist] tableFile\n".
			"\t[--quoteSpace QuoteChar] [--select 'expresssion'] [--sortBy column]\n".
			"\t[--project \"columns\"] [--postprocess 'code']\n".
			"\tExample: tableTools.pl --select 'C_firstColName =~ /^start/'\n".
			"\tExample: tableTools.pl --select '!defined(\$supplement->{".
			"\C_ColName})'\n".
			"\tExample: tableTools.pl --plist --postprocess '{ labels => PLIST }'\n".
			'',
			($0 =~ m{/?([^/]*)$}o));
		exit(!$result);
	}
	die "Table file '$ARGV[0]' doesn't exist\n" if (! -e $ARGV[0]);
	if (defined($readPlist)) {
		my $property = propertyFromString(readFile($ARGV[0]));
		$set->{list} = ref($property) eq 'HASH'? $property->{$readPlist}: $property;
		$set->{factors} = [keys %{$set->{list}[0]}];
	} else {
		$set = readHeadedTable($ARGV[0]);
	}
	# $set->{list} is garuanteed to be the table info
	# $set->{output} is usually considered to be the same, but may be hacked to
	# be different
	$set->{output} = $set->{list};
	$projection = defined($pcolumns)? [split(/\s+/, $pcolumns)]: $set->{factors};

	$supplement = readHeadedTable($supplementPath) if (defined($supplementPath));

	constraint($set, $select) if (defined($select));
	sortSet($set, $sortBy) if (defined($sortBy));

	if (defined($postprocess)) {
		$postprocess =~ s/PLIST/\$set->{list}/gos;
#print $postprocess;
		$set->{output} = eval($postprocess);
	}

	SWITCH: {
		defined($plist) && do {
			print stringFromProperty($set->{output}), "\n";
		last SWITCH };
		defined($column) && do {
			stripWhiteSpaceForColumns($set, [$column]) if ($strip);
			plistColumn($set, $column);
		last SWITCH; };
		defined($requote) && do {
			stripWhiteSpaceForColumns($set, $set->{factors});
			requote($set, $requote);
			writeHeadedTable('STDOUT', $set);
		last SWITCH; };
		do {	# default
			writeHeadedTable('STDOUT', $set, $projection);
		last SWITCH; };
	}

exit(0);
