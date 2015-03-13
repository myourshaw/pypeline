#	Set.pm
#Wed Aug 13 18:03:31 MET DST 1997

package	Set;
require	5.000;
require	Exporter;

@ISA		= qw(Exporter);

@EXPORT		= qw(&intersection &minus &product &union &pair &substitute &productJoin &join2 &makeHash &dictWithKeys &mergedHashFromHash &mergeDict2dict &arrayFromKeys &mergeDict2dictDeeply &deepCopy &valuesForKeys &readHeadedTable &readHeadedTableHandle &writeHeadedTable &productT &productTL &arrayIsEqualTo &stripWhiteSpaceForColumns &sum &max &min &scaleSetTo &dictFromDictArray &toList &definedArray &firstDef &compareArrays &inverseMap &dictIsContainedInDict &keysOfDictLevel &sortTextNumber &readUnheadedTable);

# return all keys from a dict and from dicts contained therein up to the level of $level

# for strings of the form .*[^\d]\d* sorts according to string
# then according to numerical value of \d*
sub sortTextNumber { my ($a, $b) = @_;
	my ($at, $an) = ($a =~ m{^(.*?)(\d*)$});
	my ($bt, $bn) = ($b =~ m{^(.*?)(\d*)$});
	return $at cmp $bt if ($at cmp $bt);
	return $an <=> $bn;
}

sub keysOfDictLevel { my ($d, $level) = @_;
	my ($keys, $i) = ([]);

	return $keys if ($level < 0);
	foreach $key (keys %{$d}) {
		$keys = union($keys, keysOfDictLevel($d->{$key}, $level - 1));
	}
	return union($keys, [keys %{$d}]);
}

# compare arrays numerically by means of the '==' operator

sub compareArrays { my (@arrays) = @_;
	my ($i, $j, $f);
	for ($i = 1, $f = int(@{$arrays[0]}); $i < @arrays; $i++) {
		return 1 if ($f != @{$arrays[$i]});
	}
	for ($i = 0; $i < @{$arrays[0]}; $i++) {
		for ($j = 1, $f = $arrays[0]->[$i]; $j < @arrays; $j++) {
			return 1 if ($f != $arrays[$j]->[$i]);
		}
	}
	return 0;
}
sub dictIsContainedInDict { my ($d1, $d2) = @_;
	foreach $k (keys %{$d1}) {
		return 0 if ($d1->{$k} ne $d2->{$k});
	}
	return 1;
}

sub firstDef { my @args = @_;
	foreach $a (@args) { return $a if (defined($a)); }
	return undef;
}
sub toList { my ($array) = @_;
	return @{$array};
}
sub definedArray { my ($array) = @_;
	my $ret = [];
	map { push(@{$ret}, $_) if (defined($_)) } @{$array};
	return $ret;
}

sub	max { my (@arr) = @_;
	my $max = $arr[0];
	foreach $el (@arr)
	{	$max = $el if ($el > $max);
	}
	return $max;
}
sub	min { my (@arr) = @_;
	my $min = $arr[0];
	foreach $el (@arr)
	{	$min = $el if ($el < $min);
	}
	return $min;
}
sub	sum { my ($arr) = @_;
	my $sum = 0;
	foreach $el (@{$arr})
	{	$sum += $el;
	}
	return $sum;
}
#	round by biggest decimal places
sub scaleSetTo { my ($arr, $count) = @_;
	do {
		my $sum = sum($arr);
		my ($i, $new, $rest) = (0, [], []);
		foreach $e (@{$arr})
		{
			push(@{$new}, int($e / $sum * $count));
			push(@{$rest}, [$i++, $e / $sum * $count - int($e / $sum * $count)]);
		}
		@{$rest} = sort { $b->[1] <=> $a->[1] } @{$rest};
		my $share = $count - sum($new);
		foreach $e (@{$new})
		{
			last if (!$share--);
			$new->[$e->[0]]++;
		}
		$arr = $new;
	} until (sum($arr) >= $count);
	return $arr;
}

sub intersection { my($arr1,$arr2)=@_;
	my ($ret,%keys)=[];
	foreach $i (@{$arr1}) { $keys{$i}=0; }
	foreach $i (@{$arr2}) { push(@{$ret},$i) if (defined($keys{$i})); }
	return $ret;
}
sub minus { my($minuend,$subtrahend)=@_;
	my ($ret, %keys) = ([], ());
	foreach $i (@{$subtrahend}) { $keys{$i}=0; }
	foreach $i (@{$minuend}) { push(@{$ret},$i) if (!defined($keys{$i})); }
	return $ret;
}	

sub product { my($arr1,$arr2,$eliminate)=@_;
	my ($ret,$i,$j)=[];
	foreach $i (@{$arr1})
	{ foreach $j (@{$arr2}) { push(@{$ret},$i,$j) if (!$eliminate || ($i ne $j)); } }
#	print "Product:",join(':',@{$ret}),"|",join(':',@{$arr1}),"|",join(':',@{$arr2}),"\n";
	return $ret;
}
sub productT { my($arr1,$arr2,$eliminate)=@_;
	my ($ret,$i,$j)=[];
	foreach $i (@{$arr1})
	{ foreach $j (@{$arr2}) { push(@{$ret}, [$i, $j]) if (!$eliminate || ($i ne $j)); } }
	return $ret;
}
sub productTL { my($arr1,$arr2,$eliminate)=@_;
	my ($ret,$i,$j)=[];
	foreach $i (@{$arr1})
	{ foreach $j (@{$arr2}) { push(@{$ret}, [@{$i}, $j]) if (!$eliminate || ($i ne $j)); } }
	return $ret;
}
sub productJoin { my($arr1,$str,$arr2,$eliminate)=@_;
	my ($ret,$i,$j)=[];
	foreach $i (@{$arr1})
	{ foreach $j (@{$arr2}) { push(@{$ret},$i.$str.$j) if (!$eliminate || ($i ne $j)); } }
#	print "Product:",join(':',@{$ret}),"|",join(':',@{$arr1}),"|",join(':',@{$arr2}),"\n";
	return $ret;
}
sub pair { my($arr1,$arr2)=@_;
	my ($ret,$i,$j)=[];
	for ($i=0; $i<=$#$arr1; $i++) { push(@{$ret},$arr1->[$i],$arr2->[$j]); }
	return $ret;
}
sub substitute { my($arr1,$map)=@_;
	foreach $i (@{$arr1}) { $i=$map->{$i} if (defined($map->{$i})); }
	return $arr1;
}

# join together chunks of $cnt from @arr with $str1
# join together the chunks with $str2
sub join2 { my($str1, $cnt, $str2, @arr)=@_;
	my @temp;
	while ($#arr>=0)
	{	push(@temp,join($str1,splice(@arr,0,$cnt)));
	}
	return join($str2,@temp);
}

# return an array wich is a union of all arrays referenced by $arrOfArr
# the sematics is mathematical by excluding dupblicates

sub union { my (@arrOfArr) = @_;
	my ($dict, $ret);
	foreach $arr (@arrOfArr)
	{
		foreach $el (@{$arr})
		{
			next if (defined($dict->{$el}));
			$dict->{$el} = 0;
			push(@{$ret}, $el);
		}
	}
	return $ret;
}

#	other functions

sub makeHash { my ($keys, $values, $omitKey)=@_;
	my ($hash, $i) = {};
	for ($i = 0; $i <= $#$keys; $i++)
	{	$hash->{$keys->[$i]} = $values->[$i]
			if (defined($keys->[$i]) && (!defined($omitKey) || ($keys->[$i] ne $omitKey)));
	}
	return $hash;
}
sub dictWithKeys { my($keys)=@_;
	my $hash={};
	foreach $el (@{$keys})
	{	$hash->{$el}=0;
	}
	return $hash;
}

sub mergedHashFromHash { my($hash, $merge)=@_;
	my $ret={};

	foreach $key (keys(%{$hash}))
	{	$ret->{$key}=$hash->{$key};
	}
	foreach $key (keys(%{$merge}))
	{	$ret->{$key}=$merge->{$key};
	}
	return $ret;
}
sub mergeDict2dict { my($hashS, $hashD)=@_;
	foreach $key (keys(%{$hashS}))
	{	$hashD->{$key}=$hashS->{$key};
	}
	return $hashD;
}

#	<!> not cycle proof
#	array merge is supported by appending arrays to each other
#	if arrays are to be merged the arrays var contains respective keys
#	indicating the order of merging 0: hashS after hashD; 1: hashD after hashS
#
sub mergeDict2dictDeeply { my($hashS, $hashD, $arrays)=@_;
	foreach $key (keys(%{$hashS}))
	{	if (ref($hashS->{$key}) eq 'HASH' && ref($hashD->{$key}) eq 'HASH')
		{	mergeDict2dictDeeply($hashS->{$key}, $hashD->{$key});
		} elsif (defined($arrays) && defined($arrays->{$key}) &&
			ref($hashS->{$key}) eq 'ARRAY' && ref($hashD->{$key}) eq 'ARRAY')
		{
			if (!$arrays->{$key})	# natural order
			{	push( @{$hashD->{$key}}, @{$hashS->{$key}} );
			} else {
				unshift( @{$hashD->{$key}}, @{$hashS->{$key}} );
			}
		} else { $hashD->{$key}=$hashS->{$key}; }
	}
	return $hashD;
}
#	<!> not cycle proof
sub deepCopyHash { my ($hash)=@_;
	my ($copy)={%{$hash}};
	foreach $key (keys(%{$copy})) { $copy->{$key}=deepCopy($copy->{$key}); }
	return $copy;
}
sub deepCopyArray { my ($array)=@_;
	my ($copy)=[@{$array}];
	foreach $el (@{$copy}) { $el=deepCopy($el); }
	return $copy;
}
sub deepCopy { my ($el)=@_;
	if (ref($el) eq 'HASH')
	{	return deepCopyHash($el);
	} elsif (ref($el) eq 'ARRAY')
	{	return deepCopyArray($el);
	} else { return $el; }	# not copied <N>
}

sub arrayFromKeys { my($hash, $keys, $default) = @_;
	my $arr = [];

	foreach $el (@{$keys}) {
		push(@{$arr}, defined($hash->{$el})? $hash->{$el}: $default);
	}
	return $arr;
}

sub readHeadedTableHandle { my ($fileHandle, $defaultEntry, $providedFactorList) = @_;
	my ($list, $entry, $valueList, $rawList) = ([], {});
	my $factorList = [@{$providedFactorList}];

	if (!defined($providedFactorList))
	{
		$_=<$fileHandle>, chop($_);
		chop($_) if (substr($_, -1, 1) eq "\r");	# ms-dos line separation
		# destill factors and substitue according to $map, strip surrounding space
		$factorList = [map { /^\"(.*)\"$/o? $1: $_; } split(/\t/)];
	}

	# initialize $valueList
	$valueList = { @{product($factorList, [{}])} };

	while (<$fileHandle>)
	{
		# <t> s/[\n\r]+
		chop($_);
		chop($_) if (substr($_, -1, 1) eq "\r");	# ms-dos line separation
		$rawList = [map { /^\"(.*)\"$/o? $1: $_; } split(/\t/)];
		mergeDict2dict($defaultEntry, $entry = {});
		push(@{$list}, mergeDict2dict(makeHash($factorList, $rawList), $entry));
		foreach $factor (@{$factorList}) {
			$valueList->{$factor}{$entry->{$factor}} = 0;
		}
	}
	push(@{$factorList}, sort keys %{$defaultEntry});
	my $return = { factors => $factorList, values => $valueList, data => $list};
	$return->{list} = $return->{data};	# retain compatibility
	return $return;
}

sub readHeadedTable { my ($path, $defaultEntry, $providedFactorList) = @_;

	open(STAT_INPUT, $path);
	my $sets = readHeadedTableHandle(\*STAT_INPUT, $defaultEntry, $providedFactorList);
	close(STAT_INPUT);
	return $sets;
}
sub stripWhiteSpaceForColumns { my ($set, $columns) = @_;
	foreach $row (@{$set->{list}})
	{
		foreach $column (@{$columns})
		{
			$row->{$column} =~ s{^\s*(.*?)\s*$}{$1}o;
		}
	}
}

sub readUnheadedTable { my ($path) = @_;
	my $ret = [];

	open(TABLE_INPUT, $path);
	while (<TABLE_INPUT>)
	{
		chop($_);
		push(@{$ret}, [split(/\t/)]);
	}
	close(TABLE_INPUT);

	return $ret;
}

sub writeHeadedTableHandle { my ($fileHandle, $sets, $factorList, $options) = @_;
	my $factors = defined($factorList)? $factorList: $sets->{factors};

	print $fileHandle join("\t", defined($sets->{printFactors})?
		@{$sets->{printFactors}}: @{$factors}),"\n";

	foreach $data (@{$sets->{list}})
	{	my $dmy = arrayFromKeys($data, $factors);
		if ($options->{quoteSpace})
		{
			foreach $el (@{$dmy})
			{	$el = '"'.$el.'"' if ($el =~ m{\s}o);
			}
		}
		print $fileHandle join("\t", @{$dmy}), "\n";
	}
}

sub writeHeadedTable { my ($path, $sets, $factorList, $options) = @_;
	my $outputFile;
	if ($path eq 'STDOUT') { $outputFile = \*STDOUT; }
	else { open (DATA_OUTPUT,">$path"), $outputFile = \*DATA_OUTPUT; }
		writeHeadedTableHandle($path, $sets, $factorList, $options);
	close(DATA_OUTPUT) if ($path ne 'STDOUT');
}

sub arrayIsEqualTo { my ($arr1, $arr2) = @_;
	return 0 if ($#$arr1 != $#$arr2);
	my $i;
	for ($i = 0; $i <= $#$arr1; $i++) { return 0 if ($arr1->[$i] ne $arr2->[$i]); }
	return 1;
}

# compute a somewhat inverse map:
# foreach dict take the value of some key and associate that value with the dict

sub dictFromDictArray { my ($array, $key) = @_;
	my $dict = {};

	foreach $entry (@{$array})
	{
		$dict->{$entry->{$key}} = $entry;
	}
	return $dict;
}

# compute a Map with value => key
# which of course is not unique

sub inverseMap { my ($dict) = @_;
	my $idict = {};

	foreach $k (keys %{$dict}) {
		$idict->{$dict->{$k}} = $k;
	}
	return $idict;
}

1;
