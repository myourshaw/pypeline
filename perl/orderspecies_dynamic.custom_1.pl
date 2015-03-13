#!/usr/bin/perl -w
#
# Program to order the species found in the Y column of input file, according to evolutionary distance from human.
# Also order column X in the corresponding order. Finally, print how far the mutation occurred, i.e. at what species
#

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper; 
use FileHandle;
use List::Util qw(first);
use LWP::Simple;
use HTML::TokeParser;

# run the main sub routine
&main();

# this is the main sub-routine - it needs the configured $config hash
sub main {
	print "Opening input file...\n";

	my $infile = 'orthologues.vcf';
	my $outfile = '>orthologues.sorted.vcf';

	open (OUTFILE, $outfile);	# Open the output file
	open (INFO, $infile);		# Open the input file
	my @data; 
	while (<INFO>) { 
		chomp; 
		push @data, [split "\t"]; 
	} 
	close(INFO);			# Close the file

	#count number of lines in input file
	my $numLines;
	open (FILE, $infile) or die "Can't open '$infile': $!";
	$numLines++ while (<FILE>);
	close FILE;

	print "Building phylogenetic tree...\n";

	my $species = "homo_sapiens";		# human
	#my $species = "mus_musculus";		# mouse

	# dynamically get species common and latin names
	my $speciesListCommon = [ ];
	my $speciesListLatin = [ ];
	&getEnsemblSpecies(\$speciesListCommon, \$speciesListLatin);

	# generate distance array from phylogenetic tree
	my @distArray = &generatePhyloTreeDistArray($species, $speciesListCommon, $speciesListLatin);

	print "Processing...  0%";

	my $currLine = 0;
	my $columnNumber = -1;
	foreach my $line(@data) {

		#print out what percentage of processing we're done with
		$currLine ++;
		print "\b\b\b\b";
		printf("%3d", int($currLine * 100.0 / $numLines + 0.5));
		print "%";
		
		#skip header lines and print them to outpu file as is
		if (substr(@$line[0], 0, 2) eq "##") {
			print OUTFILE join("\t", @$line)."\n";
		}
		else {
			if (substr(@$line[0], 0, 1) eq "#") {
				print OUTFILE join("\t", @$line)."\n";
				if (substr(@$line[0], 0, 6) eq "#CHROM") {

					#find out which column INFO is
					for (my $i = 0; $i < scalar(@$line); $i++) {
						my $word = @$line[$i];
						if ($word eq "INFO") {
							$columnNumber = $i;
						}
					}
				}
			}
			else {	
				#check whether we know which column the INFO column is
				if ($columnNumber eq -1) {
					print "\nERROR: INFO column not found. Make sure the line with column headers starts with '#CHROM'.\n";
					exit;
				}
				
				#print columns that we're ignoring to OUTFILE directly
				for (my $i = 0; $i < $columnNumber; $i++) {
					print OUTFILE @$line[$i]."\t";
				}	

				my $info = @$line[$columnNumber];
				if ($info eq ".") {
					print OUTFILE ".";
				}
				else {
					$info =~ m/AAs_in_Orthologues=([A-Z,\-]+);Orthologues=([A-Za-z;_.]+)/;
					my $AAsInOrthologues = $1;
					my $orthologues = $2;

					#print "-----------------\nUNORDERED:\n";
					#print "$AAsInOrthologues\n";
					#print "$orthologues\n";

					# order the orthologues according to genetic distance from my species
					&orderOrthologues($species, \$AAsInOrthologues, \$orthologues, \@distArray, $speciesListCommon, $speciesListLatin);
		
					#print "-----------------\nORDERED:\n";
					#print "$AAsInOrthologues\n";
					#print "$orthologues\n";

					print OUTFILE "AAs_in_Orthologues=$AAsInOrthologues;Orthologues=$orthologues";

				}
				#print any columns that might follow after the INFO column
				for (my $i = ($columnNumber+1); $i < scalar(@$line); $i++) {
					print OUTFILE "\t".@$line[$columnNumber];
				}

				print OUTFILE "\n";
			}
		}
	}
	close(OUTFILE);
	print "\nDone.\n";
}

sub generatePhyloTreeDistArray {
	my $species = shift;
	my $speciesListCommon = shift;
	my $speciesListLatin = shift;

	# assign arrays of species common and latin names
	#my $speciesListCommon = [ "Human", "Chimpanzee", "Gorilla", "Orangutan", "Gibbon", "Macaque", "Marmoset", "Tarsier", "Bushbaby", "Mouse_Lemur", "Tree_Shrew", "Guinea_Pig", "Kangaroo_Rat", "Mouse", "Pika", "Rabbit", "Rat", "Squirrel", "Alpaca", "Cat", "Cow", "Dog", "Dolphin", "Hedgehog", "Horse", "Megabat", "Microbat", "Panda", "Pig", "Shrew", "Armadillo", "Elephant", "Hyrax", "Lesser_hedgehog_tenrec", "Sloth", "Opossum", "Wallaby", "Platypus", "Anole_Lizard", "Chicken", "Turkey", "Zebra_Finch", "X.tropicalis", "Fugu", "Medaka", "Stickleback", "Tetraodon", "Zebrafish", "C.intestinalis", "C.savignyi", "Fruitfly", "C.elegans", "S.cerevisiae", ];
	#my $speciesListLatin = [ "Homo_sapiens", "Pan_troglodytes", "Gorilla_gorilla", "Pongo_abelii", "Nomascus_leucogenys", "Macaca_mulatta", "Callithrix_jacchus", "Tarsius_syrichta", "Otolemur_garnettii", "Microcebus_murinus", "Tupaia_belangeri", "Cavia_porcellus", "Dipodomys_ordii", "Mus_musculus", "Ochotona_princeps", "Oryctolgaus_cuniculus", "Rattus_norvegicus", "Spermophilus_tridecemlineatus", "Vicugna_pacos", "Felis_catus", "Bos_taurus", "Canis_familiaris", "Tursiops_truncatus", "Erinaceus_europaeus", "Equus_caballus", "Pteropus_vampyrus", "Myotis_lucifugus", "Ailurupoda_melanoleuca", "Sus_scrofa", "Sorex_araneus", "Dasypus_novemcinctus", "Loxodonta_africana", "Procavia_capensis", "Echinops_telfairi", "Choloepus_hoffmanni", "Monodelphis_domestica", "Macropus_eugenii", "Ornithorhynchus_anatinus", "Anolis_carolinensis", "Gallus_gallus", "Meleagris_gallopavo", "Taeniopygia_guttata", "Xenopus_tropicalis", "Takifugu_rubripes", "Oryzias_latipes", "Gasterosteus_aculeatus", "Tetraodon_nigroviridis", "Danio_rerio", "Ciona_intestinalis", "Ciona_savignyi", "Drosophila_melanogaster", "Caenorhabditis_elegans", "Saccharomyces_cerevisiae", ];

	# get index of our species
	my $speciesIndex = first { lc($speciesListLatin->[$_]) eq lc($species) } 0..scalar(@$speciesListLatin);

	# check defined species exists
	die("ERROR: Could not find species \"", $species, "\" in Ensembl species list\n") unless defined $speciesIndex;

	# use phylogenetic tree; branch lengths for Ensembl species phylogenetic tree obtained from http://tinyurl.com/ensembltree on 07/26/11
	#my $phyloTree = lc("(((((((((((((((((((((((Homo_sapiens:0.0067,Pan_troglodytes:0.006667):0.00225,Gorilla_gorilla:0.008825):0.00968,Pongo_abelii:0.018318):0.00717,Nomascus_leucogenys:0.025488):0.00717,(Macaca_mulatta:0.007853,?Papio_hamadryas:0.007637):0.029618):0.021965,Callithrix_jacchus:0.066131):0.05759,Tarsius_syrichta:0.137823):0.011062,(Microcebus_murinus:0.092749,Otolemur_garnettii:0.129725):0.035463):0.015494,Tupaia_belangeri:0.186203):0.004937,(((((Mus_musculus:0.084509,Rattus_norvegicus:0.091589):0.197773,Dipodomys_ordii:0.211609):0.022992,Cavia_porcellus:0.225629):0.01015,Spermophilus_tridecemlineatus:0.148468):0.025746,(Oryctolagus_cuniculus:0.114227,Ochotona_princeps:0.201069):0.101463):0.015313):0.020593,((((Vicugna_pacos:0.107275,(Tursiops_truncatus:0.064688,(Bos_taurus:0.061796,?Ovis_aries:0.061796):0.061796):0.025153):0.0201675,Sus_scrofa:0.079):0.0201675,((Equus_caballus:0.109397,(Felis_catus:0.098612,(Ailuropoda_melanoleuca:0.051229,Canis_familiaris:0.051229):0.051229):0.049845):0.006219,(Myotis_lucifugus:0.14254,Pteropus_vampyrus:0.113399):0.033706):0.004508):0.011671,(Erinaceus_europaeus:0.221785,Sorex_araneus:0.269562):0.056393):0.021227):0.023664,(((Loxodonta_africana:0.082242,Procavia_capensis:0.155358):0.02699,Echinops_telfairi:0.245936):0.049697,(Dasypus_novemcinctus:0.116664,Choloepus_hoffmanni:0.096357):0.053145):0.006717):0.234728,(Monodelphis_domestica:0.125686,Macropus_eugenii:0.122008):0.2151):0.071664,Ornithorhynchus_anatinus:0.456592):0.109504,((((Gallus_gallus:0.041384,Meleagris_gallopavo:0.041384):0.041384,Anas_platyrhynchos:0.082768):0.082768,Taeniopygia_guttata:0.171542):0.199223,Anolis_carolinensis:0.489241):0.105143):0.172371,Xenopus_tropicalis:0.855573):0.311354,(((Tetraodon_nigroviridis:0.224159,Takifugu_rubripes:0.203847):0.195181,(Gasterosteus_aculeatus:0.316413,Oryzias_latipes:0.48197):0.05915):0.32564,Danio_rerio:0.730752):0.147949):0.526688,?Petromyzon_marinus:0.526688),(Ciona_savignyi:0.8,Ciona_intestinalis:0.8)Cionidae:0.6)Chordata:0.2,(?Apis_mellifera:0.9,(((?Aedes_aegypti:0.25,?Culex_quinquefasciatus:0.25):0.25,?Anopheles_gambiae:0.5)Culicinae:0.2,Drosophila_melanogaster:0.8)Diptera:0.1)Endopterygota:0.7)Coelomata:0.1,Caenorhabditis_elegans:1.7)Bilateria:0.3,Saccharomyces_cerevisiae:1.9)Fungi_Metazoa_group:0.3)");
	
	# dynamically get tree with branch lengths
	my $url = 'http://tinyurl.com/ensembltree';
	my $content = get $url;
  	die "Couldn't get tree from $url" unless defined $content;
	my $phyloTree = lc($content);

	my @distanceToSpecies = ((-1) x scalar(@$speciesListLatin));
	my $specieslc;
	my $currspecieslc;

	# trim branches of tree that aren't in our species list
	my @trimmedtree = split(/([\(\):,])/, $phyloTree);
	for(my $i = 0; $i < scalar(@trimmedtree); $i++) {
		my $item = $trimmedtree[$i];
		$item =~ s/\?//;
		my $itemIndex = first { lc($speciesListLatin->[$_]) eq $item } 0..scalar(@$speciesListLatin);
		if (($item =~ m/[a-z_]+/)  && !(defined $itemIndex)) {
			#print "Item $item not found. Removing from tree.\n";
			my @replacement = ("","","");
			splice(@trimmedtree, $i, 3, @replacement);			
		}
	}
	$phyloTree = join("", @trimmedtree);
	#print $phyloTree."\n";

	# calculate distances between our species and every other species
	for(my $i = 0; $i < scalar(@$speciesListLatin); $i++) {
		
		# set distance to our own species to zero
		if ( $i == $speciesIndex) {
			$distanceToSpecies[$i] = 0;
		}
		else {
			# does our species and the current species exist in our phylogenetic tree?
			$specieslc = lc($species);
			$currspecieslc = lc($speciesListLatin->[$i]);
			if (($phyloTree =~ m/$specieslc/) && ($phyloTree =~ m/$currspecieslc/)) {
				#print "Found ".$specieslc." and ".$currspecieslc."\n";

				# make a copy of our tree
				# e.g. ((((T:0.7,F:0.3):0.1,(G:0.4,O:0.5):0.6):0.2;D:0.8):0.9;A:1.1)
				my $tree = $phyloTree;

				# trim all leaves of the tree that aren't our species or the current species
				# e.g. if our species are F and D --> ((((,F:0.3):0.1,(,):0.6):0.2;D:0.8):0.9,)
				for(my $j = 0; $j < scalar(@$speciesListLatin); $j++) {
					my $jspecies = lc($speciesListLatin->[$j]);
					if (($jspecies ne $specieslc) && ($jspecies ne $currspecieslc)) {
						$tree =~ s/$jspecies:\d.[\d]+//g;
					}
				}

				# trim all branches of the tree that aren't on the path between the two species of interest
				# e.g. if our species are F and D --> (((,F:0.3):0.1,):0.2;D:0.8)
				while ($tree =~ m/\(,\):\d.[\d]+/) {
					$tree =~ s/\(,\):\d.[\d]+//g;
				}
				while ($tree =~ m/\(,\)/) {
					$tree =~ s/\(,\)//g;
				}
				# find set of parentheses that encloses our species and discard everything outside
				# code adapted from http://www.perlmonks.org/?node_id=660316
				my @queue = ( $tree );			
				my $regex = qr/
					(			# start of bracket 1
					\(			# match an opening parentheses
						(?:               
						[^\(\)]++	# one or more non parentheses, non backtracking
							|                  
						(?1)		# recurse to parentheses 1
						)*                 
					\)			# match a closing parentheses
					)			# end of parentheses 1
					/x;

				$" = "\n\t";

				my @potentials;
				while( @queue )
				{
					my $string = shift @queue;
					my @groups = $string =~ m/$regex/g;
					#print "Found:\n\t@groups\n\n" if @groups;
					if (($string =~ m/$specieslc/) && ($string =~ m/$currspecieslc/)) {
						push @potentials, $string;
					}
					unshift @queue, map { s/^\(//; s/\)$//; $_ } @groups;
				}
				#print "Potentials: ".join("\n", @potentials)."\n";
				my $shortest = $potentials[0];
				foreach my $potential (@potentials) {
					$shortest = $potential if length($potential) < length($shortest);
				    }
				#print "Shortest: (".$shortest.")\n";

				# add together all remaining numbers in our string
				# e.g. 0.3 + 0.1 + 0.2 + 0.8 --> 1.4
				my @splittree = split(/([\(\):,])/, $shortest);
				my $sum = 0;
				foreach my $item (@splittree) {
					if ($item =~ /\d.[\d]+/) {	# if it's a number
						$sum += $item;	# branch lengths to scale
						#$sum += 0.1;	# fixed branch lengths
					}
				}
				#print "Distance between ".$speciesListCommon->[$speciesIndex]." and ".$speciesListCommon->[$i]." is: ".$sum."\n";
				$distanceToSpecies[$i] = $sum;
			}
		}
	}
	return @distanceToSpecies;
}

sub orderOrthologues {
	my $species = shift;
	my $AAsInOrthologues_ref = shift;
	my $orthologues_ref = shift;
	my $distArrayOriginal = shift;
	my $distArray = [@{$distArrayOriginal}];
	my $speciesListCommon = shift;
	my $speciesListLatin = shift;

	my $AAsInOrthologues = $$AAsInOrthologues_ref;
	my @AAsInOrthologuesArray = split(/,/, $AAsInOrthologues);
	my $orthologues = $$orthologues_ref;
	my @orthologuesArray = split(/;/, $orthologues);

	#split AAs by each character, or every n characters in the case that we have indels
	#my $numCharsPerOrthologue = length($AAsInOrthologues) / scalar(@orthologuesArray);
	#my @AAsInOrthologuesArray = ($AAsInOrthologues =~ /.{$numCharsPerOrthologue}/g);

	my @sortedAAsInOrthologues = ();
	my @sortedOrthologues = ();

	#let's NOT clean up orthologues, e.g. change X.tropicalis to Xenopus_tropicalis
	#for (my $i = 0; $i < scalar(@orthologuesArray); $i++) {
	#	if ($orthologuesArray[$i] =~ m{[A-Z]\.([a-z]*)}) {
	#		my $indexOfCorrectName = first { lc($speciesListLatin->[$_]) =~ m{$1} } 0..scalar(@$speciesListLatin);
	#		$orthologuesArray[$i] = $speciesListLatin->[$indexOfCorrectName];
	#	}
	#}

	#until we've printed all items, do...
	my $isDone = 0;
	do {
		# find smallest non-negative item in distance array
		my $shortest = 100;
		my $indexOfShortest = 0;
		for(my $i = 0; $i < scalar(@$distArray); $i++) {
			if (($distArray->[$i] != -1) && ($distArray->[$i] < $shortest)) {
				$shortest = $distArray->[$i];
				$indexOfShortest = $i;
			}
		}
		#print "Shortest distance is ".$shortest." (".$speciesListCommon->[$indexOfShortest].")... push!\n";
		if ($shortest != 0) {
			my $closestSpecies = lc($speciesListCommon->[$indexOfShortest]);

			# print corresponding species and its AA
			my $indexInOrthologues = first { lc($orthologuesArray[$_]) eq $closestSpecies } 0..scalar(@orthologuesArray);
			if (defined $indexInOrthologues) {		
				#print "Index in orthologues  = $indexInOrthologues, closest species is $closestSpecies.\n";		

				push(@sortedAAsInOrthologues, $AAsInOrthologuesArray[$indexInOrthologues]);
				push(@sortedOrthologues, $orthologuesArray[$indexInOrthologues]);

				#print "Array is now: ".join("",@sortedAAsInOrthologues)."\n";
				#print "Array is now: ".join("",@sortedOrthologues)."\n";
			}
		}

		# set item in distance array to negative
		$distArray->[$indexOfShortest] = -1;
		#print "Distance array is now: ".join(" ",@$distArray)."\n";

		$isDone = 1;
		for (my $i = 0; $i < scalar(@$distArray); $i++) {
	    		if ($distArray->[$i] != -1) {
				$isDone = 0;
			}
		}
	}
	while (!$isDone);

	$$AAsInOrthologues_ref = join(",", @sortedAAsInOrthologues);
	$$orthologues_ref = join(";", @sortedOrthologues);
}

sub getEnsemblSpecies {
	my $speciesListCommon_ref = shift;
	my $speciesListCommon = $$speciesListCommon_ref;
	my $speciesListLatin_ref = shift;
	my $speciesListLatin = $$speciesListLatin_ref;

	my $html   = get("http://uswest.ensembl.org/info/about/species.html")
		or die "Couldn't get Ensembl species names from http://uswest.ensembl.org/info/about/species.html.";
	$html =~ m{Ensembl [Ss]pecies[</a-z0-9>\s]*<table(.*?)</table>}s || die;
	my $speciestable = $1;
	while ($speciestable =~ m{<a href=.*?>(.*?)</a>[)]?<br />(.*?)</td>}g) {
		if (!($1 =~ m{preview })) { 	# skip species that are only preview
			my $commonname = $1;
			$commonname =~ tr[ ][_]d;
			my $latinname;
			my $temp = $2;
			if ($temp =~ m{<i>(.*?)</i>}) {
				$latinname = $1;
				$latinname =~ tr[ ][_]d;
				if ($commonname eq $latinname) {
					$commonname =~ m{([a-zA-Z]*?)_([a-z]*)};
					$commonname = substr($1, 0, 1).".".$2;
				}
			}
			else {
				$latinname = $commonname;
				$commonname =~ m{([a-zA-Z]*?)_([a-z]*)};
				$commonname = substr($1, 0, 1).".".$2;
			}
			push(@$speciesListCommon, $commonname);
			push(@$speciesListLatin, $latinname);
		}		
	}
	#print "@$speciesListCommon\n@$speciesListLatin\n";
}

