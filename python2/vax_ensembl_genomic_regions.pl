#!/usr/bin/perl -w

=head1 LICENSE

  Copyright (c) 2007-2013 Michael Yourshaw.  All rights reserved.

=head1 CONTACT

  Please email comments or questions to the developer at <myourshaw@ucla.edu>.

=cut

=head1 NAME

Get Genomic Regions - a script to get genomic region intervals from a gene or coordinate list

Version 1.0.β.5

by Michael Yourshaw (myourshaw@ucla.edu)
=cut

#hane disease
#--non_coding_transcripts --db_connection local -v -i "/Users/myourshaw/lab/Ataxia/gene_lists_test/Hane_disease_gene/disease_gene_GR.genes.txt"
#small test
#--db_connection local -v -i "/Users/myourshaw/lab/Ataxia/gene_lists_test/test_KCNA1.txt"
#-five_prime_UTR -three_prime_UTR --db_connection local -v -i "/Users/myourshaw/lab/dev-myourshaw/trunk/perl/get_genomic_regions_EXOSC3.txt"
#-five_prime_UTR -three_prime_UTR --db_connection local -v -i "/Users/myourshaw/lab/dev-myourshaw/trunk/perl/get_genomic_regions_test.txt"
#--include_known_genes --db_connection nelson -i "/home/myourshaw/lab/dev-myourshaw/trunk/perl/get_genomic_regions_EXOSC3.txt"
#--quiet --db_connection nelson -i "/home/myourshaw/lab/dev-myourshaw/trunk/perl/get_genomic_regions_ALLGENES.txt"
#--db_connection local -v -i "/Users/myourshaw/lab/git-myourshaw/perl/get_genomic_regions_EXOSC3.txt"
#--db_connection nelson -v -i "/home/myourshaw/lab/git-myourshaw/perl/get_genomic_regions_EXOSC3.txt" -d /scratch0/tmp/myourshaw/exome_data/diseasegene

#--quiet --db_connection nelson --five_prime_splice_region_bases 6 --three_prime_splice_region_bases 27 -i /home/myourshaw/apps/scripts/perl/get_genomic_regions_ALLGENES.txt -d /scratch0/tmp/myourshaw/exome_data/diseasegene
#--quiet --db_connection nelson --non_coding_transcripts --include_putative_genes --include_novel_genes --include_pseudogenes --five_prime_splice_region_bases 6 --three_prime_splice_region_bases 27 --five_prime_UTR --three_prime_UTR --intron -i /home/myourshaw/apps/scripts/perl/get_genomic_regions_ALLGENES.txt -d /scratch0/tmp/myourshaw/gene_intervals
#--quiet --db_connection nelson --non_coding_transcripts --include_putative_genes --include_novel_genes --include_pseudogenes --five_prime_UTR --three_prime_UTR -i /home/myourshaw/apps/scripts/perl/get_genomic_regions_ALLGENES.txt -d /scratch0/tmp/myourshaw/gene_intervals_test

#--db_connection nelson -v -i "/scratch0/tmp/myourshaw/gene_intervals_test/test_KCNA1.txt" -d /scratch0/tmp/myourshaw/gene_intervals_test

#coding_non_coding_putative_novel_pseudo_5utr_3utr_ess_5ss4_3ss13
#-- verbose --db_connection nelson --non_coding_transcripts --include_putative_genes --include_novel_genes --include_pseudogenes --five_prime_UTR --three_prime_UTR --cis_splice_site --five_prime_splice_region_bases 4 --three_prime_splice_region_bases 13 -i /home/myourshaw/apps/scripts/perl/get_genomic_regions_ALLGENES.txt -d /scratch1/tmp/myourshaw/gene_intervals_test

#default
#--verbose --db_connection nelson -i /home/myourshaw/apps/scripts/perl/get_genomic_regions_ALLGENES.txt -d /scratch1/tmp/myourshaw/gene_intervals_test

#protein_coding
#--verbose --db_connection nelson --include_biotype protein_coding -i /scratch1/tmp/myourshaw/resources/intervals/EnsemblIntervals/ALLGENES.txt -d /scratch1/tmp/myourshaw/resources/intervals/EnsemblIntervals/71/protein_coding ;
#--include_biotype protein_coding -i /scratch1/tmp/myourshaw/resources/intervals/EnsemblIntervals/ALLGENES.txt -d /scratch1/tmp/myourshaw/resources/intervals/EnsemblIntervals72/protein_coding_genes \

#all-genes
#--verbose --db_connection nelson --all_genes -d /scratch1/tmp/myourshaw/gene_intervals_test

#vax installer
#-host cortex.local -port 3306 -user ensembl -password ensembl -output_dir /scratch1/vax/75/ensembl/ensembl_intervals/all_genes --all_genes --include_putative_genes --include_novel_genes --include_pseudogenes --five_prime_UTR --three_prime_UTR --cis_splice_site --five_prime_splice_region_bases 4 --three_prime_splice_region_bases 13
#-host cortex.local -port 3306 -user ensembl -password ensembl -output_dir /scratch1/vax/vax-current/ensembl/ensembl_intervals/protein_coding_genes --all_genes --include_biotype protein_coding

use strict;
use Getopt::Long;
use Cwd;
use File::Spec;
use File::Path;
use FileHandle;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor;
use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(MAX_DISTANCE_FROM_TRANSCRIPT);
use Bio::EnsEMBL::Variation::Utils::Sequence qw(unambiguity_code);
#use Date::Calc qw(Today_and_Now Delta_YMDHMS);
use POSIX qw(ceil floor);
use List::Util qw(first max maxstr min minstr reduce shuffle sum);

my $program_header =<<'HEADER_END';
##---------------------##
## GET GENOMIC REGIONS ##
##---------------------##
##
##version 1.0.β.6
##
##© 2007-2014 Michael Yourshaw (myourshaw\@ucla.edu) all rights reserved
##
HEADER_END

my $start_seconds = time();

# configure from command line opts
my $config = &configure(scalar @ARGV);

#array to hold all genomic regions -- the union of these regions will be the output intervals
my @all_regions;
my @query_regions;
my $line_number = 0;

# run the main sub routine
&main($config);

my $end_seconds = time();

print_log($config,sprintf("program run duration: %u seconds",$end_seconds-$start_seconds));

# this is the main sub-routine - it needs the configured $config hash
sub main {
	my $config = shift;
	
	debug("Starting...");
	
	my $species = $config->{species};
	
	# get adaptors
	my $vfa = $config->{reg}->get_adaptor($species, 'variation', 'variationfeature');
	$config->{tva} = $config->{reg}->get_adaptor($species, 'variation', 'transcriptvariation');
	
	# get fake ones for species with no var DB
	if(!defined($vfa)) {
		$vfa = Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor->new_fake($species);
	}
	else {
		$vfa->db->include_failed_variations($config->{include_failed}) if defined($vfa->db) && $vfa->db->can('include_failed_variations');
	}
	
	$config->{sa} = $config->{reg}->get_adaptor($species, 'core', 'slice');
	$config->{ga} = $config->{reg}->get_adaptor($species, 'core', 'gene');
	$config->{ta} = $config->{reg}->get_adaptor($species, 'core', 'transcript');
	$config->{ea} = $config->{reg}->get_adaptor($species, 'core', 'exon');
	$config->{tra} = $config->{reg}->get_adaptor($species, 'core', 'translation');
	
	# check we got slice adaptor - can't continue without a core DB
	die("ERROR: Could not connect to core database\n") unless defined $config->{sa} and defined $config->{ga} and defined $config->{ta} and defined $config->{ea};
	
	# create a hash to hold slices so we don't get the same one twice
	my %slice_hash = ();
	my @new_vfs;
	my %vf_hash;
	
	my $transcript_cache;
	
	my $in_file_handle = $config->{in_file_handle};
    if ($config->{all_genes}) {
        debug("Considering all genes, regardless of any input queries.");
    }
    else{
        debug(sprintf("Reading queries from %s",$config->{input_file}));
    }
	debug(sprintf("Genomic region details will be written to %s",$config->{details_file}));
	debug(sprintf("Log of warnings and errors will be written to %s",$config->{log_file}));
	debug(sprintf("Genomic target intervals will be written to %s",$config->{intervals_file}));
	debug(sprintf("Bed file will be written to %s",$config->{bed_file}));
	debug(sprintf("Gene list will be written to %s",$config->{gene_file}));
	debug(sprintf("Query sizes will be written to %s",$config->{query_size_file}));
	
	# read the file
	while($config->{all_genes} || (defined($in_file_handle) && <$in_file_handle> ) ) {
        my $query;
        if ($config->{all_genes}) {
			$query = '$ALLGENES$';
            $config->{all_genes} = 0;
        }
        else{
            chomp;
            $line_number++;
            # header or comment line?
            next if /^\#/;
            $query = trim((split /\t/)[$config->{query_column}]);
        }

		#parse query for a recognizible id
		if ($query =~ /\$ALLGENES\$/) {
			my @genes = @{ $config->{ga}->fetch_all() };
			print_log($config,sprintf("considering all %u genes",scalar @genes));
			foreach my $gene(@genes){
				#last unless $gene_limit--;
				push_gene_regions($config,$query,$gene);
			}
		}
		#Chromosome position
		elsif ($query =~ /^\s*(?:chr)?([\S]+):([,\d]+)-([,\d]+)(?:\s+([01+-]{1}))?/i){
			push_chromosome_position($config, $query);
		}
		
		#Ensmbl gene stable id
		elsif ($query =~ /^\s*(ENSG\d{11})/i) {
			my $gene_id = $1;
			my $gene = $config->{ga}->fetch_by_stable_id($gene_id);
			if (defined($gene)){
				push_gene_regions($config,$query,$gene);
			}
			else{
				print_log($config, sprintf("line %u: query <%s> looks like an Ensembl gene ID but isnt; not processed", $line_number,$query));
			}
		}
		
		#Ensmbl transcript stable id
		elsif ($query =~ /^\s*(ENST\d{11})/i) {
			my $transcript_id = $1;
			my $transcript = $config->{ta}->fetch_by_stable_id($transcript_id);
			if (defined($transcript)) {
				push_transcript_regions($config,$query,$transcript);
			}
			else{
				print_log($config, sprintf("line %u: query <%s> looks like an Ensembl transcript ID but isn't; not processed", $line_number,$query));
			}
		}
		
		#Ensmbl translation (protein) stable id
		elsif ($query =~ /^\s*(ENSP\d{11})/i) {
			my $translation_id = $1;
			my $translation =$config->{tra}->fetch_by_stable_id($translation_id);
			if (defined($translation)) {
				my $transcript = $translation->transcript;
				if (defined($transcript)) {
					push_transcript_regions($config,$query,$transcript);
				}
				else{
					print_log($config, sprintf("line %u: query <%s> looks like an Ensembl protein ID but doesn't have a transcript; not processed", $line_number,$query));
				}
			}
			else {
					print_log($config, sprintf("line %u: query <%s> looks like an Ensembl protein ID but isn't; not processed", $line_number,$query));
			}
		}

		#Ensmbl exon stable id
		elsif ($query =~ /^\s*(ENSE\d{11})/i) {
			my $exon_id = $1;
			my $exon = $config->{ea}->fetch_by_stable_id($exon_id);
			#can't choose utr,CDS, splice site because we don't know the transcript
			if (defined($exon)) {
						push_region($config,
							$query,
							"", #$gene->external_name,
							"", #$gene->stable_id,
							"", #$transcript->stable_id,
							"", #$transcript->is_canonical()?"1":"0",
							"", #protein_id
							$exon->stable_id,
							"intergenic",
                            "",
                            "",
							$exon->slice->seq_region_name(),
							$exon->start,
							$exon->end,
							$exon->strand);
			}
			else{
				print_log($config, sprintf("line %u: query <%s> looks like an Ensembl exon ID but isn't; not processed", $line_number,$query));
			}
		}
		
		#GO (gene ontology) accession
		elsif ($query =~ /^\s*(GO:\d+)/i) {
			my $go_id = $1;
			my @genes = @{ $config->{ga}->fetch_all_by_GOTerm_accession($go_id) };
			if (scalar @genes){
				foreach my $gene (@genes) {
					push_gene_regions($config,$query,$gene);
				}
			}
			else{
				print_log($config, sprintf("line %u: query <%s> looks like a GO (gene ontology) accession but isnt; not processed", $line_number,$query));
			}
		}
		
		#Gene external name
		elsif ($query =~ /^\s*([^#\s]+).*#*/) {
			my $gene_id = $1;
			my @genes = @{ $config->{ga}->fetch_all_by_external_name($gene_id) };
			if (@genes) {
				my $gene = $genes[0];
				#if the external name refers to more than one gene, infor the user and pick the most likely
				if (@genes > 1) {
					foreach my $g (@genes) {
						if ($g->external_name eq $query) {
							$gene = $g;
							last;
						}
					}
					my $count_temp = @genes;
					my $external_name = $gene->external_name;
					my $stable_id = $gene->stable_id;
					print_log($config, sprintf("line %u: query <%s> identifies %u genes; using %s (%s). If you meant something else, use the stable ID, ENSG... from the following list:",
							$line_number,
							$query,
							$count_temp,
							$external_name,
							$stable_id));
					foreach my $g (@genes) {
						print_log($config, sprintf("\t%s %s:%u-%u %s (%s)",
							defined($g->external_name) ? $g->external_name : "",
							$gene->slice->seq_region_name,
							$gene->start,
							$gene->end,
							defined($g->description) ? $g->description : "",
							$g->stable_id));
					}
				}
				if ($gene->external_name ne $query) {
					print_log($config, sprintf("line %u: Ensembl's external name for query <%s> is %s (%s:%u-%u %s %s)",
						$line_number,
						$query,
						defined($gene->external_name) ? $gene->external_name : "",
						$gene->slice->seq_region_name,
						$gene->start,
						$gene->end,
						defined($gene->description) ? $gene->description : "",
						$gene->stable_id));
				}
				push_gene_regions($config,$query,$gene);
			} #if (@genes)
			else {
				print_log($config, sprintf("line %u: query <%s> looks like a gene id but isn't recognized by Ensembl", $line_number,$query));
			}
			#get_by_gene_name($1);
		} #elsif (/^\s*([^#\s]+).*#*/) gene external name
		
		#unrecognized
		else {
			print_log($config, sprintf("line %u: query <%s> has an unexpected format; not processed", $line_number,$_));
		}
	} #while(<$in_file_handle>)
	
	#close some files
    close($config->{in_file_handle}) unless !defined($config->{in_file_handle});
	close($config->{details_file_handle});
	close($config->{gene_file_handle});
	
	#create interval files
	output_union_regions($config);
    
    close $config->{intervals_file_handle} or die "can't close $config->{intervals_file_handle}";
    open(FH,">",$config->{done_file}) or die "Can't create $config->{done_file}: $!";
    close(FH);
	
	debug("Finished");

} #sub main


#get all genes (and optionally intergenic regions) from chromosome position, push regions to @all_regions, and print details
sub push_chromosome_position {
	my $config = shift;
	my $query = shift;
	my ($chr,$start,$end,$strand);
	$query =~ /^\s*(?:chr)?([\S]+):([,\d]+)-([,\d]+)(?:\s+([01+-]{1}))?/i;
	if ($1 && $2 && $3) {
		($chr,$start,$end,$strand) = ($1,$2,$3,defined($4) && ($4 eq "+" || $4 == 1) ? 1 : defined($4) && ($4 eq "-" || $4 == -1) ? -1 : 1 );
		$chr = 'MT' if $chr eq 'M';
		$start =~ s/,//g;
		$end =~ s/,//g;
		$start = $start <= $end ? $start : $end;
		$end = $end >= $start ? $end : $start;
	}
	my $slice;
	if(defined($start) && defined($end)){
		$slice = $config->{sa}->fetch_by_region('chromosome', $chr, $start, $end, $strand);
	}
	else{
		$slice = $config->{sa}->fetch_by_region('chromosome', $chr);
	}
	if ($slice) {
		my @genes = @{ $slice->get_all_Genes };
			if($config->{intergenic}){
				my $ig_start = $start;
				my $last = my $ig = 1;
				foreach my $gene (@genes) {
					#get gene not relative to slice
					$gene = $config->{ga}->fetch_by_stable_id($gene->stable_id);
					if ($gene->start > $ig_start) {
						push_region($config,
							$query,
							"", #$gene->external_name,
							"", #$gene->stable_id,
							"", #$transcript->stable_id,
							"", #$transcript->is_canonical()?"1":"0",
							"", #protein_id
							"", #$exon->stable_id,
							"intergenic",
                            "",
                            "",
							$chr,
							$ig_start,
							$gene->start - 1,
							$strand);
						$ig_start = $gene->end + 1;
						$ig++;
					}
					push_gene_regions($config,$query,$gene);
				}
				#intergenic sequence after last gene, if any
				if ($ig_start <= $end) {
						push_region($config,
							$query,
							"", #$gene->external_name,
							"", #$gene->stable_id,
							"", #$transcript->stable_id,
							"", #$transcript->is_canonical()?"1":"0",
							"", #protein_id
							"", #$exon->stable_id,
							"intergenic",
                            "",
                            "",
							$chr,
							$ig_start,
							$end,
							$strand);
				}
			}
			else {
				foreach my $gene (@genes) {
					#get gene not relative to slice
					my $gene = $config->{ga}->fetch_by_stable_id($gene->stable_id);
					push_gene_regions($config,$query,$gene);
				}
			}
		} #if ($slice)
	else {
		print_log($config, sprintf("Query <%s> isn't a region; not processed", $query));
	}
} #sub push_chromosome_position

#push gene regions to @all_regions and print details
sub push_gene_regions {
	my $config = shift;
	my $query = shift;
	my $gene = shift;
	if(
		 (keys (%{$config->{include_biotype}}) == 0 || exists(${$config->{include_biotype}}{$gene->biotype}))
		 && (keys (%{$config->{exclude_biotype}}) == 0 || !exists(${$config->{exclude_biotype}}{$gene->biotype}))
		 && (($config->{include_known_genes} && uc($gene->status) eq 'KNOWN')
		 || ($config->{include_putative_genes} && uc($gene->status) eq 'PUTATIVE')
		 || ($config->{include_novel_genes} && uc($gene->status) eq 'NOVEL'))
		 && ($gene->biotype !~ /pseudogene/i || ($gene->biotype =~ /pseudogene/i && $config->{include_pseudogenes}))
		) {
		print_gene_list($config, sprintf("%s\t%s\t%s\t%s\t%u\t%u\t%s\t%s\t%s\t%s",
			$query,
			$gene->external_name,
			$gene->stable_id,
			$gene->slice->seq_region_name,
			$gene->start,
			$gene->end,
			$gene->strand == 1 ? "+" : "-",
			$gene->biotype,
			$gene->status,
			defined($gene->description) ? $gene->description : "",
			));
		my @transcripts = @{$gene->get_all_Transcripts};
		if (@transcripts) {
			foreach my $transcript (@transcripts) {
				push_transcript_regions($config,$query,$transcript);
			} #foreach my $transcript (@transcripts)
		} #if (@transcripts)
	}
} #sub push_gene_regions

#push transcript regions to @all_regions and print details
sub push_transcript_regions {
	my $config = shift;
	my $query = shift;
	my $transcript = shift;
	my $gene = $config->{ga}->fetch_by_transcript_stable_id($transcript->stable_id);
	my $translation =$config->{tra}->fetch_by_Transcript($transcript);
	my $protein_id = defined($translation) ? $translation->{stable_id} : "" ;

	print_log($config, sprintf("line %u: query <%s>, GeneSymbol %s, GeneID %s, TranscriptID %s has no coding region",
		$line_number,
		$query,
		$gene->external_name,
		$gene->stable_id,
		$transcript->stable_id)) unless ($gene->biotype ne 'protein_coding' || (defined($transcript->coding_region_start) && defined($transcript->coding_region_end)));
	my @exons = @{ $transcript->get_all_Exons() };
	my @introns = @{ $transcript->get_all_Introns() };
	my $exon_count = @exons;
	for(my $i = 0; $i <= $#exons; $i++){
		my $exon = $exons[$i];
		push_exon_regions($config,$query,$gene,$transcript,$protein_id,$exon,$i,$#exons,$i <= $#introns ? $introns[$i] : undef);
	} #for(my $i = 0; $i <= $#exons; $i++)
}

#push exon/intron regions to @all_regions and print details
sub push_exon_regions {
	my $config = shift;
	my $query = shift;
	my $gene = shift;
	my $transcript = shift;
	my $protein_id = shift;
	my $exon = shift;
	my $exon_index = shift;
	my $exon_max_index = shift;
	my $intron = shift;
	#my $gene = $config->{ga}->fetch_by_exon_stable_id($exon->stable_id);
	my $coding_start = $transcript->coding_region_start;
	my $coding_end = $transcript->coding_region_end;
	my $transcript_has_coding = (defined($coding_start) && defined($coding_end));
	if ($config->{non_coding_transcripts} || $transcript_has_coding || $gene->biotype ne 'protein_coding') {
		#first splice region
		if (($config->{non_coding_transcripts} && $gene->biotype ne 'protein_coding' && !defined($coding_start) || !defined($coding_end))
				|| ((($config->{five_prime_UTR} && $exon->strand == 1) || ($config->{three_prime_UTR} && $exon->strand == -1)) && $exon->start < $coding_start)
				|| ($config->{CDS} && $exon->start <= $coding_end && $exon->end >= $coding_start && $exon->start >= $coding_start)
				|| ((($config->{three_prime_UTR} && $exon->strand == 1) || ($config->{five_prime_UTR} && $exon->strand == -1)) && $exon->end > $coding_end && $exon->start > $coding_end)){
			#if(!defined($config->{intron}) || !$config->{intron}){
				if ($exon->strand == 1 && $exon_index != 0 && defined($config->{three_prime_splice_region_bases}) && $config->{three_prime_splice_region_bases} > 0){
					push_region($config,
						$query,
						$gene->external_name,
						$gene->stable_id,
						$transcript->stable_id,
						$transcript->is_canonical()?"1":"0",
						$protein_id,
						$exon->stable_id,
						"intron",
                        $exon_index+1,
                        "three_prime_splice_region",
						$exon->seq_region_name,
						$exon->start-2-$config->{three_prime_splice_region_bases},
						$exon->start-3,
						$exon->strand);
				}
				if ($exon->strand == 1 && $exon_index != 0 && defined($config->{cis_splice_site}) && $config->{cis_splice_site}){
					push_region($config,
						$query,
						$gene->external_name,
						$gene->stable_id,
						$transcript->stable_id,
						$transcript->is_canonical()?"1":"0",
						$protein_id,
						$exon->stable_id,
						"intron",
                        $exon_index+1,
                        "three_prime_cis_splice_site",
						$exon->seq_region_name,
						$exon->start-2,
						$exon->start-1,
						$exon->strand);
				}
				if ($exon->strand == -1 && $exon_index != 0 && defined($config->{five_prime_splice_region_bases}) && $config->{five_prime_splice_region_bases} > 0){
					push_region($config,
						$query,
						$gene->external_name,
						$gene->stable_id,
						$transcript->stable_id,
						$transcript->is_canonical()?"1":"0",
						$protein_id,
						$exon->stable_id,
						"intron",
                        $exon_index+1,
                        "five_prime_splice_region",
						$exon->seq_region_name,
						$exon->start-2-$config->{five_prime_splice_region_bases},
						$exon->start-3,
						$exon->strand);
				}
				if ($exon->strand == -1 && $exon_index != 0 && defined($config->{cis_splice_site}) && $config->{cis_splice_site}){
					push_region($config,
						$query,
						$gene->external_name,
						$gene->stable_id,
						$transcript->stable_id,
						$transcript->is_canonical()?"1":"0",
						$protein_id,
						$exon->stable_id,
						"intron",
                        $exon_index+1,
                        "five_prime_cis_splice_site",
						$exon->seq_region_name,
						$exon->start-2,
						$exon->start-1,
						$exon->strand);
				}
			#}
		}
		#non-coding
		if (($config->{non_coding_transcripts} && !defined($coding_start) || !defined($coding_end))){
			push_region($config,
				$query,
				$gene->external_name,
				$gene->stable_id,
				$transcript->stable_id,
				$transcript->is_canonical()?"1":"0",
				$protein_id,
				$exon->stable_id,
				"exon",
                $exon_index+1,
                "noncoding",
				$exon->seq_region_name,
				$exon->start,
				$exon->end,
				$exon->strand);
		}
		#coding
		else {
			#5'utr
			if ($config->{five_prime_UTR} && (($exon->strand == 1 && $exon->start < $coding_start) || ($exon->strand == -1 && $exon->end > $coding_end))){
				my $utr_start = $exon->strand == 1 ? $exon->start : max($exon->start,$coding_end+1);
				my $utr_end = $exon->strand == 1 ? min($coding_start-1,$exon->end) : $exon->end;
				push_region($config,
					$query,
					$gene->external_name,
					$gene->stable_id,
					$transcript->stable_id,
					$transcript->is_canonical()?"1":"0",
					$protein_id,
					$exon->stable_id,
					"exon",
                    $exon_index+1,
                    "5'utr",
					$exon->seq_region_name,
					$utr_start,
					$utr_end,
					$exon->strand);
			}
			#CDS
			if ($config->{CDS} && $exon->start <= $coding_end && $exon->end >= $coding_start){
				my $cds_start = max($exon->start,$coding_start);
				my $cds_end = min($exon->end,$coding_end);
				push_region($config,
					$query,
					$gene->external_name,
					$gene->stable_id,
					$transcript->stable_id,
					$transcript->is_canonical()?"1":"0",
					$protein_id,
					$exon->stable_id,
					"exon",
                    $exon_index+1,
                    "CDS",
					$exon->seq_region_name,
					$cds_start,
					$cds_end,
					$exon->strand);
			}
			#3'utr
			if ($config->{three_prime_UTR} && (($exon->strand == 1 && $exon->end > $coding_end) || ($exon->strand == -1 && $exon->start < $coding_start))){
				my $utr_start = $exon->strand == 1 ? max($exon->start,$coding_end+1) : $exon->start;
				my $utr_end = $exon->strand == 1 ? $exon->end : min($coding_start-1,$exon->end);
				push_region($config,
					$query,
					$gene->external_name,
					$gene->stable_id,
					$transcript->stable_id,
					$transcript->is_canonical()?"1":"0",
					$protein_id,
					$exon->stable_id,
					"exon",
                    $exon_index+1,
                    "three_prime_UTR",
					$exon->seq_region_name,
					$utr_start,
					$utr_end,
					$exon->strand);
			}
		}
		#second splice region
		if (($config->{non_coding_transcripts} && !defined($coding_start) || !defined($coding_end))
				|| ((($config->{five_prime_UTR} && $exon->strand == 1) || ($config->{three_prime_UTR} && $exon->strand == -1)) && $exon->start < $coding_start && $exon->end < $coding_start)
				|| ($config->{CDS} && $exon->start <= $coding_end && $exon->end >= $coding_start && $exon->end <= $coding_end)
				|| ((($config->{three_prime_UTR} && $exon->strand == 1) || ($config->{five_prime_UTR} && $exon->strand == -1)) && $exon->end > $coding_end)){
			#if(!defined($config->{intron}) || !$config->{intron}){
				if ($exon->strand == 1 && $exon_index != $exon_max_index && defined($config->{cis_splice_site}) && $config->{cis_splice_site}){
					push_region($config,
						$query,
						$gene->external_name,
						$gene->stable_id,
						$transcript->stable_id,
						$transcript->is_canonical()?"1":"0",
						$protein_id,
						$exon->stable_id,
						"intron",
                        $exon_index+1,
                        "five_prime_cis_splice_site",
						$exon->seq_region_name,
						$exon->end+1,
						$exon->end+2,
						$exon->strand);
                }
				if ($exon->strand == 1 && $exon_index != $exon_max_index && defined($config->{five_prime_splice_region_bases}) && $config->{five_prime_splice_region_bases} > 0){
					push_region($config,
						$query,
						$gene->external_name,
						$gene->stable_id,
						$transcript->stable_id,
						$transcript->is_canonical()?"1":"0",
						$protein_id,
						$exon->stable_id,
						"intron",
                        $exon_index+1,
                        "five_prime_splice_region",
						$exon->seq_region_name,
						$exon->end+3,
						$exon->end+2+$config->{five_prime_splice_region_bases},
						$exon->strand);
				}
				if ($exon->strand == -1 && $exon_index != $exon_max_index && defined($config->{cis_splice_site}) && $config->{cis_splice_site}){
					push_region($config,
						$query,
						$gene->external_name,
						$gene->stable_id,
						$transcript->stable_id,
						$transcript->is_canonical()?"1":"0",
						$protein_id,
						$exon->stable_id,
						"intron",
                        $exon_index+1,
                        "three_prime_cis_splice_site",
						$exon->seq_region_name,
						$exon->end+1,
						$exon->end+2,
						$exon->strand);
				}
				if ($exon->strand == -1 && $exon_index != $exon_max_index && defined($config->{three_prime_splice_region_bases}) && $config->{three_prime_splice_region_bases} > 0){
					push_region($config,
						$query,
						$gene->external_name,
						$gene->stable_id,
						$transcript->stable_id,
						$transcript->is_canonical()?"1":"0",
						$protein_id,
						$exon->stable_id,
						"intron",
                        $exon_index+1,
                        "three_prime_splice_region",
						$exon->seq_region_name,
						$exon->end+3,
						$exon->end+2+$config->{three_prime_splice_region_bases},
						$exon->strand);
				}
			#}
		}
		#intron
		if (defined($intron)) {
			#entire intron
			if (defined($config->{intron}) && $config->{intron}) {
				push_region($config,
					$query,
					$gene->external_name,
					$gene->stable_id,
					$transcript->stable_id,
					$transcript->is_canonical()?"1":"0",
					$protein_id,
					$exon->stable_id,
					"intron",
                    $exon_index+1,
                    "",
					$intron->seq_region_name,
					$intron->start,
					$intron->end,
					$exon->strand);
			}
		} #if (defined($intron))
	} #if ($config->{non_coding_transcripts} || $transcript_has_coding || $gene->biotype ne 'protein_coding')
}

#push genomic region to @all_regions and print to genomic region details file 
sub push_region {
	my $config = shift;
	my $file_handle = $config->{details_file_handle};
	my ($query,$geneSymbol,$gene,$transcript,$canonical,$protein_id,$exon,$region,$region_number,$sub_region,$chr,$start,$end,$strand) = @_;
	push @all_regions, {zero_chr=>chr2zerochr($chr), start=>$start, end=>$end, strand=>$strand};
	push @query_regions, {query=>$query, zero_chr=>chr2zerochr($chr), start=>$start, end=>$end, strand=>$strand};
	my $data = sprintf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%u\t%u\t%s\n",
		$query,
		$geneSymbol,
		$gene,
		$transcript,
        $canonical,
		$protein_id,
		$exon,
		$region,
        $region_number,
        $sub_region,
		$chr,
		$start,
		$end,
		$strand == 1 ? "+" : "-");
	$data .= "\n" if substr($data, -1) ne "\n";
	print $file_handle $data;
	print $data if (defined $config->{verbose} && ! defined $config->{quiet});
}

#get union of genomic regions and print to target intervals file
sub output_union_regions {
	my $config = shift;
	my $total_base_count;
	if($#all_regions){
		my $base_count = 0;
		debug("Writing genomic target intervals to $config->{intervals_file} ...");
		my @sorted_regions = sort {$a->{zero_chr} cmp $b->{zero_chr} || $a->{start} <=> $b->{start} || $a->{end} <=> $b->{end}} @all_regions;
		my $that_chr;
		my $this_chr;
		my $that_start;
		my $this_start;
		my $that_end = -1;
		my $this_end;
		foreach my $region(@sorted_regions){
			$this_chr = $region->{zero_chr};
			$this_start = $region->{start};
			$this_end = $region->{end};
			if (defined($that_chr) && ($this_chr ne $that_chr || $this_start > $that_end+1)){
				#this region is non-overlapping or non-contiguous with that region;
				my $name = get_gene_names(zerochr2chr($that_chr), $that_start, $that_end) || 'no_genes';
				print_intervals($config, zerochr2chr($that_chr), $that_start, $that_end, '+', $name);
				$base_count += $that_end-$that_start+1;
				$that_chr = $this_chr;
				$that_start = $this_start;
				$that_end = $this_end;
			}
			else{
				$that_start = $this_start if ! defined $that_start;
				$that_end = max($that_end,$this_end);
				$that_chr = $this_chr;
			}
		}
		#last region
		my $name = get_gene_names(zerochr2chr($that_chr), $that_start, $that_end) || 'no_genes';
		print_intervals($config, zerochr2chr($that_chr), $that_start, $that_end, '+', $name);
		$base_count += $that_end-$that_start+1;
		$total_base_count = $base_count;
		print_log($config, "Total number of bases in union of all intervals: $total_base_count");
		undef @sorted_regions;
	}
	if($#query_regions){
		my $base_count = 0;
		debug("Writing query sizes to $config->{query_size_file} ...");
		my @sorted_regions = sort {$a->{query} cmp $b->{query} || $a->{zero_chr} cmp $b->{zero_chr} || $a->{start} <=> $b->{start} || $a->{end} <=> $b->{end}} @query_regions;
		my $this_query;
		my $that_query;
		my $that_chr;
		my $this_chr;
		my $that_start;
		my $this_start ;
		my $that_end = -1;
		my $this_end;
		foreach my $region(@sorted_regions){
			$this_query = $region->{query};
			$this_chr = $region->{zero_chr};
			$this_start = $region->{start};
			$this_end = $region->{end};
			if(defined($that_query) && $this_query ne $that_query){
				$base_count += $that_end-$that_start+1;
				print_query_size($config, sprintf("%s\t%u",$that_query,$base_count));
				$base_count = 0;
				$that_query = $this_query;
				$that_chr = $this_chr;
				$that_start = $this_start;
				$that_end = $this_end;
			}
			elsif (defined($that_chr) && ($this_chr ne $that_chr || $this_start > $that_end+1)){
				#this region is non-overlapping or non-contiguous with that region;
				$base_count += $that_end-$that_start+1;
				$that_query = $this_query;
				$that_chr = $this_chr;
				$that_start = $this_start;
				$that_end = $this_end;
			}
			else{
				$that_start = $this_start if ! defined $that_start;
				$that_end = max($that_end,$this_end);
				$that_chr = $this_chr;
			}
		}
		#last region
		$base_count += $that_end-$that_start+1;
		print_query_size($config, sprintf("%s\t%u",$that_query,$base_count));
		print_query_size($config, sprintf("%s\t%u","## union_all_intervals",$total_base_count));
	}
} #sub output_union_regions

sub get_gene_names {
	my $chr = shift;
	my $start = shift;
	my $end = shift;
	my @gene_names;
	my $name = '';
	my $slice;
	my @genes;
	$slice = $config->{sa}->fetch_by_region('chromosome', $chr, $start, $end);
	if ($slice) {
		my @genes = @{ $slice->get_all_Genes };
		foreach my $gene (@genes) {
			push @gene_names, $gene->external_name;
		}
	}
	if (@gene_names){
		$name = join ',', @{get_unique(\@gene_names)};
	}
	#print STDERR sprintf("chr: %s start: %s end: %s gene_name(s): %s\n",$chr,$start,$end,$name);
	return $name;
}

sub get_unique{
    my $list = shift;
    my %seen = ();
    my @uniq = grep { !$seen{$_}++ } @{$list};
    return \@uniq;
}

#print to genomic target intervals file and interval list file
sub print_intervals {
	my $config = shift;
	my $chr = shift;
	my $start = shift;
	my $end = shift;
	my $strand = shift;
	my $name = shift;
	my $file_handle = $config->{intervals_file_handle};
	print $file_handle sprintf("%s:%u-%u\n",$chr,$start,$end);
	$file_handle = $config->{interval_list_file_handle};
	print $file_handle sprintf("%s\t%u\t%u\t%s\t%s\n",$chr,$start,$end,$strand,$name);
	$file_handle = $config->{bed_file_handle};
	print $file_handle sprintf("%s\t%u\t%u\t%s\t1000\t%s\n",$chr,$start-1,$end,$name,$strand);
}

#print to gene list file
sub print_gene_list {
	my $config = shift;
	my $file_handle = $config->{gene_file_handle};
	my $data = shift;
	$data .= "\n" if substr($data, -1) ne "\n";
	print $file_handle $data;
}

#print to query size file
sub print_query_size {
	my $config = shift;
	my $file_handle = $config->{query_size_file_handle};
	my $data = shift;
	$data .= "\n" if substr($data, -1) ne "\n";
	print $file_handle $data;
}

#print to log file
sub print_log {
	my $config = shift;
	my $file_handle = $config->{log_file_handle};
	my $data = shift;
	$data .= "\n" if substr($data, -1) ne "\n";
	print $file_handle $data;
	print $data if (defined $config->{verbose});
}

# prints debug output with time to log, and to STDOUT unless quiet
sub debug {
	my $text = (@_ ? (join "", @_) : "No message");
	my $time = get_time();
	my $file_handle = $config->{log_file_handle};
	print $file_handle $time." - ".$text.($text =~ /\n$/ ? "" : "\n");
	print STDOUT $time." - ".$text.($text =~ /\n$/ ? "" : "\n") unless defined $config->{quiet};
}

##housekeeping subs
sub configure {
	my $args = shift;
	
	my $config = {};
	
	GetOptions(
		$config,
		'help|?',
		
		# input options,
		'config=s',
        'all_genes!',
		'input_file|i=s',
		'format=s',
		
		# DB options
		'species|s=s',
		'registry|r=s',
		'host=s',
		'user|u=s',
		'port=i',
		'password=s',
		'db_version=i',
		'db_connection=s',
		
		# runtime options
		'tmp_dir=s',
		
		#data options
        'cis_splice_site!',
		'five_prime_splice_region_bases=i',
		'three_prime_splice_region_bases=i',
		'upstream_regulatory_bases=i',
		'downstream_regulatory_bases=i',
		'intergenic!',
		'five_prime_UTR!',
		'CDS!',
		'three_prime_UTR!',
		'non_coding_transcripts!',
		'intron!',
		'exclude_biotype=s@',
		'include_biotype=s@',
		'include_known_genes!',
		'include_putative_genes!',
		'include_novel_genes!',
		'include_pseudogenes!',
		'query_column=s',
		
		# output options
		'output_dir|d=s',
		'intervals_file|o=s',
		'interval_list_file=s',
		'details_file=s',
		'log_file|l=s',
		'gene_file|g=s',
		'query_size_file|z',
		'verbose|v',
		'quiet',
	);
	
	# print usage message if requested or no args supplied
	if(defined($config->{help}) || !$args) {
		&usage;
		exit(0);
	}
	
	# config file?
	if(defined $config->{config}) {
		
		open CONFIG, $config->{config} or die "ERROR: Could not open config file \"".$config->{config}."\"\n";
		
		while(<CONFIG>) {
			next if /^\#/;
			my ($key, $value) = split /\s+|\=/;
			$key =~ s/^\-//g;
			$config->{$key} = $value unless defined $config->{$key};
		}
		
		close CONFIG;
	}
	
	# check file format
	if(defined $config->{input_format}) {
		die "ERROR: Unrecognised input format specified \"".$config->{input_format}."\"\n" unless $config->{input_format} =~ /pileup|vcf|guess/i;
	}
	
	# summarise options
	my $header = $program_header.<<'INTRO';
Configuration options:

INTRO
    print $header unless defined $config->{quiet};
    
    my $max_length = (sort {$a <=> $b} map {length($_)} keys %$config)[-1];
    
    foreach my $key(sort keys %$config) {
        print $key.(' ' x (($max_length - length($key)) + 4)).(ref($config->{$key}) ne 'ARRAY' ? $config->{$key} : join(',',$config->{$key}))."\n" unless defined $config->{quiet};
    }
    
    print "\n".("-" x 20)."\n\n" unless defined $config->{quiet};
		
	# set defaults
	my ($v,$d,$f);
    if (!defined($config->{input_file})) {
        ($v,$d,$f) = File::Spec->splitpath(File::Spec->join(cwd(),'ensembl_genomic_regions.txt'));
    }
    else{
        ($v,$d,$f) = File::Spec->splitpath($config->{input_file});
    }
    $config->{all_genes} ||= 0;
	$config->{db_connection} ||= 'nelson';
	$config->{output_dir} ||= $d;
	$d = $config->{output_dir};
	File::Path::mkpath($config->{output_dir});
	$config->{tmpdir} ||= '/tmp';
	$config->{format} ||= 'guess';
	$config->{query_column} ||= 0;
	$config->{CDS} ||= 1;
	$config->{include_known_genes} ||= 1;
	$config->{include_putative_genes} ||= 0;
	$config->{include_novel_genes} ||= 0;
	$config->{include_pseudogenes} ||= 0;
	$config->{cis_splice_site} ||= 1;
	$config->{five_prime_splice_region_bases} ||= 4;
	$config->{three_prime_splice_region_bases} ||= 13;
	$config->{upstream_regulatory_bases} ||= 0;
	$config->{downstream_regulatory_bases} ||= 0;
    my @file_name_info;
    #coding_non_coding_putative_novel_pseudo_5utr_3utr_ess_5ss4_3ss13
    if (defined($config->{include_biotype})) {push @file_name_info, join('_',@{$config->{include_biotype}})};
    if ($config->{include_known_genes} == 1) {push @file_name_info, "known"};
    if ($config->{include_putative_genes} == 1) {push @file_name_info, "putative"};
    if ($config->{include_novel_genes} == 1) {push @file_name_info, "novel"};
    if ($config->{include_pseudogenes} == 1) {push @file_name_info, "pseudo"};
    if ($config->{non_coding_transcripts}) {push @file_name_info, "non-coding"};
    if ($config->{CDS} == 1) {push @file_name_info, "CDS"};
    if ($config->{cis_splice_site} == 1) {push @file_name_info, "ess"};
    if ($config->{five_prime_splice_region_bases} > 0) {push @file_name_info, "5sr$config->{five_prime_splice_region_bases}"};
    if ($config->{three_prime_splice_region_bases} > 0) {push @file_name_info, "3sr$config->{three_prime_splice_region_bases}"};
    if ($config->{five_prime_UTR}) {push @file_name_info, "5utr"};
    if ($config->{three_prime_UTR}) {push @file_name_info, "3utr"};
    if ($config->{upstream_regulatory_bases} > 0) {push @file_name_info, "up$config->{upstream_regulatory_bases}"};
    if ($config->{downstream_regulatory_bases} > 0) {push @file_name_info, "down$config->{downstream_regulatory_bases}"};
    if ($config->{intergenic}) {push @file_name_info, "intergenic"};
    if ($config->{intron}) {push @file_name_info, "intron"};

	if(defined($config->{include_biotype})){
		my %biotypes = map { $_ => $_ } @{$config->{include_biotype}};
		$config->{include_biotype} = \%biotypes;
	}
	if(defined($config->{exclude_biotype})){
		my %biotypes = map { $_ => $_ } @{$config->{exclude_biotype}};
		$config->{exclude_biotype} = \%biotypes;
	}
	
	# connection settings for local Ensembl installation
	if($config->{db_connection} eq 'local'){
		$config->{species} ||= "homo_sapiens";
		$config->{host} ||= '127.0.0.1';
		$config->{port} ||= 3306;
		$config->{password} ||= 'ensembl';
		$config->{user} ||= 'ensembl';
	}
	# connection settings for Nelson lab Ensembl installation
	elsif($config->{db_connection} eq 'nelson'){
			$config->{species} ||= "homo_sapiens";
			$config->{host} ||= 'cortex.local';
			$config->{port} ||= 3306;
			$config->{password} ||= 'ensembl';
			$config->{user} ||= 'ensembl';
	}
	# connection settings for Ensembl Genomes
	elsif($config->{db_connection} eq 'genomes') {
		$config->{host} ||= 'mysql.ebi.ac.uk';
		$config->{port} ||= 4157;
		$config->{user} ||= 'anonymous';
	}
	# connection settings for main Ensembl
	else {
		$config->{species} ||= "homo_sapiens";
		$config->{host}    ||= 'ensembldb.ensembl.org';
		$config->{port}    ||= 5306;
		$config->{user}    ||= 'anonymous';
	}

    #options for file headers
    $config->{config_options_header} = join "\n", map {"## ".$_."=".$config->{$_}} sort keys %$config;
	
	# connect to databases
	$config->{reg} = &connect_to_dbs($config);
	
    #info for file headers
    $config->{config_options_header} .= "\n## Output produced at ".&get_time."\n";
	my $core_mca = $config->{reg}->get_adaptor($config->{species}, 'core', 'metacontainer');
    my $ensembl_api_version = $config->{reg}->software_version;
    my $ensembl_db_version = (defined $core_mca && $core_mca->get_schema_version ? $core_mca->get_schema_version : '?');
	$config->{config_options_header} .= "## Connected to ".$core_mca->dbc->dbname." on ".$core_mca->dbc->host."\n" if defined $core_mca;
	$config->{config_options_header} .=
		"## Using API version ".$ensembl_api_version.
		", DB version ".$ensembl_db_version."\n";

    #file names
    my $file_name_base = 'ensembl_'.$ensembl_db_version.'_'.(join '_',@file_name_info);
	$config->{intervals_file} ||= File::Spec->join($v,$d,$file_name_base.".genomic_target_intervals.txt");
	$config->{bed_file} ||= File::Spec->join($v,$d,$file_name_base.".bed");
	$config->{interval_list_file} ||= File::Spec->join($v,$d,$file_name_base.".interval_list");
	$config->{details_file} ||= File::Spec->join($v,$d,$file_name_base.".genomic_region_details.txt");
	$config->{gene_file} ||= File::Spec->join($v,$d,$file_name_base.".genes.txt");
	$config->{query_size_file} ||= File::Spec->join($v,$d,$file_name_base.".query_sizes.txt");
	$config->{log_file} ||= File::Spec->join($v,$d,$file_name_base.".log.txt");
    
    my ($vd,$dd,$fd) = File::Spec->splitpath($config->{intervals_file});
    $config->{done_file} = File::Spec->join($vd, $dd, '.'.$fd.'.done');
    if (-e $config->{done_file}) {
        print sprintf('Ensembl genomic_regions already done. To reschedule "rm %s"', $config->{done_file});
        exit;
    }

	# get input file handle
	$config->{in_file_handle} = &get_in_file_handle($config) unless $config->{all_genes};
	
	# configure intervals file
	$config->{intervals_file_handle} = &get_intervals_file_handle($config);
	
	# configure bed file
	$config->{bed_file_handle} = &get_bed_file_handle($config);
	
	# configure interval list file
	$config->{interval_list_file_handle} = &get_interval_list_file_handle($config);
	
	# configure details file
	$config->{details_file_handle} = &get_details_file_handle($config);

	# configure gene list file
	$config->{gene_file_handle} = &get_gene_file_handle($config);

	# configure query size file
	$config->{query_size_file_handle} = &get_query_size_file_handle($config);

	# configure log file
	$config->{log_file_handle} = &get_log_file_handle($config);
	
	return $config;
} #sub configure

sub usage {
	my $usage = $program_header.<<'END';
Usage:
perl get_genomic_target_intervals.pl [arguments]

Options
--help | ?             Display this message and quit
--verbose | v          Display verbose output as the script runs [default: off]
--quiet                Suppress status and warning messages [default: off]

--config               Load configuration from file. Any command line options
                       specified overwrite those in the file [default: off]
                       
--all_genes            Return all genes that pass other filters. Overrides queries in input_file. [default: off]

-i | --input_file      Input file - if input_file not specified or contains "$ALLGENES$", returns all genes that pass other filters.
                            Files may contain gene symbols, genomic regions (e.g., 9:36993071-38456365),
                            Ensembl gene, transcript, protein, or exon IDs, or GO (gene ontology) accession codes.
--query_column         Column (first = 0) of tab-separated input file that has the query (gene name, position, ...) [default: 0]
--format               Alternative input file format - unsupported
-o | --intervals_file  Output file for union of all genomic intervals [default: "<input_file | cwd>.genomic_target_intervals.txt"]
--interval_list_file  Output file for union of all genomic intervals in standard format used by picard [default: "<input_file | cwd>.interval_list"]
--details_file   			Detailed output file of each genomic region for each query [default: "<input_file | cwd>.genomic_region_details.txt"]
-g | --gene_file       List of genes [default: "<input_file | cwd>.genes.txt"]
-z | --query_size_file List of number of bases in each query [default: "<input_file | cwd>.query_sizes.txt"]
-l | --log_file        Output log file of warnings and errors [default: "<input_file | cwd>.stderr.txt"]

-s | --species         Species to use [default: "human"]
--host                 Manually define database host [default: "ensembldb.ensembl.org"]
-u | --user            Database username [default: "anonymous"]
--port                 Database port [default: 5306]
--db_connection        Sets DB connection params for
                         ensembl (main Ensembl, default)
                         genomes (Ensembl Genomes)
                         nelson (Nelson lab)
                         local (local computer)
--password             Database password [default: no password]
--genomes              Sets DB connection params for Ensembl Genomes [default: off]
--local                Sets DB connection params for local Ensembl installation [default: on]
--nelson               Sets DB connection params for Nelson lab Ensembl installation [default: off]
-r | --registry_file   Registry file to use defines DB connections [default: off]
                       Defining a registry file overrides above connection settings.
--db_version=[number]  Force script to load DBs from a specific Ensembl version. Not
                       advised due to likely incompatibilities between API and DB
--cis_splice_site           Include intronic 2 bp region bordering exon [default: on]
--five_prime_splice_region_bases        Number of bases in intron adjacent to the five_prime_cis_splice_site [default: 4]
--three_prime_splice_region_bases     Number of bases in intron adjacent to the three_prime_cis_splice_site [default: 13]
--upstream_regulatory_bases Number of bases upstream of gene to include [default: 0]
--downstream_regulatory_bases Number of bases downstream to include [default: 0]
--intergenic              Include intergenic regions of genomic regions [default: off]
--five_prime_UTR                      Include 5-prime untranslated region of exons [default: off]
--CDS                       Include CDS translated region of exons [default: on]
--three_prime_UTR                      Include 3-prime untranslated region of exons [default: off]
--non_coding_transcripts    Include non-coding transcripts of genes that are biotype protein_coding [default: off]
--intron                    Include introns [default: off]
--exclude_biotype           Exclude genes having this biotype (can be used more than once) [default: none]
--include_biotype           Include only genes having this biotype  (can be used more than once) [default: all]
--include_known_genes       Include genes with KNOWN status [default: on]
--include_putative_genes    Include genes with PUTATIVE status [default: off]
--include_novel_genes       Include genes with NOVEL status [default: off]
--include_pseudogenes       Include genes with biotype =~ /pseudogene/ [default: off]

END

	print $usage;
}


sub connect_to_dbs {
	my $config = shift;
	
	# get registry
	my $reg = 'Bio::EnsEMBL::Registry';
	
	# load DB options from registry file if given
	if(defined($config->{registry})) {
		print("Loading DB config from registry file ", $config->{registry});
		$reg->load_all($config->{registry});
	}
	
	# otherwise manually connect to DB server
	else {
		$reg->load_registry_from_db(
			-host       => $config->{host},
			-user       => $config->{user},
			-pass       => $config->{password},
			-port       => $config->{port},
			-db_version => $config->{db_version},
			-species    => $config->{species} =~ /^[a-z]+\_[a-z]+/i ? $config->{species} : undef,
			-verbose    => $config->{verbose} && ! $config->{quiet},
		);
	}
	
	#$reg->set_disconnect_when_inactive();
	
	if($config->{verbose} && ! $config->{quiet}) {
		# get a meta container adaptors to check version
		my $core_mca = $reg->get_adaptor($config->{species}, 'core', 'metacontainer');
		my $var_mca = $reg->get_adaptor($config->{species}, 'variation', 'metacontainer');
		
		if($core_mca && $var_mca) {
			print(
				"Connected to core version ", $core_mca->get_schema_version, " database ",
				"and variation version ", $var_mca->get_schema_version, " database. \n"
			);
		}
	}
	
	return $reg;
}

sub get_in_file_handle {
	my $config = shift;

	# define the filehandle to read input from
	my $in_file_handle = new FileHandle;
	
	if(defined($config->{input_file})) {
		
		# check defined input file exists
		die("ERROR: Could not find input file ", $config->{input_file}, "\n") unless -e $config->{input_file};
		
		if($config->{input_file} =~ /\.gz$/){
			$in_file_handle->open("zcat ". $config->{input_file} . " | " ) or die("ERROR: Could not read from input file ", $config->{input_file}, "\n");
		}
		else {
			$in_file_handle->open( $config->{input_file} ) or die("ERROR: Could not read from input file ", $config->{in_file}, "\n");
		}
	}
	
	# no file specified - try to read data off command line
	else {
		$in_file_handle = 'STDIN';
		print("Reading input from STDIN (or maybe you forgot to specify an input file?)...\n");
	}
	
	return $in_file_handle;
}

sub get_intervals_file_handle {
	my $config = shift;
	
	# define filehandle to write to
	my $intervals_file_handle = new FileHandle;
	$intervals_file_handle->open(">".$config->{intervals_file}) or die("ERROR: Could not write to intervals file ", $config->{intervals_file}, "\n");
	
	return $intervals_file_handle;
}

sub get_bed_file_handle {
	my $config = shift;
	
	# define filehandle to write to
	my $bed_file_handle = new FileHandle;
	$bed_file_handle->open(">".$config->{bed_file}) or die("ERROR: Could not write to bed file ", $config->{bed_file}, "\n");
	
	return $bed_file_handle;
}

sub get_interval_list_file_handle {
	my $config = shift;
	
	# define filehandle to write to
	my $interval_list_file_handle = new FileHandle;
	$interval_list_file_handle->open(">".$config->{interval_list_file}) or die("ERROR: Could not write to interval list file ", $config->{details_file}, "\n");

	my $header =<<'HEADER_END';
@HD	VN:1.0	SO:coordinate
@SQ	SN:1	LN:249250621	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:1b22b98cdeb4a9304cb5d48026a85128
@SQ	SN:2	LN:243199373	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:a0d9851da00400dec1098a9255ac712e
@SQ	SN:3	LN:198022430	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:fdfd811849cc2fadebc929bb925902e5
@SQ	SN:4	LN:191154276	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:23dccd106897542ad87d2765d28a19a1
@SQ	SN:5	LN:180915260	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:0740173db9ffd264d728f32784845cd7
@SQ	SN:6	LN:171115067	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:1d3a93a248d92a729ee764823acbbc6b
@SQ	SN:7	LN:159138663	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:618366e953d6aaad97dbe4777c29375e
@SQ	SN:8	LN:146364022	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:96f514a9929e410c6651697bded59aec
@SQ	SN:9	LN:141213431	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:3e273117f15e0a400f01055d9f393768
@SQ	SN:10	LN:135534747	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:988c28e000e84c26d552359af1ea2e1d
@SQ	SN:11	LN:135006516	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:98c59049a2df285c76ffb1c6db8f8b96
@SQ	SN:12	LN:133851895	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:51851ac0e1a115847ad36449b0015864
@SQ	SN:13	LN:115169878	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:283f8d7892baa81b510a015719ca7b0b
@SQ	SN:14	LN:107349540	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:98f3cae32b2a2e9524bc19813927542e
@SQ	SN:15	LN:102531392	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:e5645a794a8238215b2cd77acb95a078
@SQ	SN:16	LN:90354753	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:fc9b1a7b42b97a864f56b348b06095e6
@SQ	SN:17	LN:81195210	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:351f64d4f4f9ddd45b35336ad97aa6de
@SQ	SN:18	LN:78077248	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:b15d4b2d29dde9d3e4f93d1d0f2cbc9c
@SQ	SN:19	LN:59128983	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:1aacd71f30db8e561810913e0b72636d
@SQ	SN:20	LN:63025520	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:0dec9660ec1efaaf33281c0d5ea2560f
@SQ	SN:21	LN:48129895	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:2979a6085bfe28e3ad6f552f361ed74d
@SQ	SN:22	LN:51304566	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:a718acaa6135fdca8357d5bfe94211dd
@SQ	SN:X	LN:155270560	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:7e0e2e580297b7764e31dbc80c2540dd
@SQ	SN:Y	LN:59373566	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:1fa3474750af0948bdf97d5a0ee52e51
@SQ	SN:MT	LN:16569	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:c68f52674c9fb33aef52dcf399755519
@SQ	SN:GL000207.1	LN:4262	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:f3814841f1939d3ca19072d9e89f3fd7
@SQ	SN:GL000226.1	LN:15008	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:1c1b2cd1fccbc0a99b6a447fa24d1504
@SQ	SN:GL000229.1	LN:19913	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:d0f40ec87de311d8e715b52e4c7062e1
@SQ	SN:GL000231.1	LN:27386	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:ba8882ce3a1efa2080e5d29b956568a4
@SQ	SN:GL000210.1	LN:27682	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:851106a74238044126131ce2a8e5847c
@SQ	SN:GL000239.1	LN:33824	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:99795f15702caec4fa1c4e15f8a29c07
@SQ	SN:GL000235.1	LN:34474	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:118a25ca210cfbcdfb6c2ebb249f9680
@SQ	SN:GL000201.1	LN:36148	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:dfb7e7ec60ffdcb85cb359ea28454ee9
@SQ	SN:GL000247.1	LN:36422	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:7de00226bb7df1c57276ca6baabafd15
@SQ	SN:GL000245.1	LN:36651	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:89bc61960f37d94abf0df2d481ada0ec
@SQ	SN:GL000197.1	LN:37175	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:6f5efdd36643a9b8c8ccad6f2f1edc7b
@SQ	SN:GL000203.1	LN:37498	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:96358c325fe0e70bee73436e8bb14dbd
@SQ	SN:GL000246.1	LN:38154	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:e4afcd31912af9d9c2546acf1cb23af2
@SQ	SN:GL000249.1	LN:38502	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:1d78abec37c15fe29a275eb08d5af236
@SQ	SN:GL000196.1	LN:38914	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:d92206d1bb4c3b4019c43c0875c06dc0
@SQ	SN:GL000248.1	LN:39786	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:5a8e43bec9be36c7b49c84d585107776
@SQ	SN:GL000244.1	LN:39929	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:0996b4475f353ca98bacb756ac479140
@SQ	SN:GL000238.1	LN:39939	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:131b1efc3270cc838686b54e7c34b17b
@SQ	SN:GL000202.1	LN:40103	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:06cbf126247d89664a4faebad130fe9c
@SQ	SN:GL000234.1	LN:40531	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:93f998536b61a56fd0ff47322a911d4b
@SQ	SN:GL000232.1	LN:40652	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:3e06b6741061ad93a8587531307057d8
@SQ	SN:GL000206.1	LN:41001	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:43f69e423533e948bfae5ce1d45bd3f1
@SQ	SN:GL000240.1	LN:41933	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:445a86173da9f237d7bcf41c6cb8cc62
@SQ	SN:GL000236.1	LN:41934	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:fdcd739913efa1fdc64b6c0cd7016779
@SQ	SN:GL000241.1	LN:42152	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:ef4258cdc5a45c206cea8fc3e1d858cf
@SQ	SN:GL000243.1	LN:43341	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:cc34279a7e353136741c9fce79bc4396
@SQ	SN:GL000242.1	LN:43523	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:2f8694fc47576bc81b5fe9e7de0ba49e
@SQ	SN:GL000230.1	LN:43691	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:b4eb71ee878d3706246b7c1dbef69299
@SQ	SN:GL000237.1	LN:45867	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:e0c82e7751df73f4f6d0ed30cdc853c0
@SQ	SN:GL000233.1	LN:45941	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:7fed60298a8d62ff808b74b6ce820001
@SQ	SN:GL000204.1	LN:81310	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:efc49c871536fa8d79cb0a06fa739722
@SQ	SN:GL000198.1	LN:90085	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:868e7784040da90d900d2d1b667a1383
@SQ	SN:GL000208.1	LN:92689	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:aa81be49bf3fe63a79bdc6a6f279abf6
@SQ	SN:GL000191.1	LN:106433	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:d75b436f50a8214ee9c2a51d30b2c2cc
@SQ	SN:GL000227.1	LN:128374	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:a4aead23f8053f2655e468bcc6ecdceb
@SQ	SN:GL000228.1	LN:129120	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:c5a17c97e2c1a0b6a9cc5a6b064b714f
@SQ	SN:GL000214.1	LN:137718	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:46c2032c37f2ed899eb41c0473319a69
@SQ	SN:GL000221.1	LN:155397	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:3238fb74ea87ae857f9c7508d315babb
@SQ	SN:GL000209.1	LN:159169	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:f40598e2a5a6b26e84a3775e0d1e2c81
@SQ	SN:GL000218.1	LN:161147	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:1d708b54644c26c7e01c2dad5426d38c
@SQ	SN:GL000220.1	LN:161802	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:fc35de963c57bf7648429e6454f1c9db
@SQ	SN:GL000213.1	LN:164239	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:9d424fdcc98866650b58f004080a992a
@SQ	SN:GL000211.1	LN:166566	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:7daaa45c66b288847b9b32b964e623d3
@SQ	SN:GL000199.1	LN:169874	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:569af3b73522fab4b40995ae4944e78e
@SQ	SN:GL000217.1	LN:172149	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:6d243e18dea1945fb7f2517615b8f52e
@SQ	SN:GL000216.1	LN:172294	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:642a232d91c486ac339263820aef7fe0
@SQ	SN:GL000215.1	LN:172545	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:5eb3b418480ae67a997957c909375a73
@SQ	SN:GL000205.1	LN:174588	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:d22441398d99caf673e9afb9a1908ec5
@SQ	SN:GL000219.1	LN:179198	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:f977edd13bac459cb2ed4a5457dba1b3
@SQ	SN:GL000224.1	LN:179693	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:d5b2fc04f6b41b212a4198a07f450e20
@SQ	SN:GL000223.1	LN:180455	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:399dfa03bf32022ab52a846f7ca35b30
@SQ	SN:GL000195.1	LN:182896	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:5d9ec007868d517e73543b005ba48535
@SQ	SN:GL000212.1	LN:186858	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:563531689f3dbd691331fd6c5730a88b
@SQ	SN:GL000222.1	LN:186861	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:6fe9abac455169f50470f5a6b01d0f59
@SQ	SN:GL000200.1	LN:187035	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:75e4c8d17cd4addf3917d1703cacaf25
@SQ	SN:GL000193.1	LN:189789	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:dbb6e8ece0b5de29da56601613007c2a
@SQ	SN:GL000194.1	LN:191469	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:6ac8f815bf8e845bb3031b73f812c012
@SQ	SN:GL000225.1	LN:211173	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:63945c3e6962f28ffd469719a747e73c
@SQ	SN:GL000192.1	LN:547496	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:325ba9e808f669dfeee210fdd7b470ac
@SQ	SN:NC_007605	LN:171823	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:6743bd63b3ff2b5b8985d8933c53290a
@SQ	SN:hs37d5	LN:35477943	UR:file:/humgen/1kg/reference/human_g1k_v37_decoy.fasta	M5:5b6a4b3a81a2d3c134b7d14bf6ad39f1
HEADER_END

	# add headers
	print $interval_list_file_handle $header;

	return $interval_list_file_handle;
}

sub get_details_file_handle {
	my $config = shift;
	
	# define filehandle to write to
	my $details_file_handle = new FileHandle;
	$details_file_handle->open(">".$config->{details_file}) or die("ERROR: Could not write to details file ", $config->{details_file}, "\n");
	
	# make header
    my $header = $program_header.$config->{config_options_header};
#	my $time = &get_time;
#	my $core_mca = $config->{reg}->get_adaptor($config->{species}, 'core', 'metacontainer');
#	my $db_string = $core_mca->dbc->dbname." on ".$core_mca->dbc->host if defined $core_mca;
#	my $version_string =
#		"Using API version ".$config->{reg}->software_version.
#		", DB version ".(defined $core_mca && $core_mca->get_schema_version ? $core_mca->get_schema_version : '?');
#	
#	my $header = $program_header.<<HEAD;
### Output produced at $time
### Connected to $db_string
### $version_string
#HEAD
	
	# add headers
	print $details_file_handle $header;
	
	# add column headers
	print $details_file_handle join "\t", (
		'#Query',
		'GeneSymbol',
		'GeneID',
		'TranscriptID',
        'Canonical',
		'ProteinID',
		'ExonID',
		'Region',
        'RegionNumber',
        'SubRegion',
		'Chromosome',
		'Start',
		'End',
		'Strand'
	);
	
	print $details_file_handle "\n";
	
	return $details_file_handle;
}

sub get_gene_file_handle {
	my $config = shift;
	
	# define filehandle to write to
	my $gene_file_handle = new FileHandle;
	$gene_file_handle->open(">".$config->{gene_file}) or die("ERROR: Could not write to gene file ", $config->{gene_file}, "\n");
	
	# make header
    my $header = $program_header.$config->{config_options_header};
#	my $time = &get_time;
#	my $core_mca = $config->{reg}->get_adaptor($config->{species}, 'core', 'metacontainer');
#	my $db_string = $core_mca->dbc->dbname." on ".$core_mca->dbc->host if defined $core_mca;
#	my $version_string =
#		"Using API version ".$config->{reg}->software_version.
#		", DB version ".(defined $core_mca && $core_mca->get_schema_version ? $core_mca->get_schema_version : '?');
#	
#	my $header = $program_header.<<HEAD;
### Gene list produced at $time
### Connected to $db_string
### $version_string
#HEAD
	
	# add headers
	print $gene_file_handle $header;
	
	# add column headers
	print $gene_file_handle join "\t", (
		'#Query',
		'GeneSymbol',
		'GeneID',
		'Chromosome',
		'Start',
		'End',
		'Strand',
		'Biotype',
		'Status',
		'Description',
	);
	
	print $gene_file_handle "\n";
	
	return $gene_file_handle;
}
sub get_query_size_file_handle {
	my $config = shift;
	
	# define filehandle to write to
	my $query_size_file_handle = new FileHandle;
	$query_size_file_handle->open(">".$config->{query_size_file}) or die("ERROR: Could not write to query size file ", $config->{query_size_file}, "\n");
	
	# make header
    my $header = $program_header.$config->{config_options_header};
#	my $time = &get_time;
#	my $core_mca = $config->{reg}->get_adaptor($config->{species}, 'core', 'metacontainer');
#	my $db_string = $core_mca->dbc->dbname." on ".$core_mca->dbc->host if defined $core_mca;
#	my $version_string =
#		"Using API version ".$config->{reg}->software_version.
#		", DB version ".(defined $core_mca && $core_mca->get_schema_version ? $core_mca->get_schema_version : '?');
#	
#	my $header = $program_header.<<HEAD;
### Query sizes produced at $time
### Connected to $db_string
### $version_string
#HEAD
	
	# add headers
	print $query_size_file_handle $header;
	
	# add column headers
	print $query_size_file_handle join "\t", (
		'#Query',
		'BaseCount'
	);
	
	print $query_size_file_handle "\n";
	
	return $query_size_file_handle;
}
sub get_log_file_handle {
	my $config = shift;
	
	# define filehandle to write to
	my $log_file_handle = new FileHandle;
	$log_file_handle->open(">".$config->{log_file}) or die("ERROR: Could not write to log file ", $config->{log_file}, "\n");
	
	# print header
    my $header = $program_header.$config->{config_options_header};
	print $log_file_handle $header;

	return $log_file_handle;
}

# gets time
sub get_time() {
	my @time = localtime(time());

	# increment the month (Jan = 0)
	$time[4]++;

	# add leading zeroes as required
	for my $i(0..4) {
		$time[$i] = "0".$time[$i] if $time[$i] < 10;
	}

	# put the components together in a string
	my $time =
 		($time[5] + 1900)."-".
 		$time[4]."-".
 		$time[3]." ".
		$time[2].":".
		$time[1].":".
		$time[0];

	return $time;
}

# Perl trim function to remove whitespace from the start and end of the string
sub trim($) {
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}
# Left trim function to remove leading whitespace
sub ltrim($) {
	my $string = shift;
	$string =~ s/^\s+//;
	return $string;
}
# Right trim function to remove trailing whitespace
sub rtrim($) {
	my $string = shift;
	$string =~ s/\s+$//;
	return $string;
}

#standard format for start,end
sub format_coords {
	my ($start, $end) = @_;
	
	if(!defined($start)) {
		return '-';
	}
	elsif(!defined($end)) {
		return $start;
	}
	elsif($start == $end) {
		return $start;
	}
	elsif($start > $end) {
		return $end.'-'.$start;
	}
	else {
		return $start.'-'.$end;
	}
}

#standard format for chromosome position
sub slicef {
	my $slice = shift;
	return sprintf("chr%s:%u-%u[%s]",$slice->seq_region_name,$slice->start,$slice->end,$slice->strand==1?"+":"-");
} # !slicef

#convert chromosome to number
sub chr2zerochr {
    my $chr = uc(shift);
    return $chr =~ /^\d$/ ? '0'.$chr : $chr;
}
#convert number to chromosome
sub zerochr2chr {
    my $chr = shift;
    return $chr  =~ /^0\d$/ ? substr($chr,1) : $chr;
}

