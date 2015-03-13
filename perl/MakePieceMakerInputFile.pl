#!/usr/bin/perl -w
use strict;
use warnings;

#my $input_file = 'C:\Users\Michael\Desktop\gene_list.txt';
#my $input_file = 'C:\Users\Michael\Documents\Lab\Retardation\Genes\OMIM_MR+Greenwood+Raymond.txt';
my $input_file = 'C:\Users\Michael\Documents\Ensembl\ensembl.jcl';
my $output_file = 'C:\Users\Michael\Documents\Lab\Retardation\Genes\OMIM_MR+Greenwood+Raymond.piecemaker';
my $format = 'fasta';
my $include_utr;
my $include_introns;
my $roi_extra_upstream = 25;
my $roi_extra_downstream = 10;
my $padding_upstream = 1000;
my $padding_downstream= 1000;
my $warnings_only;
my $ambiguity = 'ambiguity';
CreateOUTPUTFILE($input_file,$output_file,$format,$include_utr,$include_introns,$roi_extra_upstream,$roi_extra_downstream,$padding_upstream,$padding_downstream,$warnings_only,$ambiguity);

#  Arg[0]  : string $input_file reference (first word of line is gene's external name, or a position, e.g., "chr10:40,063,307-73,063,491", rest can be anything; lines ignored if first character is #)
#  Arg[1]  : string $output_file
#  Arg[2]  : string $format fasta|piecemaker
#  Arg[3]  : int (bool) $include_utr
#  Arg[4]  : int (bool) $include_introns
#  Arg[5]  : int $roi_extra_upstream
#  Arg[6]  : int $roi_extra_downstream
#  Arg[7]  : int $padding_upstream
#  Arg[8]  : int $padding_downstream
#  Arg[9]  : int (bool) $warnings_only no download
#  Arg[10] : int (bool) $ambiguity
sub CreateOUTPUTFILE{
    my ($input_file,$output_file,$format,$include_utr,$include_introns,$roi_extra_upstream,$roi_extra_downstream,$padding_upstream,$padding_downstream,$warnings_only,$ambiguity) = @_;
    if(!defined $input_file){
        print "input file missing\n";
        help();
    }
    elsif($input_file =~ /\?/){
        help();
    }
    else{
        if(!defined $output_file){
            $output_file = $input_file.'.piecemaker'
        }
        if(!defined $padding_upstream){
            $padding_upstream = 0;
        }
        if(!defined $padding_downstream){
            $padding_downstream=$padding_upstream;
        }
        if(!defined($roi_extra_upstream)){
            $roi_extra_upstream = 0;
        }
        if(!defined($roi_extra_downstream)){
            $roi_extra_downstream = 0;
        }
        my $msg;
        my @gene_list;
        my $gene_count = 0;
        my $exon_count = 0;
        my $cds_count =0;
        my $utr_count = 0;
        my $base_count = 0;
        my $min_roi_length = 1e12;
        my $max_roi_length = 0;
        my $max_indel_length = 0;
        open GENELISTFILE,"<",$input_file
            or die "#ERROR: Couldn't open '$input_file' for reading, stopped\n";
        while(<GENELISTFILE>){
            if($_ !~ /^\s*[#$]/){
                $_ =~ /^\s*(?:chr)?[\s:]*([XY]|[0-9]{1,2})\D+([,\d]+)\D+([,\d]+).*#*/i;
                if($1 && $2 && $3){
                    
                     my $pos = $1."|".$2."|".$3;
                     $pos =~ s/,//g;
                     push @gene_list, $pos;
                }
                else{
                    $_ =~ /^\s*([^#\s]+).*#*/;
                    if($1){
                        push @gene_list, $1;
                    }
                }
             }
        }
        open OUTPUTFILE,">",$output_file
            or die "#ERROR: Couldn't open '$output_file' for writing, stopped\n";
        out (sprintf "#START: %s\n#input file: %s\n#output file: %s\n#format: %s\n#include UTR: %s\n#include introns: %s\n#extra bases in region of interest of exons: %s upstream, %s downstream\n#padding outside region of interest: %s upstream, %s downstream\n#create output file: %s\n#use ambiguity codes in sequence: %s\n",
            now(),
            $input_file,
            $output_file,
            $format,
            defined $include_utr?'yes':'no',
            defined $include_introns?'yes':'no',
            $roi_extra_upstream,
            $roi_extra_downstream,
            $padding_upstream,
            $padding_downstream,
            defined $warnings_only?'no':'yes',
            defined $ambiguity?'yes':'no');
        #connect to ensmbl database
        use Bio::EnsEMBL::Registry;
        my $registry = 'Bio::EnsEMBL::Registry';
        $registry->load_registry_from_db(-host => 'ensembldb.ensembl.org',-user => 'anonymous');
        my $gene_adaptor = $registry->get_adaptor( 'human', 'core', 'gene');
        my $exon_adaptor = $registry->get_adaptor( 'human', 'core', 'exon');
        my $slice_adaptor = $registry->get_adaptor( 'human', 'core', 'slice');
        #create records for each gene fragment
        foreach my $id(@gene_list){
            my $slice;
            my @genes;
            my $gene;
            my @exons;
            my $exon;
            #Chromosome position
            if($id =~ /^([XY]|[0-9]{1,2})\|\d+\|\d+$/){
                my ($chr,$start,$end) = split /\|/, $id;
                if($start <= $end){
                    $slice = $slice_adaptor->fetch_by_region( 'chromosome', $chr,$start,$end);
                    if($slice){
                        @genes = @{$slice->get_all_Genes};
                        if(@genes){
                            foreach $gene(@genes){
                                @exons = @{ $gene->get_all_Exons };
                                if(@exons){
                                    $gene_count++;
                                    my ($exon_count_temp,$cds_count_temp,$utr_count_temp,$base_count_temp,$min_roi_length_temp,$max_roi_length_temp) =
                                        exons(\@exons,$gene->strand,$include_utr,$padding_upstream,$padding_downstream,$roi_extra_upstream,$roi_extra_downstream);
                                    $exon_count += $exon_count_temp;
                                    $cds_count += $cds_count_temp;
                                    $utr_count += $utr_count_temp;
                                    $base_count += $base_count_temp;
                                    if($min_roi_length_temp < $min_roi_length){
                                        $min_roi_length = $min_roi_length_temp;
                                    }
                                    if($max_roi_length_temp > $max_roi_length){
                                        $max_roi_length = $max_roi_length_temp;
                                    }
                                }
                                else{
                                    out (sprintf "#ERROR: %s (%s) in chr%s:%d-%d has no exons",$gene->external_id,$gene->stable_id,$chr,$start,$end);
                                }
                            }
                        }
                        else{
                            out (sprintf "#ERROR: chr%s:%d-%d has no genes\n",$chr,$start,$end);
                        }
                    }
                    else{
                        out (sprintf "#ERROR: chr%s:%d-%d isn't there\n",$chr,$start,$end);
                    }
                }
                else{
                        out (sprintf "#ERROR: chr%s:%d-%d start greater than end\n",$chr,$start,$end);
                }
            }
            #Ensmbl gene stable id
            elsif($id =~ /ENSG\d{11}/){
                $gene = $gene_adaptor->fetch_by_stable_id($id);
                if($gene){
                    @exons = @{ $gene->get_all_Exons };
                    if(@exons){
                        $gene_count++;
                        my ($exon_count_temp,$cds_count_temp,$utr_count_temp,$base_count_temp,$min_roi_length_temp,$max_roi_length_temp) =
                            exons(\@exons,$gene->strand,$include_utr,$padding_upstream,$padding_downstream,$roi_extra_upstream,$roi_extra_downstream);
                        $exon_count += $exon_count_temp;
                        $cds_count += $cds_count_temp;
                        $utr_count += $utr_count_temp;
                        $base_count += $base_count_temp;
                        if($min_roi_length_temp < $min_roi_length){
                            $min_roi_length = $min_roi_length_temp;
                        }
                        if($max_roi_length_temp > $max_roi_length){
                            $max_roi_length = $max_roi_length_temp;
                        }
                    }
                    else{
                        out (sprintf "#ERROR: %s has no exons",$id);
                    }
                }
                else{
                    out (sprintf "#ERROR: %s isn't there\n",$id);
                }
            }
            #Ensmbl exon stable id
            elsif($id =~ /ENSE\d{11}/){
                $exon = $exon_adaptor->fetch_by_stable_id($id);
                if($exon){
                    $gene = $gene_adaptor->fetch_by_exon_stable_id($exon->stable_id);
                    $gene_count++;
                    if($include_utr || $exon->phase!=-1){
                        my $roi_length = exon($exon,0,$padding_upstream,$padding_downstream,$roi_extra_upstream,$roi_extra_downstream);
                        $exon_count++;
                        $exon->phase==-1?$utr_count++:$cds_count++;
                        $base_count += $roi_length;
                        if($roi_length < $min_roi_length){
                            $min_roi_length = $roi_length;
                        }
                        if($roi_length > $max_roi_length){
                            $max_roi_length = $roi_length;
                        }
                    }
                }
                else{
                    out (sprintf "#ERROR: %s isn't there\n",$id);                 }
                }
            #External gene name
            else{
                @genes = @{$gene_adaptor->fetch_all_by_external_name($id)};
                if(@genes){
                    $gene = $genes[0];
                    if(@genes > 1){
                        foreach my $g(@genes){
                            if($g->external_name eq $id){
                                $gene = $g;
                                last;
                            }
                        }
                        my $count_temp = @genes;
                        my $external_name = $gene->external_name;
                        my $stable_id = $gene->stable_id;
                        out (sprintf "#WARNING: %s identifies more than one gene (%d genes); using %s (%s)\n",$id,$count_temp,$external_name,$stable_id);
                        foreach my $g(@genes){
                            out (sprintf"#\t%s  chr%s:%d-%d  %s  (%s)\n",
                                 $g->external_name,$gene->slice->seq_region_name,
                                 $gene->start,
                                 $gene->end,
                                 $g->description,
                                 $g->stable_id);
                        }
                     }
                    if($gene->external_name ne $id){
                        out (sprintf "#WARNING: %s will be named %s  chr%s:%d-%d  %s  (%s)\n",
                             $id,
                             $gene->external_name,
                             $gene->slice->seq_region_name,
                             $gene->start,$gene->end,
                             $gene->description,
                             $gene->stable_id);
                    }
                    if(!$warnings_only){
                    @exons = @{ $gene->get_all_Exons };
                     if(@exons){
                        $gene_count++;
                        my ($exon_count_temp,$cds_count_temp,$utr_count_temp,$base_count_temp,$min_roi_length_temp,$max_roi_length_temp) =
                            exons(\@exons,$gene->strand,$include_utr,$padding_upstream,$padding_downstream,$roi_extra_upstream,$roi_extra_downstream);
                        $exon_count += $exon_count_temp;
                        $cds_count += $cds_count_temp;
                        $utr_count += $utr_count_temp;
                        $base_count += $base_count_temp;
                        if($min_roi_length_temp < $min_roi_length){
                            $min_roi_length = $min_roi_length_temp;
                        }
                        if($max_roi_length_temp > $max_roi_length){
                            $max_roi_length = $max_roi_length_temp;
                        }
                    }
                    else{
                        out (sprintf "#ERROR: %s has no exons",$id);
                    }
                    }
                }
                else{
                    out (sprintf "#ERROR: %s isn't there\n",$id);
                }
            }
        }
        out (sprintf "#DONE: %s\n#  %d gene%s; %d exon%s (%d CDS, %d UTR); %d base%s; ROI lengths %d-%d\n",now(),$gene_count,ss($gene_count),$exon_count,ss($exon_count),$cds_count,$utr_count,$base_count,ss($base_count),$min_roi_length,$max_roi_length);
        close(OUTPUTFILE)
            or die "#ERROR: Couldn't close '$output_file', stopped\n";
        return 1;
    } # !main sub logic ($input_file)
    
    sub out{
        my $out = $_[0];
        print $out;
        if(!$warnings_only){
        print OUTPUTFILE $out
            or die "#ERROR: Couldn't write to output file, stopped\n";
        }
    }
    
    #sort exons and write a record for each
    sub exons{
        my ($exons,$strand,$include_utr,$padding_upstream,$padding_downstream,$roi_extra_upstream,$roi_extra_downstream) = @_;
        # @exons is unordered and may have overlapping exons
        if($strand==1){
            @$exons = sort {$a->start <=> $b->start || $a->end <=> $b->end} @$exons;
        }
        else{
            @$exons = sort {$b->start <=> $a->start || $b->end <=> $a->end} @$exons;
        }
        my $exon_count=0;
        my $cds_count=0;
        my $utr_count=0;
        my $base_count=0;
        my $min_roi_length = 1e12;
        my $max_roi_length = 0;
        my $i = 1;
                
        foreach my $exon(@$exons){
            if($include_utr || $exon->phase!=-1){
                my $roi_length = exon($exon,$i,$padding_upstream,$padding_downstream,$roi_extra_upstream,$roi_extra_downstream);
                $exon_count++;
                $exon->phase==-1?$utr_count++:$cds_count++;
                $base_count += $roi_length;
                if($roi_length < $min_roi_length){
                    $min_roi_length = $roi_length;
                }
                if($roi_length > $max_roi_length){
                    $max_roi_length = $roi_length;
                }
                $i++;
            }
        }
        return ($exon_count,$cds_count,$utr_count,$base_count,$min_roi_length,$max_roi_length);
    } # !exons
    
    #write a record for each exon
    sub exon{
        my ($exon,$exon_number,$padding_upstream,$padding_downstream,$roi_extra_upstream,$roi_extra_downstream) = @_;
        my $exon_slice = $exon->slice;
        my $registry = 'Bio::EnsEMBL::Registry';
        my $gene_adaptor = $registry->get_adaptor( 'human', 'core', 'gene');
        my $gene = $gene_adaptor->fetch_by_exon_stable_id($exon->stable_id);
        my $slice_adaptor = $registry->get_adaptor( 'human', 'core', 'slice');
        my $slice = $slice_adaptor->fetch_by_region( 'chromosome', $exon_slice->seq_region_name, $exon->start-$padding_upstream-$roi_extra_upstream, $exon->end+$padding_downstream+$roi_extra_downstream, $exon->strand);
        my $sequence_with_ambiguity = fetch_sequence_with_ambiguity($exon_slice->seq_region_name, $exon->start-$padding_upstream-$roi_extra_upstream, $exon->end+$padding_downstream+$roi_extra_downstream, $exon->strand);
        my $roi_start=$exon->start-$roi_extra_upstream-$slice->start+1;
        my $roi_end=$exon->end+$roi_extra_downstream-$slice->start+1;
        my $roi_length = $roi_end - $roi_start +1;
        my $fasta_description = sprintf "%s_e%02u%s%s|%s|%s|chr%s:%u-%u|%s|%s|%s",
            $gene->external_name,
            $exon_number,
            $roi_extra_upstream?"-".$roi_extra_upstream:"",
            $roi_extra_downstream?"+".$roi_extra_downstream:"",
            $exon->phase==-1?"UTR":"CDS",
            $exon->strand==1?"+":"-",
            $exon->seq_region_name,
            $exon->start,
            $exon->end,
            $gene->stable_id,
            $exon->stable_id,
            $roi_length;
        if (length $sequence_with_ambiguity != length $slice->seq){
            out (sprintf "#WARNING: %s exon %d sequence length (including extra bases and padding) is %d but the length with ambiguity is %d\n",$gene->external_name,$exon_number,length($slice->seq),length($sequence_with_ambiguity));
            print $sequence_with_ambiguity."\n\n------------------------\n\n";
            print $slice->seq."\n\n------------------------\n";
        }
        printf "%s\t%u\t%u\t%u\t%s...%s\n",
            $fasta_description,
            $slice->start,
            $roi_start,
            $roi_end,
            substr($sequence_with_ambiguity,0,10),
            substr($sequence_with_ambiguity,-10);
        printf OUTPUTFILE "%s_e%02u\t%s\t%u\t%u\t%u\t%s\n",
            $gene->external_name,
            $exon_number,
            $fasta_description,
            $slice->start,
            $roi_start,
            $roi_end,
            $sequence_with_ambiguity;
        return $roi_length;
    } # !exon
    
    sub now{
        use Time::localtime;
        my $tm = localtime;
        return sprintf "%04d-%02d-%02d %02d:%02d:%02d",$tm->year+1900,$tm->mon+1,$tm->mday,$tm->hour,$tm->min,$tm->sec;
    } # !now
    
    sub ss{
        return $_[0] == 1?"":"s";
    } # !ss
    
    sub help{
        my $help = "usage:\n";
        $help .= "<input_file>,[<output_file>,][[fasta|piecemaker],][<include_utr>,][<include_introns>,][<padding_upstream>,][<padding_downstream>,][<roi_extra_upstream>,][<roi_extra_downstream>,][warnings_only,][ambiguity]\n";
        $help .= "input_file = path of file containg gene names, Ensembl ID's, or chromosome positions (required\n";
        $help .= "output_file = path of output file (default is <input_file>.<format>)\n";
        $help .= "format = fasta or piecemaker\n";
        $help .= "include_utr = include untranslated regions of exons\n";
        $help .= "include_introns = include introns\n";
        $help .= "roi_extra_upstream = number of extra bases to add to region of interest upstream of exons\n";
        $help .= "roi_extra_downstream = number of extra bases to add to region of interest downstream of exons\n";
        $help .= "padding_upstream = number of bases to add upstream of region of interest (needed for piecemaker)\n";
        $help .= "padding_downstream = number of bases to add downstream of region of interest (needed for piecemaker)\n";
        $help .= "warnings_only if defined don't download or create output file\n";
        $help .= "ambiguity if defined use ambiguity codes in sequence\n";
        $help .= "emitted sequence will be:\n";
        $help .= "5'[padding_upstream][roi_extra_upstream]<exon>[roi_extra_downstream][padding_downstream]3'\n";
        $help .= "or\n";
        $help .= "5'[padding_upstream]<intron>[padding_downstream]3'\n";    
        $help .= "fragment description format (used as fasta record header):\n";
        $help .= "<gene external name_<e or i><exon or intron number>-<extra upstream bases in ROI>+<extra downstream bases in ROI>|<CDS or UTR>|<strand>|chr<chromosome>:<exon start position>-<exon end position>|<gene stable name>|<exon stable name>|<ROI length>\n";
        print $help;
    } # !help

    sub fetch_sequence_with_ambiguity {
        my ($chr,$start,$end,$strand) = @_;
        use Bio::EnsEMBL::Registry;
        my $registry = 'Bio::EnsEMBL::Registry';
        $registry->load_registry_from_db(
            -host => 'ensembldb.ensembl.org',
            -user => 'anonymous'
        );
        use Bio::EnsEMBL::Variation::Utils::Sequence qw (sequence_with_ambiguity);
        my $dbCore = $registry->get_DBAdaptor( 'human', 'core');
        my $dbVar = $registry->get_DBAdaptor( 'human', 'variation' );
        my $ambiguous_slice = sequence_with_ambiguity($dbCore,$dbVar,$chr,$start,$end,$strand);
        my $ambiguous_seq = $ambiguous_slice->seq();
        return $ambiguous_seq;
    }

} # !CreateOUTPUTFILE

