#!/usr/bin/perl -w
use strict;
use warnings;

#my $input_file = 'C:\Users\Michael\Documents\Lab\Retardation\Genes\OMIM_MR+Greenwood+Raymond.txt';
my $input_file = 'C:\Users\Michael\Desktop\gene_list.txt';
my $output_file = 'C:\Users\Michael\Documents\Lab\Retardation\Genes\OMIM_MR+Greenwood+Raymond.piecemaker';
my $include_utr = 1;
my $padding = 1000;
my $padding_upstream = 1000;
my $padding_downstream= 1000;
my $roi_extra_upstream = 25;
my $roi_extra_downstream = 10;
CreatePieceMakerFile($input_file,$output_file,$include_utr,$padding_upstream,$padding_downstream,$roi_extra_upstream,$roi_extra_downstream);

#  Arg[1]      : string $input_file reference (first word of line is gene's external name, or a position, e.g., "chr10:40,063,307-73,063,491", rest can be anything; lines ignored if first character is #)
#  Arg[2]      : string $output_file
#  Arg[3]      : int (bool) $include_utr
#  Arg[4]      : int $padding_upstream
#  Arg[5]      : int $padding_downstream
#  Arg[6]      : int $roi_extra_upstream
#  Arg[7]      : int $roi_extra_downstream
sub CreatePieceMakerFile{
    my ($input_file,$output_file,$include_utr,$padding_upstream,$padding_downstream,$roi_extra_upstream,$roi_extra_downstream) = @_;
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
    my @pos_list;
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
        if($_ !~ /^\s*#/){
            $_ =~ /^\s*(?:chr)?[\s:]*([XY]|[0-9]{1,2})\D+([,\d]+)\D+([,\d]+).*#*/i;
            if($1 && $2 && $3){
                
                 my $pos = $1."|".$2."|".$3;
                 $pos =~ s/,//g;
                 push @pos_list, $pos;
            }
            else{
                $_ =~ /^\s*([^#\s]+).*#*/;
                if($1){
                    push @gene_list, $1;
                }
            }
         }
    }
    open PIECEMAKERFILE,">",$output_file
        or die "#ERROR: Couldn't open '$output_file' for writing, stopped\n";
    $msg = sprintf
"#START: %s
#input file: %s
#output file: %s
#include UTR: %s
#extra bases in region of interest: %s upstream, %s downstream
#padding outside region of interest: %s upstream, %s downstream\n",
now(),
$input_file,
$output_file,
$include_utr?'yes':'no',
$roi_extra_upstream,
$roi_extra_downstream,
$padding_upstream,
$padding_downstream;
    print $msg;
    print PIECEMAKERFILE $msg
        or die "#ERROR: Couldn't write to '$output_file', stopped\n";
    $msg = "#Fragment description format: <gene external name_<e or i><exon or intron number>-<extra upstream bases in ROI>+<extra downstream bases in ROI>|<CDS or UTR>|<gene stable name>|<exon stable name>|<strand>|chr<chromosome>:<exon start position>-<exon end position>|<ROI length>\n\n";
    print $msg;
    print PIECEMAKERFILE $msg
        or die "#ERROR: Couldn't write to '$output_file', stopped\n";
    print PIECEMAKERFILE "#Target\tFragment_Description\tFragment_Start_Position\tRegion_of_Interest_Begin_Relative_Position\tRegion_of_Interest_End_Relative_Position\tSequence\n"
        or die "#ERROR: Couldn't write to '$output_file', stopped\n";
    #connect to ensmbl database
    use Bio::EnsEMBL::Registry;
    my $registry = 'Bio::EnsEMBL::Registry';
    $registry->load_registry_from_db(-host => 'ensembldb.ensembl.org',-user => 'anonymous');
    my $gene_adaptor = $registry->get_adaptor( 'human', 'core', 'gene');
    my $exon_adaptor = $registry->get_adaptor( 'human', 'core', 'exon');
    my $slice_adaptor = $registry->get_adaptor( 'human', 'core', 'slice');
    #create records for each gene fragment
    foreach my $id(@gene_list){
        my @genes;
        my $gene;
        my $exon;
        my $slice;
        #Chromosome position
        if($id =~ /[XY]|[0-9]{1,2}\|\d+\|\d+/){
            my ($chr,$start,$end) = split /\|/, $id;
            $slice = $slice_adaptor->fetch_by_region( 'chromosome', $chr,$start,$end);
        }
        #Ensmbl gene stable id
        elsif($id =~ /ENSG\d{15}/){
            $gene = $gene_adaptor->fetch_by_gene_stable_id($id);
        }
        #Ensmbl exon stable id
        elsif($id =~ /ENSE\d{15}/)
        {
            $exon = $exon_adaptor->fetch_by_stable_id($id);
        }
        #External name
        else{
            @genes = @{$gene_adaptor->fetch_all_by_external_name($id)};
            if(@genes){
                $gene = $genes[0];
                if(@genes > 1){
                    my @stable_ids;
                    foreach my $g(@genes){
                        push @stable_ids,$g->stable_id;
                    }
                    $msg = sprintf "#WARNING: %s identifies %d genes (%s); using %s",$id,@genes,@stable_ids,$gene->stable_id;
                    print $msg;
                    print PIECEMAKERFILE $msg
                        or die "#ERROR: Couldn't write to '$output_file', stopped\n";
                }
                if($gene->external_name ne $id){
                    $msg = sprintf "#WARNING: %s will be named %s\n",$id,$gene->external_name;
                    print $msg;
                    print PIECEMAKERFILE $msg
                        or die "#ERROR: Couldn't write to '$output_file', stopped\n";
                }
            }
        }
        if($gene){
            if(gene($gene)){
                $gene_count++;
            }
        }
        elsif($exon){
            if(exon($exon)){
                $exon_count++;
            }
        }
        elsif($slice){
            if(slice($slice)){
            }
        }
        else{
            print "#ERROR: %s isn't there\n".$_;
        }
    }
    $msg = sprintf "#DONE: %s\n#  %d gene%s; %d exon%s (%d CDS, %d UTR); %d base%s; ROI lengths %d-%d\n",now(),$gene_count,ss($gene_count),$exon_count,ss($exon_count),$cds_count,$utr_count,$base_count,ss($base_count),$min_roi_length,$max_roi_length;
    print PIECEMAKERFILE $msg
        or die "#ERROR: Couldn't write to '$output_file', stopped\n";
    close(PIECEMAKERFILE)
        or die "#ERROR: Couldn't close '$output_file', stopped\n";
    print $msg;
    
    sub exon{
        my $exon = $_[0];
    }
}

sub gene{
    my $gene = $_[0];
    slice($gene);
            my $gene_chr = $gene->slice->seq_region_name;
            my @exons = @{ $gene->get_all_Exons };
            if(@exons){
                $gene_count++;
                # @exons is unordered and may have overlapping exons
                if($gene->strand==1){
                    @exons = sort {$a->start <=> $b->start || $a->end <=> $b->end} @exons;
                }
                else{
                    @exons = sort {$b->start <=> $a->start || $b->end <=> $a->end} @exons;
                }
                my $i = 1;
                foreach my $exon(@exons){
                    exon(@exons);
                    if($include_utr || $exon->phase!=-1){
                        my $slice_adaptor = $registry->get_adaptor( 'human', 'core', 'slice');
                        my $slice = $slice_adaptor->fetch_by_region( 'chromosome', $gene_chr, $exon->start-$padding_upstream-$roi_extra_upstream, $exon->end+$padding_downstream+$roi_extra_downstream, $exon->strand);
                        my $sequence_with_ambiguity = fetch_sequence_with_ambiguity($gene_chr, $exon->start-$padding_upstream-$roi_extra_upstream, $exon->end+$padding_downstream+$roi_extra_downstream, $exon->strand);
                        my $roi_start=$exon->start-$roi_extra_upstream-$slice->start+1;
                        my $roi_end=$exon->end+$roi_extra_downstream-$slice->start+1;
                        my $roi_length = $roi_end - $roi_start +1;
                        my $fasta_description = sprintf "%s_e%02u%s%s|%s|%s|%s|%s|chr%s:%u-%u|%s",$gene->external_name,$i,$roi_extra_upstream?"-".$roi_extra_upstream:"",$roi_extra_downstream?"+".$roi_extra_downstream:"",$exon->phase==-1?"UTR":"CDS",$gene->stable_id,$exon->stable_id,$exon->strand==1?"+":"-",$exon->seq_region_name,$exon->start,$exon->end,$roi_length;
                        if (length $sequence_with_ambiguity != length $slice->seq){
                            $msg = sprintf "#WARNING: %s exon %d sequence length (including extra bases and padding) is %d but the length with ambiguity is %d",$gene->external_name,$i,length($slice->seq),length($sequence_with_ambiguity);
                            print $msg;
                            print PIECEMAKERFILE $msg
                                or die "#ERROR: Couldn't write to '$output_file', stopped\n";
                        }
                        printf "%s\t%u\t%u\t%u\t%s...%s\n",$fasta_description,$slice->start,$roi_start,$roi_end,substr($sequence_with_ambiguity,0,10),substr($sequence_with_ambiguity,-10);
                        # e.g.,  "ACSL4_e01  ACSL4 108770220    1001    0   "
                        printf PIECEMAKERFILE "%s_e%02u\t%s\t%u\t%u\t%u\t%s\n",$gene->external_name,$i,$fasta_description,$slice->start,$roi_start,$roi_end,$sequence_with_ambiguity;
                        $i++;
                        $exon->phase==-1?$utr_count++:$cds_count++;
                        $exon_count = $utr_count + $cds_count;
                        if($roi_length < $min_roi_length){
                            $min_roi_length = $roi_length;
                        }
                        if($roi_length > $max_roi_length){
                            $max_roi_length = $roi_length;
                        }
                        $base_count += $roi_length;
                    }
                 }
            }
            else{
                printf "#ERROR: %s has no exons\n".$gene->external_name;
            }
}


sub slice{
    my $slice = $_[0];
}

sub now{
    use Time::localtime;
    my $tm = localtime;
    return sprintf "%04d-%02d-%02d %02d:%02d:%02d",$tm->year+1900,$tm->mon,$tm->mday,$tm->hour,$tm->min,$tm->sec;
}

sub ss{
    return $_[0] == 1?"":"s";
}
sub fetch_sequence_with_ambiguity {
    my ($chr,$start,$end,$strand) = @_;
    use Bio::EnsEMBL::Registry;
    my $registry = 'Bio::EnsEMBL::Registry';
    $registry->load_registry_from_db(
        -host => 'ensembldb.ensembl.org',
        -user => 'anonymous'
    );
    
    my $dbCore = $registry->get_DBAdaptor( 'human', 'core');
    my $dbVar = $registry->get_DBAdaptor( 'human', 'variation' );
    my $ambiguous_slice = sequence_with_ambiguity($dbCore,$dbVar,$chr,$start,$end,$strand);
    my $ambiguous_seq = $ambiguous_slice->seq();
    return $ambiguous_seq;
}

#change standard behavior so any SNP with "-" becomes "-" not ""
sub ambiguity_code {
    my $alleles = shift;
    if($alleles =~ /.*-.*/){
        return "-";
    }
    my %duplicates; #hash containing all alleles to remove duplicates
    map {$duplicates{$_}++} split /[\|\/\\]/, $alleles;
    $alleles = uc( join '', sort keys %duplicates );
    my %ambig = qw(AC M ACG V ACGT N ACT H AG R AGT D AT W CG S CGT B CT Y GT K C C A A T T G G - - ); #we will need to decide what to do with alleles like -A. Is that possible??
    return $ambig{$alleles};
}

sub ambig_code {
    my $self = shift;
    return &ambiguity_code($self->allele_string());
}

sub sequence_with_ambiguity {
    my ($dbCore,$dbVar,$chr,$start,$end,$strand) = @_;
    my $slice;
    if (ref($dbCore) ne 'Bio::EnsEMBL::DBSQL::DBAdaptor'){
	warning('You need to provide a Bio::EnsEMBL::DBSQL::DBAdaptor as a first argument');
	return $slice;
    }
    if (ref($dbVar) ne 'Bio::EnsEMBL::Variation::DBSQL::DBAdaptor'){
	warning('You need to provide a Bio::EnsEMBL::Variation::DBSQL::DBAdaptor object as second argument');
	return $slice;
    }
    my $slice_adaptor = $dbCore->get_SliceAdaptor();
    my $vf_adaptor = $dbVar->get_VariationFeatureAdaptor;
    $slice = $slice_adaptor->fetch_by_region('chromosome',$chr,$start,$end,$strand); #get the slice
    my $seq = $slice->seq;
    foreach my $vf (@{$vf_adaptor->fetch_all_by_Slice($slice)}){
	substr($seq,$vf->start-1,1,ambig_code $vf);
    }
    $slice->{'seq'} = $seq;
    return $slice;
}
