#!/usr/bin/perl -w
use warnings;
use strict;
use File::Glob;
use File::Spec;

#creates database-friendly files from SVA exports
#by normalizing column names, filling in datr for all transcripts,
#converting "-" and " " to an empty string
#converting percentages to floats, creating 2 columns for hom:het and gene[transcript]
#and concatenating multiple files of the same type

#variants indel:
#INDEL_ID,INDEL_SN,#Subjects,Chromosome,Type,Start,End,Size,Bases,Indel_quality,NumReads,NumReads_INDEL,Function,Gene[Transcript],Pathway,Gene Ontology,Relevant_refSNP,Relevant_Venter,In_Control,Hom:Het
#variants snv:
#SNV_ID,SNV_SN,#Subjects,Chromosome,Position,Allele,Ref_allele,Strand (Illumina),Avg_Consensus_score,Avg_SNP_quality,Avg_RMS_score,Avg_Read_depth,Function,Transcript,Pathway,Gene_Ontology,Exisiting_rs_#,In_HapMap,In_Illumina1M,In_Venter,In_Control,Hom:Het
#variants sv:
#SV_ID,SV_SN,#Subjects,Chromosome,Type,Start,End,Size,NumReads,Ambiguous_Start,Ambiguous_End,SV_DGV_Overlap,DGV_ID,SV_Percentage,DGV_Percentage,Function,Gene[Transcript],Pathway,Gene_Ontology,Method,In_Control
#variants indel subjects:
#INDEL_ID,INDEL_SN,Type,Chromosome,Start,End,Bases,Relevant_refSNP,Relevant_Venter,Subject_ID,Hom_status,Indel_quality,Consensus_score,RMS_score,Num_reads,Num_reads_indel,Num_reads_ref,Function,Gene[Transcript]
#variants snv subjects:
#SNV_ID,SNV_SN,Subject_ID,Alleles,Ref_allele,Illumina_Strand,Consensus_score,SNP_quality,RMS_score,Read_depth,Num_reads_SNP,Num_reads_ref,Existing_rs_#,Function,Gene[Transcript]
#variants sv subjects:
#SV_ID,SV_SN,Type,Chromosome,Start,End,Length,Subject_ID,Status,Calling_method,Function,Gene[Transcript]
my $usage = "usage: sva2db <base path/name for output files> <list of input files>\n";
if (@ARGV < 2){
    print STDERR $usage;
    exit;
}
my $outbase = shift;

my %id=(
    INDEL_ID => "INDEL",
    SNV_ID => "SNV",
    SV_ID => "SV",
);

my %header = (
    INDEL => ["analysis", "ucsc_coord","gene","transcript","transcript_index","indel_id","indel_sn","subjects_count","chromosome","type","start_indel","end_indel","size","bases","indel_quality","num_reads","num_reads_indel","function_indel","gene_transcript","pathway","gene_ontology","relevant_refsnp","relevant_venter","in_control","hom","het"],
    SNV => ["analysis", "ucsc_coord","gene","transcript","transcript_index","snv_id","snv_sn","subjects_count","chromosome","position","allele","ref_allele","illumina_strand","avg_consensus_score","avg_snp_quality","avg_rms_score","avg_read_depth","function_snv","gene_transcript","pathway","gene_ontology","exisiting_rs_num","in_hapmap","in_illumina1m","in_venter","in_control","hom","het"],
    SV => ["analysis", "ucsc_coord","gene","transcript","transcript_index","sv_id","sv_sn","subjects_count","chromosome","type","start_sv","end_sv","size","numreads","ambiguous_start","ambiguous_end","sv_dgv_overlap","dgv_id","sv_percentage","dgv_percentage","function_sv","gene_transcript","pathway","gene_ontology","calling_method","in_control"],
    INDEL_SUBJECTS => ["analysis", "ucsc_coord","gene","transcript","transcript_index","indel_id","indel_sn","type","chromosome","start_indel","end_indel","bases","relevant_refsnp","relevant_venter","subject_id","hom_status","indel_quality","consensus_score","rms_score","num_reads","num_reads_indel","num_reads_ref","function_indel","gene_transcript"],
    SNV_SUBJECTS => ["analysis", "ucsc_coord","gene","transcript","transcript_index","snv_id","snv_sn","subject_id","alleles","ref_allele","illumina_strand","consensus_score","snp_quality","rms_score","read_depth","num_reads_snp","num_reads_ref","existing_rs_num","function_snv","gene_transcript"],
    SV_SUBJECTS =>, ["analysis", "ucsc_coord","gene","transcript","transcript_index","sv_id","sv_sn","type","chromosome","start_sv","end_sv","size","subject_id","status","calling_method","function_sv","gene_transcript"]
);
my %lastcontinuationfieldindex = (
    INDEL => 11,
    SNV => 11,
    SV => 14,
    INDEL_SUBJECTS => -1,
    SNV_SUBJECTS => -1,
    SV_SUBJECTS => -1
);
my %genefieldindex = (
    INDEL => 13,
    SNV => 13,
    SV => 16,
    INDEL_SUBJECTS => 18,
    SNV_SUBJECTS => 14,
    SV_SUBJECTS => 11
);
my %fieldcount = (
    INDEL => 20,
    SNV => 22,
    SV => 21,
    INDEL_SUBJECTS => 19,
    SNV_SUBJECTS => 15,
    SV_SUBJECTS => 12
);
my %percentindices = (
    INDEL => [],
    SNV => [],
    SV => [13,14],
    INDEL_SUBJECTS => [],
    SNV_SUBJECTS => [],
    SV_SUBJECTS => []
);
my %nocopyindices = (
    INDEL => [15],
    SNV => [],
    SV => [],
    INDEL_SUBJECTS => [],
    SNV_SUBJECTS => [],
    SV_SUBJECTS => []
);
my %colonindices = (
    INDEL => [19],
    SNV => [21],
    SV => [],
    INDEL_SUBJECTS => [],
    SNV_SUBJECTS => [],
    SV_SUBJECTS => []
);
my %filehandle = (
    #INDEL => *INDEL,
    #SNV => *SNV,
    #SV => *SV,
    #INDEL_SUBJECTS => *INDELSUBJECTS,
    #SNV_SUBJECTS => *SNVSUBJECTS,
    #SV_SUBJECTS => *SVSUBJECTS
);

my %inheader = (
    INDEL => 1,
    SNV => 1,
    SV => 1,
    INDEL_SUBJECTS => 1,
    SNV_SUBJECTS => 1,
    SV_SUBJECTS => 1   
);
for my $csv(glob shift){
    if(defined $csv && -e $csv){
        my ($volume,$directories,$file) = File::Spec->splitpath($csv);
        my $csvheader = 1;
        #my @dots = split /\./,$file;
        #if ($dots[-1] eq "csv" || $dots[-1] eq "txt"){
            #pop @dots;
        #}
        #my $filter = join ".",@dots;
        my $filetype = "";
        my $filehandle;
        my $transcriptindex = 0;
        my @thatrecord = ();
        open IN,"<",$csv or die "can't open $csv $!";
        my $linecount = 0;
        print "$csv\n";
        while(my $line = <IN>){
            $linecount++;
            chomp $line;
            unless(!defined $line || length $line == 0){
                my @thisrecord = split /,/, $line;
                if($csvheader){
                    $csvheader = 0;
                    unless(exists($id{$thisrecord[0]})){
                        print STDERR "$csv has unrecognized header: [$line]\n";
                        last;
                    }
                    $filetype = $id{$thisrecord[0]};
                    $filetype .= $thisrecord[2] eq "#Subjects" ? "" : "_SUBJECTS";
                    unless(exists $filehandle{$filetype}){
                        if($filetype eq "INDEL"){
                            open INDEL,">",$outbase."_indel.txt" or die "can't open $outbase-indel.txt $!";
                            $filehandle{$filetype} = *INDEL;
                        }
                        elsif($filetype eq "SNV")
                        {
                            open SNV,">",$outbase."_snv.txt" or die "can't open $outbase-snv.txt $!";
                            $filehandle{$filetype} = *SNV;
                        }
                         elsif($filetype eq "SV")
                        {
                            open SV,">",$outbase."_sv.txt" or die "can't open $outbase-sv.txt $!";
                            $filehandle{$filetype} = *SV;
                        }
                        elsif($filetype eq "INDEL_SUBJECTS")
                        {
                            open INDELSUBJECTS,">",$outbase."_indel_subjects.txt" or die "can't open $outbase-indel_subjects.txt $!";
                            $filehandle{$filetype} = *INDELSUBJECTS;
                        }
                        elsif($filetype eq "SNV_SUBJECTS")
                        {
                            open SNVSUBJECTS,">",$outbase."_snv_subjects.txt" or die "can't open $outbase-snv_subjects.txt $!";
                            $filehandle{$filetype} = *SNVSUBJECTS;
                        }
                        elsif($filetype eq "SV_SUBJECTS")
                        {
                            open SVSUBJECTS,">",$outbase."_sv_subjects.txt" or die "can't open $outbase-sv_subjects.txt $!";
                            $filehandle{$filetype} = *SVSUBJECTS;
                        }
                        else{
                            print STDERR "unrecognized file type $thisrecord[0] in file [$csv]";
                            last;
                        }
                    }
                    $filehandle = $filehandle{$filetype};
                   if($inheader{$filetype}){
                        $inheader{$filetype} = 0;
                        printf $filehandle "%s\n", join "\t",@{$header{$filetype}};
                    }
                }
                else{
                    if(@thisrecord != $fieldcount{$filetype}){
                        print STDERR "$csv line $linecount has $#thisrecord fields instead of $fieldcount{$filetype}: [$line]\n";
                    }
                    if($thisrecord[0] =~ /\s/){
                        $transcriptindex++;
                        for(my $i=0;$i<=$#thisrecord;$i++){
                            if($thisrecord[$i] =~ /^\s+$/){
                                $thisrecord[$i] = $thatrecord[$i];
                            }
                        }
                    }
                    else{
                        $transcriptindex = 0;
                    }
                    for(my $i=0;$i<=$#thisrecord;$i++){
                        if($thisrecord[$i] eq "-"){
                            $thisrecord[$i] = "";
                        }
                    }
                    for my $i(@{$percentindices{$filetype}}){
                        $thisrecord[$i] = percent2real($thisrecord[$i]);
                    }
                    for my $i(@{$colonindices{$filetype}}){
                        $thisrecord[$i] =~ s/:/\t/;
                    }
                   my ($gene,$transcript);
                    if($thisrecord[$genefieldindex{$filetype}] =~ /(\S+)\[(\S+)]/){
                        $gene = $1;
                        $transcript = $2;
                    }
                    else{
                        $gene = $thisrecord[$genefieldindex{$filetype}];
                        $transcript = "";
                    }
                    my @ucsc=split /_/, $thisrecord[0];
                    my $ucsc_coord;
                    if($#ucsc <= 3){
                        $ucsc_coord = "chr$ucsc[0]:$ucsc[1]-$ucsc[1]";
                    }
                    else{
                        $ucsc_coord = "chr$ucsc[0]:$ucsc[1]-$ucsc[2]";
                    }
                    printf $filehandle "%s\t%s\t%s\t%s\t%u\t%s\n", ($csv,$ucsc_coord,$gene,$transcript,$transcriptindex,join "\t",@thisrecord);
                    @thatrecord = @thisrecord;
                }
            }
        }
    }
    else{
        print STDERR "undefined or non-file $csv\n";
    }
}

sub percent2real{
    my $this = shift;
    if (defined $this && substr ($this, -1) eq "%"){
        $this = substr($this,0,length($this)-1)/100.00;
    }
    return $this;
}