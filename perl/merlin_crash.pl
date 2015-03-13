#!/usr/bin/perl -w
use strict;
use warnings;

use File::Path;
use File::Spec;

#find missing partitions and merlin crashes/fatal errors

my $dir = "/Volumes/raid/merlin_files/merlin_2009-10-22";

my ($nominal_threshhold, $suggestive_threshhold, $interesting_threshhold, $significant_threshhold, $highly_significant_threshhold, $bizarre, $exclude) = (1.0, 2.2, 3.0, 3.6, 5.4, 9.999, -2);

my %codes = (
    all  => "3",
    phase1  => "1",
    phase2  => "2",
    freq => "f",
    ibd => "i",
    lod => "l",
    out => "o",
    pdf => "p",
    table => "t",
    zscore => "z",
    all_ethnic  => "a",
    non_white  => "n",
    white  => "w",
    all_gender  => "a",
    female_family  => "f",
    female_sibs  => "g",
    male_family  => "m",
    male_sibs  => "n",
    broad  => "b",
    narrow  => "n",
    all_affected_family  => "b",
    all_type  => "a",
    combined  => "c",
    combined_family  => "d",
    combined_inattentive  => "e",
    combined_inattentive_family  => "f",
    hyperactive  => "h",
    hyperactive_family  => "k",
    inattentive  => "i",
    inattentive_family  => "j"
);
my %extensions = (
    freq => "freq",
    ibd => "ibd",
    lod => "lod",
    out => "out",
    pdf => "pdf",
    table => "tbl",
    zscore => "zscore"
);
my @phases = qw(all phase1 phase2);
my @folders = qw(out);#qw(freq ibd lod out pdf table zscore);
my @ethnics = qw(all_ethnic non_white white);
my @genders = qw(all_gender female_family female_sibs male_family male_sibs);
my @adhd_lifes = qw(broad narrow);
my @adhd_types = qw(all_affected_family all_type combined combined_family combined_inattentive combined_inattentive_family hyperactive hyperactive_family inattentive inattentive_family);
my @chromosomes = ("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22","X","XY","Y");

my %nominals = ();
my %suggestives = ();
my %interestings = ();
my %significants = ();
my %highly_significants = ();

#my $lin_file = File::Spec->catfile($dir,"merlin-linear.tbl");
#open LIN,">",$lin_file;
#my $exp_file = File::Spec->catfile($dir,"merlin-exponential.tbl");
#open EXP,">",$exp_file;

#my $lin_suggestive_file = File::Spec->catfile($dir,"merlin-lin_suggestive.tbl");
#open SUG,">",$lin_suggestive_file;
##my $lin_significant_file = File::Spec->catfile($dir,"merlin--lin_significant.tbl");
##open SIG,">",$lin_significant_file;
##my $lin_highly_significant_file = File::Spec->catfile($dir,"merlin-lin_highly_significant.tbl");
##open HI,">",$lin_highly_significant_file;

#my $exp_suggestive_file = File::Spec->catfile($dir,"merlin-exp_suggestive.tbl");
#open EXSUG,">",$exp_suggestive_file;
#my $exp_significant_file = File::Spec->catfile($dir,"merlin-exp_significant.tbl");
#open EXSIG,">",$exp_significant_file;
#my $exp_highly_significant_file = File::Spec->catfile($dir,"merlin-exp_highly_significant.tbl");
#open EXHI,">",$exp_highly_significant_file;
#
#my $bizarre_file = File::Spec->catfile($dir,"merlin-bizarre.tbl");
#open BIZ,">",$bizarre_file;
#
#my $exclude_file = File::Spec->catfile($dir,"merlin-exclude.tbl");
#open EX,">",$exclude_file;

##my $marker_summary_file = File::Spec->catfile($dir,"merlin-marker_summary.tbl");
##open MARK,">",$marker_summary_file;

##my $lin_header = "ETHNIC\tGENDER\tADHD_LIFE\tADHD_TYPE\tCHR\tPOS\tMARKER\tVARIABLE\tANALYSIS\tZSCORE\tDELTA\tLOD\tPVALUE";
#my $exp_header = "ETHNIC\tGENDER\tADHD_LIFE\tADHD_TYPE\tCHR\tPOS\tMARKER\tVARIABLE\tANALYSIS\tZSCORE\tExDELTA\tExLOD\tExPVALUE";
#print LIN "$lin_header\n";
#print EXP "$exp_header\n";
#print SUG "THRESHHOLD\t$lin_header\n";
##print SIG "THRESHHOLD\t$lin_header\n";
##print HI "THRESHHOLD\t$lin_header\n";
#print EXSUG "THRESHHOLD\t$exp_header\n";
#print EXSIG "THRESHHOLD\t$exp_header\n";
#print EXHI "THRESHHOLD\t$exp_header\n";

my $missing = 0;
my $error = 0;
foreach my $phase(@phases){
    foreach my $folder(@folders){
    foreach my $chr(@chromosomes){
        foreach my $ethnic(@ethnics){
            foreach my $gender(@genders){
                foreach my $adhd_life(@adhd_lifes){
                    foreach my $adhd_type(@adhd_types){
                        my $foo = "$phase\t$folder\t$chr\t$ethnic\t$gender\t$adhd_life\t$adhd_type";
                        my $partition = "$codes{$phase}\t$codes{$folder}\t$chr\t$codes{$ethnic}\t$codes{$gender}\t$codes{$adhd_life}\t$codes{$adhd_type}";
                        (my $partition_label = $partition) =~ s/\t/-/g;
                        #print "$partition_label\n";
                        my $in_file = File::Spec->catfile($dir,"merlin_results-$phase",$folder,"merlin-$phase-chr$chr-$ethnic-$gender-$adhd_life-$adhd_type.$extensions{$folder}");
                        if (!-f $in_file){
                            print "       missing $foo\n";
                            $missing++;
                            #print STDERR "$in_file\n";
                        }
                        else{
                            open IN,"<",$in_file or die "cant't open $in_file $!\n";
                            while(<IN>)
                            {
                                if(/CRASHED/){
                                    print "CRASHED\t$foo\n";
                                    $error++;
                                    next;
                                }
                                elsif(/FATAL ERROR/){
                                    print "ERROR  \t$foo\n";
                                    $error++;
                                    next;
                                }

                            }
                        }
                        #else{
                        #    open IN,"<",$in_file or die "cant't open $in_file $!\n";
                        #    while(<IN>){
                        #        my $line = $_;
                        #        chomp;
                        #        if(/\s*([0-9XY]{1,2})\s+([e\d.-]+)\s+(\S+)\s+(\S+)\s+\[(\S+)\]\s+([e\d.-]+)\s+([e\d.-]+)\s+([e\d.-]+)\s+([e\d.-]+)\s+([e\d.-]+)\s+([e\d.-]+)\s+([e\d.-]+)/){
                        #            my ($chr,$pos,$marker,$variable,$analysis,$zscore,$delta,$lod,$pvalue,$exdelta,$exlod,$expvalue) = ($1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12);
                        #            my @lin_values = ($partition,$chr,$pos,$marker,$variable,$analysis,$zscore,$delta,$lod,$pvalue);
                        #            #my @exp_values = ($partition,$chr,$pos,$marker,$variable,$analysis,$zscore,$exdelta,$exlod,$expvalue);
                        #            #my $all_format = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n";
                        #            my $summary_format = ">=%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n";
                        #            
                        #            #printf LIN $all_format, @lin_values;
                        #            #printf EXP $all_format, @exp_values;
                        #            if($pvalue != 0.5){
                        #                #if($lod >= $suggestive_threshhold && $lod < $bizarre){
                        #                #    printf SUG $summary_format, $suggestive_threshhold,@lin_values;
                        #                #    push (@{$lin_suggestives{"$chr:$marker"}}, $partition_label);
                        #                #}
                        #                if($lod >= $significant_threshhold && $lod < $bizarre){
                        #                    printf SIG $summary_format, $significant_threshhold, @lin_values;
                        #                    push (@{$lin_significants{"$chr:$marker"}}, $partition_label);
                        #                }
                        #                if($lod >= $highly_significant_threshhold && $lod < $bizarre){
                        #                    printf HI $summary_format, $highly_significant_threshhold,@lin_values;
                        #                    push (@{$lin_highly_significants{"$chr:$marker"}}, $partition_label);
                        #                }
                        #            }
                        #            #if($expvalue != 0.5){
                        #            #    if($exlod >= $suggestive_threshhold && $exlod < $bizarre){
                        #            #        printf EXSUG $summary_format, $suggestive_threshhold,@exp_values;
                        #            #        push (@{$exp_suggestives{"$chr:$marker"}}, $partition_label);
                        #            #    }
                        #            #    if($exlod >= $significant_threshhold && $exlod < $bizarre){
                        #            #        printf EXSIG $summary_format, $significant_threshhold,@exp_values;
                        #            #        push (@{$exp_significants{"$chr:$marker"}}, $partition_label);
                        #            #    }
                        #            #    if($exlod >= $highly_significant_threshhold && $exlod < $bizarre){
                        #            #        printf EXHI $summary_format, $highly_significant_threshhold,@exp_values;
                        #            #        push (@{$exp_highly_significants{"$chr:$marker"}}, $partition_label);
                        #            #    }
                        #            #    if(($lod <= $exclude && $lod >= $bizarre) || ($exlod <= $exclude && $exlod >= $bizarre)){
                        #            #        print EX "$partition\t$line";
                        #            #    }
                        #            #    if($lod >= $bizarre || $exlod >= $bizarre || $lod <= -$bizarre || $exlod <= -$bizarre){
                        #            #        print BIZ "$partition\t$line";
                        #            #    }
                        #            #}
                        #        }
                        #    }
                        #}
                    }
                }
            }
        }
    }
    }
}

##print MARK map {"lin_suggestive\t$_\t@{$lin_suggestives{$_}}\n"} sort keys %lin_suggestives;
#print MARK map {"lin_significant\t$_\t@{$lin_significants{$_}}\n"} sort keys %lin_significants;
#print MARK map {"lin_highly_significant\t$_\t@{$lin_highly_significants{$_}}\n"} sort keys %lin_highly_significants;
##print MARK map {"exp_suggestive\t$_\t@{$exp_suggestives{$_}}\n"} sort keys %exp_suggestives;
##print MARK map {"exp_significant\t$_\t@{$exp_significants{$_}}\n"} sort keys %exp_significants;
##print MARK map {"exp_highly_significant\t$_\t@{$exp_highly_significants{$_}}\n"} sort keys %exp_highly_significants;

print "done\n$missing missing partitions\n$error crashes or fatal errors\n";
exit;
