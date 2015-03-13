#!/usr/bin/perl -w
# ©2007 Michael Yourshaw All Rights Reserved
use strict;
use warnings;
use Statistics::Descriptive;
use Bio::SeqFeature::Primer;
my %params = (
    project_id => "",
	input_format => "", #
	output_format => "", #
	forward_primer_id => "OLSF", #OLSF
	forward_primer_site => "GGCTCGCTACAGAATCAGTT", # = forward primer OLSF
	reverse_primer_id => "OLSR", #OLSR
	reverse_primer_site => "AAGGACTGGATACTCTCTGC", # = reverse complement of reverse primer OLSR
	Na_M_conc => 0.050, # [Na+] M
	oligo_M_conc => 0.00000025, # [oligo] M
    percent_formamide => 0 # percent [formamide]
);
my $input_file = "";
my $output_file = "";
my $probe_file = "";
my $coverage_file = "";
my @all_probes;
my @GCsProbe;
my @GCsOligo;
my @CenterToCenters;
#my @Tms;
#my @TmsBio;
#my @Tms_no_Na;
my @TmsProbe;
my @TmsOligo;
my @uncovereds;
my $covered = 1;
my $uncovered_regions = 0;
my %coverage;
use constant CHR_SHIFT => 1000000000;
my $usage = "usage: perl tm.pl Tiling_tdt file list [> Tm_file]\n";
unless (@ARGV) {
	print $usage;
	exit;
}
#redirect STDOUT
my $gt = "&gt;"; #komodo saves command line > as &gt;
if(index($ARGV[$#ARGV],">") == 0){
    $output_file = substr pop(@ARGV),1;
}
elsif ($ARGV[$#ARGV-1] eq ">" || $ARGV[$#ARGV-1] eq $gt) {
    $output_file = pop(@ARGV);
    pop(@ARGV);
}
elsif (rindex($ARGV[$#ARGV-1],">") == length($ARGV[$#ARGV-1])-1) {
    $output_file = pop(@ARGV);
    $ARGV[$#ARGV-1] = substr($ARGV[$#ARGV-1],0,rindex($ARGV[$#ARGV-1],">"));
}
elsif (rindex($ARGV[$#ARGV-1],$gt) == length($ARGV[$#ARGV-1])-1) {
    $output_file = pop(@ARGV);
    $ARGV[$#ARGV-1] = substr($ARGV[$#ARGV-1],0,rindex($ARGV[$#ARGV-1],$gt));
}
if($ARGV[$#ARGV] eq ">" || $ARGV[$#ARGV] eq $gt){
    pop(@ARGV);
}
if($output_file){
    open STDOUT, ">", $output_file or die "can't open $output_file: $!";
}

#read from a list of filenames on the command line
#read from SDTIN if no filenames were given
#"-" indicates STDIN
#"someprogram |" indicates the output of another program (7.14)
#read and parse input file
while (<>) {
	if ($ARGV ne $input_file) {
        print_probes();
		print STDERR sprintf("input file: %s\n", File::Spec->rel2abs($ARGV));
		$input_file = $ARGV;
	}
	chomp;
	#comment
	if (/^\s*#/) { }
	#parameter
	elsif (/^\s*\$/) {
		parameter();
	}
	#TargetID	BPStart	Sequence	ProbeLength	EndDistance
    #chrX:054085769-054085998	16	CTTACCTGCCATGAAACCAGTCCTGGCACATGTCACACTCGATCA	45	214
    elsif(/^\s*(?:chr)?[\s]*([0-9XYMT]{1,2})[\s:]+([,\d]+)\D+([,\d]+)\s+(\d+)\s+([A-Z]+)\s+(\d+)\s+(\d+)/i){
        my $chr = $1;
        my $target_start = $2;
        my $target_end = $3;
		my $this_TargetID = "chr$chr:$target_start-$target_end";
		my $this_BPStart = $4;
        my $this_Sequence = uc($5);
		my $this_ProbeLength = $6;
		my $this_EndDistance = $7;
		my $num_chr = chr2numchr($chr);
        my $start = $target_start + $this_BPStart - 1;
        my $end = $start + $this_ProbeLength - 1;
		push @all_probes, {
            chr=>$chr,
            target_start=>$target_start,
            target_end=>$target_end,
            TargetID=>$this_TargetID,
            BPStart=>$this_BPStart,
            Sequence=>$this_Sequence,
            ProbeLength=>$this_ProbeLength,
            EndDistance=>$this_EndDistance,
            num_chr=>$num_chr,
            start=>$start,
            end=>$end};
    }
}
print_probes();
printf STDERR "probe file: %s\n",$probe_file;
printf STDERR "Tm file: %s\n",$output_file;
print STDERR "DONE\n";

sub print_probes{
    my $that_TargetID = "";
    my $that_BPStart = 0;
    my $that_ProbeLength = 0;
    my $that_start = 0;
    my $center_to_center = 0;
    my @probes;
    if(@all_probes){
        unless($output_file){
            $output_file = File::Spec->rel2abs($input_file) . "_Tm.txt";
            open STDOUT, ">", $output_file or die "can't open $output_file: $!";
        }
        $probe_file = File::Spec->rel2abs($input_file) . "_probes.txt";
        open PROBEFILE, ">", $probe_file or die "can't open $probe_file: $!";
        $coverage_file = File::Spec->rel2abs($input_file) . "_coverage.txt";
        open COVERAGEFILE, ">", $coverage_file or die "can't open $coverage_file: $!";
        print "TargetID\tProbeID\tChromosome\tStart\tEnd\tBPStart\tSequence\tProbeLength\tEndDistance\tUncovered\tTmProbe\tPercentGCProbe\tOligo\tTmOligo\tPercentGCOligo\tCenterToCenter\n";
        @all_probes = sort {$a->{num_chr} <=> $b->{num_chr} || $a->{TargetID} cmp $b->{TargetID} || $a->{BPStart} <=> $b->{BPStart}} @all_probes;
        foreach my $probe(@all_probes) {
            for(my $i=$probe->{start};$i<=$probe->{end};$i++){
                $coverage{$i+CHR_SHIFT*$probe->{num_chr}}++;
            }
            if($probe->{TargetID} ne $that_TargetID){
                $covered = 1;
                $center_to_center = 0;
            }
            else{
                $center_to_center = $probe->{start} + ($probe->{ProbeLength})/2 - $that_start - ($that_ProbeLength)/2;
            }
            my $uncovered = $probe->{TargetID} eq $that_TargetID ?
                $probe->{BPStart}>($that_BPStart+$that_ProbeLength) ?
                    $probe->{BPStart}-($that_BPStart+$that_ProbeLength) :
                    0 :
                $probe->{BPStart}==1 ?
                    0 :
                    ($probe->{BPStart}-1);
            #my $this_center = $probe->{start}+($probe->{ProbeLength})/2;
            #my $that_center = $that_start+($that_ProbeLength)/2;
            if($uncovered){
                push @uncovereds,$uncovered;
                if($covered){
                    $covered = 0;
                    $uncovered_regions++;
                }
            }
            $probe->{uncovered} = $uncovered;
            $that_TargetID = $probe->{TargetID};
            $that_BPStart = $probe->{BPStart};
            $that_ProbeLength = $probe->{ProbeLength};
            $that_start = $probe->{start};
            my $TmProbe = Tm($probe->{Sequence});
            my $GCProbe = GC($probe->{Sequence});
            push @TmsProbe,$TmProbe;
            push @GCsProbe,$GCProbe;
            my $oligo = $params{forward_primer_site}.$probe->{Sequence}.$params{reverse_primer_site};
            my $TmOligo = Tm($oligo);
            my $GCOligo = GC($oligo);
            push @TmsOligo,$TmOligo;
            push @GCsOligo,$GCOligo;
            push @CenterToCenters,$center_to_center unless $center_to_center == 0;
            my $probe_id = sprintf "%schr%s_%s_%s_%s",$params{project_id}?"$params{project_id}_":"",numchr2chr($probe->{num_chr}),$probe->{target_start},$probe->{target_end},$probe->{BPStart};
            printf "%s\t%s\t%s\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%.1f\t%.0f\t%s\t%.1f\t%.0f\t%.1f\n",
                $probe->{TargetID},
                $probe_id,
                numchr2chr($probe->{num_chr}),
                $probe->{start},
                $probe->{end},
                $probe->{BPStart},
                $probe->{Sequence},
                $probe->{ProbeLength},
                $probe->{EndDistance},
                $uncovered,
                $TmProbe,
                $GCProbe*100,
                $oligo,
                $TmOligo,
                $GCOligo*100,
                $center_to_center;
            printf PROBEFILE "%s\t%s\n",$probe_id,$oligo;
        }
#        while((my $p,my $n)=each(%coverage)){
#            print STDERR "$p\t$n\n"
#        }
#        my $junk = $coverage{"101261365"};
#        print STDERR $junk;
        foreach my $position(sort {$a<=>$b} keys %coverage){
            my $c = numchr2chr(int($position/CHR_SHIFT));
            my $p = $position % CHR_SHIFT;
            my $n = $coverage{$position};
            printf COVERAGEFILE "%s\t%d\t%d\n",$c,$p,$n;
        }
        print_stats();
    }
    undef @all_probes;
    undef @GCsProbe;
    undef @GCsOligo;
    #undef @Tms;
    #undef@TmsBio;
    #undef @Tms_no_Na;
    undef @TmsProbe;
    undef @TmsOligo;
    undef @uncovereds;
}

sub print_stats{
    print "ANALYSIS\n";
    printf "[Na+] (mM)\t%.0f\n", $params{Na_M_conc}*1000;
    printf "[oligo] (nM)\t%.0f\n", $params{oligo_M_conc}*1000000000;
    printf "[formamine] (%%)\t%.1f\n",$params{percent_formamide};
    printf "Forward primer site (%s)\t%s\tlength\t%d\tTm (°C)\t%.1f\tGC%%\t%.0f\n",$params{forward_primer_id},$params{forward_primer_site},length($params{forward_primer_site}),Tm($params{forward_primer_site}),GC($params{forward_primer_site})*100;
    printf "Reverse primer site (%s)\t%s\tlength\t%d\tTm (°C)\t%.1f\tGC%%\t%.0f\n",$params{reverse_primer_id},$params{reverse_primer_site},length($params{reverse_primer_site}),Tm($params{reverse_primer_site}),GC($params{reverse_primer_site})*100;
    my $stat = Statistics::Descriptive::Full->new();
    $stat->add_data(@TmsProbe);
    print "Probes\n";
    printf "Tm mean (°C)\t%.1f\n", $stat->mean();
    printf "Tm min (°C)\t%.1f\n", $stat->min();
    printf "Tm max (°C)\t%.1f\n", $stat->max();
    printf "Tm range (max-min) (°C)\t%.1f\n", $stat->sample_range();
    printf "Tm standard deviation (°C)\t%.1f\n", $stat->standard_deviation();
    #my %f=$stat->frequency_distribution(10);
    $stat = Statistics::Descriptive::Full->new();
    $stat->add_data(@GCsProbe);
    printf "GC%% min\t%.1f\n", $stat->min()*100;
    printf "GC%% max\t%.1f\n", $stat->max()*100;
    $stat = Statistics::Descriptive::Full->new();
    $stat->add_data(@TmsOligo);
    print "Oligos\n";
    printf "Tm mean (°C)\t%.1f\n", $stat->mean();
    printf "Tm min (°C)\t%.1f\n", $stat->min();
    printf "Tm max (°C)\t%.1f\n", $stat->max();
    printf "Tm range (max-min) (°C)\t%.1f\n", $stat->sample_range();
    printf "Tm standard deviation (°C)\t%.1f\n", $stat->standard_deviation();
    #my %f=$stat->frequency_distribution(10);
    $stat = Statistics::Descriptive::Full->new();
    $stat->add_data(@GCsOligo);
    printf "GC%% min\t%.1f\n", $stat->min()*100;
    printf "GC%% max\t%.1f\n", $stat->max()*100;
    $stat = Statistics::Descriptive::Full->new();
    $stat->add_data(@CenterToCenters);
    printf "mean probe center-to-center\t%.1f\n",$stat->mean();
    printf "min probe center-to-center\t%.1f\n",$stat->min();
    printf "max probe center-to-center\t%.1f\n",$stat->max();
    printf "probe center-to-center SD\t%.1f\n",$stat->standard_deviation();
    printf "regions with at least one uncovered base\t%d\n",$uncovered_regions;
    if(@uncovereds){
        $stat = Statistics::Descriptive::Full->new();
        $stat->add_data(@uncovereds);
        printf "total uncovered bases\t%d\n",$stat->sum();
        printf "longest uncovered interval\t%d\n", $stat->max();
    }
}
sub log10 {
	return log(shift)/log(10);
}

sub Tm{
    my $seq = shift;
    my $l = length $seq;
    my $A = ($seq =~ tr/A//);
    my $C = ($seq =~ tr/C//);
    my $G = ($seq =~ tr/G//);
    my $T = ($seq =~ tr/T//);
    #Tm calculations from
    #http://www.basic.northwestern.edu/biotools/oligocalc.html
#     my $Tm_no_Na = $l<14 ?
#        2*($A+$T)+4*($G+$C) :
#        64.9+41*($G+$C-16.4)/($l);
#    push @Tms_no_Na,$Tm_no_Na;
    return $l<14 ?
        2*($A+$T)+4*($G+$C)-16.6*log10(0.050)+16.6*log10($params{Na_M_conc}) : $l<50 ?
            100.5+(41*($G+$C)/($l))-(820/($l))+16.6*log10($params{Na_M_conc}) :
            81.5+(41*($G+$C)/($l))-(500/($l))+16.6*log10($params{Na_M_conc})-0.62*$params{percent_formamide};
    #my $TmBio = Bio::SeqFeature::Primer->new(-seq=>$probe->{Sequence})->Tm(-salt=>$params{Na_M_conc}, -oligo=>$params{oligo_M_conc});
}

sub GC{
    my $seq = shift;
    my $C = ($seq =~ tr/C//);
    my $G = ($seq =~ tr/G//);
    return ($C+$G)/length $seq;
}

#
sub chr2numchr{
    my $chr = shift;
    return $chr eq "X" ? 23 : $chr eq "Y" ? 24 : ($chr eq "M" || $chr eq "MT") ? 25 : $chr;
}

sub numchr2chr{
    my $chr = shift;
    return $chr == 23 ? "X" : $chr == 24 ? "Y" : $chr == 25 ? "M" : $chr;
}

#change parameter ($key value)
sub parameter {
	if(/s*\$([^\s=>#]+){1}[\s=>#]+([^#]+)#*$/){;
		my $key = $1;
		my $value = $2;
		$value =~ s/['"]+//g;
		if(defined($value)){
			if ($value =~ /^false$/i || $value =~ /^no$/i) {
				$value = 0;
			}
			elsif ($value =~ /^true$/i || $value =~ /^yes$/i) {
				$value = 1;
			}
		}
		if (exists($params{$key})) {
			$params{$key} = $value;
		}
		else {
			printf STDERR "%s isn't a setting", $key;
		}
	}
} # !parameter
