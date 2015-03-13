#!/usr/bin/perl -w
use strict;
use File::Spec::Functions;

# uses the output from solar_bivariate.pl to produce square
# matrices of each statistic, and of RhoG/RhoE, RhoG/RhoP

my $infile = shift;
open IN,"<",$infile;
my($vol,$dir,$foo) = File::Spec->splitpath($infile);
my $header = "trait1\ttrait2\ttrait1H2r\ttrait1H2rStdErr\ttrait2H2r\ttrait2H2rStdErr\tRhoE\tRhoEp\tRhoEStdErr\tRhoG\tRhoGp0\tRhoGp1\tRhoGStdErr\tRhoP\tRhoPp0";
my @fields = split /\s/,$header,-1;
my $indata = 0;;
my $linenumber = 0;
my %traits = ();
my %RhoE = ();
my %RhoEp = ();
my %RhoG = ();
my %RhoGp0 = ();
my %RhoGp1 = ();
my %RhoP = ();
my %RhoPp0 = ();
while (<IN>){
    chomp;
    $linenumber++;
    if(length($_) > 0){
        if($indata){
            my @data = split /\s/,$_,-1;
            $#fields == $#data or die "line $linenumber has $#data elements instead of $#fields";
            my ($trait1,$trait2,$trait1H2r,$trait1H2rStdErr,$trait2H2r,$trait2H2rStdErr,$RhoE,$RhoEp,$RhoEStdErr,$RhoG,$RhoGp0,$RhoGp1,$RhoGStdErr,$RhoP,$RhoPp0) = @data;
            $traits{$trait1}++;
            $traits{$trait2}++;
            $RhoE{"$trait1.$trait2"} = $RhoE;
            $RhoEp{"$trait1.$trait2"} = $RhoEp;
            $RhoG{"$trait1.$trait2"} = $RhoG;
            $RhoGp0{"$trait1.$trait2"} = $RhoGp0;
            $RhoGp1{"$trait1.$trait2"} = $RhoGp1;
            $RhoP{"$trait1.$trait2"} = $RhoP;
            $RhoPp0{"$trait1.$trait2"} = $RhoPp0;
        }
        elsif($_ eq $header){
            $indata = 1;
        }
    }
}
my @stats = qw(RhoE RhoEp RhoG RhoGp0 RhoGp1 RhoP RhoPp0);
my @selfs = (1,0,1,0,0,1,0);
my $thisstat = 0;
my @traits = sort keys %traits;
my @datasets = ((\%RhoE,\%RhoEp,\%RhoG,\%RhoGp0,\%RhoGp1,\%RhoP,\%RhoPp0));
my $line;
my $outfile;
foreach my $stat(@datasets){
    $outfile = File::Spec->catpath($vol,$dir,"solar_bivariate.$stats[$thisstat]");
    open OUT,">",$outfile;
    $line = "\t".join ("\t",@traits)."\n";
    print OUT $line;
    print "$stats[$thisstat]$line";
    foreach my $t1(@traits){ #row
        my @re = ($t1);
        foreach my $t2(@traits){ #columns
            if(exists($stat->{"$t1.$t2"})){
                push @re,$stat->{"$t1.$t2"};
            }
            elsif(exists($stat->{"$t2.$t1"})){
                push @re,$stat->{"$t2.$t1"};
            }
            elsif($t1 eq $t2){
                push @re, $selfs[$thisstat];
            }
            else{
                push @re, "";
            }
        }
        $line = join("\t", @re)."\n";
        print OUT $line;
        print $line;
    }
    $thisstat++;
    print "\n";
    close OUT;
}

# RhoG bottom left, RhoE top right
$outfile = File::Spec->catpath($vol,$dir,"solar_bivariate.$stats[2].$stats[0]");
open OUT,">",$outfile;
$line = "\t".join ("\t",@traits)."\n";
print OUT $line;
print "$stats[2]/$stats[0]$line";
foreach my $t1(@traits){ #row
    my @re = ($t1);
    foreach my $t2(@traits){ #columns
        my $value = "";
        if($t1 gt $t2){
            if(exists($datasets[2]->{"$t1.$t2"})){
                push @re,$datasets[2]->{"$t1.$t2"};
            }
            elsif(exists($datasets[2]->{"$t2.$t1"})){
                push @re,$datasets[2]->{"$t2.$t1"};
            }
        }
        elsif($t2 gt $t1){
            if(exists($datasets[0]->{"$t1.$t2"})){
                push @re,$datasets[0]->{"$t1.$t2"};
            }
            elsif(exists($datasets[0]->{"$t2.$t1"})){
                push @re,$datasets[0]->{"$t2.$t1"};
            }

        }
        else{
            push @re, "";                
        }
    }
        $line = join("\t", @re)."\n";
        print OUT $line;
        print $line;
}
print "\n";
close OUT;

# RhoG bottom left, RhoP top right
$outfile = File::Spec->catpath($vol,$dir,"solar_bivariate.$stats[2].$stats[5]");
open OUT,">",$outfile;
$line = "\t".join ("\t",@traits)."\n";
print OUT $line;
print "$stats[2]/$stats[5]$line";
foreach my $t1(@traits){ #row
    my @re = ($t1);
    foreach my $t2(@traits){ #columns
        my $value = "";
        if($t1 gt $t2){
            if(exists($datasets[2]->{"$t1.$t2"})){
                push @re,$datasets[2]->{"$t1.$t2"};
            }
            elsif(exists($datasets[2]->{"$t2.$t1"})){
                push @re,$datasets[2]->{"$t2.$t1"};
            }
        }
        elsif($t2 gt $t1){
            if(exists($datasets[5]->{"$t1.$t2"})){
                push @re,$datasets[5]->{"$t1.$t2"};
            }
            elsif(exists($datasets[5]->{"$t2.$t1"})){
                push @re,$datasets[5]->{"$t2.$t1"};
            }

        }
        else{
            push @re, "";                
        }
    }
        $line = join("\t", @re)."\n";
        print OUT $line;
        print $line;
}
close OUT;


