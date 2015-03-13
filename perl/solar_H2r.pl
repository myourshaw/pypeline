#!/usr/bin/perl -w
use strict;
#extract H2r summary from solar stdout file
my $file = shift;
my $out = "$file.sum";
my $cov = "$file.covariates";
open IN,"<",$file;
open OUT,">",$out;
open COV,">",$cov;
my $header = "trait\tH2r\tH2rp\tsignificant_H2r\tH2r_std_error\tC2\tC2p\tsignificant_C2\tC2_std_error\tcovariate_proportion\tresidual_kurtosis\tnormal_kurtosis\tk-l_r2\tcovariates\n";
print OUT $header;
print $header;
my ($trait,$H2r,$H2rp,$significant_H2r,$H2r_std_error,$C2,$C2p,$significant_C2,$C2_std_error,$covariate_proportion,$residual_kurtosis,$normal_kurtosis,$kullback);
my @covariates = ();
my $in_summary = 0;
my $proportion_next = 0;
my $eor = 0;
while(<IN>){
    if(/\*\s+Summary of Results\s+\*/){
        $in_summary = 1;
        ($trait,$H2r,$H2rp,$significant_H2r,$H2r_std_error,$C2,$C2p,$significant_C2,$C2_std_error,$covariate_proportion,$residual_kurtosis,$normal_kurtosis,$kullback) =
        ("","","","","","","","","","","","","");
        @covariates = ();
    }
    elsif($in_summary && /Trait:\s+(\S+)\s+Individuals:/){
        $trait = $1;
    }
    elsif($in_summary && /H2r is\s([0-9.Ee\\-]+)\s+p =\s+([0-9.Ee\\-]+)\s+(.+)/){
        ($H2r,$H2rp) = ($1,$2);
        $significant_H2r = $3;# eq "Significant" ? 1 : 0;
    }
    elsif($in_summary && /H2r Std. Error:\s+([0-9.Ee\\-]+)/){
        $H2r_std_error = $1;
    }
    elsif($in_summary && /C2 is\s([0-9.Ee\\-]+)\s+p =\s+([0-9.Ee\\-]+)\s+(.+)/){
        ($C2,$C2p) = ($1,$2);
        $significant_C2 = $3;# eq "Significant" ? 1 : 0;
    }
    elsif($in_summary && /C2 Std. Error:\s+([0-9.Ee\\-]+)/){
        $C2_std_error = $1;
    }
    elsif($in_summary && /(\S+)\s+p = [0-9.Ee\\-]+\s+\(Significant\)/){#/){
        push @covariates,$1;
    }
    elsif($in_summary && /Proportion of Variance Due to All Final Covariates Is/){
        $proportion_next = 1;
    }
    elsif($in_summary && $proportion_next && /([0-9.Ee\\-]+)/){
        $covariate_proportion = $1;
        $proportion_next = 0;
    }
    elsif($in_summary && /Kullback-Leibler R-squared is ([0-9.Ee\\-]+)/){
        $kullback = $1;
        $eor = 1;
    }
    elsif($in_summary && /Residual Kurtosis is\s+([0-9.Ee\\-]+),?\s+(.+)/){
        ($residual_kurtosis,$normal_kurtosis) = ($1,$2);
        if(!defined($H2r_std_error)){
            $H2r_std_error = "";
        }
        $eor = 1;
    }
    elsif($in_summary && /An error occurred while computing residual kurtosis/){
        ($residual_kurtosis,$normal_kurtosis) = ("error","error");
        if(!defined($H2r_std_error)){
            $H2r_std_error = "";
        }
        $eor = 1;
    }
    if($eor){
        my $line = sprintf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$trait,$H2r,$H2rp,$significant_H2r,$H2r_std_error,$C2,$C2p,$significant_C2,$C2_std_error,$covariate_proportion,$residual_kurtosis,$normal_kurtosis,$kullback,join("\t",map {"covariate $trait($_)"} @covariates));
        print $line;
        print OUT $line;
        print COV map {"covariate $trait($_)\n"} @covariates;
        $in_summary = 0;
        $eor = 0;
    }
}

