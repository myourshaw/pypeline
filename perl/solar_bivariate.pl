#!/usr/bin/perl -w
use strict;
use Switch;
use File::Spec::Functions;
use File::Glob;

#extract H2r summary from solar stdout file
my $arg = $ARGV[0];#shift;
my($vol,$dir,$foo) = File::Spec->splitpath($arg);
my $out = File::Spec->catpath($vol,$dir,"solar_bivariate.sum");
open OUT,">",$out;
my $header = "trait1\ttrait2\ttrait1H2r\ttrait1H2rStdErr\ttrait2H2r\ttrait2H2rStdErr\tRhoE\tRhoEp\tRhoEStdErr\tRhoG\tRhoGp0\tRhoGp1\tRhoGStdErr\tRhoP\tRhoPp0\n";
print OUT $header;
print $header;
#my $format = "%s\t%s\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n";
my $format = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n";
my ($trait1,$trait2,$trait1H2r,$trait1H2rStdErr,$trait2H2r,$trait2H2rStdErr,$RhoE,$RhoEp,$RhoEStdErr,$RhoG,$RhoGp0,$RhoGp1,$RhoGStdErr,$RhoP,$RhoPp0);
my $in_summary = 0;
my $startSummaryRx = '\*\s+Summary of Results\s+\*';
my $TraitRx = 'Trait:\s+(\S+)\s+(\S+)\s+Individuals';
my $H2rRx1 = '';
my $H2rStdErrRx1 = '';
my $H2rRx2 = '';
my $H2rStdErrRx2 = '';
my $RhoERx = 'RhoE is\s+([0-9.Ee\\-]+)\s+p [=<]\s+([0-9.Ee\\-]+)';
my $RhoEStdErrRx = 'RhoE Std\. Error:\s+([0-9.Ee\\-]+)';
my $RhoGRx = 'RhoG is\s+([0-9.Ee\\-]+)';
my $RhoGStdErrRx = 'RhoG Std\. Error:\s+([0-9.Ee\\-]+)';
my $RhoGp0Rx = 'RhoG different from zero\s+p [=<]\s+([0-9.Ee\\-]+)';
my $RhoGp1Rx = 'RhoG different from -?1.0\s+p [=<]\s+([0-9.Ee\\-]+)';
my $RhoPRx = 'Derived Estimate of RhoP is\s+([0-9.Ee\\-]+)';
my $RhoPp0Rx = 'RhoP different from zero  p [=<]\s+([0-9.Ee\\-]+)';
my $endRx = 'Final models are named';
my @files = glob($arg);
foreach my $file(glob($arg)){
    open IN,"<",$file;
    while(<IN>){
       if(/$startSummaryRx/){
           $in_summary = 1;
           ($trait1,$trait2,$trait1H2r,$trait1H2rStdErr,$trait2H2r,$trait2H2rStdErr,$RhoE,$RhoEp,$RhoEStdErr,$RhoG,$RhoGp0,$RhoGp1,$RhoGStdErr,$RhoP)
           = ("","","","","","","","","","","","","","");
       }
       elsif($in_summary){
           if (/$TraitRx/) {
               $trait1 = $1;
               $trait2 = $2;
               $H2rRx1 = 'H2r\('.$trait1.'\) is\s+([0-9.Ee\\-]+)';
               $H2rStdErrRx1 = 'H2r\('.$trait1.'\) Std\. Error:\s+([0-9.Ee\\-]+)';
               $H2rRx2 = 'H2r\('.$trait2.'\) is\s+([0-9.Ee\\-]+)';
               $H2rStdErrRx2 = 'H2r\('.$trait2.'\) Std\. Error:\s+([0-9.Ee\\-]+)';
           }
           elsif(/$H2rRx1/){
               $trait1H2r = $1;
           }
           elsif(/$H2rStdErrRx1/){
               $trait1H2rStdErr = $1;
           }
           elsif(/$H2rRx2/){
               $trait2H2r = $1;
           }
           elsif(/$H2rStdErrRx2/){
               $trait2H2rStdErr = $1;
           }
           elsif(/$RhoERx/){
               $RhoE = $1;
               $RhoEp = $2;
           }
           elsif(/$RhoEStdErrRx/){
               $RhoEStdErr = $1;
           }
           elsif(/$RhoGRx/){
               $RhoG = $1;
           }
           elsif(/$RhoGStdErrRx/){
               $RhoGStdErr = $1;
           }
            elsif(/$RhoGp0Rx/){
               $RhoGp0 = $1;
           }
           elsif(/$RhoGp1Rx/){
               $RhoGp1 = $1;
           }
          elsif(/$RhoPRx/){
               $RhoP = $1;
          }
          elsif(/$RhoPp0Rx/){
               $RhoPp0 = $1;
          }
          elsif(/$endRx/){
               printf $format,($trait1,$trait2,$trait1H2r,$trait1H2rStdErr,$trait2H2r,$trait2H2rStdErr,$RhoE,$RhoEp,$RhoEStdErr,$RhoG,$RhoGp0,$RhoGp1,$RhoGStdErr,$RhoP,$RhoPp0);
               printf OUT $format,($trait1,$trait2,$trait1H2r,$trait1H2rStdErr,$trait2H2r,$trait2H2rStdErr,$RhoE,$RhoEp,$RhoEStdErr,$RhoG,$RhoGp0,$RhoGp1,$RhoGStdErr,$RhoP,$RhoPp0);
               $in_summary = 0;
           }
       }
   }
}


