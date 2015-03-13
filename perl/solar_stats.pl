#!/usr/bin/perl -w
use strict;
#extract statistics summary from solar stats file
my $file = shift;
my $out = "$file.sum";
open IN,"<",$file;
open OUT,">",$out;
my $header = "variable\tn\tmissing\tmean\tmin\tmax\tsd\tskewness\tkurtosis\n";
print OUT $header;
print $header;
my ($variable,$n,$p,$missing,$mean,$min,$max,$sd,$skewness,$kurtosis);
my $in_variable = 0;
my $eor = 0;
while(<IN>){
    if(/Variable:\s*(\S+)\s+Sample Size:\s*([0-9.Ee\\-]+)\s+Missing:\s*([0-9.Ee\\-]+)/){
                            #Variable:  famid  Sample Size: 2520  Missing: 0
        $in_variable = 1;
        ($variable,$n,$missing,$mean,$min,$max,$sd,$skewness,$kurtosis) =
        ($1,$2,$3,"","","","","","");
    }
    elsif($in_variable && /Mean:\s+([0-9.Ee\\-]+)/){
        $mean = $1;
    }
    elsif($in_variable && /Minimum:\s+([0-9.Ee\\-]+)/){
        $min = $1;
    }
    elsif($in_variable && /Maximum:\s+([0-9.Ee\\-]+)/){
        $max = $1;
    }
    elsif($in_variable && /Standard Deviation:\s+([0-9.Ee\\-]+)/){
        $sd = $1;
    }
    elsif($in_variable && /Skewness:\s+([0-9.Ee\\-]+)/){
        $skewness = $1;
    }
    elsif($in_variable && /Kurtosis:\s+([0-9.Ee\\-]+)/){
        $kurtosis = $1;
         $eor = 1;
   }
    if($eor){
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",($variable,$n,$missing,$mean,$min,$max,$sd,$skewness,$kurtosis);
        printf OUT "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",($variable,$n,$missing,$mean,$min,$max,$sd,$skewness,$kurtosis);
        $in_variable = 0;
        $eor = 0;
    }
}

