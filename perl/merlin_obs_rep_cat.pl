#!perl -w
use strict;
use warnings;

#concatenate merlin tables

use File::Path;
use File::Spec;
use File::Glob;
use Switch;

my $do_npl = 0;
my $do_vc = 0;
my $do_regress = 1;
my @seeds = qw(00000000 00123456);
#my @seeds = qw(00000000);
#my @datasets = qw(00123456);
my $create_files = 1;

my $dir = shift; #"/Volumes/raid/link/merlin";
print "usage: perl cat_merlin_tbl <top level dir>" unless $dir;

my %codes = (
    phase1  => "1", #phase1 markers, phase1 sample and part phase2 sample
    phase2  => "2", #phase2 markers, part phase1 and all phase2 sample
    "phase1+2"  => "3", #phase1 and phase2 markers, phase1 and phase2 sample (combined)
    "phase2-1" => "4", #phase2 markers, no phase1 and all phase2 sample
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
    inattentive_family  => "j",
    ALL => "a",
    Pairs => "p",
    QTL => "q"
);

foreach my $dataset(@seeds){
    my $catdir = File::Spec->catfile($dir,"cat-$dataset");
    File::Path->make_path($catdir) unless !$create_files;  
    my ($header,$file_regex,$line_regex,$out_format,$out_file,$glob);
    my @outfiles = ();
    
    if($do_npl){
        my $npl_all_record_count = 0;
        my $npl_pairs_record_count = 0;
        my $qtl_record_count = 0;
        my %qtls = ();
        my %npl_chr_seed = ();
        $file_regex = '\/(?:chr)?([0-9XY]{1,2}|(?:autosome)|x)-(\d{8})-(\d{5})-nonparametric\.tbl';
        #$line_regex = '\s*([0-9XY]{1,3})\s+([e\d.-]+)\s+(\S+)\s+(\S+)\s+\[(\S+)\]\s+([e\d.-]+)\s+([e\d.-]+)\s+([e\d.-]+)\s+([e\d.-]+)\s+([e\d.-]+)\s+([e\d.-]+)\s+([e\d.-]+)';
        #$header = "seed\trep\tanalysis\tvariable\tmarker\tchr\tpos\tzscore\tdelta\tlod\tpvalue\texdelta\texlod\texpvalue";
        #$out_format = "%u\t%u\t%s\t%s\t%s\t%s\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n";
        $line_regex = '\s*([0-9XY]{1,3})\s+([e\d.-]+)\s+(\S+)\s+(\S+)\s+\[(\S+)\]\s+([e\d.-]+)\s+([e\d.-]+)\s+([e\d.-]+)\s+([e\d.-]+)';
        $header = "seed\trep\tanalysis\tvariable\tmarker\tchr\tpos\tzscore\tdelta\tlod\tpvalue";
        $out_format = "%u\t%u\t%s\t%s\t%s\t%s\t%g\t%g\t%g\t%g\t%g\n";
        $out_file = File::Spec->catfile($catdir,"npl-$dataset.txt");
        $glob = "$dir/*/*-$dataset-*/*-$dataset-*-nonparametric.tbl";
        @outfiles = glob($glob);
        if (@outfiles){
            open OUT,">",$out_file or die "can't open $out_file ($!)\n"  unless !$create_files;
            print OUT "$header\n";
            foreach my $file(@outfiles){
                $file =~ /$file_regex/gi or die "malformed $file\n";
                print "$file\n";
                my $filechr = $1;
                my $seed= $2;
                my $rep = $3;
                $npl_chr_seed{"$filechr-$seed"}++;
                open IN,"<",$file or die "can't open $file ($!)\n";
                while(<IN>){
                    chomp;
                    (my $line = $_) =~ s/\tinf/\t0/g;
                    $line =~ s/\tnan/\t0/g;               
                    $line =~ s/\tna/\t0/g;               
                    if($line =~ /$line_regex/){
                        #my ($chr,$pos,$marker,$variable,$analysis,$zscore,$delta,$lod,$pvalue,$exdelta,$exlod,$expvalue) = ($chr,$2,$3,$4,$codes{$5},$6,$7,$8,$9,$10,$11,$12);
                        #my @values = ($seed,$rep,$analysis,$variable,$marker,$chr,$pos,$zscore,$delta,$lod,$pvalue,$exdelta,$exlod,$expvalue);
                        my ($chr,$pos,$marker,$variable,$analysis,$zscore,$delta,$lod,$pvalue) = (chrnumber2chr($1),$2,$3,$4,$codes{$5},$6,$7,$8,$9);
                        if ($marker =~ /^[\d.]+$/){
                            $marker = "chr$chr"."cM$marker";
                        }
                        my @values = ($seed,$rep,$analysis,$variable,$marker,$chr,$pos,$zscore,$delta,$lod,$pvalue);
                        switch ($analysis){
                            case "a"{
                                $npl_all_record_count++;
                            }
                            case "p"{
                                $npl_pairs_record_count++;
                            }
                            case "q"{
                                $qtl_record_count++;
                                $qtls{$variable}++;
                            }
                        }
                        printf OUT $out_format, @values unless !$create_files;
                    }
                    elsif(!(/^CHR\s/ || /^na\s/)){
                         print "$_\n" unless /(^CHR)|(^na)/;
                    }
                }
            }
        }
        close OUT;
        close IN;
        print map{"$dataset\tnpl\tchr$_ $npl_chr_seed{$_}\tfiles\n"} sort keys %npl_chr_seed;
        print "$dataset\t$npl_all_record_count\tnpl-all\trecords\n";
        print "$dataset\t$npl_pairs_record_count\tnpl-pairs\trecords\n";
        print "$dataset\t$qtl_record_count\tqtl\trecords\n";
        print map{"$dataset\t$qtls{$_}\tqtl\t$_\trecords\n"} sort keys %qtls;
    }
    
    if($do_vc){
        my $vc_record_count = 0;
        my %vcs = ();
        my %vc_chr_seed = ();
        $file_regex = '\/chr([0-9XY]{1,2})-(\d{8})-(\d{5})-'.$dataset.'\/merlin-vc-chr[0-9XY]{1,3}\.tbl';
        $line_regex = '\s*([0-9XY]{1,3})\s+([e\d.-]+)\s+(\S+)\s+(\S+)\s+([e\d.-]+)\s+([e\d.-]+)\s+([e\d.-]+)';
        $header = "seed\trep\tvariable\tmarker\tchr\tpos\th2\tlod\tpvalue";
        $out_format = "%u\t%u\t%s\t%s\t%s\t%g\t%g\t%g\t%g\n";
        $out_file = File::Spec->catfile($catdir,"vc-$dataset.txt");
        $glob = "$dir/chr*/chr*-*-*-$dataset/chr*-*-*-$dataset-vc-chr*.tbl";
        @outfiles = glob($glob);
        if (@outfiles){
            open OUT,">",$out_file or die "can't open $out_file ($!)\n"  unless !$create_files;
            print OUT "$header\n";
            foreach my $file(@outfiles){
                $file =~ /$file_regex/g or die "malformed $file\n";
                my $chr = "$1";
                my $seed= $2;
                my $rep = $3;
                $vc_chr_seed{"$chr-$seed"}++;
                open IN,"<",$file or die "can't open $file ($!)\n";
                while(<IN>){
                    chomp;
                    (my $line = $_) =~ s/\tinf/\t0/g;
                    $line =~ s/\tnan/\t0/g;               
                    $line =~ s/\tna/\t0/g;               
                    if($line =~ /$line_regex/){
                        my ($chr,$pos,$marker,$variable,$h2,$lod,$pvalue) = ($chr,$2,$3,$4,$5,$6,$7);
                        if ($marker =~ /^[\d.]+$/){
                            $marker = "chr$chr"."cM$marker";
                        }
                        my @values = ($seed,$rep,$variable,$marker,$chr,$pos,$h2,$lod,$pvalue);
                        $vc_record_count++;
                        $vcs{$variable}++;
                        printf OUT $out_format, @values unless !$create_files;
                    }
                    elsif(!/^CHR/ && !/^na/){
                         print "$_\n" unless /(^CHR)|(^na)/;
                    }
                }
            }
        }
        close OUT;
        close IN;
        print map{"$dataset\tvc\tchr$_ $vc_chr_seed{$_}\tfiles\n"} sort keys %vc_chr_seed;
        print "$dataset\t$vc_record_count\tvc\trecords\n";
        print map{"$dataset\t$vcs{$_}\tvc-trait\t$_ records\n"} sort keys %vcs;
    }
    
    if($do_regress){
        my $regress_record_count = 0;
        my %regresses = ();
        my %regress_chr_seed = ();
        $file_regex = '\/(?:chr)?([0-9XY]{1,2}|(?:autosome)|x)-(\d{8})-(\d{5})-regress-(\S+)\.tbl';
        $line_regex = '\s*([0-9XY]{1,3})\s+(\S+)\s+Trait:\s+(\S+)\s+([e\d.-]+)\s+([e\d.-]+)\s+([e\d.-]+)\s+([e\d.-]+)\s+([e\d.-]+)';
        $header = "seed\trep\tvariable\tmarker\tchr\th2\tsd\tinfo\tlod\tpvalue";
        $out_format = "%u\t%u\t%s\t%s\t%s\t%g\t%g\t%g\t%g\t%g\n";
        $out_file = File::Spec->catfile($catdir,"regress-$dataset.txt");
        $glob = "$dir/*/*-$dataset-*/chr*-$dataset-*-regress-model-inverseNormal.tbl";
        @outfiles = glob($glob);
        if (@outfiles){
            open OUT,">",$out_file or die "can't open $out_file ($!)\n"  unless !$create_files;
            print OUT "$header\n";
            foreach my $file(@outfiles){
                $file =~ /$file_regex/g or die "malformed $file\n";
                print "$file\n";
                my $filechr = $1 eq "999" ? "XY" : $1;
                my $seed= $2;
                my $rep = $3;
                $regress_chr_seed{"$filechr-$seed"}++;
                open IN,"<",$file or die "can't open $file ($!)\n";
                while(<IN>){
                    chomp;
                    (my $line = $_) =~ s/\tblinf/\t0/g;
                    $line =~ s/\tnan/\t0/g;               
                    $line =~ s/\tna/\t0/g;               
                    if($line =~ /$line_regex/){
                        my ($chr,$marker,$variable,$h2,$sd,$info,$lod,$pvalue) = (chrnumber2chr($1),$2,$3,$4,$5,$6,$7,$8);
                        if ($marker =~ /^[\d.]+$/){
                            $marker = "chr$chr"."cM$marker";
                        }
                        my @values = ($seed,$rep,$variable,$marker,$chr,$h2,$sd,$info,$lod,$pvalue);
                        $regress_record_count++;
                        $regresses{$variable}++;
                        printf OUT $out_format, @values unless !$create_files;
                    }
                    elsif(!/^CHR/ && !/^na/){
                         print "$_\n" unless /(^CHR)|(^na)/;
                    }
                }
            }
        }
        close OUT;
        close IN;
        print map{"$dataset\tregress\tchr$_ $regress_chr_seed{$_}\tfiles\n"} sort keys %regress_chr_seed;
        print "$dataset\t$regress_record_count\tregress\trecords\n";
        print map{"$dataset\t$regresses{$_}\tregress-trait\t$_ records\n"} sort keys %regresses;
    }
}
print "\ndone\n";
exit;
sub chrnumber2chr{
    my $chr=shift;
    return $chr =~ /^[1-9]$/ ? "0$chr" : $chr eq "999" ? "XY" : $chr;
}