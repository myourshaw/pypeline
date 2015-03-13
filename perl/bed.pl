#!/usr/bin/perl -w
use strict;
use warnings;
use XML::Simple;
  my $self = XMLin('C:\bio\GMOD\Chado\schema\chado\load\etc\load.conf',
    ForceArray  => ['token','path','file'],
    KeyAttr     => [qw(tt2 input token name file)],
    ContentKey  => '-value'
  );

my %positions;
use constant CHR_SHIFT => 1000000000;
while(<>){
    if(/chr([0-9XYMT]{1,2})\s(\d+)\s(\d+)\s([ACGT])->([ACGT])\((\d+):(\d+):([0-9.]+)%\[F:(\d+):([0-9.]+)%\|R:(\d+):([0-9.]+)%\]\)\s(\d+)\s([+-])/i){
        my $chr = $1;
        my $start = $2;
        my $end = $3;
        my $ref = uc($4);
        my $var = uc($5);
        my $reads = $6;
        my $var_reads = $7;
        my $var_percent = $8;
        my $var_fwd_reads = $9;
        my $var_fwd_percent = $10;
        my $var_rev_reads = $11;
        my $var_rev_percent = $12;
        my $score = $13;
        my $strand = $14;
        my $position = CHR_SHIFT*chr2numchr($chr)+$start;
        $positions{$position} = {
            ref => $ref,
            reads => $reads,
            score => $score,
            strand => $strand
        };
            if($var eq "A"){
                $positions{$position}->{A} = $var_reads;
                $positions{$position}->{Af} = $var_fwd_reads;
                $positions{$position}->{Ar} = $var_rev_reads;
            }
            if($var eq "C"){
                $positions{$position}->{C} = $var_reads;
                $positions{$position}->{Cf} = $var_fwd_reads;
                $positions{$position}->{Cr} = $var_rev_reads;
            }
            if($var eq "G"){
                $positions{$position}->{G} = $var_reads;
                $positions{$position}->{Gf} = $var_fwd_reads;
                $positions{$position}->{Gr} = $var_rev_reads;
            }
            if($var eq "T"){
                $positions{$position}->{T} = $var_reads;
                $positions{$position}->{Tf} = $var_fwd_reads;
                $positions{$position}->{Tr} = $var_rev_reads;
            }
    }
}
foreach my $position (keys %positions){
    unless (defined $positions{$position}->{A}){
        $positions{$position}->{A} = 0;
        $positions{$position}->{Af} = 0;
        $positions{$position}->{Ar} = 0;
    }
    unless (defined $positions{$position}->{C}){
        $positions{$position}->{C} = 0;
        $positions{$position}->{Cf} = 0;
        $positions{$position}->{Cr} = 0;
    }
    unless (defined $positions{$position}->{G}){
        $positions{$position}->{G} = 0;
        $positions{$position}->{Gf} = 0;
        $positions{$position}->{Gr} = 0;
    }
    unless (defined $positions{$position}->{T}){
        $positions{$position}->{T} = 0;
        $positions{$position}->{Tf} = 0;
        $positions{$position}->{Tr} = 0;
    }
    $positions{$position}->{ref_reads} = $positions{$position}->{reads}-$positions{$position}->{A}-$positions{$position}->{C}-$positions{$position}->{G}-$positions{$position}->{T};
}
print map { numchr2chr(int($_/CHR_SHIFT))."\t".$_ % CHR_SHIFT."\t$positions{$_}->{ref}\t".$positions{$_}->{ref_reads}."\t".$positions{$_}->{A}."\t".$positions{$_}->{C}."\t".$positions{$_}->{G}."\t".$positions{$_}->{T}."\n" } sort {$a <=> $b} keys %positions;

sub chr2numchr{
    my $chr = shift;
    return $chr eq "X" ? 23 : $chr eq "Y" ? 24 : ($chr eq "M" || $chr eq "MT") ? 25 : $chr;
}

sub numchr2chr{
    my $chr = shift;
    return $chr == 23 ? "X" : $chr == 24 ? "Y" : $chr == 25 ? "M" : $chr;
}
