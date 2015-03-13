#!/usr/bin/perl -w
use strict;
use warnings;

#BioPerl
use Bio::Seq;
use Bio::SeqIO;

my @regions;

#redirect STDOUT, which the debugger doesn't do
if ($ARGV[$#ARGV - 1] eq ">") {
	my $output_file = pop(@ARGV);
	pop(@ARGV);
	open STDOUT, ">", $output_file or die "can't open $output_file: $!";
}
while (<>) {
	chomp;
    if (/^\s*(?:chr)?[\s:]*([XYMT]|[0-9]{1,2})[\s:]+([,\d]+)\D+([,\d]+)(?:\s+([01+-]{1}))?/i){
        if ($1 && $2 && $3) {
            my $chr = $1;
            my $start = $2;
            my $end = $3;
            my $strand = defined($4)
             && $4 == "+" ? 1 : defined($4)
             && $4 == "-" ? -1 : $4;
            $start =~ s/,//g;
            $end =~ s/,//g;
            $start = $start <= $end ? $start : $end;
            $end = $end >= $start ? $end : $start;
            push @regions,sprintf("chr%s:%d-%d",$chr,$start,$end);
            #push @regions,{chr=>$chr,start=>$start,end=>$end};
        }
    }
}
@regions = sortregions(\@regions,1);
print STDERR join("\n",@regions);
print join("\n",@regions);

sub sortregions{
	#sort genomic regions by chromosome,start,end
    my ($regions,$zeros) = @_;
    #to sort chromosomes as integers
	my @zeroregions;
    foreach my $region(@$regions){
        $region =~ s/chrMT?/chr97/i;
        $region =~ s/chrX/chr98/i;
        $region =~ s/chrY/chr99/i;
		if(defined($zeros)){
			my $region =~ /chr(\d+):(\d+)-(\d+)/;
			push @zeroregions, sprintf("chr%d:%09d-%09d",$1,$2,$3);
		}
    }
	if(defined($zeros)){
		@$regions = @zeroregions;
	}
    use Sort::Key::Multi qw(iiikeysort);
    my @sorted = iiikeysort {/chr(\d+):(\d+)-(\d+)/} @$regions;
    foreach my $region(@sorted){
        $region =~ s/chr97/chrM/;
        $region =~ s/chr98/chrX/;
        $region =~ s/chr99/chrY/;
    }
    return @sorted;
}

sub sortregionslexical{
	#sort by chromosome, then by "start-end"
    my $regions = shift;
    #to sort chromosomes as integers
    foreach my $region(@$regions){
        $region =~ s/chrMT?/chr97/i;
        $region =~ s/chrX/chr98/i;
        $region =~ s/chrY/chr99/i;
    }
    #use Sort::Key qw(keysort nkeysort ikeysort);
    use Sort::Key::Multi qw(iskeysort);
    my @sorted = iskeysort {/chr(\d+):([\d-]+)/} @$regions;
    foreach my $region(@sorted){
        $region =~ s/chr97/chrM/;
        $region =~ s/chr98/chrX/;
        $region =~ s/chr99/chrY/;
    }
    return @sorted;
}
