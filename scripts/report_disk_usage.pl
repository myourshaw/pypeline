#!/usr/bin/perl
#$ -S /usr/bin/perl
use warnings;
use strict;
use Data::Dumper;

my $time = time;

my %org = (
  young => 'Young Koo - Former Visiting Scholar',
  gec_10004 => 'GEC Microarray data',
  ayang => 'Anqi Yang - Former Rotation Student',
  jdong => 'Jun Dong - Former Biostatistician',
  microarray => 'GEC Microarray data',
  corinas => 'Corina Shtir - Former Postdoc',
  yohanlee => 'Yohan Lee - Former Gradstudent',
  gec_491 => 'GEC Microarray data',
  awang => 'Anyou Wang - Guoping Lab',
  fah => 'Jarupon Fah Sathirapongsasuti - Quackenbush Lab Boston - Dept. of Biostatistics',
  gec_65534 => 'GEC Microarray data',
  traci => 'Traci Toy',
  mikyung => 'Mikyung Lee - Postdoc - Roger Lo Lab',
  reshmi => 'Reshmi Chowdhury - Programer Analyst - Albert Lai Lab',
  barrym => 'Barry Merriman',
  boconnor => 'Brian O\'Connor',
  yuan => 'Yuan Tian - Former Rotation Student',
  rwang => 'Richard Wang',
  celsius => 'Celsius data',
  nudar => 'Nitin Udar - Volunteer ?',
  'HiSeq2k.Pathology' => 'Data from Pathogy HiSeq2000 instrument',
  seqware => 'Data from Seqware (Brian\'s project)',
  aeskin => 'Ascia Eskin',
  aliz => 'Aliz Raksi',
  adamf => 'Adam Freedman - Postdoc - Novembre Lab',
  nhomer => 'Nils Homer - Former Gradstudent',
  sam => 'Sam Strom - Former Gradstudent',
  natalia => 'Natalia Shatokhina - Gradstudent - HEI - NTR Lab',
  dmitriy => 'Dmitriy Skvortsov - Former Postdoc',
  vchang => 'Vivian Chang - Fellowship Student',
  tpaige => 'Paige Taylor',
  gec_507 => 'GEC Microarray data',
  root => 'Sequence data that needs to be cleaned, Bret Harry',
  bret => 'Bret Harry',
  mclark => 'Michael James Clark - Former Gradstudent',
  biouser => 'Sequenom data',
  hlee => 'Hane Lee',
  kmsquire => 'Kevin M. Squire',
  solexa => 'Sequence data that needs to be organized, Bret Harry',
  myourshaw => 'Michael Yourshaw'
);

my %exempt = (
	      solexa => '',
	      biouser => '',
	      root => '',
	      dmitriy => '',
	      nhomer => '',
	      gec_10004 => '',
	      gec_65534 => '',
	      microarray => '',
	      gec_491 => '',
	      gec_507 => ''
	     );

###
### Get information from df.
###

my %df = (
	  archive => {
		      '/data/locus' => {},
#		      '/data/vault' => {},
		      '/data/storage-1-00' => {},
		      '/data/storage-1-01' => {},
		      '/data/storage-1-02' => {},
		      '/data/storage-1-03' => {},
		      '/data/storage-1-04' => {}
		     },
	  scratch => {
		      '/scratch0' => {},
		      '/scratch1' => {}
		     }
	 );

sub get_free_bytes {
  my $target = shift @_;
  unless(-d $target) {
    return 0;
  }

  my @output = `df -B1 $target`;
  my @ans = split(/\s+/, $output[1]);

  if(not defined $ans[3]) {
    my @ans = split(/\s+/, $output[2]); ## Line begins with space :(
    if(not defined $ans[3]) {
      return 0;
    } else {
      return $ans[3];
    }
  } else {
    return $ans[3];
  }
}

sub get_max_bytes {
  my $target = shift @_;
  unless(-d $target) {
    return 0;
  }

  my @output = `df -B1 $target`;
  my @ans = split(/\s+/, $output[1]);

  if(not defined $ans[3]) {
    my @ans = split(/\s+/, $output[2]); ## Line begins with space :(
    if(not defined $ans[1]) { ### diff &get_free_bytes();
      return 0;
    } else {
      return $ans[1]; ### diff &get_free_bytes();
    }
  } else {
    return $ans[1]; ### diff &get_free_bytes();
  }
}


###
### Figure out which systems are >95%, for global quota.
###

for my $t (keys %df) {
  for my $k (keys %{$df{$t}}) {
    my $free = &get_free_bytes($k);
    my $max = &get_max_bytes($k);
    my $volume = $k;
    $volume =~ s/^\/(data\/)?//;

    $df{$t}->{$k}->{free} = $free;
    $df{$t}->{$k}->{max} = $max;
    $df{$t}->{$k}->{volume} = $volume;

  }
}

my %total;
my %total_by_type;
my %total_by_volume;

my @files;
push @files, glob "/home/bret/bin/mybs/data/disk_usage_xfs_*txt";
push @files, glob "/home/bret/bin/mybs/data/disk_usage_lustre_*txt";

for my $usage_file ( @files ) {
  next if $usage_file =~ /vault/; ### ignore this volume.
  $usage_file =~ /disk_usage_(?:xfs|lustre)_([^\.]+).txt/;
  my $volume = $1;
  $volume =~ s/_/-/g;
  my $type = "archive";
  $type = "scratch" if $usage_file =~ /lustre/;

  my $last_line;
  open(F,$usage_file) or die $!;
  while(<F>) {
    $last_line = $_;
  }
  close(F) or die $!;
  chomp $last_line;
  $last_line =~ s/(\d+)://;
  if($time - $1 > 7200) { # Usage data is over 2 hours old.
    print "WARNING OLD DATA: $usage_file\n";
  }
  my @data = split(/;/, $last_line);

  while(my $datum = shift @data) {
    my ($user, $value) = split /,/, $datum;
    if($user =~ /^\d+$/) {
      my $uid = $user;
      $user = getpwuid($user);
      $user = $uid if not defined $user;
    }
    $total_by_volume{$volume}->{$user} += $value;
    $total_by_type{$type}->{$user} += $value;
    $total{$user} += $value;
  }
}

#&print_top_users();

for my $t (keys %df) {
  for my $k (keys %{$df{$t}}) {
    my $free = $df{$t}->{$k}->{free};
    my $max  = $df{$t}->{$k}->{max};
    my $volume = $df{$t}->{$k}->{volume};
    my $percent = int(100 - ($free / $max * 100));
    my $fmt = "%11.11s - %-12.12s %3.3s%%\n";
    my $quota;
    my %threshhold =
      (
       archive => {
		   hard => '95',
		   soft => '85',
		   reset => '80',
		  },
       scratch => {
		   hard => '90',
		   soft => '75',
		   reset => '70',
		  }
      );

    if($percent >= $threshhold{$t}->{hard}) {
      printf($fmt,"Hard Quota",$volume,$percent);
      $quota = "hard";
    } elsif($percent >= $threshhold{$t}->{soft}) {
      printf($fmt,"Soft Quota",$volume,$percent);
      $quota = "soft";
    } elsif($percent <= $threshhold{$t}->{reset}) {
      $quota = "reset";
    } else {
      $quota = "keep";
    }

    if($quota eq "hard" or $quota eq "soft") {
      if(&has_greedy_users($volume)) {
	#&print_greedy_users($volume);
	&print_top_consumers($volume);
      } else {
	&print_top_consumers($volume);
      }
    } elsif ($quota eq "reset") {
      # printf($fmt,"Reset Quota",$volume,$percent);
      # print "\n";
    } elsif ($quota eq "keep") {
      # printf($fmt,"Keep Quota",$volume,$percent);
      # print "\n";
    }
  }
}

#exit; ### NO ASCII REPORT BY DEFAULT
{
  my %usage;
  my $total;
  print "BIG ASCII REPORT\n";
  for my $type (keys %df) {
    for my $volume (keys %{$df{$type}}) {
#      print "  $volume\n";
      my $summary ="$volume/big.ascii.files.summary.txt";
      ### Need to do a stat here  and skip over this file if it's old.
      if(-e $summary) {
	open(F,$summary) or die $!;
	while(<F>) {
	  die unless $_ =~ /^(\d+) (\d+)/;
#	  print "$1 $2\n";
	  $usage{$1} += $2;
	  $total += $2;
	}
	close(F) or die $!;
      }
    }
  }

  for my $uid ( sort {
    return $usage{$b} <=> $usage{$a};
  } keys %usage) {
    my $user = getpwuid($uid);
    $user = $uid if not defined $user;
    $user = "solid" if $user eq "501";
    printf("%-16s %12s\n",$user,&get_filesize_str($usage{$uid}));
  }	   

  printf("%-16s %12s\n","total",&get_filesize_str($total));
#  print "$total\n";
#  print Dumper \%usage;
}

sub get_max_by_volume {
  my $volume = shift @_;
  for my $t (keys %df) {
    for my $k (keys %{$df{$t}}) {
      if($df{$t}->{$k}->{volume} eq $volume) {
	return $df{$t}->{$k}->{max};
      }
    }
  }
  die "Did not find max for volume: $volume\n";
}

sub print_top_consumers {
  my $volume = shift @_;
  my $first = 1;
  my $report;

  for my $user ( sort {$total_by_volume{$volume}->{$b} <=> $total_by_volume{$volume}->{$a}}
		 keys %{$total_by_volume{$volume}}) {
    my $used = $total_by_volume{$volume}->{$user};
    my $percent = int( $used / &get_max_by_volume($volume) * 100);


    if($percent < 5) {
      next;
    }

    if(exists $exempt{$user}) {
      next;
    }
    if($first) {

      $first = 0;
      print "\n";
      print "                       --- TOP USERS ---\n";
    }
    my $fmt = "                %10.10s %3.3s%%  %s\n";
    printf($fmt,&get_filesize_str($used),$percent,$user);
  }
  print "\n" unless $first;
}

### Wish I had something more sophisticated here...
sub has_greedy_users {
  my $volume = shift @_;

  my @sorted_user;
  for my $user ( sort {$total_by_volume{$volume}->{$b} <=> $total_by_volume{$volume}->{$a}}
		 keys %{$total_by_volume{$volume}}) {
    if(exists $exempt{$user}) {
      next;
    }
    push @sorted_user, $user;
  }

  for(my $i=0; $i < scalar @sorted_user; $i++) {
    my $user_a = $total_by_volume{$volume}->{$sorted_user[$i]};
    my $user_b = 0;
    my $user_c = 0;

    if(exists $sorted_user[$i+1]) {
      $user_b = $total_by_volume{$volume}->{$sorted_user[$i+1]};
    }

    if(exists $sorted_user[$i+2]) {
      $user_c = $total_by_volume{$volume}->{$sorted_user[$i+2]};
    }

    if($user_a > ($user_b + $user_c)) {
      return 1;
    } else {
      return 0;
    }
  }


}

sub print_greedy_users {
  my $volume = shift @_;
  my $first = 1;
  my @sorted_user;
  for my $user ( sort {$total_by_volume{$volume}->{$b} <=> $total_by_volume{$volume}->{$a}}
		 keys %{$total_by_volume{$volume}}) {
    if(exists $exempt{$user}) {
      next;
    }

    push @sorted_user, $user;

  }

  for(my $i=0; $i < scalar @sorted_user; $i++) {
    my $user_a = $total_by_volume{$volume}->{$sorted_user[$i]};
    my $user_b = 0;
    my $user_c = 0;

    if(exists $sorted_user[$i+1]) {
      $user_b = $total_by_volume{$volume}->{$sorted_user[$i+1]};
    }

    if(exists $sorted_user[$i+2]) {
      $user_c = $total_by_volume{$volume}->{$sorted_user[$i+2]};
    }

    if($user_a > ($user_b + $user_c)) {
      my $user = $sorted_user[$i];
      my $used = $total_by_volume{$volume}->{$user};
      my $percent = int( $used / &get_max_by_volume($volume) * 100);

#      my $fmt = "                %10.10s %3.3s%%  %s\n";
      my $fmt = "                %10.10s %3.3s%%  %s\n";
      if($first) {
	print "\n";
	print "                    --- GREEDY USERS ---\n";
#	printf($fmt,"","","GREEDY USERS");
	$first = 0;
      }

      printf($fmt,&get_filesize_str($used),$percent,$user);
    } else {
      last; ### There are no greedy users, buy more storage.
    }
  }
  print "\n" unless $first;
}

sub get_filesize_str {
  my $size = shift; 
  return "" if not defined $size;

    if ($size > 1099511627776)  #   TiB: 1024 GiB
    {
        return sprintf("%.2f T", $size / 1099511627776);
    }
    elsif ($size > 1073741824)  #   GiB: 1024 MiB
    {
        return sprintf("%.2f G", $size / 1073741824);
    }
    elsif ($size > 1048576)     #   MiB: 1024 KiB
    {
        return sprintf("%.2f M", $size / 1048576);
    }
    elsif ($size > 1024)        #   KiB: 1024 B
    {
        return sprintf("%.2f K", $size / 1024);
    }
    else                        #   bytes
    {
        return sprintf("%.2f bytes", $size);
    }
  return "";
}

sub print_top_users {
  my %max = (
	     scratch => 0,
	     archive => 0
	    );
  my %used = (
	    scratch => 0,
	    archive => 0
	   );

  for my $t (keys %max) {
    for my $v (keys %{$df{$t}}) {
      $max{$t} += $df{$t}->{$v}->{max};
      $used{$t} += ($df{$t}->{$v}->{max} - $df{$t}->{$v}->{free});
    }
  }

  my $fmt = "%10.10s %4.4s %10.10s %4.4s   %s\n";

  my $cutoff = 0.75;
  my $first = 1;
  for my $user (sort {$total{$a} <=> $total{$b}} keys %total) {

    if($used{scratch} / $max{scratch} > $cutoff) {
    } elsif( $used{archive} / $max{archive} > $cutoff) {
    } else {
      last; ### The filesystem usage is not an issue.
    }

    if($first) {
      my $total_scratch = int($used{scratch} / $max{scratch} * 100);
      $total_scratch .= "%";
      my $total_archive = int($used{archive} / $max{archive} * 100);
      $total_archive .= "%";

      print "Global Usage - In excess of 75%\n\n";
      printf($fmt, "scratch",$total_scratch, "archive",$total_archive, "\n");

      $first = 0;
    }

    my $percent_scratch;
    if(exists $total_by_type{scratch}->{$user}) {
      $percent_scratch = int($total_by_type{scratch}->{$user} / $max{scratch} * 100 );
    }

    my $percent_archive;
    if(exists $total_by_type{archive}->{$user}) {
      $percent_archive = int($total_by_type{archive}->{$user} / $max{archive} * 100 );
    }

    if (defined $percent_scratch and $percent_scratch > 1) {
    } elsif( defined $percent_archive and $percent_archive > 1) {
    } else {
      next; ### This user is below 1% usage.
    }

    if(exists $exempt{$user}) {
      next; ### This user is exempt.
    }

    if(not defined $percent_scratch) {
      $percent_scratch = "";
    } else {
      $percent_scratch .= "%";
      $percent_scratch = "" if $percent_scratch eq "0%";
    }

    if(not defined $percent_archive) {
      $percent_archive = "";
    } else {
      $percent_archive .= "%";
      $percent_archive = "" if $percent_archive eq "0%";
    }

    ### Ignore usage less than 1%.
    printf($fmt,
	   &get_filesize_str($total_by_type{scratch}->{$user}),
	   $percent_scratch,
	   &get_filesize_str($total_by_type{archive}->{$user}),
	   $percent_archive,
	   $org{$user} ? $org{$user} : $user);
  }
  print "\n";
}


