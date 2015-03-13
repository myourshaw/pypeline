=head1 LICENSE

 Copyright (c) 2011-2012 Michael Yourshaw.  All rights reserved.                                                                      

=head1 CONTACT                                                                                                       

 Michael Yourshaw <myourshaw@ucla.edu>
    
=cut

=head1 NAME

 vw

=head1 SYNOPSIS

 mv vw.pm ~/.vep/Plugins
 perl variant_effect_predictor.pl -i variations.vcf --plugin vw[,cortex.local,3306,vw,vw,mysql,vw]

=head1 DESCRIPTION

 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that
 provides a connection to the vw databaase:
 
=head1 PARAMETERS

    host (default: cortex.local),port (default: 3306),user (default: vw),password (default: vw),platform (default: mysql),database (default: vw)

=cut

package vw;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);
use DBI;
our $vw_conn;
our $vw_database;

sub new {
    my $class = shift;
    my $self = $class->SUPER::new(@_);
    
    my ($host,$port,$user,$password,$platform,$database) = @{$self->{params}};
    $host ||= 'cortex.local';
    $port ||= 3306;
    $user ||= 'vw';
    $password ||= 'vw';
    $platform ||= 'mysql';
    $database ||= 'vw';
    $vw_database = $database;
    $self->{vw_database} = $vw_database;
    
    my $dsn = "dbi:$platform:$database:$host:$port";
    my $conn = DBI->connect($dsn, $user, $password)
        or die "Unable to connect: $DBI::errstr\n";
    $vw_conn = $conn;
    $self->{conn} = $conn;
    
    return $self;
}

sub version {
    return '2.4';
}

sub feature_types {
    return [];
}

sub get_header_info {
    return {};
}

sub run {
    return {};
}


1;

