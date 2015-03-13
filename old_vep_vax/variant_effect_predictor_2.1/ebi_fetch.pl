#FILE: ebi_fetch.pl
#AUTH: Paul Stothard (paul.stothard@gmail.com)
#DATE: October 7, 2010
#VERS: 1.0

use warnings;
use strict;
use Getopt::Long;
use LWP::Simple;
use URI::Escape;

use LWP::UserAgent;
use HTTP::Request::Common;

my %options = (
    ids          => undef,
    database     => undef,
    output       => undef,
    format       => 'xml',
    style        => 'raw',
    verbose      => undef,
    url          => 'http://www.ebi.ac.uk/cgi-bin/dbfetch?',
    max_attempts => 10,
    retries      => 0
);

GetOptions(
    'i=s@{,}'   => \$options{ids},
    'd=s'       => \$options{database},
    'o=s'       => \$options{output},
    'f=s'       => \$options{format},
    's=s'       => \$options{style},
    'verbose|v' => \$options{verbose}
);

if (   ( !defined( $options{ids} ) )
    || ( !defined( $options{database} ) )
    || ( !defined( $options{format} ) )
    || ( !defined( $options{output} ) )
    || ( !defined( $options{style} ) ) )
{
    die( print_usage() . "\n" );
}

#build query
#see http://www.ebi.ac.uk/cgi-bin/dbfetch
#http://www.ebi.ac.uk/cgi-bin/dbfetch?db=EMBL&id=J00231,HSFOS,ROD894,LOP242600&format=fasta&style=text

my $request
    = $options{url} . 'db='
    . uri_escape( $options{database} ) . '&id='
    . uri_escape( join( ',', @{ $options{ids} } ) )
    . '&format='
    . uri_escape( $options{format} )
    . '&style='
    . uri_escape( $options{style} );

message( $options{verbose},
    "Retrieving records using the request '$request'.\n" );

my $result = get($request);

while ( !defined($result) ) {
    message( $options{verbose},
        "Request results could not be parsed. Trying again.\n" );
    sleep(10);
    if ( $options{retries} > $options{max_attempts} ) {
#MY
        #die("Too many failures--giving up search.");
        open( my $OUTFILE, ">" . $options{output} )
            or die("Cannot open file '$options{output}' : $!");
        close($OUTFILE) or die("Cannot close file : $!");
        return;
#MY
    }

    $result = get($request);
    $options{retries}++;
}

open( my $OUTFILE, ">" . $options{output} )
    or die("Cannot open file '$options{output}' : $!");
print( $OUTFILE $result );
close($OUTFILE) or die("Cannot close file : $!");

message( $options{verbose}, "Results written to '$options{output}'.\n" );

sub message {
    my $verbose = shift;
    my $message = shift;
    if ($verbose) {
        print $message;
    }
}

sub print_usage {
    print <<BLOCK;
USAGE: perl ebi_fetch.pl [-arguments]
 -i [STRINGS] : IDs of records to retrieve (Required).
 -d [STRING]  : database containing the records (Required).
 -o [FILE]    : the output file to create (Required).
 -f [STRING]  : the output format (Optional; default is 'xml').
 -s [STRING]  : the output style (Optional; default is 'raw').
 -v           : provide progress messages (Optional).

perl ebi_fetch.pl -i J00231 HSFOS ROD894 LOP242600 -d EMBL \
-f fasta -o sequences.fasta
BLOCK
}
