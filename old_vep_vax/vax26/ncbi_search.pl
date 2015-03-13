#!/usr/bin/perl
#FILE: ncbi_search.pl
#AUTH: Paul Stothard (paul.stothard@gmail.com)
#DATE: August 20, 2008
#VERS: 1.0

use warnings;
use strict;
use Getopt::Long;
use LWP::Simple;
use URI::Escape;

use LWP::UserAgent;
use HTTP::Request::Common;

my %param = (
    query      => undef,
    outputFile => undef,
    database   => undef,
    returnType => undef,
    maxRecords => undef,
    format     => undef,
    verbose    => undef,
    url        => 'http://www.ncbi.nlm.nih.gov/entrez/eutils',
    retries    => 0,
    maxRetries => 10
);

Getopt::Long::Configure('bundling');
GetOptions(
    'q|query=s'       => \$param{query},
    'o|output_file=s' => \$param{outputFile},
    'd|database=s'    => \$param{database},
    'r|return_type=s' => \$param{returnType},
    'm|max_records=i' => \$param{maxRecords},
    'verbose|v'       => \$param{verbose}
);

if (   !( defined( $param{query} ) )
    or !( defined( $param{outputFile} ) )
    or !( defined( $param{database} ) )
    or !( defined( $param{returnType} ) ) )
{
    die( print_usage() . "\n" );
}

$param{returnType} = lc( $param{returnType} );

$param{query} = uri_escape( $param{query} );

_doSearch(%param);

sub _doSearch {
    my %param = @_;

    my $esearch = "$param{url}/esearch.fcgi?db=$param{database}"
        . "&retmax=1&usehistory=y&term=$param{query}";
    my $esearch_result = get($esearch);

    while (
        ( !defined($esearch_result) )
        || (!(  $esearch_result
                =~ m/<Count>(\d+)<\/Count>.*<QueryKey>(\d+)<\/QueryKey>.*<WebEnv>(\S+)<\/WebEnv>/s
            )
        )
        )
    {
        message( $param{verbose},
            "ESearch results could not be parsed. Resubmitting query.\n" );
        sleep(10);
        if ( $param{retries} >= $param{maxRetries} ) {
#MY
            #die("Too many failures--giving up search.");
            return;
#MY
        }

        $esearch_result = get($esearch);
        $param{retries}++;
    }

    $param{retries} = 0;

    $esearch_result
        =~ m/<Count>(\d+)<\/Count>.*<QueryKey>(\d+)<\/QueryKey>.*<WebEnv>(\S+)<\/WebEnv>/s;

    my $count     = $1;
    my $query_key = $2;
    my $web_env   = $3;

    if ( defined( $param{maxRecords} ) ) {
        if ( $count > $param{maxRecords} ) {
            message( $param{verbose},
                "Retrieving $param{maxRecords} records out of $count available records.\n"
            );
            $count = $param{maxRecords};
        }
        else {
            message( $param{verbose},
                "Retrieving $count records out of $count available records.\n"
            );
        }
    }
    else {
        message( $param{verbose},
            "Retrieving $count records out of $count available records.\n" );
    }

    my $retmax = 500;
    if ( $retmax > $count ) {
        $retmax = $count;
    }

    open( my $OUTFILE, ">" . $param{outputFile} )
        or die("Error: Cannot open $param{outputFile} : $!");

    for (
        my $retstart = 0;
        $retstart < $count;
        $retstart = $retstart + $retmax
        )
    {
        message( $param{verbose},
                  "Downloading records "
                . ( $retstart + 1 ) . " to "
                . ( $retstart + $retmax )
                . "\n" );
        my $efetch
            = "$param{url}/efetch.fcgi?rettype=$param{returnType}&retmode=text&retstart=$retstart&retmax=$retmax&db=$param{database}&query_key=$query_key&WebEnv=$web_env";
        my $efetch_result = get($efetch);

        while ( !defined($efetch_result) ) {
            message( $param{verbose},
                "EFetch results could not be parsed. Resubmitting query.\n" );
            sleep(10);
            if ( $param{retries} >= $param{maxRetries} ) {
                #die("Too many failures--giving up search.");
#MY
                close($OUTFILE) or die("Error: Cannot close $param{outputFile} file: $!");
                return;
#MY
            }

            $efetch_result = get($efetch);
            $param{retries}++;
        }

        print( $OUTFILE $efetch_result );

        unless ( $param{maxRecords} == 1 ) {
            sleep(3);
        }
    }

    close($OUTFILE) or die("Error: Cannot close $param{outputFile} file: $!");
}

sub message {
    my $verbose = shift;
    my $message = shift;
    if ($verbose) {
        print $message;
    }
}

sub print_usage {
    print <<BLOCK;
USAGE: perl ncbi_search.pl [-arguments]
 -q [STRING]     : raw query text (Required).
 -o [FILE]       : output file to create (Required).
 -d [STRING]     : name of the NCBI database to search, such as 'nucleotide', 
                   'protein', or 'gene' (Required).
 -r [STRING]     : the type of information requested. For sequences 'fasta' is
                   often used. The accepted formats vary depending on the 
                   database being queried (Required).
 -m [INTEGER]    : the maximum number of records to return (Optional; default
                   is to return all matches satisfying the query).
 -v              : provide progress messages (Optional).

perl ncbi_search.pl -q 'dysphagia AND homo sapiens[ORGN]' \
-o results.txt -d pubmed -r uilist -m 100
BLOCK
}
