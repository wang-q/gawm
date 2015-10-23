#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw(HelpMessage);
use Config::Tiny;
use FindBin;
use YAML qw(Dump Load DumpFile LoadFile);

use MongoDB;
$MongoDB::BSON::looks_like_number = 1;
use MongoDB::OID;

use AlignDB::Stopwatch;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->read("$FindBin::RealBin/config.ini");

# record ARGV and Config
my $stopwatch = AlignDB::Stopwatch->new(
    program_name => $0,
    program_argv => [@ARGV],
    program_conf => $Config,
);

=head1 NAME

gsw_csv.pl - Write gsws to a csv file

=head1 SYNOPSIS

    perl gsw_csv.pl [options]
      Options:
        --help      -?          brief help message
        --server        STR     MongoDB server IP/Domain name
        --port          INT     MongoDB server port
        --db        -d  STR     database name
        --output    -o  STR     output filename, default is [db.gsw.csv]

=cut

GetOptions(
    'help|?' => sub { HelpMessage(0) },
    'server=s' => \( my $server = $Config->{database}{server} ),
    'port=i'   => \( my $port   = $Config->{database}{port} ),
    'db|d=s'   => \( my $dbname = $Config->{database}{db} ),
    'output|o=i' => \my $outfile,
) or HelpMessage(1);

$outfile = "$dbname.gsw.csv" unless $outfile;

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
$stopwatch->start_message("Write gsws for $dbname...");

my $db = MongoDB::MongoClient->new(
    host          => $server,
    port          => $port,
    query_timeout => -1,
)->get_database($dbname);

#----------------------------------------------------------#
# Write
#----------------------------------------------------------#

## Original SQL
#select
#    gsw_distance distance,
#    gsw_wave_length wave_length,
#    gsw_amplitude amplitude,
#    gsw_trough_gc trough_gc,
#    gsw_gradient gradient,
#    window_average_gc window_gc,
#    gsw_cv window_cv,
#    gsw_intra_cv intra_cv,
#    window_indel window_indel,
#    if(window_indel > 0, 1, 0) has_indel
#from
#    gsw g
#        inner join
#    window w ON g.window_id = w.window_id
#where
#    gsw_amplitude >= 0.1;

open my $fh, '>', $outfile;
my @headers = qw{
    distance_to_trough distance_to_crest window_gc window_cv trough_gc crest_gc
    gradient value flag
};

print {$fh} join( ",", @headers ), "\n";

my $coll   = $db->get_collection('gsw');
my $cursor = $coll->find;

#my $cursor = $coll->find( { chr_name => 'chr1' } );
while ( my $object = $cursor->next ) {
    print {$fh} join( ",",
        $object->{distance},  $object->{distance_crest},
        $object->{gc},        $object->{gc_cv},
        $object->{trough_gc}, $object->{crest_gc},
        $object->{gradient},  $object->{bed_count},
        $object->{bed_count} > 0 ? 1 : 0, ),
        "\n";
}

$stopwatch->end_message;

exit;
