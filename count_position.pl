#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Config::Tiny;
use FindBin;
use YAML::Syck;

use MCE;
use MongoDB;
$MongoDB::BSON::looks_like_number = 1;
use MongoDB::OID;

use AlignDB::IntSpan;
use AlignDB::Window;
use AlignDB::Stopwatch;
use App::RL::Common;

use lib "$FindBin::RealBin/lib";
use MyUtil qw(check_coll);

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->read("$FindBin::RealBin/config.ini");

# record ARGV and Config
my $stopwatch = AlignDB::Stopwatch->new();

=head1 NAME

count_position.pl - Add position files and count intersections.

=head1 SYNOPSIS

    perl count_position.pl [options]
      Options:
        --help      -?          brief help message
        --server        STR     MongoDB server IP/Domain name
        --port          INT     MongoDB server port
        --db        -d  STR     database name
        --file      -f  @STR    bed files
        --run       -r  STR     all, insert or count, default is [all]
        --nochr                 remove 'chr' from chr_name
        --parallel      INT     run in parallel mode, [1]
        --batch         INT     alignments processed in one child process, [10]

    perl count_bed.pl -d S288c --run insert -f ~/Scripts/alignDB/ofg/spo11/spo11_hot.bed --batch 1 --parallel 1
    perl count_bed.pl -d S288c --run count --batch 1 --parallel 4

=cut

GetOptions(
    'help|?' => sub { Getopt::Long::HelpMessage(0) },
    'server=s'   => \( my $server       = $Config->{database}{server} ),
    'port=i'     => \( my $port         = $Config->{database}{port} ),
    'db|d=s'     => \( my $dbname       = $Config->{database}{db} ),
    'file|f=s'   => \my @files,
    'run|r=s'    => \( my $run          = "all" ),
    'parallel=i' => \( my $parallel     = 1 ),
    'batch=i'    => \( my $batch_number = 10 ),
) or Getopt::Long::HelpMessage(1);

#----------------------------------------------------------#
# Init objects
#----------------------------------------------------------#
$stopwatch->start_message("Count positions of $dbname...");

$stopwatch->block_message( "Total position files are " . scalar(@files) ) if @files;

#----------------------------------------------------------#
# workers
#----------------------------------------------------------#
my $worker_insert = sub {
    my ( $self, $chunk_ref, $chunk_id ) = @_;
    my $file = $chunk_ref->[0];

    my $wid = MCE->wid;

    my $inner_watch = AlignDB::Stopwatch->new;
    $inner_watch->block_message("Process task [$chunk_id] by worker #$wid");

    print "Reading file [$file]\n";

    # wait forever for responses
    #@type MongoDB::Database
    my $db = MongoDB::MongoClient->new(
        host          => $server,
        port          => $port,
        query_timeout => -1,
    )->get_database($dbname);

    my $coll_align = $db->get_collection('align');

    my @beds;
    open my $data_fh, '<', $file;
    while ( my $string = <$data_fh> ) {
        next unless defined $string;
        chomp $string;

        my $info = App::RL::Common::decode_header($string);
        next unless defined $info->{chr};

        my $align = $coll_align->find_one(
            {   'chr.name'  => $info->{chr},
                'chr.start' => { '$lte' => $info->{start} },
                'chr.end'   => { '$gte' => $info->{end} }
            }
        );

        if ( !$align ) {
            print "    Can't locate an align for $string\n";
            next;
        }
        push @beds,
            {
            align => { _id => $align->{_id}, },
            chr   => {
                name    => $info->{chr},
                start   => $info->{start},
                end     => $info->{end},
                runlist => AlignDB::IntSpan->new->add_pair( $info->{start}, $info->{end} )->runlist,
            },
            };
    }
    close $data_fh;

    print "Inserting file [$file]\n";

    #@type MongoDB::Collection
    my $coll_bed = $db->get_collection('bed');
    while ( scalar @beds ) {
        my @batching = splice @beds, 0, 10000;
        $coll_bed->insert_many( \@batching );
    }
    print "Insert done.\n";
};

my $worker_count = sub {
    my ( $self, $chunk_ref, $chunk_id ) = @_;
    my @jobs = @{$chunk_ref};

    my $wid = MCE->wid;

    my $inner_watch = AlignDB::Stopwatch->new;
    $inner_watch->block_message("Process task [$chunk_id] by worker #$wid");

    # wait forever for responses
    #@type MongoDB::Database
    my $db = MongoDB::MongoClient->new(
        host          => $server,
        port          => $port,
        query_timeout => -1,
    )->get_database($dbname);

    my $coll_align = $db->get_collection('align');
    my $coll_gsw   = $db->get_collection('gsw');
    my $coll_ofgsw = $db->get_collection('ofgsw');
    my $coll_bed   = $db->get_collection('bed');

    for my $job (@jobs) {
        my $align = $coll_align->find_one( { _id => $job->{_id} } );
        if ( !defined $align ) {
            printf "Can't find align for %s\n", $job->{_id};
            next;
        }

        printf "Process align %s:%s-%s\n", $align->{chr}{name}, $align->{chr}{start},
            $align->{chr}{end};

        # all beds in this align
        my @beds = $coll_bed->find( { 'align._id' => $align->{_id} } )->all;
        next unless @beds;
        printf "    %d beds in this align\n", scalar @beds;
        my $bed_chr_set = AlignDB::IntSpan->new;
        for (@beds) {
            $bed_chr_set->add_runlist( $_->{chr}{runlist} );
        }

        # all gsws in this align
        my @gsws = $coll_gsw->find( { 'align._id' => $align->{_id} } )->all;
        my %gsw_count_of;
        for my $gsw (@gsws) {
            next if $bed_chr_set->intersect( $gsw->{chr}{runlist} )->is_empty;

            my $count = count_bed_in_sw( $coll_bed, $align->{_id}, $gsw );

            if ($count) {
                $gsw_count_of{ $gsw->{_id} } = $count;
            }
            else {
                printf "gsw %s matching wrong\n", $gsw->{chr}{runlist};
            }
        }
        printf "    %d gsws in this align, %d overlapping with beds\n", scalar @gsws,
            scalar keys %gsw_count_of;

        # all ofgsw in this align
        my @ofgsws = $coll_ofgsw->find( { 'align._id' => $align->{_id} } )->all;
        my %ofgsw_count_of;
        for my $ofgsw (@ofgsws) {
            next if $bed_chr_set->intersect( $ofgsw->{chr}{runlist} )->is_empty;

            my $count = count_bed_in_sw( $coll_bed, $align->{_id}, $ofgsw );

            if ($count) {
                $ofgsw_count_of{ $ofgsw->{_id} } = $count;
            }
            else {
                printf "ofgsw %s matching wrong\n", $ofgsw->{chr}{runlist};
            }
        }
        printf "    %d ofgsws in this align, %d overlapping with beds\n", scalar @ofgsws,
            scalar keys %ofgsw_count_of;

        for my $key ( keys %gsw_count_of ) {
            $coll_gsw->update_one(
                { _id    => MongoDB::OID->new( value => $key ) },
                { '$set' => { bed_count              => $gsw_count_of{$key}, } },
            );
        }
        for my $key ( keys %ofgsw_count_of ) {
            $coll_ofgsw->update_one(
                { _id    => MongoDB::OID->new( value => $key ) },
                { '$set' => { bed_count              => $ofgsw_count_of{$key}, } },
            );
        }
    }
};

#----------------------------------------------------------#
# start update
#----------------------------------------------------------#
if ( $run eq "all" or $run eq "insert" ) {

    #@type MongoDB::Database
    my $db = MongoDB::MongoClient->new(
        host          => $server,
        port          => $port,
        query_timeout => -1,
    )->get_database($dbname);

    #@type MongoDB::Collection
    my $coll = $db->get_collection('bed');
    $coll->drop;

    my $mce = MCE->new( max_workers => $parallel, );
    $mce->foreach( \@files, $worker_insert, );    # foreach implies chunk_size => 1.

    #@type MongoDB::IndexView
    my $indexes = $coll->indexes;
    $indexes->create_one( [ 'align._id' => 1 ] );
    $indexes->create_one( [ chr_name => 1, chr_start => 1, chr_end => 1 ] );

    $stopwatch->block_message( check_coll( $db, 'bed', '_id' ) );
}

if ( $run eq "all" or $run eq "count" ) {
    #@type MongoDB::Database
    my $db = MongoDB::MongoClient->new(
        host          => $server,
        port          => $port,
        query_timeout => -1,
    )->get_database($dbname);

    my @jobs = $db->get_collection('align')->find->fields( { _id => 1 } )->all;
    $stopwatch->block_message("Total align: " . scalar(@jobs));

    my $mce = MCE->new( max_workers => $parallel, chunk_size => $batch_number, );
    $mce->forchunk( \@jobs, $worker_count, );

    $stopwatch->block_message( check_coll( $db, 'gsw',   'bed_count' ) );
    $stopwatch->block_message( check_coll( $db, 'ofgsw', 'bed_count' ) );
}

$stopwatch->end_message;

exit;

sub count_bed_in_sw {

    #@type MongoDB::Collection
    my $coll     = shift;
    my $align_id = shift;
    my $sw       = shift;

    my $count = $coll->find(
        {   'align._id' => $align_id,
            '$or'       => [

                # bed    |----|
                # sw  |----|
                {   'chr.start' => {
                        '$gte' => $sw->{chr}{start},
                        '$lte' => $sw->{chr}{end},
                    }
                },

                # bed |----|
                # sw    |----|
                {   'chr.end' => {
                        '$gte' => $sw->{chr}{start},
                        '$lte' => $sw->{chr}{end},
                    }
                },

                # bed |--------|
                # sw    |----|
                {   'chr.start' => { '$lte' => $sw->{chr}{start}, },
                    'chr.end'   => { '$gte' => $sw->{chr}{end}, }
                },

                # bed   |----|
                # sw  |--------|
                {   'chr.start' => { '$gte' => $sw->{chr}{start}, },
                    'chr.end'   => { '$lte' => $sw->{chr}{end}, }
                },
            ]
        }
    )->count;

    return $count;
}

__END__
