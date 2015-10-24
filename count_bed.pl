#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw(HelpMessage);
use Config::Tiny;
use FindBin;
use YAML qw(Dump Load DumpFile LoadFile);

use MCE;
use MongoDB;
$MongoDB::BSON::looks_like_number = 1;
use MongoDB::OID;

use Roman;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use List::MoreUtils qw(any all uniq zip);

use AlignDB::IntSpan;
use AlignDB::Window;
use AlignDB::Stopwatch;

use lib "$FindBin::RealBin/lib";
use MyUtil qw(check_coll);

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

count_bed.pl - Add bed files to bed and count intersections.

=head1 SYNOPSIS

    perl count_bed.pl [options]
      Options:
        --help      -?          brief help message
        --server        STR     MongoDB server IP/Domain name
        --port          INT     MongoDB server port
        --db        -d  STR     database name
        --file      -f  @STR    bed files
        --run       -r  STR     all, insert or count, default is [all]
        --nochr                 remove 'chr' from chr_name
        --parallel      INT     run in parallel mode, [1]
        --batch         INT     number of alignments processed in one child
                                process, [10]

    perl count_bed.pl -d S288c --run insert -f ~/Scripts/alignDB/ofg/spo11/spo11_hot.bed --batch 1 --parallel 1
    perl count_bed.pl -d S288c --run count --batch 1 --parallel 4

=cut

GetOptions(
    'help|?' => sub { HelpMessage(0) },
    'server=s'   => \( my $server       = $Config->{database}{server} ),
    'port=i'     => \( my $port         = $Config->{database}{port} ),
    'db|d=s'     => \( my $dbname       = $Config->{database}{db} ),
    'f|file=s'   => \my @files,
    'r|run=s'    => \( my $run          = "all" ),
    'nochr'      => \my $nochr,
    'parallel=i' => \( my $parallel     = 1 ),
    'batch=i'    => \( my $batch_number = 10 ),
) or HelpMessage(1);

#----------------------------------------------------------#
# Init objects
#----------------------------------------------------------#
$stopwatch->start_message("Count bed of $dbname...");

$stopwatch->block_message( "Total bed files are " . scalar(@files) ) if @files;

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
        my ( $chr, $start, $end )
            = ( split /\t/, $string )[ 0, 1, 2 ];
        next unless $chr =~ /^\w+$/;
        if ( !$nochr ) {
            $chr =~ s/chr0?//i;
            $chr = "chr$chr";
        }
        next unless $start =~ /^\d+$/;
        next unless $end =~ /^\d+$/;

        if ( $start > $end ) {
            ( $start, $end ) = ( $end, $start );
        }

        my $align = $coll_align->find_one(
            {   'chr.name'  => $chr,
                'chr.start' => { '$lte' => $start },
                'chr.end'   => { '$gte' => $end }
            }
        );
        if ( !$align ) {
            print "    Can't locate an align for $chr:$start-$end\n";
            next;
        }
        push @beds,
            {
            align => { _id => $align->{_id}, },
            chr   => {
                name    => $chr,
                start   => $start,
                end     => $end,
                runlist => AlignDB::IntSpan->new->add_range( $start, $end )->runlist,
            },
            };
    }
    close $data_fh;

    print "Inserting file [$file]\n";
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
    my $db = MongoDB::MongoClient->new(
        host          => $server,
        port          => $port,
        query_timeout => -1,
    )->get_database($dbname);
    my $coll = $db->get_collection('bed');
    $coll->drop;

    my $mce = MCE->new( max_workers => $parallel, );
    $mce->foreach( \@files, $worker_insert, );    # foreach implies chunk_size => 1.

    my $indexes = $coll->indexes;
    $indexes->create_one( [ 'align._id' => 1 ] );
    $indexes->create_one( [ chr_name => 1, chr_start => 1, chr_end => 1 ] );

    $stopwatch->block_message( check_coll( $db, 'bed', '_id' ) );
}

if ( $run eq "all" or $run eq "count" ) {
    my $db = MongoDB::MongoClient->new(
        host          => $server,
        port          => $port,
        query_timeout => -1,
    )->get_database($dbname);

    my @jobs = $db->get_collection('align')->find->fields( { _id => 1 } )->all;

    my $mce = MCE->new( max_workers => $parallel, chunk_size => $batch_number, );
    $mce->forchunk( \@jobs, $worker_count, );

    $stopwatch->block_message( check_coll( $db, 'gsw',   'bed_count' ) );
    $stopwatch->block_message( check_coll( $db, 'ofgsw', 'bed_count' ) );
}

$stopwatch->end_message;

exit;

sub count_bed_in_sw {
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
