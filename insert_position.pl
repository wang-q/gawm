#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Config::Tiny;
use FindBin;
use YAML::Syck;

use MCE;
use MCE::Flow Sereal => 1;

use MongoDB;
$MongoDB::BSON::looks_like_number = 1;
use MongoDB::OID;

use List::MoreUtils;

use AlignDB::IntSpan;
use AlignDB::Stopwatch;
use AlignDB::Window;

use App::RL::Common;
use App::Fasops::Common;

use lib "$FindBin::RealBin/lib";
use MyUtil qw(process_message check_coll);

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

insert_position.pl - Add position files to ofg and generate ofgsw.

=head1 SYNOPSIS

    perl insert_position.pl [options]
      Options:
        --help      -?          brief help message
        --server        STR     MongoDB server IP/Domain name
        --port          INT     MongoDB server port
        --db        -d  STR     database name
        --file      -f  @STR    bed files
        --tag       -f  @STR    bed tags
        --type      -f  @STR    bed types
        --style         STR     center_intact or center, default is [center_intact]
        --parallel      INT     run in parallel mode, [1]
        --batch         INT     alignments processed in one child process, [10]

    mongo S288c --eval "db.dropDatabase();"
    perl gen_mg.pl -d S288c -n S288c --dir ~/data/alignment/yeast_combine/S288C  --parallel 1
    perl insert_bed.pl -d S288c --tag hot --type hot -f ~/Scripts/alignDB/ofg/spo11/spo11_hot.bed --batch 10 --parallel 1

=cut

GetOptions(
    'help|?' => sub { Getopt::Long::HelpMessage(0) },
    'server=s'   => \( my $server       = $Config->{database}{server} ),
    'port=i'     => \( my $port         = $Config->{database}{port} ),
    'db|d=s'     => \( my $dbname       = $Config->{database}{db} ),
    'file|f=s'   => \my @files,
    'tag=s'      => \my @tags,
    'type=s'     => \my @types,
    'style=s'    => \( my $style        = "center_intact" ),
    'parallel=i' => \( my $parallel     = 1 ),
    'batch=i'    => \( my $batch_number = 10 ),
) or Getopt::Long::HelpMessage(1);

#----------------------------------------------------------#
# Init
#----------------------------------------------------------#
$stopwatch->start_message("Insert positions to $dbname...");

#----------------------------------------------------------#
# workers
#----------------------------------------------------------#
my $worker_insert = sub {
    my ( $self, $chunk_ref, $chunk_id ) = @_;
    my $job = $chunk_ref->[0];

    my $wid = MCE->wid;

    my $inner_watch = AlignDB::Stopwatch->new;
    $inner_watch->block_message("Process task [$chunk_id] by worker #$wid");

    my ( $file, $tag, $type ) = @$job;
    print "Reading file [$file]\n";

    # wait forever for responses
    #@type MongoDB::Database
    my $db = MongoDB::MongoClient->new(
        host          => $server,
        port          => $port,
        query_timeout => -1,
    )->get_database($dbname);

    #@type MongoDB::Collection
    my $coll_align = $db->get_collection('align');

    my @data;
    open my $data_fh, '<', $file;
    while ( my $string = <$data_fh> ) {
        next unless defined $string;
        chomp $string;

        my $info = App::RL::Common::decode_header($string);
        next unless defined $info->{chr};
        $info->{tag}  = $tag;
        $info->{type} = $type;

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
        else {
            my $length          = $info->{end} - $info->{start} + 1;
            my $ofg_align_start = $info->{start} - $align->{chr}{start} + 1;
            my $ofg_align_end   = $info->{end} - $align->{chr}{start} + 1;
            my $ofg_seq         = substr $align->{seq}, $ofg_align_start - 1, $length;
            my $ofg_gc          = App::Fasops::Common::calc_gc_ratio( [$ofg_seq] );
            push @data,
                {
                align => {
                    _id     => $align->{_id},
                    start   => $ofg_align_start,
                    end     => $ofg_align_end,
                    runlist => AlignDB::IntSpan->new->add_pair( $ofg_align_start, $ofg_align_end )
                        ->runlist,
                },
                chr => {
                    name  => $info->{chr},
                    start => $info->{start},
                    end   => $info->{end},
                    runlist =>
                        AlignDB::IntSpan->new->add_pair( $info->{start}, $info->{end} )->runlist,
                },
                length => $length,
                gc     => $ofg_gc,
                tag    => $tag,
                type   => $type,
                };
        }
    }
    close $data_fh;

    print "Inserting file [$file]\n";

    #@type MongoDB::Collection
    my $coll_ofg = $db->get_collection('ofg');
    while ( scalar @data ) {
        my @batching = splice @data, 0, 10000;
        $coll_ofg->insert_many( \@batching );
    }
    print "Insert done.\n";
};

my $worker_sw = sub {
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

    #@type MongoDB::Collection
    my $coll_ofg = $db->get_collection('ofg');

    #@type MongoDB::Collection
    my $coll_ofgsw = $db->get_collection('ofgsw');

    for my $job (@jobs) {
        my $align = process_message( $db, $job->{_id} );
        next unless $align;

        my @align_ofgs = $coll_ofg->find( { 'align._id' => $align->{_id} } )->all;
        if ( @align_ofgs == 0 ) {
            warn "No ofgs in this align\n";
            next;
        }
        printf "    Find %d ofgs in this align\n", scalar @align_ofgs;

        my $align_set = AlignDB::IntSpan->new->add_pair( 1, $align->{length} );

        #----------------------------#
        # ofgsw
        #----------------------------#
        my $window_maker = AlignDB::Window->new(
            sw_size          => 100,
            max_out_distance => 20,
            max_in_distance  => 20,
        );

        for my $ofg (@align_ofgs) {
            my @rsws;
            if ( $style eq 'center_intact' ) {
                @rsws = $window_maker->center_intact_window(
                    $align_set,
                    $ofg->{align}{start},
                    $ofg->{align}{end}
                );
            }
            elsif ( $style eq 'center' ) {
                @rsws = $window_maker->center_window(
                    $align_set,
                    $ofg->{align}{start},
                    $ofg->{align}{end}
                );
            }

            my @ofgsws;
            for my $rsw (@rsws) {
                my $ofgsw = {
                    chr => {
                        name  => $align->{chr}{name},
                        start => $rsw->{set}->min + $align->{chr}{start} - 1,
                        end   => $rsw->{set}->max + $align->{chr}{start} - 1,
                    },
                    align => {
                        _id   => $align->{_id},
                        start => $rsw->{set}->min,
                        end   => $rsw->{set}->max,
                    },
                    ofg => {
                        _id      => $ofg->{_id},
                        tag      => $ofg->{tag},
                        type     => $ofg->{type},
                        distance => $rsw->{distance},
                    },
                    type   => $rsw->{type},
                    length => $rsw->{set}->size,
                };
                $ofgsw->{chr}{runlist}
                    = AlignDB::IntSpan->new->add_pair( $ofgsw->{chr}{start}, $ofgsw->{chr}{end} )
                    ->runlist;
                $ofgsw->{align}{runlist}
                    = AlignDB::IntSpan->new->add_pair( $ofgsw->{align}{start},
                    $ofgsw->{align}{end} )->runlist;

                # pre allocate
                $ofgsw->{bed_count} = 0;
                my $ofgsw_seq = substr $align->{seq}, $rsw->{set}->min - 1, $ofgsw->{length};
                $ofgsw->{gc} = {
                    gc => App::Fasops::Common::calc_gc_ratio( [$ofgsw_seq] ),
                    mean => 0.0,
                    cv   => 0.0,
                    std  => 0.0,
                };

                push @ofgsws, $ofgsw;
            }
            $coll_ofgsw->insert_many( \@ofgsws );
        }
    }
};

#----------------------------------------------------------#
# start update
#----------------------------------------------------------#
# insert position files
{
    #@type MongoDB::Database
    my $db = MongoDB::MongoClient->new(
        host          => $server,
        port          => $port,
        query_timeout => -1,
    )->get_database($dbname);
    my $name = "ofg";

    #@type MongoDB::Collection
    my $coll = $db->get_collection($name);
    $coll->drop;

    my @args = List::MoreUtils::PP::mesh @files, @tags, @types;
    my @jobs;
    while ( scalar @args ) {
        my @batching = splice @args, 0, 3;
        push @jobs, [@batching];
    }

    my $mce = MCE->new( max_workers => $parallel, );
    $mce->foreach( \@jobs, $worker_insert, );    # foreach implies chunk_size => 1.

    $stopwatch->block_message("Indexing $name");

    #@type MongoDB::IndexView
    my $indexes = $coll->indexes;
    $indexes->create_one( [ 'align._id' => 1 ] );
    $indexes->create_one( [ 'chr.name'  => 1, 'chr.start' => 1, 'chr.end' => 1 ] );
    $indexes->create_one( [ tag         => 1 ] );
    $indexes->create_one( [ type        => 1 ] );

    $stopwatch->block_message( check_coll( $db, $name, '_id' ) );
}

# ofgsw
{
    #@type MongoDB::Database
    my $db = MongoDB::MongoClient->new(
        host          => $server,
        port          => $port,
        query_timeout => -1,
    )->get_database($dbname);

    # get aligns
    my @jobs = $db->get_collection('align')->find->fields( { _id => 1 } )->all;

    # insert ofgsw
    my $name = "ofgsw";

    #@type MongoDB::Collection
    my $coll = $db->get_collection($name);
    $coll->drop;

    my $mce = MCE->new( max_workers => $parallel, chunk_size => $batch_number, );
    $mce->forchunk( \@jobs, $worker_sw, );

    # indexing
    $stopwatch->block_message("Indexing $name");

    #@type MongoDB::IndexView
    my $indexes = $coll->indexes;
    $indexes->create_one( [ 'align._id'    => 1 ] );
    $indexes->create_one( [ 'ofg._id'      => 1 ] );
    $indexes->create_one( [ 'chr.name'     => 1, 'chr.start' => 1, 'chr.end' => 1 ] );
    $indexes->create_one( [ type           => 1 ] );
    $indexes->create_one( [ 'ofg.distance' => 1 ] );
    $indexes->create_one( { 'ofg.tag'  => 1 } );
    $indexes->create_one( { 'ofg.type' => 1 } );

    $stopwatch->block_message( check_coll( $db, $name, '_id' ) );
}

$stopwatch->end_message;

exit;

__END__
