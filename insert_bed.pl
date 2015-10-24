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

use AlignDB::GC;
use AlignDB::IntSpan;
use AlignDB::Window;
use AlignDB::Stopwatch;

use lib "$FindBin::RealBin/lib";
use MyUtil qw(process_message check_coll center_resize calc_gc_ratio);

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

# AlignDB::GC options
my $wave_window_size  = $Config->{gc}{wave_window_size};
my $wave_window_step  = $Config->{gc}{wave_window_step};
my $vicinal_size      = $Config->{gc}{vicinal_size};
my $fall_range        = $Config->{gc}{fall_range};
my $gsw_size          = $Config->{gc}{gsw_size};
my $stat_segment_size = $Config->{gc}{stat_segment_size};
my $stat_window_size  = $Config->{gc}{stat_window_size};
my $stat_window_step  = $Config->{gc}{stat_window_step};

=cut

=head1 NAME

insert_bed.pl - Add bed files to ofg and generate ofgsw.

=head1 SYNOPSIS

    perl insert_bed.pl [options]
      Options:
        --help      -?          brief help message
        --server        STR     MongoDB server IP/Domain name
        --port          INT     MongoDB server port
        --db        -d  STR     database name
        --file      -f  @STR    bed files
        --tag       -f  @STR    bed tags
        --type      -f  @STR    bed types
        --style         STR     center_intact or center, default is [center_intact]
        --nochr                 remove 'chr' from chr_name
        --parallel      INT     run in parallel mode, [1]
        --batch         INT     number of alignments processed in one child
                                process, [10]

    mongo S288c --eval "db.dropDatabase();"
    perl gen_mg.pl -d S288c -n S288c --dir ~/data/alignment/yeast_combine/S288C  --parallel 1
    perl insert_bed.pl -d S288c --tag hot --type hot -f ~/Scripts/alignDB/ofg/spo11/spo11_hot.bed --batch 10 --parallel 1

=cut

GetOptions(
    'help|?' => sub { HelpMessage(0) },
    'server=s'   => \( my $server       = $Config->{database}{server} ),
    'port=i'     => \( my $port         = $Config->{database}{port} ),
    'db|d=s'     => \( my $dbname       = $Config->{database}{db} ),
    'file|f=s'   => \my @files,
    'tag=s'      => \my @tags,
    'type=s'     => \my @types,
    'style=s'    => \( my $style        = "center_intact" ),
    'nochr'      => \my $nochr,
    'parallel=i' => \( my $parallel     = 1 ),
    'batch=i'    => \( my $batch_number = 10 ),
) or HelpMessage(1);

#----------------------------------------------------------#
# Init
#----------------------------------------------------------#
$stopwatch->start_message("Insert bed to $dbname...");

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
        else {
            my $length          = $end - $start + 1;
            my $ofg_align_start = $start - $align->{chr}{start} + 1;
            my $ofg_align_end   = $end - $align->{chr}{start} + 1;
            my $ofg_seq         = substr $align->{seq}, $ofg_align_start - 1, $length;
            my $ofg_gc          = calc_gc_ratio($ofg_seq);
            push @beds,
                {
                align => {
                    _id     => $align->{_id},
                    start   => $ofg_align_start,
                    end     => $ofg_align_end,
                    runlist => AlignDB::IntSpan->new->add_range( $ofg_align_start, $ofg_align_end )
                        ->runlist,
                },
                chr => {
                    name    => $chr,
                    start   => $start,
                    end     => $end,
                    runlist => AlignDB::IntSpan->new->add_range( $start, $end )->runlist,
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
    my $coll_ofg = $db->get_collection('ofg');
    while ( scalar @beds ) {
        my @batching = splice @beds, 0, 10000;
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
    my $db = MongoDB::MongoClient->new(
        host          => $server,
        port          => $port,
        query_timeout => -1,
    )->get_database($dbname);

    my $coll_ofg   = $db->get_collection('ofg');
    my $coll_ofgsw = $db->get_collection('ofgsw');

    # AlignDB::GC
    my $obj = AlignDB::GC->new(
        wave_window_size => $wave_window_size,
        wave_window_step => $wave_window_step,
        vicinal_size     => $vicinal_size,
        fall_range       => $fall_range,
        gsw_size         => $gsw_size,
        stat_window_size => $stat_window_size,
        stat_window_step => $stat_window_step,
        skip_mdcw        => 1,
    );

    for my $job (@jobs) {
        my $align = process_message( $db, $job->{_id} );
        next unless $align;

        my @align_ofgs = $coll_ofg->find( { 'align._id' => $align->{_id} } )->all;
        if ( @align_ofgs == 0 ) {
            warn "No ofgs in this align\n";
            next;
        }
        printf "    Find %d ofgs in this align\n", scalar @align_ofgs;

        my $align_set = AlignDB::IntSpan->new->add_range( 1, $align->{length} );

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
                    = AlignDB::IntSpan->new->add_range( $ofgsw->{chr}{start}, $ofgsw->{chr}{end} )
                    ->runlist;
                $ofgsw->{align}{runlist}
                    = AlignDB::IntSpan->new->add_range( $ofgsw->{align}{start},
                    $ofgsw->{align}{end} )->runlist;

                # pre allocate
                $ofgsw->{bed_count} = 0;
                my $ofgsw_seq = substr $align->{seq}, $rsw->{set}->min - 1, $ofgsw->{length};
                $ofgsw->{gc} = {
                    gc   => calc_gc_ratio($ofgsw_seq),
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
# insert bed files
{
    my $db = MongoDB::MongoClient->new(
        host          => $server,
        port          => $port,
        query_timeout => -1,
    )->get_database($dbname);
    my $name = "ofg";
    my $coll = $db->get_collection($name);
    $coll->drop;

    my @args = zip @files, @tags, @types;
    my @jobs;
    while ( scalar @args ) {
        my @batching = splice @args, 0, 3;
        push @jobs, [@batching];
    }

    my $mce = MCE->new( max_workers => $parallel, );
    $mce->foreach( \@jobs, $worker_insert, );    # foreach implies chunk_size => 1.

    $stopwatch->block_message("Indexing $name");
    my $indexes = $coll->indexes;
    $indexes->create_one( [ 'align._id' => 1 ] );
    $indexes->create_one( [ 'chr.name'  => 1, 'chr.start' => 1, 'chr.end' => 1 ] );
    $indexes->create_one( [ tag         => 1 ] );
    $indexes->create_one( [ type        => 1 ] );

    $stopwatch->block_message( check_coll( $db, $name, '_id' ) );
}

# ofgsw
{
    my $db = MongoDB::MongoClient->new(
        host          => $server,
        port          => $port,
        query_timeout => -1,
    )->get_database($dbname);

    # get aligns
    my @jobs = $db->get_collection('align')->find->fields( { _id => 1 } )->all;

    # insert ofgsw
    my $name = "ofgsw";
    my $coll = $db->get_collection($name);
    $coll->drop;

    my $mce = MCE->new( max_workers => $parallel, chunk_size => $batch_number, );
    $mce->forchunk( \@jobs, $worker_sw, );

    # indexing
    $stopwatch->block_message("Indexing $name");
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
