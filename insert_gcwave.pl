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

use AlignDB::GC;
use AlignDB::IntSpan;
use AlignDB::Stopwatch;

use lib "$FindBin::RealBin/lib";
use MyUtil qw(center_resize check_coll calc_gc_ratio);

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

=head1 NAME

insert_gcwave.pl - Add GC ralated tables to alignDB

=head1 SYNOPSIS

    perl insert_gcwave.pl [options]
      Options:
        --help      -?          brief help message
        --server        STR     MongoDB server IP/Domain name
        --port          INT     MongoDB server port
        --db        -d  STR     database name
        --parallel      INT     run in parallel mode, [1]
        --batch         INT     number of alignments processed in one child
                                process, [10]

    # S288c
    mongo S288c --eval "db.dropDatabase();"
    
    # Runtime 5 seconds.
    perl gen_mg.pl -d S288c -n S288c --dir ~/data/alignment/yeast_combine/S288C  --parallel 1
    
    # Runtime 38 seconds.
    perl insert_gcwave.pl -d S288c --batch 1 --parallel 4
    
    # Runtime 21 seconds.
    perl update_sw_cv.pl -d S288c --batch 1 --parallel 4

=cut

GetOptions(
    'help|?' => sub { HelpMessage(0) },
    'server=s'   => \( my $server       = $Config->{database}{server} ),
    'port=i'     => \( my $port         = $Config->{database}{port} ),
    'db|d=s'     => \( my $dbname       = $Config->{database}{db} ),
    'parallel=i' => \( my $parallel     = 1 ),
    'batch=i'    => \( my $batch_number = 10 ),
) or HelpMessage(1);

#----------------------------------------------------------#
# Init
#----------------------------------------------------------#
$stopwatch->start_message("Update GC tables of $dbname...");

# retrieve all _id from align
my @jobs;
{
    my $db = MongoDB::MongoClient->new(
        host          => $server,
        port          => $port,
        query_timeout => -1,
    )->get_database($dbname);

    @jobs = $db->get_collection('align')->find->fields( { _id => 1 } )->all;
    printf "There are [%d] aligns totally.\n", scalar @jobs;

    $db->get_collection('gce')->drop;
    $db->get_collection('gsw')->drop;
}

#----------------------------------------------------------#
# Worker
#----------------------------------------------------------#
my $worker = sub {
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

    # mocking AlignDB::GC
    my $obj_gc = AlignDB::GC->new(
        wave_window_size => $wave_window_size,
        wave_window_step => $wave_window_step,
        vicinal_size     => $vicinal_size,
        fall_range       => $fall_range,
        gsw_size         => $gsw_size,
        stat_window_size => $stat_window_size,
        stat_window_step => $stat_window_step,
        skip_mdcw        => 1,
    );

    my $coll_align = $db->get_collection('align');
    for my $job (@jobs) {
        my $align = $coll_align->find_one( { _id => $job->{_id} } );
        if ( !defined $align ) {
            printf "Can't find align for %s\n", $job->{_id};
            next;
        }

        printf "Process align %s:%s-%s\n", $align->{chr}{name}, $align->{chr}{start},
            $align->{chr}{end};

        print "    Insert gc extremes\n";
        insert_gce( $obj_gc, $db, $align );

        print "    Insert gc sliding windows\n";
        insert_gsw( $obj_gc, $db, $align );
    }

    $inner_watch->block_message( "Task [$chunk_id] has been processed.", "duration" );

    return;
};

#----------------------------------------------------------#
# start update
#----------------------------------------------------------#
my $mce = MCE->new( max_workers => $parallel, chunk_size => $batch_number, );
$mce->forchunk( \@jobs, $worker, );

#----------------------------#
# index and check
#----------------------------#
{
    my $db = MongoDB::MongoClient->new(
        host          => $server,
        port          => $port,
        query_timeout => -1,
    )->get_database($dbname);

    {
        my $name = "gce";
        $stopwatch->block_message("Indexing $name");
        my $coll    = $db->get_collection($name);
        my $indexes = $coll->indexes;
        $indexes->create_one( [ 'align._id' => 1 ] );
        $indexes->create_one( [ prev_id     => 1 ] );
        $indexes->create_one( [ "chr.name"  => 1, "chr.start" => 1, "chr.end" => 1 ] );
        $indexes->create_one( [ type        => 1 ] );
    }

    {
        my $name = "gsw";
        $stopwatch->block_message("Indexing $name");
        my $coll    = $db->get_collection($name);
        my $indexes = $coll->indexes;
        $indexes->create_one( [ 'align._id'          => 1 ] );
        $indexes->create_one( [ 'gce._id'            => 1 ] );
        $indexes->create_one( [ 'gce.prev_id'        => 1 ] );
        $indexes->create_one( [ "chr.name"           => 1, "chr.start" => 1, "chr.end" => 1 ] );
        $indexes->create_one( [ type                 => 1 ] );
        $indexes->create_one( [ 'gce.distance'       => 1 ] );
        $indexes->create_one( [ 'gce.distance_crest' => 1 ] );
        $indexes->create_one( [ 'gce.gradient'       => 1 ] );
    }

    $stopwatch->block_message( check_coll( $db, 'gce', '_id' ) );
    $stopwatch->block_message( check_coll( $db, 'gsw', 'gc.cv' ) );
}

$stopwatch->end_message;

exit;

#----------------------------------------------------------#
# Subroutines
#----------------------------------------------------------#
sub insert_gce {
    my $obj_gc = shift;
    my $db     = shift;
    my $align  = shift;

    my $comparable_set = AlignDB::IntSpan->new( "1-" . $align->{length} );
    my @gce_site;
    my @slidings = $obj_gc->gc_wave( [ $align->{seq} ], $comparable_set );
    for my $s (@slidings) {
        my $flag = $s->{high_low_flag};

        # rename high_low_flag to type
        if ( $flag eq 'T' or $flag eq 'C' ) {
            $s->{type} = $flag;
            delete $s->{high_low_flag};
            push @gce_site, $s;
        }
    }

    # get extreme sliding windows' sizes
    my $windows_size = $obj_gc->wave_window_size;
    my $half_length  = int( $windows_size / 2 );

    # left
    my $prev_gce_middle_right_idx = 1;
    my $prev_gce_gc               = $slidings[0]->{gc};
    for ( my $i = 0; $i < scalar @gce_site; $i++ ) {

        # wave_length
        my $gce_set             = $gce_site[$i]->{set};
        my $gce_middle_left     = $gce_set->at($half_length);
        my $gce_middle_left_idx = $comparable_set->index($gce_middle_left);
        my $left_wave_length    = $gce_middle_left_idx - $prev_gce_middle_right_idx + 1;
        $prev_gce_middle_right_idx = $gce_middle_left_idx + 1;
        $gce_site[$i]->{left_wave_length} = $left_wave_length;

        # amplitude
        my $gce_gc         = $gce_site[$i]->{gc};
        my $left_amplitude = abs( $gce_gc - $prev_gce_gc );
        $gce_site[$i]->{left_amplitude} = $left_amplitude;
        $prev_gce_gc = $gce_gc;
    }

    # right
    my $next_gce_middle_left_idx = $comparable_set->cardinality;
    my $next_gce_gc              = $slidings[-1]->{gc};
    for ( my $i = scalar @gce_site - 1; $i >= 0; $i-- ) {

        # wave_length
        my $gce_set              = $gce_site[$i]->{set};
        my $gce_middle_right     = $gce_set->at( $half_length + 1 );
        my $gce_middle_right_idx = $comparable_set->index($gce_middle_right);
        my $right_wave_length    = $next_gce_middle_left_idx - $gce_middle_right_idx + 1;
        $next_gce_middle_left_idx = $gce_middle_right_idx - 1;
        $gce_site[$i]->{right_wave_length} = $right_wave_length;

        # amplitude
        my $gce_gc          = $gce_site[$i]->{gc};
        my $right_amplitude = abs( $gce_gc - $next_gce_gc );
        $gce_site[$i]->{right_amplitude} = $right_amplitude;
        $prev_gce_gc = $gce_gc;
    }

    my $coll        = $db->get_collection('gce');
    my $prev_gce_id = undef;
    for my $exs (@gce_site) {
        my $chr_start = $exs->{set}->min + $align->{chr}{start} - 1;
        my $chr_end   = $exs->{set}->max + $align->{chr}{start} - 1;
        $exs->{chr} = {
            common_name => $align->{chr}{common_name},
            name        => $align->{chr}{name},
            start       => $chr_start,
            end         => $chr_end,
            strand      => $align->{chr}{strand},
            runlist     => $chr_start . '-' . $chr_end,
        };
        $exs->{align} = {
            _id     => $align->{_id},
            start   => $exs->{set}->min,
            end     => $exs->{set}->max,
            runlist => $exs->{set}->runlist,
        };
        $exs->{length} = $exs->{set}->cardinality;
        delete $exs->{set};

        $exs->{prev_id} = $prev_gce_id;
        $prev_gce_id = $coll->insert_one($exs)->inserted_id;
    }

    return;
}

sub insert_gsw {
    my $obj_gc = shift;
    my $db     = shift;
    my $align  = shift;

    my $comparable_set = AlignDB::IntSpan->new( "1-" . $align->{length} );

    my $coll_gce = $db->get_collection('gce');
    my $coll_gsw = $db->get_collection('gsw');

    # get gc sliding windows' sizes
    my $gsw_size = $obj_gc->gsw_size;

    my @gces = $coll_gce->find( { "align._id" => $align->{_id} } )->all;
    for my $ex (@gces) {

        # bypass the first gce
        next unless $ex->{prev_id};

        # get attrs of the two gces
        my $ex_type             = $ex->{type};
        my $ex_gc               = $ex->{gc};
        my $ex_left_amplitude   = $ex->{left_amplitude};
        my $ex_left_wave_length = $ex->{left_wave_length};
        my $ex_set              = AlignDB::IntSpan->new( $ex->{align}{runlist} );

        my $prev_ex     = $coll_gce->find_one( { _id => $ex->{prev_id} } );
        my $prev_ex_gc  = $prev_ex->{gc};
        my $prev_ex_set = AlignDB::IntSpan->new( $prev_ex->{align}{runlist} );

        # determining $gsw_density here, which is different from isw_density
        my $half_length = int( $ex_set->cardinality / 2 );
        my $gsw_density = int( ( $ex_left_wave_length - $half_length ) / $gsw_size );

        # wave length, amplitude, trough_gc and gradient
        my $gsw_wave_length = $ex_left_wave_length;
        my $gsw_amplitude   = $ex_left_amplitude;
        my $gsw_trough_gc   = $ex_type eq 'T' ? $ex_gc : $prev_ex_gc;
        my $gsw_crest_gc    = $ex_type eq 'T' ? $prev_ex_gc : $ex_gc;
        my $gsw_gradient    = int( $gsw_amplitude / $ex_left_wave_length / 0.00001 );

        # determining $gsw_type here, ascend and descent
        my $gsw_type;
        my @gsw_windows;
        if ( $ex_type eq 'T' ) {    # push trough to gsw
            $gsw_type = 'D';        # descend, left of trough£¬right of crest
            push @gsw_windows,
                {
                type => 'T',
                set  => $ex_set,
                gce  => {
                    distance       => 0,
                    distance_crest => $gsw_density + 1,
                },
                };
        }
        elsif ( $ex_type eq 'C' ) {    # crest
            $gsw_type = 'A';           # ascend, right of trough, left of crest
            push @gsw_windows,
                {
                type => 'C',
                set  => $ex_set,
                gce  => {
                    distance       => $gsw_density + 1,
                    distance_crest => 0,
                },
                };
        }
        else {
            warn "gce_type [$ex_type] error\n";
        }

        {    # More windows will be submitted in the following section

            # distance start counting from trough
            # $sw_start and $sw_end are both index of $comprarable_set
            # $gsw_distance is from 1 to $gsw_density
            # window 0 is trough
            # ..., D2, D1, T0, A1, A2, ...
            my ( $sw_start, $sw_end );
            if ( $gsw_type eq 'A' ) {
                $sw_start = $comparable_set->index( $prev_ex_set->max ) + 1;
                $sw_end   = $sw_start + $gsw_size - 1;
            }
            elsif ( $gsw_type eq 'D' ) {
                $sw_end   = $comparable_set->index( $ex_set->min ) - 1;
                $sw_start = $sw_end - $gsw_size + 1;
            }

            for my $gsw_distance ( 1 .. $gsw_density ) {
                my $gsw_set = $comparable_set->slice( $sw_start, $sw_end );

                push @gsw_windows,
                    {
                    type => $gsw_type,
                    set  => $gsw_set,
                    gce  => {
                        distance       => $gsw_distance,
                        distance_crest => $gsw_density - $gsw_distance + 1,
                    }
                    };

                if ( $gsw_type eq 'A' ) {
                    $sw_start = $sw_end + 1;
                    $sw_end   = $sw_start + $gsw_size - 1;
                }
                elsif ( $gsw_type eq 'D' ) {
                    $sw_end   = $sw_start - 1;
                    $sw_start = $sw_end - $gsw_size + 1;
                }
            }

            for my $gsw (@gsw_windows) {
                my $chr_start = $gsw->{set}->min + $align->{chr}{start} - 1;
                my $chr_end   = $gsw->{set}->max + $align->{chr}{start} - 1;

                $gsw->{chr} = {
                    common_name => $align->{chr}{common_name},
                    name        => $align->{chr}{name},
                    start       => $chr_start,
                    end         => $chr_end,
                    strand      => $align->{chr}{strand},
                    runlist     => $chr_start . '-' . $chr_end,
                };
                $gsw->{align} = {
                    _id     => $align->{_id},
                    start   => $gsw->{set}->min,
                    end     => $gsw->{set}->max,
                    runlist => $gsw->{set}->runlist,
                };

                $gsw->{length} = $gsw->{set}->cardinality;

                # wave properties
                $gsw->{gce} = {
                    _id         => $ex->{_id},
                    prev_id     => $prev_ex->{_id},
                    wave_length => $gsw_wave_length,
                    amplitude   => $gsw_amplitude,
                    trough_gc   => $gsw_trough_gc,
                    crest_gc    => $gsw_crest_gc,
                    gradient    => $gsw_gradient,
                    %{ $gsw->{gce} },
                };

                # pre allocate
                $gsw->{bed_count} = 0;
                my $gsw_seq = substr $align->{seq}, $gsw->{set}->min - 1, $gsw->{length};
                $gsw->{gc} = {
                    gc   => calc_gc_ratio($gsw_seq),
                    mean => 0.0,
                    cv   => 0.0,
                    std  => 0.0,
                };

                delete $gsw->{set};
            }
            $coll_gsw->insert_many( \@gsw_windows );
        }
    }

    return;
}

__END__

