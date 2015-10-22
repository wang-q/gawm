#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw(HelpMessage);
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use MCE;
use MongoDB;
$MongoDB::BSON::looks_like_number = 1;
use MongoDB::OID;

use AlignDB::GC;
use AlignDB::IntSpan;
use AlignDB::Stopwatch;

use FindBin;
use lib "$FindBin::RealBin/lib";
use MyUtil qw(center_resize check_coll calc_gc_ratio);

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->read("$FindBin::Bin/config.ini");

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

    @jobs = $db->get_collection('align')->find->all;
    printf "There are [%d] aligns totally.\n", scalar @jobs;

    $db->get_collection('extreme')->drop;
    $db->get_collection('gsw')->drop;
}

#----------------------------------------------------------#
# Worker
#----------------------------------------------------------#
my $worker = sub {
    my ( $self, $chunk_ref, $chunk_id ) = @_;
    my @aligns = @{$chunk_ref};

    my $wid = MCE->wid;

    my $inner_watch = AlignDB::Stopwatch->new;
    $inner_watch->block_message("Process task [$chunk_id] by worker #$wid");

    # wait forever for responses
    my $db = MongoDB::MongoClient->new(
        host          => $server,
        port          => $port,
        query_timeout => -1,
    )->get_database($dbname);

    my $coll_seq = $db->get_collection('sequence');

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

    for my $align (@aligns) {
        printf "Process align %s:%s-%s\n", $align->{chr_name}, $align->{chr_start}, $align->{chr_end};

        my $seq = $coll_seq->find_one( { _id => $align->{seq_id} } )->{seq};

        print "    Insert gc extremes\n";
        insert_extreme( $obj_gc, $db, $align, $seq );

        print "    Insert gc sliding windows\n";
        insert_gsw( $obj_gc, $db, $align, $seq );
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
        my $name = "extreme";
        $stopwatch->block_message("Indexing $name");
        my $coll    = $db->get_collection($name);
        my $indexes = $coll->indexes;
        $indexes->create_one( [ align_id => 1 ] );
        $indexes->create_one( [ prev_id  => 1 ] );
        $indexes->create_one( [ chr_name => 1, chr_start => 1, chr_end => 1 ] );
        $indexes->create_one( [ type     => 1 ] );
    }

    {
        my $name = "gsw";
        $stopwatch->block_message("Indexing $name");
        my $coll    = $db->get_collection($name);
        my $indexes = $coll->indexes;
        $indexes->create_one( [ align_id        => 1 ] );
        $indexes->create_one( [ extreme_id      => 1 ] );
        $indexes->create_one( [ prev_extreme_id => 1 ] );
        $indexes->create_one( [ chr_name        => 1, chr_start => 1, chr_end => 1 ] );
        $indexes->create_one( [ type            => 1 ] );
        $indexes->create_one( [ distance        => 1 ] );
        $indexes->create_one( [ distance_crest  => 1 ] );
        $indexes->create_one( [ gradient        => 1 ] );
    }

    $stopwatch->block_message( check_coll( $db, 'extreme', '_id' ) );
    $stopwatch->block_message( check_coll( $db, 'gsw',     '_id' ) );
}

$stopwatch->end_message;

exit;

#----------------------------------------------------------#
# Subroutines
#----------------------------------------------------------#
sub insert_extreme {
    my $obj_gc = shift;
    my $db     = shift;
    my $align  = shift;
    my $seq    = shift;

    my $comparable_set = AlignDB::IntSpan->new( "1-" . $align->{length} );
    my @extreme_site;
    my @slidings = $obj_gc->gc_wave( [$seq], $comparable_set );
    for my $s (@slidings) {
        my $flag = $s->{high_low_flag};

        # rename high_low_flag to type
        if ( $flag eq 'T' or $flag eq 'C' ) {
            $s->{type} = $flag;
            delete $s->{high_low_flag};
            push @extreme_site, $s;
        }
    }

    # get extreme sliding windows' sizes
    my $windows_size = $obj_gc->wave_window_size;
    my $half_length  = int( $windows_size / 2 );

    # left
    my $prev_extreme_middle_right_idx = 1;
    my $prev_extreme_gc               = $slidings[0]->{gc};
    for ( my $i = 0; $i < scalar @extreme_site; $i++ ) {

        # wave_length
        my $extreme_set             = $extreme_site[$i]->{set};
        my $extreme_middle_left     = $extreme_set->at($half_length);
        my $extreme_middle_left_idx = $comparable_set->index($extreme_middle_left);
        my $left_wave_length        = $extreme_middle_left_idx - $prev_extreme_middle_right_idx + 1;
        $prev_extreme_middle_right_idx = $extreme_middle_left_idx + 1;
        $extreme_site[$i]->{left_wave_length} = $left_wave_length;

        # amplitude
        my $extreme_gc     = $extreme_site[$i]->{gc};
        my $left_amplitude = abs( $extreme_gc - $prev_extreme_gc );
        $extreme_site[$i]->{left_amplitude} = $left_amplitude;
        $prev_extreme_gc = $extreme_gc;
    }

    # right
    my $next_extreme_middle_left_idx = $comparable_set->cardinality;
    my $next_extreme_gc              = $slidings[-1]->{gc};
    for ( my $i = scalar @extreme_site - 1; $i >= 0; $i-- ) {

        # wave_length
        my $extreme_set              = $extreme_site[$i]->{set};
        my $extreme_middle_right     = $extreme_set->at( $half_length + 1 );
        my $extreme_middle_right_idx = $comparable_set->index($extreme_middle_right);
        my $right_wave_length = $next_extreme_middle_left_idx - $extreme_middle_right_idx + 1;
        $next_extreme_middle_left_idx = $extreme_middle_right_idx - 1;
        $extreme_site[$i]->{right_wave_length} = $right_wave_length;

        # amplitude
        my $extreme_gc      = $extreme_site[$i]->{gc};
        my $right_amplitude = abs( $extreme_gc - $next_extreme_gc );
        $extreme_site[$i]->{right_amplitude} = $right_amplitude;
        $prev_extreme_gc = $extreme_gc;
    }

    my $coll_extreme    = $db->get_collection('extreme');
    my $prev_extreme_id = undef;
    for my $ex (@extreme_site) {
        $ex->{chr_name}      = $align->{chr_name};
        $ex->{chr_start}     = $ex->{set}->min + $align->{chr_start} - 1;
        $ex->{chr_end}       = $ex->{set}->max + $align->{chr_start} - 1;
        $ex->{length}        = $ex->{set}->cardinality;
        $ex->{runlist}       = $ex->{chr_start} . '-' . $ex->{chr_end};
        $ex->{align_runlist} = $ex->{set}->runlist;
        delete $ex->{set};

        $ex->{align_id} = $align->{_id};

        $ex->{prev_id} = $prev_extreme_id;
        $prev_extreme_id = $coll_extreme->insert_one($ex)->inserted_id;
    }

    return;
}

sub insert_gsw {
    my $obj_gc = shift;
    my $db     = shift;
    my $align  = shift;
    my $seq    = shift;

    my $comparable_set = AlignDB::IntSpan->new( "1-" . $align->{length} );

    my $coll_extreme = $db->get_collection('extreme');
    my $coll_gsw     = $db->get_collection('gsw');

    # get gc sliding windows' sizes
    my $gsw_size = $obj_gc->gsw_size;

    my @extremes = $coll_extreme->find( { align_id => $align->{_id} } )->all;
    for my $ex (@extremes) {

        # bypass the first extreme
        next unless $ex->{prev_id};

        # get attrs of the two extremes
        my $ex_type             = $ex->{type};
        my $ex_gc               = $ex->{gc};
        my $ex_left_amplitude   = $ex->{left_amplitude};
        my $ex_left_wave_length = $ex->{left_wave_length};
        my $ex_set              = AlignDB::IntSpan->new( $ex->{align_runlist} );

        my $prev_ex     = $coll_extreme->find_one( { _id => $ex->{prev_id} } );
        my $prev_ex_gc  = $prev_ex->{gc};
        my $prev_ex_set = AlignDB::IntSpan->new( $prev_ex->{align_runlist} );

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
                type           => 'T',
                set            => $ex_set,
                distance       => 0,
                distance_crest => $gsw_density + 1,
                };
        }
        elsif ( $ex_type eq 'C' ) {    # crest
            $gsw_type = 'A';           # ascend, right of trough, left of crest
            push @gsw_windows,
                {
                type           => 'C',
                set            => $ex_set,
                distance       => $gsw_density + 1,
                distance_crest => 0,
                };
        }
        else {
            warn "extreme_type [$ex_type] error\n";
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
                    type           => $gsw_type,
                    set            => $gsw_set,
                    distance       => $gsw_distance,
                    distance_crest => $gsw_density - $gsw_distance + 1,
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
                $gsw->{chr_name}      = $align->{chr_name};
                $gsw->{chr_start}     = $gsw->{set}->min + $align->{chr_start} - 1;
                $gsw->{chr_end}       = $gsw->{set}->max + $align->{chr_start} - 1;
                $gsw->{length}        = $gsw->{set}->cardinality;
                $gsw->{runlist}       = $gsw->{chr_start} . '-' . $gsw->{chr_end};
                $gsw->{align_runlist} = $gsw->{set}->runlist;

                my $gsw_seq = substr $seq, $gsw->{set}->min - 1, $gsw->{length};
                $gsw->{gc} = calc_gc_ratio($gsw_seq);

                # ids
                $gsw->{align_id}        = $align->{_id};
                $gsw->{extreme_id}      = $ex->{_id};
                $gsw->{prev_extreme_id} = $prev_ex->{_id};

                # wave properties
                $gsw->{wave_length} = $gsw_wave_length;
                $gsw->{amplitude}   = $gsw_amplitude;
                $gsw->{trough_gc}   = $gsw_trough_gc;
                $gsw->{crest_gc}    = $gsw_crest_gc;
                $gsw->{gradient}    = $gsw_gradient;

                # pre allocate
                $gsw->{bed_count} = 0;

                # gsw cv
                my $resize_set = center_resize( $gsw->{set}, $comparable_set, $stat_segment_size );
                if ( !$resize_set ) {
                    print "    Can't resize window!\n";
                    $gsw->{gc_mean} = undef;
                    $gsw->{gc_cv}   = undef;
                    $gsw->{gc_std}  = undef;
                }
                else {
                    my ( $gc_mean, $gc_std, $gc_cv, $gc_mdcw )
                        = $obj_gc->segment_gc_stat( [$seq], $resize_set );
                    $gsw->{gc_mean} = $gc_mean;
                    $gsw->{gc_cv}   = $gc_cv;
                    $gsw->{gc_std}  = $gc_std;
                }

                delete $gsw->{set};
            }
            $coll_gsw->insert_many( \@gsw_windows );
        }
    }

    return;
}

__END__

