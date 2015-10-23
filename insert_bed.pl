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
    my $coll_seq   = $db->get_collection('sequence');

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
            {   chr_name  => $chr,
                chr_start => { '$lte' => $start },
                chr_end   => { '$gte' => $end }
            }
        );
        if ( !$align ) {
            print "    Can't locate an align for $chr:$start-$end\n";
            next;
        }
        else {
            my $seq             = $coll_seq->find_one( { _id => $align->{seq_id} } )->{seq};
            my $length          = $end - $start + 1;
            my $ofg_align_start = $start - $align->{chr_start} + 1;
            my $ofg_align_end   = $end - $align->{chr_start} + 1;
            my $ofg_seq         = substr $seq, $ofg_align_start - 1, $length;
            my $ofg_gc          = calc_gc_ratio($ofg_seq);
            push @beds,
                {
                align_id    => $align->{_id},
                chr_name    => $chr,
                chr_start   => $start,
                chr_end     => $end,
                length      => $length,
                runlist     => AlignDB::IntSpan->new("$start-$end")->runlist,
                align_start => $ofg_align_start,
                align_end   => $ofg_align_end,
                gc          => $ofg_gc,
                tag         => $tag,
                type        => $type,
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

    my $coll_seq   = $db->get_collection('sequence');
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

    for my $align (@aligns) {
        printf "Process align %s:%s-%s\n", $align->{chr_name}, $align->{chr_start},
            $align->{chr_end};

        my @align_ofgs
            = $coll_ofg->find( { align_id => $align->{_id} } )->all;
        if ( @align_ofgs == 0 ) {
            warn "No ofgs in this align\n";
            next;
        }
        printf "    Find %d ofgs in this align\n", scalar @align_ofgs;

        my $seq = $coll_seq->find_one( { _id => $align->{seq_id} } )->{seq};
        my $align_set = AlignDB::IntSpan->new( "1-" . $align->{length} );

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
                @rsws = $window_maker->center_intact_window( $align_set,
                    $ofg->{align_start}, $ofg->{align_end} );
            }
            elsif ( $style eq 'center' ) {
                @rsws = $window_maker->center_window( $align_set,
                    $ofg->{align_start}, $ofg->{align_end} );
            }

            my @ofgsws;
            for my $rsw (@rsws) {
                my $ofgsw = {
                    chr_name => $align->{chr_name},
                    align_id => $align->{_id},
                    ofg_id   => $ofg->{_id},
                    type     => $rsw->{type},
                    distance => $rsw->{distance},
                };
                $ofgsw->{length}    = $rsw->{set}->size;
                $ofgsw->{chr_start} = $rsw->{set}->min + $align->{chr_start} - 1;
                $ofgsw->{chr_end}   = $rsw->{set}->max + $align->{chr_start} - 1;

                my $ofgsw_seq = substr $seq, $rsw->{set}->min - 1, $ofgsw->{length};
                $ofgsw->{gc} = calc_gc_ratio($ofgsw_seq);

                $ofgsw->{bed_count} = 0;
                $ofgsw->{ofg}       = {
                    tag  => $ofg->{tag},
                    type => $ofg->{type},
                };

                # gsw cv
                my $resize_set = center_resize( $rsw->{set}, $align_set, $stat_segment_size );
                if ( !$resize_set ) {
                    print "    Can't resize window!\n";
                    $ofgsw->{gc_mean} = undef;
                    $ofgsw->{gc_cv}   = undef;
                    $ofgsw->{gc_std}  = undef;
                }
                else {
                    my ( $gc_mean, $gc_std, $gc_cv, $gc_mdcw )
                        = $obj->segment_gc_stat( [$seq], $resize_set );
                    $ofgsw->{gc_mean} = $gc_mean;
                    $ofgsw->{gc_cv}   = $gc_cv;
                    $ofgsw->{gc_std}  = $gc_std;
                }

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
    $indexes->create_one( [ align_id => 1 ] );
    $indexes->create_one( [ chr_name => 1, chr_start => 1, chr_end => 1 ] );
    $indexes->create_one( [ tag      => 1 ] );
    $indexes->create_one( [ type     => 1 ] );

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
    my $coll_align = $db->get_collection('align');
    my @jobs       = $coll_align->find->all;

    # insert ofgsw
    my $name = "ofgsw";
    my $coll = $db->get_collection($name);
    $coll->drop;

    my $mce = MCE->new( max_workers => $parallel, chunk_size => $batch_number, );
    $mce->forchunk( \@jobs, $worker_sw, );

    # indexing
    $stopwatch->block_message("Indexing $name");
    my $indexes = $coll->indexes;
    $indexes->create_one( [ align_id => 1 ] );
    $indexes->create_one( [ ofg_id   => 1 ] );
    $indexes->create_one( [ chr_name => 1, chr_start => 1, chr_end => 1 ] );
    $indexes->create_one( [ type     => 1 ] );
    $indexes->create_one( [ distance => 1 ] );
    $indexes->create_one( { 'ofg.tag'  => 1 } );
    $indexes->create_one( { 'ofg.type' => 1 } );

    $stopwatch->block_message( check_coll( $db, $name, '_id' ) );
}

$stopwatch->end_message;

exit;

__END__
