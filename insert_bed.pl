#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use Roman;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use List::MoreUtils qw(any all uniq zip);

use MongoDB;
$MongoDB::BSON::looks_like_number = 1;
$MongoDB::BSON::utf8_flag_on      = 0;
use MongoDB::OID;

use MCE;

use AlignDB::GC;
use AlignDB::IntSpan;
use AlignDB::Window;
use AlignDB::Stopwatch;
use AlignDB::Util qw(:all);

use FindBin;
use lib "$FindBin::Bin/lib";
use MyUtil qw(center_resize);

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

# Database init values
my $server = $Config->{database}{server};
my $port   = $Config->{database}{port};
my $dbname = $Config->{database}{db};

my @tags;
my @types;
my $style = "center_intact";    # or center

my @files;

# AlignDB::GC options
my $wave_window_size  = $Config->{gc}{wave_window_size};
my $wave_window_step  = $Config->{gc}{wave_window_step};
my $vicinal_size      = $Config->{gc}{vicinal_size};
my $fall_range        = $Config->{gc}{fall_range};
my $gsw_size          = $Config->{gc}{gsw_size};
my $stat_segment_size = 500;
my $stat_window_size  = $Config->{gc}{stat_window_size};
my $stat_window_step  = $Config->{gc}{stat_window_step};

# run in parallel mode
my $parallel     = 1;
my $batch_number = 10;

my $nochr;

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'     => \$help,
    'man'        => \$man,
    's|server=s' => \$server,
    'P|port=i'   => \$port,
    'd|db=s'     => \$dbname,
    'f|file=s'   => \@files,
    'tag=s'      => \@tags,
    'type=s'     => \@types,
    'style=s'    => \$style,
    'parallel=i' => \$parallel,
    'batch=i'    => \$batch_number,
    'nochr'      => \$nochr,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# Init objects
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
    $inner_watch->block_message("* Process task [$chunk_id] by worker #$wid");

    my ( $file, $tag, $type ) = @$job;
    print "Reading file [$file]\n";

    # wait forever for responses
    my $mongo = MongoDB::MongoClient->new(
        host          => $server,
        port          => $port,
        query_timeout => -1,
    );
    my $db = $mongo->get_database($dbname);

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
            my $seq = $coll_seq->find_one( { _id => $align->{seq_id} } )->{seq};
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
        $coll_ofg->batch_insert( \@batching, { safe => 1 } );
    }
    print "Insert done.\n";
};

my $worker_sw = sub {
    my ( $self, $chunk_ref, $chunk_id ) = @_;
    my @aligns = @{$chunk_ref};

    my $wid = MCE->wid;

    my $inner_watch = AlignDB::Stopwatch->new;
    $inner_watch->block_message("* Process task [$chunk_id] by worker #$wid");

    # wait forever for responses
    my $mongo = MongoDB::MongoClient->new(
        host          => $server,
        port          => $port,
        query_timeout => -1,
    );
    my $db = $mongo->get_database($dbname);

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
        my $chr_name  = $align->{chr_name};
        my $chr_start = $align->{chr_start};
        my $chr_end   = $align->{chr_end};
        printf "Process align %s:%s-%s\n", $chr_name, $chr_start, $chr_end;

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
                    chr_name => $chr_name,
                    align_id => $align->{_id},
                    ofg_id   => $ofg->{_id},
                    type     => $rsw->{type},
                    distance => $rsw->{distance},
                };
                $ofgsw->{length}    = $rsw->{set}->size;
                $ofgsw->{chr_start} = $rsw->{set}->min + $chr_start - 1;
                $ofgsw->{chr_end}   = $rsw->{set}->max + $chr_start - 1;

                my $ofgsw_seq = substr $seq, $rsw->{set}->min - 1,
                    $ofgsw->{length};
                $ofgsw->{gc} = calc_gc_ratio($ofgsw_seq);

                $ofgsw->{bed_count} = 0;
                $ofgsw->{ofg}       = {
                    tag  => $ofg->{tag},
                    type => $ofg->{type},
                };

                # gsw cv
                my $resize_set
                    = center_resize( $rsw->{set}, $align_set,
                    $stat_segment_size );
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
            $coll_ofgsw->batch_insert( \@ofgsws, { safe => 1 } );
        }
    }
};

#----------------------------------------------------------#
# start update
#----------------------------------------------------------#
# insert bed files
{
    my $mongo = MongoDB::MongoClient->new(
        host          => $server,
        port          => $port,
        query_timeout => -1,
    );
    my $db       = $mongo->get_database($dbname);
    my $coll_ofg = $db->get_collection('ofg');
    $coll_ofg->drop;

    my @args = zip @files, @tags, @types;
    my @jobs;
    while ( scalar @args ) {
        my @batching = splice @args, 0, 3;
        push @jobs, [@batching];
    }

    my $mce = MCE->new( max_workers => $parallel, );
    $mce->foreach( \@jobs, $worker_insert, ); # foreach implies chunk_size => 1.

    $stopwatch->block_message("Add index to collection ofg");
    $coll_ofg->ensure_index( { 'align_id'  => 1 } );
    $coll_ofg->ensure_index( { 'chr_name'  => 1 } );
    $coll_ofg->ensure_index( { 'chr_start' => 1 } );
    $coll_ofg->ensure_index( { 'chr_end'   => 1 } );
    $coll_ofg->ensure_index( { 'tag'       => 1 } );
    $coll_ofg->ensure_index( { 'type'      => 1 } );
}

# ofgsw
{
    my $mongo = MongoDB::MongoClient->new(
        host          => $server,
        port          => $port,
        query_timeout => -1,
    );
    my $db = $mongo->get_database($dbname);

    my $coll = $db->get_collection('align');

    my @jobs = $coll->find->all;

    my $coll_ofgsw = $db->get_collection('ofgsw');
    $coll_ofgsw->drop;
    $coll_ofgsw->ensure_index( { 'align_id'  => 1 } );
    $coll_ofgsw->ensure_index( { 'ofg_id'    => 1 } );
    $coll_ofgsw->ensure_index( { 'chr_name'  => 1 } );
    $coll_ofgsw->ensure_index( { 'chr_start' => 1 } );
    $coll_ofgsw->ensure_index( { 'chr_end'   => 1 } );
    $coll_ofgsw->ensure_index( { 'type'      => 1 } );
    $coll_ofgsw->ensure_index( { 'distance'  => 1 } );
    $coll_ofgsw->ensure_index( { 'ofg.tag'   => 1 } );
    $coll_ofgsw->ensure_index( { 'ofg.type'  => 1 } );

    my $mce
        = MCE->new( max_workers => $parallel, chunk_size => $batch_number, );
    $mce->forchunk( \@jobs, $worker_sw, );

=head1 shell indexing
    ~/share/mongodb/bin/mongo Human_GC --eval "db.ofgsw.ensureIndex( { align_id: 1 }, {background: true} );" 
    ~/share/mongodb/bin/mongo Human_GC --eval "db.ofgsw.ensureIndex( { ofg_id: 1 }, {background: true} );" 
    ~/share/mongodb/bin/mongo Human_GC --eval "db.ofgsw.ensureIndex( { chr_name: 1 }, {background: true} );" 
    ~/share/mongodb/bin/mongo Human_GC --eval "db.ofgsw.ensureIndex( { chr_start: 1 }, {background: true} );" 
    ~/share/mongodb/bin/mongo Human_GC --eval "db.ofgsw.ensureIndex( { chr_end: 1 }, {background: true} );" 
    ~/share/mongodb/bin/mongo Human_GC --eval "db.ofgsw.ensureIndex( { type: 1 }, {background: true} );" 
    ~/share/mongodb/bin/mongo Human_GC --eval "db.ofgsw.ensureIndex( { distance: 1 }, {background: true} );" 
    ~/share/mongodb/bin/mongo Human_GC --eval "db.ofgsw.ensureIndex( { 'ofg.tag': 1 }, {background: true} );" 
    ~/share/mongodb/bin/mongo Human_GC --eval "db.ofgsw.ensureIndex( { 'ofg.type': 1 }, {background: true} );"
=cut    

}

$stopwatch->end_message;

exit;

__END__

=head1 NAME

insert_bed.pl - Add bed files to ofg and generate ofgsw.

=head1 SYNOPSIS

    mongo S288c --eval "db.dropDatabase();"
    perl gen_mg.pl -d S288c -n S288c --dir ~/data/alignment/yeast_combine/S288C  --parallel 1
    perl insert_bed.pl -d alignDB --tag hot --type hot -f ~/Scripts/alignDB/ofg/spo11/spo11_hot.bed --batch 10 --parallel 1

=cut