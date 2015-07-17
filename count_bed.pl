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
use AlignDB::IntSpan;
use AlignDB::Window;
use AlignDB::Stopwatch;

use FindBin;
use lib "$FindBin::Bin/lib";
use MyUtil qw(check_coll);

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->new;
$Config = Config::Tiny->read("$FindBin::Bin/config.ini");

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

my @files;

my $run = "all";    # insert/count

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
    'r|run=s'    => \$run,
    'parallel=i' => \$parallel,
    'batch=i'    => \$batch_number,
    'nochr'      => \$nochr,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# Init objects
#----------------------------------------------------------#
$stopwatch->start_message("Count bed of $dbname...");

#----------------------------------------------------------#
# workers
#----------------------------------------------------------#
my $worker_insert = sub {
    my ( $self, $chunk_ref, $chunk_id ) = @_;
    my $file = $chunk_ref->[0];

    my $wid = MCE->wid;

    my $inner_watch = AlignDB::Stopwatch->new;
    $inner_watch->block_message("* Process task [$chunk_id] by worker #$wid");

    print "Reading file [$file]\n";

    # wait forever for responses
    my $mongo = MongoDB::MongoClient->new(
        host          => $server,
        port          => $port,
        query_timeout => -1,
    );
    my $db = $mongo->get_database($dbname);

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
            {   chr_name  => $chr,
                chr_start => { '$lte' => $start },
                chr_end   => { '$gte' => $end }
            }
        );
        if ( !$align ) {
            print "    Can't locate an align for $chr:$start-$end\n";
            next;
        }
        push @beds,
            {
            align_id  => $align->{_id},
            chr_name  => $chr,
            chr_start => $start,
            chr_end   => $end,
            };
    }
    close $data_fh;

    print "Inserting file [$file]\n";
    my $coll_bed = $db->get_collection('bed');
    while ( scalar @beds ) {
        my @batching = splice @beds, 0, 10000;
        $coll_bed->batch_insert( \@batching, { safe => 1 } );
    }
    print "Insert done.\n";
};

my $worker_count = sub {
    my $job    = shift;
    my @aligns = @$job;

    # wait forever for responses
    my $mongo = MongoDB::MongoClient->new(
        host          => $server,
        port          => $port,
        query_timeout => -1,
    );
    my $db = $mongo->get_database($dbname);

    my $coll_gsw   = $db->get_collection('gsw');
    my $coll_ofgsw = $db->get_collection('ofgsw');
    my $coll_bed   = $db->get_collection('bed');
    for my $align (@aligns) {
        my $chr_name  = $align->{chr_name};
        my $chr_start = $align->{chr_start};
        my $chr_end   = $align->{chr_end};
        printf "Process align %s:%s-%s\n", $chr_name, $chr_start, $chr_end;

        # all beds in this align
        my @beds = $coll_bed->find( { align_id => $align->{_id} } )->all;
        next unless @beds;
        printf "    %d beds in this align\n", scalar @beds;
        my $bed_chr_set = AlignDB::IntSpan->new;
        for (@beds) {
            $bed_chr_set->add_range( $_->{chr_start}, $_->{chr_end} );
        }

        # all gsws in this align
        my @gsws = $coll_gsw->find( { align_id => $align->{_id} } )->all;
        printf "    %d gsws in this align\n", scalar @gsws;
        my %gsw_count_of;
        for my $gsw (@gsws) {
            next if $bed_chr_set->intersect( $gsw->{runlist} )->is_empty;

            my $count = count_bed_in_sw( $coll_bed, $align->{_id}, $gsw );

            if ($count) {
                $gsw_count_of{ $gsw->{_id} } = $count;
            }
            else {
                printf "gsw %s matching wrong\n", $gsw->{runlist};
            }
        }

        # all gsws in this align
        my @ofgsws = $coll_ofgsw->find( { align_id => $align->{_id} } )->all;
        printf "    %d ofgsws in this align\n", scalar @ofgsws;
        my %ofgsw_count_of;
        for my $ofgsw (@ofgsws) {
            my $ofgsw_set = AlignDB::IntSpan->new(
                $ofgsw->{chr_start} . '-' . $ofgsw->{chr_end} );
            next if $bed_chr_set->intersect($ofgsw_set)->is_empty;

            my $count = count_bed_in_sw( $coll_bed, $align->{_id}, $ofgsw );

            if ($count) {
                $ofgsw_count_of{ $ofgsw->{_id} } = $count;
            }
            else {
                printf "ofgsw %s matching wrong\n", $ofgsw->{runlist};
            }
        }

        for my $key ( keys %gsw_count_of ) {
            $coll_gsw->update(
                { _id => MongoDB::OID->new( value => $key ) },
                { '$set' => { bed_count => $gsw_count_of{$key}, } },
                { safe   => 1 },
            );
        }
        for my $key ( keys %ofgsw_count_of ) {
            $coll_ofgsw->update(
                { _id => MongoDB::OID->new( value => $key ) },
                { '$set' => { bed_count => $ofgsw_count_of{$key}, } },
                { safe   => 1 },
            );
        }
    }
};

#----------------------------------------------------------#
# start update
#----------------------------------------------------------#
if ( $run eq "all" or $run eq "insert" ) {
    my $mongo = MongoDB::MongoClient->new(
        host          => $server,
        port          => $port,
        query_timeout => -1,
    );
    my $db       = $mongo->get_database($dbname);
    my $coll_bed = $db->get_collection('bed');

    $coll_bed->drop;
    $coll_bed->ensure_index( { 'align_id'  => 1 } );
    $coll_bed->ensure_index( { 'chr_name'  => 1 } );
    $coll_bed->ensure_index( { 'chr_start' => 1 } );
    $coll_bed->ensure_index( { 'chr_end'   => 1 } );

    my $mce = MCE->new( max_workers => $parallel, );
    $mce->foreach( \@files, $worker_insert, )
        ;    # foreach implies chunk_size => 1.

}

if ( $run eq "all" or $run eq "count" ) {
    my $mongo = MongoDB::MongoClient->new(
        host          => $server,
        port          => $port,
        query_timeout => -1,
    );
    my $db = $mongo->get_database($dbname);

    my $coll   = $db->get_collection('align');
    my @aligns = $coll->find->all;

    my $mce
        = MCE->new( max_workers => $parallel, chunk_size => $batch_number, );
    $mce->forchunk( \@aligns, $worker_count, );

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
        {   'align_id' => $align_id,
            '$or'      => [

                # bed    |----|
                # sw  |----|
                {   'chr_start' => {
                        '$gte' => $sw->{chr_start},
                        '$lte' => $sw->{chr_end},
                    }
                },

                # bed |----|
                # sw    |----|
                {   'chr_end' => {
                        '$gte' => $sw->{chr_start},
                        '$lte' => $sw->{chr_end},
                    }
                },

                # bed |--------|
                # sw    |----|
                {   'chr_start' => { '$lte' => $sw->{chr_start}, },
                    'chr_end'   => { '$gte' => $sw->{chr_end}, }
                },

                # bed   |----|
                # sw  |--------|
                {   'chr_start' => { '$gte' => $sw->{chr_start}, },
                    'chr_end'   => { '$lte' => $sw->{chr_end}, }
                },
            ]
        }
    )->count;

    return $count;
}

__END__

=head1 NAME

count_bed.pl - Add bed files to bed and count intersections.

=head1 SYNOPSIS

    perl count_bed.pl -d alignDB -f ~/Scripts/alignDB/ofg/spo11/spo11_hot.bed --batch 1 --parallel 1

=cut
