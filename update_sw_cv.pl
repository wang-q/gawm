#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use MongoDB;
$MongoDB::BSON::looks_like_number = 1;
$MongoDB::BSON::utf8_flag_on      = 0;
use MongoDB::OID;

use MCE;

use AlignDB::GC;
use AlignDB::Stopwatch;
use AlignDB::Util qw(:all);

use FindBin;
use lib "$FindBin::Bin/lib";
use MyUtil qw(center_resize check_coll);

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

my $stat_segment_size = 500;
my $stat_window_size  = 100;
my $stat_window_step  = 100;

# run in parallel mode
my $parallel = 1;

# number of alignments process in one child process
my $batch_number = 10;

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'     => \$help,
    'man'        => \$man,
    's|server=s' => \$server,
    'P|port=i'   => \$port,
    'd|db=s'     => \$dbname,
    'parallel=i' => \$parallel,
    'batch=i'    => \$batch_number,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
$stopwatch->start_message("Update sliding cv of $dbname...");

my @jobs;
{
    my $mongo = MongoDB::MongoClient->new(
        host          => $server,
        port          => $port,
        query_timeout => -1,
    );
    my $db = $mongo->get_database($dbname);

    my $coll = $db->get_collection('align');
    @jobs = $coll->find->all;
    printf "There are %d aligns totally.\n", scalar @jobs;

    # pre-allocate disk space
    print "Init disk space\n";
    $db->get_collection('ofgsw')->update(
        { gc_cv => { '$exists' => 0, } },
        {   '$set' => {
                gc_mean => undef,
                gc_cv   => undef,
                gc_std  => undef,
            }
        },
        { multiple => 1 },
    );

    $db->get_collection('gsw')->update(
        { gc_cv => { '$exists' => 0, } },
        {   '$set' => {
                gc_mean => undef,
                gc_cv   => undef,
                gc_std  => undef,
            }
        },
        { multiple => 1 },
    );
}

#----------------------------------------------------------#
# start update
#----------------------------------------------------------#
my $worker = sub {
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
    my $coll_ofgsw = $db->get_collection('ofgsw');
    my $coll_gsw   = $db->get_collection('gsw');

    # AlignDB::GC
    my $obj = AlignDB::GC->new(
        stat_window_size => $stat_window_size,
        stat_window_step => $stat_window_step,
        skip_mdcw        => 1,
    );

    for my $align (@aligns) {
        my $chr_name  = $align->{chr_name};
        my $chr_start = $align->{chr_start};
        my $chr_end   = $align->{chr_end};
        printf "Process align %s:%s-%s\n", $chr_name, $chr_start, $chr_end;

        my $seq = $coll_seq->find_one( { _id => $align->{seq_id} } )->{seq};
        my $align_set = AlignDB::IntSpan->new( "1-" . $align->{length} );

        #----------------------------#
        # ofgsw
        #----------------------------#
        my @ofgsws = $coll_ofgsw->find( { align_id => $align->{_id} } )->all;
        printf "    Updating %d ofgsws\n", scalar @ofgsws;
        my %stat_ofgsw_of;

        #----------------------------#
        # gsw
        #----------------------------#
        my @gsws = $coll_gsw->find( { align_id => $align->{_id} } )->all;
        printf "    Updating %d gsws\n", scalar @gsws;
        my %stat_gsw_of;

        #----------------------------#
        # calc
        #----------------------------#
        for my $ofgsw (@ofgsws) {
            my $window_set = AlignDB::IntSpan->new(
                $ofgsw->{align_start} . '-' . $ofgsw->{align_end} );
            my $resize_set
                = center_resize( $window_set, $align_set, $stat_segment_size );

            if ( !$resize_set ) {
                print "    Can't resize window!\n";
                next;
            }

            my ( $gc_mean, $gc_std, $gc_cv, $gc_mdcw )
                = $obj->segment_gc_stat( [$seq], $resize_set );

            $stat_ofgsw_of{ $ofgsw->{_id} } = {
                gc_mean => $gc_mean,
                gc_cv   => $gc_cv,
                gc_std  => $gc_std,
            };
        }
        for my $gsw (@gsws) {
            my $window_set = AlignDB::IntSpan->new( $gsw->{align_runlist} );
            my $resize_set
                = center_resize( $window_set, $align_set, $stat_segment_size );

            if ( !$resize_set ) {
                print "    Can't resize window!\n";
                next;
            }

            my ( $gc_mean, $gc_std, $gc_cv, $gc_mdcw )
                = $obj->segment_gc_stat( [$seq], $resize_set );

            $stat_gsw_of{ $gsw->{_id} } = {
                gc_mean => $gc_mean,
                gc_cv   => $gc_cv,
                gc_std  => $gc_std,
            };
        }

        #----------------------------#
        # update
        #----------------------------#
        # MongoDB::OID would be overloaded to string when as hash key
        for my $key ( keys %stat_ofgsw_of ) {
            $coll_ofgsw->update(
                { _id    => MongoDB::OID->new( value => $key ) },
                { '$set' => $stat_ofgsw_of{$key}, },
                { safe   => 1 },
            );
        }
        for my $key ( keys %stat_gsw_of ) {
            $coll_gsw->update(
                { _id    => MongoDB::OID->new( value => $key ) },
                { '$set' => $stat_gsw_of{$key}, },
                { safe   => 1 },
            );
        }
    }
};

#----------------------------------------------------------#
# start update
#----------------------------------------------------------#
my $mce = MCE->new( max_workers => $parallel, chunk_size => $batch_number, );
$mce->forchunk( \@jobs, $worker, );

#----------------------------#
# check
#----------------------------#
{
    my $mongo = MongoDB::MongoClient->new(
        host          => $server,
        port          => $port,
        query_timeout => -1,
    );
    my $db = $mongo->get_database($dbname);
    $stopwatch->block_message( check_coll( $db, 'ofgsw', 'gc_cv' ) );
    $stopwatch->block_message( check_coll( $db, 'gsw',   'gc_cv' ) );

}

$stopwatch->end_message;

exit;

__END__

=head1 NAME

update_sw_cv.pl - CV for ofgsw and gsw

=head1 SYNOPSIS

    perl update_sw_cv.pl -d alignDB --batch 10 --parallel 1

=cut
