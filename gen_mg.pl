#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Config::Tiny;
use FindBin;
use YAML::Syck;

use File::Find::Rule;
use MCE;
use Path::Tiny;

use MongoDB;
$MongoDB::BSON::looks_like_number = 1;
use MongoDB::OID;

use AlignDB::IntSpan;
use AlignDB::Stopwatch;

use App::RL::Common;
use App::Fasops::Common;

use lib "$FindBin::RealBin/lib";
use MyUtil qw(check_coll);

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->read("$FindBin::RealBin/config.ini");

# record ARGV and Config
my $stopwatch = AlignDB::Stopwatch->new();

=head1 NAME

gen_mg.pl - Generate database from fasta files

=head1 SYNOPSIS

    perl gen_mg.pl --dir <directory> [options]
      Options:
        --help      -?          brief help message
        --server        STR     MongoDB server IP/Domain name
        --port          INT     MongoDB server port
        --db        -d  STR     database name
        --dir           STR     genome files' directory
        --size      -s  STR     chr.sizes
        --name      -n  STR     common name
        --length        INT     MongoDB has write lock, so we break genome into
                                pieces and make perl busy, [500_000]
        --fill          INT     fill gaps smaller than this value, [50]
        --min           INT     skip pieces smaller than this value, [5000]
        --parallel      INT     run in parallel mode, [1]

    cd ~/Scripts/gawm
    mongo S288c --eval "db.dropDatabase();"
    perl gen_mg.pl -d S288c --dir ~/data/alignment/yeast_combine/S288C --name S288c --parallel 1

    cd ~/Scripts/alignDB
    perl init/init_alignDB.pl -d S288Cvsself
    perl init/gen_alignDB_genome.pl -d S288Cvsself -t "4932,S288C" --dir /home/wangq/data/alignment/yeast65/S288C/  --parallel 4
    perl init/insert_gc.pl -d S288Cvsself --parallel 4

=cut

GetOptions(
    'help|?' => sub { Getopt::Long::HelpMessage(0) },
    'server=s'   => \( my $server           = $Config->{database}{server} ),
    'port=i'     => \( my $port             = $Config->{database}{port} ),
    'db|d=s'     => \( my $dbname           = $Config->{database}{db} ),
    'dir=s'      => \my $dir,
    'size|s=s'   => \my $size_file,
    'name|n=s'   => \my $common_name,
    'length=i'   => \( my $truncated_length = 500_000 ),
    'fill=i'     => \( my $fill             = 50 ),
    'min=i'      => \( my $min_length       = 5000 ),
    'parallel=i' => \( my $parallel         = 1 ),
) or Getopt::Long::HelpMessage(1);

# die unless we got the mandatory argument
Getopt::Long::HelpMessage(1) unless defined $dir;

#----------------------------------------------------------#
# Init
#----------------------------------------------------------#
$stopwatch->start_message("Generate database from fasta files");

if ( !defined $size_file ) {
    if ( path($dir)->is_file ) {
        if ( path($dir)->parent->child('chr.sizes')->is_file ) {
            $size_file = path($dir)->parent->child('chr.sizes')->absolute->stringify;
            printf "--size set to $size_file\n";
        }
    }
    elsif ( path( $dir, 'chr.sizes' )->is_file ) {
        $size_file = path( $dir, 'chr.sizes' )->absolute->stringify;
        printf "--size set to $size_file\n";
    }
    else {
        die "--size chr.sizes is needed\n";
    }
}
elsif ( !-e $size_file ) {
    die "--size chr.sizes doesn't exist\n";
}

if ( !defined $common_name ) {
    $common_name = path($dir)->basename;
}

{
    $stopwatch->block_message("Init database $dbname");
    my $client = MongoDB::MongoClient->new(
        host          => $server,
        port          => $port,
        query_timeout => -1,
    );
    $client->ns("$dbname.chr")->drop;
    $client->ns("$dbname.align")->drop;
}

#----------------------------------------------------------#
# Search for size and fasta files
#----------------------------------------------------------#
$stopwatch->block_message("Insert to chr");

my @files = sort File::Find::Rule->file->name( '*.fa', '*.fas', '*.fasta' )->in($dir);
printf "    Total .fa Files: [%s]\n", scalar @files;

{
    my $client = MongoDB::MongoClient->new(
        host          => $server,
        port          => $port,
        query_timeout => -1,
    );

    #@type MongoDB::Collection
    my $coll_chr = $client->ns("$dbname.chr");

    my @chrs;
    my $length_of = App::RL::Common::read_sizes($size_file);
    for my $key ( keys %{$length_of} ) {
        push @chrs,
            {
            name        => $key,
            length      => $length_of->{$key},
            common_name => $common_name,
            };
    }
    $coll_chr->insert_many( \@chrs );

    printf "    There are [%s] documents in collection chromosome\n", $coll_chr->count;
}

#----------------------------------------------------------#
# worker
#----------------------------------------------------------#
my $worker = sub {
    my ( $self, $chunk_ref, $chunk_id ) = @_;
    my $infile = $chunk_ref->[0];
    my $wid    = MCE->wid;

    my $inner    = AlignDB::Stopwatch->new;
    my $basename = path($infile)->basename;
    $inner->block_message("Process task [$chunk_id] by worker #$wid. [$basename]");

    #@type MongoDB::Database
    my $db = MongoDB::MongoClient->new(
        host          => $server,
        port          => $port,
        query_timeout => -1,
    )->get_database($dbname);

    my $seq_of = App::Fasops::Common::read_fasta($infile);

    for my $chr_name ( keys %{$seq_of} ) {
        print "    Processing chromosome $chr_name\n";

        my $chr_seq    = $seq_of->{$chr_name};
        my $chr_length = length $chr_seq;

        # find chromosome OID
        #@type MongoDB::Collection
        my $coll_chr = $db->get_collection('chr');
        my $chr_id = $coll_chr->find_one( { 'common_name' => $common_name, 'name' => $chr_name } );
        return unless $chr_id;
        $chr_id = $chr_id->{_id};

        my $ambiguous_set = AlignDB::IntSpan->new;
        for ( my $pos = 0; $pos < $chr_length; $pos++ ) {
            my $base = substr $chr_seq, $pos, 1;
            if ( $base =~ /[^ACGT-]/i ) {
                $ambiguous_set->add( $pos + 1 );
            }
        }
        printf "    Ambiguous region for %s:\n        %s\n", $chr_name, $ambiguous_set->runlist;

        my $valid_set = AlignDB::IntSpan->new("1-$chr_length");
        $valid_set->subtract($ambiguous_set);
        $valid_set = $valid_set->fill( $fill - 1 );
        printf "    Valid region for %s:\n        %s\n", $chr_name, $valid_set->runlist;

        my @regions;    # ([start, end], [start, end], ...)
        for my $set ( $valid_set->sets ) {
            my $size = $set->size;
            next if $size < $min_length;

            my @set_regions;
            my $pos = $set->min;
            my $max = $set->max;
            while ( $max - $pos + 1 > $truncated_length ) {
                push @set_regions, [ $pos, $pos + $truncated_length - 1 ];
                $pos += $truncated_length;
            }
            if ( scalar @set_regions > 0 ) {
                $set_regions[-1]->[1] = $max;
            }
            else {
                @set_regions = ( [ $pos, $max ] );
            }
            push @regions, @set_regions;
        }

        # insert to collection align
        #@type MongoDB::Collection
        my $coll_align = $db->get_collection('align');
        for my $region (@regions) {
            my ( $start, $end ) = @{$region};
            my $length  = $end - $start + 1;
            my $runlist = AlignDB::IntSpan->new("$start-$end")->runlist;
            my $seq     = substr $chr_seq, $start - 1, $length;

            my $data = {
                chr => {
                    common_name => $common_name,
                    _id         => $chr_id,
                    name        => $chr_name,
                    start       => $start,
                    end         => $end,
                    strand      => '+',
                    runlist     => $runlist,
                },
                length => $length,
                seq    => $seq,
            };
            $coll_align->insert_one($data);
        }
    }

    $inner->block_message( "$infile has been processed.", "duration" );

    return;
};

#----------------------------------------------------------#
# start
#----------------------------------------------------------#
my $mce = MCE->new( max_workers => $parallel, );
$mce->foreach( \@files, $worker, );    # foreach implies chunk_size => 1.

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
        my $name = "chr";
        $stopwatch->block_message("Indexing $name");
        my $coll    = $db->get_collection($name);
        my $indexes = $coll->indexes;
        $indexes->create_one( [ common_name => 1, name => 1 ], { unique => 1 } );
    }

    {
        my $name = "align";
        $stopwatch->block_message("Indexing $name");
        my $coll    = $db->get_collection($name);
        my $indexes = $coll->indexes;
        $indexes->create_one( [ "chr.name" => 1, "chr.start" => 1, "chr.end" => 1 ] );
    }

    $stopwatch->block_message( check_coll( $db, 'chr',   '_id' ) );
    $stopwatch->block_message( check_coll( $db, 'align', 'chr._id' ) );
}

$stopwatch->end_message( "All files have been processed.", "duration" );

exit;

__END__
