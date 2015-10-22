#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw(HelpMessage);
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use MCE;
use File::Find::Rule;
use Path::Tiny;
use Text::CSV_XS;

use MongoDB;
$MongoDB::BSON::looks_like_number = 1;
use MongoDB::OID;

use AlignDB::IntSpan;
use AlignDB::Stopwatch;
use AlignDB::Util qw(:all);

use FindBin;

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
    'help|?' => sub { HelpMessage(0) },
    'server=s' => \( my $server           = $Config->{database}{server} ),
    'port=i'   => \( my $port             = $Config->{database}{port} ),
    'db|d=s'     => \( my $dbname           = $Config->{database}{db} ),
    'dir=s'      => \( my $dir ),
    'size|s=s'   => \( my $size_file ),
    'name|n=s'   => \( my $common_name ),
    'length=i'   => \( my $truncated_length = 500_000 ),
    'fill=i'     => \( my $fill             = 50 ),
    'min=i'      => \( my $min_length       = 5000 ),
    'parallel=i' => \( my $parallel         = 1 ),
) or HelpMessage(1);

# die unless we got the mandatory argument
HelpMessage(1) unless defined $dir;

#----------------------------------------------------------#
# Init
#----------------------------------------------------------#
$stopwatch->start_message("Generate database from fasta files");

{
    $stopwatch->block_message("Init database $dbname");
    my $db = MongoDB::MongoClient->new(
        host          => $server,
        port          => $port,
        query_timeout => -1,
    )->get_database($dbname);

    my $coll_chr = $db->get_collection('chromosome');
    $coll_chr->drop;
    $coll_chr->ensure_index( { 'chr_name' => 1 } );

    my $coll_seq = $db->get_collection('sequence');
    $coll_seq->drop;
    $coll_seq->ensure_index( { 'chr_name'  => 1 } );
    $coll_seq->ensure_index( { 'chr_start' => 1 } );
    $coll_seq->ensure_index( { 'chr_end'   => 1 } );
    $coll_seq->ensure_index( { 'level'     => 1 } );

    my $coll_align = $db->get_collection('align');
    $coll_align->drop;
    $coll_align->ensure_index( { 'chr_name'  => 1 } );
    $coll_align->ensure_index( { 'chr_start' => 1 } );
    $coll_align->ensure_index( { 'chr_end'   => 1 } );
    $coll_align->ensure_index( { 'seq_id'    => 1 } );
}

#----------------------------------------------------------#
# Search for size and fasta files
#----------------------------------------------------------#
if ( !defined $size_file ) {
    if ( path( $dir, 'chr.sizes' )->is_file ) {
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

my @files = sort File::Find::Rule->file->name( '*.fa', '*.fas', '*.fasta' )->in($dir);
printf "\n----Total .fa Files: %4s----\n\n", scalar @files;

{
    my $mongo = MongoDB::MongoClient->new(
        host          => $server,
        port          => $port,
        query_timeout => -1,
    );
    my $db = $mongo->get_database($dbname);

    my $coll_chr = $db->get_collection('chromosome');

    my $csv = Text::CSV_XS->new( { binary => 1, eol => "\n", sep => "\t" } );
    open my $csv_fh, "<", $size_file;
    my @fields = qw{chr_name chr_length};

    my @chrs;
    while ( my $row = $csv->getline($csv_fh) ) {
        my $data = {};
        for my $i ( 0 .. $#fields ) {
            $data->{ $fields[$i] } = $row->[$i];
        }
        $data->{common_name} = $common_name;

        push @chrs, $data;
    }
    close $csv_fh;

    $coll_chr->batch_insert( \@chrs, { safe => 1 }, );

    print "There are @{[$coll_chr->count]} documents in collection chromosome\n";
}

#----------------------------------------------------------#
# worker
#----------------------------------------------------------#
my $worker = sub {
    my ( $self, $chunk_ref, $chunk_id ) = @_;
    my $infile = $chunk_ref->[0];

    my $wid = MCE->wid;

    my $inner_watch = AlignDB::Stopwatch->new;
    $inner_watch->block_message("* Process task [$chunk_id] by worker #$wid\n* File [$infile]...");

    my $db = MongoDB::MongoClient->new(
        host          => $server,
        port          => $port,
        query_timeout => -1,
    )->get_database($dbname);

    my ( $seq_of, $seq_names ) = read_fasta($infile);
    my $chr_name   = $seq_names->[0];
    my $chr_seq    = $seq_of->{ $seq_names->[0] };
    my $chr_length = length $chr_seq;

    # find chromosome OID
    my $coll_chr = $db->get_collection('chromosome');
    my $chr_id = $coll_chr->find_one( { 'common_name' => $common_name, 'chr_name' => $chr_name } );
    return unless $chr_id;
    $chr_id = $chr_id->{_id};

    my $ambiguous_set = AlignDB::IntSpan->new;
    for ( my $pos = 0; $pos < $chr_length; $pos++ ) {
        my $base = substr $chr_seq, $pos, 1;
        if ( $base =~ /[^ACGT-]/i ) {
            $ambiguous_set->add( $pos + 1 );
        }
    }

    print "Ambiguous chromosome region for $chr_name:\n    " . $ambiguous_set->runlist . "\n";

    my $valid_set = AlignDB::IntSpan->new("1-$chr_length");
    $valid_set->subtract($ambiguous_set);
    $valid_set = $valid_set->fill( $fill - 1 );   

    print "Valid chromosome region for $chr_name:\n    " . $valid_set->runlist . "\n";

    my @regions;                                   # ([start, end], [start, end], ...)
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

    # collections
    my $coll_seq   = $db->get_collection('sequence');
    my $coll_align = $db->get_collection('align');

    for my $region (@regions) {
        my ( $start, $end ) = @{$region};
        my $length  = $end - $start + 1;
        my $runlist = AlignDB::IntSpan->new("$start-$end")->runlist;
        my $seq     = substr $chr_seq, $start - 1, $length;

        my $data_seq = {
            common_name => $common_name,
            chr_id      => $chr_id,
            chr_name    => $chr_name,
            chr_start   => $start,
            chr_end     => $end,
            chr_strand  => '+',
            seq         => $seq,
            length      => $length,
            runlist     => $runlist,
            level       => 1,              # top level
        };
        my $seq_id = $coll_seq->insert($data_seq);

        my $data_align = {
            chr_name   => $chr_name,
            chr_start  => $start,
            chr_end    => $end,
            chr_strand => '+',
            length     => $length,
            runlist    => $runlist,
            seq_id     => $seq_id,
        };
        $coll_align->insert($data_align);
    }

    $inner_watch->block_message( "$infile has been processed.", "duration" );

    return;
};

#----------------------------------------------------------#
# start
#----------------------------------------------------------#
my $mce = MCE->new( max_workers => $parallel, );
$mce->foreach( \@files, $worker, );    # foreach implies chunk_size => 1.

$stopwatch->end_message( "All files have been processed.", "duration" );

exit;

__END__
