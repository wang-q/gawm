#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw(HelpMessage);
use Config::Tiny;
use FindBin;
use YAML qw(Dump Load DumpFile LoadFile);

use MongoDB;
$MongoDB::BSON::looks_like_number = 1;
use MongoDB::OID;

use AlignDB::IntSpan;
use AlignDB::Stopwatch;
use AlignDB::ToXLSX;

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

stat_mg.pl - Do stats 

=head1 SYNOPSIS

    perl stat_mg.pl [options]
      Options:
        --help      -?          brief help message
        --server        STR     MongoDB server IP/Domain name
        --port          INT     MongoDB server port
        --db        -d  STR     database name
        --output    -o  STR     output filename, default is [db.mg.xlsx]
        --by            STR     tag, type or tt, default is [tag]
        --replace       STR=STR replace strings in axis names
        --index                 add an index sheet

=cut

GetOptions(
    'help|?' => sub { HelpMessage(0) },
    'server=s' => \( my $server = $Config->{database}{server} ),
    'port=i'   => \( my $port   = $Config->{database}{port} ),
    'db|d=s'   => \( my $dbname = $Config->{database}{db} ),
    'output|o=s' => \my $outfile,
    'by=s'       => \( my $by = "tag" ),
    'replace=s'  => \( my %replace ),
    'index'      => \( my $add_index_sheet, ),
) or HelpMessage(1);

$outfile = "$dbname.mg.xlsx" unless $outfile;

#----------------------------------------------------------#
# Init section
#----------------------------------------------------------#
$stopwatch->start_message("Do stat for $dbname...");

my $write_obj = AlignDB::ToXLSX->new(
    outfile => $outfile,
    mocking => 1,
    replace => \%replace,
);

my $db = MongoDB::MongoClient->new(
    host          => $server,
    port          => $port,
    query_timeout => -1,
)->get_database($dbname);

#----------------------------------------------------------#
# worksheet -- distance_to_trough
#----------------------------------------------------------#
my $distance_to_trough = sub {
    my $sheet_name = 'distance_to_trough';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    my $coll = $db->get_collection('gsw');
    my $exists = $coll->find( { "gc.cv" => { '$exists' => 1 } } )->count;
    if ( !$exists ) {
        print "    gsw.gc.cv doesn't exist\n";
        print "    Skip sheet $sheet_name\n";
        return;
    }

    my @names = qw{_id AVG_gc AVG_cv AVG_bed COUNT};
    my $data = [ [], [], [], [], [] ];

    {    # write header
        ( $sheet_row, $sheet_col ) = ( 0, 0 );
        my %option = (
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
            header    => \@names,
        );
        ( $sheet, $sheet_row ) = $write_obj->write_header_direct( $sheet_name, \%option );
    }

    {    # write data
        my @docs = $coll->aggregate(
            [   { '$match' => { 'gce.distance' => { '$lte' => 15 } } },
                {   '$group' => {
                        $names[0] => '$gce.distance',
                        $names[1] => { '$avg' => '$gc.gc' },
                        $names[2] => { '$avg' => '$gc.cv' },
                        $names[3] => { '$avg' => '$bed_count' },
                        $names[4] => { '$sum' => 1 },
                    }
                },
                { '$sort' => { $names[0] => 1 } },
            ]
        )->all;
        for my $row (@docs) {
            for my $i ( 0 .. $#names ) {
                push @{ $data->[$i] }, $row->{ $names[$i] };
            }
        }
        $sheet->write( $sheet_row, 0, $data, $write_obj->format->{NORMAL} );
    }

    {    # chart
        my %option = (
            x_column    => 0,
            y_column    => 1,
            first_row   => 1,
            last_row    => 16,
            x_max_scale => 15,
            y_data      => $data->[1],
            x_title     => "Distance to GC troughs",
            y_title     => "Window GC",
            top         => 1,
            left        => 6,
        );
        $write_obj->draw_y( $sheet, \%option );

        $option{y_column} = 2;
        $option{y_title}  = "Window CV";
        $option{y_data}   = $data->[2];
        $option{top} += 18;
        $write_obj->draw_y( $sheet, \%option );

        $option{y_column} = 3;
        $option{y_title}  = "BED count";
        $option{y_data}   = $data->[3];
        $option{top} += 18;
        $write_obj->draw_y( $sheet, \%option );

        $option{y_column}  = 1;
        $option{y_title}   = "Window GC";
        $option{y_data}    = $data->[1];
        $option{y2_column} = 2;
        $option{y2_data}   = $data->[2];
        $option{y2_title}  = "Window CV";
        $option{top}       = 1;
        $option{left}      = 12;
        $write_obj->draw_2y( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- distance_to_crest
#----------------------------------------------------------#
my $distance_to_crest = sub {
    my $sheet_name = 'distance_to_crest';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    my $coll = $db->get_collection('gsw');
    my $exists = $coll->find( { 'gc.cv' => { '$exists' => 1 } } )->count;
    if ( !$exists ) {
        print "    gsw.gc.cv doesn't exist\n";
        print "    Skip sheet $sheet_name\n";
        return;
    }

    my @names = qw{_id AVG_gc AVG_cv AVG_bed COUNT};
    my $data = [ [], [], [], [], [] ];

    {    # write header
        ( $sheet_row, $sheet_col ) = ( 0, 0 );
        my %option = (
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
            header    => \@names,
        );
        ( $sheet, $sheet_row ) = $write_obj->write_header_direct( $sheet_name, \%option );
    }

    {    # write data
        my @docs = $coll->aggregate(
            [   { '$match' => { 'gce.distance_crest' => { '$lte' => 15 } } },
                {   '$group' => {
                        $names[0] => '$gce.distance_crest',
                        $names[1] => { '$avg' => '$gc.gc' },
                        $names[2] => { '$avg' => '$gc.cv' },
                        $names[3] => { '$avg' => '$bed_count' },
                        $names[4] => { '$sum' => 1 },
                    }
                },
                { '$sort' => { $names[0] => 1 } },
            ]
        )->all;
        for my $row (@docs) {
            for my $i ( 0 .. $#names ) {
                push @{ $data->[$i] }, $row->{ $names[$i] };
            }
        }
        $sheet->write( $sheet_row, 0, $data, $write_obj->format->{NORMAL} );
    }

    {    # chart
        my %option = (
            x_column    => 0,
            y_column    => 1,
            first_row   => 1,
            last_row    => 16,
            x_max_scale => 15,
            y_data      => $data->[1],
            x_title     => "Distance to GC crests",
            y_title     => "Window GC",
            top         => 1,
            left        => 6,
        );
        $write_obj->draw_y( $sheet, \%option );

        $option{y_column} = 2;
        $option{y_title}  = "Window CV";
        $option{y_data}   = $data->[2];
        $option{top} += 18;
        $write_obj->draw_y( $sheet, \%option );

        $option{y_column} = 3;
        $option{y_title}  = "BED count";
        $option{y_data}   = $data->[3];
        $option{top} += 18;
        $write_obj->draw_y( $sheet, \%option );

        $option{y_column}  = 1;
        $option{y_title}   = "Window GC";
        $option{y_data}    = $data->[1];
        $option{y2_column} = 2;
        $option{y2_data}   = $data->[2];
        $option{y2_title}  = "Window CV";
        $option{top}       = 1;
        $option{left}      = 12;
        $write_obj->draw_2y( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- gradient
#----------------------------------------------------------#
my $gradient = sub {
    my $sheet_name = 'gradient';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    my $coll = $db->get_collection('gsw');
    my $exists = $coll->find( { 'gc.cv' => { '$exists' => 1 } } )->count;
    if ( !$exists ) {
        print "    gsw.gc.cv doesn't exist\n";
        print "    Skip sheet $sheet_name\n";
        return;
    }

    my @names = qw{_id AVG_gc AVG_cv AVG_bed COUNT};
    my $data = [ [], [], [], [], [] ];

    {    # write header
        ( $sheet_row, $sheet_col ) = ( 0, 0 );
        my %option = (
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
            header    => \@names,
        );
        ( $sheet, $sheet_row ) = $write_obj->write_header_direct( $sheet_name, \%option );
    }

    {    # write data
        my @docs = $coll->aggregate(
            [   { '$match' => { 'gce.gradient' => { '$gte' => 1, '$lte' => 40 } } },
                {   '$group' => {
                        $names[0] => '$gce.gradient',
                        $names[1] => { '$avg' => '$gc.gc' },
                        $names[2] => { '$avg' => '$gc.cv' },
                        $names[3] => { '$avg' => '$bed_count' },
                        $names[4] => { '$sum' => 1 },
                    }
                },
                { '$sort' => { $names[0] => 1 } },
            ]
        )->all;
        for my $row (@docs) {
            for my $i ( 0 .. $#names ) {
                push @{ $data->[$i] }, $row->{ $names[$i] };
            }
        }
        $sheet->write( $sheet_row, 0, $data, $write_obj->format->{NORMAL} );
    }

    {    # chart
        my %option = (
            x_column    => 0,
            y_column    => 1,
            first_row   => 1,
            last_row    => 41,
            x_max_scale => 40,
            y_data      => $data->[1],
            x_title     => "Gradient",
            y_title     => "Window GC",
            top         => 1,
            left        => 6,
        );
        $write_obj->draw_y( $sheet, \%option );

        $option{y_column} = 2;
        $option{y_title}  = "Window CV";
        $option{y_data}   = $data->[2];
        $option{top} += 18;
        $write_obj->draw_y( $sheet, \%option );

        $option{y_column} = 3;
        $option{y_title}  = "BED count";
        $option{y_data}   = $data->[3];
        $option{top} += 18;
        $write_obj->draw_y( $sheet, \%option );

        $option{y_column}  = 1;
        $option{y_title}   = "Window GC";
        $option{y_data}    = $data->[1];
        $option{y2_column} = 2;
        $option{y2_data}   = $data->[2];
        $option{y2_title}  = "Window CV";
        $option{top}       = 1;
        $option{left}      = 12;
        $write_obj->draw_2y( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- ofg_all
#----------------------------------------------------------#
my $ofg_all = sub {
    my $sheet_name = 'ofg_all';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    my $coll = $db->get_collection('ofgsw');
    my $exists = $coll->find( { 'gc.cv' => { '$exists' => 1 } } )->count;
    if ( !$exists ) {
        print "    ofgsw.gc.cv doesn't exist\n";
        print "    Skip sheet $sheet_name\n";
        return;
    }

    my @names = qw{_id AVG_gc AVG_cv AVG_bed COUNT};
    my $data = [ [], [], [], [], [] ];

    {    # write header
        ( $sheet_row, $sheet_col ) = ( 0, 0 );
        my %option = (
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
            header    => \@names,
        );
        ( $sheet, $sheet_row ) = $write_obj->write_header_direct( $sheet_name, \%option );
    }

    {    # write data
        my @docs = $coll->aggregate(
            [   { '$match' => { 'ofg.distance' => { '$lte' => 15 } } },
                {   '$group' => {
                        $names[0] => '$ofg.distance',
                        $names[1] => { '$avg' => '$gc.gc' },
                        $names[2] => { '$avg' => '$gc.cv' },
                        $names[3] => { '$avg' => '$bed_count' },
                        $names[4] => { '$sum' => 1 },
                    }
                },
                { '$sort' => { $names[0] => 1 } },
            ]
        )->all;
        for my $row (@docs) {
            for my $i ( 0 .. $#names ) {
                push @{ $data->[$i] }, $row->{ $names[$i] };
            }
        }
        $sheet->write( $sheet_row, 0, $data, $write_obj->format->{NORMAL} );
    }

    {    # chart
        my %option = (
            x_column    => 0,
            y_column    => 1,
            first_row   => 1,
            last_row    => 16,
            x_max_scale => 15,
            y_data      => $data->[1],
            x_title     => "Distance to ofg",
            y_title     => "Window GC",
            top         => 1,
            left        => 6,
        );
        $write_obj->draw_y( $sheet, \%option );

        $option{y_column} = 2;
        $option{y_title}  = "Window CV";
        $option{y_data}   = $data->[2];
        $option{top} += 18;
        $write_obj->draw_y( $sheet, \%option );

        $option{y_column} = 3;
        $option{y_title}  = "BED count";
        $option{y_data}   = $data->[3];
        $option{top} += 18;
        $write_obj->draw_y( $sheet, \%option );

        $option{y_column}  = 1;
        $option{y_title}   = "Window GC";
        $option{y_data}    = $data->[1];
        $option{y2_column} = 2;
        $option{y2_data}   = $data->[2];
        $option{y2_title}  = "Window CV";
        $option{top}       = 1;
        $option{left}      = 12;
        $write_obj->draw_2y( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

my $ofg_tag_type = sub {

    my $coll = $db->get_collection('ofgsw');
    my $exists = $coll->find( { 'gc.cv' => { '$exists' => 1 } } )->count;
    if ( !$exists ) {
        print "    ofgsw.gc.cv doesn't exist\n";
        print "    Skip sheets ofg_tag_type\n";
        return;
    }

    my $write_sheet = sub {
        my ( $by, $bind ) = @_;

        my $sheet_name;
        if ( $by eq "tag" ) {
            $sheet_name = "ofg_tag_$bind";
        }
        elsif ( $by eq "type" ) {
            $sheet_name = "ofg_type_$bind";
        }
        elsif ( $by eq "tt" ) {
            $sheet_name = "ofg_tt_$bind";
        }

        # length limit of excel sheet names
        $sheet_name = substr $sheet_name, 0, 31;
        my $sheet;
        my ( $sheet_row, $sheet_col );

        my @names = qw{_id AVG_gc AVG_cv AVG_bed COUNT};
        my $data = [ [], [], [], [], [] ];

        {    # write header
            ( $sheet_row, $sheet_col ) = ( 0, 0 );
            my %option = (
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
                header    => \@names,
            );
            ( $sheet, $sheet_row ) = $write_obj->write_header_direct( $sheet_name, \%option );
        }

        {    # write data

            my $condition;
            if ( $by eq "tag" ) {
                $condition = { 'ofg.distance' => { '$lte' => 15 }, 'ofg.tag' => $bind, };
            }
            elsif ( $by eq "type" ) {
                $condition = { 'ofg.distance' => { '$lte' => 15 }, 'ofg.type' => $bind, };
            }
            elsif ( $by eq "tt" ) {
                my ( $tag, $type ) = split /\-/, $bind;
                $condition = {
                    'ofg.distance' => { '$lte' => 20 },
                    'ofg.tag'      => $tag,
                    'ofg.type'     => $type,
                };
            }

            my @docs = $coll->aggregate(
                [   { '$match' => $condition },
                    {   '$group' => {
                            $names[0] => '$ofg.distance',
                            $names[1] => { '$avg' => '$gc.gc' },
                            $names[2] => { '$avg' => '$gc.cv' },
                            $names[3] => { '$avg' => '$bed_count' },
                            $names[4] => { '$sum' => 1 },
                        }
                    },
                    { '$sort' => { $names[0] => 1 } },
                ]
            )->all;
            for my $row (@docs) {
                for my $i ( 0 .. $#names ) {
                    push @{ $data->[$i] }, $row->{ $names[$i] };
                }
            }
            $sheet->write( $sheet_row, 0, $data, $write_obj->format->{NORMAL} );
        }

        {    # chart
            my %option = (
                x_column    => 0,
                y_column    => 1,
                first_row   => 1,
                last_row    => 16,
                x_max_scale => 15,
                y_data      => $data->[1],
                x_title     => "Distance to ofg",
                y_title     => "Window GC",
                top         => 1,
                left        => 6,
            );
            $write_obj->draw_y( $sheet, \%option );

            $option{y_column} = 2;
            $option{y_title}  = "Window CV";
            $option{y_data}   = $data->[2];
            $option{top} += 18;
            $write_obj->draw_y( $sheet, \%option );

            $option{y_column} = 3;
            $option{y_title}  = "BED count";
            $option{y_data}   = $data->[3];
            $option{top} += 18;
            $write_obj->draw_y( $sheet, \%option );

            $option{y_column}  = 1;
            $option{y_title}   = "Window GC";
            $option{y_data}    = $data->[1];
            $option{y2_column} = 2;
            $option{y2_data}   = $data->[2];
            $option{y2_title}  = "Window CV";
            $option{top}       = 1;
            $option{left}      = 12;
            $write_obj->draw_2y( $sheet, \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    };

    my $ary_ref;
    if ( $by eq "tag" ) {
        $ary_ref = get_tags($db);
    }
    elsif ( $by eq "type" ) {
        $ary_ref = get_types($db);
    }
    elsif ( $by eq "tt" ) {
        $ary_ref = get_tts($db);
    }

    for ( @{$ary_ref} ) {
        $write_sheet->( $by, $_ );
    }
};

{
    &$distance_to_trough;
    &$distance_to_crest;
    &$gradient;
    &$ofg_all;
    &$ofg_tag_type;
}

if ($add_index_sheet) {
    $write_obj->add_index_sheet;
            print "Sheet [INDEX] has been generated.\n";

}

$stopwatch->end_message;
exit;

sub get_tags {
    my $db = shift;

    my $result = $db->run_command(
        [   "distinct" => "ofg",
            "key"      => "tag",
            "query"    => {},
        ]
    );
    my @values = sort @{ $result->{values} };

    return \@values;
}

sub get_types {
    my $db = shift;

    my $result = $db->run_command(
        [   "distinct" => "ofg",
            "key"      => "type",
            "query"    => {},
        ]
    );
    my @values = sort @{ $result->{values} };

    return \@values;
}

sub get_tts {
    my $db = shift;

    my $coll = $db->get_collection('ofg');

    my $result
        = $coll->aggregate( [ { '$group' => { "_id" => { type => '$type', tag => '$tag' } } } ] );

    my @values;
    for ( @{$result} ) {
        my $hash_ref = $_->{_id};
        push @values, $hash_ref->{tag} . '-' . $hash_ref->{type};
    }
    @values = sort @values;

    return \@values;
}

