#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw(HelpMessage);
use Config::Tiny;
use FindBin;
use YAML qw(Dump Load DumpFile LoadFile);

use AlignDB::Excel;
use AlignDB::Stopwatch;

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

chart_mg.pl - Use Win32::OLE to automate Excel chart

=head1 SYNOPSIS

    perl chart_mg.pl [options]
      Options:
        --help            brief help message
        --man             full documentation
        --infile          input file name (full path)
        --outfile         output file name

=cut

GetOptions(
    'help|?'      => sub { HelpMessage(0) },
    'infile|i=s'  => \( my $infile ),
    'outfile|o=s' => \( my $outfile ),
    'jc=s'              => \( my $jc_correction   = $Config->{stat}{jc_correction} ),
    'time_stamp=s'      => \( my $time_stamp      = $Config->{stat}{time_stamp} ),
    'add_index_sheet=s' => \( my $add_index_sheet = $Config->{stat}{add_index_sheet} ),
    'replace=s'         => \( my %replace ),
    'add_trend=s'       => \( my $add_trend ),
) or HelpMessage(1);

# die unless we got the mandatory argument
HelpMessage(1) unless defined $infile;

#----------------------------------------------------------#
# Init section
#----------------------------------------------------------#
$stopwatch->start_message("Processing $infile...");

my $excel_obj;
if ($outfile) {
    $excel_obj = AlignDB::Excel->new(
        infile  => $infile,
        outfile => $outfile,
        replace => \%replace,
    );
}
else {
    $excel_obj = AlignDB::Excel->new(
        infile  => $infile,
        replace => \%replace,
    );
    $outfile = $excel_obj->outfile;
}

$excel_obj->jc_correction if $jc_correction;

#----------------------------------------------------------#
# draw charts section
#----------------------------------------------------------#
{

    #----------------------------#
    # worksheet -- distance_to_trough
    #----------------------------#
    my $sheet_name = 'distance_to_trough';
    my %option     = (
        chart_serial => 1,
        x_column     => 1,
        y_column     => 2,
        first_row    => 2,
        last_row     => 17,
        x_max_scale  => 15,
        x_title      => "Distance to GC troughs",
        y_title      => "Window GC",
        Height       => 200,
        Width        => 260,
        Top          => 14.25,
        Left         => 550,
    );
    $option{x_scale_unit} = 5;
    $excel_obj->draw_y( $sheet_name, \%option );

    # chart 2
    $option{chart_serial}++;
    $option{y_column} = 3;
    $option{y_title}  = "Window CV";
    $option{Top} += $option{Height} + 14.25;
    $excel_obj->draw_y( $sheet_name, \%option );

    # chart 3
    $option{chart_serial}++;
    $option{y_column} = 4;
    $option{y_title}  = "BED count";
    $option{Top} += $option{Height} + 14.25;
    $excel_obj->draw_y( $sheet_name, \%option );

    #----------------------------#
    # worksheet -- distance_to_crest
    #----------------------------#
    $sheet_name           = 'distance_to_crest';
    $option{chart_serial} = 1;
    $option{y_column}     = 2;
    $option{x_title}      = "Distance to GC peaks";
    $option{y_title}      = "Window GC";
    $option{Top}          = 14.25;
    $option{x_scale_unit} = 5;
    $excel_obj->draw_y( $sheet_name, \%option );

    # chart 2
    $option{chart_serial}++;
    $option{y_column} = 3;
    $option{y_title}  = "Window CV";
    $option{Top} += $option{Height} + 14.25;
    $excel_obj->draw_y( $sheet_name, \%option );

    # chart 3
    $option{chart_serial}++;
    $option{y_column} = 4;
    $option{y_title}  = "BED count";
    $option{Top} += $option{Height} + 14.25;
    $excel_obj->draw_y( $sheet_name, \%option );

    #----------------------------#
    # worksheet -- distance_to_crest
    #----------------------------#
    $sheet_name           = 'gradient';
    $option{chart_serial} = 1;
    $option{y_column}     = 2;
    $option{last_row}     = 32;
    $option{x_max_scale}  = 30;
    $option{x_title}      = "Gradient";
    $option{y_title}      = "Window GC";
    $option{Top}          = 14.25;
    $option{x_scale_unit} = 5;
    $excel_obj->draw_y( $sheet_name, \%option );

    # chart 2
    $option{chart_serial}++;
    $option{y_column} = 3;
    $option{y_title}  = "Window CV";
    $option{Top} += $option{Height} + 14.25;
    $excel_obj->draw_y( $sheet_name, \%option );

    # chart 3
    $option{chart_serial}++;
    $option{y_column} = 4;
    $option{y_title}  = "BED count";
    $option{Top} += $option{Height} + 14.25;
    $excel_obj->draw_y( $sheet_name, \%option );
}

my @sheet_names = @{ $excel_obj->sheet_names };

{

    #----------------------------#
    # worksheet
    #----------------------------#
    my @sheets = grep {/^ofg/} @sheet_names;

    for (@sheets) {
        my $sheet_name = $_;

        my %option = (
            chart_serial => 1,
            x_column     => 1,
            first_row    => 2,
            last_row     => 17,
            x_max_scale  => 15,
            x_title      => "Distance to ofg",
            Height       => 200,
            Width        => 260,
            Top          => 14.25,
            Left         => 550,
        );
        $option{x_scale_unit} = 5;
        $option{y_column}     = 2;
        $option{y_title}      = "GC proportion";
        $option{y2_column}    = 3;
        $option{y2_title}     = "Window CV";
        $excel_obj->draw_2y( $sheet_name, \%option );
        delete $option{y2_column};
        delete $option{y2_title};

        # chart 2
        $option{chart_serial}++;
        $option{y_column} = 2;
        $option{y_title}  = "GC proportion";
        $option{Top} += $option{Height} + 14.25;
        $excel_obj->draw_y( $sheet_name, \%option );

        # chart 3
        $option{chart_serial}++;
        $option{y_column} = 3;
        $option{y_title}  = "Window CV";
        $option{Top} += $option{Height} + 14.25;
        $excel_obj->draw_y( $sheet_name, \%option );

        # chart 4
        $option{chart_serial}++;
        $option{y_column} = 4;
        $option{y_title}  = "BED count";
        $option{Top} += $option{Height} + 14.25;
        $excel_obj->draw_y( $sheet_name, \%option );
    }
}

#----------------------------------------------------------#
# POST Processing
#----------------------------------------------------------#
# add time stamp to "summary" sheet
$excel_obj->time_stamp("summary") if $time_stamp;

# add an index sheet
$excel_obj->add_index_sheet if $add_index_sheet;

print "$outfile has been generated.\n";

$stopwatch->end_message;
exit;

__END__
