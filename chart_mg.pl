#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw(HelpMessage);
use Config::Tiny;
use FindBin;
use YAML qw(Dump Load DumpFile LoadFile);

require AlignDB::Excel;
use AlignDB::Stopwatch;

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

chart_mg.pl - Use Win32::OLE to automate Excel chart

=head1 SYNOPSIS

    perl chart_mg.pl -i <file> [options]
      Options:
        --help      -?          brief help message
        --infile    -i  STR     input file name (full path)
        --outfile   -o  STR     output file name
        --time                  time stamp
        --index                 add an index sheet

=cut

GetOptions(
    'help|?'      => sub { HelpMessage(0) },
    'infile|i=s'  => \my $infile,
    'outfile|o=s' => \my $outfile,
    'jc=s'    => \( my $jc_correction   = $Config->{stat}{jc_correction} ),
    'time=s'  => \( my $time_stamp      = $Config->{stat}{time_stamp} ),
    'index=s' => \( my $add_index_sheet = $Config->{stat}{add_index_sheet} ),
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
    );
}
else {
    $excel_obj = AlignDB::Excel->new( infile => $infile, );
    $outfile = $excel_obj->outfile;
}

$excel_obj->jc_correction if $jc_correction;

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
