package MyUtil;
use strict;
use warnings;
use autodie;

use Carp;

use AlignDB::IntSpan;

use base 'Exporter';
use vars qw(@ISA @EXPORT_OK %EXPORT_TAGS);
@ISA = qw(Exporter);

%EXPORT_TAGS = (
    all => [
        qw{
            check_coll process_message
            },
    ],
);

@EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

sub check_coll {
    my $db    = shift;
    my $name  = shift;
    my $field = shift;

    my $coll = $db->get_collection($name);

    my $total      = $coll->find->count;
    my $exists     = $coll->find( { $field => { '$exists' => 1 } } )->count;
    my $non_exists = $coll->find( { $field => { '$exists' => 0 } } )->count;

    return "For collection [$name], check field [$field]:\n"
        . "    Total $total\n    Exists $exists\n    Non exists $non_exists\n";
}

sub process_message {
    my $db       = shift;
    my $align_id = shift;

    my $coll_align = $db->get_collection('align');
    my $align = $coll_align->find_one( { _id => $align_id } );
    if ( !defined $align ) {
        printf "Can't find align for %s\n", $align_id;
        return;
    }

    printf "Process align %s(%s):%s-%s\n", $align->{chr}{name}, $align->{chr}{strand},
        $align->{chr}{start},
        $align->{chr}{end};

    return $align;
}

1;
