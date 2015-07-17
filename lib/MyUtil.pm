package MyUtil;
use strict;
use warnings;
use autodie;

use File::Spec;
use File::HomeDir;

use List::MoreUtils qw(firstidx uniq zip);
use List::Util qw(first max maxstr min minstr reduce shuffle sum);

use base 'Exporter';
use vars qw(@ISA @EXPORT_OK %EXPORT_TAGS);
@ISA = qw(Exporter);

%EXPORT_TAGS = (
    all => [
        qw{
            replace_home center_resize check_coll
            },
    ],
);

@EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

sub replace_home {
    my $path = shift;

    if ( $path =~ /^\~\// ) {
        $path =~ s/^\~\///;
        $path = File::Spec->catdir( File::HomeDir->my_home, $path );
    }

    return $path;
}

sub center_resize {
    my $old_set    = shift;
    my $parent_set = shift;
    my $resize     = shift;

    # find the middles of old_set
    my $half_size           = int( $old_set->size / 2 );
    my $midleft             = $old_set->at($half_size);
    my $midright            = $old_set->at( $half_size + 1 );
    my $midleft_parent_idx  = $parent_set->index($midleft);
    my $midright_parent_idx = $parent_set->index($midright);

    return unless $midleft_parent_idx and $midright_parent_idx;

    # map to parent
    my $parent_size  = $parent_set->size;
    my $half_resize  = int( $resize / 2 );
    my $new_left_idx = $midleft_parent_idx - $half_resize + 1;
    $new_left_idx = 1 if $new_left_idx < 1;
    my $new_right_idx = $midright_parent_idx + $half_resize - 1;
    $new_right_idx = $parent_size if $new_right_idx > $parent_size;

    my $new_set = $parent_set->slice( $new_left_idx, $new_right_idx );

    return $new_set;
}

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


1;
