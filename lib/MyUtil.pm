package MyUtil;
use strict;
use warnings;
use autodie;

use Carp;
use Path::Tiny;
use List::MoreUtils qw(firstidx uniq zip);
use List::Util qw(first max maxstr min minstr reduce shuffle sum);

use base 'Exporter';
use vars qw(@ISA @EXPORT_OK %EXPORT_TAGS);
@ISA = qw(Exporter);

%EXPORT_TAGS = (
    all => [
        qw{
            center_resize check_coll read_fasta calc_gc_ratio
            },
    ],
);

@EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

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

sub read_fasta {
    my $filename = shift;

    my @lines = path($filename)->lines;

    my @seq_names;
    my %seqs;
    for my $line (@lines) {
        if ( $line =~ /^\>[\w:-]+/ ) {
            $line =~ s/\>//;
            chomp $line;
            push @seq_names, $line;
            $seqs{$line} = '';
        }
        elsif ( $line =~ /^[\w-]+/ ) {
            $line =~ s/[^\w-]//g;
            chomp $line;
            my $seq_name = $seq_names[-1];
            $seqs{$seq_name} .= $line;
        }
        else {    # Blank line, do nothing
        }
    }

    return ( \%seqs, \@seq_names );
}

sub calc_gc_ratio {
    my @seqs = @_;

    for my $seq (@seqs) {
        _ref2str( \$seq );
    }

    my @ratios;
    for my $seq (@seqs) {

        # Count all four bases
        my $a_count = $seq =~ tr/Aa/Aa/;
        my $g_count = $seq =~ tr/Gg/Gg/;
        my $c_count = $seq =~ tr/Cc/Cc/;
        my $t_count = $seq =~ tr/Tt/Tt/;

        my $four_count = $a_count + $g_count + $c_count + $t_count;
        my $gc_count   = $g_count + $c_count;

        if ( $four_count == 0 ) {
            next;
        }
        else {
            my $gc_ratio = $gc_count / $four_count;
            push @ratios, $gc_ratio;
        }
    }

    return mean(@ratios);
}

# in situ convert reference of string to string
# For the sake of efficiency, the return value should be discarded
sub _ref2str {
    my $ref = shift;

    if ( ref $ref eq "REF" ) {
        $$ref = $$$ref;    # this is very weird, but it works
    }

    unless ( ref $ref eq "SCALAR" ) {
        carp "Wrong parameter passed\n";
    }

    return $ref;
}

sub mean {
    @_ = grep { defined $_ } @_;
    return unless @_;
    return $_[0] unless @_ > 1;
    return sum(@_) / scalar(@_);
}

1;
