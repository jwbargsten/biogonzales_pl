package Bio::Gonzales::Stat::Util;
use warnings;
use strict;
use Carp;
use Statistics::Descriptive;
use POSIX qw/ceil/;

use 5.010;

use base 'Exporter';
our ( @EXPORT, @EXPORT_OK, %EXPORT_TAGS );
# VERSION

@EXPORT      = qw();
%EXPORT_TAGS = ();
@EXPORT_OK   = qw(hist_text);

sub hist_text {
  my $v    = shift;
  my $stat = Statistics::Descriptive::Full->new();
  $stat->add_data(@$v);
  my $iqr = $stat->quantile(3) - $stat->quantile(1);
  $iqr = 1 unless ($iqr);
  my $nclass = ceil( $stat->sample_range / ( ( $stat->count()**( -1 / 3 ) ) * 2 * $iqr ) );

  my $fd = $stat->frequency_distribution_ref( $nclass + 1 );

  my $max_len = 60;
  my $count   = $stat->count;

  my @keys = sort { $a <=> $b } keys %$fd;

  my @res;
  my $longest = -1;
  for my $i ( 0 .. $#keys ) {
    my $last = $i == 0 ? "-inf" : $keys[ $i - 1 ];
    my $curr = $keys[$i];
    my $kcnt = $fd->{$curr};
    my $len  = sprintf( "%.0f", ( $kcnt * $max_len / $count ) );
    push @res, [ "($last,$curr]", "#" x $len ];
    $longest = length( $res[-1][0] ) if ( length( $res[-1][0] ) > $longest );
  }

  my $res;
  for my $r (@res) {
    $res .= sprintf( "%${longest}s %s\n", @$r );
  }
  return $res;
}
