package Bio::Gonzales::Stat::Util;
use warnings;
use strict;
use Carp;
use Statistics::Descriptive;
use Number::Format qw/format_number/;
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
  my $c    = shift;
  my $stat = Statistics::Descriptive::Full->new();
  $stat->add_data(@$v);
  my $nclass = $c->{breaks} || nclass_sturges($stat);

  my $fd = $stat->frequency_distribution_ref( $nclass + 1 );

  my $max_len = 60;
  my $count   = $stat->count;

  my @keys = sort { $a <=> $b } keys %$fd;

  my @res;
  my $longest_interval = -1;
  my $longest_cnt      = -1;
  for my $i ( 0 .. $#keys ) {
    my $last = $i == 0 ? "-inf" : $keys[ $i - 1 ];
    my $curr = $keys[$i];
    my $kcnt = $fd->{$curr};
    next if ( $kcnt == 0 && $c->{skip_empty} );
    my $len = sprintf( "%.0f", ( $kcnt * $max_len / $count ) );
    push @res, [ sprintf( "%.2f", $curr ), "$kcnt", "#" x $len ];
    $longest_interval = length( $res[-1][0] ) if ( length( $res[-1][0] ) > $longest_interval );
    $longest_cnt      = length( $res[-1][1] ) if ( length( $res[-1][1] ) > $longest_cnt );
  }

  my $res;
  $res .= "median: " . format_number( $stat->median ) . "\n";
  $res .= "count: " . format_number( $stat->count ) . "\n";
  $res .= "\n";
  for my $r (@res) {
    $res .= sprintf( "<= %${longest_interval}s (%${longest_cnt}s)  %s\n", @$r );
  }
  return $res;
}

sub nclass_fd {
  my ($stat) = @_;
  my $iqr = $stat->quantile(3) - $stat->quantile(1);
  $iqr = 1 unless ($iqr);
  return ceil( $stat->sample_range / ( ( $stat->count()**( -1 / 3 ) ) * 2 * $iqr ) );
}

sub nclass_sturges {
  my ($stat) = @_;
  return ceil(log($stat->count)/log(2) + 1);
}
