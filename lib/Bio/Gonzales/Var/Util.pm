package Bio::Gonzales::Var::Util;

use warnings;
use strict;
use Carp;

use 5.010;

use Exporter 'import';

our $VERSION = 0.01_01;

our %EXPORT_TAGS = ( 'all' => ['geno2haplo' ] );
our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

sub geno2haplo {
  my $genotypes = shift;
  # check if also coverage, etc. is part of the genotype then
  # split the genotypes into haplotypes
  my @haplotypes;
  my $phased = 1;
  for my $g_raw (@$genotypes) {
    my $g = index( $g_raw, ':' ) >= 0 ? substr( $g_raw, 0, index( $g_raw, ':' ) ) : $g_raw;
    # we need to find only one genotype of x/y to set phased to false
    $phased &&= not index( $g, '|' ) < 0;
    push @haplotypes, split /[|\/]/, $g;
  }
  return (\@haplotypes, $phased);
}


1;
