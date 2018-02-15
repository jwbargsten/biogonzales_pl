package Bio::Gonzales::Var::Util;

use warnings;
use strict;
use Carp;

use 5.010;

use Exporter 'import';

use List::MoreUtils qw/uniq/;
use List::Util qw/max/;

our $VERSION = 0.01_01;

our %EXPORT_TAGS = ( 'all' => [qw/geno2haplo renumber_genotypes merge_alleles geno_get_gt_simple/ ] );
our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

sub geno2haplo {
  my $genotypes = shift;
  my $ploidy = shift;

  die "no ploidy defined for geno2haplo" unless(defined($ploidy));
  # check if also coverage, etc. is part of the genotype then
  # split the genotypes into haplotypes
  my $phased = 1;
  my @haplotypes;
  for my $g_raw (@$genotypes) {
    my $g = index( $g_raw, ':' ) >= 0 ? substr( $g_raw, 0, index( $g_raw, ':' ) ) : $g_raw;
    # we need to find only one genotype of x/y to set phased to false
    $phased &&= not index( $g, '|' ) < 0;
    my @h = split /[|\/]/, $g;
    @h = ('.')x$ploidy if(@h == 1 && $h[0] eq '.');
    die "ploidy mismatch in geno2haplo" if(@h != $ploidy);
    push @haplotypes, @h;
  }
  return ( \@haplotypes, $phased );
}

sub renumber_genotypes {
  my ( $map , $genotypes, ) = @_;
  my @renumbered;
  for my $g_raw (@$genotypes) {
    my $idx = index( $g_raw, ':' );
    my $g = $idx >= 0 ? substr( $g_raw, 0, $idx ) : $g_raw;
    my @g_split = split /([|\/])/, $g;
    for ( my $i = 0; $i < @g_split; $i += 2 ) {
      $g_split[$i] = $map->[ $g_split[$i] ] if($g_split[$i] ne '.');
    }
    if($idx < 0) {
      $g_raw = join '', @g_split;
    } else {
      substr( $g_raw, 0, $idx, join( '', @g_split));
    }
  }
  return $genotypes;
}

sub merge_alleles {
  my ( $ref_alleles, $alleles ) = @_;

  my $i = 0;
  my %ra = map { $_ => $i++ } @$ref_alleles;

  my @map;
  my @merged_alleles = @$ref_alleles;
  my $allele_idx     = @$ref_alleles;
  for ( my $idx = 0; $idx < @$alleles; $idx++ ) {
    if ( defined $ra{ $alleles->[$idx] } ) {
      $map[$idx] = $ra{ $alleles->[$idx] };
    } else {
      $map[$idx] = $allele_idx++;
      push @merged_alleles, $alleles->[$idx];
    }
  }
  return ( \@merged_alleles, \@map );
}

sub geno_get_gt_simple {
  my $genotypes = shift;
  my @res = map { index( $_, ':' ) >= 0 ? substr( $_, 0, index( $_, ':' ) ) : $_ } @$genotypes;
  return \@res;
}

sub geno_get_gt {
  my ( $gt, $smp_idcs, $gt_idx ) = @_;
  $gt_idx //= 0;

  my @gt;
  $smp_idcs = [ 0 .. (@$gt -1 )] unless($smp_idcs);
  for my $i (@$smp_idcs) {
    my $call = $gt->[$i];

    my @fields = split /:/, $call;
    if ( index( $fields[$gt_idx], '.' ) >= 0 ) {
      push @gt, '.';
    } else {
      push @gt, $fields[$gt_idx];
    }
  }
  return \@gt;
}

sub geno_count_gt {
  my $gt = shift;

  my %cnt;
  map { $cnt{$_}++ } @$gt;
  return \%cnt;
}

sub geno_nhet {
  my $gt = shift;

  my $nhet = 0;
  for my $g (@$gt) {
    my @called_alleles = split /[\/|]/, $g;
    $nhet++ if ( uniq(@called_alleles) > 1 );

  }
  return $nhet;
}

sub geno_is_poly {
  my $gt = shift;
  my @valid_gt = grep { $_ ne '.' } @$gt;
  return 0 unless (@valid_gt);
  return uniq(@valid_gt) == 1 ? 0 : 1;
}

sub geno_split_hom_het {
  my $gt = shift;
  my @hom;
  my @het;
  my $nmissing = 0;

  for my $g (@$gt) {
    if ( $g eq '.' ) {
      $nmissing++;
      next;
    }
    my @called_alleles = split /[\/|]/, $g;
    if ( uniq(@called_alleles) > 1 ) {
      push @het, $g;
    } else {
      push @hom, $g;
    }
  }
  return ( \@hom, \@het, $nmissing );
}

sub geno_nmissing {
  my $gt = shift;
  return scalar grep { $_ eq '.' } @$gt;
}

sub geno_major_alleles {
  my $cnt = shift;

  my @valid         = grep { $_ ne '.' } keys %$cnt;
  my $cnt_major     = max @{$cnt}{@valid};
  my @major_alleles = grep { $cnt->{$_} == $cnt_major && $_ ne '.' } @valid;
  return \@major_alleles;
}


1;
