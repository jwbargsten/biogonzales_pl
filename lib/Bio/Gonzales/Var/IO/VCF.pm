package EZ::VarNext::VcfReader;

use Mouse;

use warnings;
use strict;
use Carp;

use 5.010;

our $VERSION = 0.01_01;

with 'Bio::Gonzales::Util::Role::FileIO';

has meta       => ( is => 'rw', default => sub { {} } );
has sample_ids => ( is => 'rw', default => sub { [] } );

# stay consistent with GFF3 io
sub pragmas { shift->meta(@_) }

sub BUILD {
  my ($self) = @_;

  $self->_parse_header if ( $self->mode eq '<' );
}

sub _parse_header {
  my $self = shift;
  my $fhi  = $self->_fhi;

  my @sample_ids;
  my %meta;
  my $l;
  while ( defined( $l = $fhi->() ) ) {
    next if ( !$l || $l =~ /^\s*$/ );
    #looks like the header is over!
    last unless $l =~ /^\#/;
    if ( $l =~ /^\s*#CHROM/ ) {

      ( undef, undef, undef, undef, undef, undef, undef, undef, undef, @sample_ids ) = split /\t/, $l;
    } elsif ( $l =~ s/^##// ) {
      my ( $k, $v ) = split /=/, $l, 2;
      $meta{$k} = $v;
    } else {
      next;
    }
  }
  push @{ $self->_cached_records }, $l;

  $self->meta( \%meta );
  $self->sample_ids( \@sample_ids );

  return;
}

sub next_var {
  my ($self) = @_;

  my $fhi = $self->_fhi;

  my $l;
  while ( defined( $l = $fhi->() ) ) {
    if ( $l =~ /^\#/ || $l =~ /^\s*$/ ) {
      next;
    } else {
      last;
    }
  }
  return unless $l;

  my ( $chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @variants ) = split /\t/, $l;
  return {
    seq_id  => $chr,
    pos     => $pos,
    var_id  => $id,
    alleles => [ $ref, split(/,/, $alt) ],
    qual   => $qual,
    filter => $filter,
    info   => $info,
    format => $format,
    genotypes    => \@variants,
  };
}

__PACKAGE__->meta->make_immutable();
