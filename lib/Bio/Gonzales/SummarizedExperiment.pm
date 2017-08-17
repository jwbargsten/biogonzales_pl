package Bio::Gonzales::SummarizedExperiment;

use warnings;
use strict;
use utf8;
use Carp;

use v5.11;
use IO::Handle ();

use List::MoreUtils qw/firstidx indexes/;
# Imports hier

use Moo;
use namespace::clean;

our $VERSION = 0.01_01;

has [qw/assay col_data row_data row_names col_names/] => ( is => 'rw', default => [] );

sub data   { shift->assay(@_) }
sub header { shift->col_names(@_) }

sub _idx {
  my ( $self, $m, $name ) = @_;
  firstidx { $_ eq $name } @{ $self->$m };
}

sub row_idx { shift->_idx( 'row_names', @_ ) }

sub col_idx { shift->_idx( 'col_names', @_ ) }

sub header_idx { shift->col_idx(@_) }

sub _idx_match {
  my ( $self, $m, $rex ) = @_;
  [ indexes { $_ =~ /$rex/ } @{ $self->$m } ];
}

sub row_idx_match { shift->_idx_match( 'row_names', @_ ) }

sub col_idx_match { shift->_idx_match( 'col_names', @_ ) }

sub header_idx_match {
  shift->col_idx_match(@_);
}

sub add_col {
  my ( $self, $col, $name ) = @_;

  my @names;
  push @names, $name if ( defined($name) );

  my @data = map { [$_] } @$col;

  return $self->cbind( \@data, \@names );
}

sub cbind {
  my ( $self, $data, $names ) = @_;

  push @{ $self->col_names }, @$names if ( $names && @$names );

  my $assay = $self->assay;

  if ( ref $data eq 'CODE' ) {
    for ( my $i = 0; $i < @$assay; $i++ ) {
      push @{ $assay->[$i] }, $data->( $self, $assay->[$i] );
    }
  } elsif ( ref $data eq 'ARRAY' ) {
    die "number of rows differ" unless ( @$data == @$assay );
    for ( my $i = 0; $i < @$assay; $i++ ) {
      push @{ $assay->[$i] }, @{ $data->[$i] };
    }
  } else {
    die "no code or array";
  }

  return $self;
}

sub add_cols {
  shift->cbind(@_);
}

sub rbind {
  my ( $self, $rows, $names ) = @_;

  push @{ $self->row_names }, @$names if ( $names && @$names );

  push @{ $self->assay }, @$rows;

  return $self;
}

sub add_rows { shift->rbind(@_) }

sub aggregate_by_idcs {
  my ( $self, $idcs, $code, $col_names ) = @_;

  my $assay = $self->assay;
  my %row_groups;
  my @key_names = @{ $self->col_names }[@$idcs];
  for ( my $i = 0; $i < @$assay; $i++ ) {
    my @key = @{ $assay->[$i] }[@$idcs];
    my $key = join( $;, @key );
    $row_groups{$key} //= { idcs => [], rows => [], key => \@key, key_names => \@key_names };
    push @{ $row_groups{$key}{idcs} }, $i;
    push @{ $row_groups{$key}{rows} }, $assay->[$i];
  }

  my @agg_assay;
  my @agg_row_names;
  for my $v ( values %row_groups ) {
    my ( $row, $row_name )
      = $code->( $self, $v->{key}, $v->{rows}, { names => $v->{key_names}, row_idcs => $v->{idcs} } );
    push @agg_assay, $row;
    push @agg_row_names, $row_name if ( defined($row_name) );
  }

  return __PACKAGE__->new(
    assay     => \@agg_assay,
    col_names => $col_names // [],
    row_names => \@agg_row_names
  );
}

sub aggregate_by_names {
  my ( $self, $names, $code, $col_names ) = @_;
  my @idcs;
  for my $n (@$names) {
    push @idcs, $self->col_idx($n);
  }

  return $self->aggregate_by_idcs( \@idcs, $code, $col_names );
}

1;
