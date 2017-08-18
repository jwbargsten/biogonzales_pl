package Bio::Gonzales::SummarizedExperiment;

use warnings;
use strict;
use utf8;
use Carp;

use v5.11;
use IO::Handle ();

use List::MoreUtils qw/firstidx indexes uniq any/;
use Bio::Gonzales::Matrix::IO qw/mspew/;
use Algorithm::Loops qw/MapCarU/;

use Data::Dumper;

use Clone;

# Imports hier

use Moo;
use namespace::clean;

our $VERSION  = 0.01_01;
our $NA_VALUE = 'NA';

has [qw/assay col_data row_data row_names col_names/] => ( is => 'rw', default => sub { [] } );

sub data   { shift->assay(@_) }
sub header { shift->col_names(@_) }

sub _idx {
  my ( $self, $m, $name ) = @_;
  unless ($name) {
    carp "$name not found in column names";
    return -1;
  }
  firstidx { $_ eq $name } @{ $self->$m };
}

sub new_from_mslurp {
  my ( $class, $m, $cn, $rn ) = @_;
  $cn //= [];
  $rn //= [];
  return __PACKAGE__->new( assay => $m, col_names => $cn, row_names => $rn );
}

sub spew_assay {
  my $self  = shift;
  my $src   = shift;
  my $param = shift // {};
  my %c     = %$param;
  $c{header} = $self->col_names if ( $param->{header} || $param->{col_names} );
  $c{row_names} = $self->row_names if ( $param->{row_names} );
  mspew( $src, $self->assay, \%c );
  return $self;
}

sub sort {
  my ( $self, $code ) = @_;

  my $assay     = $self->assay;
  my $row_names = $self->row_names;
  my $row_data  = $self->row_data;

  my $nrow = $self->nrow;
  my @idcs = 0 .. ( $nrow - 1 );
  @idcs = sort { $code->( $assay->[$a], $assay->[$b] ) } @idcs;
  $self->assay(     [ @{$assay}[@idcs] ] );
  $self->row_names( [ @{$row_names}[@idcs] ] ) if (@$row_names);
  $self->row_data(  [ @{$row_data}[@idcs] ] ) if (@$row_data);

  return $self;
}

sub row_idx { shift->_idx( 'row_names', @_ ) }

sub col_idx { shift->_idx( 'col_names', @_ ) }

sub header_idx { shift->col_idx(@_) }

sub transpose {
  my $self = shift;

 my @assay_t = MapCarU { [@_] } @{$self->{assay}};
 my @row_data_t = MapCarU { [@_] } @{$self->{row_data}};
 my @col_data_t = MapCarU { [@_] } @{$self->{col_data}};


  return __PACKAGE__->new(
    assay     => \@assay_t,
    row_names => Clone::clone( $self->col_names ),
    row_data  => \@row_data_t,
    col_names => Clone::clone( $self->row_names ),
    col_data  => \@col_data_t,
  );
}

sub make_consistent {

}

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
  my ( $self, $assay_col, $name, $col_data ) = @_;

  my @names;
  push @names, $name if ( defined($name) );

  $col_data = [ map { [$_] } @$col_data ] if ( $col_data && @$col_data );

  my $data;
  if ( ref $assay_col eq 'CODE' ) {
    $data = $assay_col;
  } else {
    $data = [ map { [$_] } @$assay_col ];
  }

  return $self->cbind( $data, \@names, $col_data );
}

sub cbind {
  my ( $self, $data, $names, $col_data_n ) = @_;

  $col_data_n //= [];
  push @{ $self->col_names }, @$names if ( $names && @$names );

  my $assay = $self->assay;
  my $ncol  = $self->ncol;

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

  my $col_data = $self->col_data;
  # the assay is already updated, so ncol will represent the new number
  my $ncol_added = $self->ncol - $ncol;
  die "col data dims differ" if ( @$col_data && @$col_data_n && @$col_data != @$col_data_n );

  my $col_data_ncol = @$col_data > @$col_data_n ? @$col_data : @$col_data_n;
  for ( my $i = 0; $i < $col_data_ncol; $i++ ) {
    $col_data->[$i] //= [ ($NA_VALUE) x $ncol ];
    push @{ $col_data->[$i] }, @{ $col_data_n->[$i] // [ ($NA_VALUE) x $ncol_added ] };
  }

  return $self;
}

sub add_cols {
  shift->cbind(@_);
}

sub group {
  my ( $self, $idcs ) = @_;

  my $assay     = $self->assay;
  my $row_names = $self->row_names;
  my $row_data  = $self->row_data;
  my %groups;
  my @key_names = @{ $self->col_names }[@$idcs];
  for ( my $i = 0; $i < @$assay; $i++ ) {
    my @key = @{ $assay->[$i] }[@$idcs];
    my $key = join( $;, @key );
    $groups{$key}
      //= { idcs => [], rows => [], key => \@key, key_names => \@key_names, row_names => [], row_data => [] };
    push @{ $groups{$key}{idcs} },      $i;
    push @{ $groups{$key}{rows} },      $assay->[$i];
    push @{ $groups{$key}{row_names} }, $row_names->[$i] if ( $row_names && @$row_names );
    push @{ $groups{$key}{row_data} },  $row_data->[$i] if ( $row_data && @$row_data );
  }
  return \%groups;
}

sub ncol {
  my $assay = shift->assay;
  return unless ( $assay && @$assay );

  return scalar @{ $assay->[0] };
}

sub nrow { scalar @{ shift->assay }; }

sub subset {
  my ( $self, $code ) = @_;

  my $assay     = $self->assay;
  my $row_names = $self->row_names;
  my $row_data  = $self->row_data;

  my @row_names_new;
  my @assay_new;
  my @row_data_new;
  for ( my $i = 0; $i < @$assay; $i++ ) {
    if ( $code->( $self, $assay->[$i], $row_names->[$i], $row_data->[$i] ) ) {
      push @assay_new,     Clone::clone( $assay->[$i] );
      push @row_names_new, Clone::clone( $row_names->[$i] ) if ( $row_names && @$row_names );
      push @row_data_new,  Clone::clone( $row_data->[$i] ) if ( $row_data && @$row_data );
    }
  }
  return __PACKAGE__->new(
    assay     => \@assay_new,
    row_names => \@row_names_new,
    row_data  => \@row_data_new,
    col_names => Clone::clone( $self->col_names ),
    col_data  => Clone::clone( $self->col_data )
  );
}

sub _invert_idcs {
  my ( $n, $idcs ) = @_;

  my @inv;
  my %m = map { $_ => 1 } @$idcs;
  for ( my $i = 0; $i < $n; $i++ ) {
    push @inv, $i unless ( $m{$i} );
  }
  return \@inv;
}

sub merge {
  my ( $se_x, $se_y, $param ) = @_;

  my %param = ( join => 'inner', %{ $param // {} } );
  my $by_x = $param{by_x} // $param{by};
  my $by_y = $param{by_y} // $param{by};
  die "join needs same amount of rows on both sets" unless ( @$by_x == @$by_y );
  my $idcs_x = $se_x->col_names_to_idcs($by_x);
  my $idcs_y = $se_y->col_names_to_idcs($by_y);

  my $col_data_x = $se_x->col_data;
  my $col_data_y = $se_y->col_data;

  die "col data has different number of rows"
    if ( @$col_data_x && @$col_data_y && @$col_data_x != @$col_data_y );

  my $assay_x = $se_x->assay;
  my $assay_y = $se_y->assay;

  my $groups_x = $se_x->group($idcs_x);
  my $groups_y = $se_y->group($idcs_y);

  my $ncol_common = @$by_x;
  my $ncol_only_x = $se_x->ncol - $ncol_common;
  my $ncol_only_y = $se_y->ncol - $ncol_common;
  my $ncol_total  = $ncol_only_x + $ncol_only_y + $ncol_common;

  my @row_names_new;
  my @assay_new;
  my @row_data_new;
  my @keys = uniq( keys(%$groups_x), keys(%$groups_y) );
  my $inv_by_y = _invert_idcs( $se_y->ncol, $idcs_y );
  for my $k (@keys) {
    my $data_x = $groups_x->{$k};
    my $data_y = $groups_y->{$k};

    next if ( $param{join} eq 'inner' && !( $data_x && $data_y ) );
    if ( $param{join} eq 'left' || $param{join} eq 'full' ) {
      next if ( $param{join} eq 'left' && !$data_x );
      $data_y //= {
        idcs      => [-1],
        rows      => [ [ ($NA_VALUE) x $se_y->ncol ] ],
        key       => $data_x->{key},
        key_names => $data_x->{key_names},
        row_names => [],
        row_data  => []
      };
    }
    if ( $param{join} eq 'right' || $param{join} eq 'full' ) {
      next if ( $param{join} eq 'right' && !$data_y );
      my @row = ( ($NA_VALUE) x $se_x->ncol );
      @row[@$idcs_x] = @{ $data_y->{key} };
      $data_x //= {
        idcs      => [-1],
        rows      => [ \@row ],
        key       => $data_y->{key},
        key_names => $data_y->{key_names},
        row_names => [],
        row_data  => []
      };
    }

    for ( my $i = 0; $i < @{ $data_x->{rows} }; $i++ ) {
      for ( my $j = 0; $j < @{ $data_y->{rows} }; $j++ ) {

        my @row = ( @{ $data_x->{rows}[$i] }, @{ $data_y->{rows}[$j] }[@$inv_by_y] );
        push @assay_new,     \@row;
        push @row_names_new, $data_x->{row_names}[$i]
          if ( $data_x->{row_names} && @{ $data_x->{row_names} } );
        push @row_data_new, Clone::clone( $data_x->{row_data}[$i] )
          if ( $data_x->{row_data} && @{ $data_x->{row_data} } );
      }
    }
  }

  my @col_names = ( @{ $se_x->col_names }, @{ $se_y->col_names }[@$inv_by_y] );
  my @col_data;

  if ( @$col_data_x || @$col_data_y ) {
    for ( my $i = 0; $i < $ncol_total; $i++ ) {
      my $cd_x = $col_data_x->[$i] // [ ($NA_VALUE) x $ncol_only_x ];
      my $cd_y = $col_data_y->[$i] // [ ($NA_VALUE) x $ncol_only_y ];
      push @col_data, [ @$cd_x, @{$cd_y}[@$inv_by_y] ];
    }
  }
  return __PACKAGE__->new(
    assay     => \@assay_new,
    col_data  => \@col_data,
    col_names => \@col_names,
    row_names => \@row_names_new,
    row_data  => \@row_data_new
  );

}

sub clone {
  my $self = shift;
  return __PACKAGE__->new(
    assay     => Clone::clone( $self->assay ),
    col_data  => Clone::clone( $self->col_data ),
    col_names => Clone::clone( $self->col_names ),
    row_names => Clone::clone( $self->row_names ),
    row_data  => Clone::clone( $self->row_data ),
  );
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

  my $row_groups = $self->group($idcs);

  my @agg_assay;
  my @agg_row_names;
  for my $v ( values %$row_groups ) {
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

sub col_names_to_idcs {
  my ( $self, $names ) = @_;
  my @idcs = map { $self->col_idx($_) } @$names;
  die "could not find all idcs " . jon( ", ", @$names ) if ( any { $_ < 0 } @idcs );
  return \@idcs;
}

sub aggregate_by_names {
  my ( $self, $names, $code, $col_names ) = @_;
  my $idcs = $self->col_names_to_idcs($names);

  return $self->aggregate_by_idcs( $idcs, $code, $col_names );
}

sub col_idx_map {
  my $i = 0;
  return { map { $_ => $i++ } @{ shift->col_names } };
}

sub col_rename {
  my ( $self, $old, $new ) = @_;

  my $idx = $self->col_idx($old);
  die if ( $idx < 0 );
  $self->col_names->[$idx] = $new;
  return $self;
}

sub row_apply {
  my ( $self, $code ) = @_;

  my @res;
  my $assay = $self->assay;
  for ( my $i = 0; $i < @$assay; $i++ ) {
    push @res, $code->( $self, $assay->[$i] );
  }
  return @res;
}

sub slice_by_idcs {
  my ( $self, $idcs ) = @_;
  #FIXME implement cloning

  my @assay_new = map { [ @{$_}[@$idcs] ] } @{ $self->assay };
  my @new_colnames;
  @new_colnames = @{ $self->col_names }[@$idcs] if ( $self->has_col_names );

  my @new_coldata;
  @new_coldata = map { [ @{$_}[@$idcs] ] } @{ $self->col_data } if ( $self->has_col_data );

  return __PACKAGE__->new(
    assay     => \@assay_new,
    row_names => Clone::clone( $self->row_names ),
    row_data  => Clone::clone( $self->row_data ),
    col_names => \@new_colnames,
    col_data  => \@new_coldata,
  );
}

sub has_col_names {
  my $c = shift->col_names;
  return $c && @$c;
}

sub has_col_data {
  my $c = shift->col_data;
  return $c && @$c;
}

sub has_row_data {
  my $c = shift->row_data;
  return $c && @$c;
}

sub has_row_names {
  my $c = shift->row_names;
  return $c && @$c;
}

sub extract_col_by_idx {
  my ( $self, $idx ) = @_;

  my $assay = $self->assay;
  return [ map { $_->[$idx] } @$assay ];
}

sub extract_col_by_name {
  my $self = shift;
  return $self->extract_col_by_idx( $self->col_idx(@_) );
}

sub slice_by_names {
  my ( $self, $names ) = @_;
  my $idcs = $self->col_names_to_idcs($names);
  return $self->slice_by_idcs($idcs);
}

1;
