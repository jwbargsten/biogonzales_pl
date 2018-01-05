package Bio::Gonzales::SummarizedExperiment;

use warnings;
use strict;
use utf8;
use Carp;

use v5.11;
use IO::Handle ();

use List::MoreUtils qw/firstidx indexes any/;
use List::Util qw/max/;
use Bio::Gonzales::Matrix::IO qw/mspew mslurp/;
use Algorithm::Loops qw/MapCarU/;

use Bio::Gonzales::Matrix::Util;

use Data::Dumper;
use JSON::XS;

use Clone;

# Imports hier

use Moo;
use namespace::clean;

our $VERSION  = 0.01_01;
our $NA_VALUE = 'NA';

has [qw/assay col_data row_data row_names col_names row_data_names col_data_names meta_data/] =>
  ( is => 'rw', default => sub { [] } );

has na_value => ( is => 'rw', default => $NA_VALUE );

sub data { shift->assay(@_) }

sub header { shift->col_names(@_) }

sub slurp_assay {
  my $class = shift;

  my ( $m, $cn, $rn ) = mslurp(@_);
  $cn //= [];
  $rn //= [];
  return $class->new( assay => $m, col_names => $cn, row_names => $rn );
}

sub spew_assay {
  my $self  = shift;
  my $src   = shift;
  my $param = shift // {};
  my %c     = %$param;
  $c{header}    = $self->col_names if ( $param->{header} || $param->{col_names} );
  $c{row_names} = $self->row_names if ( $param->{row_names} );
  $c{col_data}  = $self->col_data  if ( $param->{col_data} );
  $c{na_value} = $self->na_value;
  mspew( $src, $self->assay, \%c );
  return $self;
}

sub _reorder {
  my $self = shift;
  my $idcs = shift;

  my $assay          = $self->assay;
  my $row_names      = $self->row_names;
  my $row_data       = $self->row_data;
  my $row_data_names = $self->row_data_names;

  $self->assay(          [ @{$assay}[@$idcs] ] );
  $self->row_names(      [ @{$row_names}[@$idcs] ] ) if (@$row_names);
  $self->row_data(       [ @{$row_data}[@$idcs] ] ) if (@$row_data);
  $self->row_data_names( [ @{$row_data_names}[@$idcs] ] ) if (@$row_data_names);

  return $self;
}

sub sort {
  my $self = shift;
  my $cb   = shift;

  my $assay = $self->assay;
  my $nrow  = $self->nrow;
  my @idcs  = 0 .. ( $nrow - 1 );
  @idcs = sort { $cb->( $assay->[$a], $assay->[$b], $a, $b ) } @idcs;

  return $self->_reorder( \@idcs, @_ );
}

sub shuffle {
  my $self = shift;

  my $nrow = $self->nrow;
  my @idcs = 0 .. ( $nrow - 1 );
  @idcs = List::Util::shuffle(@idcs);

  return $self->_reorder( \@idcs, @_ );
}

sub _idx {
  my ( $self, $m, $name ) = @_;
  unless ($name) {
    return -1;
  }
  firstidx { $_ eq $name } @{ $self->$m };
}

sub row_idx {
  my ( $self, $name ) = @_;
  return -1 unless ($name);
  return firstidx { $_ eq $name } @{ $self->row_names };
}

sub col_idx {
  my ( $self, $name ) = @_;
  return -1 unless ($name);
  return firstidx { $_ eq $name } @{ $self->col_names };
}

sub transpose {
  my $self = shift;

  my @assay_t    = MapCarU { [@_] } @{ $self->{assay} };
  my @row_data_t = MapCarU { [@_] } @{ $self->{row_data} };
  my @col_data_t = MapCarU { [@_] } @{ $self->{col_data} };

  return __PACKAGE__->new(
    assay          => \@assay_t,
    row_names      => Clone::clone( $self->col_names ),
    row_data       => \@row_data_t,
    col_names      => Clone::clone( $self->row_names ),
    col_data       => \@col_data_t,
    col_data_names => Clone::clone( $self->row_data_names ),
    row_data_names => Clone::clone( $self->col_data_names ),
  );
}

sub make_consistent { die 'function not implemented, yet'; }

sub _idx_grep {
  my ( $self, $names, $cb ) = ( shift, shift, shift );

  return [ indexes { $_ =~ $cb } @$names ] if ref $cb eq 'Regexp';
  return [ indexes { $cb->($_) } @$names ];
}

sub row_idx_grep {
  my $self = shift;

  return $self->_idx_grep( $self->row_names, @_ );
}

sub col_idx_grep {
  my $self = shift;

  return $self->_idx_grep( $self->col_names, @_ );
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
  my $na_value = $self->na_value;

  $col_data_n //= [];
  push @{ $self->col_names }, @$names if ( $names && @$names );

  my $assay = $self->assay;
  my $ncol  = $self->ncol;

  if ( ref $data eq 'CODE' ) {
    for ( my $i = 0; $i < @$assay; $i++ ) {
      local $_ = $assay->[$i];
      push @{ $assay->[$i] }, $data->( $assay->[$i], $i );
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
    $col_data->[$i] //= [ ($na_value) x $ncol ];
    push @{ $col_data->[$i] }, @{ $col_data_n->[$i] // [ ($na_value) x $ncol_added ] };
  }

  return $self;
}

sub add_cols {
  die 'function not implemented, yet';
}

sub group {
  return shift->group_by_idcs(@_);
}

sub ncol {
  my $assay = shift->assay;
  return unless ( $assay && @$assay );

  return scalar @{ $assay->[0] };
}

sub nrow { scalar @{ shift->assay }; }

sub as_hash {
  my $self = shift;
  my %data;
  for my $entry (qw/assay col_data row_data row_names col_names row_data_names col_data_names meta_data/) {
    $data{$entry} = $self->$entry;
  }
  return \%data;
}

sub encode_as_json {
  my $self = shift;
  my $js   = JSON::XS->new->utf8->allow_nonref->indent(1);    #->canonical(1);
  return $js->encode( $self->as_hash );
}

sub json_spew {
  my ( $self, $f ) = @_;
  open my $fh, '>', $f or die "Can't open filehandle: $!";
  print $fh $self->encode_as_json;
  close $fh;
}

sub subset {
  my ( $self, $cb ) = @_;

  my $assay          = $self->assay;
  my $row_names      = $self->row_names;
  my $row_data       = $self->row_data;
  my $row_data_names = $self->row_data_names;

  my $idcs;
  if ( ref $cb eq 'CODE' ) {
    for ( my $i = 0; $i < @$assay; $i++ ) {
      local $_ = $assay->[$i];
      push @$idcs, $i if ( $cb->( $assay->[$i], $i ) );
    }
  } elsif ( ref $cb eq 'ARRAY' ) {
    $idcs = $cb;
  }

  my @row_names_new;
  my @assay_new;
  my @row_data_new;

  for my $i (@$idcs) {
    push @assay_new,     Clone::clone( $assay->[$i] );
    push @row_names_new, Clone::clone( $row_names->[$i] ) if ( $row_names && @$row_names );
    push @row_data_new,  Clone::clone( $row_data->[$i] ) if ( $row_data && @$row_data );
  }

  return __PACKAGE__->new(
    assay          => \@assay_new,
    row_names      => \@row_names_new,
    row_data       => \@row_data_new,
    row_data_names => Clone::clone( $self->row_data_names ),
    col_names      => Clone::clone( $self->col_names ),
    col_data       => Clone::clone( $self->col_data ),
    col_data_names => Clone::clone( $self->col_data_names ),
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

#TODO names to idcs as wantarray?
sub merge {
  my ( $se_x, $se_y, $param ) = @_;

  my $na_value = $se_x->na_value;

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
  my @keys = List::MoreUtils::uniq( keys(%$groups_x), keys(%$groups_y) );
  my $inv_by_y = _invert_idcs( $se_y->ncol, $idcs_y );
  for my $k (@keys) {
    my $data_x = $groups_x->{$k};
    my $data_y = $groups_y->{$k};

    next if ( $param{join} eq 'inner' && !( $data_x && $data_y ) );
    if ( $param{join} eq 'left' || $param{join} eq 'full' ) {
      next if ( $param{join} eq 'left' && !$data_x );
      $data_y //= {
        idcs      => [-1],
        rows      => [ [ ($na_value) x $se_y->ncol ] ],
        key       => $data_x->{key},
        key_names => $data_x->{key_names},
        row_names => [],
        row_data  => []
      };
    }
    if ( $param{join} eq 'right' || $param{join} eq 'full' ) {
      next if ( $param{join} eq 'right' && !$data_y );
      my @row = ( ($na_value) x $se_x->ncol );
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
    my $col_data_nrow = max( ( scalar @$col_data_x ), ( scalar @$col_data_y ) );
    for ( my $i = 0; $i < $col_data_nrow; $i++ ) {
      my $cd_x = $col_data_x->[$i] // [ ($na_value) x $ncol_only_x ];
      my $cd_y = $col_data_y->[$i] // [ ($na_value) x $ncol_only_y ];
      push @col_data, [ @$cd_x, @{$cd_y}[@$inv_by_y] ];
    }
  }
  return __PACKAGE__->new(
    assay => \@assay_new,

    row_data => \@row_data_new,
    col_data => \@col_data,

    row_data_names => Clone::clone( $se_x->row_data_names ),
    col_data_names => Clone::clone( $se_x->col_data_names ),

    row_names => \@row_names_new,
    col_names => \@col_names,
  );
}

sub inconsistencies {
  my $self = shift;

  # TODO
  # check if assay is rectangular
  # check if row data is rectangular and has the

}

sub _is_rectangular_matrix {
  my $aoa = shift;

  return unless ( ref $aoa eq 'ARRAY' );
  my $rlen;
  for ( my $i = 0; $i < @$aoa; $i++ ) {
    my $row = $aoa->[$i];
    return unless ( $row && ref $row eq 'ARRAY' );
    $rlen = @$row unless ( defined($rlen) );
    return unless ( @$row == $rlen );
  }
  return 1;
}

sub clone {
  my $self = shift;
  return __PACKAGE__->new(
    assay          => Clone::clone( $self->assay ),
    col_data       => Clone::clone( $self->col_data ),
    col_names      => Clone::clone( $self->col_names ),
    row_names      => Clone::clone( $self->row_names ),
    row_data       => Clone::clone( $self->row_data ),
    row_data_names => Clone::clone( $self->row_data_names ),
    col_data_names => Clone::clone( $self->col_data_names ),
    meta_data      => Clone::clone( $self->meta_data ),
  );
}

sub rbind {
  my $self           = shift;
  my $row_elems      = shift;
  my $names          = shift;
  my $row_data_elems = shift;

  my $col_names      = $self->col_names;
  my $row_data_names = $self->row_data_names;

  my @rows;
  for my $o (@$row_elems) {
    if ( ref $o eq 'ARRAY' ) {
      push @rows, $o;
    } else {
      push @rows, [ @{$o}{@$col_names} ];
    }
  }
  my @row_data;
  for my $o (@$row_data_elems) {
    if ( ref $o eq 'ARRAY' ) {
      push @row_data, $o;
    } else {
      push @row_data, [ @{$o}{@$col_names} ];
    }
  }

  return $self->_rbind( \@rows, $names, \@row_data );
}

sub _rbind {
  my ( $self, $rows, $names, $row_data ) = @_;

  my $nrow = $self->nrow;
  # FIXME check if all input params have the same length

  push @{ $self->assay }, @$rows;

  if ( $names && @$names ) {
    $self->row_names->[ $nrow - 1 ] //= undef if ( $nrow > 0 );
    push @{ $self->row_names }, @$names;
  }

  if ( $row_data && @$row_data ) {
    $self->row_data->[ $nrow - 1 ] //= [] if ( $nrow > 0 );
    push @{ $self->row_data }, @$row_data;
  }

  return $self;
}

sub dim {
  my $self = shift;
  return ( $self->nrow, $self->ncol );
}

sub _max_dim {
  my $aoa = shift;

  return unless ( ref $aoa eq 'ARRAY' );
  my $max_ncol;
  for ( my $i = 0; $i < @$aoa; $i++ ) {
    my $row = $aoa->[$i];
    return unless ( $row && ref $row eq 'ARRAY' );
    $max_ncol = @$row if ( !defined($max_ncol) || @$row > $max_ncol );
  }
  return unless ( defined $max_ncol );

  return [ scalar(@$aoa), $max_ncol ];
}

sub _na_fill_2d {
  my $data     = shift;
  my $dim      = shift // [ 0, 0 ];
  my $na_value = shift // $NA_VALUE;

  return unless ( $data && ref $data eq 'ARRAY' );
  my $dim_data = _max_dim($data);
  return unless ($dim);
  my $nrow = max( $dim_data->[0], $dim->[0] );
  my $ncol = max( $dim_data->[1], $dim->[1] );

  for ( my $i = 0; $i < $nrow; $i++ ) {
    $data->[$i] //= [];
    next if ( @{ $data->[$i] } == $ncol );
    for ( my $j = 0; $j < $ncol; $j++ ) {
      $data->[$i][$j] //= $na_value;
    }
  }
  return $data;
}

sub _na_fill_1d {
  my $data     = shift;
  my $dim      = shift;
  my $na_value = shift // $NA_VALUE;

  return unless ( $data && ref $data eq 'ARRAY' );
  my $len = max( $dim, @$data );
  return unless ($len);

  for ( my $i = 0; $i < $len; $i++ ) {
    $data->[$i] //= $na_value;
  }
  return $data;
}

sub add_rows { shift->rbind(@_) }

sub aggregate {
  return shift->aggregate_by_idcs(@_);
}

sub aggregate_by_idcs {
  my ( $self, $idcs, $cb, $col_names ) = @_;

  my $row_groups = $self->group($idcs);

  my @agg_assay;
  my @agg_row_names;
  my @agg_row_data;
  my @agg_row_data_names;
  for my $v ( values %$row_groups ) {
    local $_ = $v;
    my ( $row, $row_name, $row_data, $row_data_name )
      = $cb->( $v->{key}, $v->{rows}, $v->{idcs} );
    push @agg_assay,          $row           if ( defined($row) );
    push @agg_row_names,      $row_name      if ( defined($row_name) );
    push @agg_row_data,       $row_data      if ( defined($row_data) );
    push @agg_row_data_names, $row_data_name if ( defined($row_data_name) );
  }

  return __PACKAGE__->new(
    assay          => \@agg_assay,
    col_names      => $col_names // [],
    row_names      => \@agg_row_names,
    row_data_names => \@agg_row_data_names,
    row_data       => \@agg_row_data,
  );
}

sub col_names_to_idcs {
  my $self  = shift;
  my @names = @_;
  return unless (@names);
  @names = @{ $names[0] } if ( @names == 1 && ref $names[0] eq 'ARRAY' );

  my @idcs = map { $self->col_idx($_) } @names;
  die "could not find all idcs " . jon( ", ", @names ) if ( any { $_ < 0 } @idcs );
  return \@idcs;
}

sub aggregate_by_names {
  my ( $self, $names, $cb, $col_names ) = @_;
  my $idcs = $self->col_names_to_idcs($names);

  return $self->aggregate_by_idcs( $idcs, $cb, $col_names );
}

sub group_by_idcs {
  my ( $self, $idcs, $args ) = @_;

  my $assay          = $self->assay;
  my $row_names      = $self->row_names;
  my $row_data       = $self->row_data;
  my $row_data_names = $self->row_data_names;
  my %groups;
  my @key_names = @{ $self->col_names }[@$idcs];
  for ( my $i = 0; $i < @$assay; $i++ ) {
    my @key = @{ $assay->[$i] }[@$idcs];
    my $key = join( $;, @key );
    $groups{$key} //= {
      idcs           => [],
      rows           => [],
      key            => \@key,
      key_names      => \@key_names,
      row_names      => [],
      row_data       => [],
      row_data_names => $row_data_names
    };
    push @{ $groups{$key}{idcs} },      $i;
    push @{ $groups{$key}{rows} },      $assay->[$i];
    push @{ $groups{$key}{row_names} }, $row_names->[$i] if ( $row_names && @$row_names );
    push @{ $groups{$key}{row_data} },  $row_data->[$i] if ( $row_data && @$row_data );
  }
  return \%groups;
}

#if ( $uniq && !defined($vidx) ) {
#$map{$k} = 1;
#} elsif ( not defined $vidx ) {
#$map{$k}++;
#} elsif ($uniq) {
#confess "strict mode: two times the same key $k" if ( $is_strict && defined( $map{$k} ) );
#$map{$k} = ( ref $vidx ? [ @{$r}[@$vidx] ] : ( $vidx eq 'all' ? $r : $r->[$vidx] ) );
#} else {
#$map{$k} //= [];
#push @{ $map{$k} }, ( ref $vidx ? [ @{$r}[@$vidx] ] : ( $vidx eq 'all' ? $r : $r->[$vidx] ) );
#}
#}
#}
#return \%map;
#}

sub group_by_names {
  my $self  = shift;
  my $names = shift;

  my $idcs = $self->col_names_to_idcs($names);

  return $self->group( $idcs, @_ );
}

sub names_to_idcs {
  return shift->col_names_to_idcs(@_);
}

sub c2i {
  my $i = 0;
  my %I = ( map { $_ => $i++ } @{ shift->col_names } );

  return unless (%I);
  return wantarray ? %I : \%I;
}

sub col_idx_map { shift->c2i }

sub row_idx_map {
  my $i = 0;
  my %I = ( map { $_ => $i++ } @{ shift->row_names } );

  return unless (%I);
  return wantarray ? %I : \%I;
}

sub col_rename {
  my ( $self, $old, $new ) = @_;

  my $idx = $self->col_idx($old);
  die if ( $idx < 0 );
  $self->col_names->[$idx] = $new;
  return $self;
}

sub row_apply {
  my ( $self, $cb ) = @_;

  my @res;
  my $assay = $self->assay;
  for ( my $i = 0; $i < @$assay; $i++ ) {
    local $_ = $assay->[$i];
    push @res, $cb->( $assay->[$i], $i );
  }
  return \@res;
}

sub element_apply {
  my ( $self, $cb ) = @_;

  my @res;
  my $assay = $self->assay;
  for ( my $i = 0; $i < @$assay; $i++ ) {
    my $j = 0;
    my @row_res = map { $cb->( $_, $i, $j++ ) } @{ $assay->[$i] // [] };

    push @res, \@row_res;
  }
  return \@res;
}

sub col_apply {
  my ( $self, $cb ) = @_;

  my @res;
  my @assay_t = MapCarU { [@_] } @{ $self->{assay} };

  for ( my $i = 0; $i < @assay_t; $i++ ) {
    local $_ = $assay_t[$i];
    push @res, $cb->( $assay_t[$i], $i );
  }
  return \@res;
}

sub apply {
  my ( $self, $dir, $cb, @args ) = @_;

  if ( $dir eq 'r' || $dir == 1 ) {
    return $self->row_apply( $cb, @args );
  } elsif ( $dir eq 'c' || $dir == 2 ) {
    return $self->col_apply( $cb, @args );
  } elsif ( $dir eq 'rc' || $dir eq 'cr' || $dir == 3 ) {
    return $self->element_apply( $cb, @args );
  }
}

sub slice_by_idcs {
  my ( $self, $idcs ) = @_;

  my @assay_new = map { [ @{$_}[@$idcs] ] } @{ $self->assay };
  my @new_colnames;
  @new_colnames = @{ $self->col_names }[@$idcs] if ( $self->has_col_names );

  my @new_coldata;
  @new_coldata = map { [ @{$_}[@$idcs] ] } @{ $self->col_data } if ( $self->has_col_data );

  return __PACKAGE__->new(
    assay          => \@assay_new,
    row_names      => Clone::clone( $self->row_names ),
    row_data       => Clone::clone( $self->row_data ),
    col_names      => \@new_colnames,
    col_data       => \@new_coldata,
    row_data_names => Clone::clone( $self->row_data_names ),
    col_data_names => Clone::clone( $self->col_data_names ),
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

sub each {
  my $self = shift;
  $self->apply( 1, @_ );
  return $self;
}

sub uniq {
  my $self = shift;

  my %seen;
  return $self->subset(
    sub {
      return if ( $seen{ join $;, @$_ }++ );
      return 1;
    }
  );
}

# sub grep -> subset
# from Mojo::Collection
# first
# last

1;

__END__

=head1 NAME

Bio::Gonzales::SummarizedExperiment - represent experimental matrix-like data (assay) with features and sample info

=head1 SYNOPSIS


=head1 DESCRIPTION

L<http://bioconductor.org/packages/devel/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html>

=head1 ATTRIBUTES

=head2 assay

    my $assay = $se->assay;

Return the assay of the summarized experiment.

=head2 col_data

    my $col_data = $se->col_data;
    $se->col_data(\@col_data);

=head2 row_data

=head2 row_names

=head2 col_names

=head2 row_data_names

=head2 col_data_names

=head2 meta_data

=head1 METHODS

=head2 data

    my $assay = $se->data;

A alias for assay.

=head2 add_col

=head2 add_cols

=head2 add_rows

=head2 aggregate

=head2 C<< $se = $se->aggregate_by_idcs(\@idcs, sub { ... }, \@col_names)

The callback gets passed the grouping keys, rows and row indices. C<$_> is set to the
group has that comes from the (internally used) C<< $se->group >> function.

    sub {
      my ($key, $rows, $row_idcs) = @_;
      my $group = $_;
    }

=head2 C<< $se = $se->aggregate_by_names(\@names, sub { ... }, \@col_names)

=head2 apply

=head2 as_hash

=head2 cbind

=head2 clone

=head2 col_apply

=head2 col_idx

=head2 col_idx_map

    my $I = $se->col_idx_map;
    my %I = $se->col_idx_map;

Returns a hash that maps the column names to their column index. col_idx_map is context
sensitve and returns a hash in list context and a hash reference in scalar context.

=head2 col_idx_match

=head2 col_names_to_idcs

=head2 col_rename


=head2 dim

=head2 each

=head2 extract_col_by_idx

=head2 extract_col_by_name

=head2 group

=head2 group_by_idcs

=head2 group_by_names

=head2 has_col_data

=head2 has_col_names
=head2 has_row_data
=head2 has_row_names
=head2 header
=head2 header_idx
=head2 header_idx_match
=head2 inconsistencies
=head2 json_spew
=head2 make_consistent
=head2 merge

Merge two SummarizedExperiment objects.

    use Bio::Gonzales::SummarizedExperiment;
    my $se_x = Bio::Gonzales::SummarizedExperiment->new(
      assay => [ [ 1, "homer", "simpson" ], [ 2, "bart", "simpson" ], [ 3, "lisa simpson" ] ],
      col_names => [qw(user_id first_name surname)]
    );
    
    my $se_y = Bio::Gonzales::SummarizedExperiment->new(
      assay => [ [ 1, 120 ], [ 2, 20 ] ],
      col_names => [qw(user_id weight_kg)]
    );
    
    # inner join by default
    my $merged_se = $se_x->merge($se_y, { by => [ 'user_id' ] });
    
    # user_id first_name surname weight_kg
    # 1       homer      simpson 120
    # 2       bart       simpson 20
    # Lisa is missing, because the $se_y lacks weight information.

=head2 names_to_idcs
=head2 ncol
=head2 nrow
=head2 rbind
=head2 row_apply
=head2 row_idx
=head2 row_idx_map
=head2 row_idx_match
=head2 shuffle

=head2 slice_by_idcs

    $se->slice_by_idcs(\@idcs);
    $se->slice_by_idcs([0,5,13]);

Extract a column-"slice" from the summarized experiment. The indices select the columns.

=head2 slice_by_names

=head2 slurp_assay
    
    my $se = Bio::Gonzales::SummarizedExperiment->slurp_assay($source, \%params);
    my $se = Bio::Gonzales::SummarizedExperiment->slurp_assay("data.csv", { header => 1, sep => ';' });

Create a new summarized experiment from matrix/tabular data.

=head2 sort
=head2 spew_assay
=head2 subset
=head2 encode_as_json
=head2 transpose
=head2 uniq

=head1 LIMITATIONS

=head1 NOTES

By convention,

=over 4

=item * constructor or function arguments ending in C<?> are optional

=item * methods ending in C<!> will modify the object it is called on

=back

=head1 SEE ALSO

=head1 AUTHOR

jw bargsten, C<< <jwb at cpan dot org> >>

=cut


