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

has [qw/assay col_data row_data row_names col_names/] => ( is => 'rw', default => []);

sub data { shift->assay(@_) };
sub header { shift->col_names(@_) };

sub _idx {
  my ($self, $m, $name) = @_;
  firstidx { $_ eq $name } @{$self->$m};
}

sub row_idx { shift->_idx('row_names', @_) }

sub col_idx { shift->_idx('col_names', @_) }

sub header_idx { shift->col_idx(@_) }


sub _idx_match {
  my ($self, $m, $rex) = @_;
  [ indexes { $_ =~ /$rex/ } @{$self->$m} ]
}

sub row_idx_match { shift->_idx_match('row_names', @_) }

sub col_idx_match { shift->_idx_match('col_names', @_) }

sub header_idx_match {
  shift->col_idx_match(@_);
}


sub add_cols {
  my ($self, $cols, $names) = @_;

  push @{$self->col_names}, $name if(defined($name));

  my $assay = $self->assay;
  die unless(@$col == @$assay);

  # FIXME col kein eine fkt sein
  for(my $i = 0; $i < @$assay; $i++) {
    push @{$assay->[$i]}, $col->[$i];
  }

  return $self;
}

sub add_rows {
  my ($self, $rows, $names) = @_;

  push @{$self->row_names}, $name if(defined($name));

  push @{$self->assay}, $row;
  return $self;
}

sub row_ref_apply {

}

sub col_apply {

}

sub aggregate_by_idcs {
  my ($self, $idcs, $fun, $col_names) = @_;

  sub {
    my ( $data,  $row_idcs,$col_idcs, $row_names, $col_names) = @_;
    return ($row, $row_name);
  }
  
  return 

}

sub aggregate_by_names {

}

1;
