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

sub row_idx {
  my ($self, $name) = @_;

  firstidx { $_ eq $name } @{$self->row_names};
}

sub col_idx {
  my ($self, $name) = @_;

  firstidx { $_ eq $name } @{$self->col_names};
}

sub header_idx { shift->col_idx(@_) }

sub row_idx_match {
  my ($self, $rex) = @_;

  [ indexes { $_ =~ /$rex/ } @{$self->row_names} ]
}

sub col_idx_match {
  my ($self, $rex) = @_;

  [ indexes { $_ =~ /$rex/ } @{$self->col_names} ]
}

sub header_idx_match {
  shift->col_idx_match(@_);
}

1;
