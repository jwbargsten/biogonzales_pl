#!/usr/bin/env perl
# created on 2015-01-22

use warnings;
use strict;
use 5.010;

use Bio::Gonzales::Matrix::IO qw(:DEFAULT xlsx_slurp);

my $in_f   = shift;
my $prefix = shift;

my $d = xlsx_slurp($in_f);
my $i = 1;
for my $s (@$d) {
  if ( @$s > 0 ) {
    mspew( "$prefix$i.tsv", $s );
  }
  $i++;
}
