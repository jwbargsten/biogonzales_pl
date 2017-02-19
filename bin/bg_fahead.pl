#!/usr/bin/env perl
# created on 2016-09-09

use warnings;
use strict;
use 5.010;
use IO::Handle ();

use Bio::Gonzales::Seq::Filter qw/clean_dna_seq/;

my $fait = faiterate(\*STDIN);

use Bio::Gonzales::Seq;

my $fait = faiterate(\*STDIN);
while ( my $so = $fait->() ) {
  clean_dna_seq($so);
  print $so;
}

