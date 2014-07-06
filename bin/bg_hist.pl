#!/usr/bin/env perl
# created on 2014-07-06

use warnings;
use strict;
use 5.010;
use Bio::Gonzales::Stat::Util qw/hist_text/;

my @values = map { chomp; $_ } <STDIN>;

print hist_text(\@values);
