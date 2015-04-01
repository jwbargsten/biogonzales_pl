use warnings;
use 5.010;
use strict;

use Test::More;
use Data::Dumper;

BEGIN { use_ok('Bio::Gonzales::Project::Functions'); }

my @a = (0..10);
my @b;
gonz_iterate(\@a, sub { $b[$_[0]] = $_[1] * $_[1] });

my @b_ref = map { $_*$_ } @a;
is_deeply(\@b, \@b_ref);
done_testing();
