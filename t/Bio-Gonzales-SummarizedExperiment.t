use warnings;
use v5.11;
use strict;

use IO::Handle ();
use Test::More;
use Data::Dumper;
use File::Temp qw/tempfile tempdir/;
use File::Spec::Functions;
use File::Copy;
use File::Path qw/make_path/;
use Test::Exception;

BEGIN { use_ok('Bio::Gonzales::SummarizedExperiment'); }

my $se = Bio::Gonzales::SummarizedExperiment->new(
  assay => [
    [qw/a b c/],    #
    [ 1, 2, 3 ],    #
    [qw/d e f/],
  ],
  row_names => [qw/rn1 rn2 rn3/],
  col_names => [qw/cn1 cn2 cn3/],
  row_data  => [
    [qw/r1d1 r1d2 r1d3/],    #
    [qw/r2d1 r2d2 r2d3/],
    [qw/r3d1 r3d2 r3d3/],
  ],
  col_data => [
    [qw/c1d1 c2d1 c3d1/],    #
    [qw/c1d2 c2d2 c3d2/]
  ],
);

$se->add_col( [qw/u v w/], 'cn4' );

is_deeply( $se->col_names, [qw/cn1 cn2 cn3 cn4/] );
is_deeply( $se->assay, [ [qw/a b c u/], [ 1, 2, 3, 'v' ], [qw/d e f w/], ] );
is( $se->col_idx('cn4'), 3 );
is($se->col_data->[0][3], 'NA');
is($se->col_data->[1][3], 'NA');

# MERGE
my $se1 = Bio::Gonzales::SummarizedExperiment->new(
  assay => [
    [qw/a b c/],    #
    [ 1, 2, 3 ],    #
    [qw/d e f/],
  ],
  row_names => [qw/rn1 rn2 rn3/],
  col_names => [qw/cn1 cn2 cn3/],
  row_data  => [
    [qw/r1d1 r1d2 r1d3/],    #
    [qw/r2d1 r2d2 r2d3/],
    [qw/r3d1 r3d2 r3d3/],
  ],
  col_data => [
    [qw/c1d1 c2d1 c3d1/],    #
    [qw/c1d2 c2d2 c3d2/]
  ],
);

my $se2 = Bio::Gonzales::SummarizedExperiment->new(
  assay => [
    [ 2, 2, 3 ],    #
    [qw/d e f/],
    [qw/a b c/],    #
  ],
  row_names => [qw/rn1 rn2 rn3/],
  col_names => [qw/ckn1 ckn2 ckn3/],
  row_data  => [
    [qw/r1d1 r1d2 r1d3/],    #
    [qw/r2d1 r2d2 r2d3/],
    [qw/r3d1 r3d2 r3d3/],
  ],
  col_data => [
    [qw/c1d1 c2d1 c3d1/],    #
    [qw/c1d2 c2d2 c3d2/]
  ],
);

my $sem;
$sem = $se1->merge($se2, { by_x => [qw/cn1/], by_y => [qw/ckn1/] })->sort(sub { $_[0][0] cmp $_[1][0] });

is_deeply($sem->assay, [ [ qw/a b c b c/ ], [ qw/d e f e f/ ] ]);

$sem = $se1->merge($se2, { join => 'left', by_x => [qw/cn1/], by_y => [qw/ckn1/] })->sort(sub { $_[0][0] cmp $_[1][0] });
is_deeply($sem->assay, [ [qw/1 2 3 NA NA/], [qw/a b c b c/ ], [ qw/d e f e f/ ] ]);

$sem = $se1->merge($se2, { join => 'right', by_x => [qw/cn1/], by_y => [qw/ckn1/] })->sort(sub { $_[0][0] cmp $_[1][0] });
is_deeply($sem->assay, [ [qw/2 NA NA 2 3/], [qw/a b c b c/ ], [ qw/d e f e f/ ] ]);

dies_ok { $se1->merge($se2, { join => 'left', by_x => [qw/cn1 cn2/], by_y => [qw/ckn1/] })};

$sem = $se1->merge($se2, { join => 'left', by_x => [qw/cn1 cn3/], by_y => [qw/ckn1 ckn3/] })->sort(sub { $_[0][0] cmp $_[1][0] });
is_deeply($sem->assay, [ [qw/1 2 3 NA/], [qw/a b c b/ ], [ qw/d e f e/ ] ]);
is_deeply($sem->col_names, [qw/cn1 cn2 cn3 ckn2/]);

done_testing();

