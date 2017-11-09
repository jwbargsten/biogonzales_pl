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

ok( Bio::Gonzales::SummarizedExperiment::_is_rectangular_matrix(  [ [],  [], [] ] ) );
ok( !Bio::Gonzales::SummarizedExperiment::_is_rectangular_matrix( [ [1], [], [] ] ) );
ok( !Bio::Gonzales::SummarizedExperiment::_is_rectangular_matrix( [ 1,   [], [] ] ) );

my $se_base = Bio::Gonzales::SummarizedExperiment->new(
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
  col_data_names => [qw/cdn1 cdn2/],
  row_data_names => [qw/rdn1 rdn2 rdn3/],
);

my $se_t = $se_base->clone;

$se_t = $se_base->transpose;

is_deeply( $se_t->col_data_names, [qw/rdn1 rdn2 rdn3/] );
is_deeply( $se_t->assay,          [ [qw/a 1 d/], [qw/b 2 e/], [qw/c 3 f/] ] );
is_deeply( $se_t->row_names,      [qw/cn1 cn2 cn3/] );
is_deeply( $se_t->col_data,       [ [qw/c1d1 c1d2/], [qw/c2d1 c2d2/], [qw/c3d1 c3d2/] ] );

my $se = $se_base->clone;

$se->add_col( [qw/u v w/], 'cn4' );

is_deeply( $se->col_names, [qw/cn1 cn2 cn3 cn4/] );
is_deeply( $se->assay, [ [qw/a b c u/], [ 1, 2, 3, 'v' ], [qw/d e f w/], ] );
is( $se->col_idx('cn4'),   3 );
is( $se->col_data->[0][3], 'NA' );
is( $se->col_data->[1][3], 'NA' );

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
    [ 2, 2, 3 ],             #
    [qw/d e f/],
    [qw/a b c/],             #
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

is_deeply(
  $se2->subset( [ 0, 1 ] )->assay,
  [
    [ 2, 2, 3 ],             #
    [qw/d e f/],
  ]
);

is_deeply( $se2->subset( sub { $_->[0] eq '2' } )->assay, [ [ 2, 2, 3 ], ] );
is_deeply( $se2->subset( sub { $_->[0] eq '2' } )->row_data,  [ [qw/r1d1 r1d2 r1d3/] ] );
is_deeply( $se2->subset( sub { $_->[0] eq '2' } )->row_names, [qw/rn1/] );

my $sem;
$sem = $se1->merge( $se2, { by_x => [qw/cn1/], by_y => [qw/ckn1/] } )->sort( sub { $_[0][0] cmp $_[1][0] } );

is_deeply( $sem->assay, [ [qw/a b c b c/], [qw/d e f e f/] ] );

$sem = $se1->merge( $se2, { join => 'left', by_x => [qw/cn1/], by_y => [qw/ckn1/] } )
  ->sort( sub { $_[0][0] cmp $_[1][0] } );
is_deeply( $sem->assay, [ [qw/1 2 3 NA NA/], [qw/a b c b c/], [qw/d e f e f/] ] );

$sem = $se1->merge( $se2, { join => 'right', by_x => [qw/cn1/], by_y => [qw/ckn1/] } )
  ->sort( sub { $_[0][0] cmp $_[1][0] } );
is_deeply( $sem->assay, [ [qw/2 NA NA 2 3/], [qw/a b c b c/], [qw/d e f e f/] ] );

dies_ok { $se1->merge( $se2, { join => 'left', by_x => [qw/cn1 cn2/], by_y => [qw/ckn1/] } ) };

$sem = $se1->merge( $se2, { join => 'left', by_x => [qw/cn1 cn3/], by_y => [qw/ckn1 ckn3/] } )
  ->sort( sub { $_[0][0] cmp $_[1][0] } );
is_deeply( $sem->assay, [ [qw/1 2 3 NA/], [qw/a b c b/], [qw/d e f e/] ] );
is_deeply( $sem->col_names, [qw/cn1 cn2 cn3 ckn2/] );

my $se_slice = $se_base->slice_by_names( [qw/cn1 cn3/] );
is_deeply(
  $se_slice->assay,
  [
    [qw/a c/],    #
    [ 1, 3 ],     #
    [qw/d f/],
  ]
);

is_deeply( $se_slice->col_names, [qw/cn1 cn3/] );

my @data;
@data = ( [qw/1 2 3/], [qw/2/], [qw/3 2/] );
Bio::Gonzales::SummarizedExperiment::_na_fill_2d( \@data );
is_deeply( \@data, [ [qw/1 2 3/], [qw/2 NA NA /], [qw/3 2 NA/] ] );

@data = ( [qw/1 2 3/], [qw/2/], [qw/3 2/] );
Bio::Gonzales::SummarizedExperiment::_na_fill_2d( \@data, [ 5, 6 ] );
is_deeply( \@data,
  [ [qw/1 2 3 NA NA NA/], [qw/2 NA NA NA NA NA/], [qw/3 2 NA NA NA NA/], [ ('NA') x 6 ], [ ('NA') x 6 ] ] );

is_deeply(
  $se_slice->col_data,
  [
    [qw/c1d1 c3d1/],    #
    [qw/c1d2 c3d2/]
  ]
);

my $se_elem_apply = $se_base->clone;

my $i = 0;
my $res = $se_elem_apply->element_apply( sub { s/[a-z]/u/g; return $i++ } );

is_deeply( $se_elem_apply->assay, [ [ 'u', 'u', 'u' ], [ 1, 2, 3 ], [ 'u', 'u', 'u' ] ] );
is_deeply( $res,                  [ [ 0,   1,   2 ],   [ 3, 4, 5 ], [ 6,   7,   8 ] ] );

done_testing();
