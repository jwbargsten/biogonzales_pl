use warnings;
use v5.11;
use strict;

use IO::Handle ();
use Test2::V0;
use Data::Dumper;
use File::Temp qw/tempfile tempdir/;
use File::Spec::Functions;
use File::Copy;
use File::Path qw/make_path/;
use Bio::Gonzales::Util::Cerial;
use Path::Tiny;

use Bio::Gonzales::Project;
$ENV{ANALYSIS_VERSION} = '2017-08-08_test';

{
  no warnings 'redefine';

  sub Bio::Gonzales::Project::path_to {
    shift;
    return catfile( "/tmp", @_ );
  }
}

my $structure = {
  a => "__av__/wurst",
  b => [qw(~/a b __an__)],
  c => { a => '__data__' },
  d => 1,
};

my ( $tmp_fh, $tmp_f ) = tempfile();

yspew( $tmp_fh, $structure );
$tmp_fh->flush;

my $bgp = Bio::Gonzales::Project->new( config_file => $tmp_f );

is(
  $bgp->config,
  {
    a => "2017-08-08_test/wurst",
    b => [ "$ENV{HOME}/a", 'b', path(".")->realpath->basename ],
    c => { a => '/tmp/data' },
    d => 1,
  }
);
done_testing();
