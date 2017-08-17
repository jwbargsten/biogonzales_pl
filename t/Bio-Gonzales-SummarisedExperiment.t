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

BEGIN { use_ok("Bio::Gonzales::SummarizedExperiment"); }

done_testing();
