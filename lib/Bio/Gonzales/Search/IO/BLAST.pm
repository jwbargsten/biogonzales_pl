package Bio::Gonzales::Search::IO::BLAST;

use warnings;
use strict;
use Carp;

use 5.010;

use base 'Exporter';
use Bio::Gonzales::Util::File qw/basename regex_glob is_archive/;
use List::Util qw/min max/;
use File::Temp qw/tempdir tempfile mktemp/;
use Bio::Gonzales::Seq::IO qw(faiterate faspew);
use Capture::Tiny qw/capture_merged/;
use Path::Tiny;
use File::Spec::Functions qw/catfile/;
use String::ShellQuote;

use Params::Validate qw/validate/;
our ( @EXPORT, @EXPORT_OK, %EXPORT_TAGS );
# VERSION

@EXPORT      = qw();
%EXPORT_TAGS = ();
@EXPORT_OK   = qw(makeblastdb);

sub makeblastdb {
  my %c = validate(
    @_,
    {
      seq_file     => 1,
      title        => 0,
      parse_seqids => { default => 0 },
      hash_index   => { default => 1 },
      alphabet     => 1,
      db_prefix    => 0,
      wd           => 0,
    }
  );

  $c{wd} ||= './';
  my $unlink;
  my $seqf     = shell_quote( $c{seq_file} );
  my $basename = basename( $c{seq_file} );
  if ( my $type = is_archive($seqf) ) {
    say STDERR "$seqf is an archive, extracting first ...";

    my $tmp_f = mktemp( catfile( $c{wd}, 'tempXXXXXX' ) );

    if ( $type eq 'gzip' or $type eq 'bgzf' ) {
      system("gzip -dc <$seqf >$tmp_f") == 0 or die "system failed: $?";
    } elsif ( $type eq 'bzip2' ) {
      system("bzip2 -dc <$seqf >$tmp_f") == 0 or die "system failed: $?";
    } else {
      confess("archive type $type not supported");
    }

    $unlink = 1;
    $seqf   = $tmp_f;
    say STDERR "extraction finished. making blast DB";
    $basename = basename($basename);    # remove 2nd extension, e.g. a.fa.gz -> a.fa -> a
  }

  my @cmd = 'makeblastdb';
  push @cmd, '-in',    $seqf;
  push @cmd, '-title', $basename;
  push @cmd, '-parse_seqids' if ( $c{parse_seqids} );
  push @cmd, '-hash_index'   if ( $c{hash_index} );

  if    ( $c{alphabet} =~ /^(?:a|p)/ )   { push @cmd, '-dbtype', 'prot' }
  elsif ( $c{alphabet} =~ /^(?:n|d|r)/ ) { push @cmd, '-dbtype', 'nucl' }

  $c{db_prefix} //= $basename;
  $c{db_prefix} .= '.bdb';
  my $db_name = File::Spec->catfile( $c{wd}, $c{db_prefix} );
  push @cmd, '-out', $db_name;

  my @existing_db_files = regex_glob( $c{wd}, qr/^\Q$c{db_prefix}.\E[np]\w\w$/ );
  if ( @existing_db_files > 0 ) {
    my $oldest_db_file_age = min( map { ( stat $_ )[9] } @existing_db_files );
    my $seq_file_age = ( stat $c{seq_file} )[9];

    unless ( $seq_file_age > $oldest_db_file_age ) {
      #sequence file is older than db, so do noting
      say STDERR "sequence file $c{seq_file} is older than $db_name";
      say STDERR "Skipping blast db creation";
      return $db_name;
    }
  }

  say STDERR "Creating blast db:";
  say STDERR join " ", @cmd;

  my $merged = capture_merged { system(@cmd) };

  say STDERR $merged;

  unlink $seqf if ($unlink);
  return $db_name;
}

1;
__END__

=head1 NAME

Bio::Gonzales::Search::IO::BLAST

=head1 SYNOPSIS

    Bio::Gonzales::Search::IO::BLAST qw(makeblastdb instantblast)

=head1 DESCRIPTION

=head1 OPTIONS

=head1 SUBROUTINES

=over 4

=item B<< $db_location = makeblastdb(\%config) >>

Creates a blast database with the config options supplied. Config options are:

    %config = (
        'seq_file!'  => undef,
        title        => 'basename of seq_file',
        parse_seqids => 1,
        hash_index   => 1,
        'alphabet!'  => 'n(ucleotide)? || p(rotein)? || d(na)? || a(a)?',
        db_prefix    => 'basename of seq_file.bdb',
        wd           => './',
    );

Options with C<!> are required.

=back

=head1 SEE ALSO

=head1 AUTHOR

jw bargsten, C<< <joachim.bargsten at wur.nl> >>

=cut
