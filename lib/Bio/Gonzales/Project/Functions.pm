package Bio::Gonzales::Project::Functions;

use warnings;
use strict;
use Carp;

use 5.010;

use File::Spec::Functions qw/catfile/;
use Bio::Gonzales::Project;
use Carp;
use Bio::Gonzales::Util::Cerial;
use Parallel::ForkManager;
use Bio::Gonzales::Util qw/ask sys_pipe_fatal sys_fmt/;

use base 'Exporter';
our (@EXPORT, @EXPORT_OK, %EXPORT_TAGS);
# VERSION

@EXPORT = qw(
  catfile
  nfi
  analysis_version
  path_to
  analysis_path
  gonzlog
  gonzconf
  iof
  gonzc
  gonzl
  gonz_iterate
  gonzsys
  analysis_name
  gsx
  gonzsys_pipe
  gsxp
  yes_or_die
);
%EXPORT_TAGS = ();
@EXPORT_OK   = qw();

sub _bgp {
  state $bgp = Bio::Gonzales::Project->new();
}

sub analysis_version { _bgp->analysis_version(@_) }
sub path_to          { _bgp->path_to(@_) }
sub nfi              { _bgp->nfi(@_) }
sub iof              { _bgp->conf(@_) }
sub gonzconf         { _bgp->conf(@_) }
sub gonzc            { _bgp->conf(@_) }
sub analysis_path    { _bgp->analysis_path(@_) }
sub analysis_name    { _bgp->analysis_name(@_) }

sub gonzlog { confess "deprecated call syntax, use gonzlog->info" if (@_ > 0 && $_[0]); _bgp->log() }
sub gonzl   { confess "deprecated call syntax, use gonzl->info"   if (@_ > 0 && $_[0]); _bgp->log() }

sub gonzsys {
  _bgp->log->info("(exec) " . join(" ", @_));
  system(@_) == 0 or _bgp->log->fatal_confess("system failed: $?");
}

*gsx = \&gonzsys;

sub gonzsys_pipe {
  _bgp->log->info("(sh) " . sys_fmt(@_));
  sys_pipe_fatal(@_);
}

*gsxp = \&gonzsys_pipe;

sub yes_or_die {
  my $q = shift;
  unless (ask($q)) {
    _bgp->log->fatal_croak("answered no to \"$q\"");
  }
}

sub gonz_iterate {
  my ($src, $code, $conf) = @_;
  $conf->{processes} //= 4;
  my $data;
  my $ref_type = ref($src);
  if (!$ref_type || ($ref_type ne 'ARRAY' && $ref_type ne 'HASH')) {
    $data = jslurp($src);
  } else {
    $data = $src;
  }

  if ($conf->{test}) {
    $code = sub { say jfreeze(\@_); return };
  }

  my $pm = Parallel::ForkManager->new($conf->{processes});

  my @result_all;
  $pm->run_on_finish(
    sub {
      my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $res) = @_;

      if (defined($res) && @$res > 0) {
        push @result_all, $res;
      }
    }
  );

  if (ref($data) eq 'ARRAY') {
    for (my $i = 0; $i < @$data; $i++) {
      $pm->start and next;     # do the fork
      my $res = $code->($i, $data->[$i]);
      $pm->finish(0, $res);    # do the exit in the child process
    }
    $pm->wait_all_children;
  } elsif (ref($data) eq 'HASH') {
    for my $k (keys %$data) {
      $pm->start and next;     # do the fork
      my $res = $code->($k, $data->{$k});
      $pm->finish(0, $res);    # do the exit in the child process
    }
    $pm->wait_all_children;

  }
  return \@result_all;
}
1;

__END__

=head1 NAME

Bio::Gonzales::AV - analysis project utils

=head1 SYNOPSIS

    use Bio::Gonzales::AV
      qw(
        catfile
        nfi
        analysis_version
        path_to
        analysis_path
        gonzlog
        gonzconf
        iof
        gonzc
        gonzl
        gonz_iterate
        gonzsys
        analysis_name
        gsx
        gonzsys_pipe
        gsxp
        yes_or_die
      );

=head1 SUBROUTINES

=over 4

=item B<< msg(@stuff) >>

say C<@stuff> to C<STDERR>.

=item B<< path_to($filename) >>

Locate the root of the project and prepend it to the C<$filename>.

=item B<< iof() >>

get access to the IO files config file. Use like

    my $protein_files = iof()->{protein_files}

=item B<< nfi($filename) >>

Prepend the current analysis version diretory to the filename.


=item B<< catfile($path, $file) >>

make them whole again...

=back

=head1 SEE ALSO

=head1 AUTHOR

jw bargsten, C<< <joachim.bargsten at wur.nl> >>

=cut
