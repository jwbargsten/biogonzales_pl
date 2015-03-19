package Bio::Gonzales::Project::Functions;

use warnings;
use strict;
use Carp;

use 5.010;

use File::Spec::Functions qw/catfile/;
use Bio::Gonzales::Project;
use Carp;
use Bio::Gonzales::Util::Cerial;

use base 'Exporter';
our ( @EXPORT, @EXPORT_OK, %EXPORT_TAGS );
# VERSION

@EXPORT
  = qw(catfile nfi analysis_version path_to analysis_path gonzlog gonzconf iof $GONZLOG gonzc gonzl gonz_iterate);
%EXPORT_TAGS = ();
@EXPORT_OK   = qw();

my $bgp = Bio::Gonzales::Project->new();

our $GONZLOG = $bgp->log;

sub analysis_version { $bgp->analysis_version(@_) }
sub path_to          { $bgp->path_to(@_) }
sub gonzlog          { $bgp->log() }
sub gonzl            { $bgp->log() }
sub nfi              { $bgp->nfi(@_) }
sub iof              { $bgp->conf(@_) }
sub gonzconf         { $bgp->conf(@_) }
sub gonzc            { $bgp->conf(@_) }
sub analysis_path    { $bgp->analysis_path(@_) }

sub gonz_iterate {
  my ( $src, $code ) = @_;
  my $data;
  my $ref_type = ref($src);
  if ( !$ref_type || ( $ref_type ne 'ARRAY' && $ref_type ne 'HASH' ) ) {
    $data = jslurp($src);
  } else {
    $data = $src;
  }

  my @res;
  if ( ref($data) eq 'ARRAY' ) {
    for ( my $i = 0; $i < @$data; $i++ ) {
      push @res, $code->( $i, $data->[$i] );
    }

  } elsif ( ref($data) eq 'HASH' ) {
    for my $k ( keys %$data ) {
      push @res, $code->( $k, $data->{$k} );
    }

  }
  return \@res;
}
1;

__END__

=head1 NAME

Bio::Gonzales::AV - analysis project utils

=head1 SYNOPSIS

    use Bio::Gonzales::AV qw(catfile nfi $ANALYSIS_VERSION iof path_to analysis_path msg error debug);

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
