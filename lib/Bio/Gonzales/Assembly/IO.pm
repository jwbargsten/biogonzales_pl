package Bio::Gonzales::Assembly::IO;

use warnings;
use strict;
use Carp;

use 5.010;

use List::MoreUtils qw/zip/;
use Bio::Gonzales::Seq::IO qw/fahash/;
use Bio::Gonzales::Matrix::IO qw(mslurp);
use Bio::Gonzales::Util::File qw/slurpc/;

use base 'Exporter';
our ( @EXPORT, @EXPORT_OK, %EXPORT_TAGS );
# VERSION

@EXPORT      = qw();
%EXPORT_TAGS = ();
@EXPORT_OK   = qw(agpslurp agp2fasta);

our @AGP_COLUMN_NAMES = qw/
  object
  object_beg
  object_end
  part_number
  component_type
  component_id
  component_beg
  component_end
  orientation/;

sub INFO { say STDERR @_; }

sub agpslurp {
  my ($file) = @_;

  my @lines = slurpc( $file);

  my @agp;
  for my $l (@lines) {
    my @a = split /\t/, $l;

    #add last field in case somebody forgot to add it...
    push @a, '' if ( @a == 8 );

    #but if sth is really going wrong, break
    confess "error in agp file $file\n$l" unless ( @a == 9 );

    if    ( $a[8] eq '-' ) { $a[8] = -1 }
    elsif ( $a[8] eq '+' ) { $a[8] = 1 }
    else                   { $a[8] = 0 }

    push @agp, +{ zip @AGP_COLUMN_NAMES, @a };
  }

  return \@agp;
}

sub agp2fasta {
  my ( $agp, $seq, $out ) = @_;

  INFO("reading scf_seqs");
  my %scf_seqs = map { INFO("  $_"); ( fahash($_) ) } @$seq;

  INFO("reading agp data");
  my %agp_data;
  for my $agp_file (@$agp) {
    INFO("  $agp_file");
    my $data = agpslurp($agp_file);
    for my $e (@$data) {
      $agp_data{ $e->{object} } //= [];
      push @{ $agp_data{ $e->{object} } }, $e;

    }
  }

  INFO("processing agp data");
  while ( my ( $chr_id, $objs ) = each %agp_data ) {
    INFO("processing $chr_id");
    $objs = [ sort { $a->{part_number} <=> $b->{part_number} } @$objs ];

    open my $out_fh, '>', $out or confess "Can't open filehandle: $!";
    my $last_obj_end;
    say $out_fh ">$chr_id";
    for my $o (@$objs) {
      die "error" if ( $last_obj_end && ( $last_obj_end + 1 != $o->{object_beg} ) );
      if ( $o->{component_type} eq 'W' ) {
        INFO( "processing scf " . $o->{component_id} );
        unless ( exists( $scf_seqs{ $o->{component_id} } ) ) {
          die "could not find " . Dumper $o;
        }

        print $out_fh $scf_seqs{ $o->{component_id} }
          ->subseq( [ @{$o}{qw/component_beg component_end orientation/} ] )->seq;
      } elsif ( $o->{component_id} eq 'N' ) {

        my $len = $o->{object_end} - $o->{object_beg} + 1;
        INFO( "processing gap of length " . $len );
        print $out_fh 'N' x $len;
      } elsif ( $o->{component_id} eq 'U' ) {
        INFO("processing gap of UNKNOWN length 100");
        print $out_fh 'N' x 100

      } else {
        die "unknown component type" . Dumper $o;
      }
    }
    $out_fh->close;
  }

}

1;

__END__

=head1 NAME

Bio::Gonzales::Assembly::IO - assembly related stuff

=head1 SYNOPSIS

    use Bio::Gonzales::Assembly::IO qw(agpslurp);

=head1 DESCRIPTION

=head1 OPTIONS

=head1 SUBROUTINES

=head1 METHODS

=head1 SEE ALSO

=head1 AUTHOR

jw bargsten, C<< <joachim.bargsten at wur.nl> >>

=cut
