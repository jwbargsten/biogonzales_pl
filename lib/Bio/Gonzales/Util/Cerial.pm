package Bio::Gonzales::Util::Cerial;

use warnings;
use strict;
use Carp;

use v5.11;

use Exporter 'import';

use Bio::Gonzales::Util::File qw/open_on_demand/;
use Bio::Gonzales::Util qw/deep_value flatten/;

use Try::Tiny;
use YAML::XS;
use JSON::XS;
use Data::Dumper;
use Storable qw/nstore_fd fd_retrieve/;

# VERSION

our %EXPORT_TAGS = (
  'all' => [
    qw(
      ndjson_iterate ndjson_hash ndjson_slurp ndjson_spew
      ndjson_freeze ndjson_thaw
      ythaw yfreeze yslurp yspew
      jthaw jfreeze jslurp jspew
      stoslurp stospew
      )
  ],
  std => [
    qw(
      ythaw yfreeze yslurp yspew
      jthaw jfreeze jslurp jspew
      stoslurp stospew
      )
  ]
);
our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );
our @EXPORT    = ( @{ $EXPORT_TAGS{'std'} } );

BEGIN {
  *yfreeze = \&YAML::XS::Dump;
  *ythaw   = \&YAML::XS::Load;
  #*jfreeze = \&JSON::XS::encode_json;
  *jthaw = \&JSON::XS::decode_json;
}

our $JSON = JSON::XS->new->indent(1)->utf8->allow_nonref;

sub jfreeze {
  my $r;
  my @d = @_;
  try {
    $r = $JSON->encode(@d);
  }
  catch {
    confess Dumper \@d;
  };

}

sub _spew {
  my $dest = shift;
  my $data = shift;

  my ( $fh, $was_open ) = open_on_demand( $dest, '>' );
  binmode $fh, ':utf8' unless ( ref $fh eq 'IO::Zlib' );
  local $/ = "\n";

  print $fh $data;
  $fh->close unless $was_open;
}

sub _slurp {
  my $src = shift;
  my ( $fh, $was_open ) = open_on_demand( $src, '<' );
  binmode $fh, ':utf8' unless ( ref $fh eq 'IO::Zlib' );
  local $/ = "\n";

  my $data = do { local $/; <$fh> };

  $fh->close unless $was_open;
  return $data;
}

sub yslurp { return ythaw( _slurp(shift) ) }
sub jslurp { return jthaw( _slurp(shift) ) }
sub yspew  { my $file = shift; _spew( $file, yfreeze( $_[0] ) ) }
sub jspew  { my $file = shift; _spew( $file, jfreeze( $_[0] ) ) }

sub stospew {
  my $dest = shift;
  my $data = shift;

  my ( $fh, $was_open ) = open_on_demand( $dest, '>' );
  nstore_fd( $data, $fh );
  $fh->close unless ($was_open);
}

sub stoslurp {
  my $src = shift;
  my ( $fh, $was_open ) = open_on_demand( $src, '<' );
  my $data = fd_retrieve($fh);
  $fh->close unless $was_open;
  return $data;
}

sub ndjson_freeze {
  my $entries = shift;
  return unless (@$entries);
  state $js = JSON::XS->new->utf8->allow_nonref;
  return join( "\n", ( map { $js->encode($_) } @$entries ) ) . "\n";
}

sub ndjson_thaw {
  my $data = shift;

  return unless ($data);
  my @entries = split /\n/, $data;

  return unless (@entries);

  state $js = JSON::XS->new->utf8->allow_nonref;

  return [ map { $js->decode($_) } @entries ];
}

sub ndjson_hash {
  my $keys  = shift;
  my $files = shift;
  my $cfg   = shift;

  my %res;
  my $it = ndjson_iterate($files);
  while ( my $elem = $it->() ) {
    my $val = deep_value( $elem, $keys );
    if ( $cfg->{uniq} ) {
      die $val . " already exists" if ( $res{$val} );
      $res{$val} = $elem;
    } else {
      $res{$val} //= [];
      push @{ $res{$val} }, $elem;
    }
  }
  return \%res;
}

sub ndjson_slurp {
  my $it = ndjson_iterate(@_);

  my @res;
  while ( defined( my $elem = $it->() ) ) {
    push @res, $elem;
  }
  return \@res;
}

sub ndjson_spew {
  my ( $dest, $elems, $c ) = @_;

  my $js = JSON::XS->new->utf8->allow_nonref;
  $js = $js->canonical(1) if ( $c->{canonical} );
  my ( $fh, $fh_was_open ) = open_on_demand( $dest, '>' );

  try {
    for my $elem (@$elems) {
      say $fh $js->encode($elem);
    }
  }
  catch {
    confess "could not spew: $_";
  };

  $fh->close unless ($fh_was_open);
  return;
}

sub ndjson_iterate {
  my @srcs = flatten(@_);

  my $i = 0;
  my ( $fh, $fh_was_open ) = open_on_demand( $srcs[$i], '<' );

  return sub {
    while ( $i < @srcs ) {
      while ( my $record = <$fh> ) {
        next if ( !$record || $record =~ /^\s*$/ );
        my $data;
        try {
          $data = decode_json($record);
        }
        catch {
          warn "caught error: $_ JSON string >$record<";
        };
        confess "no valid data in record" unless ($data);

        return $data;
      }

      $fh->close unless ($fh_was_open);
      if ( ++$i >= @srcs ) {
        # return if next file does not exist
        return;
      }

      # open next file
      ( $fh, $fh_was_open ) = open_on_demand( $srcs[$i], '<' );
    }
    return;
  };

}

__END__

=head1 NAME

Bio::Gonzales::Util::Cerial - convenience functions for yaml and json IO

=head1 SYNOPSIS

    use Bio::Gonzales::Util::Cerial;

    # YAML IO
    my $yaml_string = yfreeze(\%data);
    my $data = ythaw($yaml_string);

    yspew($filename, \%data);
    my $data = yslurp($filename);

    # JSON IO
    my $json_string = jfreeze(\%data);
    my $data = jthaw($json_string);

    jspew($filename, \%data);
    my $data = jslurp($filename);


=head1 DESCRIPTION
    
=item B<< $yaml_string = yfreeze(\%data) >>

Serialize data structure as yaml string

=item B<< $data = ythaw($yaml_string) >>

UNserialize data structure from yaml string

=item B<< yspew($filename, \%data) >>

Serialize data structure as yaml string to a file

=item B<< my $data = yslurp($filename) >>

UNserialize data structure from yaml file

=item B<< my $json_string = jfreeze(\%data) >>

Serialize data structure as json string

=item B<< my $data = jthaw($json_string) >>

UNserialize data structure from json string

=item B<< jspew($filename, \%data) >>

Serialize data structure as json string to a file

=item B<< my $data = jslurp($filename) >>

UNserialize data structure from json file

=head1 EXPORT

The following functions are exported by default

    ythaw
    yfreeze
    yslurp
    yspew

    jthaw
    jfreeze
    jslurp
    jspew

=cut
