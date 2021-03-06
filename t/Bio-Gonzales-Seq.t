use warnings;
use Test::More;
use Data::Dumper;

BEGIN { use_ok("Bio::Gonzales::Seq"); }

my $seq = Bio::Gonzales::Seq->new(seq => 'aggct', id => 'test');
is($seq->revcom->seq, 'agcct');

my $seq2 = Bio::Gonzales::Seq->new(seq => 'ag.gct', id => 'test');

is($seq2->ungapped_seq, 'aggct');
done_testing();

