# Before 'make install' is performed this script should be runnable with
# 'make test'. After 'make install' it should work as 'perl SeqAn.t'

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use strict;
use warnings;

use Test::More tests => 8;
BEGIN { use_ok('PEPE::SeqAn') };

#########################

# Insert your test code below, the Test::More module is use()ed here so read
# its man page ( perldoc Test::More ) for help writing this test script.

my $seqan = new PEPE::SeqAn();

isa_ok($seqan,'PEPE::SeqAn');

my ($score1, $aln12, $aln11) = split(';', $seqan->Align2Seq("TAAGTNNGAGGAATTGAAAAAAAAATA", "TTACCGATAAACGTAGAGACAAAGTCGGACGGGGTTGAAAAAAAATAGAAGAGACCAAAAAGTAGAGAAGAAGAGACGTAGAGaCAGTAACGGTAGAGAgAGaAAAACCC"));

cmp_ok($score1, '==', 37, "Check alignment 1 score");
cmp_ok($aln11,'eq','--------------------TAAGTNNGA-GGAATTGAAAAAAAAATA---------------------------------------------------------------', "Check aligned sequence 1 in alignment 1");
cmp_ok($aln12,'eq','TTACCGATAAACGTAGAGACAAAGTCGGACGGGGTTG-AAAAAAAATAGAAGAGACCAAAAAGTAGAGAAGAAGAGACGTAGAGACAGTAACGGTAGAGAGAGAAAAACCC', "Check aligned sequence 2 in alignment 1");

my ($score2, $aln22, $aln21) = split(';', $seqan->Align2Seq("AAAGTNNGAGGAATTGAAAAAATA", "TTACCGATAAACGTAGAGACAAAGTCGGACGGGGTTGAAAAAAAATAGAAGAGACCAAAAAGTAGAGAAGAAGAGACGTAGAGaCAGTAACGGTAGAGAgAGaAAAACCC"));

cmp_ok($score2, '==', 35, "Check alignment 2 score");
cmp_ok($aln21,'eq','--------------------AAAGTNNGA-GGAATTGAAAAAATA-----------------------------------------------------------------', "Check aligned sequence 1 in alignment 2");
cmp_ok($aln22,'eq','TTACCGATAAACGTAGAGACAAAGTCGGACGGGGTTGAAAAAAAATAGAAGAGACCAAAAAGTAGAGAAGAAGAGACGTAGAGACAGTAACGGTAGAGAGAGAAAAACCC', "Check aligned sequence 2 in alignment 2");
