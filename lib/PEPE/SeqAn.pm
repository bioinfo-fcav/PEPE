package PEPE::SeqAn;

use XSLoader;
use vars qw($VERSION @ISA);
 
BEGIN {
   @ISA = qw( );

   $VERSION = '0.01';
 
   # Put Perl code used in the BOOT: section here
 
   XSLoader::load __PACKAGE__, $VERSION;
}
 
# Put Perl code used in onBOOT() function here; calls to XSUBs are
# prototype-checked.
 
onBOOT;
 
# Put Perl initialization code assuming that XS is initialized here



