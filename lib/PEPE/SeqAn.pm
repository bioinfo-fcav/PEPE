package PEPE::SeqAn;

use 5.022001;
use strict;
use warnings;

require Exporter;

use vars qw($AUTOLOAD);

our @ISA = qw(Exporter);

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# This allows declaration	use PEPE::SeqAn ':all';
# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK
# will save memory.
our %EXPORT_TAGS = ( 'all' => [ qw(
	
) ] );

our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

our @EXPORT = qw(
	
);

our $VERSION = '0.01';

sub new {
	my ( $caller, @args ) = @_;

	my $self = bless {}, ref($caller) || $caller;

	require XSLoader;
	XSLoader::load('SeqAn', $VERSION);
	
	$self->{'SeqAn'} = new SeqAn();
	
	return $self;
}

# Preloaded methods go here.

sub AUTOLOAD {
	my ($self, @value) = @_;
	my ($name) = $AUTOLOAD=~/([^:]+)$/;
	return $self->{'SeqAn'}->$name(@value);
}

1;
__END__
# Below is stub documentation for your module. You'd better edit it!

=head1 NAME

PEPE::SeqAn - Perl extension for blah blah blah

=head1 SYNOPSIS

  use PEPE::SeqAn;
  blah blah blah

=head1 DESCRIPTION

Stub documentation for PEPE::SeqAn, created by h2xs. It looks like the
author of the extension was negligent enough to leave the stub
unedited.

Blah blah blah.

=head2 EXPORT

None by default.



=head1 SEE ALSO

Mention other useful documentation such as the documentation of
related modules or operating system documentation (such as man pages
in UNIX), or any relevant external documentation such as RFCs or
standards.

If you have a mailing list set up for your module, mention it here.

If you have a web site set up for your module, mention it here.

=head1 AUTHOR

Daniel Guariz Pinheiro, E<lt>dgpinheiro@E<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2018 by Daniel Guariz Pinheiro

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.22.1 or,
at your option, any later version of Perl 5 you may have available.


=cut
