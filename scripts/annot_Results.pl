#!/usr/bin/perl
#
#              INGLÊS/ENGLISH
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  http://www.gnu.org/copyleft/gpl.html
#
#
#             PORTUGUÊS/PORTUGUESE
#  Este programa é distribuído na expectativa de ser útil aos seus
#  usuários, porém NÃO TEM NENHUMA GARANTIA, EXPLÍCITAS OU IMPLÍCITAS,
#  COMERCIAIS OU DE ATENDIMENTO A UMA DETERMINADA FINALIDADE.  Consulte
#  a Licença Pública Geral GNU para maiores detalhes.
#  http://www.gnu.org/copyleft/gpl.html
#
#  Copyright (C) 2018  Universidade Estadual Paulista "Júlio de Mesquita Filho"
#
#  Universidade Estadual Paulista "Júlio de Mesquita Filho"
#
#  Daniel Guariz Pinheiro
#  dgpinheiro@gmail.com
#  http://zulu.fmrp.usp.br/bioinfo
#
# $Id$

=head1 NAME

=head1 SYNOPSIS

=head1 ABSTRACT

=head1 DESCRIPTION
    
    Arguments:

        -h/--help   Help
        -l/--level  Log level [Default: FATAL] 
            OFF
            FATAL
            ERROR
            WARN
            INFO
            DEBUG
            TRACE
            ALL

=head1 AUTHOR

Daniel Guariz Pinheiro E<lt>dgpinheiro@gmail.comE<gt>

Copyright (c) 2012 Universidade de São Paulo

=head1 LICENSE

GNU General Public License

http://www.gnu.org/copyleft/gpl.html


=cut

use strict;
use warnings;
use Getopt::Long;

use File::Basename;

use vars qw/$LOGGER/;
use FileHandle;


INIT {
    use Log::Log4perl qw/:easy/;
    Log::Log4perl->easy_init($FATAL);
    $LOGGER = Log::Log4perl->get_logger($0);
}

my ($level, $infile, $genefile, $deffile);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
	    "i|infile=s"=>\$infile,
	    "g|genefile=s"=>\$genefile,
	    "d|deffile=s"=>\$deffile
    ) or &Usage();


if ($level) {
    my %LEVEL = (   
    'OFF'   =>$OFF,
    'FATAL' =>$FATAL,
    'ERROR' =>$ERROR,
    'WARN'  =>$WARN,
    'INFO'  =>$INFO,
    'DEBUG' =>$DEBUG,
    'TRACE' =>$TRACE,
    'ALL'   =>$ALL);
    $LOGGER->logdie("Wrong log level ($level). Choose one of: ".join(', ', keys %LEVEL)) unless (exists $LEVEL{$level});
    Log::Log4perl->easy_init($LEVEL{$level});
}

$LOGGER->logdie("Missing gene name file") if (!defined $genefile);
$LOGGER->logdie("Wrong gene name file ($genefile)") if (! -e $genefile);

my %symbol;
open(GENE, "<", $genefile) or $LOGGER->logdie($!);
while(<GENE>) {
	chomp;
	my ($gene, $name) = split(/\t/, $_);
	$gene=~s/\.\d+//;
	$symbol{ $gene } = $name;
}
close(GENE);

$LOGGER->logdie("Missing defline file") if (!defined $deffile);
$LOGGER->logdie("Wrong defline file ($deffile)") if (! -e $deffile);

my %definition;
open(DEF, "<", $deffile) or $LOGGER->logdie($!);
while(<DEF>) {
	chomp;
	my ($gene, undef, $defline) = split(/\t/, $_);
	$gene=~s/\.\d+$//;
	$defline=~s/defLine\s+//;
	$definition{ $gene } = $defline;
}
close(DEF);

$LOGGER->logdie("Missing input counts file") if (!defined $infile);
$LOGGER->logdie("Wrong input counts file ($infile)") if (! -e $infile);

open(IN, "<", $infile) or $LOGGER->logdie($!);
my $first=<IN>;
chomp($first);
my @header = split(/\t/, $first);
my @sample = @header[1..$#header];
my %sample;
@sample{ @sample } = ();

print join("\t",$header[0], 'Symbol', 'Description', @sample),"\n";

while(<IN>) {
	chomp;
	my %line;
	@line{ @header } = split(/\t/, $_);
	print join("\t", $line{ $header[0] }, $symbol{$line{$header[0]}}||'', $definition{$line{$header[0]}}||'', @line{ @sample }),"\n";
}
close(IN);

my %data;


# Subroutines

sub Usage {
    my ($msg) = @_;
    my $USAGE = <<"END_USAGE";
Daniel Guariz Pinheiro (dgpinheiro\@gmail.com)
(c)2012 Universidade de São Paulo

Usage

        $0	[-h/--help] [-l/--level <LEVEL>]

Argument(s)

	-h	--help		Help
	-l	--level		Log level [Default: FATAL]
	-i	--infile	Input sample counts file
	-g	--genefile	Phytozome gene name file
	-d	--deffile	Phytozome defline file

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
    print STDERR $USAGE;
    exit(1);
}

