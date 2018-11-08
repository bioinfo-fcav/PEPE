#!/usr/bin/env perl
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

Copyright (c) 2018 Universidade Estadual Paulista "Júlio de Mesquita Filho"

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

my ($level, $indir, $sample, $suffix, $outfile);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
	    "i|indir=s"=>\$indir,
	    "s|sample=s"=>\$sample,
	    "o|outfile=s"=>\$outfile,
	    "x|suffix=s"=>\$suffix
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

$LOGGER->logdie("Missing input directory") if (!defined $indir);
$LOGGER->logdie("Wrong input directory ($indir)") if (! -e $indir);

$LOGGER->logdie("Missing sample name") if (! defined $sample);

$LOGGER->logdie("Missing suffix") if (! defined $suffix);

$LOGGER->logdie("Missing output file") if (!defined $outfile);


my %outfh;

$outfh{'counts'} = FileHandle->new(">".$outfile);
autoflush { $outfh{'counts'} } 1;
$LOGGER->logdie($!) if (! defined $outfh{'counts'});

my @gene;
foreach my $idfile (glob("$indir/$sample$suffix")) {
	my $bnfile = basename($idfile);
	$LOGGER->info("Processing [$sample] $bnfile ...");
	open(IN, "<", $idfile) or $LOGGER->logdie($!);
	while(<IN>) {
		chomp;
		my @F = split(/\t/, $_);
		push(@gene, $F[1]);
	}
	close(IN);
}

my %count;
foreach my $g (sort {$a cmp $b} @gene) {
	$count{$g} = 0 unless (exists $count{$g});
	$count{$g}++;
}

foreach my $g (sort { $count{$b} <=> $count{$a} } keys %count) {
	print { $outfh{'counts'} } $g,"\t",$count{$g},"\n";
}

foreach my $s (keys %outfh) {
	$outfh{$s}->close();
}

# Subroutines

sub Usage {
    my ($msg) = @_;
    my $USAGE = <<"END_USAGE";
Daniel Guariz Pinheiro (dgpinheiro\@gmail.com)
(c)2018 Universidade Estadual Paulista "Júlio de Mesquita Filho"

Usage

        $0	[-h/--help] [-l/--level <LEVEL>]

Argument(s)

	-h	--help		Help
	-l	--level		Log level [Default: FATAL]
	-i	--indir		Input directory
	-s	--sample	Sample name
	-x	--suffix	Suffix
	-o	--outfile	Output file

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
    print STDERR $USAGE;
    exit(1);
}

