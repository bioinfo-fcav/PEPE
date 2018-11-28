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

my ($level, %defines, %groups);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
	    "d|define=s"=>\%defines,
	    "g|group=s"=>\%groups
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

my %gene;
my %sum;
my @sample;
foreach my $g (keys %groups) {
	
	#print $g,"\n";
	
	$groups{$g}=~s/\s//g;
	
	foreach my $k ( split(',', $groups{$g}) ) {
		if (! exists $defines{$k}) {
			$LOGGER->logdie("Missing definition for element ($k) of group $g");
		} else {
			if (! -e $defines{$k}) {
				$LOGGER->logdie("Wrong definition file ($defines{$k}) for element ($k) of group ($g)");
			}
			
			push(@sample, $k);
			
			#print "\t",$k,"\t",$defines{$k},"\n";
			
			open(IN, "<", $defines{$k}) or $LOGGER->logdie($!);

			while(<IN>) {
				chomp;
				my ($geneid, $count) = split(/\t/, $_);
				$gene{$geneid}->{$k} = $count;
				$sum{$geneid}+=$count;
			}
			close(IN);
		}
	}
}


print join("\t", 'ID', 'SUM', @sample),"\n";

foreach my $geneid (sort { $sum{ $b } <=> $sum{ $a } } keys %sum) {
	my @v;
	foreach my $k (@sample) {
		push(@v, $gene{$geneid}->{$k}||0);	
	}

	print $geneid,"\t",$sum{$geneid},"\t",join("\t", @v),"\n";
}

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
	-d	--define	Defines input counts file
	-g	--group		Groups of samples

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
    print STDERR $USAGE;
    exit(1);
}

