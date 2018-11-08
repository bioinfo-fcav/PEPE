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

my ($level, $infile, %groups, $min, $outgroup, $minsamples, $mingroups);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
	    "i|infile=s"=>\$infile,
	    "g|group=s"=>\%groups,
	    "m|min=i"=>\$min,
	    "o|outgroup=s"=>\$outgroup,
	    "s|minsamples=i"=>\$minsamples,
	    "u|mingroups=i"=>\$mingroups
    ) or &Usage();

$min||=1;
$minsamples||=2;
$mingroups||=2;

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


$LOGGER->logdie("Missing input counts file") if (!defined $infile);
$LOGGER->logdie("Wrong input counts file ($infile)") if (! -e $infile);

my %data;

open(IN, "<", $infile) or $LOGGER->logdie($!);
my $first=<IN>;
chomp($first);
my @header = split(/\t/, $first);
my @sample = @header[2..$#header];
my %sample;
@sample{ @sample } = ();

foreach my $g (keys %groups) {
	$groups{$g}=~s/\s//g;
	my @smp = split(',', $groups{$g});
	$groups{$g} = [];
	foreach my $k ( @smp ) {
		if (! exists $sample{ $k }) {
			$LOGGER->logdie("Missing sample ($k) on input file");
		} else {
			push(@{ $groups{$g} }, $k);
		}
	}
}

$LOGGER->logdie("Missing outgroup") if (!defined $outgroup);
$LOGGER->logdie("Wrong outgroup ($outgroup). Not defined in any -g/--group") if (! exists $groups{ $outgroup });

while(<IN>) {
	chomp;
	my %line;
	@line{ @header } = split(/\t/, $_);
	foreach my $g (keys %groups) {
		$data{ $line{'ID'} }->{ $g } = 0;
		foreach my $s ( @{ $groups{$g} } ) {
			if ($line{ $s }>=$min) {
				$data{ $line{'ID'} }->{ $g }++;
			}
		}
#		print "\t",$data{ $line{'ID'} }->{ $g },"\n";
	}
}

close(IN);

my @main;
foreach my $g (keys %groups) {
	if ($g ne $outgroup) {
		push(@main, $g);
	}
}

print join("\t",'ID', @main, $outgroup),"\n";
foreach my $id (keys %data) {
	my $mainsc = 0;
	my $maingc = 0;
	foreach my $g (@main) {
		$mainsc+=$data{$id}->{$g};
		$maingc++ if ($data{$id}->{$g} > 0);
	}
	if (	( $mainsc >= $minsamples )&&
		( (! exists $data{$id}->{$outgroup} )||($data{$id}->{$outgroup}==0) )&&
		( $maingc >= $mingroups )
	   ) {
		print join("\t", $id,@{$data{$id}}{@main,$outgroup}),"\n";
	}
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
	-i	--infile	Input sample counts file
	-g	--group		Groups of sample
	-m	--min		Minimum evidence counts to be considered present in a sample [Deafult: 1]
	-s	--minsamples	Minimum samples in the main group (the groups except outgroup) [Default: 2]
	-u	--mingroups	Minimum groups with evidence accounted (* see -m) [Default:2]
	-o	--outgroup	Outgroup

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
    print STDERR $USAGE;
    exit(1);
}

