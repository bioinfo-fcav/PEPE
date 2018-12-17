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

use Bio::SeqIO;

use vars qw/$LOGGER/;
use FileHandle;

use File::Temp qw/ tempfile tempdir /;

INIT {
    use Log::Log4perl qw/:easy/;
    Log::Log4perl->easy_init($FATAL);
    $LOGGER = Log::Log4perl->get_logger($0);
}

my ($level, $infile, $minsize, $maxsize, $id);

#Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
	    "i|infile=s"=>\$infile,
            "n|minsize=i"=>\$minsize,
            "m|maxsize=i"=>\$maxsize,
            "d|id=s"=>\$id
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

$id||='UNKNOWN';

my @data;
$minsize||=0;
$maxsize||=1000000000;


if ($infile) {
	$LOGGER->logdie("Wrong input genomeCoverageBed output file ($infile)") if (! -e $infile);
	open(IN, "<", $infile) or $LOGGER->logdie($!);
	while(<IN>) {
		chomp;
		my @F=split(/\t/, $_);
		push(@data, $F[2]);
	}
	close(IN);
} elsif (! -t STDIN) {
	while(<STDIN>) {
		chomp;
		push(@data, $_);
	}
}

my $proteinsize=scalar(@data);
if ( $proteinsize > 0 ) {

	# Determine mean. Or use Statistics::Descriptive.
	my $sum;
	$sum += $_ for @data;
	my $mean = $sum / @data;

	# Make a pass over the data to find contiguous runs of values
	# that are either less than or greater than the mean. Also
	# keep track of the mins and maxes within those groups.
	my $group = -1;
	my $gt_mean_prev = '';
	my @mins_maxs;
	my $i = -1;

	for my $d (@data){
	    $i ++;
	    my $gt_mean = $d > $mean ? 1 : 0;

	    unless ($gt_mean eq $gt_mean_prev){
		$gt_mean_prev = $gt_mean;
		$group ++;
		$mins_maxs[$group] = $d;
	    }

	    if ($gt_mean){
		$mins_maxs[$group] = $d if $d > $mins_maxs[$group];
	    }
	    else {
		$mins_maxs[$group] = $d if $d < $mins_maxs[$group];
	    }

	    $d = {
		i       => $i,
		val     => $d,
		group   => $group,
		gt_mean => $gt_mean,
	    };
	}

	# A fun picture.
	my @coord;
	my $j = 0;
	my $set;
	for my $d (@data){
		if (defined $set) {
			if (!$d->{gt_mean}) {
				$coord[$j]->[1] = $d->{i}-1;
				$set = undef;
				$j++;
			}
		} else {
			if ($d->{gt_mean}) {
				$coord[$j] = [$d->{i}, $#data];
				$set = 1;
			}
		}
	#    printf
	#        "%6.1f  %2d  %1s  %1d  %3s  %s\n",
	#        $d->{val},
	#        $d->{i},
	#        $d->{gt_mean} ? '+' : '-',
	#        $d->{group},
	#        $d->{val} == $mins_maxs[$d->{group}] ? '==>' : '',
	#        '.' x ($d->{val} / 2)
	#    ;
	}

	foreach my $c (0..$#coord) {
		my $len = $coord[$c]->[1]-$coord[$c]->[0];
		if (($len>=$minsize)&&($len<=$maxsize)) {
			print $id,"\t",join("\t", @{$coord[$c]}[0,1]),"\t",$len,"\t",$proteinsize,"\n";
		}
	}
}

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
	-i	--infile	Input file (genomeCoverageBed output file) or STDIN with coverage
	-n	--minsize	Minimum size [Default: 0]
	-m	--maxsize	Maximum size [Default: 1000000000]
	-d	--id		Identifier [Deafaut: UNKNOWN]

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
    print STDERR $USAGE;
    exit(1);
}
