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
use Term::ProgressBar;
use FileHandle;


INIT {
    use Log::Log4perl qw/:easy/;
    Log::Log4perl->easy_init($FATAL);
    $LOGGER = Log::Log4perl->get_logger($0);
}

my ($level, $infile, $outdir, $suffix, $sample);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
	    "i|infile=s"=>\$infile,
	    "o|outdir=s"=>\$outdir,
	    "s|sample=s"=>\$sample,
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

$LOGGER->logdie("Missing input alignment file") if (!defined $infile);
$LOGGER->logdie("Wrong input alignment file ($infile)") if (! -e $infile);

$LOGGER->logdie("Missing output directory") if (! defined $outdir);
$LOGGER->logdie("Wrong output directory ($outdir)") if (! -e $outdir);

$LOGGER->logdie("Missing sample name") if (! defined $sample);

$suffix||='.txt';

my %outfh;

$outfh{'id'} = FileHandle->new(">".$outdir.'/'.$sample.'.id'.$suffix);
autoflush { $outfh{'id'} } 1;
$LOGGER->logdie($!) if (! defined $outfh{'id'});

$outfh{'notid'} = FileHandle->new(">".$outdir.'/'.$sample.'.notid'.$suffix);
autoflush { $outfh{'notid'} } 1;
$LOGGER->logdie($!) if (! defined $outfh{'notid'});

my $linecountexe = 'echo "$(wc -l '.$infile.' | sed \'s/ .*$//\')"';
my $linecount = `$linecountexe`;
$linecount+=0;

print "Loading ($linecount) alignments [$sample] ...\n";

my $progress = Term::ProgressBar->new($linecount);

open(IN, "<", $infile) or $LOGGER->logdie($!);

my $i = 0;

my %read;
my %bestqframeerr;
my %bestscore;
while(<IN>) {
chomp;
	$i++;
	if ( ($i%10) == 0 ) {
		$progress->update($i);
	}
	
	my ($qseqid, $qlen, $sseqid, $bitscore, $qcovhsp, $qframe, $sstart, $send) = split(/\t/, $_);

	#$sseqid =~s/\.\d+$//;
	my $r = 1;
	if ($qseqid=~/\/(\d+)$/) {
		$r = $1;
		$qseqid=~s/\/(\d+)$//;
	}
	
	if ($r == 1) {
		if ($qframe != 1) {
			$bestqframeerr{ $qseqid } = $sseqid;
		}
	}	

	if (exists $bestqframeerr{ $qseqid }) {
		next;
	}
	
	if (! exists $read{$qseqid}->{$r}) {
		$read{$qseqid}->{$r}->{$sseqid} = [$bitscore, $sstart, $send];
		$bestscore{$qseqid}->{$r} = $bitscore;
	} else {
		if ( $bitscore >= ($bestscore{$qseqid}->{$r}-(0.00*$bestscore{$qseqid}->{$r})) ) {
			$read{$qseqid}->{$r}->{$sseqid} = [$bitscore, $sstart, $send];
		}
	}
}
close(IN);
$progress->update($i);

my $readcount = scalar(keys %read);
print "Processing (".$readcount.") alignments [$sample] ...\n";

$progress = Term::ProgressBar->new($readcount);

$i=0;
foreach my $qseqid (keys %read) {
	$i++;
	if ( ($i%10) == 0 ) {
		$progress->update($i);
	}
	my @subject;
	if (scalar(keys %{ $read{$qseqid} }) == 1) {
		my ($r) = %{ $read{$qseqid} };
		# Discard match if only read 2 had previously selected
		next if ($r == 2);
		
		if (scalar(keys %{ $read{$qseqid}->{$r} }) == 1) {
			my ($sseqid) = keys %{ $read{$qseqid}->{$r} };
			#OLD: M03855:91:000000000-B4GJ9:1:2108:9708:19434	AT3G48990	75	128
			#print { $outfh{'id'} } $qseqid,"\t",$sseqid,"\t",join("\t",@{ $read{$qseqid}->{$r}->{$sseqid}}[1,2]),"\n";
			#BED: AT3G48990.1	75	128	M03855:91:000000000-B4GJ9:1:2108:9708:19434	0	+
			print { $outfh{'id'} } $sseqid,"\t",join("\t",@{ $read{$qseqid}->{$r}->{$sseqid}}[1,2]),"\t",$qseqid,"\t",0,"\t",'+',"\n";
		} else {
			foreach my $sseqid (keys %{ $read{$qseqid}->{$r} }) {
				#print { $outfh{'notid'} } $qseqid,"\t",$sseqid,"\t",join("\t",@{ $read{$qseqid}->{$r}->{$sseqid}}[1,2]),"\n";
				print { $outfh{'notid'} } $sseqid,"\t",join("\t",@{ $read{$qseqid}->{$r}->{$sseqid}}[1,2]),"\t",$qseqid,"\t",0,"\t",'+',"\n";
			}
		}
	} else {
		my @match;
		foreach my $sseqid1 (keys  %{ $read{$qseqid}->{1} }) {
			foreach my $sseqid2 (keys %{ $read{$qseqid}->{2} }) {
				if ($sseqid1 eq $sseqid2) {
					push(@match, $sseqid1);
				}
			}
		}
		if (scalar(@match)==1) {
			my $sseqid = $match[0];
			my ($min_start, $max_end) = (1000000,0);
			foreach my $r (1,2) {
				if ($read{$qseqid}->{$r}->{$sseqid}->[1] < $min_start) {
					$min_start = $read{$qseqid}->{$r}->{$sseqid}->[1];
				}
				if ($read{$qseqid}->{$r}->{$sseqid}->[2] > $max_end) {
					$max_end = $read{$qseqid}->{$r}->{$sseqid}->[2];
				}
			}
			#print { $outfh{'id'} } $qseqid,"\t",$sseqid,"\t",join("\t",$min_start,$max_end),"\n";
			print { $outfh{'id'} } $sseqid,"\t",join("\t",$min_start,$max_end),"\t",$qseqid,"\t",0,"\t",'+',"\n";
		} else {
			foreach my $sseqid ( @match ) {
				my ($min_start, $max_end) = (1000000,0);
				foreach my $r (1,2) {
					if ($read{$qseqid}->{$r}->{$sseqid}->[1] < $min_start) {
						$min_start = $read{$qseqid}->{$r}->{$sseqid}->[1];
					}
					if ($read{$qseqid}->{$r}->{$sseqid}->[2] > $max_end) {
						$max_end = $read{$qseqid}->{$r}->{$sseqid}->[2];
					}
				}
				
				#print { $outfh{'notid'} } $qseqid,"\t",$sseqid,"\t",join("\t",$min_start,$max_end),"\n";
				print { $outfh{'notid'} } $sseqid,"\t",join("\t",$min_start,$max_end),"\t",$qseqid,"\t",0,"\t",'+',"\n";
			}
		}
	}
}
$progress->update($i);

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
	-i	--infile	Input alignment file
	-o	--outdir	Output directory
	-s	--sample	Sample name
	-x	--suffix	Suffix added to output file [Default: .txt]

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
    print STDERR $USAGE;
    exit(1);
}
