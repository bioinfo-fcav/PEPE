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

use FindBin qw/$Bin/;

INIT {
    use Log::Log4perl qw/:easy/;
    Log::Log4perl->easy_init($FATAL);
    $LOGGER = Log::Log4perl->get_logger($0);
}

my ($level, $infile, $outdir, $protfile);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
	    "i|infile=s"=>\$infile,
	    "o|outdir=s"=>\$outdir,
	    "p|protfile=s"=>\$protfile
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

$LOGGER->logdie("Missing input protein sequence fasta file") if (!defined $protfile);
$LOGGER->logdie("Wrong input protein sequence fasta file ($protfile)") if (! -e $protfile);

my %data;
if ($infile) {
	$LOGGER->logdie("Wrong input file with peaks ($infile)") if (! -e $infile);
	open(IN, "<", $infile) or $LOGGER->logdie($!);
	while(<IN>) {
		chomp;
		my @F=split(/\t/, $_);
		push(@{ $data{$F[0]} }, \@F);
	}
	close(IN);
} elsif (! -t STDIN) {
	while(<STDIN>) {
		chomp;
		my @F=split(/\t/, $_);
		push(@{ $data{$F[0]} }, \@F);
	}
}

$LOGGER->logdie("Missing output dir") if (!defined $outdir);
$LOGGER->logdie("Wrong output dir ($outdir)") if (! -d $outdir);

my %protein;

my $seqin = Bio::SeqIO->new(-file=>$protfile, -format=>'FASTA');

while(my $seq = $seqin->next_seq()) {
	if ( exists $data{ $seq->display_id() } ) {
		my $sseq = $seq->seq();
		$sseq=~s/\*//g;
		$seq->seq($sseq);
		$protein{ $seq->display_id() } = $seq;
		#print $seq->display_id(),"\t",$seq->length(),"\n";
	}
}

my %outfh;

my $tmpdir = tempdir('pepeXXXX', CLEANUP=>1, DIR=>'./');

foreach my $protid (keys %data) {
	#print $protid,"\n";
	
	my $seqout = Bio::SeqIO->new(-file=>'>'.$tmpdir.'/'.$protid.'.fa', -format=>'FASTA');
	$seqout->write_seq( $protein{ $protid } );
	my $coordsparam='';
	
	foreach my $ar_data (@{ $data{$protid} }) {
		$coordsparam.=' -c '.$ar_data->[1].' '.$ar_data->[2];
	}	
	
	$coordsparam||='';
	
	my $gcbedcmd = "$Bin/hydroplot.py -i $tmpdir/$protid.fa  -s hw $coordsparam -w 9 -o $outdir/$protid.png";

	system($gcbedcmd) == 0 or $LOGGER->logdie("System call ($gcbedcmd) failed");
	
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
	-i	--infile	Input file with peaks
	-o	--outdir	Output directory
	-p	--protfile	Protein fasta (.fa) file

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
    print STDERR $USAGE;
    exit(1);
}

