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

my ($level, $indir, $samplegroup, $suffix, $outfile, $protfile, $protlist, $minpeaksize, $maxpeaksize);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
	    "i|indir=s"=>\$indir,
	    "g|group=s"=>\$samplegroup,
	    "o|outfile=s"=>\$outfile,
	    "x|suffix=s"=>\$suffix,
	    "p|protfile=s"=>\$protfile,
            "t|protlist=s"=>\$protlist,
            "n|minpeaksize=i"=>\$minpeaksize,
            "m|maxpeaksize=i"=>\$maxpeaksize
    ) or &Usage();


$minpeaksize||=0;
$maxpeaksize||=1000000000;

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

$LOGGER->logdie("Missing input directory with protein x read mapping bed files") if (!defined $indir);
$LOGGER->logdie("Wrong input directory with protein x read mapping bed files ($indir)") if (! -e $indir);

$LOGGER->logdie("Missing input protein sequence fasta file") if (!defined $protfile);
$LOGGER->logdie("Wrong input protein sequence fasta file ($protfile)") if (! -e $protfile);

my %selprot;

if ($protlist) {
	$LOGGER->logdie("Wrong input protein list ($protlist)") if (! -e $protlist);
	open(IN, "<", $protlist) or $LOGGER->logdie($!);
	while(<IN>) {
		chomp;
		my @F=split(/\t/, $_);
		$selprot{$F[0]} = undef;
	}
	close(IN);
} elsif (! -t STDIN) {
	while(<STDIN>) {
		chomp;
		$selprot{$_} = undef;
	}
}

$LOGGER->logdie("Missing sample group name list [ prefix(es) ]") if (! defined $samplegroup);

$LOGGER->logdie("Missing suffix *including .bed") if (! defined $suffix);

$LOGGER->logdie("Missing output file") if (!defined $outfile);


my %protein;

my $seqin = Bio::SeqIO->new(-file=>$protfile, -format=>'FASTA');

while(my $seq = $seqin->next_seq()) {
	my $seqid=$seq->display_id();
	my $seqidwv=$seqid;
	$seqidwv=~s/\.\d+$//;
	if ( (scalar(keys %selprot) == 0) || (exists $selprot{ $seqid }) || ( exists $selprot{ $seqidwv } ) ) {
		
		my $sseq=$seq->seq();
		$sseq=~s/\*//g;
		$seq->seq($sseq);
		
		$protein{ $seq->display_id() } = $seq->length();
		#print $seq->display_id(),"\t",$seq->length(),"\n";
	}
}

my %outfh;

$outfh{'cov'} = FileHandle->new(">".$outfile);
autoflush { $outfh{'cov'} } 1;
$LOGGER->logdie($!) if (! defined $outfh{'cov'});

my $tmpdir = tempdir('pepeXXXX', CLEANUP=>1, DIR=>'./');

my %data;
$samplegroup=~s/\s+//g;

foreach my $sg (split(/,/, $samplegroup)) {
	foreach my $idfile (glob("$indir/$sg$suffix")) {
		my $bnfile = basename($idfile);
		$LOGGER->info("Loading [$sg] $bnfile ...");
		#my ($bnfileid) = $bnfile=~/^([^\.]+)/;
		#print $bnfileid,"\n";
		open(IN, "<", $idfile) or $LOGGER->logdie($!);
		while(<IN>) {
			chomp;
			my @F=split(/\t/, $_);
			if (exists $protein{ $F[0] }) {
				push(@{ $data{$F[0]} }, \@F);
			}
		}
		close(IN);
	}
}

foreach my $protid (keys %data) {
	#print $protid,"\n";
	
	my ($fh, $filename) = tempfile( DIR => $tmpdir, CLEANUP=>0 );
	
	foreach my $armap (sort { $a->[1] <=> $b->[1] } @{ $data{$protid} }) { 
		print { $fh } join("\t",@{ $armap }),"\n";
	}
	$fh->close();
	
	my ($gfh, $gfilename) = tempfile( DIR => $tmpdir, CLEANUP=>0 );
	print { $gfh } $protid."\t".$protein{$protid};
	$gfh->close();

	my $gcbedcmd = 'genomeCoverageBed -i '.$filename.' -g '.$gfilename.' -d > '."$tmpdir/$protid.cov";
	system($gcbedcmd) == 0 or $LOGGER->logdie("System call ($gcbedcmd) failed");
	
	print { $outfh{'cov'} } `$Bin/peaksfinder.pl -i $tmpdir/$protid.cov -d $protid -n $minpeaksize -m $maxpeaksize`;
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
	-i	--indir		Input directory with .bed files
	-s	--group		Sample group name (Prefix)
	-x	--suffix	Suffix (including .bed)
	-o	--outfile	Output file
	-p	--protfile	Protein fasta (.fa) file
	-t	--protlist	Protein list file
	-n	--minpeaksize	Minimum peak size [Default: 0]
	-m	--maxpeaksize	Maximum peak size [Default: 1000000000]

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
    print STDERR $USAGE;
    exit(1);
}

