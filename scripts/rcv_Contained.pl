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

use Bio::SeqIO;
use PEPE::SeqAn;

INIT {
    use Log::Log4perl qw/:easy/;
    Log::Log4perl->easy_init($FATAL);
    $LOGGER = Log::Log4perl->get_logger($0);
}

use constant DEFAULT_MIN_ALIGNMENT_LENGTH=>15;
use constant DEFAULT_SCORE_FOR_MATCH=>2;
use constant DEFAULT_SCORE_FOR_MISMATCH=>-4;
use constant DEFAULT_SCORE_FOR_GAP_OPEN=>-6;
use constant DEFAULT_SCORE_FOR_GAP_EXTENSION=>-6;

my ($level, $ctmfile, $outdir, $vectorfile, $suffix, $sample, $minalnlen, $score_for_match, $score_for_mismatch, $score_for_gapopen, $score_for_gapextension);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
	    "i|infile=s"=>\$ctmfile,
	    "o|outdir=s"=>\$outdir,
	    "v|vectorfile=s"=>\$vectorfile,
	    "x|suffix=s"=>\$suffix,
            "s|sample=s"=>\$sample,
            "m|minalnlen=i"=>\$minalnlen,
            "sma|match=i"=>\$score_for_match,
            "smi|mismatch=i"=>\$score_for_mismatch,
            "sgo|gapopen=i"=>\$score_for_gapopen,
            "sge|gapextension=i"=>\$score_for_gapextension
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
            
$score_for_match||=DEFAULT_SCORE_FOR_MATCH;
$score_for_mismatch||=DEFAULT_SCORE_FOR_MISMATCH;
$score_for_gapopen||=DEFAULT_SCORE_FOR_GAP_OPEN;
$score_for_gapextension||=DEFAULT_SCORE_FOR_GAP_EXTENSION;
$minalnlen||=DEFAULT_MIN_ALIGNMENT_LENGTH;

$LOGGER->logdie("Missing input fastq file") if (!defined $ctmfile);
$LOGGER->logdie("Wrong input fastq file ($ctmfile)") if (! -e $ctmfile);

$LOGGER->logdie("Missing output directory") if (! defined $outdir);
$LOGGER->logdie("Wrong output directory ($outdir)") if (! -e $outdir);

$LOGGER->logdie("Missing vector fata file") if (! defined $vectorfile);
$LOGGER->logdie("Wrong vector fasta file ($vectorfile)") if (! -e $vectorfile);

$LOGGER->logdie("Missing sample name") if (! defined $sample);


$suffix||='.fastq';



my $seqin = Bio::SeqIO->new(-file=>$vectorfile, -format=>'FASTA');

my $vectorseq = $seqin->next_seq();

$|=1;

my $bn = basename($ctmfile, '.fastq');

my %outfh;

$outfh{'decontaminated'} = FileHandle->new(">".$outdir.'/'.$sample.'.decontaminated'.$suffix);
autoflush { $outfh{'decontaminated'} } 1;
$LOGGER->logdie($!) if (! defined $outfh{'decontaminated'});


my $seqcountexe = 'echo "$(wc -l '.$ctmfile.' | sed \'s/ .*$//\')/4" | bc -l';
my $seqcount = `$seqcountexe`;
$seqcount+=0;
print "Loading ($seqcount) fragments [$bn] ...\n";

my $progress = Term::ProgressBar->new($seqcount);

open(IN, "<", $ctmfile) or $LOGGER->logdie($!);

my $i = 0;

my @read;
while(<IN>) {
	chomp;
	if ( ($. % 4) == 1 ) {
		$read[$i]->[0] = [ split(/\t/, $_) ];
	} elsif ( ($. % 4) == 2) {
		$read[$i]->[1] = [ split(//, $_) ];
	} elsif ( ($. % 4) == 3) {
	  	next;
	} elsif ( ($. % 4) == 0 ) {
	 	$read[$i]->[2] = [ split(//, $_) ];
		$i++;
		if ( ($i%10) == 0 ) {
			$progress->update($i);
		}
		
	} else {
		$LOGGER->logdie("Error found: $. mod 4 cannot be distinct from 1, 2, 3 or 0 !")
	}
}
close(IN);
$progress->update($i);

print "Processing ($seqcount) fragments [$bn] ...\n";

$progress = Term::ProgressBar->new($seqcount);

my $seqan = new PEPE::SeqAn();

$i=0;
foreach my $ar_read (@read) {
	if ( ($i%10) == 0 ) {
		$progress->update($i);
	}
	$i++;
	#print join("\t", @{ $ar_read->[0] }),"\n",join("", @{ $ar_read->[1] }),"\n",
	#join(" ", @{ $ar_read->[2] }),"\n",
	#join(" ", map { &toq($_) } @{ $ar_read->[2] }),"\n";
	my $seq = join('', @{ $ar_read->[1] });
	my $qual = join('', @{ $ar_read->[2] });
	my $id = $ar_read->[0]->[0];
	# Vector sequence, Read sequene, Local Alignment, Match score, Mismatch score, Extension gap score, Open gap score
	if (length($seq)>(2*$minalnlen)) {
		my ($scorer, $aln1r, $aln2r, $begaln2r, $endaln2r) = split(';', $seqan->Align2Seq($seq, $vectorseq->seq(), 1, $score_for_match, $score_for_mismatch, $score_for_gapextension, $score_for_gapopen));
		if ( (length($aln1r) >= $minalnlen) ) {
	#		print $scorer,"\n";
	#		print $aln1r,"\n";
	#		print $aln2r,"\n";
	#		print $begaln2r,"\n";
	#		print $endaln2r,"\n";
	#		print '>'.length($seq),"\n";
	#		print "$id\n$seq\n";#,substr($seq,$start,($end-$start)),"\n";

			print { $outfh{'decontaminated'} } $id,"\t$aln2r\n".substr($seq, 0, $begaln2r-1),"\n","+\n",substr($qual, 0, $begaln2r-1),"\n";
		}
	#	} else {
	#		print ':::',($scorer/&bestscore($aln1r,1)),"\n";
	#		print ':::',$scorer,"\n";
	#		print ':::',$aln1r,"\n";
	#		print ':::',$aln2r,"\n";
	#		print ':::',$begaln2r,"\n";
	#		print ':::',$endaln2r,"\n";
	#		print ':::>'.length($seq),"\n";
	#	}
	}
	$progress->update($i);
}


# Subroutines

sub Usage {
    my ($msg) = @_;
	
	$score_for_match=DEFAULT_SCORE_FOR_MATCH;
	$score_for_mismatch=DEFAULT_SCORE_FOR_MISMATCH;
	$score_for_gapopen=DEFAULT_SCORE_FOR_GAP_OPEN;
	$score_for_gapextension=DEFAULT_SCORE_FOR_GAP_EXTENSION;
	$minalnlen=DEFAULT_MIN_ALIGNMENT_LENGTH;

    my $USAGE = <<"END_USAGE";
Daniel Guariz Pinheiro (dgpinheiro\@gmail.com)
(c)2018 Universidade Estadual Paulista "Júlio de Mesquita Filho"

Usage

        $0	[-h/--help] [-l/--level <LEVEL>]

Argument(s)

	-h	--help		Help
	-l	--level		Log level [Default: FATAL]
	-i	--infile	Contaminated input fastq file
	-o	--outdir	Output directory
	-v	--vectorfile	Vector fasta file
	-s	--sample	Sample name
	-x	--suffix	Suffix added to output file [Default .fastq]
	-m	--minalnlen	Minimum alignment length to consider a vector contamination [Default: $minalnlen]
	-sma	--match		Score for match [Default: $score_for_match]
	-smi	--mismatch	Score for mismatch [Default: $score_for_mismatch]
	-sgo	--gapopen	Score for gap open [Default: $score_for_gapopen]
	-sge	--gapextension	Score for gap extension [Default: $score_for_gapextension]

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
    print STDERR $USAGE;
    exit(1);
}

sub simpleAlign {
	my ($seq1, $seq2) = @_;
	my @s1 = split(//, $seq1);
	my @s2 = split(//, $seq2);
	my $nmm = 0;
	for (my $i=0; $i<=$#s1; $i++) {
		if ($s1[$i] ne $s2[$i]) {
			$nmm++ if (($s1[$i] ne 'N') && ($s2[$i] ne 'N'));
		}
	}
	return($nmm, $seq1, $seq2);
}

sub toq {
	my ($qchr) = @_;
	return (ord($qchr)-33);
}

sub bestscore {
	my ($seq, $match) = @_;
	
	return(length($seq)*$match);
}
