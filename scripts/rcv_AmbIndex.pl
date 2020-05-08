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

use constant DEFAULT_MODEL_ADAPTER_SEQUENCE_3PRIME=>'CCTCTAAACGGGTCTTGAGGGGTT';



use PEPE::SeqAn;

INIT {
    use Log::Log4perl qw/:easy/;
    Log::Log4perl->easy_init($FATAL);
    $LOGGER = Log::Log4perl->get_logger($0);
}

my ($level, $ambfile, $outdir, $barcodesfile, $suffix, $model_adapter_sequence_3prime);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
	    "a|ambfile=s"=>\$ambfile,
	    "o|outdir=s"=>\$outdir,
	    "b|barcodes=s"=>\$barcodesfile,
            "m3p|model3pseq=s"=>\$model_adapter_sequence_3prime,
	    "s|suffix=s"=>\$suffix
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

$model_adapter_sequence_3prime||=DEFAULT_MODEL_ADAPTER_SEQUENCE_3PRIME;

$LOGGER->logdie("Missing input fastq file") if (!defined $ambfile);
$LOGGER->logdie("Wrong input fastq file ($ambfile)") if (! -e $ambfile);

$LOGGER->logdie("Missing output directory") if (! defined $outdir);
$LOGGER->logdie("Wrong output directory ($outdir)") if (! -e $outdir);

$LOGGER->logdie("Missing barcodes file") if (! defined $barcodesfile);
$LOGGER->logdie("Wrong barcodes file ($barcodesfile)") if (! -e $barcodesfile);

my %sample;

$suffix||='.fastq';

open(BC, "<", $barcodesfile) or $LOGGER->logdie($!);

while(<BC>) {
	chomp;
	my ($bcseq, $sample_name) = split(/\t/, $_);
	$sample_name=~s/^\s*//;
	$sample_name=~s/\s*$//;

	$sample{$bcseq} = {	'name'=>$sample_name,
				'afh'=> FileHandle->new(">".$outdir.'/'.$sample_name.'.recovered'.$suffix)
			  };
	autoflush { $sample{$bcseq}->{'afh'} } 1;
	$LOGGER->logdie($!) if (! defined $sample{$bcseq}->{'afh'});
}

my @bc = keys %sample;

close(BC);

$|=1;

my $bn = basename($ambfile, '.fastq');

my $seqcountexe = 'echo "$(wc -l '.$ambfile.' | sed \'s/ .*$//\')/4" | bc -l';
my $seqcount = `$seqcountexe`;
$seqcount+=0;
print "Loading ($seqcount) fragments [$bn] ...\n";

my $progress = Term::ProgressBar->new($seqcount);

open(IN, "<", $ambfile) or $LOGGER->logdie($!);

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
	my $tag = $ar_read->[0]->[1];
	$seq =~ /($tag)/g;
	my $p = pos($seq);
	if (defined $p) {
		my $start = $p-length($tag);
		my $k = 0;
		my %phred;
		my $min = 100;
		my @index = ();
		for (my $j=($start+2); $j<=($start+4); $j++) { 
			my $q = &toq($ar_read->[2]->[$j]);
			if ($q < $min) {
				$min = $q;
			}
			push(@index, $ar_read->[1]->[$j]);
			$phred{ $q }->{ $k } = undef;
			$k++;
			#print "\t",$j,"\t",$ar_read->[1]->[$j],"\t",&toq($ar_read->[2]->[$j]),"\n";
		}
		#print $tag,"\n",substr($seq, $start+2, 3),"\n\n";
		if (scalar(keys %{$phred{$min}}) == 1) {
			my ($x) = keys %{ $phred{$min} };
			$index[$x] = 'N';

			my $idx = join('', @index);
			my @sorder;
			my %mm;
			foreach my $kseq (@bc) {
				(@{ $mm{$kseq} }) = &simpleAlign($idx, $kseq);
			}
			@sorder = sort { $mm{$a}->[0] <=> $mm{$b}->[0] } keys %mm;
			my $s = $sorder[0];
			if ((!defined $sorder[1])||($mm{$s}->[0] < $mm{ $sorder[1] }->[0])) {
				#print "************** NEW SAMPLE **************\t",$sample{$s}->{'name'},"\n";
				my $start = $p;
				my $end = length($seq);
				my ($scorer, $aln1r, $aln2r) = split(';', $seqan->Align2Seq($model_adapter_sequence_3prime, $seq));
				while ($aln2r=~/^-/g) {
					$aln2r=~s/^-//;
					$aln1r=~s/^.//;
				}
				while ($aln2r=~/-$/g) {
					$aln2r=~s/-$//;
					$aln1r=~s/.$//;
				}
				$aln1r=~s/-/\./g;

				if ($seq=~/($aln1r)/g) {
					$end = pos($seq)-length($model_adapter_sequence_3prime);
				}
				print { $sample{ $s }->{'afh'} } $id,"\t$tag\t$ar_read->[0]->[2]\t$s\n",substr($seq,$start,($end-$start)),"\n","+\n",substr($qual,$start,($end-$start)),"\n";
			}
		}
	} else {
		$LOGGER->logdie("Not found match ($id) tag ($tag) position on $seq");
	}
}
$progress->update($i);

foreach my $s (keys %sample) {
	$sample{$s}->{'afh'}->close();
}

# Subroutines

sub Usage {
    my ($msg) = @_;

	$model_adapter_sequence_3prime=DEFAULT_MODEL_ADAPTER_SEQUENCE_3PRIME;

    my $USAGE = <<"END_USAGE";
Daniel Guariz Pinheiro (dgpinheiro\@gmail.com)
(c)2018 Universidade Estadual Paulista "Júlio de Mesquita Filho"

Usage

        $0	[-h/--help] [-l/--level <LEVEL>]

Argument(s)

	-h	--help		Help
	-l	--level		Log level [Default: FATAL]
	-a	--ambfile	Ambiguous fastq file
	-o	--outdir	Output directory
	-b	--barcodes	Barcodes file
	-m3p	--model3pseq	Model adapter sequence 3' [Default: $model_adapter_sequence_3prime]
	-s	--suffix	Suffix added to output file [Default .fastq]

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
