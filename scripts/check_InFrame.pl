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
use Bio::Seq;
use Bio::Tools::CodonTable;

use vars qw/$LOGGER/;
use Term::ProgressBar;
use FileHandle;

use SeqAn;

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

$LOGGER->logdie("Missing input fastq file") if (!defined $infile);
$LOGGER->logdie("Wrong input fastq file ($infile)") if (! -e $infile);

$LOGGER->logdie("Missing output directory") if (! defined $outdir);
$LOGGER->logdie("Wrong output directory ($outdir)") if (! -e $outdir);

$LOGGER->logdie("Missing sample name") if (! defined $sample);

$suffix||='.fastq';

my %outfh;

$outfh{'inframe'} = FileHandle->new(">".$outdir.'/'.$sample.'.inframe'.$suffix);
autoflush { $outfh{'inframe'} } 1;
$LOGGER->logdie($!) if (! defined $outfh{'inframe'});

$outfh{'notinframe'} = FileHandle->new(">".$outdir.'/'.$sample.'.notinframe'.$suffix);
autoflush { $outfh{'notinframe'} } 1;
$LOGGER->logdie($!) if (! defined $outfh{'notinframe'});

$outfh{'wrongdna'} = FileHandle->new(">".$outdir.'/'.$sample.'.wrongdna'.$suffix);
autoflush { $outfh{'wrongdna'} } 1;
$LOGGER->logdie($!) if (! defined $outfh{'wrongdna'});


my $bn = basename($infile, '.fastq');

my $seqcountexe = 'echo "$(wc -l '.$infile.' | sed \'s/ .*$//\')/4" | bc -l';
my $seqcount = `$seqcountexe`;
$seqcount+=0;
print "Loading ($seqcount) fragments [$bn] ...\n";

my $progress = Term::ProgressBar->new($seqcount);

open(IN, "<", $infile) or $LOGGER->logdie($!);

my $i = 0;

my @read;
while(<IN>) {
	chomp;
	if ( ($. % 4) == 1 ) {
		$_=~s/^@//;
		$read[$i]->[0] = [ split(/\t/, $_) ];
		$read[$i]->[0]->[0]=~s/ .*//;
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

my $seqan = new SeqAn();

$i=0;
foreach my $ar_read (@read) {
	if ( ($i%10) == 0 ) {
		$progress->update($i);
	}
	$i++;
	my $fullseq = join('', @{ $ar_read->[1] });
	my $fullqual = join('', @{ $ar_read->[2] });
	my $id = $ar_read->[0]->[0];
	my $tag = $ar_read->[0]->[1];
	
	my @dual_seq  = split(/NNNNN/, $fullseq);
	my @dual_qual = split(/!!!!!/, $fullqual);
	
	my $c=undef;
	if (scalar(@dual_seq)>1) {
		$c = 1;
	}
	for (my $j=0; $j<=$#dual_seq; $j++) {
		my $seq = $dual_seq[$j];
		my $qual = $dual_qual[$j];
		
		my $seqo = Bio::Seq->new(-display_id => $id.(($c) ? '/'.$c++ : ''),
					 -seq => $seq,
					 -alphabet=>'dna'
					);

		next unless ($seqo->length()>12);

		my @aa;
		foreach my $frame (0,1,2) {
			my $prot_obj = $seqo->translate(-frame=>$frame);
			my $aaseq = $prot_obj->seq();
			if ($aaseq=~/(MLGDPNS?)/g) {
				my $start = pos($aaseq);
				my $end = length($aaseq);
				if ($aaseq=~/(KLAAA)/g) {
					$end = pos($aaseq)-5;
				}
				#print "#",$aaseq,"\n";
				#print ">",substr($aaseq,$start,$end-$start),"\n";
				push(@aa, [$aaseq,$start,$end,$end-$start]);
				#print { $outfh{ 'inframe' } } '@',$seqo->display_id(),"\t*** NOT FOUND MLGDPN! ***\n",$seq,"\n","+\n",$qual,"\n";
			}
		}
		if (scalar(@aa)) {
			my ($aaseq) = @{ [ sort { $b->[3] <=> $a->[3] } @aa ]->[0] };
			#print ">>>>",$aaseq,"\n";
			my ($scorer, $aln1r, $aln2r) = split(';', $seqan->Align2Seq('ATGCTCGGGGATCCGAATTCN',$seq));
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
				my $start = pos($seq);
				if ($start) {
					my $end = length($seq);
					if ($seq=~/(AAGCTTGCGGCCGCA)/g) {
						$end = pos($seq)-15;
					} else {
						#AAGCTTGCGGCCGCA
	#					($scorer, $aln1r, $aln2r) = split(';', $seqan->Align2Seq('AAGCTTGCGGCCGCACTCGAG',$seq));
	#					while ($aln2r=~/^-/g) {
	#						$aln2r=~s/^-//;
	#						$aln1r=~s/^.//;
	#					}
	#					while ($aln2r=~/-$/g) {
	#						$aln2r=~s/-$//;
	#						$aln1r=~s/.$//;
	#					}
	#					$aln1r=~s/-/\./g;
	#					if ($seq=~/($aln1r)/g) {
	#						print ">>>>>$aln1r\n";
	#					}
					}
					print { $outfh{ 'inframe' } } '@',$seqo->display_id(),"\t$aaseq\n",substr($seq,$start,$end-$start),"\n","+\n",substr($qual,$start,$end-$start),"\n";
				}
			} else {
				print { $outfh{ 'wrongdna' } } '@',$seqo->display_id(),"\n",$seq,"\n","+\n",$qual,"\n";
			}
		} else {
			my $start = 0;
			my $end = length($seq);
			my @aa;
			foreach my $frame (0,1,2) {
				my $prot_obj = $seqo->translate(-frame=>$frame);
				my $aaseq = $prot_obj->seq();
				if ( $aaseq=~/(KLAAALE\*)/g) {
					$end = pos($aaseq)-8;
					if ($seq=~/AAGCTTGCGGCCGCA/) {
						push(@aa, [$aaseq,$start,$end,$end-$start]);
					}
				} 
			}
			if (scalar(@aa)) {
				my ($aaseq) = @{ [ sort { $b->[3] <=> $a->[3] } @aa ]->[0] };
				$start = 0;
				if ($seq=~/(AAGCTTGCGGCCGCA)/g){
					$end = pos($seq)-15;
					print { $outfh{ 'inframe' } } '@',$seqo->display_id(),"\t$aaseq\n",substr($seq,$start,$end-$start),"\n","+\n",substr($qual,$start,$end-$start),"\n";
				} else {
					print { $outfh{ 'wrongdna' } } '@',$seqo->display_id(),"\n",$seq,"\n","+\n",$qual,"\n";
				}
			}
			else {
				print { $outfh{ 'notinframe' } } '@',$seqo->display_id(),"\t*** NOT FOUND MLGDPN OR KLAAALE! ***\n",$seq,"\n","+\n",$qual,"\n";
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
	-i	--infile	Input fastq file
	-o	--outdir	Output directory
	-s	--sample	Sample name
	-x	--suffix	Suffix added to output file [Default .fastq]

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
    print STDERR $USAGE;
    exit(1);
}

