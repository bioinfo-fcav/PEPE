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
use SeqAn;
use FileHandle;

INIT {
    use Log::Log4perl qw/:easy/;
    Log::Log4perl->easy_init($FATAL);
    $LOGGER = Log::Log4perl->get_logger($0);
}

my ($level, $r1file, $r2file, $outdir, $barcodesfile, $maxmismatches);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
	    "r1|r1file=s"=>\$r1file,
	    "r2|r2file=s"=>\$r2file,
	    "o|outdir=s"=>\$outdir,
	    "b|barcodes=s"=>\$barcodesfile,
	    "mm|maxmismatches=s"=>\$maxmismatches
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

$LOGGER->logdie("Missing R1 and/or R2 input file(s)") if (!(defined $r1file)&&(defined $r2file));
$LOGGER->logdie("Wrong R1 input file ($r1file)") if (! -e $r1file);
$LOGGER->logdie("Wrong R2 input file ($r2file)") if (! -e $r2file);

$LOGGER->logdie("Missing output directory") if (! defined $outdir);
$LOGGER->logdie("Wrong output directory ($outdir)") if (! -e $outdir);

$LOGGER->logdie("Missing barcodes file") if (! defined $barcodesfile);
$LOGGER->logdie("Wrong barcodes file ($barcodesfile)") if (! -e $barcodesfile);

my %sample;

open(BC, "<", $barcodesfile) or $LOGGER->logdie($!);

while(<BC>) {
	chomp;
	my ($bcseq, $sample_name) = split(/\t/, $_);
	$sample_name=~s/^\s*//;
	$sample_name=~s/\s*$//;

	$sample{$bcseq} = {	'name'=>$sample_name,
				'afh'=> FileHandle->new(">".$outdir.'/'.$sample_name.'.assembled.fastq'),
				'ufh'=> FileHandle->new(">".$outdir.'/'.$sample_name.'.unassembled.fastq'),
			  };
	autoflush { $sample{$bcseq}->{'afh'} } 1;
	autoflush { $sample{$bcseq}->{'ufh'} } 1;
	$LOGGER->logdie($!) if (! defined $sample{$bcseq}->{'afh'});
}

my @bc = keys %sample;

close(BC);

$|=1;
{
	my ($bcseq, $sample_name) = ('UNID', 'UNIDENTIFIED');
	$sample{$bcseq} = {	'name'=>$sample_name,
				'afh'=> FileHandle->new(">".$outdir.'/'.$sample_name.'.assembled.fastq'),
				'ufh'=> FileHandle->new(">".$outdir.'/'.$sample_name.'.unassembled.fastq') 
			  };
	
	$LOGGER->logdie($!) if (! defined $sample{$bcseq}->{'afh'});
	autoflush { $sample{$bcseq}->{'afh'} } 1;
	autoflush { $sample{$bcseq}->{'ufh'} } 1;
	($bcseq, $sample_name) = ('WTMP', 'WRONGTEMPLATE');
	$sample{$bcseq} = {	'name'=>$sample_name,
				'afh'=> FileHandle->new(">".$outdir.'/'.$sample_name.'.assembled.fastq'),
				'ufh'=> FileHandle->new(">".$outdir.'/'.$sample_name.'.unassembled.fastq') 
			  };
	
	$LOGGER->logdie($!) if (! defined $sample{$bcseq}->{'afh'});
	autoflush { $sample{$bcseq}->{'afh'} } 1;
	autoflush { $sample{$bcseq}->{'ufh'} } 1;
	
	($bcseq, $sample_name) = ('AMBI', 'AMBIGUOUS');
	$sample{$bcseq} = {	'name'=>$sample_name,
				'afh'=> FileHandle->new(">".$outdir.'/'.$sample_name.'.assembled.fastq'),
				'ufh'=> FileHandle->new(">".$outdir.'/'.$sample_name.'.unassembled.fastq') 
			  };
	
	$LOGGER->logdie($!) if (! defined $sample{$bcseq}->{'afh'});
	autoflush { $sample{$bcseq}->{'afh'} } 1;
	autoflush { $sample{$bcseq}->{'ufh'} } 1;
	
	($bcseq, $sample_name) = ('REFU', 'REFUSED');
	$sample{$bcseq} = {	'name'=>$sample_name,
				'afh'=> FileHandle->new(">".$outdir.'/'.$sample_name.'.assembled.fastq'),
				'ufh'=> FileHandle->new(">".$outdir.'/'.$sample_name.'.unassembled.fastq') 
			  };
	
	$LOGGER->logdie($!) if (! defined $sample{$bcseq}->{'afh'});
	autoflush { $sample{$bcseq}->{'afh'} } 1;
	autoflush { $sample{$bcseq}->{'ufh'} } 1;
}

my $bn = basename($r1file, '.fastq');
$bn=~s/_R1_/_/;

$maxmismatches||=1;
$LOGGER->logdie("Maximum number of mismatches is 2") if ($maxmismatches>2);

my $seqan = new SeqAn();

if (! -e $outdir.'/'.$bn.'.assembled.fastq') {
	print "Running PEAR to assemble fragments ...\n";
	`pear --threads 20  --quality-threshold 20 --min-overlap 10 -f $r1file -r $r2file -o $outdir/$bn --keep-original`;
	`paste $outdir/$bn.unassembled.forward.fastq $outdir/$bn.unassembled.reverse.fastq > $outdir/$bn.unassembled.txt`;
}


my $seqcountexe = 'echo "$(wc -l '.$outdir.'/'.$bn.'.assembled.fastq | sed \'s/ .*$//\')/4" | bc -l';
my $seqcount = `$seqcountexe`;
$seqcount+=0;
print "Processing assembled ($seqcount) fragments ...\n";

my $progress = Term::ProgressBar->new($seqcount);

open(IN, "<", $outdir.'/'.$bn.'.assembled.fastq') or $LOGGER->logdie($!);

my $id;
my $getqual_start;
my $getqual_end;
my $getqual_sample;
my $seqc = 0;


while(<IN>) {
	chomp;
	&printToFile($seqan, 'afh', $_, $., \$id, \%sample, \$seqc, \$getqual_start, \$getqual_end, \$getqual_sample, $progress);
}
close(IN);
$progress->update($seqc);

$seqcountexe = 'echo "$(wc -l '.$outdir.'/'.$bn.'.unassembled.txt | sed \'s/ .*$//\')/4" | bc -l';
$seqcount = `$seqcountexe`;
$seqcount+=0;
print "Processing unassembled ($seqcount) fragments ...\n";

$progress = Term::ProgressBar->new($seqcount);

open(IN, "<", $outdir.'/'.$bn.'.unassembled.txt') or $LOGGER->logdie($!);

$id = undef;
$getqual_start = undef;
$getqual_end = undef;
$getqual_sample = undef;
$seqc = 0;
while(<IN>) {
	chomp;
	
	if ( ($. % 4) == 2 ) {
		$_=~s/\s/NNNNN/g;
	} elsif ( ($. % 4) == 0 ) {
		$_=~s/\s/!!!!!/g;
	}
	
	&printToFile($seqan, 'ufh', $_, $., \$id, \%sample, \$seqc, \$getqual_start, \$getqual_end, \$getqual_sample, $progress);
}
close(IN);
$progress->update($seqc);

foreach my $s (keys %sample) {
	$sample{$s}->{'afh'}->close();
	$sample{$s}->{'ufh'}->close();
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
	-r1	--r1file	Read 1 (forward) file
	-r2	--r2file	Read 2 (reverse) file
	-o	--outdir	Output directory
	-b	--barcodes	Barcodes file
	-mm	--maxmismatches	Maximum number of mismatches [Default: 1]

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
    print STDERR $USAGE;
    exit(1);
}

sub printToFile {
	my ($seqan, $ifh, $l, $c, $sr_id, $hr_sample, $sr_seqc, $sr_getqual_start, $sr_getqual_end, $sr_getqual_sample, $progress) = @_;
	
	if ( ($c % 4) == 1 ) {
		(${$sr_id}) = $l=~/^@(\S+(?: \S+)?)/;
		${$sr_getqual_start} = undef;
		${$sr_getqual_end} = undef;
		${$sr_getqual_sample} = undef;
		${$sr_seqc}++;
		if ( (${$sr_seqc}%1000) == 0 ) {
			$progress->update(${$sr_seqc});
		}
	} elsif ( ($. % 4) == 2 ) {
		$LOGGER->logdie("Missing ID ($.): $_") unless (defined ${$sr_id});
		if (length($l) > 60) {
			my ($score, $aln1, $aln2) = split(';', $seqan->Align2Seq('tcNNNGGAGCTGTCGTATTCCAGTCAGG',$l));
			while ($aln2=~/^-/g) {
				$aln2=~s/^-//;
				$aln1=~s/^.//;
			}
			while ($aln2=~/-$/g) {
				$aln2=~s/-$//;
				$aln1=~s/.$//;
			}
			$aln1=~s/-//g;
			
			if ($score >= 30) {
				my %mismatch;
				my $selseq = '';
				if ($aln2=~/^TC-NNN/) {
					$aln2=~s/^TC-NNN/TCNNN-/;
				}
				while($aln2=~/(N)/g) {
					
					my $p = pos($aln2)-1;
					$selseq.=substr($aln1, $p, 1);
				}
				if ($selseq!~/^[ACGTN]+$/) {
					${$sr_getqual_sample} = 'WTMP';
					${$sr_getqual_start} = 0;
					${$sr_getqual_end} = length($l);
					print { $hr_sample->{ ${$sr_getqual_sample} }->{$ifh} } '@',${$sr_id},"\t$aln1\n",substr($l,${$sr_getqual_start},(${$sr_getqual_end}-${$sr_getqual_start})),"\n";
				} else {
					my @sorder;
					if (exists $hr_sample->{$selseq}) {
						@{ $mismatch{$selseq} } = (0, $selseq, $selseq);
						@sorder = ($selseq);
					} else {
						foreach my $kseq (@bc) {
							(@{ $mismatch{$kseq} })= &simpleAlign($selseq, $kseq);
						}
						
						@sorder = sort { $mismatch{$a}->[0] <=> $mismatch{$b}->[0] } keys %mismatch;
					}
					
					my $s = $sorder[0];
					if ($mismatch{$s}->[0] <= $maxmismatches) {
						if ((!defined $sorder[1])||($mismatch{$s}->[0] < $mismatch{ $sorder[1] }->[0])) {
							#print $s,"\t",$hr_sample->{$s}->{'name'},"\t",$mismatch{$s}->[0],"\t",$mismatch{$s}->[1],"\t",$mismatch{$s}->[2],"\n";
							my $start = 0;
							my $end = length($l);
							if ($l =~/($aln1)/g) {
								$start = pos($l);
							}
							my ($scorer, $aln1r, $aln2r) = split(';', $seqan->Align2Seq('CCTCTAAACGGGTCT',$l));
							while ($aln2r=~/^-/g) {
								$aln2r=~s/^-//;
								$aln1r=~s/^.//;
							}
							while ($aln2r=~/-$/g) {
								$aln2r=~s/-$//;
								$aln1r=~s/.$//;
							}
							$aln1r=~s/-//g;
							
							if ($l=~/($aln1r)/g) {
								$end = pos($l)-15; # CCTCTAAACGGGTCT (15 bases)
							}
							${$sr_getqual_sample} = $s;
							${$sr_getqual_start} = $start;
							${$sr_getqual_end} = $end;

							print { $hr_sample->{ ${$sr_getqual_sample} }->{$ifh} } '@',${$sr_id},"\t$aln1\n",substr($l,${$sr_getqual_start},(${$sr_getqual_end}-${$sr_getqual_start})),"\n";
						} else {
							${$sr_getqual_sample} = 'AMBI';
							${$sr_getqual_start} = 0;
							${$sr_getqual_end} = length($l);
							my @amb;
							foreach my $o (@sorder) {
								push(@amb, $o) if ($mismatch{$o}->[0] == $mismatch{$s}->[0]);
							}
							print { $hr_sample->{ ${$sr_getqual_sample} }->{$ifh} } '@',${$sr_id},"\t$aln1\t".join(';',@amb)."\n",substr($l,${$sr_getqual_start},(${$sr_getqual_end}-${$sr_getqual_start})),"\n";
						}
					} else {
						${$sr_getqual_sample} = 'UNID';
						${$sr_getqual_start} = 0;
						${$sr_getqual_end} = length($l);
						print { $hr_sample->{ ${$sr_getqual_sample} }->{$ifh} } '@',${$sr_id},"\t$aln1\n",substr($l,${$sr_getqual_start},(${$sr_getqual_end}-${$sr_getqual_start})),"\n";
					}
				}
			} else {
				${$sr_getqual_sample} = 'WTMP';
				${$sr_getqual_start} = 0;
				${$sr_getqual_end} = length($l);
				print { $hr_sample->{ ${$sr_getqual_sample} }->{$ifh} } '@',${$sr_id},"\t$aln1\n",substr($l,${$sr_getqual_start},(${$sr_getqual_end}-${$sr_getqual_start})),"\n";
			}
		} else {
			${$sr_getqual_sample} = 'REFU';
			${$sr_getqual_start} = 0;
			${$sr_getqual_end} = length($l);
			print { $hr_sample->{ ${$sr_getqual_sample} }->{$ifh} } '@',${$sr_id},"\t\n",substr($l,${$sr_getqual_start},(${$sr_getqual_end}-${$sr_getqual_start})),"\n";

		}
		
	} elsif ( ($. % 4) == 0 ) {
		if (defined ${$sr_getqual_sample}) {
			print { $hr_sample->{ ${$sr_getqual_sample} }->{$ifh} } '+',"\n",substr($l,${$sr_getqual_start},(${$sr_getqual_end}-${$sr_getqual_start})),"\n";
		}
		#exit;
	}
}


sub simpleAlign {
	my ($seq1, $seq2) = @_;
	my @s1 = split(//, $seq1);
	my @s2 = split(//, $seq2);
	my $nmm = 0;
	for (my $i=0; $i<=$#s1; $i++) {
		if ($s1[$i] ne $s2[$i]) {
			$nmm++;
		}
	}
	return($nmm, $seq1, $seq2);
}
