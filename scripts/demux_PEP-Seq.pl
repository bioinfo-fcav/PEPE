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
#  Copyright (C) 2020  Universidade Estadual Paulista "Júlio de Mesquita Filho"
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

Copyright (c) 2020 Universidade Estadual Paulista "Júlio de Mesquita Filho"

=head1 LICENSE

GNU General Public License

http://www.gnu.org/copyleft/gpl.html


=cut

use strict;
use warnings;
use Getopt::Long;

use File::Basename;
use File::Temp qw/ tempdir /;

use Bio::SeqIO;
use Bio::Seq::Quality;
use Bio::Seq;

use vars qw/$LOGGER/;
use Term::ProgressBar;
use PEPE::SeqAn;
use FileHandle;

INIT {
    use Log::Log4perl qw/:easy/;
    Log::Log4perl->easy_init($FATAL);
    $LOGGER = Log::Log4perl->get_logger($0);
}

# In MODEL_ADAPTER_SEQUENCE_5PRIME we need to have 3 consecutive Ns
use constant DEFAULT_MODEL_ADAPTER_SEQUENCE_5PRIME=>'tcNNNGGAGCTGTCGTATTCCAGTCAGG';
use constant DEFAULT_MODEL_ADAPTER_SEQUENCE_3PRIME=>'CCTCTAAACGGGTCTTGAGGGGTT';
use constant DEFAULT_MIN_ALN_SCORE=>30;
use constant DEFAULT_MIN_LENGTH=>60;
use constant DEFAULT_MAX_MISMATCHES=>1;
use constant DEFAULT_MERGER=>'pear';
use constant DEFAULT_THREADS_NUMBER => 1;
use constant DEFAULT_DIAMOND_QUERY_COVER => 100;
use constant DEFAULT_DIAMOND_ID => 90;
use constant DEFAULT_DIAMOND_EVALUE => 0.001;
#tcNNNGGAGCTGTCGTATTCCAGTCAGGTGTGATGCTCGGGGATCCGAATT (51)
use constant DEFAULT_TRIM_FWD_FOR_BLASTX=> 55;
#AACCCCTCAAGACCCGTTTAGAGGCCCCAAGGGGTTAACTAGTTACTCGAGTGCGGCCGCAAGCT (65)
use constant DEFAULT_TRIM_REV_FOR_BLASTX=> 70;
use constant DEFAULT_MIN_BASES_TO_KEEP_FOR_BLASTX=>30;
use constant DEFAULT_MIN_OVERLAP=>10;
use constant DEFAULT_QUAL_THRESHOLD=>10;
use constant DEFAULT_MERGE_MAX_DIFFS=>2;

my (	$level, 
	$r1file, 
	$r2file, 
	$outdir, 
	$barcodesfile, 
	$maxmismatches, 
	$minlength,
	$minalnscore,
	$merger, 
	$dmnd, 
	$nthreads, 
	$noalign, 
	$diamond_id,
	$diamond_evalue,
	$diamond_query_cover,
	$only_one_end,
	$nbases,
	$model_adapter_sequence_5prime,
	$model_adapter_sequence_3prime,
	$trim_fwd_for_blastx,
	$trim_rev_for_blastx,
	$overwrite,
	$minoverlap,
	$qualthreshold,
	$mergemaxdiffs,
	$pear_memory
	);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions(	"h|?|help" => sub { &Usage(); },
		"l|level=s"=> \$level,
		"r1|r1file=s"=>\$r1file,
		"r2|r2file=s"=>\$r2file,
		"o|outdir=s"=>\$outdir,
		"b|barcodes=s"=>\$barcodesfile,
		"m|merger=s"=>\$merger,
		"d|dmnd=s"=>\$dmnd,
		"mm|maxmismatches=s"=>\$maxmismatches,
		"ml|minlength=i"=>\$minlength,
		"ms|minalnscore=i"=>\$minalnscore,
		"t|threads=i"=>\$nthreads,
		"y|memory=s"=>\$pear_memory,
		"ooe|only-one-end"=>\$only_one_end,
		"nb|nbases"=>\$nbases,
		"did|diamond-id=f"=>\$diamond_id,
		"dev|diamond-evalue=f"=>\$diamond_evalue,
		"dqc|diamond-query-cover=f"=>\$diamond_query_cover,
		"m5p|model5pseq=s"=>\$model_adapter_sequence_5prime,
		"m3p|model3pseq=s"=>\$model_adapter_sequence_3prime,
		"tfr|trimfwdread=i"=>\$trim_fwd_for_blastx,
		"trr|trimrevread=i"=>\$trim_rev_for_blastx,
		"w|overwrite"=>\$overwrite,
		"mo|minoverlap=i"=>\$minoverlap,
		"qt|qualthreshold=i"=>\$qualthreshold,
		"md|mergemaxdiffs=i"=>\$mergemaxdiffs
    ) or &Usage();

$minoverlap||=DEFAULT_MIN_OVERLAP;
$trim_fwd_for_blastx||=DEFAULT_TRIM_FWD_FOR_BLASTX;
$trim_rev_for_blastx||=DEFAULT_TRIM_REV_FOR_BLASTX;
$model_adapter_sequence_5prime||=DEFAULT_MODEL_ADAPTER_SEQUENCE_5PRIME;
$model_adapter_sequence_3prime||=DEFAULT_MODEL_ADAPTER_SEQUENCE_3PRIME;
$minlength||=DEFAULT_MIN_LENGTH;
$minalnscore||=DEFAULT_MIN_ALN_SCORE;
$maxmismatches||=DEFAULT_MAX_MISMATCHES;
$merger||=DEFAULT_MERGER;
$nthreads||=DEFAULT_THREADS_NUMBER;
$diamond_id||=DEFAULT_DIAMOND_ID;
$diamond_evalue||=DEFAULT_DIAMOND_EVALUE;
$diamond_query_cover||=DEFAULT_DIAMOND_QUERY_COVER;
$qualthreshold||=DEFAULT_QUAL_THRESHOLD;
$mergemaxdiffs||=DEFAULT_MERGE_MAX_DIFFS;

unless ((defined $only_one_end)||(defined $nbases)) {
	$LOGGER->logdie("Missing DIAMOND index or choose -ooe/--only-one-end option") if (! defined $dmnd);
	$dmnd=~s/\.dmnd//;
	$LOGGER->logdie("Wrong DIAMOND index. DIAMOND index file ($dmnd.dmnd) doesn't exist.") if (! -e $dmnd.'.dmnd');
} else {
	if ((defined $only_one_end)&&(defined $nbases)) {
		$LOGGER->logdie("It is not possible to set -nb/--nbases together with -ooe/--only-one-end");
	}
}

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

$LOGGER->info("Running with $nthreads thread(s)");

open(BC, "<", $barcodesfile) or $LOGGER->logdie($!);

while(<BC>) {
	chomp;
	my ($bcseq, $sample_name) = split(/\t/, $_);
	$sample_name=~s/^\s*//;
	$sample_name=~s/\s*$//;

	if (-e $outdir.'/'.$sample_name.'.assembled.fastq') {
		unless ($overwrite) {
			$LOGGER->logdie("File exists (".$outdir.'/'.$sample_name.'.assembled.fastq'."). Please use parameter -w/--overwrite to overwrite.");
			$LOGGER->logdie("File exists (".$outdir.'/'.$sample_name.'.unassembled.fastq'."). Please use parameter -w/--overwrite to overwrite.");
		}
	}

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
$bn=~s/_R1(_?)/$1/;

$LOGGER->logdie("Maximum number of mismatches is 2") if ($maxmismatches>2);

my $seqan = new PEPE::SeqAn();

if (! -e $outdir.'/'.$bn.'.assembled.fastq') {
	print "Running ".uc($merger)." to assemble fragments ...\n";
	if ($merger eq 'pear') {
		my $memory="";
		if ($pear_memory) {
			$memory="--memory $pear_memory";
		}
		system("pear --threads $nthreads --quality-threshold $qualthreshold --min-overlap $minoverlap -f $r1file -r $r2file -o $outdir/$bn --keep-original $memory > $outdir/$bn.PEAR.log.out.txt 2> $outdir/$bn.PEAR.log.err.txt");
	} elsif ($merger eq 'vsearch') {
		system("vsearch --threads $nthreads --fastq_minovlen $minoverlap --fastq_allowmergestagger --fastq_truncqual $qualthreshold --fastq_maxdiffs $mergemaxdiffs --fastq_mergepairs $r1file --reverse $r2file --fastqout $outdir/$bn.assembled.fastq --fastqout_notmerged_fwd $outdir/$bn.unassembled.forward.fastq --fastqout_notmerged_rev $outdir/$bn.unassembled.reverse.fastq > $outdir/$bn.vsearch.log.out.txt 2> $outdir/$bn.vsearch.log.err.txt");
		
	} elsif ($merger eq 'usearch') {
		system("usearch -threads $nthreads -fastq_minovlen $minoverlap -fastq_allowmergestagger -fastq_truncqual $qualthreshold -fastq_maxdiffs $mergemaxdiffs -fastq_mergepairs $r1file -reverse $r2file -fastqout $outdir/$bn.assembled.fastq -fastqout_notmerged_fwd $outdir/$bn.unassembled.forward.fastq -fastqout_notmerged_rev $outdir/$bn.unassembled.reverse.fastq > $outdir/$bn.usearch.log.out.txt 2> $outdir/$bn.usearch.log.err.txt");
	} else {
		$LOGGER->logdie("Program to merge ($merger) not recognized! Please choose PEAR or usearch");
	}
	`paste $outdir/$bn.unassembled.forward.fastq $outdir/$bn.unassembled.reverse.fastq > $outdir/$bn.unassembled.txt`;
}


my $seqcountexe = 'echo "$(wc -l '.$outdir.'/'.$bn.'.assembled.fastq | sed \'s/ .*$//\')/4" | bc -l';
my $seqcount = `$seqcountexe`;
$seqcount+=0;
print "Processing assembled ($seqcount) fragments ...\n";

my $progress = Term::ProgressBar->new($seqcount);

if ( -e $outdir.'/'.$bn.'.assembled.fastq' ) {
	open(IN, "<", $outdir.'/'.$bn.'.assembled.fastq') or $LOGGER->logdie($!);

	my $id;
	my $getqual_start;
	my $getqual_end;
	my $getqual_sample;
	my $seqc = 0;


	while(<IN>) {
		chomp;
		&printToFile($seqan, 'afh', $_, '', $., \$id, \%sample, \$seqc, \$getqual_start, \$getqual_end, \$getqual_sample, $progress, undef, undef);
	}
	close(IN);
	$progress->update($seqc);
}

$seqcountexe = 'echo "$(wc -l '.$outdir.'/'.$bn.'.unassembled.txt | sed \'s/ .*$//\')/4" | bc -l';
$seqcount = `$seqcountexe`;
$seqcount+=0;
print "Processing unassembled ($seqcount) fragments ...\n";

my $hr_align;
unless ((defined $only_one_end)||(defined $nbases)) {
	print "	BLASTx alignment step for unassembled fragments ...\n";
	
	&trim("$outdir/$bn.unassembled.forward.fastq","$outdir/$bn.unassembled.forward.blastx.fastq", $trim_fwd_for_blastx, DEFAULT_MIN_BASES_TO_KEEP_FOR_BLASTX);
	&trim("$outdir/$bn.unassembled.reverse.fastq","$outdir/$bn.unassembled.reverse.blastx.fastq", $trim_rev_for_blastx, DEFAULT_MIN_BASES_TO_KEEP_FOR_BLASTX);

	$hr_align = &blastx("$outdir/$bn.unassembled.forward.blastx.fastq", "$outdir/$bn.unassembled.reverse.blastx.fastq", \$dmnd, $outdir);
}

$progress = Term::ProgressBar->new($seqcount);


if ( -e $outdir.'/'.$bn.'.unassembled.txt' ) {
	open(IN, "<", $outdir.'/'.$bn.'.unassembled.txt') or $LOGGER->logdie($!);

	my $id = undef;
	my $getqual_start = undef;
	my $getqual_end = undef;
	my $getqual_sample = undef;
	my $seqc = 0;
	my $lastseq1='';
	my $lastseqn='';
	my $lastseq2='';
	my $lastsign='';
	my $nline=0;
	my %scfseq;
	while(<IN>) {
		chomp;
		$nline++;

		my ($seq1, $seqn, $seq2, $sign);
		
		($seq1, $seq2) = split(/\t/, $_);
		
		#&printToFileDMND($seqan, 'ufh', $seq1,$seq2, $nline, \$id, \%sample, \$seqc, \$getqual_start, \$getqual_end, \$getqual_sample, $progress,\$dmnd,\%scfseq);
		&printToFile($seqan, 'ufh', $seq1, $seq2, $nline, \$id, \%sample, \$seqc, \$getqual_start, \$getqual_end, \$getqual_sample, $progress, $hr_align, \%scfseq);
	}
	close(IN);
	$progress->update($seqc);
}

foreach my $s (keys %sample) {
	$sample{$s}->{'afh'}->close() if ($sample{$s}->{'afh'});
	$sample{$s}->{'ufh'}->close() if ($sample{$s}->{'ufh'});
}

# Subroutines

sub Usage {
    my ($msg) = @_;

	$minlength=DEFAULT_MIN_LENGTH;
	$minalnscore=DEFAULT_MIN_ALN_SCORE;
	$maxmismatches=DEFAULT_MAX_MISMATCHES;
	$merger=DEFAULT_MERGER;
	$nthreads=DEFAULT_THREADS_NUMBER;
	$diamond_id=DEFAULT_DIAMOND_ID;
	$diamond_evalue=DEFAULT_DIAMOND_EVALUE;
	$diamond_query_cover=DEFAULT_DIAMOND_QUERY_COVER;
	$model_adapter_sequence_5prime=DEFAULT_MODEL_ADAPTER_SEQUENCE_5PRIME;
	$model_adapter_sequence_3prime=DEFAULT_MODEL_ADAPTER_SEQUENCE_3PRIME;
	$trim_fwd_for_blastx=DEFAULT_TRIM_FWD_FOR_BLASTX;
	$trim_rev_for_blastx=DEFAULT_TRIM_REV_FOR_BLASTX;
	$minoverlap=DEFAULT_MIN_OVERLAP;
	$qualthreshold=DEFAULT_QUAL_THRESHOLD;
	$mergemaxdiffs=DEFAULT_MERGE_MAX_DIFFS;
	
    my $USAGE = <<"END_USAGE";
Daniel Guariz Pinheiro (dgpinheiro\@gmail.com)
(c)2020 Universidade Estadual Paulista "Júlio de Mesquita Filho"

Usage

        $0	[-h/--help] [-l/--level <LEVEL>]

Argument(s)

	-h	--help			Help
	-t	--threads		Number of threads [Default: $nthreads]
	-y	--memory		Amount of memory to be used by pear
	-m	--merger		Program to merge R1 and R2 (pear, usearch or vsearch) [Default: $merger]
	-mo	--minoverlap		Minimum overlap required for merging the two reads [Default: $minoverlap]
	-qt	--qualthreshold		Quality threshold (--quality-threshold) for pear / Truncation by quality (--fastq_truncqual) for usearch/vsearch [Default: $qualthreshold]
	-md	--mergemaxdiffs		Merger (usearch/vsearch) max differences in overlap [Default: $mergemaxdiffs]
	-l	--level			Log level [Default: FATAL]
	-ooe	--only-one-end		Choose only the read 1 (R1) of the paired-ends for unassembled fragments
	-nb	--nbases		Concatenate unassembled R1 and R2 splitted by a number of N bases
	-d	--dmnd			DIAMOND index
	-r1	--r1file		Read 1 (forward) file
	-r2	--r2file		Read 2 (reverse) file
	-o	--outdir		Output directory
	-b	--barcodes		Barcodes file
	-mm	--maxmismatches		Maximum number of mismatches [Default: $maxmismatches]
	-ms	--minalnscore		Minimum SeqAn alignment score [Default: $minalnscore]
	-ml	--minlength		Minimum read size [Default: $minlength]
	-n	--noalign		No DIAMOND alignment
	-did	--diamond-id		DIAMOND --id [Default: $diamond_id]
	-dev	--diamond-evalue	DIAMOND --evalue [Default: $diamond_evalue]
	-dqc	--diamond-query-cover	DIAMOND --query-cover [Default: $diamond_query_cover]
	-m5p	--model5pseq		Model adapter sequence 5' [Default: $model_adapter_sequence_5prime]
	-m3p	--model3pseq		Model adapter sequence 3' [Default: $model_adapter_sequence_3prime]
	-tfr	--trimfwdread		Number of bases to trim in forward reads for blastx [Default: $trim_fwd_for_blastx]
	-trr	--trimrevread		Number of bases to trim in reverse reads for blastx [Default: $trim_rev_for_blastx]
	-w	--overwrite		Overwrite previous results

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
    print STDERR $USAGE;
    exit(1);
}

sub printToFile {
	my ($seqan, $ifh, $s1, $s2, $c, $sr_id, $hr_sample, $sr_seqc, $sr_getqual_start, $sr_getqual_end, $sr_getqual_sample, $progress, $hr_align, $hr_scfseq) = @_;
	if ( ($c % 4) == 1 ) {
		(${$sr_id}) = $s1=~/^@(\S+(?: \S+)?)/;
		${$sr_getqual_start} = undef;
		${$sr_getqual_end} = undef;
		${$sr_getqual_sample} = undef;
		${$sr_seqc}++;
		if ( (${$sr_seqc}%1000) == 0 ) {
			$progress->update(${$sr_seqc});
		}
	} elsif ( ($c % 4) == 2 ) {
		$LOGGER->logdie("Missing ID ($c): $_") unless (defined ${$sr_id});
		if (length($s1) >= $minlength) {
			my ($score, $aln1, $aln2) = split(';', $seqan->Align2Seq($model_adapter_sequence_5prime, $s1));
			
			while ($aln2=~/^-/g) {
				$aln2=~s/^-//;
				$aln1=~s/^.//;
			}
			while ($aln2=~/-$/g) {
				$aln2=~s/-$//;
				$aln1=~s/.$//;
			}
			$aln1=~s/-//g;
			
			if ($score >= $minalnscore) {
				my %mismatch;
				my $selseq = '';

				if ($aln2=~/^[ACGTN]+-NNN/) {
					$aln2=~s/^([ACGTN+])(-+)NNN/$1NNN$2/;
				}
				while($aln2=~/(N)/g) {
					
					my $p = pos($aln2)-1;
					$selseq.=substr($aln1, $p, 1);
				}
				if ($selseq!~/^[ACGTN]+$/) {
					${$sr_getqual_sample} = 'WTMP';
					${$sr_getqual_start} = 0;
					${$sr_getqual_end} = length($s1);
					print { $hr_sample->{ ${$sr_getqual_sample} }->{$ifh} } '@',${$sr_id},"\t$aln1\n",substr($s1,${$sr_getqual_start},(${$sr_getqual_end}-${$sr_getqual_start})),"\n";
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
						my $l;
						
						if ($ifh eq 'ufh') {
							my $sid = ${$sr_id};
							$sid=~s/\s+.*$//;

							if ($only_one_end) {
								$l=$s1;
							} else {
								my ($seq1, $seqn, $seq2, $sign, $sbjctid)=&scaffoldSeq(\$s1, \$s2, $hr_align->{ $sid }, $nbases);
						
								if ($sbjctid) {

									#$LOGGER->logdie("Not defined \$seq2 for $sid ! ($seq1;$seqn)") unless (defined $seq2);
									
									$l=$seq1.$seqn.$seq2;
									$hr_scfseq->{ ${$sr_id} } = { 'seq1'=>$seq1,
												      'seq2'=>$seq2,
												      'seqn'=>$seqn
												    };
								} else {
									$l=$s1;
								}
							}
								
						} else {
							$l=$s1;
						}
						
						if ((!defined $sorder[1])||($mismatch{$s}->[0] < $mismatch{ $sorder[1] }->[0])) {
							#print $s,"\t",$hr_sample->{$s}->{'name'},"\t",$mismatch{$s}->[0],"\t",$mismatch{$s}->[1],"\t",$mismatch{$s}->[2],"\n";
							my $start = 0;
							my $end = length($l);
							if ($l =~/($aln1)/g) {
								$start = pos($l);
							}
							my ($scorer, $aln1r, $aln2r) = split(';', $seqan->Align2Seq($model_adapter_sequence_3prime, $l));
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
								$end = pos($l)-length($model_adapter_sequence_3prime);
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
						${$sr_getqual_end} = length($s1);
						print { $hr_sample->{ ${$sr_getqual_sample} }->{$ifh} } '@',${$sr_id},"\t$aln1\n",substr($s1,${$sr_getqual_start},(${$sr_getqual_end}-${$sr_getqual_start})),"\n";
					}
				}
			} else {
				${$sr_getqual_sample} = 'WTMP';
				${$sr_getqual_start} = 0;
				${$sr_getqual_end} = length($s1);
				print { $hr_sample->{ ${$sr_getqual_sample} }->{$ifh} } '@',${$sr_id},"\t$aln1\n",substr($s1,${$sr_getqual_start},(${$sr_getqual_end}-${$sr_getqual_start})),"\n";
			}
		} else {
			${$sr_getqual_sample} = 'REFU';
			${$sr_getqual_start} = 0;
			${$sr_getqual_end} = length($s1);
			print { $hr_sample->{ ${$sr_getqual_sample} }->{$ifh} } '@',${$sr_id},"\t\n",substr($s1,${$sr_getqual_start},(${$sr_getqual_end}-${$sr_getqual_start})),"\n";

		}
		
	} elsif ( ($c % 4) == 0 ) {
		if (defined ${$sr_getqual_sample}) {
			my $l;
			if ($ifh eq 'ufh') {
				if ( $hr_scfseq->{ ${$sr_id} } ) {
					
					$LOGGER->logdie("Not defined \$s2 for ${$sr_id} !!!") unless (defined $s2);
					
					my $seq1=substr($s1, 0, length($hr_scfseq->{ ${$sr_id} }->{'seq1'}));
					my $seqn="!" x length($hr_scfseq->{ ${$sr_id} }->{'seqn'});
					my $seq2=reverse(substr($s2, 0, length( $hr_scfseq->{ ${$sr_id} }->{'seq2'} )));
									
					$LOGGER->logdie("Not defined \$seq2 for ${$sr_id} !!!") unless (defined $seq2);
					
					$l=$seq1.$seqn.$seq2;
				} else {
					$l=$s1;
				}
					
			} else {
				$l=$s1;
			}
			
			print { $hr_sample->{ ${$sr_getqual_sample} }->{$ifh} } '+',"\n",substr($l,${$sr_getqual_start},(${$sr_getqual_end}-${$sr_getqual_start})),"\n";

			

		} else {
			$LOGGER->logdie("Not defined quality coordinates to print fastq file for ${$sr_id}");
		}
	}
}

#sub printToFileDMND {
#	my ($seqan, $ifh, $s1, $s2, $c, $sr_id, $hr_sample, $sr_seqc, $sr_getqual_start, $sr_getqual_end, $sr_getqual_sample, $progress, $rdmnd,$hr_scfseq) = @_;
#	if ( ($c % 4) == 1 ) {
#		(${$sr_id}) = $s1=~/^@(\S+(?: \S+)?)/;
#		${$sr_getqual_start} = undef;
#		${$sr_getqual_end} = undef;
#		${$sr_getqual_sample} = undef;
#		${$sr_seqc}++;
#		if ( (${$sr_seqc}%1000) == 0 ) {
#			$progress->update(${$sr_seqc});
#		}
#	} elsif ( ($c % 4) == 2 ) {
#		$LOGGER->logdie("Missing ID ($c): $_") unless (defined ${$sr_id});
#		if (length($s1) >= $minlength) {
#			my ($score, $aln1, $aln2) = split(';', $seqan->Align2Seq('tcNNNGGAGCTGTCGTATTCCAGTCAGG',$s1));
#			
#			while ($aln2=~/^-/g) {
#				$aln2=~s/^-//;
#				$aln1=~s/^.//;
#			}
#			while ($aln2=~/-$/g) {
#				$aln2=~s/-$//;
#				$aln1=~s/.$//;
#			}
#			$aln1=~s/-//g;
#			
#			if ($score >= DEFAULT_MIN_ALN_SCORE) {
#				my %mismatch;
#				my $selseq = '';
#				if ($aln2=~/^TC-NNN/) {
#					$aln2=~s/^TC-NNN/TCNNN-/;
#				}
#				while($aln2=~/(N)/g) {
#					
#					my $p = pos($aln2)-1;
#					$selseq.=substr($aln1, $p, 1);
#				}
#				if ($selseq!~/^[ACGTN]+$/) {
#					${$sr_getqual_sample} = 'WTMP';
#					${$sr_getqual_start} = 0;
#					${$sr_getqual_end} = length($s1);
#					print { $hr_sample->{ ${$sr_getqual_sample} }->{$ifh} } '@',${$sr_id},"\t$aln1\n",substr($s1,${$sr_getqual_start},(${$sr_getqual_end}-${$sr_getqual_start})),"\n";
#				} else {
#					my @sorder;
#					if (exists $hr_sample->{$selseq}) {
#						@{ $mismatch{$selseq} } = (0, $selseq, $selseq);
#						@sorder = ($selseq);
#					} else {
#						foreach my $kseq (@bc) {
#							(@{ $mismatch{$kseq} })= &simpleAlign($selseq, $kseq);
#						}
#						
#						@sorder = sort { $mismatch{$a}->[0] <=> $mismatch{$b}->[0] } keys %mismatch;
#					}
#					
#					my $s = $sorder[0];
#					if ($mismatch{$s}->[0] <= $maxmismatches) {
#						my $l;
#						
#						if ($ifh eq 'ufh') {
#							my ($seq1, $seqn, $seq2, $sign, $sbjctid)=&scaffoldSeqWithDiamond(\$s1, \$s2, $rdmnd);
#						
#							if ($sbjctid) {
#								$l=$seq1.$seqn.$seq2;
#								$hr_scfseq->{ ${$sr_id} } = { 'seq1'=>$seq1,
#											      'seq2'=>$seq2,
#											      'seqn'=>$seqn
#											    };
#							} else {
#								$l=$s1;
#							}
#								
#						} else {
#							$l=$s1;
#						}
#						
#						if ((!defined $sorder[1])||($mismatch{$s}->[0] < $mismatch{ $sorder[1] }->[0])) {
#							#print $s,"\t",$hr_sample->{$s}->{'name'},"\t",$mismatch{$s}->[0],"\t",$mismatch{$s}->[1],"\t",$mismatch{$s}->[2],"\n";
#							my $start = 0;
#							my $end = length($l);
#							if ($l =~/($aln1)/g) {
#								$start = pos($l);
#							}
#							my ($scorer, $aln1r, $aln2r) = split(';', $seqan->Align2Seq('CCTCTAAACGGGTCT',$l));
#							while ($aln2r=~/^-/g) {
#								$aln2r=~s/^-//;
#								$aln1r=~s/^.//;
#							}
#							while ($aln2r=~/-$/g) {
#								$aln2r=~s/-$//;
#								$aln1r=~s/.$//;
#							}
#							$aln1r=~s/-//g;
#							
#							if ($l=~/($aln1r)/g) {
#								$end = pos($l)-15; # CCTCTAAACGGGTCT (15 bases)
#							}
#							${$sr_getqual_sample} = $s;
#							${$sr_getqual_start} = $start;
#							${$sr_getqual_end} = $end;
#
#							print { $hr_sample->{ ${$sr_getqual_sample} }->{$ifh} } '@',${$sr_id},"\t$aln1\n",substr($l,${$sr_getqual_start},(${$sr_getqual_end}-${$sr_getqual_start})),"\n";
#						} else {
#							${$sr_getqual_sample} = 'AMBI';
#							${$sr_getqual_start} = 0;
#							${$sr_getqual_end} = length($l);
#							my @amb;
#							foreach my $o (@sorder) {
#								push(@amb, $o) if ($mismatch{$o}->[0] == $mismatch{$s}->[0]);
#							}
#							print { $hr_sample->{ ${$sr_getqual_sample} }->{$ifh} } '@',${$sr_id},"\t$aln1\t".join(';',@amb)."\n",substr($l,${$sr_getqual_start},(${$sr_getqual_end}-${$sr_getqual_start})),"\n";
#						}
#					} else {
#						${$sr_getqual_sample} = 'UNID';
#						${$sr_getqual_start} = 0;
#						${$sr_getqual_end} = length($s1);
#						print { $hr_sample->{ ${$sr_getqual_sample} }->{$ifh} } '@',${$sr_id},"\t$aln1\n",substr($s1,${$sr_getqual_start},(${$sr_getqual_end}-${$sr_getqual_start})),"\n";
#					}
#				}
#			} else {
#				${$sr_getqual_sample} = 'WTMP';
#				${$sr_getqual_start} = 0;
#				${$sr_getqual_end} = length($s1);
#				print { $hr_sample->{ ${$sr_getqual_sample} }->{$ifh} } '@',${$sr_id},"\t$aln1\n",substr($s1,${$sr_getqual_start},(${$sr_getqual_end}-${$sr_getqual_start})),"\n";
#			}
#		} else {
#			${$sr_getqual_sample} = 'REFU';
#			${$sr_getqual_start} = 0;
#			${$sr_getqual_end} = length($s1);
#			print { $hr_sample->{ ${$sr_getqual_sample} }->{$ifh} } '@',${$sr_id},"\t\n",substr($s1,${$sr_getqual_start},(${$sr_getqual_end}-${$sr_getqual_start})),"\n";
#
#		}
#		
#	} elsif ( ($c % 4) == 0 ) {
#		if (defined ${$sr_getqual_sample}) {
#			my $l;
#			if ($ifh eq 'ufh') {
#				if ( $hr_scfseq->{ ${$sr_id} } ) {
#					my ($seq1, $seq2, $seqn);
#					
#					$seq1=substr($s1, 0, length($hr_scfseq->{ ${$sr_id} }->{'seq1'}));
#					$seqn="!" x length($hr_scfseq->{ ${$sr_id} }->{'seqn'});
#					$seq2=reverse(substr($s2, 0, length( $hr_scfseq->{ ${$sr_id} }->{'seq2'} )));
#					
#					$l=$seq1.$seqn.$seq2;
#				} else {
#					$l=$s1;
#				}
#					
#			} else {
#				$l=$s1;
#			}
#			
#			print { $hr_sample->{ ${$sr_getqual_sample} }->{$ifh} } '+',"\n",substr($l,${$sr_getqual_start},(${$sr_getqual_end}-${$sr_getqual_start})),"\n";
#
#			
#
#		} else {
#			$LOGGER->logdie("Not defined quality coordinates to print fastq file for ${$sr_id}");
#		}
#	}
#}


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

#sub scaffoldSeqWithDiamond {
#	my ($rr1, $rr2, $rdmnd) = @_;
#	my $tmpdir = tempdir( CLEANUP => 1 );
#	
#	if (length(${$rr1})>=57) {
#		my $seqr1 = Bio::Seq->new(	-display_id => 'R1',
#		                         	-seq => ${$rr1});
#		open(R1, ">", "$tmpdir/R1.fa") or $LOGGER->logdie($!);
#		print R1 ">R1\n${$rr1}";
#		close(R1);
#		`diamond blastx --threads $nthreads --db ${$rdmnd} --out $tmpdir/R1.txt --outfmt 6 qseqid qlen sseqid slen qstart qend sstart send bitscore length --query $tmpdir/R1.fa --query-cover 0.6 --more-sensitive --top 1 &>/dev/null`;
#		open(OUT, "<", "$tmpdir/R1.txt") or $LOGGER->logdie($!);
#		my %match;
#		while(<OUT>) {
#			chomp;
#			my ($qseqid, $qlen, $sseqid, $slen, $qstart, $qend, $sstart, $send, $bitscore, $length) = split(/\t/, $_);
#			$LOGGER->logwarn("Found more than one match of R1 with $sseqid") if (exists $match{$sseqid});
#			$match{$sseqid} = {	'subject' => [$sstart, $send],
#						'query'   => [$qstart, $qend]
#						};
#		}
#		close(OUT);
#		if (length(${$rr2})>=57) {
#			my $seqr2 = Bio::Seq->new(	-display_id => 'R2',
#							-seq => ${$rr2});
#			open(R2, ">", "$tmpdir/R2.fa") or $LOGGER->logdie($!);
#			print R2 ">R2\n${$rr2}";
#			close(R2);
#			`diamond blastx --threads $nthreads --db ${$rdmnd} --out $tmpdir/R2.txt --outfmt 6 qseqid qlen sseqid slen qstart qend sstart send bitscore length --query $tmpdir/R2.fa --query-cover 0.6 --more-sensitive --top 1 &>/dev/null`;
#
#			open(OUT, "<", "$tmpdir/R2.txt") or $LOGGER->logdie($!);
#			while(<OUT>) {
#				chomp;
#				my ($qseqid, $qlen, $sseqid, $slen, $qstart, $qend, $sstart, $send, $bitscore, $length) = split(/\t/, $_);
#				if (exists $match{$sseqid}) {
#					if (($match{$sseqid}->{'query'}->[0] < $match{$sseqid}->{'query'}->[1]) &&
#					    ($qstart > $qend) ) {
#
#						my $rseq1=$seqr1->trunc(1, $match{$sseqid}->{'query'}->[1])->seq();
#						my $rseq2=$seqr2->trunc(1, $qstart)->revcom()->seq();
#						my $n = (($sstart-$match{$sseqid}->{'subject'}->[1]-1)*3);
#						my $rseqn="";
#						if ($n>0) {
#							$rseqn="N"x$n;
#						} else {
#							$rseq2 = substr($rseq2,abs($n),length($rseq2));
#						}
#						return($rseq1,$rseqn,$rseq2,'+',$sseqid);
#					}
#					elsif (($match{$sseqid}->{'query'}->[0] > $match{$sseqid}->{'query'}->[1]) &&
#						 ($qstart < $qend) ) {
#						
#						my $rseq1=$seqr1->trunc(1, $match{$sseqid}->{'query'}->[0])->seq();
#						my $rseq2=$seqr2->trunc(1, $qend)->revcom()->seq();
#						
#						my $n = (($match{$sseqid}->{'subject'}->[0]-$send-1)*3);
#						my $rseqn="";
#						if ($n>0) {
#							$rseqn="N"x$n;
#						} else {
#							$rseq2 = substr($rseq2,abs($n),length($rseq2));
#						}
#						
#						return($rseq1,$rseqn,$rseq2,'-',$sseqid);
#					}
#				}
#			}
#			close(OUT);
#			
#		}
#	}
#	return(${$rr1}, '', '', '+', '');
#}


sub scaffoldSeq {
	my ($rr1, $rr2, $hr_seq_align, $pnbases) = @_;
	
	
	my %match;
		
	my $seqr1 = Bio::Seq->new(	-display_id => 'R1',
					-seq => ${$rr1});
	my $seqr2 = Bio::Seq->new(	-display_id => 'R2',
					-seq => ${$rr2});
	
	if ($hr_seq_align) {
		my $ar_1 = $hr_seq_align->{'R1'};
		my $ar_2 = $hr_seq_align->{'R2'};

		if (defined $ar_1) {
			foreach my $aln1 (@{$ar_1}) {
				chomp($aln1);
				my ($qseqid, $qlen, $sseqid, $slen, $qstart, $qend, $sstart, $send, $bitscore, $length) = split(/\t/, $aln1);
				if (exists $match{$sseqid}) {
					next;
					#$LOGGER->logwarn("Found more than one match of R1 ($qseqid) with $sseqid");
				} else {
					$match{$sseqid} = {	'subject' => [$sstart, $send],
								'query'   => [$qstart+$trim_fwd_for_blastx, $qend+$trim_fwd_for_blastx]
							};
				}
			}
		}

		if (defined $ar_2) {
			foreach my $aln2 (@{$ar_2}) {
				chomp($aln2);
				my ($qseqid, $qlen, $sseqid, $slen, $qstart, $qend, $sstart, $send, $bitscore, $length) = split(/\t/, $aln2);
				$qstart+=$trim_rev_for_blastx;
				$qend+=$trim_rev_for_blastx;
				if (exists $match{$sseqid}) {
					if (($match{$sseqid}->{'query'}->[0] < $match{$sseqid}->{'query'}->[1]) &&
					    ($qstart > $qend) ) {

						my $rseq1=$seqr1->trunc(1, $match{$sseqid}->{'query'}->[1])->seq();
						my $rseq2=$seqr2->trunc(1, $qstart)->revcom()->seq();
						my $n = (($sstart-$match{$sseqid}->{'subject'}->[1]-1)*3);
						my $rseqn="";
						if ($n>0) {
							$rseqn="N"x$n;
						} else {
							$rseq2 = substr($rseq2,abs($n),length($rseq2));
						}
						
						return($rseq1,$rseqn,$rseq2,'+',$sseqid);
						
					} elsif (($match{$sseqid}->{'query'}->[0] > $match{$sseqid}->{'query'}->[1]) &&
						 ($qstart < $qend) ) {
						
						my $rseq1=$seqr1->trunc(1, $match{$sseqid}->{'query'}->[0])->seq();
						my $rseq2=$seqr2->trunc(1, $qend)->revcom()->seq();
						
						my $n = (($match{$sseqid}->{'subject'}->[0]-$send-1)*3);

						if (abs($n)>length($rseq2)) {
							$n=length($rseq2);
						}
						my $rseqn="";
						if ($n>0) {
							$rseqn="N"x$n;
						} else {
							$rseq2 = substr($rseq2,abs($n),length($rseq2));
						}
						
						return($rseq1,$rseqn,$rseq2,'-',$sseqid);
					}
				}
			}
		}
		
	} elsif (defined $pnbases) {

		my $rseq1=$seqr1->seq();
		my $rseq2=$seqr2->revcom()->seq();
		my $rseqn="N"x$pnbases;
		
		return($rseq1, $rseqn, $rseq2, '+', '');
	} 
	
	return(${$rr1}, '', '', '+', '');
}

sub blastx {
	my ($rf1, $rf2, $rdmnd, $outdir) = @_;

	my $blastdir = tempdir( DIR => $outdir, CLEANUP => 1 );	
	
	`diamond blastx --threads $nthreads --db ${$rdmnd} --out $blastdir/BLASTX_R1.txt --outfmt 6 qseqid qlen sseqid slen qstart qend sstart send bitscore length --query $rf1 --evalue $diamond_evalue --id $diamond_id --query-cover $diamond_query_cover --more-sensitive --top 1 &>/dev/null`;
	
	`diamond blastx --threads $nthreads --db ${$rdmnd} --out $blastdir/BLASTX_R2.txt --outfmt 6 qseqid qlen sseqid slen qstart qend sstart send bitscore length --query $rf2 --evalue $diamond_evalue --id $diamond_id --query-cover $diamond_query_cover --more-sensitive --top 1 &>/dev/null`;
	
	my %BLASTX;
	
	open(IN, "<", "$blastdir/BLASTX_R1.txt") or $LOGGER->logdie($!);
	while(<IN>) {
		chomp;
		my (@f) = split(/\t/, $_);
		push(@{ $BLASTX{$f[0]}->{'R1'} }, $_);
	}
	close(IN);
	
	open(IN, "<", "$blastdir/BLASTX_R2.txt") or $LOGGER->logdie($!);
	while(<IN>) {
		chomp;
		my (@f) = split(/\t/, $_);
		push(@{ $BLASTX{$f[0]}->{'R2'} }, $_);
	}
	close(IN);
	
	return(\%BLASTX);
}

sub trim {
	my ($fin, $fout, $tsize, $msize) = @_;

	my $seqout = Bio::SeqIO->new(-file=>'>'.$fout, -format=>'FASTQ');
	
	open(IN, "<", $fin) or $LOGGER->logdie($!);
	my $seqid;
	my $seqstr;
	my $seqqual;
	while(<IN>) {
		my $lt = $.%4;
		chomp($_);
		if ($lt == 1) {
			$seqid=$_;
			$seqid=~s/^@//;
		} elsif ($lt == 2) {
			$seqstr=$_;
		} elsif ($lt == 3) {
			#
		} elsif ($lt == 0) {
			$seqqual=$_;

			if ($seqstr) {
				my $seq = Bio::Seq::Quality->new(-id=>$seqid,
								 -seq=>$seqstr,
								 -raw_quality=>$seqqual);
				my $len = $seq->length();
				if ($len >= ($tsize+$msize)) {
					$seqout->write_seq( $seq->trunc($tsize+1, $len) );
				}
			}
			
		} else {
			$LOGGER->logdie("Error: This line type ($lt) is not possible.");
		}
	}
}
