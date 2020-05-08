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

use lib '../blib';

use PEPE::SeqAn;


INIT {
    use Log::Log4perl qw/:easy/;
    Log::Log4perl->easy_init($FATAL);
    $LOGGER = Log::Log4perl->get_logger($0);
}

use constant DEFAULT_MODEL_ORF_START=>'ATGCTCGGGGATCCGAATTNN';
use constant DEFAULT_VECTOR_ORF_END=>'AAGCTTGCGGCCGCACTCGAG';
use constant DEFAULT_VECTOR_ORF_INNER=>'CCTGCAGGGATATCCCGGGAGCTCGTCGAC';
use constant DEFAULT_MIN_SCORE_VECTOR_ORF_INNER=>40;
use constant DEFAULT_MIN_SEQ_LENGTH=>45;

my %stop;
@stop{'TAG','TAA', 'TGA' } = ("amber", "ochre", "opal");

my ($level, $infile, $outdir, $suffix, $sample, $minseqlen, $model_orf_start, $vector_orf_end, $vector_orf_inner, $minscore_vector_orf_inner);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
	    "i|infile=s"=>\$infile,
	    "o|outdir=s"=>\$outdir,
	    "s|sample=s"=>\$sample,
	    "x|suffix=s"=>\$suffix,
            "m|minseqlen=s"=>\$minseqlen,
            "mos|modorfstart=s"=>\$model_orf_start,
            "voe|vectororfend=s"=>\$vector_orf_end,
            "voi|vectororfinner=s"=>\$vector_orf_inner,
            "msv|minscorevoi=s"=>\$minscore_vector_orf_inner
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

$minseqlen||=DEFAULT_MIN_SEQ_LENGTH;
$model_orf_start||=DEFAULT_MODEL_ORF_START;
$vector_orf_end||=DEFAULT_VECTOR_ORF_END;
$vector_orf_inner||=DEFAULT_VECTOR_ORF_INNER;
$minscore_vector_orf_inner||=DEFAULT_MIN_SCORE_VECTOR_ORF_INNER;

my $seqo_model_orf_start = Bio::Seq->new(-display_id => 'MODEL_ORF_START',
					 -seq => $model_orf_start,
					 -alphabet=>'dna'
					);

my $model_cds_start=$seqo_model_orf_start->translate(-frame=>0)->seq();
$model_cds_start=~s/X/\./g;

my $MODEL_CDS_START_RE=qr/($model_cds_start)/;
my $VECTOR_ORF_END_RE=qr/($vector_orf_end)/;

$LOGGER->logdie("Missing input fastq file") if (!defined $infile);
$LOGGER->logdie("Wrong input fastq file ($infile)") if (! -e $infile);

$LOGGER->logdie("Missing output directory") if (! defined $outdir);
$LOGGER->logdie("Wrong output directory ($outdir)") if (! -e $outdir);

$LOGGER->logdie("Missing sample name") if (! defined $sample);

$suffix||='.fastq';
							


#my $seqanx = new PEPE::SeqAn();
#my ($score_vector_orf_inner) = split(';', $seqanx->Align2Seq($vector_orf_inner, 'CCTATAGTAAGGTCCCGGGAGCTCGTCGAC','L')); 
#print $score_vector_orf_inner,"\n";
#exit;

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

$outfh{'smallseq'} = FileHandle->new(">".$outdir.'/'.$sample.'.smallseq'.$suffix);
autoflush { $outfh{'smallseq'} } 1;
$LOGGER->logdie($!) if (! defined $outfh{'smallseq'});

$outfh{'phage'} = FileHandle->new(">".$outdir.'/'.$sample.'.phage'.$suffix);
autoflush { $outfh{'phage'} } 1;
$LOGGER->logdie($!) if (! defined $outfh{'phage'});

$outfh{'empty'} = FileHandle->new(">".$outdir.'/'.$sample.'.empty'.$suffix);
autoflush { $outfh{'empty'} } 1;
$LOGGER->logdie($!) if (! defined $outfh{'empty'});

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

my $seqan = new PEPE::SeqAn();

$i=0;
foreach my $ar_read (@read) {
	
	my $status=undef;
	
	if ( ($i%10) == 0 ) {
		$progress->update($i);
	}
	$i++;
	my $seq = join('', @{ $ar_read->[1] });
	my $qual = join('', @{ $ar_read->[2] });
	my $id = $ar_read->[0]->[0];
	my $tag = $ar_read->[0]->[1];
	
	my $seqo = Bio::Seq->new(-display_id => $id,
				 -seq => $seq,
				 -alphabet=>'dna'
				);

	if ($seqo->length() < $minseqlen) {

		$status='smallseq';
		print { $outfh{ $status } } '@',$id,"\n",$seq,"\n","+\n",$qual,"\n";

	} else {
		
		my @aa;
		foreach my $frame (0,1,2) {
			my $prot_obj = $seqo->translate(-frame=>$frame);
			my $aaseq = $prot_obj->seq();
			if ($aaseq =~ /$MODEL_CDS_START_RE/g) {
				my $start = pos($aaseq);
				my $end = length($aaseq);
				#print "#",$aaseq,"\n";
				#print ">",substr($aaseq,$start,$end-$start),"\n";
				push(@aa, [$aaseq,$start,$end,$end-$start]);
				#print { $outfh{ 'inframe' } } '@',$seqo->display_id(),"\t*** NOT FOUND $model_cds_start! ***\n",$seq,"\n","+\n",$qual,"\n";
			}
		}
		
		if (scalar(@aa)) {
			my ($aaseq) = @{ [ sort { $b->[3] <=> $a->[3] } @aa ]->[0] };
			#print ">>>>",$aaseq,"\n";
			my ($scorer, $aln1r, $aln2r) = split(';', $seqan->Align2Seq($model_orf_start, $seq));
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
					$status='inframe';
					my $end = length($seq);
					
                                        if ($seq=~/$VECTOR_ORF_END_RE/g) {
						$end = undef;
						$end = pos($seq)-length($vector_orf_end);
						$status='phage';
						
                                                my ($score_vector_orf_inner) = split(';', $seqan->Align2Seq($vector_orf_inner, substr($seq,$start,$end-$start))); 
						if ($score_vector_orf_inner >= $minscore_vector_orf_inner) {
							$status='empty';
						}
						
                                        } else { 
                                               ($scorer, $aln1r, $aln2r) = split(';', $seqan->Align2Seq($vector_orf_end, $seq));
					       if ($scorer>30) {
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
							       $end = undef;
							       $end = pos($seq)-length($aln1r);
							       $status='phage';
	
								my ($score_vector_orf_inner) = split(';', $seqan->Align2Seq($vector_orf_inner, substr($seq,$start,$end-$start))); 
								if ($score_vector_orf_inner >= $minscore_vector_orf_inner) {
									$status='empty';
								}
							}
					       } 
					}
					
					if ($aaseq=~/\*/) {
						my $tmpseq = substr($seq,$start,$end-$start);
						#print $seqo->display_id(),"\t$aaseq\n";
						#print $seq,"\n";
						my $ncodon=1;
						foreach my $codon ( $tmpseq=~/.{3}/g ) {
							#print $ncodon,"\t",$codon,"\t",(($stop{$codon}) ? $stop{$codon} : ''),"\n";
							last if ($stop{$codon});
							$ncodon++;
						}
						$end=(($ncodon-1)*3)+$start;
						#print ">$start $end\n";
					}

					#print $seqo->display_id(),">>>$start\t$end\n";
					#print $seq,"\n";
					#print substr($seq,$start,$end-$start),"\n";
					print { $outfh{ $status } } '@',$seqo->display_id(),"\t$aaseq\n",substr($seq,$start,$end-$start),"\n","+\n",substr($qual,$start,$end-$start),"\n";
				} else {
					$LOGGER->logdie("NOT FOUND start !!!");
				}
			} else {
				$status='wrongdna';
				print { $outfh{ 'wrongdna' } } '@',$seqo->display_id(),"\t$aaseq\n",$seq,"\n","+\n",$qual,"\n";
			}
		} else {
			$status='notinframe';
			print { $outfh{ 'notinframe' } } '@',$seqo->display_id(),"\t*** NOT FOUND $model_cds_start ***\n",$seq,"\n","+\n",$qual,"\n";
		}

	}
	
	$LOGGER->info("Sequence ($id) status ($status)");
	
}
$progress->update($i);


foreach my $s (keys %outfh) {
	$outfh{$s}->close();
}

# Subroutines

sub Usage {
    my ($msg) = @_;
	
	$minseqlen = DEFAULT_MIN_SEQ_LENGTH;
	$model_orf_start=DEFAULT_MODEL_ORF_START;
	$vector_orf_end=DEFAULT_VECTOR_ORF_END;
	$vector_orf_inner=DEFAULT_VECTOR_ORF_INNER;
	$minscore_vector_orf_inner=DEFAULT_MIN_SCORE_VECTOR_ORF_INNER;
	
    my $USAGE = <<"END_USAGE";
Daniel Guariz Pinheiro (dgpinheiro\@gmail.com)
(c)2018 Universidade Estadual Paulista "Júlio de Mesquita Filho"

Usage

        $0	[-h/--help] [-l/--level <LEVEL>]

Argument(s)

	-h	--help			Help
	-l	--level			Log level [Default: FATAL]
	-i	--infile		Input fastq file
	-o	--outdir		Output directory
	-s	--sample		Sample name
	-x	--suffix		Suffix added to output file [Default .fastq]
	-m	--minseqlen		Minimum sequence length [Default: $minseqlen]
        -mos	--modorfstart		Model ORF start sequence (including start codon) [Default: $model_orf_start]
        -voe	--vectororfend		Vector ORF end sequence (including stop codon) [Default: $vector_orf_end]
        -voi	--vectororfinner	Vector ORF inner sequence (closed vector - sequence between -mos e -moe) [Default: $vector_orf_inner]
        -msv	--minscorevoi		Minimum alignment score to recognize closed vector [Default: $minscore_vector_orf_inner]

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
    print STDERR $USAGE;
    exit(1);
}

