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
use FileHandle;

INIT {
    use Log::Log4perl qw/:easy/;
    Log::Log4perl->easy_init($FATAL);
    $LOGGER = Log::Log4perl->get_logger($0);
}

my ($level, $indir, $sample, $suffix, $outfile, $countprotein, $finalfile, $anninfofile, $deflinefile, $bcfile);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
	    "i|indir=s"=>\$indir,
            "b|bcfile=s"=>\$bcfile,
	    "o|outfile=s"=>\$outfile,
	    "f|final=s"=>\$finalfile,
            "a|anninfo=s"=>\$anninfofile,
            "d|defline=s"=>\$deflinefile
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

$LOGGER->logdie("Missing input directory") if (!defined $indir);
$LOGGER->logdie("Wrong input directory ($indir)") if (! -e $indir);

$LOGGER->logdie("Missing output file") if (!defined $outfile);

$LOGGER->logdie("Missing final output file") if (!defined $finalfile);
$LOGGER->logdie("Wrong final output file ($finalfile)") if (! -e $finalfile);

$LOGGER->logdie("Missing barcodes file") if (!defined $bcfile);
$LOGGER->logdie("Wrong barcodes file ($bcfile)") if (! -e $bcfile);

$LOGGER->logdie("Missing annotation_info file") if (!defined $anninfofile);
$LOGGER->logdie("Wrong annotation_info file ($anninfofile)") if (! -e $anninfofile);

$LOGGER->logdie("Missing defline file") if (!defined $deflinefile);
$LOGGER->logdie("Wrong defline file ($deflinefile)") if (! -e $deflinefile);

open(BC, "<", $bcfile) or $LOGGER->logdie($!);
my %biosample;
while(<BC>) {
	chomp;
	my ($bcseq, $sample_name, $sample_group) = split(/\t/, $_);
	$sample_name=~s/^\s*//;
	$sample_name=~s/\s*$//;
	$biosample{$sample_group}->{$sample_name} = $bcseq;
}

close(BC);

my %outfh;

$outfh{'counts'} = FileHandle->new(">".$outfile);
autoflush { $outfh{'counts'} } 1;
$LOGGER->logdie($!) if (! defined $outfh{'counts'});

my %fastqdata;
my %fastqheader;

my @biosample_order = map {  sort keys %{ $biosample{$_} }   } sort keys %biosample;

#print join(";", @biosample_order),"\n";

my @ftype=(
'assembled',
'unassembled',
'recovered.assembled',
'recovered.unassembled',
'inframe.assembled',
'inframe.unassembled',
'phage.assembled',
'phage.unassembled',
'decontaminated.assembled',
'decontaminated.unassembled',
'smallseq.assembled',
'smallseq.unassembled',
'empty.assembled',
'empty.unassembled',
'notinframe.assembled',
'notinframe.unassembled',
'recovered.assembled',
'recovered.unassembled',
'inframe.recovered.assembled',
'inframe.recovered.unassembled',
'phage.recovered.assembled',
'phage.recovered.unassembled',
'decontaminated.recovered.assembled',
'decontaminated.recovered.unassembled',
'smallseq.recovered.assembled',
'smallseq.recovered.unassembled',
'empty.recovered.assembled',
'empty.recovered.unassembled',
'notinframe.recovered.assembled',
'notinframe.recovered.unassembled'
);

my @bsft_order;
foreach my $bs (@biosample_order) {
	foreach my $ft (@ftype) {
		#print  $bs.'.'.$ft,"\n";
		push(@bsft_order, $bs.'.'.$ft);
	}
}
my @extra = (	'atropos_final',
		'REFUSED',
		'WRONGTEMPLATE',
		'UNIDENTIFIED',
		'AMBIGUOUS',
		);

my @e_order;
foreach my $e (@extra) {
	foreach my $ft ('assembled','unassembled') {
		push(@e_order, $e.'.'.$ft);
	}
}


my %info;
my @info_order = (	'INPUT_READ_PAIRS',
			'atropos_final',
			@e_order,
			@bsft_order
		);

my %disregarded;
@disregarded{	
		'atropos_final.unassembled.forward.blastx',
		'atropos_final.unassembled.forward',
		'atropos_final.unassembled.reverse',
		'atropos_final.unassembled.reverse.blastx' } = ();


@info{@info_order} = ();

open(IN, "<", "$indir/atropos/atropos_insert.log.out.txt") or $LOGGER->logdie($!);
while(<IN>) {
	chomp;
	if ($_=~/^Total read pairs processed:\s+(\S+)/) {
		$info{"INPUT_READ_PAIRS"}=$1;
		last;
	}
}
close(IN);

unless ($info{"INPUT_READ_PAIRS"}) {
	$LOGGER->logdie("Wrong execution of ATROPOS!");
} else {
	$info{"INPUT_READ_PAIRS"}=~s/,//g;
	$info{"INPUT_READ_PAIRS"}+=0;
} 


foreach my $idfile (glob("$indir/atropos/*_R1*.atropos_final.fastq"), glob("$indir/*.fastq")) {
	my $bnfile = basename($idfile,'.fastq');
	$LOGGER->info("Processing FASTQ $bnfile ...");

	open(IN, "<", $idfile) or $LOGGER->logdie($!);
	while(<IN>) {
		chomp;
		#print $.,"\t".($.%4),"\n";
		if ( ($.%4)==1 ) {
			#print $_,"\n";
			$fastqdata{ $bnfile } = 0 unless (exists $fastqdata{ $bnfile });
			$fastqdata{ $bnfile }++;
		}
	}
	close(IN);
}


foreach my $bnfile (keys %fastqdata) {
	my ($sample, $suffix)=$bnfile=~/^([^\.]+)\.?(.*)/;
	$suffix||='';
	if (exists $info{$bnfile}) {
		$info{$bnfile} = $fastqdata{$bnfile};
	} elsif (exists $info{$suffix}) {
		$info{$suffix} = $fastqdata{$bnfile};
	} elsif (exists $disregarded{$suffix}) {
		# Nothing to do
	} else {
		$LOGGER->logdie("Not recognized: $bnfile");
	}
	
}

my @analysis_order;
foreach my $tp ('Nreads', 'Ngenes', 'Nproteins') {
	foreach my $cmp ('id.inframe', 'notid.inframe', 'id.decontaminated', 'notid.decontaminated', 'id_wrongframe.inframe') {
		foreach my $asm ('assembled', 'recovered.assembled', 'unassembled', 'recovered.unassembled') {
			push(@analysis_order, $tp.'.'.$cmp.'.'.$asm);
		}
	}
}

foreach my $id (@biosample_order) {
	foreach my $an (@analysis_order) {
		$info{ $id.'.'.$an } = undef;
	}
}


my %bedheader;
my %beddata;
foreach my $bedfile (glob("$indir/*.bed")) {
	my $bnfile = basename($bedfile,'.bed');
	my ($sample, $suffix)=$bnfile=~/^([^\.]+)\.?(.*)/;
	$suffix||='';
	$LOGGER->info("Processing BED [$sample] $bnfile ...");

	my %analysis;
	my %gene;
	my %prot;
	my %read;
	open(IN, "<", $bedfile) or $LOGGER->logdie($!);
	while(<IN>) {
		chomp;
		my @F=split(/\t/, $_);
		$analysis{'Nhits'} = 0 unless (exists $analysis{'Nhits'});
		$analysis{'Nhits'}++;
		$prot{$F[0]}=undef;
		my ($g) = $F[0]=~/^([^\.]+)/;
		$gene{$g}=undef;
		$read{$F[3]}=undef;
	}
	close(IN);

	$analysis{'Nproteins'} = scalar(keys %prot);
	$analysis{'Ngenes'} = scalar(keys %gene);
	$analysis{'Nreads'} = scalar(keys %read);
	foreach my $anls (keys %analysis) {

		if (exists $info{ $sample.'.'.$anls.'.'.$suffix }) {
			$info{ $sample.'.'.$anls.'.'.$suffix } = $analysis{$anls}||0;
		}
		
		$beddata{ $sample }->{ $anls.'.'.$suffix } = 0 unless (exists $beddata{ $sample }->{ $anls.'.'.$suffix });
		$beddata{ $sample }->{ $anls.'.'.$suffix } = $analysis{$anls};
		#print ">>>>",$sample,"\t",$anls.'.'.$suffix,"\t",$analysis{$anls},"\n";
		$bedheader{ $anls.'.'.$suffix } = undef;
	}
}


my @fastqh = (sort { $a cmp $b } keys %fastqheader);
my @bedh = (sort { $a cmp $b } keys %bedheader);

my %finalprot;
open(IN, "<", $finalfile) or $LOGGER->logdie($!);
while(<IN>) {
	chomp;
	my @F = split(/\t/, $_);
	$finalprot{$F[0]} = undef;
}
close(IN);

my %defline;
{
	my %check_protein = %finalprot;
	open(IN, "<", $deflinefile) or $LOGGER->logdie($!);
	while(<IN>) {
		chomp;
		my @F=split(/\t/, $_);
		if (exists $finalprot{ $F[0] } ) {
			delete($check_protein{ $F[0] });
			$finalprot{ $F[0] }=0;
			$finalprot{ $F[0] }=1 if ($F[1]=~/hypothetical protein/);
		}
	}
	close(IN);

	if (scalar(keys %check_protein)>0) {
		$LOGGER->logdie("Not found all proteins in defline file ($deflinefile)");
	}
}

my %finalgene;
{
	my %check_protein = %finalprot;

	open(IN, "<", $anninfofile) or $LOGGER->logdie($!);
	my $ann_info_header = <IN>;
	chomp($ann_info_header);
	my @ann_info = split(/\t/, $ann_info_header);
	while(<IN>) {
		chomp;
		my %ann;
		@ann{ @ann_info } = split(/\t/, $_);

		if (exists $finalprot{ $ann{'peptideName'}  }) {
			delete($check_protein{ $ann{'peptideName'} });
			unless (exists $finalgene{ $ann{'locusName'} }) {
				$finalgene{ $ann{'locusName'} }=0;
			}
			$finalgene{ $ann{'locusName'} }=1 if ($finalprot{ $ann{'peptideName'}  });
		}
	}
	close(IN);
	
	if (scalar(keys %check_protein)>0) {
		$LOGGER->logdie("Not found all proteins in annotation_info file ($anninfofile)");
	}
}

print { $outfh{'counts'} } '#','Input read pairs .........................................: ', $info{ 'INPUT_READ_PAIRS' },"\n";
print { $outfh{'counts'} } '#','Processed read pairs (atropos_final) .....................: ', $info{ 'atropos_final' },"\n";
print { $outfh{'counts'} } '#','Number of genes passing filter............................: '.scalar(keys %finalgene)."\n";
print { $outfh{'counts'} } '#','Number of genes passing filter ("hypothetical protein")...: '.&count(\%finalgene)."\n";
print { $outfh{'counts'} } '#','Number of proteins passing filter.........................: '.scalar(keys %finalprot)."\n";
print { $outfh{'counts'} } '#','Number of proteins passing filter ("hypothetical protein"): '.&count(\%finalprot)."\n";
print { $outfh{'counts'} } '#',join("\t", 'ID', @ftype, @analysis_order),"\n";
foreach my $id (@extra, @biosample_order) {
	print { $outfh{'counts'} } $id;
	foreach my $ft (@ftype, @analysis_order) {
		print { $outfh{'counts'} } "\t",($info{ $id.'.'.$ft }||0);
	}
	print { $outfh{'counts'} } "\n";
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
	-i	--indir		Input directory
	-f	--final		Final results peaks file ( _peaks.txt )
	-a	--anninfo	Protein annotation information (.annotation_info.txt)
	-d	--defline	Protein definition file (.defline.txt)
	-o	--outfile	Output file

END_USAGE
    print STDERR "\nERR: $msg\n\n" if $msg;
    print STDERR qq[$0  ] . q[$Revision$] . qq[\n];
    print STDERR $USAGE;
    exit(1);
}

sub count {
	my ($hr_g) = @_;
	my $r=0;
	if ($hr_g) {
		foreach my $k (keys %{ $hr_g }) {
			if ($hr_g->{$k}) {
				$r++;
			}
		}
	}
	return $r;
}

sub sum {
	my ($hr_g) = @_;
	my $r=0;
	if ($hr_g) {
		foreach my $k (keys %{ $hr_g }) {
			if ($hr_g->{$k}) {
				$r+=$hr_g->{$k};
			}
		}
	}
	return $r;
}
