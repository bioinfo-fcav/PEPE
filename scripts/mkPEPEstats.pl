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

my ($level, $indir, $sample, $suffix, $outfile, $countprotein, $finalfile, $protdef);

Usage("Too few arguments") if $#ARGV < 0;
GetOptions( "h|?|help" => sub { &Usage(); },
            "l|level=s"=> \$level,
	    "i|indir=s"=>\$indir,
	    "o|outfile=s"=>\$outfile,
	    "f|final=s"=>\$finalfile,
            "p|protdef=s"=>\$protdef
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

my %outfh;

$outfh{'counts'} = FileHandle->new(">".$outfile);
autoflush { $outfh{'counts'} } 1;
$LOGGER->logdie($!) if (! defined $outfh{'counts'});

my %fastqdata;
my %fastqheader;
foreach my $idfile (glob("$indir/*.fastq")) {
	my $bnfile = basename($idfile,'.fastq');
	my ($sample, $suffix)=$bnfile=~/^([^\.]+)\.?(.*)/;
	$suffix||='';
	$LOGGER->info("Processing FASTQ [$sample] $bnfile ...");
	open(IN, "<", $idfile) or $LOGGER->logdie($!);
	while(<IN>) {
		chomp;
		#print $.,"\t".($.%4),"\n";
		if ( ($.%4)==1 ) {
			#print $_,"\n";
			$fastqdata{ $sample }->{ $suffix } = 0 unless (exists $fastqdata{ $sample }->{ $suffix });
			$fastqdata{ $sample }->{ $suffix }++;
			$fastqheader{ $suffix } = undef;
		}
	}
	close(IN);
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
		my ($g) = $F[0]=~s/\.\d+$//;
		$gene{$g}=undef;
		$read{$F[3]}=undef;
	}
	close(IN);
	$analysis{'Nproteins'} = scalar(keys %prot);
	$analysis{'Ngenes'} = scalar(keys %gene);
	$analysis{'Nreads'} = scalar(keys %read);
	foreach my $analysis (keys %analysis) {
		$beddata{ $sample }->{ $analysis.'.'.$suffix } = 0 unless (exists $beddata{ $sample }->{ $analysis.'.'.$suffix });
		$beddata{ $sample }->{ $analysis.'.'.$suffix }++;
		$bedheader{ $analysis.'.'.$suffix } = undef;
	}
}

my @fastqh = (sort { $a cmp $b } keys %fastqheader);
my @bedh = (sort { $a cmp $b } keys %bedheader);


my %finalgene;
open(IN, "<", $finalfile) or $LOGGER->logdie($!);
while(<IN>) {
	chomp;
	next if ($.==1);
	my @F = split(/\t/, $_);
	$finalgene{$F[0]}=0;
	$finalgene{$F[0]}=1 if ($F[2]=~/hypothetical protein/);
}
close(IN);

my %finalprot;
open(IN, "<", $protdef) or $LOGGER->logdie($!);
while(<IN>) {
	chomp;
	my @F = split(/\t/, $_);
	my ($g) = $F[0]=~/^([^\.]+)/;
	if (exists $finalgene{$g}) {
		$finalprot{$F[0]}=0;
		$finalprot{$F[0]}=1 if ($finalgene{$g});
	}
}
close(IN);

print { $outfh{'counts'} } '# Number of genes passing filter: '.scalar(keys %finalgene)."\n";
print { $outfh{'counts'} } '# Number of genes passing filter ("hypothetical protein"): '.&count(\%finalgene)."\n";
print { $outfh{'counts'} } '# Number of proteins passing filter: '.scalar(keys %finalprot)."\n";
print { $outfh{'counts'} } '# Number of proteins passing filter ("hypothetical protein"): '.&count(\%finalprot)."\n";
print { $outfh{'counts'} } join("\t", '#ID', @fastqh, @bedh),"\n";
foreach my $sample (sort { $a cmp $b } (keys %fastqdata)) {	
	print { $outfh{'counts'} } join("\t", $sample, map {$_||0} (@{$fastqdata{$sample}}{@fastqh}, @{$beddata{$sample}}{@bedh})),"\n";
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
	-f	--final		Final results file
	-p	--protdef	Protein definition file
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
