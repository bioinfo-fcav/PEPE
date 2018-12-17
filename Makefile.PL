use 5.022001;
use ExtUtils::MakeMaker;

my $CC = 'g++';

# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    NAME              => 'PEPE::SeqAn',
    VERSION_FROM      => 'lib/PEPE/SeqAn.pm', # finds $VERSION, requires EU::MM from perl >= 5.5
    PREREQ_PM         => {	'Term::ProgressBar'=>'0',
				'Log::Log4perl'=>'0',
				'Bio::Tools::CodonTable'=>'0',
				'Bio::Seq'=>'0',
				'Bio::SeqIO'=>'0'
			}, # e.g., Module::Name => 1.1
    ABSTRACT_FROM     => 'lib/PEPE/SeqAn.pm', # retrieve abstract from module
    AUTHOR            => 'Daniel Guariz Pinheiro <dgpinheiro@gmail.com>',
    #LICENSE           => 'perl',
    #Value must be from legacy list of licenses here
    #http://search.cpan.org/perldoc?Module%3A%3ABuild%3A%3AAPI
    LIBS              => [''], # e.g., '-lm'
    DEFINE            => '', # e.g., '-DHAVE_SOMETHING'
    INC               => '-I.', # e.g., '-I. -I/usr/include/other'
    CC                => $CC,
    LD                => '$(CC)',
    XSOPT             => '-C++',
    TYPEMAPS          => ['perlobject.map']
	# Un-comment this if you add C files to link with later:
    # OBJECT            => '$(O_FILES)', # link all the C files too
);
