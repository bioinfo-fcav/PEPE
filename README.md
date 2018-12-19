# PEPE

## PEP-Seq Explorer

### Dependencies:

	Tools and libraries

		DIAMOND - alignment tool for short DNA sequencing reads [v.0.9.22] [https://ab.inf.uni-tuebingen.de/software/diamond]

		PEAR - Paired-End reAd mergeR [v.0.9.8] [https://sco.h-its.org/exelixis/web/software/pear/ ]

		SeqAn - The Library for Sequence Analysis [v1.4.2] [https://github.com/seqan/seqan]

		Perl libraries

			Term::ProgressBar	[v.2.22] [https://metacpan.org/release/Term-ProgressBar]

			Log::Log4perl		[v.1.49] [https://metacpan.org/release/Log-Log4perl]

			BioPerl			[v.1.007002] [https://metacpan.org/release/BioPerl]

	Files:
	
		Athaliana_447_Araport11.protein.fa
		Athaliana_447_Araport11.geneName.txt
		Athaliana_447_Araport11.defline.txt
	

### Usage:

 Before you run the following command, please, format the Arabidopsis proteome file (Athaliana_447_Araport11.protein.fa) with DIAMOND build to
get the index (Athaliana_447_Araport11.dmnd).

```
 cd ./scripts/

 ./run.sh <R1 file path> <R2 file path> <barcode-sample file path> <output directory path> <outgroup sample name> <proteome diamond index (.dmnd) file path> <Phytozome gene Name file path> <Phytozome defline file path>
```

### Example:

```
 time ./run.sh ../example/data/INPUT_R1.fastq ../example/data/INPUT_R2.fastq ../example/barcodes.txt ../example/output/ BSA ../../Atha/refs/DIAMOND/Athaliana_447_Araport11.dmnd ../../Atha/refs/Athaliana_447_Araport11.geneName.txt ../../Atha/refs/Athaliana_447_Araport11.defline.txt
```

### Flowchart diagram:

![Flowchart diagram](https://raw.githubusercontent.com/bioinfo-fcav/PEPE/master/DiagramPEP-Seq-v4.png)

