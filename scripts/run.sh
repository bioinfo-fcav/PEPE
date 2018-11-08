#!/bin/bash

r1file=$1
r2file=$2

bcfile=$3
outdir=$4

outsample=$5

dmndfile=$6

genenamefile=$7
deflinefile=$8

if [ ! ${r1file} ] || [ ! ${r2file} ] || [ ! ${bcfile} ] || [ ! ${outdir} ] || [ ! ${outsample} ] || [ ! ${genenamefile} ] || [ ! ${deflinefile} ]; then
	echo "Please, put the arguments in this order: <R1 file path> <R2 file path> <barcode-sample file path> <output directory path> <outgroup sample name> <proteome diamond index (.dmnd) file path> <Phytozome gene Name file path> <Phytozome defline file path>"
	exit;
fi

if [ ! ${r1file} ]; then
        echo "Missing R1 file"
        exit
else   
        if [ ! -e ${r1file} ]; then
                echo "Wrong R1 file (${r1file})"
                exit
        else   
                if [[ ! ${r1file} =~ _R1.fastq$ ]]; then
                        if [[ ${r1file} =~ _R1.fastq.gz ]]; then
                                echo "Uncompressing ${r1file} ..."
                                gunzip ${r1file}
                                r1file=`echo ${r1file} | sed 's/.gz$//'`
                        else   
                                echo "Wrong name for R1 file (${r1file})"
                                exit
                     fi
                fi
        fi
fi

if [ ! ${r2file} ]; then
        echo "Missing R2 file"
        exit
else   
        if [ ! -e ${r2file} ]; then
                echo "Wrong R2 file (${r2file})"
                exit
        else   
                if [[ ! ${r2file} =~ _R2.fastq$ ]]; then
                        if [[ ${r2file} =~ _R2.fastq.gz ]]; then
                                echo "Uncompressing ${r2file} ..."
                                gunzip ${r2file}
                                r2file=`echo ${r2file} | sed 's/.gz$//'`
                        else   
                                echo "Wrong name for R2 file (${r2file})"
                                exit
                     fi
                fi
        fi
fi


if [ ! -e ${bcfile} ]; then
	echo "Wrong BC file (${bcfile})"
	exit
fi

if [ ! -d ${outdir} ]; then
	echo "Wrong output directory (${outdir})"
	exit
fi


if [ $(grep -c -P "\b${outsample}[0-9]*" ${bcfile}) == 0 ]; then
	echo "Not found ${outsample} in ${bcfile}"
	exit
fi

if [ ! ${dmndfile} ]; then
        echo "Missing diamond index file (.dmnd)"
        exit
else   
        if [ ! -e ${dmndfile} ]; then
                echo "Wrong diamond index file (${dmndfile})"
                exit
        else   
                if [[ ! ${dmndfile} =~ .dmnd$ ]]; then
			echo "Wrong diamond index file (${dmndfile})"
			exit
		fi
	fi
fi

if [ ! -e ${genenamefile} ]; then
	echo "Wrong Phytozome gene Name file (${genenamefile})"
	exit
fi

if [ ! -e ${deflinefile} ]; then
	echo "Wrong Phytozome defline file (${deflinefile})"
	exit
fi

./demux_PEP-Seq.pl -r1 ${r1file} -r2 ${r2file} -o ${outdir} -b ${bcfile}

./rcv_AmbIndex.pl -b ${bcfile} -o ${outdir} -a ${outdir}/AMBIGUOUS.assembled.fastq -s .assembled.fastq
./rcv_AmbIndex.pl -b ${bcfile} -o ${outdir} -a ${outdir}/AMBIGUOUS.unassembled.fastq -s .unassembled.fastq


declare -A GROUP
biogroups=()
counts_file_param=()

for s in $(cut -f 2 barcodes.txt); do
	
	echo "Processing sample $s ...";
	
	./check_InFrame.pl -i ${outdir}/${s}.assembled.fastq -o ${outdir}/ -s ${s} -x .assembled.fastq
	./check_InFrame.pl -i ${outdir}/${s}.unassembled.fastq -o ${outdir}/ -s ${s} -x .unassembled.fastq
	./check_InFrame.pl -i ${outdir}/${s}.recovered.assembled.fastq -o ${outdir}/ -s ${s} -x .recovered.assembled.fastq
	./check_InFrame.pl -i ${outdir}/${s}.recovered.unassembled.fastq -o ${outdir}/ -s ${s} -x .recovered.unassembled.fastq
	
	diamond blastx --more-sensitive --strand plus --top 10 --query-cover 80 --threads 20 --db ${dmndfile} \
	--query ./${outdir}/${s}.inframe.assembled.fastq --out ./${outdir}/${s}.inframe.assembled.diamond.daa --outfmt 100

	diamond blastx --more-sensitive --strand plus --top 10 --query-cover 80 --threads 20 --db ${dmndfile} \
	--query ./${outdir}/${s}.inframe.unassembled.fastq --out ./${outdir}/${s}.inframe.unassembled.diamond.daa --outfmt 100

	diamond view --daa ${outdir}/${s}.inframe.assembled.diamond.daa --outfmt 6 qseqid qlen sseqid bitscore qcovhsp qframe \
	--out ./${outdir}/${s}.inframe.assembled.diamond.txt

	diamond view --daa ${outdir}/${s}.inframe.unassembled.diamond.daa --outfmt 6 qseqid qlen sseqid bitscore qcovhsp qframe \
	--out ./${outdir}/${s}.inframe.unassembled.diamond.txt

	diamond blastx --more-sensitive --strand plus --top 10 --query-cover 80 --threads 20 --db ${dmndfile} \
	--query ./${outdir}/${s}.inframe.recovered.assembled.fastq --out ./${outdir}/${s}.inframe.recovered.assembled.diamond.daa --outfmt 100
	diamond blastx --more-sensitive --strand plus --top 10 --query-cover 80 --threads 20 --db ${dmndfile} --query ./${outdir}/${s}.inframe.recovered.unassembled.fastq --out ./${outdir}/${s}.inframe.recovered.unassembled.diamond.daa --outfmt 100
	diamond view --daa ${outdir}/${s}.inframe.recovered.assembled.diamond.daa --outfmt 6 qseqid qlen sseqid bitscore qcovhsp qframe --out ./${outdir}/${s}.inframe.recovered.assembled.diamond.txt
	diamond view --daa ${outdir}/${s}.inframe.recovered.unassembled.diamond.daa --outfmt 6 qseqid qlen sseqid bitscore qcovhsp qframe --out ./${outdir}/${s}.inframe.recovered.unassembled.diamond.txt

	./find_Protein.pl -s ${s} -o ${outdir}/ -i ${outdir}/${s}.inframe.assembled.diamond.txt -x .assembled.txt
	./find_Protein.pl -s ${s} -o ${outdir}/ -i ${outdir}/${s}.inframe.unassembled.diamond.txt -x .unassembled.txt
	./find_Protein.pl -s ${s} -o ${outdir}/ -i ${outdir}/${s}.inframe.recovered.assembled.diamond.txt -x .recovered.assembled.txt
	./find_Protein.pl -s ${s} -o ${outdir}/ -i ${outdir}/${s}.inframe.recovered.unassembled.diamond.txt -x .recovered.unassembled.txt
	
	./count_Reads.pl -i ${outdir} -s ${s} -x *id*assembled.txt -o ${outdir}/${s}.counts.txt -l INFO

	counts_file_param=(${counts_file_param[@]} "-d ${s}=${outdir}/${s}.counts.txt")
	bgroup=$(echo ${s} | sed 's/[0-9]\+$//')
	biogroups=($( printf "%s\n" ${biogroups[@]} ${bgroup} | sort -u ))
	GROUP[${bgroup}]="${GROUP[${bgroup}]} ${s}"
done

groups_param=()
for K in "${biogroups[@]}"; do
	groups_param=(${groups_param[@]} "-g ${K}=$(echo ${GROUP[${K}]} | sed 's/ /,/')")
done

./join_Results.pl ${counts_file_param[*]} ${groups_param[*]} > ${outdir}/ReadCountsMatrix.txt

./filter_Results.pl ${groups_param} -i ${outdir}/ReadCountsMatrix.txt -o ${refsample} -u 2 -s 3 -m 1 > ${outdir}/FilteredReadCountsMatrix_2_3_1.txt

./annot_Results.pl -i ${outdir}/FilteredReadCountsMatrix_2_3_1.txt -g ${genenamefile} -d ${deflinefile} > ${outdir}/AnnotFilteredReadCountsMatrix_2_3_1.txt

