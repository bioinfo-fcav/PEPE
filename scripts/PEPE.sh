#!/bin/bash

configfile=${1}

if [ ! ${configfile} ]; then
	echo "Missing configuration file"
	exit
fi

if [ ! -e ${configfile} ]; then
	echo "Not found configuration file (${configfile})"
	exit
else
	dos2unix ${configfile}
fi

echo "Using ${configfile} as configuration file!"

configfile_secured=$(mktemp)

# checagem se o arquivo contem algo indesejavel (que não segue o padrao esperado para arquivos de configuracao chave=valor)
if egrep -q -v '^[^ ]*=[^;]*' "${configfile}"; then
  echo "Config file is unclean, cleaning it..." >&2
  # filtragem do arquivo original para um novo arquivo
  egrep '^[^ ]*=[^;&]*'  "${configfile}" > "${configfile_secured}"
  configfile="${configfile_secured}"
fi

# carrega o arquivo original ou sua versao filtrada
source "${configfile}"

# remove arquivo de configuração temporario
rm ${configfile_secured}

scriptsdir=`dirname $0`;

if [ ${2} ]; then
	# overwrite config file num_threads
	num_threads=${2}
fi

echo "Executing PEPE with ${num_threads} threads ..."


if [ ! ${r1file} ]; then
        echo "Missing R1 file"
        exit
else   
        if [ ! -e ${r1file} ]; then
                echo "Wrong R1 file (${r1file})"
                exit
        else   
                if [[ ! ${r1file} =~ _R1(_[0-9]+)?.fastq$ ]]; then
                        if [[ ${r1file} =~ _R1(_[0-9]+)?.fastq.gz ]]; then
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
                if [[ ! ${r2file} =~ _R2(_[0-9]+)?.fastq$ ]]; then
                        if [[ ${r2file} =~ _R2(_[0-9]+)?.fastq.gz ]]; then
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
	echo "Wrong output directory (${outdir}). Please create the output directory before running this script ($0)."
	exit
fi


if [ $(grep -c -P "\b${outgroup}$" ${bcfile}) == 0 ]; then
	echo "Not found ${outgroup} in ${bcfile}"
	exit
fi

declare -A TOEXCLUDE
if [ ${excludegroups} ]; then
	IFS="," read -ra NAMES <<< "${excludegroups}"; 

	for egroup in ${NAMES[@]}; do 
		if [ $(grep -c -P "\b${egroup}$" ${bcfile}) == 0 ]; then
			echo "Not found ${egroup} in ${bcfile} to exclude from analysis"
			exit
		else
			TOEXCLUDE[${egroup}]=1
		fi
	done
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


bnr1=`basename ${r1file} .fastq`
bnr2=`basename ${r2file} .fastq`


if [ ! -e "${outdir}/atropos/${bnr1}.atropos_final.fastq" ] && 
   [ ! -e "${outdir}/atropos/${bnr2}.atropos_final.fastq" ]; then


	echo "Trimming step ..."
	
	if [ ${atropos_pyenv} ]; then
		eval "$(pyenv init -)"
		eval "$(pyenv virtualenv-init -)"
		pyenv activate ${atropos_pyenv}
	fi
	
	mkdir -p ${outdir}/atropos
	
	CHECK_ATROPOS=$(atropos --version 2>&1)
	if [ ! ${CHECK_ATROPOS}=~1 ]; then
		echo "Not found atropos ($(atropos))"
		rm -fr ${outdir}/atropos/
		exit
	fi

	if [ ! -e "${outdir}/atropos/${bnr1}.atropos_insert.fastq" ] && 
	   [ ! -e "${outdir}/atropos/${bnr2}.atropos_insert.fastq" ]; then

		atropos trim --aligner insert \
			     -e ${atropos_insert_e} \
			     -n ${atropos_insert_n} \
			     -m ${atropos_insert_m} \
			     --op-order ${atropos_insert_op_order} \
			     --match-read-wildcards \
			     -O ${atropos_insert_O} \
			     -q ${atropos_insert_q} \
			     -T ${num_threads} \
			     --correct-mismatches conservative \
			     --pair-filter any \
			     -a ${atropos_a} \
			     -A ${atropos_A} \
			     -o ${outdir}/atropos/${bnr1}.atropos_insert.fastq \
			     -p ${outdir}/atropos/${bnr2}.atropos_insert.fastq \
			     -pe1 ${r1file} \
			     -pe2 ${r2file} \
			     --untrimmed-output        ${outdir}/atropos/${bnr1}.atropos_untrimmed.fastq \
			     --untrimmed-paired-output ${outdir}/atropos/${bnr2}.atropos_untrimmed.fastq \
			       > ${outdir}/atropos/atropos_insert.log.out.txt \
			      2> ${outdir}/atropos/atropos_insert.log.err.txt
	else
		echo "::: Atropos (insert mode) step appears to have been performed!"
	fi
	
	if 	[ -e "${outdir}/atropos/${bnr1}.atropos_untrimmed.fastq" ] && 
		[ -e "${outdir}/atropos/${bnr2}.atropos_untrimmed.fastq" ]; then
	
		# o parâmetro --pair-filter any é obrigatório pois caso contrário o pear apresentará problemas
		atropos trim    --aligner adapter \
				-e ${atropos_adapter_e} \
				-n ${atropos_adapter_n} \
				-m ${atropos_adapter_m} \
		     		--op-order ${atropos_adapter_op_order} \
				--match-read-wildcards  \
				-O ${atropos_adapter_O} \
				-q ${atropos_adapter_q} \
				--pair-filter any \
				-pe1 ${outdir}/atropos/${bnr1}.atropos_untrimmed.fastq \
				-pe2 ${outdir}/atropos/${bnr2}.atropos_untrimmed.fastq \
				-a ${atropos_a} \
				-A ${atropos_A} \
				-g ${atropos_g} \
				-G ${atropos_G} \
				-T ${num_threads} \
				-o  ${outdir}/atropos/${bnr1}.atropos_adapter.fastq \
				-p  ${outdir}/atropos/${bnr2}.atropos_adapter.fastq \
				 > ${outdir}/atropos/atropos_adapter.log.out.txt \
				2> ${outdir}/atropos/atropos_adapter.log.err.txt

	else
		touch ${outdir}/atropos/${bnr1}.atropos_adapter.fastq
		touch ${outdir}/atropos/${bnr2}.atropos_adapter.fastq
	fi

	cat     ${outdir}/atropos/${bnr1}.atropos_insert.fastq \
		${outdir}/atropos/${bnr1}.atropos_adapter.fastq \
	   > ${outdir}/atropos/${bnr1}.atropos_final.fastq

	cat     ${outdir}/atropos/${bnr2}.atropos_insert.fastq \
		${outdir}/atropos/${bnr2}.atropos_adapter.fastq \
	   > ${outdir}/atropos/${bnr2}.atropos_final.fastq

	rm -f ${outdir}/atropos/${bnr1}.atropos_insert.fastq ${outdir}/atropos/${bnr1}.atropos_adapter.fastq
	rm -f ${outdir}/atropos/${bnr2}.atropos_insert.fastq ${outdir}/atropos/${bnr2}.atropos_adapter.fastq

	if [ ${atropos_pyenv} ]; then
		pyenv deactivate atropos
	fi

fi

infilebn=`basename ${outdir}/atropos/${bnr1}.atropos_final .fastq | sed  's/_R1\(_\?\)/\\1/'`

if 	[ ! -e ${outdir}/${infilebn}.assembled.fastq ] && 
	[ ! -e ${outdir}/${infilebn}.unassembled.forward.fastq ] && 
	[ ! -e ${outdir}/${infilebn}.unassembled.reverse.fastq ]; then

	echo "::: Executing Demultiplex step ..."

	if [ ${demux_ooe} ]; then
		${scriptsdir}/demux_PEP-Seq.pl 	-m ${demux_reads_merger} \
						-y ${demux_reads_merger_memory} \
						-mo ${demux_reads_merger_min_overlap} \
						-mm ${demux_maxmismatches} \
						-qt ${demux_quality_threshold} \
						-ml ${minlength} \
						-md ${demux_maxdiffs:="0"} \
						-ms ${demux_minalnscore} \
						-m5p ${model_adapter_sequence_5prime} \
						-m3p ${model_adapter_sequence_3prime} \
						-ooe \
						-r1 ${outdir}/atropos/${bnr1}.atropos_final.fastq \
						-r2 ${outdir}/atropos/${bnr2}.atropos_final.fastq \
						-o ${outdir} -b ${bcfile} -t ${num_threads} 
	else
		if [ ${demux_nbases} ]; then
			${scriptsdir}/demux_PEP-Seq.pl 	-m ${demux_reads_merger} \
							-y ${demux_reads_merger_memory} \
							-mm ${demux_maxmismatches} \
							-ml ${minlength} \
							-ms ${demux_minalnscore} \
							-mm ${demux_maxmismatches} \
							-nb ${demux_nbases} \
							-m5p ${model_adapter_sequence_5prime} \
							-m3p ${model_adapter_sequence_3prime} \
							-r1 ${outdir}/atropos/${bnr1}.atropos_final.fastq \
							-r2 ${outdir}/atropos/${bnr2}.atropos_final.fastq \
							-o ${outdir} -b ${bcfile} -t ${num_threads}
			
		else
			${scriptsdir}/demux_PEP-Seq.pl 	-m ${demux_reads_merger} \
							-y ${demux_reads_merger_memory} \
							-d ${dmndfile} \
							-mm ${demux_maxmismatches} \
							-ml ${minlength} \
							-ms ${demux_minalnscore} \
							-mm ${demux_maxmismatches} \
							-m5p ${model_adapter_sequence_5prime} \
							-m3p ${model_adapter_sequence_3prime} \
							-tfr ${demux_trim_fwd_for_blastx} \
							-trr ${demux_trim_rev_for_blastx} \
							-dqc ${demux_diamond_qcover} \
							-dev ${demux_diamond_evalue} \
							-did ${demux_diamond_id} \
							-r1 ${outdir}/atropos/${bnr1}.atropos_final.fastq \
							-r2 ${outdir}/atropos/${bnr2}.atropos_final.fastq \
							-o ${outdir} -b ${bcfile} -t ${num_threads}

		fi
	fi

else
	echo "::: Demultiplex step appears to have been performed!"
fi


if [ ${recovery_step} ]; then


        recovered_file=()

        for rf in `ls ${outdir}/${s}.recovered*.fastq 2>/dev/null`; do
                recovered_file=(${recovered_file[@]} ${rf})
        done

	if [ ${#recovered_file[@]} -eq 0 ]; then

		${scriptsdir}/rcv_AmbIndex.pl 	-b ${bcfile} \
						-o ${outdir} \
						-a ${outdir}/AMBIGUOUS.assembled.fastq \
						-m3p ${model_adapter_sequence_3prime} \
						-s .assembled.fastq

		${scriptsdir}/rcv_AmbIndex.pl	-b ${bcfile} \
						-o ${outdir} \
						-a ${outdir}/AMBIGUOUS.unassembled.fastq \
						-m3p ${model_adapter_sequence_3prime} \
						-s .unassembled.fastq
	else
		echo "::: Recovery step appears to have been performed!"
	fi
else
	echo "::: Recovery step will not be performed!"
fi

declare -A GROUP
counts_file_param=()

for s in $(cut -f 2 ${bcfile}); do
	
	echo "Processing sample $s ...";
	ftype=("assembled" "recovered.assembled" "unassembled" "recovered.unassembled") 
	for ft in ${ftype[@]}; do

		ft_file=()
	        for ftf in `ls ${outdir}/${s}.${ft}.fastq 2>/dev/null`; do
        	        ft_file=(${ft_file[@]} ${ftf})
	        done
		
		if [ ${#ft_file[@]} -gt 0 ]; then

			echo -e "\t>>> Considering ${ft} dataset:"
			
			inframe_file=()

			for iff in `ls ${outdir}/${s}.inframe.${ft}.fastq 2>/dev/null`; do
				inframe_file=(${inframe_file[@]} ${iff})
			done

			if [ ${#inframe_file[@]} -eq 0 ]; then
				echo "Searching for in-frame translation sequences..."
				${scriptsdir}/check_InFrame.pl	-i ${outdir}/${s}.${ft}.fastq \
								-m ${minlength} \
 	      							-mos ${model_orf_start} \
							        -voe ${vector_orf_end} \
       	 							-voi ${vector_orf_inner} \
 								-msv ${minscore_vector_orf_inner} \
								-o ${outdir} \
								-s ${s} \
								-x .${ft}.fastq 
			else
				echo -e "\t::: Check in-frame translation step appears to have been performed!"
			fi
			
			if [ ${vectorfile} ]; then
				if [ ! -e "${vectorfile}" ]; then
					echo "Wrong vector file (${vectorfile})"
					exit;
				fi
				
				decontaminated_file=()

				for dcf in `ls ${outdir}/${s}.decontaminated.${ft}.fastq 2>/dev/null`; do
					decontaminated_file=(${decontaminated_file[@]} ${dcf})
				done
				
				if [ ${#decontaminated_file[@]} -eq 0 ]; then
					echo -e "\tDecontamination step ..."
					${scriptsdir}/rcv_Contained.pl 	-v ${vectorfile} \
									-i ${outdir}/${s}.phage.${ft}.fastq \
									-m   ${minalnlen} \
									-sma ${score_for_match} \
									-smi ${score_for_mismatch} \
									-sgo ${score_for_gap_open} \
									-sge ${score_for_gap_extension} \
									-o ${outdir} -s ${s} -x .${ft}.fastq
				else
					echo -e "\t::: Decontamination step appears to have been performed!"
				fi
			fi
			
			
			alnlist=("inframe" "decontaminated") 
			for al in ${alnlist[@]}; do
				if [ -e "${outdir}/${s}.${al}.${ft}.fastq" ]; then
					echo -e "\tFound ${outdir}/${s}.${al}.${ft}.fastq"
					if [ ! -e "${outdir}/${s}.${al}.${ft}.diamond.daa" ]; then
						echo -e "\t\tdiamond blastx againts proteome ..."
						diamond blastx	--more-sensitive \
								--evalue ${diamond_evalue} \
								--strand plus \
								--top 1 \
								--id ${diamond_id} \
								--query-cover ${diamond_qcover} \
								--threads ${num_threads} \
								--db ${dmndfile} \
								--query ${outdir}/${s}.${al}.${ft}.fastq \
								--out ${outdir}/${s}.${al}.${ft}.diamond.daa \
								--outfmt 100 \
								2> /dev/null > ${outdir}/${s}.${al}.${ft}.diamond.log.txt
					fi
					
					if [ ! -e "${outdir}/${s}.${al}.${ft}.diamond.txt" ]; then
						echo -e "\t\tdiamond view"
						diamond view	--daa ${outdir}/${s}.${al}.${ft}.diamond.daa \
								--outfmt 6 qseqid qlen sseqid bitscore qcovhsp qframe sstart send \
								--out ${outdir}/${s}.${al}.${ft}.diamond.txt \
								2> /dev/null > ${outdir}/${s}.${al}.${ft}.diamond.daa.log.txt
					fi

					if [ ! -e "${outdir}/${s}.id.${al}.${ft}.bed" ]; then
						echo -e "\tSearching for valid protein matches (${outdir}/${s}.${al}.${ft}.bed) ..."
						${scriptsdir}/find_Protein.pl	-s ${s} \
										-o ${outdir}/ \
										-i ${outdir}/${s}.${al}.${ft}.diamond.txt \
										-x .${al}.${ft}.bed \
										2> ${outdir}/${s}.${al}.${ft}.diamond.log.txt
					else
						echo -e "\t::: Searching for valid protein matches step appears to have been performed!"
					fi
				else
					echo -e "\tWarning! Not found ${outdir}/${s}.${al}.${ft}.fastq"
				fi
			done
			
		fi
		
	done
	
	eval "${scriptsdir}/count_Reads.pl -i ${outdir} -s ${s} -x *id.*.bed -o ${outdir}/${s}.gene.counts.txt -l INFO"
	eval "${scriptsdir}/count_Reads.pl -p -i ${outdir} -s ${s} -x *id.*.bed -o ${outdir}/${s}.protein.counts.txt -l INFO"
	ln -f -s $(readlink -f ${outdir}/${s}.protein.counts.txt) ${outdir}/${s}.counts.txt
	
	counts_file_param=(${counts_file_param[@]} "-d ${s}=\"${outdir}/${s}.counts.txt\"")
	biogroup=$(grep "${s}" ${bcfile} | cut -f 3)
	
	if [ ! -v TOEXCLUDE[${biogroup}] ]; then
		GROUP[${biogroup}]="${GROUP[${biogroup}]} ${s}"
	else
		echo "Exclude ${biogroup} = ${s} from the following analysis ..."
	fi

done

groups_param=()
for K in "${!GROUP[@]}"; do
	groups_param=(${groups_param[@]} "-g ${K}=\"$(echo ${GROUP[${K}]} | sed 's/ /,/g')\"")
done

echo "Joining Results [ReadCountsMatrix.txt] ..."

eval "./join_Results.pl ${counts_file_param[*]} ${groups_param[*]} > ${outdir}/ReadCountsMatrix.txt"

filtered_file="${outdir}/FilteredReadCountsMatrix_${filter_mingroups}_${filter_minsamples}_${filter_min}.txt"
annotated_filtered_file="${outdir}/AnnotFilteredReadCountsMatrix_${filter_mingroups}_${filter_minsamples}_${filter_min}.txt"
annotated_filtered_peaks_file="${outdir}/AnnotFilteredReadCountsMatrix_${filter_mingroups}_${filter_minsamples}_${filter_min}_peaks.txt"

echo "Filtering Results [FilteredReadCountsMatrix_${filter_mingroups}_${filter_minsamples}_${filter_min}.txt] ..."

eval "./filter_Results.pl ${groups_param[*]} -i ${outdir}/ReadCountsMatrix.txt -o ${outgroup} -u ${filter_mingroups} -s ${filter_minsamples} -m ${filter_min} > ${filtered_file}"

echo "Annotating Results [AnnotFilteredReadCountsMatrix_${filter_mingroups}_${filter_minsamples}_${filter_min}.txt] ..."

eval "./annot_Results.pl -i ${filtered_file} -g ${genenamefile} -d ${deflinefile} > ${annotated_filtered_file}"

if [ ${fastafile} ]; then

	samps=`echo ${!GROUP[@]} | sed "s/${outgroup}//" | sed 's/^ \+//' | sed 's/ \+$//' | sed 's/ \+/,/g'`
	
	
	echo "Creating coverage plots ..."
	
	mkdir -p ${outdir}/covplot
	
	echo "cut -f 1 ${annotated_filtered_file} | grep -v '^ID' | ./plotCoverage.pl -i ${outdir}/ -g ${samps} -x *id.*.bed -p ${fastafile} -o ${outdir}/covplot"

	echo "Finding peaks [AnnotFilteredReadCountsMatrix_${filter_mingroups}_${filter_minsamples}_${filter_min}_peaks.txt] ..."
	
	eval "cut -f 1 ${annotated_filtered_file} | grep -v '^ID' | ./findPeaks.pl -i ${outdir}/ -g ${samps} -x *id.*.bed -p ${fastafile} -o ${annotated_filtered_peaks_file} -l INFO"

	echo "Plot Hydrophobic/Hydrophilic indexes ..."

	mkdir -p ${outdir}/hydroplot/
	
	eval "./plotHydro.pl -p ${fastafile} -i ${annotated_filtered_peaks_file} -o ${outdir}/hydroplot/"
fi


./mkPEPEstats.pl -i ${outdir} 	-f ${annotated_filtered_peaks_file} \
				-o ${outdir}/STATISTICS.txt \
				-a ${anninfofile} \
				-d ${deflinefile} \
				-b ${bcfile}

