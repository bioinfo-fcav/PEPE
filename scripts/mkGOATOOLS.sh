#!/bin/bash

annttype=$1

if [ ! ${annttype} ]; then
	echo "Missing annotation type: phytozome or agrigo"
	exit
else
	if [ ! ${annttype} == "phytozome" ] && [ ! ${annttype} == "agrigo" ]; then
		echo "Please choose annotation type: phytozome or agrigo"
		exit
	fi
fi

anntfile=$2

resfile=$3

outdir=$4

if [ ! ${anntfile} ]; then
	echo "Missing <Phytozome genome>.annotation_info.txt file."
	exit
else
	if [ ! -e ${anntfile} ]; then
		echo "Wrong <Phytozome genome>.annotation_info.txt file (${anntfile})"
		exit
	fi
fi

if [ ! ${resfile} ]; then
	echo "Missing PEPE final result file."
	exit
else
	if [ ! -e ${resfile} ]; then
		echo "Wrong PEPE final result file (${resfile})"
		exit
	fi
fi

if [ ! ${outdir} ]; then
	echo "Missing output directory."
	exit
else
	if [ ! -d ${outdir} ]; then
		echo "Wrong output directory (${outdir})"
		exit
	fi
fi

runid=${RANDOM}
runoutdir=${outdir}/pepego.${runid}
mkdir -p ${runoutdir}

echo "Running ID ${runid} ..."

cut -f 1 ${resfile} | grep -v '^ID' | sort -u > ${runoutdir}/study

if [ ${annttype} == "phytozome" ]; then
	cut -f 2,10 ${anntfile} | grep -v '^locusname' | sed 's/,/;/g' | sort -u > ${runoutdir}/association
elif [ ${annttype} == "agrigo" ]; then
	grep '^TAIRgenemodel' ${anntfile} | perl -F"\t" -lane ' INIT {our %data;} my ($g)=$F[1]=~/^([^\.]+)\./; $data{$g}->{$F[2]}=undef; END { foreach my $g (keys %data) { print $g,"\t",join(";", keys %{ $data{$g} }); }  }' > ${runoutdir}/association
else
	echo "Wrong annotation type: ${annttype}"
	exit
fi

cut -f 1 ${runoutdir}/association | sort -u > ${runoutdir}/population


GO_OBO_FILE=go-basic.obo
GO_OBO_DOWNLOAD=http://geneontology.org/ontology/${GO_OBO_FILE}

wget --timestamping ${GO_OBO_DOWNLOAD} 2> /dev/null

find_enrichment.py ${runoutdir}/study ${runoutdir}/population ${runoutdir}/association --outfile=${outdir}/goea_all.xlsx,${outdir}/goea_all.tsv --pvalcalc=fisher_scipy_stats --obo=${GO_OBO_FILE} --pval=1 --alpha=0.05


rm -fr ${runoutdir}
