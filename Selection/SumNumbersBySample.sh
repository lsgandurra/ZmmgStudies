#!/usr/local/bin/bash
# Script to sum Numbers by sample
# Written by O. Bondu (August 2010)

CMSSW_release="CMSSW_3_9_7_v2"

syntax="${0} {SelectionVersion}"

if [[ -z ${1} ]]
then
	echo ${syntax}
	exit 1
fi
SelectionVersion=${1}

echo "*** SelectionCutsNumbersSummed_${SelectionVersion}.txt ***"
echo -e "*** RESULTS ***" > SelectionCutsNumbersSummed_${SelectionVersion}.txt
echo -e "Sample \t0\tPthatFilter\tCSA07ID\tZJETVETO\t1.a\t1.b\t1.c\t1.d\t1.e\t2.a\t2.b\t2.c\t3\t4\t5\t6\t7\t8\t9\t10\tSelected" >> SelectionCutsNumbersSummed_${SelectionVersion}.txt
echo "" >> SelectionCutsNumbersSummed_${SelectionVersion}.txt


for sample in `echo "Run2010 G_Pt_ TTJets_TuneZ2_7TeV-madgraph-tauola FSR_DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia nonFSR_DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia WJetsToLNu_TuneZ2_7TeV-madgraph-tauola G[0-4]Jet"`
do
	allcuts=""
#	cat SelectionCutsNumbers_${SelectionVersion}.txt | grep ${sample}
	cols=`cat SelectionCutsNumbers_${SelectionVersion}.txt | grep ${sample} | head -n +1 | wc -w`
	for column in `seq 2 ${cols}`
	do
		if [[ ${sample} = "part" ]]
		then
			currentSum=`cat SelectionCutsNumbers_${SelectionVersion}.txt | grep ${sample} | awk 'BEGIN {SUM=0} {SUM+=$'"${column}"'} END {printf "%4i", SUM}'`
		else
			currentSum=`cat SelectionCutsNumbers_${SelectionVersion}.txt | grep ${sample} | awk 'BEGIN {SUM=0} {SUM+=$'"${column}"'} END {printf "%4.5f", SUM}'`
		fi
#		echo "${column} ${currentSum}"
		allcuts=`echo "${allcuts}\t${currentSum}"`
	done
	echo -e "${sample}\t${allcuts}" >> SelectionCutsNumbersSummed_${SelectionVersion}.txt
done

exit 0

