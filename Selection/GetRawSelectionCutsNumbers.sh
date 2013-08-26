#!/usr/local/bin/bash
# Script to fetch cut informations
# Created by Olivier Bondu for CMSSW_1_6_12
# Modified for 336 (January 2010)
# Modified for 361p4 (August 2010)

CMSSW_release="CMSSW_3_9_7_v2"

syntax="${0} {SelectionVersion}"

if [[ -z ${1} ]]
then
	echo ${syntax}
	exit 1
fi
SelectionVersion=${1}
SELECTEDDIR=/sps/cms/obondu/${CMSSW_release}/src/Zmumugamma/Selection/Selected/${SelectionVersion}

echo "*** RawSelectionCutsNumbers_${SelectionVersion}.txt ***"

rm RawSelectionCutsNumbers_${SelectionVersion}.txt
echo -e "*** RAW RESULTS ***" > RawSelectionCutsNumbers_${SelectionVersion}.txt
echo -e "Sample \t0\tPthatFilter\tCSA07ID\tZJETVETO\t1.a\t1.b\t1.c\t1.d\t1.e\t2.a\t2.b\t2.c\t3\t4\t5\t6\t7\t8\t9\t10\tSelected" >> RawSelectionCutsNumbers_${SelectionVersion}.txt
echo "" >> RawSelectionCutsNumbers_${SelectionVersion}.txt

for sample in `'ls' -l -r ${SELECTEDDIR} | grep drw | awk '{print $9}' | grep -v OLD`
do
	for file in `'ls' -r ${SELECTEDDIR}/${sample} | grep ${sample}_ | grep _${SelectionVersion}.out`
	do
		SAMPLEandNb=`echo ${file} | cut -d . -f -1 | rev | cut -d _ -f 2- | rev`
	  total=`grep "nBeforeAllCuts=" ${SELECTEDDIR}/${sample}/${file} | awk '{print $2}'`
#		Pthat=`grep "nAfterCutPthatFilter=" ${SELECTEDDIR}/${sample}/${file} | awk '{print $2}'`
#	  csa07id=`grep "nAfterCutCSA07ID=" ${SELECTEDDIR}/${sample}/${file} | awk '{print $2}'`
#	  zjetveto=`grep "nAfterCutZJETVETO=" ${SELECTEDDIR}/${sample}/${file} | awk '{print $2}'`
	  allcuts=""
#	  allcuts="${allcuts}\t${total}\t${Pthat}\t${csa07id}\t${zjetveto}"
	  allcuts="${allcuts}\t${total}"
	  for cutvalue in `grep "nAfterCut" ${SELECTEDDIR}/${sample}/${file} | awk '{print $2}'`
	  do
      allcuts="${allcuts}\t${cutvalue}"
	  done
#	  if [ "${SAMPLEandNb}" != "7TeV_Zmumugamma" ]
#	  then
#			if [ "${SAMPLEandNb}" != "7TeV_Zmumugamma_10GeV" ]
#			then
#			  cutvalue=`echo "scale=1 ; $cutvalue/5.0" |bc -l`
#			  allcuts="${allcuts}\t${cutvalue}"
#			fi
#	  fi
		# echo -e "${SAMPLEandNb}\t( ${NbOfFiles} )${allcuts}"
		echo -e  "${SAMPLEandNb} \t${allcuts}" >> TEMPRESULTS_${SelectionVersion}.txt
		total=""
#		Pthat=""
#		csa07id=""
		allcuts=""
		cutvalue=""																								
		done # end of loop over sample files
done # end of loop over samples


for sample in `'ls' -l -r ${SELECTEDDIR} | grep drw | awk '{print $9}' | grep -v OLD`
do
	cols=`cat TEMPRESULTS_${SelectionVersion}.txt | grep "${sample}_0" | wc -w`
  SubFiles=`cat TEMPRESULTS_${SelectionVersion}.txt | grep "${sample}_" | wc -l`
#FIXME  NbOfFiles=`cat TEMPRESULTS_${SelectionVersion}.txt | grep "${sample}" | awk 'BEGIN {NbOfFiles=0} {NbOfFiles+=$3} END {print NbOfFiles}'`
  allcuts=""
#	echo "$cols $SubFiles $NbOfFiles"
  for currentColumn in `seq 2 $cols`
  do
		currentSum=`cat TEMPRESULTS_${SelectionVersion}.txt | grep "${sample}_" | awk 'BEGIN {SUM=0} {SUM+=$'"${currentColumn}"'} END {print SUM}'`
		allcuts="${allcuts}\t${currentSum}"
  done
#	if [[ `echo ${sample} | grep 7TeV_Zmumugamma` ]]
#	then
#		if [[ `echo ${sample} | grep 7TeV_Zmumugamma_10GeV` ]]
#		then
#			cat TEMPRESULTS_${SelectionVersion}.txt | grep -v "${sample}" > TEMPRESULTSBIS_${SelectionVersion}.txt
#		else
#			cat TEMPRESULTS_${SelectionVersion}.txt | grep -v "${sample}" > TEMPRESULTSBIS_${SelectionVersion}.txt
#			cat TEMPRESULTS_${SelectionVersion}.txt | grep "${sample}_10GeV" >> TEMPRESULTSBIS_${SelectionVersion}.txt
#		fi
#	else
		cat TEMPRESULTS_${SelectionVersion}.txt | grep -v "${sample}_" > TEMPRESULTSBIS_${SelectionVersion}.txt
#	fi
#FIXME  echo -e "${sample}\t(${NbOfFiles})${allcuts}" >> TEMPRESULTSBIS_${SelectionVersion}.txt
  echo -e "${sample}\t${allcuts}" >> TEMPRESULTSBIS_${SelectionVersion}.txt
  mv TEMPRESULTSBIS_${SelectionVersion}.txt TEMPRESULTS_${SelectionVersion}.txt
done # end of loop over samples

cat TEMPRESULTS_${SelectionVersion}.txt | sort >> RawSelectionCutsNumbers_${SelectionVersion}.txt

rm TEMPRESULTS_${SelectionVersion}.txt

exit 0


