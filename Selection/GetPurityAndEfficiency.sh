#!/bin/bash
# gets the selection efficiency and purity out of the weighted raw numbers. TODO: merge with Sum
# Written by O. Bondu (October 2010)

CMSSW_release="CMSSW_3_9_7_v2"

syntax="${0} {selection version}"
if [[ -z ${1} ]]
then
	echo "ERROR: syntax is ${syntax}"
	exit 2
fi

SelectionVersion=${1}
SELECTEDDIR=/sps/cms/obondu/${CMSSW_release}/src/Zmumugamma/Selection/Selected/${SelectionVersion}
echo -e "*** PurityEfficiency_${SelectionVersion}.txt ***"
rm  PurityEfficiency_${SelectionVersion}.txt

IntegratedLuminosity="36.145000992"

G_Pt_0to15_TuneZ2_7TeV_pythia6="84200000.0"
G_Pt_15to30_TuneZ2_7TeV_pythia6="171700.0"
TTJets_TuneZ2_7TeV_madgraph_tauola="121.0"
DYToMuMu_M_20_CT10_TuneZ2_7TeV_powheg_pythia="1614.0"
WJetsToLNu_TuneZ2_7TeV_madgraph_tauola="24640.0"
G4Jets_Pt_60to120_TuneZ2_7TeV_alpgen="32.08"
G4Jets_Pt_300to5000_TuneZ2_7TeV_alpgen="0.25"
G4Jets_Pt_240to300_TuneZ2_7TeV_alpgen="0.42"
G4Jets_Pt_20to60_TuneZ2_7TeV_alpgen="148.30"
G4Jets_Pt_180to240_TuneZ2_7TeV_alpgen="1.40"
G4Jets_Pt_120to180_TuneZ2_7TeV_alpgen="5.73"
G3Jets_Pt_60to120_TuneZ2_7TeV_alpgen="118.70"
G3Jets_Pt_300to5000_TuneZ2_7TeV_alpgen="0.47"
G3Jets_Pt_240to300_TuneZ2_7TeV_alpgen="0.76"
G3Jets_Pt_20to60_TuneZ2_7TeV_alpgen="794.70"
G3Jets_Pt_180to240_TuneZ2_7TeV_alpgen="3.09"
G3Jets_Pt_120to180_TuneZ2_7TeV_alpgen="15.27"
G2Jets_Pt_60to120_TuneZ2_7TeV_alpgen="425.80"
G2Jets_Pt_300to5000_TuneZ2_7TeV_alpgen="0.72"
G2Jets_Pt_240to300_TuneZ2_7TeV_alpgen="1.45"
G2Jets_Pt_20to60_TuneZ2_7TeV_alpgen="4080.00"
G2Jets_Pt_180to240_TuneZ2_7TeV_alpgen="5.94"
G2Jets_Pt_120to180_TuneZ2_7TeV_alpgen="35.56"
QCD_Pt_20_MuEnrichedPt_15_TuneZ2_7TeV_pythia6="84679.3"

InitialNumberG_Pt_0to15_TuneZ2_7TeV_pythia6="1057100"
InitialNumberG_Pt_15to30_TuneZ2_7TeV_pythia6="1025840"
InitialNumberTTJets_TuneZ2_7TeV_madgraph_tauola="1164732"
InitialNumberDYToMuMu_M_20_CT10_TuneZ2_7TeV_powheg_pythia="1998931"
InitialNumberWJetsToLNu_TuneZ2_7TeV_madgraph_tauola="15123740"
InitialNumberG4Jets_Pt_60to120_TuneZ2_7TeV_alpgen="333214"
InitialNumberG4Jets_Pt_300to5000_TuneZ2_7TeV_alpgen="336836"
InitialNumberG4Jets_Pt_240to300_TuneZ2_7TeV_alpgen="333241"
InitialNumberG4Jets_Pt_20to60_TuneZ2_7TeV_alpgen="335546"
InitialNumberG4Jets_Pt_180to240_TuneZ2_7TeV_alpgen="331020"
InitialNumberG4Jets_Pt_120to180_TuneZ2_7TeV_alpgen="328253"
InitialNumberG3Jets_Pt_60to120_TuneZ2_7TeV_alpgen="333681"
InitialNumberG3Jets_Pt_300to5000_TuneZ2_7TeV_alpgen="311653"
InitialNumberG3Jets_Pt_240to300_TuneZ2_7TeV_alpgen="339273"
InitialNumberG3Jets_Pt_20to60_TuneZ2_7TeV_alpgen="322557"
InitialNumberG3Jets_Pt_180to240_TuneZ2_7TeV_alpgen="324607"
InitialNumberG3Jets_Pt_120to180_TuneZ2_7TeV_alpgen="331691"
InitialNumberG2Jets_Pt_60to120_TuneZ2_7TeV_alpgen="338909"
InitialNumberG2Jets_Pt_300to5000_TuneZ2_7TeV_alpgen="333630"
InitialNumberG2Jets_Pt_240to300_TuneZ2_7TeV_alpgen="331339"
InitialNumberG2Jets_Pt_20to60_TuneZ2_7TeV_alpgen="1782042"
InitialNumberG2Jets_Pt_180to240_TuneZ2_7TeV_alpgen="329962"
InitialNumberG2Jets_Pt_120to180_TuneZ2_7TeV_alpgen="330992"
InitialNumberQCD_Pt_20_MuEnrichedPt_15_TuneZ2_7TeV_pythia6="28979866"

echo -e "*** Raw Yields and errors\n">> PurityEfficiency_${SelectionVersion}.txt
# ************
# First: if there is exactly zero RAW event, overestimate the RAW number of event to 1, to get an overestimate of the purity
# ************

cols=`head -n 4 RawSelectionCutsNumbers_${SelectionVersion}.txt | tail -n 1 | wc -w`
for sample in `tac RawSelectionCutsNumbers_${SelectionVersion}.txt | awk '{print $1}' | head -n -3 | grep -v Run2010 | grep -v G[1-4]Jet`
do
  line=`grep -w "${sample}" RawSelectionCutsNumbers_${SelectionVersion}.txt`
  newline=""
  newerrorline=""
  cuterror=""
  cut=""
  allcuts=""
  allcuterrors=""
	iszero="0"
#  for currentColumn in `seq 2 $cols`
	for currentColumn in "${cols}"
  do
    if [[ ${sample} = "G_Pt_0to15_TuneZ2_7TeV_pythia6" ]]
    then
			currentNumber=`echo ${line} | awk '{print $'"${currentColumn}"'}'`
			if [[ "${currentNumber}" == "0" ]]
			then
				currentNumber="1"
				iszero="1"
			fi
      cut=`echo ${line} | awk '{printf "%4.10f", ('"${G_Pt_0to15_TuneZ2_7TeV_pythia6}"' / '"${InitialNumberG_Pt_0to15_TuneZ2_7TeV_pythia6}"' * '"${currentNumber}"' * '"${IntegratedLuminosity}"')}'`
      cuterror=`echo ${line} | awk '{printf "%4.10f", ('"${G_Pt_0to15_TuneZ2_7TeV_pythia6}"' / '"${InitialNumberG_Pt_0to15_TuneZ2_7TeV_pythia6}"' * sqrt( '"${currentNumber}"' ) * '"${IntegratedLuminosity}"')}'`
    elif [[ ${sample} = "G_Pt_15to30_TuneZ2_7TeV_pythia6" ]]
    then
			currentNumber=`echo ${line} | awk '{print $'"${currentColumn}"'}'`
			if [[ "${currentNumber}" == "0" ]]
			then
				currentNumber="1"
				iszero="1"
			fi
      cut=`echo ${line} | awk '{printf "%4.10f", ('"${G_Pt_15to30_TuneZ2_7TeV_pythia6}"' / '"${InitialNumberG_Pt_15to30_TuneZ2_7TeV_pythia6}"' * '"${currentNumber}"' * '"${IntegratedLuminosity}"')}'`
      cuterror=`echo ${line} | awk '{printf "%4.10f", ('"${G_Pt_15to30_TuneZ2_7TeV_pythia6}"' / '"${InitialNumberG_Pt_15to30_TuneZ2_7TeV_pythia6}"' * sqrt( '"${currentNumber}"' ) * '"${IntegratedLuminosity}"')}'`
    elif [[ ${sample} = "TTJets_TuneZ2_7TeV-madgraph-tauola" ]]
    then
			currentNumber=`echo ${line} | awk '{print $'"${currentColumn}"'}'`
			if [[ "${currentNumber}" == "0" ]]
			then
				currentNumber="1"
				iszero="1"
			fi
      cut=`echo ${line} | awk '{printf "%4.10f", ('"${TTJets_TuneZ2_7TeV_madgraph_tauola}"' / '"${InitialNumberTTJets_TuneZ2_7TeV_madgraph_tauola}"' * '"${currentNumber}"' * '"${IntegratedLuminosity}"')}'`
      cuterror=`echo ${line} | awk '{printf "%4.10f", ('"${TTJets_TuneZ2_7TeV_madgraph_tauola}"' / '"${InitialNumberTTJets_TuneZ2_7TeV_madgraph_tauola}"' * sqrt( '"${currentNumber}"' ) * '"${IntegratedLuminosity}"')}'`
    elif [[ ${sample} = "FSR_DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia" ]]
    then
			currentNumber=`echo ${line} | awk '{print $'"${currentColumn}"'}'`
			if [[ "${currentNumber}" == "0" ]]
			then
				currentNumber="1"
				iszero="1"
			fi
      cut=`echo ${line} | awk '{printf "%4.10f", ('"${DYToMuMu_M_20_CT10_TuneZ2_7TeV_powheg_pythia}"' / '"${InitialNumberDYToMuMu_M_20_CT10_TuneZ2_7TeV_powheg_pythia}"' * '"${currentNumber}"' * '"${IntegratedLuminosity}"')}'`
      cuterror=`echo ${line} | awk '{printf "%4.10f", ('"${DYToMuMu_M_20_CT10_TuneZ2_7TeV_powheg_pythia}"' / '"${InitialNumberDYToMuMu_M_20_CT10_TuneZ2_7TeV_powheg_pythia}"' * sqrt( '"${currentNumber}"' ) * '"${IntegratedLuminosity}"')}'`
    elif [[ ${sample} = "nonFSR_DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia" ]]
    then
			currentNumber=`echo ${line} | awk '{print $'"${currentColumn}"'}'`
			if [[ "${currentNumber}" == "0" ]]
			then
				currentNumber="1"
				iszero="1"
			fi
      cut=`echo ${line} | awk '{printf "%4.10f", ('"${DYToMuMu_M_20_CT10_TuneZ2_7TeV_powheg_pythia}"' / '"${InitialNumberDYToMuMu_M_20_CT10_TuneZ2_7TeV_powheg_pythia}"' * '"${currentNumber}"' * '"${IntegratedLuminosity}"')}'`
      cuterror=`echo ${line} | awk '{printf "%4.10f", ('"${DYToMuMu_M_20_CT10_TuneZ2_7TeV_powheg_pythia}"' / '"${InitialNumberDYToMuMu_M_20_CT10_TuneZ2_7TeV_powheg_pythia}"' * sqrt( '"${currentNumber}"' ) * '"${IntegratedLuminosity}"')}'`
    elif [[ ${sample} = "WJetsToLNu_TuneZ2_7TeV-madgraph-tauola" ]]
    then
			currentNumber=`echo ${line} | awk '{print $'"${currentColumn}"'}'`
			if [[ "${currentNumber}" == "0" ]]
			then
				currentNumber="1"
				iszero="1"
			fi
      cut=`echo ${line} | awk '{printf "%4.10f", ('"${WJetsToLNu_TuneZ2_7TeV_madgraph_tauola}"' / '"${InitialNumberWJetsToLNu_TuneZ2_7TeV_madgraph_tauola}"' * '"${currentNumber}"' * '"${IntegratedLuminosity}"')}'`
      cuterror=`echo ${line} | awk '{printf "%4.10f", ('"${WJetsToLNu_TuneZ2_7TeV_madgraph_tauola}"' / '"${InitialNumberWJetsToLNu_TuneZ2_7TeV_madgraph_tauola}"' * sqrt( '"${currentNumber}"' ) * '"${IntegratedLuminosity}"')}'`
    elif [[ ${sample} = "QCD_Pt-20_MuEnrichedPt-15_TuneZ2_7TeV-pythia6" ]]
    then
			currentNumber=`echo ${line} | awk '{print $'"${currentColumn}"'}'`
			if [[ "${currentNumber}" == "0" ]]
			then
				currentNumber="1"
				iszero="1"
			fi
      cut=`echo ${line} | awk '{printf "%4.10f", ('"${QCD_Pt_20_MuEnrichedPt_15_TuneZ2_7TeV_pythia6}"' / '"${InitialNumberQCD_Pt_20_MuEnrichedPt_15_TuneZ2_7TeV_pythia6}"' * '"${currentNumber}"' * '"${IntegratedLuminosity}"')}'`
      cuterror=`echo ${line} | awk '{printf "%4.10f", ('"${QCD_Pt_20_MuEnrichedPt_15_TuneZ2_7TeV_pythia6}"' / '"${InitialNumberQCD_Pt_20_MuEnrichedPt_15_TuneZ2_7TeV_pythia6}"' * sqrt( '"${currentNumber}"' ) * '"${IntegratedLuminosity}"')}'`
#    elif [[ ${sample} = "SAMPLE" ]]
#    then
#			currentNumber=`echo ${line} | awk '{print $'"${currentColumn}"'}'`
#			if [[ "${currentNumber}" == "0" ]]
#			then
#				currentNumber="1"
#				iszero="1"
#			fi
#      cut=`echo ${line} | awk '{printf "%4.10f", ('"${SAMPLE}"' / '"${InitialNumberSAMPLE}"' * '"${currentNumber}"' * '"${IntegratedLuminosity}"')}'`
#      cuterror=`echo ${line} | awk '{printf "%4.10f", ('"${SAMPLE}"' / '"${InitialNumberSAMPLE}"' * sqrt( '"${currentNumber}"' ) * '"${IntegratedLuminosity}"')}'`
    else
			currentNumber=`echo ${line} | awk '{print $'"${currentColumn}"'}'`
			if [[ "${currentNumber}" == "0" ]]
			then
				currentNumber="1"
				iszero="1"
			fi
      cut=`echo ${line} | awk '{print '"${currentNumber}"'}'`
      cuterror=`echo ${line} | awk '{print sqrt( '"${currentNumber}"' )}'`
    fi
    allcuts="${allcuts}\t${cut}"
    allcuterrors="${allcuterrors}\t${cuterror}"
  done
	if [[ "${iszero}" == "1" ]]
	then
		newline="FLAGGED_${sample}\t${allcuts}"
		newerrorline="ERROR_FLAGGED_${sample}\t${allcuterrors}"
	else
	  newline="${sample}\t${allcuts}"
	  newerrorline="ERROR_${sample}\t${allcuterrors}"
	fi
echo -e ${newline} >> PurityEfficiency_${SelectionVersion}_TEMP.txt
echo -e ${newerrorline} >> PurityEfficiency_${SelectionVersion}_TEMP.txt
echo -e ${newline} >> PurityEfficiency_${SelectionVersion}.txt
echo -e ${newerrorline} >> PurityEfficiency_${SelectionVersion}.txt
done
echo -e "\n\n" >> PurityEfficiency_${SelectionVersion}.txt


echo -e "*** Summed yields\n" >> PurityEfficiency_${SelectionVersion}.txt
# ************
# Sum samples
# ************
for sample in `echo "G_Pt_ TTJets FSR_DYToMuMu nonFSR_DYToMuMu WJets QCD"`
do
	allcuts=""
	if [[ "${sample}" = "FSR_DYToMuMu" ]]
	then
		cols=`cat PurityEfficiency_${SelectionVersion}_TEMP.txt | grep ${sample} | grep -v nonFSR | grep -v ERROR | head -n +1 | wc -w`
	else
		cols=`cat PurityEfficiency_${SelectionVersion}_TEMP.txt | grep ${sample} | grep -v ERROR | head -n +1 | wc -w`
	fi
	for column in ${cols}
	do
		if [[ "${sample}" = "FSR_DYToMuMu" ]]
		then
			currentSum=`cat PurityEfficiency_${SelectionVersion}_TEMP.txt | grep ${sample} | grep -v nonFSR | grep -v ERROR | awk 'BEGIN {SUM=0} {SUM+=$'"${column}"'} END {printf "%4.10f", SUM}'`
			cat PurityEfficiency_${SelectionVersion}_TEMP.txt | grep ${sample} | grep -v nonFSR | grep -v ERROR | grep "FLAGGED" > /dev/null
			if [[ "$?" == "0" ]]
			then
				currentMax=`cat PurityEfficiency_${SelectionVersion}_TEMP.txt | grep ${sample} | grep -v nonFSR | grep -v ERROR | awk 'BEGIN {MAX=0} {MAX=(MAX > $2 ? MAX : $2 ) } END {print MAX}'`
				allcuts=`echo "${allcuts}\t MORE_THAN ${currentMax}"`
			else
				allcuts=`echo "${allcuts}\t${currentSum}"`
			fi
		else
			currentSum=`cat PurityEfficiency_${SelectionVersion}_TEMP.txt | grep ${sample} | grep -v ERROR | awk 'BEGIN {SUM=0} {SUM+=$'"${column}"'} END {printf "%4.10f", SUM}'`
			cat PurityEfficiency_${SelectionVersion}_TEMP.txt | grep ${sample} | grep -v ERROR | grep "FLAGGED" > /dev/null
			if [[ "$?" == "0" ]]
			then
				currentMax=`cat PurityEfficiency_${SelectionVersion}_TEMP.txt | grep ${sample} | grep -v ERROR | awk 'BEGIN {MAX=0} {MAX=(MAX > $2 ? MAX : $2 ) } END {print MAX}'`
				allcuts=`echo "${allcuts}\t MORE_THAN ${currentMax}"`
			else
				allcuts=`echo "${allcuts}\t${currentSum}"`
			fi
		fi
#		cat PurityEfficiency_${SelectionVersion}_TEMP.txt | grep ${sample} | grep -v ERROR | grep "FLAGGED" > /dev/null
	done
echo -e "${sample}\t${allcuts}" >> PurityEfficiency_${SelectionVersion}_TEMP2.txt
echo -e "${sample}\t${allcuts}" >> PurityEfficiency_${SelectionVersion}.txt
done
echo -e "\n\n" >> PurityEfficiency_${SelectionVersion}.txt


echo -e "*** Root mean square of errors\n" >> PurityEfficiency_${SelectionVersion}.txt
# ************
# Sum errors (root mean square)
# ************
for sample in `echo "G_Pt_ TTJets FSR_DYToMuMu nonFSR_DYToMuMu WJets QCD"`
do
	allcuterrors=""
	if [[ "${sample}" = "FSR_DYToMuMu" ]]
	then
		cols=`cat PurityEfficiency_${SelectionVersion}_TEMP.txt | grep ${sample} | grep -v nonFSR | grep ERROR | head -n +1 | wc -w`
	else
		cols=`cat PurityEfficiency_${SelectionVersion}_TEMP.txt | grep ${sample} | grep ERROR | head -n +1 | wc -w`
	fi
	for column in ${cols}
	do
		if [[ "${sample}" = "FSR_DYToMuMu" ]]
		then
			currentSum=`cat PurityEfficiency_${SelectionVersion}_TEMP.txt | grep ${sample} | grep -v nonFSR | grep ERROR | awk 'BEGIN {SUM=0} {SUM += $'"${column}"' * $'"${column}"' } END {printf "%4.10f", sqrt( SUM )}'`
	#		allcuterrors=`echo "${allcuterrors}\t${currentSum}"`
			cat PurityEfficiency_${SelectionVersion}_TEMP.txt | grep ${sample} | grep -v nonFSR | grep ERROR | grep "FLAGGED" > /dev/null
			if [[ "$?" == "0" ]]
			then
				currentMax=`cat PurityEfficiency_${SelectionVersion}_TEMP.txt | grep ${sample} | grep -v nonFSR | grep ERROR | awk 'BEGIN {MAX=0} {MAX=(MAX > $2 ? MAX : $2 ) } END {print MAX}'`
				allcuterrors=`echo "${allcuterrors}\t${currentMax}"`
			else
				allcuterrors=`echo "${allcuterrors}\t${currentSum}"`
			fi
		else
			currentSum=`cat PurityEfficiency_${SelectionVersion}_TEMP.txt | grep ${sample} | grep ERROR | awk 'BEGIN {SUM=0} {SUM += $'"${column}"' * $'"${column}"' } END {printf "%4.10f", sqrt( SUM )}'`
	#		allcuterrors=`echo "${allcuterrors}\t${currentSum}"`
			cat PurityEfficiency_${SelectionVersion}_TEMP.txt | grep ${sample} | grep ERROR | grep "FLAGGED" > /dev/null
			if [[ "$?" == "0" ]]
			then
				currentMax=`cat PurityEfficiency_${SelectionVersion}_TEMP.txt | grep ${sample} | grep ERROR | awk 'BEGIN {MAX=0} {MAX=(MAX > $2 ? MAX : $2 ) } END {print MAX}'`
				allcuterrors=`echo "${allcuterrors}\t${currentMax}"`
			else
				allcuterrors=`echo "${allcuterrors}\t${currentSum}"`
			fi
		fi
	done
echo -e "ERROR_${sample}\t${allcuterrors}" >> PurityEfficiency_${SelectionVersion}_TEMP3.txt
echo -e "ERROR_${sample}\t${allcuterrors}" >> PurityEfficiency_${SelectionVersion}.txt
done
echo -e "\n\n" >> PurityEfficiency_${SelectionVersion}.txt








# ************
# COMPUTE PURITY
# ************
signal="FSR_DYToMuMu"
BGList="G_Pt_ TTJets nonFSR_DYToMuMu WJets QCD"
isExact="1"
allBGs="0.0"
allBGerrors="0.0"
for background in `echo "${BGList}"`
do
	if [[ ${isExact} == "1" ]]
	then
		cat PurityEfficiency_${SelectionVersion}_TEMP2.txt | grep "MORE_THAN" | grep ${background} > /dev/null
		if [[ $? == "0" ]]
		then
			isExact="0"
		fi
	fi
	cols=`cat PurityEfficiency_${SelectionVersion}_TEMP2.txt | grep ${background} | wc -w`
	errorcols=`cat PurityEfficiency_${SelectionVersion}_TEMP3.txt | grep ERROR | grep ${background} | wc -w`
	currentBG=`cat PurityEfficiency_${SelectionVersion}_TEMP2.txt | grep ${background} | awk '{print $'"${cols}"'}'`
	currentBGerror=`cat PurityEfficiency_${SelectionVersion}_TEMP3.txt | grep ${background} | awk '{print $'"${errorcols}"'}'`
	allBGs=`echo "${allBGs} + ${currentBG}"`
	allBGerrors=`echo "${allBGerrors} + ${currentBGerror} * ${currentBGerror}"`
done
echo "SIGNAL = ${signal}" >> PurityEfficiency_${SelectionVersion}.txt
cols=`cat PurityEfficiency_${SelectionVersion}_TEMP2.txt | grep ${signal} | grep -v nonFSR | wc -w`
errorcols=`cat PurityEfficiency_${SelectionVersion}_TEMP3.txt | grep ERROR | grep ${signal} | grep -v nonFSR | wc -w`
signalYield=`cat PurityEfficiency_${SelectionVersion}_TEMP2.txt | grep ${signal} | grep -v nonFSR | awk '{print $'"${cols}"'}'`
signalYieldError=`cat PurityEfficiency_${SelectionVersion}_TEMP3.txt | grep ERROR | grep ${signal} | grep -v nonFSR | awk '{print $'"${cols}"'}'`
echo "signal yield = ${signalYield} +- ${signalYieldError}" >> PurityEfficiency_${SelectionVersion}.txt
BGYield=`awk 'BEGIN {printf "%4.10f", ( '"${allBGs}"' )}'`
BGYieldError=`awk 'BEGIN {printf "%4.10f", ( sqrt( '"${allBGerrors}"' ) )}'`
echo "BG =  ${BGList}" >> PurityEfficiency_${SelectionVersion}.txt
echo "BG yield : ${BGYield} +- ${BGYieldError}" >> PurityEfficiency_${SelectionVersion}.txt
TotalYield=`awk 'BEGIN {printf "%4.10f", ( '"${BGYield}"' +  '"${signalYield}"' )}'`
TotalYieldError=`awk 'BEGIN {printf "%4.10f", ( sqrt( '"${BGYieldError}"' * '"${BGYieldError}"' +  '"${signalYieldError}"' * '"${signalYieldError}"' ) )}'`
echo "Total Yield = ${TotalYield} +- ${TotalYieldError}"  >> PurityEfficiency_${SelectionVersion}.txt
PURITY=`awk 'BEGIN { printf "%4.10f", ('"${signalYield}"') / ('"${signalYield}"' + '"${allBGs}"')}'`
PURITYError=`awk 'BEGIN { printf "%4.10f", ('"${PURITY}"' * ( '"${signalYieldError}"' / '"${signalYield}"' + '"${TotalYieldError}"' / '"${TotalYield}"' ) ) }'`
if [[ "${isExact}" == "0" ]]
then
	echo "PURITY < ${PURITY} +- ${PURITYError}" >> PurityEfficiency_${SelectionVersion}.txt
else
	echo "PURITY = ${PURITY} +- ${PURITYError}" >> PurityEfficiency_${SelectionVersion}.txt
fi



echo "" >> PurityEfficiency_${SelectionVersion}.txt

# ************
# COMPUTE PURITY
# ************
signal="FSR_DYToMuMu"
BGList="TTJets nonFSR_DYToMuMu WJets QCD"
isExact="1"
allBGs="0.0"
allBGerrors="0.0"
for background in `echo "${BGList}"`
do
	if [[ ${isExact} == "1" ]]
	then
		cat PurityEfficiency_${SelectionVersion}_TEMP2.txt | grep "MORE_THAN" | grep ${background} > /dev/null
		if [[ $? == "0" ]]
		then
			isExact="0"
		fi
	fi
	cols=`cat PurityEfficiency_${SelectionVersion}_TEMP2.txt | grep ${background} | wc -w`
	errorcols=`cat PurityEfficiency_${SelectionVersion}_TEMP3.txt | grep ERROR | grep ${background} | wc -w`
	currentBG=`cat PurityEfficiency_${SelectionVersion}_TEMP2.txt | grep ${background} | awk '{print $'"${cols}"'}'`
	currentBGerror=`cat PurityEfficiency_${SelectionVersion}_TEMP3.txt | grep ${background} | awk '{print $'"${errorcols}"'}'`
	allBGs=`echo "${allBGs} + ${currentBG}"`
	allBGerrors=`echo "${allBGerrors} + ${currentBGerror} * ${currentBGerror}"`
done
echo "SIGNAL = ${signal}" >> PurityEfficiency_${SelectionVersion}.txt
cols=`cat PurityEfficiency_${SelectionVersion}_TEMP2.txt | grep ${signal} | grep -v nonFSR | wc -w`
errorcols=`cat PurityEfficiency_${SelectionVersion}_TEMP3.txt | grep ERROR | grep ${signal} | grep -v nonFSR | wc -w`
signalYield=`cat PurityEfficiency_${SelectionVersion}_TEMP2.txt | grep ${signal} | grep -v nonFSR | awk '{print $'"${cols}"'}'`
signalYieldError=`cat PurityEfficiency_${SelectionVersion}_TEMP3.txt | grep ERROR | grep ${signal} | grep -v nonFSR | awk '{print $'"${cols}"'}'`
echo "signal yield = ${signalYield} +- ${signalYieldError}" >> PurityEfficiency_${SelectionVersion}.txt
BGYield=`awk 'BEGIN {printf "%4.10f", ( '"${allBGs}"' )}'`
BGYieldError=`awk 'BEGIN {printf "%4.10f", ( sqrt( '"${allBGerrors}"' ) )}'`
echo "BG =  ${BGList}" >> PurityEfficiency_${SelectionVersion}.txt
echo "BG yield : ${BGYield} +- ${BGYieldError}" >> PurityEfficiency_${SelectionVersion}.txt
TotalYield=`awk 'BEGIN {printf "%4.10f", ( '"${BGYield}"' +  '"${signalYield}"' )}'`
TotalYieldError=`awk 'BEGIN {printf "%4.10f", ( sqrt( '"${BGYieldError}"' * '"${BGYieldError}"' +  '"${signalYieldError}"' * '"${signalYieldError}"' ) )}'`
echo "Total Yield = ${TotalYield} +- ${TotalYieldError}"  >> PurityEfficiency_${SelectionVersion}.txt
PURITY=`awk 'BEGIN { printf "%4.10f", ('"${signalYield}"') / ('"${signalYield}"' + '"${allBGs}"')}'`
PURITYError=`awk 'BEGIN { printf "%4.10f", ('"${PURITY}"' * ( '"${signalYieldError}"' / '"${signalYield}"' + '"${TotalYieldError}"' / '"${TotalYield}"' ) ) }'`
if [[ "${isExact}" == "0" ]]
then
	echo "PURITY < ${PURITY} +- ${PURITYError}" >> PurityEfficiency_${SelectionVersion}.txt
else
	echo "PURITY = ${PURITY} +- ${PURITYError}" >> PurityEfficiency_${SelectionVersion}.txt
fi



echo "" >> PurityEfficiency_${SelectionVersion}.txt

# ************
# COMPUTE SIGNAL EFFICIENCY
# ************
line=`cat SelectionCutsNumbersSummed_${SelectionVersion}.txt | grep ${signal} | grep -v nonFSR `
cols=`echo ${line} | wc -w`
EFFICIENCY=`echo ${line} | awk '{printf "%4.10f", ( $'"${cols}"' / $5 )}'`
echo "EFFICIENCY = ${EFFICIENCY}" >> PurityEfficiency_${SelectionVersion}.txt
#echo "EFFICIENCY = ${EFFICIENCY}"
echo "" >> PurityEfficiency_${SelectionVersion}.txt
echo "" >> PurityEfficiency_${SelectionVersion}.txt

rm PurityEfficiency_${SelectionVersion}_TEMP.txt
rm PurityEfficiency_${SelectionVersion}_TEMP2.txt
rm PurityEfficiency_${SelectionVersion}_TEMP3.txt
exit 0
