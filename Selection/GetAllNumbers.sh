#!/bin/bash
# Small script to run all GetRawSelectionCutsNumbers.sh in one shot
# Written by Olivier Bondu (March 2010)

#for selection in `echo "hadEt noMuIso-hadEt sumPt noMuIso-sumPt"`
#for selection in `echo "hadEt-lowDeltaRmin hadEt-noDeltaRmin-relaxedMuEta hadEt-noDeltaRmin-relaxedpT hadEt-noDeltaRmin-tightedPtMu"`
#for selection in `echo "hadEt-lowDeltaRmin"`
#for selection in `echo "hadEt-noDeltaRmin-relaxedMuEta"`
for selection in `echo "hadEt-noDeltaRmin-relaxedMuEta hadEt-noDeltaRmin-relaxedMuEta-relaxedMMGv2"`
#for selection in `echo "hadEt-noDeltaRmin-relaxedMuEta-relaxedMMG hadEt-noDeltaRmin-relaxedMuEta-relaxedMMGv2"`
#for selection in `echo "hadEt-noDeltaRmin-relaxedpT"`
#for selection in `echo "hadEt-noDeltaRmin-tightedPtMu"`
#for selection in `echo "hadEt-noDeltaRmin-relaxedpT-looseWindow"`
#for selection in `echo "hadEt-noDeltaRmin-singleTightedPtMu-0.950EScale"`
do
	./GetRawSelectionCutsNumbers.sh ${selection}
	./ConvertRawNumbersIntoNumbers.sh ${selection}
	./SumNumbersBySample.sh ${selection}
	./GetSelectionNumbers.sh ${selection}
#	echo "HaddAllRootFiles.sh ${selection}"
	eval `./HaddAllRootFiles.sh ${selection}`
	./GetPurityAndEfficiency.sh ${selection}
	cd LaTeXTables/
	eval `pdflatex SelectionCutsNumbersEfficiencies_${selection}.tex`
	cd ..
	echo -e "\n\n"
done

exit 0

