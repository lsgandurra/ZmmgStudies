#!/bin/bash
# Small script to run all GetRawSelectionCutsNumbers.sh in one shot
# Written by Olivier Bondu (March 2010)

#for selection in `echo "hadEt noMuIso-hadEt sumPt noMuIso-sumPt"`
for selection in `echo "noMuIso-hadEt-noDeltaRmin"`
do
	./GetRawSelectionCutsNumbers.sh ${selection}
	./ConvertRawNumbersIntoNumbers.sh ${selection}
	./SumNumbersBySample.sh ${selection}
	./GetSelectionNumbers.sh ${selection}
	./HaddAllRootFiles.sh ${selection}
	cd LaTeXTables/
	eval `pdflatex SelectionCutsNumbersEfficiencies_${selection}.tex`
done

exit 0

