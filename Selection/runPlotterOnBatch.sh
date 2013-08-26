#!/bin/bash
# Script to launch plots on batch
# Written by O. Bondu (August 2010)

time=828600
queue=T
memory="4096MB"

#for selection in `echo "hadEt-noDeltaRmin-relaxedMuEta hadEt-noDeltaRmin-relaxedpT hadEt-noDeltaRmin-tightedPtMu hadEt-lowDeltaRmin"`
for selection in `echo "hadEt-noDeltaRmin-relaxedpT-looseWindow"`
do
	if [[ ! -d Plots_${selection} ]]
	then
		mkdir Plots_${selection}
	fi

	cd Plots_${selection}
	
	if [[ ! -d gif ]]
	then
		mkdir gif eps pdf
	fi

# copy files to current directory
	cp ../setTDRStyle.C .
	cp ../rootlogon.C .
	cp ../DrawDataMC.C .
	cp ../DrawDataMC.h .
	cp ../plotDataMC_TDR_miniTree.C .
	cp ../Makefile .
	cp ../batch_template.sh batch_plots_${selection}.sh

# pre-batch file
	chmod u+x batch_plots_${selection}.sh
	sed -i -e "s/TIME/${time}/g" batch_plots_${selection}.sh
	sed -i -e "s/QUEUE/${queue}/g" batch_plots_${selection}.sh
	sed -i -e "s/MEMORY/${memory}/g" batch_plots_${selection}.sh
	sed -i -e "s,EXEDIR,Zmumugamma/Selection/Plots_${selection},g" batch_plots_${selection}.sh

	for cut in `echo "isBeforeAllCuts isAfterCut1a isAfterCut1b isAfterCut1c isAfterCut1d isAfterCut1e isAfterCut2a isAfterCut2b isAfterCut2c isAfterCut3 isAfterCut4 isAfterCut5 isAfterCut6 isAfterCut7 isAfterCut8"`
#	for cut in `echo "isAfterCut1a"`
#for cut in `echo "isBeforeAllCuts"`
	do
	# batch file
		name=`echo "plot_${selection}_${cut}"`
		outlog=`echo "plot_${selection}_${cut}"`
		errlog=`echo "plot_${selection}_${cut}"`
		macro=`echo "plotDataMC_${selection}_${cut}"`

		cp batch_plots_${selection}.sh batch_plots_${selection}_${cut}.sh
		sed -i -e "s/NAME/${name}/g" batch_plots_${selection}_${cut}.sh
		sed -i -e "s/OUTLOG/${outlog}/g" batch_plots_${selection}_${cut}.sh
		sed -i -e "s/ERRLOG/${errlog}/g" batch_plots_${selection}_${cut}.sh
		sed -i -e "s/MACRO/${macro}/g" batch_plots_${selection}_${cut}.sh

	# makefile
		cp Makefile Makefile_${selection}_${cut}
		sed -i -e "s/plotDataMC_TEST/${macro}/g" Makefile_${selection}_${cut}
		sed -i -e "s/plotDataMC_TDR_miniTree/plotDataMC_${selection}_${cut}/g" Makefile_${selection}_${cut}

	# macro
		cp plotDataMC_TDR_miniTree.C plotDataMC_${selection}_${cut}.C
		sed -i -e "/selection = /s/\".*$/\"${selection}\";/g" plotDataMC_${selection}_${cut}.C
		sed -i -e "/${cut}/,/;/s_//__g" plotDataMC_${selection}_${cut}.C
		sed -i -e "s,Selected/,../Selected/,g" plotDataMC_${selection}_${cut}.C

	# compile
		make -f Makefile_${selection}_${cut}

	# submit
		qsub batch_plots_${selection}_${cut}.sh

	# execute
#		./plotDataMC_${selection}_${cut} #&

	done

#cleaning
	rm batch_plots_${selection}.sh
	cd ..
done

exit 0

