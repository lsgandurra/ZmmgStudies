#./bin/bash
# Small script to run samples on batch
# Written by Olivier Bondu (March 2010)

DATE=`date +%F`

for lumi in `echo "May10 Promptv4 July05 Aug05 Oct03 2011A 2011B 2011"`
#for lumi in `echo "2011A"`
do
	NAME=`echo "S6_${lumi}_40_80"`
	LOC_REP=`echo "${DATE}-Plots_${NAME}"`
	mkdir ${LOC_REP}
	for img in `echo "C png gif pdf eps"`
	do
  	mkdir ${LOC_REP}/${img}
	done
	sed -e "s/NAME/${NAME}/g" -e "s/LUMI/${lumi}/g" batch_GE_plotter_template.sh > batch_GE_plotter_${NAME}.sh
	qsub batch_GE_plotter_${NAME}.sh ${lumi} ${LOC_REP}
	mv batch_GE_plotter_${NAME}.sh stored_batch/
done

return 0;

