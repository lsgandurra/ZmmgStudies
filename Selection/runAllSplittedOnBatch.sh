#./bin/bash
# Small script to run samples on batch
# Written by Olivier Bondu (March 2010)

for zgamma in `seq 1 2`
#for zgamma in `seq 2 2`
do
#	for lumi in `echo "May10 Promptv4 July05 Aug05 Oct03 2011A 2011B 2011"`
#	for lumi in `echo "2011A 2011B 2011"`
	for lumi in `echo "2011"`
	do
		low[1]=38
		high[1]=82
		#
		low[2]=36
		high[2]=84
		#
		low[3]=38
		high[3]=84
		#
		low[4]=38
		high[4]=86
		#
		low[5]=38
		high[5]=88
#		for idx in `seq 1 5`
#		for low in `echo "36 44"`
#		for low in `echo "38 40 42"`
#		for low in `echo "40 42"`
		for low in `echo "40"`
		do
			high=`echo "${low} + 40" | bc -ql`
#			low=${low[${idx}]}
#			high=${high[${idx}]}
			if [[ "${zgamma}" == "1" ]]
			then
				zg_disp="FSR"
			fi
			if [[ "${zgamma}" == "2" ]]
			then
				zg_disp="nonFSR"
			fi
			name=`echo "v10_${zg_disp}_DYToMuMu_S6_${lumi}_${low}_${high}"`
#			name=`echo "v10_${zg_disp}_DYToMuMu_S6_${lumi}_${low}"`
			echo ${name}
#
			sed -e "s/NAME/${name}/g" batch_GE_selection_template.sh > batch_GE_selection_${name}.sh
			qsub batch_GE_selection_${name}.sh Reduc_v03_DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_PU_S6_v2 ${name} ${zgamma} ${lumi} PU_S6 ${low} ${high} 
			mv batch_GE_selection_${name}.sh stored_batch/
#
#			qdel ${name}
#
#			hadd miniTree_${name}_partALL.root miniTree_${name}_part[0-9]*root
#			mv miniTree_${name}_part[0-9]*root stored_miniTree/
#			mv ${name}_part[0-9]*err stored_errput/
#			mv ${name}_part[0-9]*out stored_output/
#			mv ${name}.o[0-9]* stored_logGE/
		done
	done
done

return 1;

for lumi in `echo "2011"`
#for lumi in `echo "2011A 2011B 2011"`
#for lumi in `echo "May10 Promptv4 July05 Aug05 Oct03 2011A 2011B 2011"`
do
	low[1]=38
	high[1]=82
	#
	low[2]=36
	high[2]=84
	#
	low[3]=38
	high[3]=84
	#
	low[4]=38
	high[4]=86
	#
	low[5]=38
	high[5]=88
	for idx in `seq 1 5`
#	for low in `echo "38 40 42"`
#	for low in `echo "40"`
	do
		low=${low[${idx}]}
		high=${high[${idx}]}
#		hadd miniTree_v10_DYToMuMu_S6_${lumi}_${low}_partALL.root miniTree_v10_FSR_DYToMuMu_S6_${lumi}_${low}_partALL.root miniTree_v10_nonFSR_DYToMuMu_S6_${lumi}_${low}_partALL.root
		hadd miniTree_v10_DYToMuMu_S6_${lumi}_${low}_${high}_partALL.root miniTree_v10_FSR_DYToMuMu_S6_${lumi}_${low}_${high}_partALL.root miniTree_v10_nonFSR_DYToMuMu_S6_${lumi}_${low}_${high}_partALL.root
	done
done



return 0;


