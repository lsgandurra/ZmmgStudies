#./bin/bash

directory=Results_v6_RecoEnergy
for injectedResolution in '0' 
do
	directory=/sps/cms/sgandurr/CMSSW_5_3_7_RECO_5_3_3_v4/src/Energy_scale_extraction/Results_v6_RecoEnergy/InjectedResolution_${injectedResolution}Percent
	rm ${directory}/MC/Selected_Fits/FitFunctionSystematics/Summary_fitFunctionSystematics.txt
	rm ${directory}/data/Selected_Fits/FitFunctionSystematics/Summary_fitFunctionSystematics.txt
	echo "rm ${directory}/MC/Selected_Fits/FitFunctionSystematics/Summary_fitFunctionSystematics.txt"
	echo "rm ${directory}/data/Selected_Fits/FitFunctionSystematics/Summary_fitFunctionSystematics.txt"
	
	for dataType in 'data' 'MC' ##`seq 0 1`
	do
		for fitVariable in 'mmg_s' 
	        do
			for eta in 'Barrel' 'Endcaps'
			do
	        		for r9 in 'low' 'high' 'all'
	        		do
					##./fitFunctionSystematics.exe ${directory} ${dataType} ${fitVariable} ${eta} ${r9} voigtian cruijff	
					qsub batchJob_fitFunctionSystematics.sh ${directory} ${dataType} ${fitVariable} ${eta} ${r9} voigtian cruijff
				done
			done
		done
	done
done
##qsub batchJob_fitFunctionSystematics.sh /sps/cms/sgandurr/CMSSW_5_3_7_RECO_5_3_3_v4/src/Energy_scale_extraction/Results_v6_RecoEnergy/InjectedResolution_0Percent MC mmg_s Endcaps low voigtian cruijff
