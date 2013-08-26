#./bin/bash
for dataType in 'data' 'MC'
do
	for fitVariable in 'mmg_s' 
	do
	        for fitPercentage in `seq 60 100`
		do
			for injectedResolution in '0' ##'0.5' '1' '1.5' '2' '2.5' '3' '3.5' '4' '4.5'
			do

				echo "dataType = ${dataType}, fitVariable = ${fitVariable}, fitPercentage = ${fitPercentage}" 	
				qsub batchJob_closureTestZee.sh input Results_v6_Zee_Regressionb ${dataType} ${fitVariable} Photon_Et LimitesAllPtOneBin.txt ${fitPercentage} ${injectedResolution} 0 Zee_Regression
				qsub batchJob_closureTestZee.sh input Results_v6_Zee_RecoEnergyb ${dataType} ${fitVariable} Photon_Et LimitesAllPtOneBin.txt ${fitPercentage} ${injectedResolution} 0 Zee_RecoEnergy	
			done 
		done
        done
done

##qsub batchJob_closureTest.sh input Results_v6_autoClosureTest data mmg_s Photon_Et LimitesAllPtOneBin.txt 90 0
##./SFits.exe input Results_v1 data mmg_s Photon_Et LimitesAllPtOneBin.txt Barrel all voigtian 95

