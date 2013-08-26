#./bin/bash
for dataType in 'data' 'MC'
do
	for fitVariable in 'mmg_s' 'mmg_s_true'
	do
		if [ "$dataType" = "data" ] && [ "$fitVariable" = "mmg_s_true" ]
                then
			continue
                fi

	        for fitPercentage in `seq 60 100`
		do
			for injectedResolution in '0' ##'0.5' '1' '1.5' '2' '2.5' '3' '3.5' '4' '4.5'
			do

				echo "dataType = ${dataType}, fitVariable = ${fitVariable}, fitPercentage = ${fitPercentage}" 	
				qsub batchJob.sh input Results_v6_RecoEnergy ${dataType} ${fitVariable} Photon_Et LimitesAllPtOneBin.txt ${fitPercentage} ${injectedResolution}
			done 
		done
        done
done


##qsub batchJob.sh input Results_v6 data mmg_s Photon_Et LimitesAllPtOneBin.txt 90 0
##./SFits.exe input Results_v1 data mmg_s Photon_Et LimitesAllPtOneBin.txt Barrel all voigtian 95

