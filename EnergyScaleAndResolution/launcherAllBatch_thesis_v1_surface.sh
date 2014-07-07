#./bin/bash
directory=Results_thesis_v1
for dataType in 'data' 'MC'
do
	for fitVariable in 'mmg_s_true' 'mmg_s_MZ_Surface'
	do
		if [ "$dataType" = "data" ] && [ "$fitVariable" = "mmg_s_true" ]
                then
			continue
                fi

	        for fitPercentage in `seq 60 100`
		do
			for injectedResolution in '0' ##'0.5' '1' '1.5' '2' '2.5' '3' '3.5' '4' '4.5'
			do
				for analysisVersion in '2012'
				do
					if [ "$analysisVersion" = "2011" ]
					then 
						directory=Results_2011_thesis_v1_recoEnergy_surface			
					fi
					
					if [ "$analysisVersion" = "2012" ]
                                        then 
                                                directory=Results_2012_thesis_v1_recoEnergy_surface_10GeV        
                                        fi	
			
					echo "dataType = ${dataType}, fitVariable = ${fitVariable}, fitPercentage = ${fitPercentage}" 	
					qsub batchJob.sh input ${directory}  ${dataType} ${fitVariable} Photon_Et LimitesAllPtFiveBin.txt ${fitPercentage} ${injectedResolution} ${analysisVersion} -l os=sl6
				done
			done 
		done
        done
done


##qsub batchJob.sh input Results_v6 data mmg_s Photon_Et LimitesAllPtOneBin.txt 90 0
##./SFits.exe input Results_v1 data mmg_s Photon_Et LimitesAllPtOneBin.txt Barrel all voigtian 95

