#./bin/bash
fitVariable=mmg_s
fitPercentage=100
injectedResolution=0
analysisVersion=2011
for dataType in 'data' 'MC'
do
	for toy in `seq 1 100` 
	do

		echo "dataType = ${dataType}, fitVariable = ${fitVariable}, fitPercentage = ${fitPercentage}, injectedResolution = ${injectedResolution}, toy = ${toy}" 	
		qsub batchJob_muSys.sh input Results_2011_thesis_v1_MuSys ${dataType} ${fitVariable} Photon_Et LimitesAllPtFiveBin.txt ${fitPercentage} ${injectedResolution} ${toy} ${analysisVersion}
	done 
done


