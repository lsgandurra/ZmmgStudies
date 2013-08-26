#./bin/bash
fitVariable=mmg_s
fitPercentage=100
injectedResolution=0.00
for dataType in 'data' 'MC'
do
	for toy in `seq 1 100` 
	do

		echo "dataType = ${dataType}, fitVariable = ${fitVariable}, fitPercentage = ${fitPercentage}, injectedResolution = ${injectedResolution}, toy = ${toy}" 	
		qsub batchJob_muSys.sh input Results_v6_RecoEnergy_MuSys ${dataType} ${fitVariable} Photon_Et LimitesAllPtOneBin.txt ${fitPercentage} ${injectedResolution} ${toy}
	done 
done


