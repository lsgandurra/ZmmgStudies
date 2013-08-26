#./bin/bash

directoryName=Results_v4
fitFunction=voigtian
rm -rf ${directoryName}/CombinedGraphs/
for dataType in 'data' 'MC'
do
	for fitVariable in 'mmg_s' ##'mmg_s_true'
        do
                if [ "$dataType" = "data" ] && [ "$fitVariable" = "mmg_s_true" ]
                then
                        continue
                fi
		for eta in 'Barrel' 'Endcaps' 
		do
			for r9 in 'low' 'high' 'all'
        		do
				./Mean_vs_injectedResolution.exe ${directoryName} ${dataType} ${fitVariable} ${eta} ${r9} ${fitFunction}
			done
		done
	done
done
