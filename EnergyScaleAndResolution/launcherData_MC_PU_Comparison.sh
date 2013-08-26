#./bin/bash

cutVariable=Photon_Et
directoryName=Data_MC_PU_Comparison
for cutVariableValue in '0' '25'
do
	for eta in 'Barrel' 'Endcaps' 'all'
	do
		for r9 in 'low' 'high' 'all'
	        do
			if [ "$eta" = "all" ] && [ "$r9" = "low" ]
			then 
				continue
			fi

			if [ "$eta" = "all" ] && [ "$r9" = "high" ]
                        then 
                                continue
                        fi	

			for norm in 'integral' 
			do
				for xVariable in 'nVertices' 'weight_pileUp'
				do  
				
					for log in '0' '1'
					do
						./Data_MC_PU_Comparison.exe ${directoryName} ${eta} ${r9} ${xVariable} ${log} ${norm} ${cutVariable} ${cutVariableValue}	
					done
				done 
			done
		done
	done
done

