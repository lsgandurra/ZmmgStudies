#./bin/bash

	for dataType in 'data' 'MC' 
	do
		for fitVariable in 'mmg_s' 
	        do
			for eta in 'Barrel' 'Endcaps'
			do
	        		for r9 in 'low' 'high' 'all'
	        		do
					##./ProfilePlotter.exe SVsNVerticesProfile ${dataType} ${eta} ${r9} nVertices ${fitVariable} Limites_nVertices_allR9_all_data_10bins.txt
					./ProfilePlotter2.exe SVsNVerticesProfile_combined_v2 ${dataType} ${eta} ${r9} nVertices ${fitVariable} Limites_nVertices_allR9_all_data_5bins.txt
				done
			done
		done
	done
