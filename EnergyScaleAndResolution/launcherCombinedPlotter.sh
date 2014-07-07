#./bin/bash

##rm -rf Results_v6_autoClosureTest/MC/Selected_Fits/
##rm -rf Results_v6_autoClosureTest/data/Selected_Fits/

##directory=Results_v6_Zee_Regressionb
directory=Results_22Jan_v1
for injectedResolution in '0' ##'0.5' '1' '1.5' '2' '2.5' '3' '3.5' '4' 
do
	##directory=Results_v6_Zee_Regressionb/InjectedResolution_${injectedResolution}Percent
	directory=Results_22Jan_v1/InjectedResolution_${injectedResolution}Percent
	echo "${directory}"
	##rm -rf ${directory}/MC/Selected_Fits/
	##rm -rf ${directory}/data/Selected_Fits/
	##echo "rm -rf ${directory}/MC/Selected_Fits/"
	##echo "rm -rf ${directory}/data/Selected_Fits/"
	
	for dataType in 'data' 'MC' 
	do
		for fitVariable in 'mmg_s' ##'mmg_s' ##'mmg_s_true'
	        do
			if [ "$dataType" = "data" ] && [ "$fitVariable" = "mmg_s_true" ]
	                then
	                        continue
	                fi
			
			for eta in 'Barrel' 'Endcaps'
			do
	        		for r9 in 'low' 'high' 'all'
	        		do
	               			for fitFunction in 'voigtian' ##'cruijff'
	               			do	
						./PvaluesPlotter.exe ${directory} ${dataType} ${fitVariable} Photon_Et ${eta} ${r9} ${fitFunction} 60 70
						./PvaluesPlotter.exe ${directory} ${dataType} ${fitVariable} Photon_Et ${eta} ${r9} ${fitFunction} 70 80
						./PvaluesPlotter.exe ${directory} ${dataType} ${fitVariable} Photon_Et ${eta} ${r9} ${fitFunction} 80 90
						./PvaluesPlotter.exe ${directory} ${dataType} ${fitVariable} Photon_Et ${eta} ${r9} ${fitFunction} 90 101	
						./Pvalues_vs_percent.exe ${directory} ${dataType} ${fitVariable} Photon_Et ${eta} ${r9} ${fitFunction} 60 100	
						./Mean_vs_percent.exe ${directory} ${dataType} ${fitVariable} Photon_Et ${eta} ${r9} ${fitFunction} 60 100
						./Fit_Selection.exe ${directory} ${dataType} ${fitVariable} ${eta} ${r9} ${fitFunction} 60 100
					done
					##if [ "$fitVariable" = "mmg_s" ] 
					##then 
					##	./fitFunctionSystematics.exe ${directory} ${dataType} ${fitVariable} ${eta} ${r9} voigtian cruijff	
					##fi
					##./ProfilePlotter.exe SVsNVerticesProfile ${dataType} ${eta} ${r9} nVertices ${fitVariable} Limites_nVertices_allR9_all_data_10bins.txt
				done
			done
		done
	done
done
