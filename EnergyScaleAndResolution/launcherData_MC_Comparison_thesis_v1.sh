#./bin/bash

##rm -rf Data_MC_Pt_Comparison/
##rm -rf Data_MC_r9_Comparison/
##rm -rf Data_MC_Mmumugamma_Comparison/

cutVariable=Photon_Et
directoryName=Data_MC_Comparison
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

			for norm in 'lumi' 'integral' ##'lumi2' 'integral2'
			do
				if [ "$norm" = "lumi" ] || [ "$norm" = "integral" ]
	                	then
					directoryName=Data_MC_Comparison_thesis_v1
				fi	

		
				if [ "$norm" = "lumi2" ] || [ "$norm" = "integral2" ]
	                        then
	                                directoryName=Data_MC_Comparison_thesis_v1
	                        fi	
	
				##for xVariable in 'Photon_Et' 'Mmumugamma' 'Photon_r9' 'Photon_SC_Eta' 'nVertices' 'MuonM_Eta' 'MuonP_Eta'
				##do  
				
				##	for log in '0' '1'
				##	do
						##./Data_MC_Var_Comparison.exe ${directoryName} ${eta} ${r9} ${xVariable} ${log} ${norm} ${cutVariable} ${cutVariableValue}	
				for analysisVersion in '2012' ##'2011' '2012' 
				do

						if [ "$analysisVersion" = "2011" ] 
                                                then
                                                        directoryName=Data_MC_Comparison_2011_thesis_v1
                                                fi
      
                                                if [ "$analysisVersion" = "2012" ] 
                                                then
                                                        directoryName=Data_MC_Comparison_2012_thesis_v1_MmumugammaTest
                                                fi

						echo "${directoryName} ${eta} ${r9} ${norm} ${cutVariable} ${cutVariableValue}" 
						qsub batchJob_Data_MC_Comparison.sh ${directoryName} ${eta} ${r9} ${norm} ${cutVariable} ${cutVariableValue} ${analysisVersion} -l os=sl6

				##	done
				done 
			done
		done
	done
done
##./Data_MC_Var_Comparison.exe Data_MC_Comparison_testNewPU all Barrel nVertices 0 lumi Photon_Et 25
