#./bin/bash

##rm -rf Data_MC_Pt_Comparison/
##rm -rf Data_MC_r9_Comparison/
##rm -rf Data_MC_Mmumugamma_Comparison/

##cutVariable=iRunID4
cutVariable=NoCut
directoryName=RD1_S10_Comparison_before_muons_selection_v7_nbSuperClusterSup25
##directoryName=RD1_S10_Comparison_after_dimuons_selection
for cutVariableValue in '0'
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
	
				##for xVariable in 'Photon_r9_4' 'Photon_Eta_4' 'Photon_Phi_4' 'Photon_SC_Eta_4' 'Photon_SC_Phi_4' 'PhotonMC_Eta_4' 'PhotonMC_Phi_4' 'Photon_SC_E_4' 'Photon_SC_Eraw_4' 'Photon_correctedE_4' 'Photon_E_4' 'PhotonMC_E_4' 'PhotonMC_Rconv_4'  
				##for xVariable in 'muon_Pt_2' 'muon_E_2' 'muon_Eta_2' 'muon_Phi_2' 'muon_DeltaR_2'
				##for xVariable in 'Photon_Et_2' 'diMuonPt_2' 'Photon_SC_E_2' 	
				##for xVariable in 'ParticleMC_E_2-Photon_E_2' 'ParticleMC_E_2'
				for xVariable in 'Photon_SC_E_2'
				do  
				
					for log in '0' '1'
					do
						./RD1_S10_Comparison.exe ${directoryName} ${cutVariable} ${cutVariableValue} ${xVariable} ${eta} ${r9} ${log} 
						##./RD1_S10_Comparison_2.exe ${directoryName} ${cutVariable} ${cutVariableValue} ${xVariable} ${eta} ${r9} ${log}
					done
				done 
		done
	done
done
##./Data_MC_Var_Comparison.exe Data_MC_Comparison_testNewPU all Barrel nVertices 0 lumi Photon_Et 25
