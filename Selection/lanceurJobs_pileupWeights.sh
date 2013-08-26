#./bin/bash

for sample in 'DYToMuMu_Summer12_NewMuonID' 'TTJets_Summer12_S7_NewMuonID' 'WJetsToLNu_Summer12_S10_NewMuonID'  
do


	##qsub batch_GE_PileUp.sh ${sample}  
	##hadd miniTree_pileUp_${sample}_partALL.root miniTree_pileUp_${sample}_part[0-9]*root
	mv miniTree_pileUp_${sample}_part[0-9]*root stored_miniTree/
        mv pileUp_${sample}_part[0-9]*err stored_errput/
        mv pileUp_${sample}_part[0-9]*out stored_output/

done

