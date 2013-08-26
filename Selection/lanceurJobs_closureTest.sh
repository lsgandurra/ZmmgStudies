#./bin/bash

pileup_set=PU_S10
lumi_set=2012
isZgamma=0
lowMuMuCut=35
highMuMuCut=9999
itoy=0
for sample in 'Run2012A_13Jul2012_v1_NewMuonID' 'Run2012A_recover_06Aug2012_v1_NewMuonID' 'Run2012B_13Jul2012_v4_NewMuonID' 'Run2012C-24Aug2012-v1_NewMuonID' 'Run2012C-EcalRecover_11Dec2012-v1_NewMuonID' 'Run2012C_PromptReco_v2_NewMuonID' 'Run2012D_PromptReco_v1_NewMuonID' 
do

		if [ "$sample" = "Run2012A_13Jul2012_v1_NewMuonID" ] || [ "$sample" = "Run2012A_recover_06Aug2012_v1_NewMuonID" ] || [ "$sample" = "Run2012B_13Jul2012_v4_NewMuonID" ] || [ "$sample" = "Run2012C-24Aug2012-v1_NewMuonID" ] || [ "$sample" = "Run2012C-EcalRecover_11Dec2012-v1_NewMuonID" ] || [ "$sample" = "Run2012C_PromptReco_v2_NewMuonID" ]
		then
			pileup_set=PU_S10
			lumi_set=2012ABC
			isZgamma=0	
		fi

		if [ "$sample" = "Run2012D_PromptReco_v1_NewMuonID" ]
		then
			pileup_set=PU_S10
                        lumi_set=2012D
                        isZgamma=0
		fi
		
                if [ "$sample" = "DYToMuMu_Summer12_NewMuonID_FSR" ]
		then
                        pileup_set=PU_S10
			lumi_set=2012       
                        isZgamma=1
			sample=DYToMuMu_Summer12_NewMuonID
                fi

		if [ "$sample" = "DYToMuMu_Summer12_NewMuonID_nonFSR" ]
                then
                        pileup_set=PU_S10  
			lumi_set=2012         
                        isZgamma=2
			sample=DYToMuMu_Summer12_NewMuonID
                fi	
		
		if [ "$sample" = "TTJets_Summer12_S7_NewMuonID" ]
                then
                        pileup_set=PU_S7
			lumi_set=2012
                        isZgamma=3
                fi

		if [ "$sample" = "WJetsToLNu_Summer12_S10_NewMuonID" ]
                then
                        pileup_set=PU_S10
			lumi_set=2012
                        isZgamma=3
                fi

                for inj_resolution in '0' ##'0.5' '1' '1.5' '2' '2.5' '3' '3.5' '4'
		do
			for extraScale in '0.999909' '0.988684' '0.994057' '0.994712' '0.978923' '0.987336'
			do	
		
				echo "sample = $sample, pileup_set = ${pileup_set}, lumi_set = ${lumi_set}, isZgamma = ${isZgamma}, inj_resolution = ${inj_resolution}, extraScale = ${extraScale}"
				

				##qsub batchJob.sh ${sample} ${sample}_NewSelection_${isZgamma}_injRe${inj_resolution}_${extraScale}_v6 ${isZgamma} ${lumi_set} ${pileup_set} ${lowMuMuCut} ${highMuMuCut} MITregression ${extraScale} 3 0.0 ${inj_resolution} ${itoy} -l os=sl5 
				##hadd miniTree_${sample}_NewSelection_${isZgamma}_injRe${inj_resolution}_${extraScale}_v6_partALL.root miniTree_${sample}_NewSelection_${isZgamma}_injRe${inj_resolution}_${extraScale}_v6_part[0-9]*root
                        	##echo "miniTree_${sample}_${Zgamma2}_${SetOfCorrections2}_Triangle_v1_part[0-9].root"
				mv miniTree_${sample}_NewSelection_${isZgamma}_injRe${inj_resolution}_${extraScale}_v6_part[0-9]*root stored_miniTree/
                        	mv ${sample}_NewSelection_${isZgamma}_injRe${inj_resolution}_${extraScale}_v6_part[0-9]*err stored_errput/
                        	mv ${sample}_NewSelection_${isZgamma}_injRe${inj_resolution}_${extraScale}_v6_part[0-9]*out stored_output/
		
			done
		done
done

for sample in 'DYToMuMu_Summer12_NewMuonID_FSR' 'DYToMuMu_Summer12_NewMuonID_nonFSR' 'TTJets_Summer12_S7_NewMuonID' 'WJetsToLNu_Summer12_S10_NewMuonID'
do

		if [ "$sample" = "Run2012A_13Jul2012_v1_NewMuonID" ] || [ "$sample" = "Run2012A_recover_06Aug2012_v1_NewMuonID" ] || [ "$sample" = "Run2012B_13Jul2012_v4_NewMuonID" ] || [ "$sample" = "Run2012C-24Aug2012-v1_NewMuonID" ] || [ "$sample" = "Run2012C-EcalRecover_11Dec2012-v1_NewMuonID" ] || [ "$sample" = "Run2012C_PromptReco_v2_NewMuonID" ]
		then
			pileup_set=PU_S10
			lumi_set=2012ABC
			isZgamma=0	
		fi

		if [ "$sample" = "Run2012D_PromptReco_v1_NewMuonID" ]
		then
			pileup_set=PU_S10
                        lumi_set=2012D
                        isZgamma=0
		fi
		
                if [ "$sample" = "DYToMuMu_Summer12_NewMuonID_FSR" ]
		then
                        pileup_set=PU_S10
			lumi_set=2012       
                        isZgamma=1
			sample=DYToMuMu_Summer12_NewMuonID
                fi

		if [ "$sample" = "DYToMuMu_Summer12_NewMuonID_nonFSR" ]
                then
                        pileup_set=PU_S10  
			lumi_set=2012         
                        isZgamma=2
			sample=DYToMuMu_Summer12_NewMuonID
                fi	
		
		if [ "$sample" = "TTJets_Summer12_S7_NewMuonID" ]
                then
                        pileup_set=PU_S7
			lumi_set=2012
                        isZgamma=3
                fi

		if [ "$sample" = "WJetsToLNu_Summer12_S10_NewMuonID" ]
                then
                        pileup_set=PU_S10
			lumi_set=2012
                        isZgamma=3
                fi

                for inj_resolution in '0' ##'0.5' '1' '1.5' '2' '2.5' '3' '3.5' '4'
		do
			for extraScale in '0.997083' '0.99926' '0.998738' '0.989083' '0.995002' '0.991969'
			do	
		
				echo "sample = $sample, pileup_set = ${pileup_set}, lumi_set = ${lumi_set}, isZgamma = ${isZgamma}, inj_resolution = ${inj_resolution}, extraScale = ${extraScale}"
				

				##qsub batchJob.sh ${sample} ${sample}_NewSelection_${isZgamma}_injRe${inj_resolution}_${extraScale}_v6 ${isZgamma} ${lumi_set} ${pileup_set} ${lowMuMuCut} ${highMuMuCut} MITregression ${extraScale} 3 0.0 ${inj_resolution} ${itoy} -l os=sl5 
				##hadd miniTree_${sample}_NewSelection_${isZgamma}_injRe${inj_resolution}_${extraScale}_v6_partALL.root miniTree_${sample}_NewSelection_${isZgamma}_injRe${inj_resolution}_${extraScale}_v6_part[0-9]*root
                        	##echo "miniTree_${sample}_${Zgamma2}_${SetOfCorrections2}_Triangle_v1_part[0-9].root"
				mv miniTree_${sample}_NewSelection_${isZgamma}_injRe${inj_resolution}_${extraScale}_v6_part[0-9]*root stored_miniTree/
                        	mv ${sample}_NewSelection_${isZgamma}_injRe${inj_resolution}_${extraScale}_v6_part[0-9]*err stored_errput/
                        	mv ${sample}_NewSelection_${isZgamma}_injRe${inj_resolution}_${extraScale}_v6_part[0-9]*out stored_output/
		
			done
		done
done

