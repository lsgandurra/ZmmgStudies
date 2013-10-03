#./bin/bash

pileup_set=PU_S6
isZgamma=0
lowMuMuCut=40
highMuMuCut=80
lumiset=2011
itoy=0
scCorrection=1.0 ##FIXME MITregression
for sample in '2011A_03Oct2011V1ReReco_toto_v2' '2011A_05Jul2011ReReco_toto_v2' '2011A_PromptSkimV5ReReco_toto_v2' '2011B_PromptSkimV1ReReco_toto_v2' 'DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Fall11-PU_S6_START42_V14B-v1_September12_FSR' 'DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Fall11-PU_S6_START42_V14B-v1_September12_nonFSR' 'TTJets_TuneZ2_7TeV-madgraph-tauola' 'WJetsToLNu_TuneZ2_7TeV-madgraph-tauola'
do

		if [ "$sample" = "2011A_03Oct2011V1ReReco_toto_v2" ] || [ "$sample" = "2011A_05Jul2011ReReco_toto_v2" ] || [ "$sample" = "2011A_PromptSkimV5ReReco_toto_v2" ] 
		then
			pileup_set=PU_S6
			isZgamma=0
			lumiset=2011A	
		fi

		if [ "$sample" = "2011B_PromptSkimV1ReReco_toto_v2" ] 
                then
                        pileup_set=PU_S6
                        isZgamma=0   
			lumiset=2011B  
                fi	
		
                if [ "$sample" = "DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Fall11-PU_S6_START42_V14B-v1_September12_FSR" ]
		then
                        pileup_set=PU_S6       
                        isZgamma=1
			sample=DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Fall11-PU_S6_START42_V14B-v1_September12
                	lumiset=2011
		fi

		if [ "$sample" = "DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Fall11-PU_S6_START42_V14B-v1_September12_nonFSR" ]
                then
                        pileup_set=PU_S6           
                        isZgamma=2
			sample=DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Fall11-PU_S6_START42_V14B-v1_September12
                	lumiset=2011
		fi	
		
		if [ "$sample" = "TTJets_TuneZ2_7TeV-madgraph-tauola" ]
                then
                        pileup_set=PU_S6
                        isZgamma=3
			lumiset=2011
                fi

		if [ "$sample" = "WJetsToLNu_TuneZ2_7TeV-madgraph-tauola" ]
                then
                        pileup_set=PU_S6
                        isZgamma=3
			lumiset=2011
                fi
		for inj_resolution in '0' ##'0.5' '1' '1.5' '2' '2.5' '3' '3.5' '4'
		do
				echo "sample = $sample, pileup_set = ${pileup_set}, isZgamma = ${isZgamma}, lumiset=${lumiset}, inj_resolution=${inj_resolution}"
				

				qsub batchJob.sh ${sample} ${sample}_${isZgamma}_injRe${inj_resolution}_v1 ${isZgamma} ${lumiset} ${pileup_set} ${lowMuMuCut} ${highMuMuCut} ${scCorrection} 1.0 3 0.0 ${inj_resolution} ${itoy} 
				##hadd miniTree_${sample}_${isZgamma}_injRe${inj_resolution}_v1_partALL.root miniTree_${sample}_${isZgamma}_injRe${inj_resolution}_v1_part[0-9]*root
                        	##echo "miniTree_${sample}_${Zgamma2}_${SetOfCorrections2}_Triangle_v1_part[0-9].root"
				##mv miniTree_${sample}_${isZgamma}_injRe${inj_resolution}_v1_part[0-9]*root stored_miniTree/
                        	##mv ${sample}_${isZgamma}_injRe${inj_resolution}_v1_part[0-9]*err stored_errput/
                        	##mv ${sample}_${isZgamma}_injRe${inj_resolution}_v1_part[0-9]*out stored_output/

		done
done

