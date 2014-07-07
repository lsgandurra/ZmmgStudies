#./bin/bash

pileup_set=RD1
lumi_set=2012_22Jan
isZgamma=0
lowMuMuCut=35
highMuMuCut=9999
itoy=0
scCorrection=MITregression5
analysisVersion=2012
##for sample in 'Run2012A_13Jul2012_v1_NewMuonID' 'Run2012A_recover_06Aug2012_v1_NewMuonID' 'Run2012B_13Jul2012_v4_NewMuonID' 'Run2012C-24Aug2012-v1_NewMuonID' 'Run2012C-EcalRecover_11Dec2012-v1_NewMuonID' 'Run2012C_PromptReco_v2_NewMuonID' 'Run2012D_PromptReco_v1_NewMuonID' 'DYToMuMu_Summer12_NewMuonID_FSR' 'DYToMuMu_Summer12_NewMuonID_nonFSR' 'TTJets_Summer12_S7_NewMuonID' 'WJetsToLNu_Summer12_S10_NewMuonID'  
##for sample in 'totouples_Run2012A_22Jan2013_v1_November2013' 'totouples_parked_Run2012B_22Jan2013_v1_November2013' 'totouples_parked_Run2012C_22Jan2013_v1_November2013' 'totouples_parked_Run2012D_22Jan2013_v1_November2013' 'totouples_DYToMuMu_Summer12_November2013_FSR' 'totouples_DYToMuMu_Summer12_November2013_nonFSR' 'totouples_WJetsToLNu_Summer12_S10_November2013' 'totouples_TTJets_SemiLeptMGDecays_Summer12_November2013' 'totouples_TTJets_FullLeptMGDecays_Summer12_November2013_v2'  
##for sample in 'totouples_DYToMuMu_Summer12_February2014_noskim_FSR' 'totouples_DYToMuMu_Summer12_February2014_noskim_nonFSR'
for sample in 'totouples_DYToMuMu_Summer12_November2013_FSR' 'totouples_DYToMuMu_Summer12_November2013_nonFSR' 'totouples_TTJets_FullLeptMGDecays_Summer12_November2013_v2' 'totouples_TTJets_SemiLeptMGDecays_Summer12_November2013' 'totouples_TTJets_HadronicMGDecays_Summer12_November2013_v2_temp218'  
##for sample in 'totouples_DYToMuMu_Summer12_February2014_noskim_FSR'
##for sample in 'totouples_TTJets_FullLeptMGDecays_Summer12_November2013_v2' 'totouples_TTJets_SemiLeptMGDecays_Summer12_November2013' 'totouples_TTJets_HadronicMGDecays_Summer12_November2013_v2_temp218'
##for sample in 'totouples_DYToMuMu_Summer12_February2014_noskim_FSR' 'totouples_Run2012A_22Jan2013_v1_November2013' 'totouples_parked_Run2012B_22Jan2013_v1_November2013' 'totouples_parked_Run2012C_22Jan2013_v1_November2013' 'totouples_parked_Run2012D_22Jan2013_v1_November2013'
do

		if [ "$sample" = "totouples_Run2012A_22Jan2013_v1_November2013" ]
                then
                        pileup_set=RD1
                        lumi_set=2012A_22Jan
                        isZgamma=0
                fi

		if [ "$sample" = "totouples_parked_Run2012B_22Jan2013_v1_November2013" ]
                then
                        pileup_set=RD1
                        lumi_set=2012B_22Jan
                        isZgamma=0
                fi

		if [ "$sample" = "totouples_parked_Run2012C_22Jan2013_v1_November2013" ]
                then
                        pileup_set=RD1
                        lumi_set=2012C_22Jan
                        isZgamma=0
                fi

		if [ "$sample" = "totouples_parked_Run2012D_22Jan2013_v1_November2013" ]
		then
			pileup_set=RD1
                        lumi_set=2012D_22Jan
                        isZgamma=0
		fi
		
                if [ "$sample" = "totouples_DYToMuMu_Summer12_February2014_noskim_FSR" ]
		then
                        pileup_set=RD1
			lumi_set=2012_22Jan      
                        isZgamma=1
			sample=totouples_DYToMuMu_Summer12_February2014_noskim
			scCorrection=MITregression5
			analysisVersion=2012
                fi

		if [ "$sample" = "totouples_DYToMuMu_Summer12_February2014_noskim_nonFSR" ]
                then
                        pileup_set=RD1
			lumi_set=2012_22Jan        
                        isZgamma=2
			sample=totouples_DYToMuMu_Summer12_February2014_noskim
			scCorrection=MITregression5
			analysisVersion=2012
                fi	
	
		if [ "$sample" = "totouples_DYToMuMu_Summer12_November2013_FSR" ]
                then
                        pileup_set=RD1
                        lumi_set=2012_22Jan    
                        isZgamma=1
                        sample=totouples_DYToMuMu_Summer12_November2013
                        scCorrection=MITregression5
                        analysisVersion=2012
                fi

                if [ "$sample" = "totouples_DYToMuMu_Summer12_November2013_nonFSR" ]
                then
                        pileup_set=RD1
                        lumi_set=2012_22Jan    
                        isZgamma=2
                        sample=totouples_DYToMuMu_Summer12_November2013
                        scCorrection=MITregression5
                        analysisVersion=2012
                fi

	
		if [ "$sample" = "DYToMuMu_Summer12_PU_S10_noskim_FSR" ]
                then
                        pileup_set=PU_S10
                        lumi_set=2012_22Jan    
                        isZgamma=1
                        sample=DYToMuMu_Summer12_PU_S10_noskim
			scCorrection=MITregressionOld
			analysisVersion=2012
                fi

                if [ "$sample" = "DYToMuMu_Summer12_PU_S10_noskim_nonFSR" ]
                then
                        pileup_set=PU_S10
                        lumi_set=2012_22Jan    
                        isZgamma=2
                        sample=DYToMuMu_Summer12_PU_S10_noskim
			scCorrection=MITregressionOld
			analysisVersion=2012
                fi 	
	

		if [ "$sample" = "totouples_TTJets_SemiLeptMGDecays_Summer12_November2013" ] || [ "$sample" = "totouples_TTJets_FullLeptMGDecays_Summer12_November2013_v2" ] || [ "$sample" = "totouples_TTJets_HadronicMGDecays_Summer12_November2013_v2_temp218" ]
                then
                        pileup_set=RD1
                        lumi_set=2012_22Jan
                        isZgamma=3
			scCorrection=MITregression5
			analysisVersion=2012
                fi


		if [ "$sample" = "totouples_WJetsToLNu_Summer12_S10_November2013" ]
                then
                        pileup_set=PU_S10
			lumi_set=2012_22Jan
                        isZgamma=3
                fi

		if [ "$sample" = "TTJets_Summer12_S7_NewMuonID" ]
                then
                        pileup_set=PU_S7
                        lumi_set=2012_22Jan
                        isZgamma=3
			scCorrection=MITregressionOld
                fi

                for inj_resolution in '0' ##'0.5' '1' '1.5' '2' '2.5' '3' '3.5' '4'
		do

				echo "sample = $sample, pileup_set = ${pileup_set}, lumi_set = ${lumi_set}, isZgamma = ${isZgamma}, inj_resolution = ${inj_resolution}"
				

				##qsub batchJob.sh ${sample} ${sample}_${isZgamma}_RD1_thesis_v1 ${isZgamma} ${lumi_set} ${pileup_set} ${lowMuMuCut} ${highMuMuCut} ${scCorrection} 1.0 3 0.0 ${inj_resolution} ${itoy} ${analysisVersion} -l os=sl6 
				##hadd miniTree_${sample}_${isZgamma}_RD1_thesis_v1_partALL.root miniTree_${sample}_${isZgamma}_RD1_thesis_v1_part[0-9]*root
				##hadd miniFriend2_${sample}_${isZgamma}_RD1_thesis_v1_partALL.root miniFriend2_${sample}_${isZgamma}_RD1_thesis_v1_part[0-9]*root	
				##hadd miniFriend3_${sample}_${isZgamma}_RD1_thesis_v1_partALL.root miniFriend3_${sample}_${isZgamma}_RD1_thesis_v1_part[0-9]*root
				##hadd miniFriend4_${sample}_${isZgamma}_RD1_thesis_v1_partALL.root miniFriend4_${sample}_${isZgamma}_RD1_thesis_v1_part[0-9]*root
				mv miniTree_${sample}_${isZgamma}_RD1_thesis_v1_part[0-9]*root stored_miniTree/
                        	##mv miniFriend2_${sample}_${isZgamma}_RD1_thesis_v1_part[0-9]*root stored_miniTree/
				##mv miniFriend3_${sample}_${isZgamma}_RD1_thesis_v1_part[0-9]*root stored_miniTree/
				##mv miniFriend4_${sample}_${isZgamma}_RD1_thesis_v1_part[0-9]*root stored_miniTree/
				mv ${sample}_${isZgamma}_RD1_thesis_v1_part[0-9]*err stored_errput/
                        	mv ${sample}_${isZgamma}_RD1_thesis_v1_part[0-9]*out stored_output/
		done

done

