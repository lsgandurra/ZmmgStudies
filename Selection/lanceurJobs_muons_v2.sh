#./bin/bash

pileup_set=PU_S10
lumi_set=2012
isZgamma=0
lowMuMuCut=35
highMuMuCut=9999
itoy=0
analysisVersion=2011
scCorrection=MITregression2
##for sample in 'totouples_Run2012A_22Jan2013_v1_November2013' 'totouples_parked_Run2012B_22Jan2013_v1_November2013' 'totouples_parked_Run2012C_22Jan2013_v1_November2013' 'totouples_parked_Run2012D_22Jan2013_v1_November2013' 'totouple_DYToMuMu_Summer12_S10_reg5' 'totouples_WJetsToLNu_Summer12_S10_November2013' 'totouple_TTJets_Summer12_S10_reg5' 'totouples_Run2011A_12Oct2013_v1' 'totouples_Run2011B_12Oct2013_v1' 'totouple_DYToMuMu_Summer11_S13_v2' 'totouple_WJetsToLNu_Summer11_S13' 'totouple_TTJets_Summer11_S13'
##for sample in 'totouple_DYToMuMu_Summer12_S10_reg5_noSkim_v1' 'totouple_TTJets_Summer12_S10_reg5_noSkim_v1' 'totouples_Run2012A_22Jan2013_v1_noSkim_v1' 'totouples_parked_Run2012B_22Jan2013_v1_noSkim_v1' 'totouples_parked_Run2012C_22Jan2013_v1_noSkim_v1' 'totouples_parked_Run2012D_22Jan2013_v1_noSkim_v1' 'totouples_WJetsToLNu_Summer12_S10_noSkim_v1'
##for sample in 'totouples_parked_Run2012B_22Jan2013_v1_noSkim_v1' 'totouples_parked_Run2012C_22Jan2013_v1_noSkim_v1' 'totouples_parked_Run2012D_22Jan2013_v1_noSkim_v1' 'totouples_WJetsToLNu_Summer12_S10_noSkim_v1'
##for sample in 'totouples_parked_Run2012B_22Jan2013_v1_noSkim_v1_part1' 'totouples_parked_Run2012B_22Jan2013_v1_noSkim_v1_part2'
##for sample in 'totouples_parked_Run2012C_22Jan2013_v1_noSkim_v1_part1' 'totouples_parked_Run2012C_22Jan2013_v1_noSkim_v1_part2' 'totouples_parked_Run2012D_22Jan2013_v1_noSkim_v1_part1' 'totouples_parked_Run2012D_22Jan2013_v1_noSkim_v1_part2' 'totouples_WJetsToLNu_Summer12_S10_noSkim_v1_part1' 'totouples_WJetsToLNu_Summer12_S10_noSkim_v1_part2'
##for sample in 'totouples_WJetsToLNu_Summer12_S10_noSkim_v1_part1' 'totouples_WJetsToLNu_Summer12_S10_noSkim_v1_part2' 'totouples_parked_Run2012B_22Jan2013_v1_noSkim_v1_part1' 'totouples_parked_Run2012B_22Jan2013_v1_noSkim_v1_part2'
##for sample in 'totouples_parked_Run2012C_22Jan2013_v1_noSkim_v1_part1' 'totouples_parked_Run2012C_22Jan2013_v1_noSkim_v1_part2' 'totouples_parked_Run2012C_22Jan2013_v1_noSkim_v1_part3' 'totouples_parked_Run2012D_22Jan2013_v1_noSkim_v1_part1' 'totouples_parked_Run2012D_22Jan2013_v1_noSkim_v1_part2' 'totouples_parked_Run2012D_22Jan2013_v1_noSkim_v1_part3'
##for sample in 'totouple_DYToMuMu_Summer12_S10_reg5_noSkim_v1' 'totouple_TTJets_Summer12_S10_reg5_noSkim_v1' 'totouples_WJetsToLNu_Summer12_S10_noSkim_v1_part1' 'totouples_WJetsToLNu_Summer12_S10_noSkim_v1_part2' 'totouples_WJetsToLNu_Summer12_S10_noSkim_v1_part3' 'totouples_Run2012A_22Jan2013_v1_noSkim_v1' 'totouples_parked_Run2012B_22Jan2013_v1_noSkim_v1_part1' 'totouples_parked_Run2012B_22Jan2013_v1_noSkim_v1_part2' 'totouples_parked_Run2012C_22Jan2013_v1_noSkim_v1_part1' 'totouples_parked_Run2012C_22Jan2013_v1_noSkim_v1_part2' 'totouples_parked_Run2012C_22Jan2013_v1_noSkim_v1_part3' 'totouples_parked_Run2012D_22Jan2013_v1_noSkim_v1_part1' 'totouples_parked_Run2012D_22Jan2013_v1_noSkim_v1_part2' 'totouples_parked_Run2012D_22Jan2013_v1_noSkim_v1_part3'
for sample in 'totouple_DYToMuMu_Summer12_S10_reg5_noSkim_v1' 'totouples_Run2012A_22Jan2013_v1_noSkim_v1' 'totouples_parked_Run2012B_22Jan2013_v1_noSkim_v1_part1' 'totouples_parked_Run2012B_22Jan2013_v1_noSkim_v1_part2' 'totouples_parked_Run2012C_22Jan2013_v1_noSkim_v1_part1' 'totouples_parked_Run2012C_22Jan2013_v1_noSkim_v1_part2' 'totouples_parked_Run2012C_22Jan2013_v1_noSkim_v1_part3' 'totouples_parked_Run2012D_22Jan2013_v1_noSkim_v1_part1' 'totouples_parked_Run2012D_22Jan2013_v1_noSkim_v1_part2' 'totouples_parked_Run2012D_22Jan2013_v1_noSkim_v1_part3'
##for sample in 'totouples_WJetsToLNu_Summer12_S10_noSkim_v1_part1' 'totouples_WJetsToLNu_Summer12_S10_noSkim_v1_part2' 'totouples_WJetsToLNu_Summer12_S10_noSkim_v1_part3'
##for sample in 'totouples_parked_Run2012B_22Jan2013_v1_noSkim_v1_part1' 'totouples_parked_Run2012B_22Jan2013_v1_noSkim_v1_part2' 'totouples_parked_Run2012C_22Jan2013_v1_noSkim_v1_part1' 'totouples_parked_Run2012C_22Jan2013_v1_noSkim_v1_part2'
do

		if [ "$sample" = "totouples_Run2012A_22Jan2013_v1_noSkim_v1" ]
		then
			pileup_set=PU_S10
                        lumi_set=2012A_22Jan
                        isZgamma=0
			analysisVersion=2012
			scCorrection=MITregression5
		fi
		
		if [ "$sample" = "totouples_parked_Run2012B_22Jan2013_v1_noSkim_v1" ] || [ "$sample" = "totouples_parked_Run2012B_22Jan2013_v1_noSkim_v1_part1" ] || [ "$sample" = "totouples_parked_Run2012B_22Jan2013_v1_noSkim_v1_part2" ]
                then
                        pileup_set=PU_S10
                        lumi_set=2012B_22Jan
                        isZgamma=0
                        analysisVersion=2012
			scCorrection=MITregression5
                fi

		if [ "$sample" = "totouples_parked_Run2012C_22Jan2013_v1_noSkim_v1" ] || [ "$sample" = "totouples_parked_Run2012C_22Jan2013_v1_noSkim_v1_part1" ] || [ "$sample" = "totouples_parked_Run2012C_22Jan2013_v1_noSkim_v1_part2" ] || [ "$sample" = "totouples_parked_Run2012C_22Jan2013_v1_noSkim_v1_part3" ]
                then
                        pileup_set=PU_S10
                        lumi_set=2012C_22Jan
                        isZgamma=0
                        analysisVersion=2012
                        scCorrection=MITregression5
                fi

		if [ "$sample" = "totouples_parked_Run2012D_22Jan2013_v1_noSkim_v1" ] || [ "$sample" = "totouples_parked_Run2012D_22Jan2013_v1_noSkim_v1_part1" ] || [ "$sample" = "totouples_parked_Run2012D_22Jan2013_v1_noSkim_v1_part2" ] || [ "$sample" = "totouples_parked_Run2012D_22Jan2013_v1_noSkim_v1_part3" ]
                then
                        pileup_set=PU_S10
                        lumi_set=2012D_22Jan
                        isZgamma=0
                        analysisVersion=2012
                        scCorrection=MITregression5
                fi

	

                if [ "$sample" = "totouple_DYToMuMu_Summer12_S10_reg5_noSkim_v1" ] 
		then
                        pileup_set=PU_S10
			lumi_set=2012_22Jan       
                        isZgamma=1
			analysisVersion=2012
			scCorrection=MITregression5
                fi
		
		if [ "$sample" = "totouples_WJetsToLNu_Summer12_S10_noSkim_v1" ] || [ "$sample" = "totouple_TTJets_Summer12_S10_reg5_noSkim_v1" ] || [ "$sample" = "totouples_WJetsToLNu_Summer12_S10_noSkim_v1_part1" ] || [ "$sample" = "totouples_WJetsToLNu_Summer12_S10_noSkim_v1_part2" ] || [ "$sample" = "totouples_WJetsToLNu_Summer12_S10_noSkim_v1_part3" ]
                then
                        pileup_set=PU_S10
                        lumi_set=2012_22Jan
                        isZgamma=3
                        analysisVersion=2012
                        scCorrection=MITregression5
                fi


		if [ "$sample" = "totouples_Run2011A_12Oct2013_v1" ]
                then
                        pileup_set=PU_S13
                        lumi_set=2011_12OctA
                        isZgamma=0
                        analysisVersion=2011
                        scCorrection=MITregression2
                fi

                if [ "$sample" = "totouples_Run2011B_12Oct2013_v1" ]
                then
                        pileup_set=PU_S13
                        lumi_set=2011_12OctB
                        isZgamma=0
                        analysisVersion=2011
                        scCorrection=MITregression2
                fi



                if [ "$sample" = "totouple_DYToMuMu_Summer11_S13_v2" ] 
                then
                        pileup_set=PU_S13
                        lumi_set=2011_12OctB
                        isZgamma=1
                        analysisVersion=2011
                        scCorrection=MITregression2
                fi

		if [ "$sample" = "totouple_TTJets_Summer11_S13" ] || [ "$sample" = "totouple_WJetsToLNu_Summer11_S13" ]
                then
                        pileup_set=PU_S13
                        lumi_set=2011_12OctB
                        isZgamma=3
                        analysisVersion=2011
                        scCorrection=MITregression2
                fi	
		
                for inj_resolution in '0' ##'0.5' '1' '1.5' '2' '2.5' '3' '3.5' '4'
		do

				echo "sample = $sample, pileup_set = ${pileup_set}, lumi_set = ${lumi_set}, isZgamma = ${isZgamma}, inj_resolution = ${inj_resolution}"
				
				##qsub batchJob_muons.sh ${sample} muons_${sample}_June_v1_reduced ${isZgamma} ${lumi_set} ${pileup_set} ${lowMuMuCut} ${highMuMuCut} ${scCorrection} 1.0 3 0.0 ${inj_resolution} ${itoy} ${analysisVersion} -l os=sl6 
				##hadd miniTree_muons_${sample}_June_v1_reduced_partALL.root miniTree_muons_${sample}_June_v1_reduced_part[0-9]*root
                        	##echo "miniTree_${sample}_${Zgamma2}_${SetOfCorrections2}_Triangle_v1_part[0-9].root"
				mv miniTree_muons_${sample}_June_v1_reduced_part[0-9]*root stored_miniTree/
                        	mv muons_${sample}_June_v1_reduced_part[0-9]*err stored_errput/
                        	mv muons_${sample}_June_v1_reduced_part[0-9]*out stored_output/
		done

done

