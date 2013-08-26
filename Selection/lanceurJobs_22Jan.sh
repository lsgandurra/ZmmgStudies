#./bin/bash

##CHECK LUMI DATA + PILEUP MC !!!!!!!!!!!!!!!!!!!!!
pileup_set=PU_S10
lumi_set=2012
isZgamma=0
lowMuMuCut=35
highMuMuCut=9999
itoy=0
scCorrection=MITregression ##FIXME MITregression
##for sample in 'parked_Run2012A_22Jan2013_v1_3' 'parked_Run2012B_22Jan2013_v1' 'parked_Run2012C_22Jan2013_v1' 'Run2012A_22Jan2013_v1' ##'parked_Run2012D_22Jan2013_v1' 'DYToMuMu_Summer12_2013_v2_FSR' 'DYToMuMu_Summer12_2013_v2_nonFSR'  
for sample in 'parked_Run2012D_22Jan2013_v1' 'DYToMuMu_Summer12_2013_v2_FSR' 'DYToMuMu_Summer12_2013_v2_nonFSR'  
do

		if [ "$sample" = "parked_Run2012A_22Jan2013_v1_3" ] || [ "$sample" = "parked_Run2012B_22Jan2013_v1" ] || [ "$sample" = "parked_Run2012C_22Jan2013_v1" ] || [ "$sample" = "Run2012A_22Jan2013_v1" ] 
		then
			pileup_set=PU_S10
			lumi_set=2012ABC
			isZgamma=0	
		fi

		if [ "$sample" = "parked_Run2012D_22Jan2013_v1" ]
		then
			pileup_set=PU_S10
                        lumi_set=2012D
                        isZgamma=0
		fi
		
                if [ "$sample" = "DYToMuMu_Summer12_2013_v2_FSR" ]
		then
                        pileup_set=PU_S10 ##FIXME
			lumi_set=2012       
                        isZgamma=1
			sample=DYToMuMu_Summer12_2013_v2
                fi

		if [ "$sample" = "DYToMuMu_Summer12_2013_v2_nonFSR" ]
                then
                        pileup_set=PU_S10  ##FIXME
			lumi_set=2012         
                        isZgamma=2
			sample=DYToMuMu_Summer12_2013_v2
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

				echo "sample = $sample, pileup_set = ${pileup_set}, lumi_set = ${lumi_set}, isZgamma = ${isZgamma}, inj_resolution = ${inj_resolution}"
				

				##qsub batchJob.sh ${sample} ${sample}_NewSelection_${isZgamma}_injRe${inj_resolution}_v1_22Jan2013 ${isZgamma} ${lumi_set} ${pileup_set} ${lowMuMuCut} ${highMuMuCut} ${scCorrection} 1.0 3 0.0 ${inj_resolution} ${itoy} -l os=sl5 
				##hadd miniTree_${sample}_NewSelection_${isZgamma}_injRe${inj_resolution}_v1_22Jan2013_partALL.root miniTree_${sample}_NewSelection_${isZgamma}_injRe${inj_resolution}_v1_22Jan2013_part[0-9]*root
                        	##echo "miniTree_${sample}_${Zgamma2}_${SetOfCorrections2}_Triangle_v1_part[0-9].root"
				mv miniTree_${sample}_NewSelection_${isZgamma}_injRe${inj_resolution}_v1_22Jan2013_part[0-9]*root stored_miniTree/
                        	mv ${sample}_NewSelection_${isZgamma}_injRe${inj_resolution}_v1_22Jan2013_part[0-9]*err stored_errput/
                        	mv ${sample}_NewSelection_${isZgamma}_injRe${inj_resolution}_v1_22Jan2013_part[0-9]*out stored_output/
		done

done

