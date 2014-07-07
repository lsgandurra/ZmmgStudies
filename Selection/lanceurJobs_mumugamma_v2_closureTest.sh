#./bin/bash

pileup_set=PU_S10
lumi_set=2012
isZgamma=0
lowMuMuCut=35
highMuMuCut=9999
##itoy=0
analysisVersion=2011
scCorrection=MITregression2
##2012 Data
for sample in 'totouples_Run2012A_22Jan2013_v1_November2013' 'totouples_parked_Run2012B_22Jan2013_v1_November2013' 'totouples_parked_Run2012C_22Jan2013_v1_November2013' 'totouples_parked_Run2012D_22Jan2013_v1_November2013'
do

		if [ "$sample" = "totouples_Run2012A_22Jan2013_v1_November2013" ]
		then
			pileup_set=PU_S10
                        lumi_set=2012A_22Jan
                        isZgamma=0
			analysisVersion=2012
			scCorrection=MITregression5
		fi
		
		if [ "$sample" = "totouples_parked_Run2012B_22Jan2013_v1_November2013" ]
                then
                        pileup_set=PU_S10
                        lumi_set=2012B_22Jan
                        isZgamma=0
                        analysisVersion=2012
			scCorrection=MITregression5
                fi

		if [ "$sample" = "totouples_parked_Run2012C_22Jan2013_v1_November2013" ]
                then
                        pileup_set=PU_S10
                        lumi_set=2012C_22Jan
                        isZgamma=0
                        analysisVersion=2012
                        scCorrection=MITregression5
                fi

		if [ "$sample" = "totouples_parked_Run2012D_22Jan2013_v1_November2013" ]
                then
                        pileup_set=PU_S10
                        lumi_set=2012D_22Jan
                        isZgamma=0
                        analysisVersion=2012
                        scCorrection=MITregression5
                fi

	

                if [ "$sample" = "totouple_DYToMuMu_Summer12_S10_reg5_FSR" ] 
		then
                        pileup_set=PU_S10
			lumi_set=2012_22Jan       
                        isZgamma=1
			analysisVersion=2012
			scCorrection=MITregression5
			sample=totouple_DYToMuMu_Summer12_S10_reg5
                fi

		if [ "$sample" = "totouple_DYToMuMu_Summer12_S10_reg5_nonFSR" ]
                then
                        pileup_set=PU_S10
                        lumi_set=2012_22Jan
                        isZgamma=2
                        analysisVersion=2012
                        scCorrection=MITregression5
                        sample=totouple_DYToMuMu_Summer12_S10_reg5
                fi
		
		if [ "$sample" = "totouples_WJetsToLNu_Summer12_S10_November2013" ] || [ "$sample" = "totouple_TTJets_Summer12_S10_reg5" ]
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



                if [ "$sample" = "totouple_DYToMuMu_Summer11_S13_v2_FSR" ] 
                then
                        pileup_set=PU_S13
                        lumi_set=2011_12Oct
                        isZgamma=1
                        analysisVersion=2011
                        scCorrection=MITregression2
                	sample=totouple_DYToMuMu_Summer11_S13_v2
		fi
	
		if [ "$sample" = "totouple_DYToMuMu_Summer11_S13_v2_nonFSR" ]
                then
                        pileup_set=PU_S13
                        lumi_set=2011_12Oct
                        isZgamma=2
                        analysisVersion=2011
                        scCorrection=MITregression2
                        sample=totouple_DYToMuMu_Summer11_S13_v2
                fi
		

		if [ "$sample" = "totouple_TTJets_Summer11_S13" ] || [ "$sample" = "totouple_WJetsToLNu_Summer11_S13" ]
                then
                        pileup_set=PU_S13
                        lumi_set=2011_12Oct
                        isZgamma=3
                        analysisVersion=2011
                        scCorrection=MITregression2
                fi	
		
                for inj_resolution in '0' ##'0.5' '1' '1.5' '2' '2.5' '3' '3.5' '4'
		do

			for itoy in '0' ##`seq 1 100` modif title with MuSys_itoy${itoy}
			do
				for extraScale in '0.998546' '0.989881' '0.994203' '0.995068' '0.982388' '0.988869'
				do	
			
				echo "sample = $sample, pileup_set = ${pileup_set}, lumi_set = ${lumi_set}, isZgamma = ${isZgamma}, inj_resolution = ${inj_resolution}, analysisVersion = ${analysisVersion}"
		
					##qsub batchJob.sh ${sample} ${sample}_${isZgamma}_extraScale${extraScale}_thesis_v1 ${isZgamma} ${lumi_set} ${pileup_set} ${lowMuMuCut} ${highMuMuCut} ${scCorrection} ${extraScale} 3 0.0 ${inj_resolution} ${itoy} ${analysisVersion} -l os=sl6 
                                	##hadd miniTree_${sample}_${isZgamma}_extraScale${extraScale}_thesis_v1_partALL.root miniTree_${sample}_${isZgamma}_extraScale${extraScale}_thesis_v1_part[0-9]*root
                                	mv miniTree_${sample}_${isZgamma}_extraScale${extraScale}_thesis_v1_part[0-9]*root stored_miniTree/
                                	mv ${sample}_${isZgamma}_extraScale${extraScale}_thesis_v1_part[0-9]*err stored_errput/
					mv ${sample}_${isZgamma}_extraScale${extraScale}_thesis_v1_part[0-9]*out stored_output/
				done
			done

		done

done


##2012 MC
for sample in 'totouple_DYToMuMu_Summer12_S10_reg5_FSR' 'totouple_DYToMuMu_Summer12_S10_reg5_nonFSR' 'totouples_WJetsToLNu_Summer12_S10_November2013' 'totouple_TTJets_Summer12_S10_reg5'
do

		if [ "$sample" = "totouples_Run2012A_22Jan2013_v1_November2013" ]
		then
			pileup_set=PU_S10
                        lumi_set=2012A_22Jan
                        isZgamma=0
			analysisVersion=2012
			scCorrection=MITregression5
		fi
		
		if [ "$sample" = "totouples_parked_Run2012B_22Jan2013_v1_November2013" ]
                then
                        pileup_set=PU_S10
                        lumi_set=2012B_22Jan
                        isZgamma=0
                        analysisVersion=2012
			scCorrection=MITregression5
                fi

		if [ "$sample" = "totouples_parked_Run2012C_22Jan2013_v1_November2013" ]
                then
                        pileup_set=PU_S10
                        lumi_set=2012C_22Jan
                        isZgamma=0
                        analysisVersion=2012
                        scCorrection=MITregression5
                fi

		if [ "$sample" = "totouples_parked_Run2012D_22Jan2013_v1_November2013" ]
                then
                        pileup_set=PU_S10
                        lumi_set=2012D_22Jan
                        isZgamma=0
                        analysisVersion=2012
                        scCorrection=MITregression5
                fi

	

                if [ "$sample" = "totouple_DYToMuMu_Summer12_S10_reg5_FSR" ] 
		then
                        pileup_set=PU_S10
			lumi_set=2012_22Jan       
                        isZgamma=1
			analysisVersion=2012
			scCorrection=MITregression5
			sample=totouple_DYToMuMu_Summer12_S10_reg5
                fi

		if [ "$sample" = "totouple_DYToMuMu_Summer12_S10_reg5_nonFSR" ]
                then
                        pileup_set=PU_S10
                        lumi_set=2012_22Jan
                        isZgamma=2
                        analysisVersion=2012
                        scCorrection=MITregression5
                        sample=totouple_DYToMuMu_Summer12_S10_reg5
                fi
		
		if [ "$sample" = "totouples_WJetsToLNu_Summer12_S10_November2013" ] || [ "$sample" = "totouple_TTJets_Summer12_S10_reg5" ]
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



                if [ "$sample" = "totouple_DYToMuMu_Summer11_S13_v2_FSR" ] 
                then
                        pileup_set=PU_S13
                        lumi_set=2011_12Oct
                        isZgamma=1
                        analysisVersion=2011
                        scCorrection=MITregression2
                	sample=totouple_DYToMuMu_Summer11_S13_v2
		fi
	
		if [ "$sample" = "totouple_DYToMuMu_Summer11_S13_v2_nonFSR" ]
                then
                        pileup_set=PU_S13
                        lumi_set=2011_12Oct
                        isZgamma=2
                        analysisVersion=2011
                        scCorrection=MITregression2
                        sample=totouple_DYToMuMu_Summer11_S13_v2
                fi
		

		if [ "$sample" = "totouple_TTJets_Summer11_S13" ] || [ "$sample" = "totouple_WJetsToLNu_Summer11_S13" ]
                then
                        pileup_set=PU_S13
                        lumi_set=2011_12Oct
                        isZgamma=3
                        analysisVersion=2011
                        scCorrection=MITregression2
                fi	
		
                for inj_resolution in '0' ##'0.5' '1' '1.5' '2' '2.5' '3' '3.5' '4'
		do

			for itoy in '0' ##`seq 1 100` modif title with MuSys_itoy${itoy}
			do
				for extraScale in '0.9984' '1.00245' '1.00034' '0.98276' '0.99439' '0.988379'
				do	
			
				echo "sample = $sample, pileup_set = ${pileup_set}, lumi_set = ${lumi_set}, isZgamma = ${isZgamma}, inj_resolution = ${inj_resolution}, analysisVersion = ${analysisVersion}"
		
					##qsub batchJob.sh ${sample} ${sample}_${isZgamma}_extraScale${extraScale}_thesis_v1 ${isZgamma} ${lumi_set} ${pileup_set} ${lowMuMuCut} ${highMuMuCut} ${scCorrection} ${extraScale} 3 0.0 ${inj_resolution} ${itoy} ${analysisVersion} -l os=sl6 
                                	##hadd miniTree_${sample}_${isZgamma}_extraScale${extraScale}_thesis_v1_partALL.root miniTree_${sample}_${isZgamma}_extraScale${extraScale}_thesis_v1_part[0-9]*root
                                	mv miniTree_${sample}_${isZgamma}_extraScale${extraScale}_thesis_v1_part[0-9]*root stored_miniTree/
                                	mv ${sample}_${isZgamma}_extraScale${extraScale}_thesis_v1_part[0-9]*err stored_errput/
					mv ${sample}_${isZgamma}_extraScale${extraScale}_thesis_v1_part[0-9]*out stored_output/
				done
			done

		done

done

##2011 Data
for sample in 'totouples_Run2011A_12Oct2013_v1' 'totouples_Run2011B_12Oct2013_v1' 
do

		if [ "$sample" = "totouples_Run2012A_22Jan2013_v1_November2013" ]
		then
			pileup_set=PU_S10
                        lumi_set=2012A_22Jan
                        isZgamma=0
			analysisVersion=2012
			scCorrection=MITregression5
		fi
		
		if [ "$sample" = "totouples_parked_Run2012B_22Jan2013_v1_November2013" ]
                then
                        pileup_set=PU_S10
                        lumi_set=2012B_22Jan
                        isZgamma=0
                        analysisVersion=2012
			scCorrection=MITregression5
                fi

		if [ "$sample" = "totouples_parked_Run2012C_22Jan2013_v1_November2013" ]
                then
                        pileup_set=PU_S10
                        lumi_set=2012C_22Jan
                        isZgamma=0
                        analysisVersion=2012
                        scCorrection=MITregression5
                fi

		if [ "$sample" = "totouples_parked_Run2012D_22Jan2013_v1_November2013" ]
                then
                        pileup_set=PU_S10
                        lumi_set=2012D_22Jan
                        isZgamma=0
                        analysisVersion=2012
                        scCorrection=MITregression5
                fi

	

                if [ "$sample" = "totouple_DYToMuMu_Summer12_S10_reg5_FSR" ] 
		then
                        pileup_set=PU_S10
			lumi_set=2012_22Jan       
                        isZgamma=1
			analysisVersion=2012
			scCorrection=MITregression5
			sample=totouple_DYToMuMu_Summer12_S10_reg5
                fi

		if [ "$sample" = "totouple_DYToMuMu_Summer12_S10_reg5_nonFSR" ]
                then
                        pileup_set=PU_S10
                        lumi_set=2012_22Jan
                        isZgamma=2
                        analysisVersion=2012
                        scCorrection=MITregression5
                        sample=totouple_DYToMuMu_Summer12_S10_reg5
                fi
		
		if [ "$sample" = "totouples_WJetsToLNu_Summer12_S10_November2013" ] || [ "$sample" = "totouple_TTJets_Summer12_S10_reg5" ]
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



                if [ "$sample" = "totouple_DYToMuMu_Summer11_S13_v2_FSR" ] 
                then
                        pileup_set=PU_S13
                        lumi_set=2011_12Oct
                        isZgamma=1
                        analysisVersion=2011
                        scCorrection=MITregression2
                	sample=totouple_DYToMuMu_Summer11_S13_v2
		fi
	
		if [ "$sample" = "totouple_DYToMuMu_Summer11_S13_v2_nonFSR" ]
                then
                        pileup_set=PU_S13
                        lumi_set=2011_12Oct
                        isZgamma=2
                        analysisVersion=2011
                        scCorrection=MITregression2
                        sample=totouple_DYToMuMu_Summer11_S13_v2
                fi
		

		if [ "$sample" = "totouple_TTJets_Summer11_S13" ] || [ "$sample" = "totouple_WJetsToLNu_Summer11_S13" ]
                then
                        pileup_set=PU_S13
                        lumi_set=2011_12Oct
                        isZgamma=3
                        analysisVersion=2011
                        scCorrection=MITregression2
                fi	
		
                for inj_resolution in '0' ##'0.5' '1' '1.5' '2' '2.5' '3' '3.5' '4'
		do

			for itoy in '0' ##`seq 1 100` modif title with MuSys_itoy${itoy}
			do
				for extraScale in '0.977904' '0.982256' '0.980437' '0.993921' '1.00056' '0.996481'
				do	
			
				echo "sample = $sample, pileup_set = ${pileup_set}, lumi_set = ${lumi_set}, isZgamma = ${isZgamma}, inj_resolution = ${inj_resolution}, analysisVersion = ${analysisVersion}"
		
					##qsub batchJob.sh ${sample} ${sample}_${isZgamma}_extraScale${extraScale}_thesis_v1 ${isZgamma} ${lumi_set} ${pileup_set} ${lowMuMuCut} ${highMuMuCut} ${scCorrection} ${extraScale} 3 0.0 ${inj_resolution} ${itoy} ${analysisVersion} -l os=sl6 
                                	##hadd miniTree_${sample}_${isZgamma}_extraScale${extraScale}_thesis_v1_partALL.root miniTree_${sample}_${isZgamma}_extraScale${extraScale}_thesis_v1_part[0-9]*root
                                	mv miniTree_${sample}_${isZgamma}_extraScale${extraScale}_thesis_v1_part[0-9]*root stored_miniTree/
                                	mv ${sample}_${isZgamma}_extraScale${extraScale}_thesis_v1_part[0-9]*err stored_errput/
					mv ${sample}_${isZgamma}_extraScale${extraScale}_thesis_v1_part[0-9]*out stored_output/
				done
			done

		done

done

##2011 MC
for sample in 'totouple_DYToMuMu_Summer11_S13_v2_FSR' 'totouple_DYToMuMu_Summer11_S13_v2_nonFSR' 'totouple_WJetsToLNu_Summer11_S13' 'totouple_TTJets_Summer11_S13'
do

		if [ "$sample" = "totouples_Run2012A_22Jan2013_v1_November2013" ]
		then
			pileup_set=PU_S10
                        lumi_set=2012A_22Jan
                        isZgamma=0
			analysisVersion=2012
			scCorrection=MITregression5
		fi
		
		if [ "$sample" = "totouples_parked_Run2012B_22Jan2013_v1_November2013" ]
                then
                        pileup_set=PU_S10
                        lumi_set=2012B_22Jan
                        isZgamma=0
                        analysisVersion=2012
			scCorrection=MITregression5
                fi

		if [ "$sample" = "totouples_parked_Run2012C_22Jan2013_v1_November2013" ]
                then
                        pileup_set=PU_S10
                        lumi_set=2012C_22Jan
                        isZgamma=0
                        analysisVersion=2012
                        scCorrection=MITregression5
                fi

		if [ "$sample" = "totouples_parked_Run2012D_22Jan2013_v1_November2013" ]
                then
                        pileup_set=PU_S10
                        lumi_set=2012D_22Jan
                        isZgamma=0
                        analysisVersion=2012
                        scCorrection=MITregression5
                fi

	

                if [ "$sample" = "totouple_DYToMuMu_Summer12_S10_reg5_FSR" ] 
		then
                        pileup_set=PU_S10
			lumi_set=2012_22Jan       
                        isZgamma=1
			analysisVersion=2012
			scCorrection=MITregression5
			sample=totouple_DYToMuMu_Summer12_S10_reg5
                fi

		if [ "$sample" = "totouple_DYToMuMu_Summer12_S10_reg5_nonFSR" ]
                then
                        pileup_set=PU_S10
                        lumi_set=2012_22Jan
                        isZgamma=2
                        analysisVersion=2012
                        scCorrection=MITregression5
                        sample=totouple_DYToMuMu_Summer12_S10_reg5
                fi
		
		if [ "$sample" = "totouples_WJetsToLNu_Summer12_S10_November2013" ] || [ "$sample" = "totouple_TTJets_Summer12_S10_reg5" ]
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



                if [ "$sample" = "totouple_DYToMuMu_Summer11_S13_v2_FSR" ] 
                then
                        pileup_set=PU_S13
                        lumi_set=2011_12Oct
                        isZgamma=1
                        analysisVersion=2011
                        scCorrection=MITregression2
                	sample=totouple_DYToMuMu_Summer11_S13_v2
		fi
	
		if [ "$sample" = "totouple_DYToMuMu_Summer11_S13_v2_nonFSR" ]
                then
                        pileup_set=PU_S13
                        lumi_set=2011_12Oct
                        isZgamma=2
                        analysisVersion=2011
                        scCorrection=MITregression2
                        sample=totouple_DYToMuMu_Summer11_S13_v2
                fi
		

		if [ "$sample" = "totouple_TTJets_Summer11_S13" ] || [ "$sample" = "totouple_WJetsToLNu_Summer11_S13" ]
                then
                        pileup_set=PU_S13
                        lumi_set=2011_12Oct
                        isZgamma=3
                        analysisVersion=2011
                        scCorrection=MITregression2
                fi	
		
                for inj_resolution in '0' ##'0.5' '1' '1.5' '2' '2.5' '3' '3.5' '4'
		do

			for itoy in '0' ##`seq 1 100` modif title with MuSys_itoy${itoy}
			do
				for extraScale in '0.9849' '0.996171' '0.991103' '0.999223' '1.00251' '1.00114' 
				do	
			
				echo "sample = $sample, pileup_set = ${pileup_set}, lumi_set = ${lumi_set}, isZgamma = ${isZgamma}, inj_resolution = ${inj_resolution}, analysisVersion = ${analysisVersion}"
		
					##qsub batchJob.sh ${sample} ${sample}_${isZgamma}_extraScale${extraScale}_thesis_v1 ${isZgamma} ${lumi_set} ${pileup_set} ${lowMuMuCut} ${highMuMuCut} ${scCorrection} ${extraScale} 3 0.0 ${inj_resolution} ${itoy} ${analysisVersion} -l os=sl6 
                                	##hadd miniTree_${sample}_${isZgamma}_extraScale${extraScale}_thesis_v1_partALL.root miniTree_${sample}_${isZgamma}_extraScale${extraScale}_thesis_v1_part[0-9]*root
                                	mv miniTree_${sample}_${isZgamma}_extraScale${extraScale}_thesis_v1_part[0-9]*root stored_miniTree/
                                	mv ${sample}_${isZgamma}_extraScale${extraScale}_thesis_v1_part[0-9]*err stored_errput/
					mv ${sample}_${isZgamma}_extraScale${extraScale}_thesis_v1_part[0-9]*out stored_output/
				done
			done

		done

done


