#./bin/bash

rm -rf Results_v6_MuSys_RecoEnergy/InjectedResolution_0Percent/MC/MuonSystematics/
rm -rf Results_v6_MuSys_RecoEnergy/InjectedResolution_0Percent/data/MuonSystematics/

for dataType in 'MC' 'data'
do
	for eta in 'Barrel' 'Endcaps' 
	do
		for r9 in 'low' 'high' 'all'
        	do
			./muonsSystematics.exe Results_v6_RecoEnergy_MuSys ${dataType} mmg_s ${eta} ${r9} voigtian
		done
	done
done
