#./bin/bash

rm -rf ${1}/InjectedResolution_0Percent/MC/MuonSystematics/
rm -rf ${1}/InjectedResolution_0Percent/data/MuonSystematics/

for dataType in 'MC' 'data'
do
	for eta in 'Barrel' 'Endcaps' 
	do
		for r9 in 'low' 'high' 'all'
        	do
			./muonsSystematics.exe ${1} ${dataType} mmg_s ${eta} ${r9} voigtian
		done
	done
done
