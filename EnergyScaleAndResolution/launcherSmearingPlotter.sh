#./bin/bash

for eta in 'Barrel_1' 'Barrel_2' 'Endcaps_1' 'Endcaps_2' 
do
	for r9 in 'low' 'high' 'all'
        do
		./smearingPlotter.exe smearingFactors ${eta} ${r9} 25 
	done
done
