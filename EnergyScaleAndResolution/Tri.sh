#./bin/bash

echo "argument 1 = ${1}"
echo "argument 2 = ${2}"

for i in `seq 0 ${2}` 
do 
 	echo Mmumu_${i}.txt >> ${1}liste.txt 
done


