#./bin/bash
rm $1.txt
root -l -b -q "EventsList.C(\"$1\")" >> temp.txt
cat  temp.txt | awk '{print $4":"$6":"$8}' | tail -n +6 | head -n -3 >> $1.txt
rm temp.txt
