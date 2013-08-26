#! /usr/local/bin/bash -l

ijob=2
if [ ! -e "${ijob}.done" ]
then
	echo "${ijob}.done n'existe pas"
	touch SAMPLE_${ijob}.fail	
fi
