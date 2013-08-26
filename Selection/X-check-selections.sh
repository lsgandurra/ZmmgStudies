#!/bin/bash
# Parse the errlog files of data for a given selection to X-check what is happening of concurrent selection events
# Written by O. Bondu (September 2010)

syntax="${0} {selection} {Guy}"
if [[ -z ${1} ]]
then
	echo ${syntax}
	exit 1
else
	if [[ -z ${2} ]]
	then
		echo ${syntax}
		exit 2
	fi
fi

selection=${1}
guy=${2}
echo "*** X-check_${selection}_${guy}.txt ***"
rm X-check_${selection}_${guy}.txt
for event in `cat eventlist_${guy}.txt | awk '{print $3}'`
#for event in `cat eventlist_${guy}.txt | awk '{print $2}'`
do
	lumi=`grep ${event} eventlist_${guy}.txt | awk '{print $2}'`
	run=`grep ${event} eventlist_${guy}.txt | awk '{print $1}'`
	echo "looking for ( ${run} , ${lumi} , ${event} )"
	grep -B 10 "( ${run} , ${lumi} , ${event} )" Selected/${selection}/part*/part*err >> X-check_${selection}_${guy}.txt
#	grep -B 10 " , ${event} )" Selected/${selection}/part*/part*err >> X-check_${selection}_${guy}.txt
	echo "" >> X-check_${selection}_${guy}.txt
done


exit 0
