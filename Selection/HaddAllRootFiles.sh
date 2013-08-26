#!/bin/bash
# script to hadd all root files into big ones
# Written by O. Bondu (August 2010)

syntax="Syntax: ${0} {cutVersion}"
if [[ -z ${1} ]]
then
  echo ${syntax}
  exit 1
fi

cutVersion=${1}

for sample in `'ls' -l Selected/${cutVersion} | grep drw | awk '{print $9}' | grep -v OLD`
do
	echo "Processing ${cutVersion} ${sample}"
#	summedRootFile=`'ls' Selected/${cutVersion}/${sample} | grep _0.root | sed -e "s/_0/_ALL/g"`
	summedRootFile=`'ls' Selected/${cutVersion}/${sample} | grep _0_ | grep .root | sed -e "s/_0/_ALL/1"`
	if [[ -e Selected/${cutVersion}/${sample}/${summedRootFile} ]]
	then
		echo "Removing Selected/${cutVersion}/${sample}/${summedRootFile}"
		rm Selected/${cutVersion}/${sample}/${summedRootFile}
	fi
	fileList=""
	for file in `'ls' Selected/${cutVersion}/${sample} | grep root | grep -v minimini`
	do
		fileList=`echo "${fileList} Selected/${cutVersion}/${sample}/${file}"`
	done
	hadd Selected/${cutVersion}/${sample}/${summedRootFile} ${fileList}
done

exit 0;

fileList=""
for data in `'ls' Selected/${cutVersion} | grep Run2010`
do
	summedRootFile="miniTree_DATA_ALL.root"
	if [[ -e Selected/${cutVersion}/${summedRootFile} ]]
	then
		echo "Removing Selected/${cutVersion}/${summedRootFile}"
    rm Selected/${cutVersion}/${summedRootFile}
	fi
  for file in `'ls' Selected/${cutVersion}/${data} | grep ALL | grep root`
  do
    fileList=`echo "${fileList} Selected/${cutVersion}/${data}/${file}"`
  done
done
hadd Selected/${cutVersion}/${summedRootFile} ${fileList}


exit 0
