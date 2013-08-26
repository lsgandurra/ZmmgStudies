#!/bin/bash
# script to merge / sort everything if a job completed successfully
# Written by O. Bondu (March 2012)

syntax="${0} (job number)"
if [[ -z ${1} ]]
then
	echo "${syntax}"
	return 234
fi

ijob=${1}
njobs="100"
#itoybegin="1"; itoyend="25"
itoybegin="26"; itoyend="50"
#itoybegin="51"; itoyend="75"
#itoybegin="76"; itoyend="99"
ntoys=`echo "${itoyend} - ${itoybegin} + 1" | bc -ql`
#ntoys="24"
let ntot="${njobs}*${ntoys}"
echo "ijob= ${ijob}"
echo "njobs= ${njobs}"
echo "ntoys= ${ntoys}"
echo "ntot= ${ntot}"
echo "itoybegin= ${itoybegin}"
echo "itoyend= ${itoyend}"


# STOP if job still running on the batch
qstat | grep ${ijob}
if [[ "$?" == "0" ]]
then
	echo "job still running on the batch"
	return 22
fi

name=`'grep' "argv\[2\]=" Selection.o${ijob}.1 | awk '{print $2}' | head -n 1`
echo "Treating job: ${ijob} corresponding to sample ${name}"

nfinished=`'grep' "TOTALnbEventsAfterMuMuGammaID\[7\]=" Selection.o${ijob}.*| wc -l`

# STOP if jobs did not finish
echo "Number of finished jobs= ${nfinished}"
if [[ "${nfinished}" -lt "${ntot}" ]]
then
	echo "storing GE logs in stored_logGE/"
	mv Selection.o${ijob}.* stored_logGE/
	echo "stopping, GE failed"
	return 1
elif [[ "${nfinished}" -gt "${ntot}" ]]
then
	echo "more finished jobs than expected, please investigate"
	return 2
fi

ntrees="0"
# STOP if jobs finished but minitrees not copied
for itoy in `seq -w ${itoybegin} ${itoyend}`
do
	ntemp=`'ls' miniTree_${name}_part[0-9]*_toy${itoy}.root | wc -l`
	ntrees=`echo "${ntrees} + ${ntemp}" | bc -ql`
done

#ntrees=`'ls' miniTree_${name}_part[0-9]*.root | wc -l`
echo "ntrees= ${ntrees}"
if [[ "${ntrees}" != "${ntot}" ]]
then
	echo "erh, ${njobs} successful jobs (hence ntot= ${ntot}) but only ${ntrees} miniTrees ? Something wrong here"
	n_m_one=`echo "${njobs} -1" | bc -ql`
	for ifile in `seq 0 ${n_m_one}`
	do
		for itoy in `seq -w ${itoybegin} ${itoyend}`
		do
			if [[ ! -e "miniTree_${name}_part${ifile}_toy${itoy}.root" ]]
			then
				echo "File miniTree_${name}_part${ifile}_toy${itoy}.root is missing!!"
			fi
		done
	done
	mv Selection.o${ijob}.* stored_logGE/
	echo "stopping, GE failed"
	return 1
fi



echo "storing GE logs in stored_logGE/"
mv Selection.o${ijob}.* stored_logGE/

echo "merging trees into miniTree_${name}_toyITOY_partALL.root"
for itoy in `seq -w ${itoybegin} ${itoyend}`
do
	echo "storing outputs in stored_output/"
	mv ${name}_part*_toy${itoy}*out stored_output/
	echo "storing errputs in stored_errput/"
	mv ${name}_part*_toy${itoy}*err stored_errput/
	echo "###### MERGING TREES FOR itoy= ${itoy} #####"
	if [[ -e miniTree_${name}_toy${itoy}_partALL.root ]]
	then
		rm miniTree_${name}_toy${itoy}_partALL.root
	fi
	hadd -f miniTree_${name}_toy${itoy}_partALL.root miniTree_${name}_part[0-9]*_toy${itoy}.root
	echo "storing trees into stored_miniTree/"
	mv miniTree_${name}_part[0-9]*_toy${itoy}.root stored_miniTree/
done


return 0
