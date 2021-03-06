#! /usr/local/bin/bash -l
#$ -l ct=40000        ###time in seconds
#$ -P P_cmsf         
#$ -l vmem=4G
#$ -l fsize=30G
#$ -q long
#$ -l sps=1
###$ -l hpss=1
#$ -N EnergyScale_February2012_ 
### Merge the stdout et stderr in a single file
#$ -j y
### fichiers .e et .o copied to current working directory
#$ -cwd
###$ -m be
### set array job indices 'min-max:interval'
####$ -t 1-16


syntax="${0} {parameter}"
#if [[ -z ${6} ]]
if [[ -z ${2} ]]
then
        echo ${syntax}
        exit 1
fi

###ijob=`echo "${SGE_TASK_ID} - 1" | bc -ql`

echo "USER=${USER}"

# LOAD CORRECT ENVIRONMENT VARIABLES FROM SPS
echo "LOAD CORRECT ENVIRONMENT VARIABLES FROM SPS"
###export HOMEDIR=/afs/in2p3.fr/home/o/obondu
###source ${HOMEDIR}/428v2.sh
export HOMEDIR=/afs/in2p3.fr/home/s/sgandurr
source ${HOMEDIR}/5311p3_RECO_5_3_11_v1.sh
SPSDIR=`pwd`
WORKDIR=${TMPDIR}

echo "USER=${USER}"

# CHECK THE ENVIRONMENT VARIABLES
echo "CHECK THE ENVIRONMENT VARIABLES"
echo "ROOTSYS :"
echo ${ROOTSYS}
echo ""

echo "USER=${USER}"

cd ${TMPDIR}/
echo "pwd; ls -als"
pwd; ls -als
echo ""


# COPY HEADER FILES TO WORKER
echo "COPY HEADER FILES TO WORKER"
#mkdir ${TMPDIR}/interface
#cp ${SPSDIR}/Toto/IpnTreeProducer/interface/*h ${TMPDIR}/interface/
if [[ ! -e ${TMPDIR}/interface ]]
then
  mkdir ${TMPDIR}/interface
        cp ${SPSDIR}/Toto/IpnTreeProducer/interface/*h ${TMPDIR}/interface/
fi

echo "USER=${USER}"

# COPY IpnTree LIB FILE TO WORKER
#mkdir ${TMPDIR}/lib
#cp ${SPSDIR}/Toto/IpnTreeProducer/src/libToto.so ${TMPDIR}/lib/
echo "COPY IpnTree LIB FILE TO WORKER"
if [[ ! -e ${TMPDIR}/lib ]]
then
        mkdir ${TMPDIR}/lib
        cp ${SPSDIR}/Toto/IpnTreeProducer/src/libToto.so ${TMPDIR}/lib/
fi

echo "USER=${USER}"

# ADD CURRENT DIRECTORY AND LIB TO LIRARY PATH
LD_LIBRARY_PATH=`echo "${LD_LIBRARY_PATH}:${WORKDIR}/lib:${WORKDIR}"`
echo "LD_LIBRARY_PATH"
echo ${LD_LIBRARY_PATH}
echo ""

echo "USER=${USER}"


# COPY EXECUTABLE TO WORKER
echo "COPY EXECUTABLE TO WORKER"
cp ${SPSDIR}/ZmmgStudies/EnergyScaleAndResolution/SFits.exe ${TMPDIR}/
cp ${SPSDIR}/ZmmgStudies/EnergyScaleAndResolution/*.txt ${TMPDIR}/
cp /afs/in2p3.fr/home/s/sgandurr/loadRoot.sh ${TMPDIR}/ 
##cp /sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/cvs_developpment/Selection_NewMuID/miniTree*v6*partALL.root ${TMPDIR}/
cp ${SPSDIR}/ZmmgStudies/Selection/ClosureTest_thesis_v1/miniTree_*thesis_v1*.root ${TMPDIR}/

echo "LOAD GOOD ROOT VERSION"
source loadRoot.sh


echo "pwd; ls -als"
pwd; ls -als
echo ""

echo "USER=${USER}"

# EXECUTE JOB
echo "EXECUTE JOB"
cd ${TMPDIR}/


for eta in 'Barrel' 'Endcaps'
do
	for r9 in 'low' 'high' 'all'
	do
		for fitFunction in 'voigtian' 'cruijff'
		do
			if [ "$4" = "mmg_s" ] && [ "$fitFunction" = "cruijff" ]
                	then
                        	continue
                	fi
		
			if [ "$4" = "mmg_s_true" ] && [ "$fitFunction" = "voigtian" ]
                        then
                                continue
                        fi	
			

			echo "${eta}, ${r9} r9, ${fitFunction} :" 
			./SFits.exe ${1} ${2} ${3} ${4} ${5} ${6} ${eta} ${r9} ${fitFunction} ${7} ${8} 0 1 ${9}>> sortie_${1}_${2}_${3}_${4}_${5}_${6}_${eta}_${r9}_${fitFunction}_${7}_${8}_${9}.out 2> sortie_${1}_${2}_${3}_${4}_${5}_${6}_${eta}_${r9}_${fitFunction}_${7}_${8}_${9}.err
		done
	done
done






echo "pwd; ls -als"
pwd; ls -als
echo ""

# GET BACK OUTPUT FILES TO SPS
echo "GET BACK OUTPUT FILES TO SPS AND REMOVE THEM FROM DISTANT DIR"
mkdir -p ${SPSDIR}/ZmmgStudies/EnergyScaleAndResolution/logFiles/
mv ${TMPDIR}/sortie_${1}_${2}_${3}_${4}_${5}_${6}_*_*_*_${7}_${8}_${9}.err ${SPSDIR}/ZmmgStudies/EnergyScaleAndResolution/logFiles/
mv ${TMPDIR}/sortie_${1}_${2}_${3}_${4}_${5}_${6}_*_*_*_${7}_${8}_${9}.out ${SPSDIR}/ZmmgStudies/EnergyScaleAndResolution/logFiles/
cp -r ${TMPDIR}/${2}/ ${SPSDIR}/ZmmgStudies/EnergyScaleAndResolution/
rm -rf ${2}/
rm ${TMPDIR}/SFits.exe
rm ${TMPDIR}/setTDRStyle.C
rm ${TMPDIR}/*.txt
rm ${TMPDIR}/*.root

#"cd ${SPSDIR}/Toto/IpnTreeProducer/OlivierMiniTrees/"
#cd ${SPSDIR}/Toto/IpnTreeProducer/OlivierMiniTrees/
#echo "pwd; ls -als"
#pwd; ls -als
#echo ""

exit 0


