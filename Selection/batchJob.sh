#! /usr/local/bin/bash -l
#$ -l ct=40000
#$ -P P_cmsf
#$ -l vmem=4G
#$ -l fsize=29G
#$ -q long
#$ -l sps=1
#$ -l dcache=1
###$ -l hpss=1
#$ -N Selection_NewMuID_
### Merge the stdout et stderr in a single file
#$ -j y
### fichiers .e et .o copied to current working directory
#$ -cwd
###$ -m be
### set array job indices 'min-max:interval'
#$ -t 1-20

syntax="${0} {parameter}"
#if [[ -z ${6} ]]
if [[ -z ${2} ]]
then
	echo ${syntax}
	exit 1
fi

ijob=`echo "${SGE_TASK_ID} - 1" | bc -ql`

echo "USER=${USER}"

# LOAD CORRECT ENVIRONMENT VARIABLES FROM SPS
echo "LOAD CORRECT ENVIRONMENT VARIABLES FROM SPS"
###export HOMEDIR=/afs/in2p3.fr/home/o/obondu
###source ${HOMEDIR}/428v2.sh
export HOMEDIR=/afs/in2p3.fr/home/s/sgandurr
source ${HOMEDIR}/537_RECO_5_3_3_v4.sh
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
#cp ${SPSDIR}/UserCode/IpnTreeProducer/interface/*h ${TMPDIR}/interface/
if [[ ! -e ${TMPDIR}/interface ]]
then
  mkdir ${TMPDIR}/interface
	cp ${SPSDIR}/UserCode/IpnTreeProducer/interface/*h ${TMPDIR}/interface/
fi

echo "USER=${USER}"

# COPY IpnTree LIB FILE TO WORKER
#mkdir ${TMPDIR}/lib
#cp ${SPSDIR}/UserCode/IpnTreeProducer/src/libToto.so ${TMPDIR}/lib/
echo "COPY IpnTree LIB FILE TO WORKER"
if [[ ! -e ${TMPDIR}/lib ]]
then
	mkdir ${TMPDIR}/lib
	cp ${SPSDIR}/UserCode/IpnTreeProducer/src/libToto.so ${TMPDIR}/lib/
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
###cp ${SPSDIR}/Selection_July2013/Selection_miniTree.exe ${TMPDIR}/
cp ${SPSDIR}/Selection_July2013/Selection_miniTree.exe ${TMPDIR}/
cp ${SPSDIR}/Selection_July2013/*.C ${TMPDIR}/
cp ${SPSDIR}/Selection_July2013/*.h ${TMPDIR}/
cp ${SPSDIR}/Selection_July2013/*.txt ${TMPDIR}/
cp ${SPSDIR}/Selection_July2013/*.dat ${TMPDIR}/
##cp -r /sps/cms/obondu/CMSSW_4_2_8__RECO_4_2_8_v2/src/Zmumugamma/TotoSamples/DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_PU_S6_v2/ ${TMPDIR}/
##cp /sps/cms/sgandurr/CMSSW_5_3_6/src/UserCode/IpnTreeProducer/ListZmumugamma/listFiles_* ${TMPDIR}/
cp /sps/cms/sgandurr/CMSSW_5_3_7_RECO_5_3_3_v4/src/Selection_July2013/listFiles_* ${TMPDIR}/

echo "pwd; ls -als"
pwd; ls -als
echo ""

echo "USER=${USER}"

# EXECUTE JOB
echo "EXECUTE JOB"
cd ${TMPDIR}/
###./Selection_miniTree.exe ${1} ${2} ${3} ${4} ${5} ${6} 2> ${2}.err | tee ${2}.out
##./Selection_miniTree.exe ${1} ${2} 20 ${ijob} ${3} 2> ${2}_part${ijob}.err | tee ${2}_part${ijob}.out
##./Selection_miniTree.exe ${1} ${2} 100 ${ijob} ${3} 2011 PU_S6 ${4} ${5} ${6} ${7} 2> ${2}_part${ijob}.err | tee ${2}_part${ijob}.out
./Selection_miniTree.exe ${1} ${2} 20 ${ijob} ${3} ${4} ${5} ${6} ${7} ${8} ${9} ${10} ${11} ${12} ${13} 2> ${2}_part${ijob}.err | tee ${2}_part${ijob}.out


#LOOK IF EVERYTHING IS OK
if [ ! -e "miniTree_${ijob}.done" ]
then
        echo "${ijob}.done doesn't exist"
        touch ${2}_${ijob}.fail   
	mv ${2}_${ijob}.fail ${SPSDIR}/Selection_July2013/ 
fi

echo "pwd; ls -als"
pwd; ls -als
echo ""

# GET BACK OUTPUT FILES TO SPS
echo "GET BACK OUTPUT FILES TO SPS AND REMOVE THEM FROM DISTANT DIR"
mv ${TMPDIR}/miniTree_*root ${SPSDIR}/Selection_July2013/
mv ${TMPDIR}/${2}_part${ijob}.out ${SPSDIR}/Selection_July2013/
mv ${TMPDIR}/${2}_part${ijob}.err ${SPSDIR}/Selection_July2013/
rm ${TMPDIR}/Selection_miniTree.exe
rm ${TMPDIR}/*.C
rm ${TMPDIR}/*.h
rm ${TMPDIR}/*.txt
##rm -rf ${TMPDIR}/DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_PU_S6_v2/
rm ${TMPDIR}/listFiles_* 

#"cd ${SPSDIR}/Selection_July2013/"
#cd ${SPSDIR}/Selection_July2013/
#echo "pwd; ls -als"
#pwd; ls -als
#echo ""

##./Selection_miniTree.exe DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Fall11-PU_S6_START42_V14B-v1 DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Fall11-PU_S6_START42_V14B-v1_Triangle_vtest 9999 -1 1 2011 PU_S6 0 9999 1.0 1.0 3

##./Selection_miniTree.exe DYToMuMu_Summer12_NewMuonID testRegression 9999 -1 1 2012 PU_S10 0 9999 MITRegression 1.0 3

exit 0


