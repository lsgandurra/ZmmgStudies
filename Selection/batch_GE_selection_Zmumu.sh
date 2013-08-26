#! /usr/local/bin/bash -l
#$ -l ct=2000
#$ -P P_cmsf
#$ -l vmem=3G
#$ -l fsize=1G
#$ -q medium
#$ -l sps=1
#$ -l dcache=1
###$ -l hpss=1
#$ -N MuonsSelection
### Merge the stdout et stderr in a single file
#$ -j y
### fichiers .e et .o copied to current working directory
#$ -cwd
###$ -m be
### set array job indices 'min-max:interval'
#$ -t 1-30

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
export HOMEDIR=/afs/in2p3.fr/home/o/obondu
source ${HOMEDIR}/428v2.sh
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
cp ${SPSDIR}/Zmumugamma/Selection/rochcor_v2.h ${TMPDIR}/
cp ${SPSDIR}/Zmumugamma/Selection/Muons_v10.exe ${TMPDIR}/
cp ${SPSDIR}/Zmumugamma/Selection/listFiles_${1} ${TMPDIR}/

echo "pwd; ls -als"
pwd; ls -als
echo ""

echo "USER=${USER}"

# EXECUTE JOB
echo "EXECUTE JOB"
cd ${TMPDIR}/
./Muons_v10.exe ${1} ${2} 30 ${ijob} ${3} ${4} ${5} ${6} ${7} ${8} 2> ${2}_part${ijob}.err | tee ${2}_part${ijob}.out
#./Muons_v10.exe ${1} ${2} ${ijob} 3 Nov03 PU_S6 40 80 START42_V11 2> ${2}_part${ijob}.err | tee ${2}_part${ijob}.out

echo "pwd; ls -als"
pwd; ls -als
echo ""

# GET BACK OUTPUT FILES TO SPS
echo "GET BACK OUTPUT FILES TO SPS AND REMOVE THEM FROM DISTANT DIR"
mv ${TMPDIR}/miniTreeMuons_${2}_part${ijob}*root ${SPSDIR}/Zmumugamma/Selection/
mv ${TMPDIR}/${2}_part${ijob}*out ${SPSDIR}/Zmumugamma/Selection/
mv ${TMPDIR}/${2}_part${ijob}*err ${SPSDIR}/Zmumugamma/Selection/
rm ${TMPDIR}/Muons_v10.exe

#"cd ${SPSDIR}/Zmumugamma/Selection/"
#cd ${SPSDIR}/Zmumugamma/Selection/
#echo "pwd; ls -als"
#pwd; ls -als
#echo ""



exit 0
