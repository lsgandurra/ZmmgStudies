#! /usr/local/bin/bash -l
#$ -l ct=40000
#$ -P P_cmsf
#$ -l vmem=4G
#$ -l fsize=29G
#$ -q long
#$ -l sps=1
#$ -l dcache=1
###$ -l hpss=1
#$ -N Selection_November2013_
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
##ijob=85

echo "USER=${USER}"

# LOAD CORRECT ENVIRONMENT VARIABLES FROM SPS
echo "LOAD CORRECT ENVIRONMENT VARIABLES FROM SPS"
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
if [[ ! -e ${TMPDIR}/interface ]]
then
  mkdir ${TMPDIR}/interface
	cp ${SPSDIR}/Toto/IpnTreeProducer/interface/*h ${TMPDIR}/interface/
fi

echo "USER=${USER}"

# COPY IpnTree LIB FILE TO WORKER
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
cp ${SPSDIR}/ZmmgStudies/Selection/Selection_miniTree.exe ${TMPDIR}/
cp ${SPSDIR}/ZmmgStudies/Selection/*.C ${TMPDIR}/
cp ${SPSDIR}/ZmmgStudies/Selection/*.h ${TMPDIR}/
cp ${SPSDIR}/ZmmgStudies/Selection/*.txt ${TMPDIR}/
cp ${SPSDIR}/ZmmgStudies/Selection/*.dat ${TMPDIR}/
cp ${SPSDIR}/ZmmgStudies/Selection/listFiles_* ${TMPDIR}/

echo "pwd; ls -als"
pwd; ls -als
echo ""

echo "USER=${USER}"

# EXECUTE JOB
echo "EXECUTE JOB"
cd ${TMPDIR}/
./Selection_miniTree.exe ${1} ${2} 20 ${ijob} ${3} ${4} ${5} ${6} ${7} ${8} ${9} ${10} ${11} ${12} ${13} ${14} 2> ${2}_part${ijob}.err | tee ${2}_part${ijob}.out


#LOOK IF EVERYTHING IS OK
if [ ! -e "miniTree_${ijob}.done" ]
then
        echo "${ijob}.done doesn't exist"
        touch ${2}_${ijob}.fail   
	mv ${2}_${ijob}.fail ${SPSDIR}/ZmmgStudies/Selection/ 
fi

echo "pwd; ls -als"
pwd; ls -als
echo ""

# GET BACK OUTPUT FILES TO SPS
echo "GET BACK OUTPUT FILES TO SPS AND REMOVE THEM FROM DISTANT DIR"
mv ${TMPDIR}/miniTree_*root ${SPSDIR}/ZmmgStudies/Selection/
##mv ${TMPDIR}/miniFriend2_*root ${SPSDIR}/ZmmgStudies/Selection/
##mv ${TMPDIR}/miniFriend3_*root ${SPSDIR}/ZmmgStudies/Selection/
##mv ${TMPDIR}/miniFriend4_*root ${SPSDIR}/ZmmgStudies/Selection/
mv ${TMPDIR}/${2}_part${ijob}.out ${SPSDIR}/ZmmgStudies/Selection/
mv ${TMPDIR}/${2}_part${ijob}.err ${SPSDIR}/ZmmgStudies/Selection/
rm ${TMPDIR}/miniFriend*root
rm ${TMPDIR}/Selection_miniTree.exe
rm ${TMPDIR}/*.C
rm ${TMPDIR}/*.h
rm ${TMPDIR}/*.txt
rm ${TMPDIR}/listFiles_* 

#"cd ${SPSDIR}/ZmmgStudies/Selection/"
#cd ${SPSDIR}/ZmmgStudies/Selection/
#echo "pwd; ls -als"
#pwd; ls -als
#echo ""

##./Selection_miniTree.exe DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Fall11-PU_S6_START42_V14B-v1 DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Fall11-PU_S6_START42_V14B-v1_Triangle_vtest 9999 -1 1 2011 PU_S6 0 9999 1.0 1.0 3

##./Selection_miniTree.exe DYToMuMu_Summer12_NewMuonID testRegression 9999 -1 1 2012 PU_S10 0 9999 MITRegression 1.0 3

exit 0


