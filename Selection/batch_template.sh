#!/usr/local/bin/bash
#PBS -l platform=LINUX,u_sps_cmsf     # Plateforme d'execution
#PBS -l T=TIME              # Nombre d'unite normalisee (consommation cpu)
#PBS -l scratch=SCRATCHSPACE
#PBS -q QUEUE
#PBS -l M=MEMORY
########PBS -N NAME               # Job Name
#PBS -o OUTLOG.out
#PBS -e ERRLOG.err

# LOAD CORRECT ENVIRONMENT VARIABLES FROM SPS
echo "LOAD CORRECT ENVIRONMENT VARIABLES FROM SPS"
export HOMEDIR=/afs/in2p3.fr/home/o/obondu
source ${HOMEDIR}/397v2.sh
SPSDIR=`pwd`
WORKDIR=${TMPBATCH}

# CHECK THE ENVIRONMENT VARIABLES
echo "CHECK THE ENVIRONMENT VARIABLES"
echo "ROOTSYS :"
echo ${ROOTSYS}
echo ""

cd ${TMPBATCH}/
echo "pwd; ls -als"
pwd; ls -als
echo ""

# COPY HEADER FILES TO WORKER
echo "COPY HEADER FILES TO WORKER"
if [[ ! -e ${TMPBATCH}/interface ]]
then
  mkdir ${TMPBATCH}/interface
  cp ${SPSDIR}/UserCode/IpnTreeProducer/interface/*h ${TMPBATCH}/interface/
fi

# COPY IpnTree LIB FILE TO WORKER
echo "COPY IpnTree LIB FILE TO WORKER"
if [[ ! -e ${TMPBATCH}/lib ]]
then
  mkdir ${TMPBATCH}/lib
  cp ${SPSDIR}/UserCode/IpnTreeProducer/src/libToto.so ${TMPBATCH}/lib/
	cp ${SPSDIR}/UserCode/IpnTreeProducer/src/libToto.so ${TMPBATCH}/
fi

# ADD CURRENT DIRECTORY AND LIB TO LIRARY PATH
LD_LIBRARY_PATH=`echo "${LD_LIBRARY_PATH}:${WORKDIR}/lib:${WORKDIR}"`
echo "LD_LIBRARY_PATH"
echo ${LD_LIBRARY_PATH}
echo ""

# COPY EXECUTABLE TO WORKER
echo "COPY EXECUTABLE TO WORKER"
#cp ${SPSDIR}/Zmumugamma/Selection/LOCATION/MACRO ${TMPBATCH}/
cp EXEDIR/MACRO ${TMPBATCH}/

echo "pwd; ls -als"
pwd; ls -als
echo ""

# EXECUTE JOB
echo "EXECUTE JOB"
cd ${TMPBATCH}/
./MACRO

echo "pwd; ls -als"
pwd; ls -als
echo ""

# GET BACK OUTPUT FILES TO SPS
echo "GET BACK OUTPUT FILES TO SPS AND REMOVE THEM FROM DISTANT DIR"
mv ${TMPBATCH}/miniTree*root EXEDIR/
mv ${TMPBATCH}/SAMPLEPART*out EXEDIR/
mv ${TMPBATCH}/SAMPLEPART*err EXEDIR/
rm ${TMPBATCH}/MACRO

echo "cd EXEDIR/"
cd EXEDIR/

echo "pwd; ls -als"
pwd; ls -als
echo ""

exit 0

