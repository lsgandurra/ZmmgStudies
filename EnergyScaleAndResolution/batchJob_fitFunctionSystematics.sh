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
cp ${SPSDIR}/Energy_scale_extraction/fitFunctionSystematics.exe ${TMPDIR}/
cp /afs/in2p3.fr/home/s/sgandurr/loadRoot.sh ${TMPDIR}/ 


echo "LOAD GOOD ROOT VERSION"
source loadRoot.sh


echo "pwd; ls -als"
pwd; ls -als
echo ""

echo "USER=${USER}"

# EXECUTE JOB
echo "EXECUTE JOB"
cd ${TMPDIR}/

echo "arguments : ${1}, ${2}, ${3}, ${4}, ${5}, ${6}." 

./fitFunctionSystematics.exe ${1} ${2} ${3} ${4} ${5} ${6} >> sortie_${2}_${3}_${4}_${5}_${6}.out 2> sortie_${2}_${3}_${4}_${5}_${6}.err

echo "pwd; ls -als"
pwd; ls -als
echo ""

# GET BACK OUTPUT FILES TO SPS
echo "GET BACK OUTPUT FILES TO SPS AND REMOVE THEM FROM DISTANT DIR"
mkdir -p ${SPSDIR}/Energy_scale_extraction/logFiles/
mv ${TMPDIR}/sortie_${2}_${3}_${4}_${5}_${6}.err ${SPSDIR}/Energy_scale_extraction/logFiles/
mv ${TMPDIR}/sortie_${2}_${3}_${4}_${5}_${6}.out ${SPSDIR}/Energy_scale_extraction/logFiles/
rm ${TMPDIR}/fitFunctionSystematics.exe
rm ${TMPDIR}/setTDRStyle.C

exit 0


