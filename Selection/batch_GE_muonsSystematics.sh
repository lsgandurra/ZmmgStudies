#! /usr/local/bin/bash -l
#$ -l ct=18000
#$ -P P_cmsf
#$ -l vmem=4G
#$ -l fsize=5G
#$ -q long
#$ -l sps=1
#$ -l dcache=1
###$ -l hpss=1
#$ -N muonsSystematics
### Merge the stdout et stderr in a single file
#$ -j y
### fichiers .e et .o copied to current working directory
#$ -cwd
#$ -m a
### set array job indices 'min-max:interval'
#$ -t 1-6

syntax="${0} {sample} {data/mc}"
#if [[ -z ${6} ]]
if [[ -z ${2} ]]
then
	echo ${syntax}
	exit 1
fi

sample=${1}
data=${2}
icatmin=`echo "${SGE_TASK_ID} - 1" | bc -ql`
icatmax="${SGE_TASK_ID}"

echo "Running job with the following parameters:"
echo -e "sample= \t${sample}"
echo -e "icatmin= \t${icatmin}"
echo -e "icatmax= \t${icatmax}"
echo -e "data= ${data}"

echo "USER=${USER}"

# LOAD CORRECT ENVIRONMENT VARIABLES FROM SPS
echo "LOAD CORRECT ENVIRONMENT VARIABLES FROM SPS"
export HOMEDIR=/afs/in2p3.fr/home/o/obondu
source ${HOMEDIR}/428v4.sh
SPSDIR=`pwd`
WORKDIR=${TMPDIR}
source ${HOMEDIR}/load_ROOT.sh

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

if [[ ! -e ${TMPDIR}/gif ]]
then
  mkdir ${TMPDIR}/gif
  mkdir ${TMPDIR}/pdf
  mkdir ${TMPDIR}/C
fi

echo "USER=${USER}"

# COPY EXECUTABLE TO WORKER
echo "COPY EXECUTABLE TO WORKER"
cp ${SPSDIR}/Zmumugamma/Selection/Muons_Systematics_v04.exe ${TMPDIR}/

echo "pwd; ls -als"
pwd; ls -als
echo ""

echo "USER=${USER}"

# EXECUTE JOB
echo "EXECUTE JOB"
cd ${TMPDIR}/
itoyMin="1"
itoyMax="25"
echo -e "itoyMin= ${itoyMin}\titoyMax= ${itoyMax}"
./Muons_Systematics_v04.exe ${sample}_${icatmin}_${icatmax}_${itoyMin}_${itoyMax} ${icatmin} ${icatmax} ${itoyMin} ${itoyMax} ${data} 2> ${sample}_${icatmin}_${icatmax}_${itoyMin}_${itoyMax}.err > ${sample}_${icatmin}_${icatmax}_${itoyMin}_${itoyMax}.out

itoyMin="25"
itoyMax="50"
echo -e "itoyMin= ${itoyMin}\titoyMax= ${itoyMax}"
./Muons_Systematics_v04.exe ${sample}_${icatmin}_${icatmax}_${itoyMin}_${itoyMax} ${icatmin} ${icatmax} ${itoyMin} ${itoyMax} ${data} 2> ${sample}_${icatmin}_${icatmax}_${itoyMin}_${itoyMax}.err > ${sample}_${icatmin}_${icatmax}_${itoyMin}_${itoyMax}.out

itoyMin="50"
itoyMax="75"
echo -e "itoyMin= ${itoyMin}\titoyMax= ${itoyMax}"
./Muons_Systematics_v04.exe ${sample}_${icatmin}_${icatmax}_${itoyMin}_${itoyMax} ${icatmin} ${icatmax} ${itoyMin} ${itoyMax} ${data} 2> ${sample}_${icatmin}_${icatmax}_${itoyMin}_${itoyMax}.err > ${sample}_${icatmin}_${icatmax}_${itoyMin}_${itoyMax}.out

itoyMin="75"
itoyMax="100"
echo -e "itoyMin= ${itoyMin}\titoyMax= ${itoyMax}"
./Muons_Systematics_v04.exe ${sample}_${icatmin}_${icatmax}_${itoyMin}_${itoyMax} ${icatmin} ${icatmax} ${itoyMin} ${itoyMax} 2> ${sample}_${icatmin}_${icatmax}_${itoyMin}_${itoyMax}.err > ${sample}_${icatmin}_${icatmax}_${itoyMin}_${itoyMax}.out

echo "pwd; ls -als"
pwd; ls -als
echo ""

# GET BACK OUTPUT FILES TO SPS
echo "GET BACK OUTPUT FILES TO SPS AND REMOVE THEM FROM DISTANT DIR"
mv ${TMPDIR}/outfile_${sample}_${icatmin}_${icatmax}*.root ${SPSDIR}/Zmumugamma/Selection/
mv ${TMPDIR}/${sample}_${icatmin}_${icatmax}*.err ${SPSDIR}/Zmumugamma/Selection/
mv ${TMPDIR}/${sample}_${icatmin}_${icatmax}*.out ${SPSDIR}/Zmumugamma/Selection/
mv ${TMPDIR}/plots_${sample}_${icatmin}_${icatmax}*_EB ${SPSDIR}/Zmumugamma/Selection/
mv ${TMPDIR}/plots_${sample}_${icatmin}_${icatmax}*_EE ${SPSDIR}/Zmumugamma/Selection/
mv ${TMPDIR}/gif/*gif ${SPSDIR}/Zmumugamma/Selection/gif/
mv ${TMPDIR}/pdf/*pdf ${SPSDIR}/Zmumugamma/Selection/pdf/
mv ${TMPDIR}/C/*C ${SPSDIR}/Zmumugamma/Selection/C/
rm ${TMPDIR}/Muons_Systematics_v04.exe

exit 0
