#! /usr/local/bin/bash -l
#$ -l ct=2000
#$ -P P_cmsf
#$ -l vmem=2G
#$ -l fsize=1G
#$ -q medium
#$ -l sps=1
#$ -N NAME
### Merge the stdout et stderr in a single file
#$ -j y
### fichiers .e et .o copied to current working directory
#$ -cwd
###$ -m be
### set array job indices 'min-max:interval'
#$ -t 1-6

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
cp ${SPSDIR}/Zmumugamma/Selection/DrawDataMC.h ${TMPDIR}/
cp ${SPSDIR}/Zmumugamma/Selection/setTDRStyle.C ${TMPDIR}/

echo "USER=${USER}"

# COPY miniTrees TO WORKER
echo "COPY miniTrees TO WORKER"
for sample in `echo "Data FSR_DYToMuMu nonFSR_DYToMuMu TTJets WJetsToLNu QCDMu"`
do
	miniTree=`cat ${SPSDIR}/Zmumugamma/Selection/plotDataMC_TDR_miniTree.C | grep -v "//" | grep "string ${sample}" | cut -d \" -f 2 | sed -e "s/%s/LUMI/g"`
	echo "Copying ${miniTree}"
	cp ${SPSDIR}/Zmumugamma/Selection/${miniTree} ${TMPDIR}/
done

for img in `echo "C png gif pdf eps"`
do
	mkdir ${TMPDIR}/${img}
done

echo "USER=${USER}"

# ADD CURRENT DIRECTORY AND LIB TO LIRARY PATH
LD_LIBRARY_PATH=`echo "${LD_LIBRARY_PATH}:${WORKDIR}/lib:${WORKDIR}"`
echo "LD_LIBRARY_PATH"
echo ${LD_LIBRARY_PATH}
echo ""

echo "USER=${USER}"

# COPY EXECUTABLE TO WORKER
echo "COPY EXECUTABLE TO WORKER"
cp ${SPSDIR}/Zmumugamma/Selection/PLOT.exe ${TMPDIR}/

echo "pwd; ls -als"
pwd; ls -als
echo ""

echo "USER=${USER}"

# EXECUTE JOB
echo "EXECUTE JOB"
cd ${TMPDIR}/
PLOT.exe ${1} ${ijob}

echo "pwd; ls -als"
pwd; ls -als
echo ""

# GET BACK OUTPUT FILES TO SPS
echo "GET BACK OUTPUT FILES TO SPS AND REMOVE THEM FROM DISTANT DIR"
for img in `echo "C png gif pdf eps"`
do
	cp -r ${TMPDIR}/${img} ${SPSDIR}/Zmumugamma/Selection/${2}
	rm -r ${TMPDIR}/${img}
done
for sample in `echo "Data FSR_DYToMuMu nonFSR_DYToMuMu TTJets WJetsToLNu QCDMu"`
do
  miniTree=`cat ${SPSDIR}/Zmumugamma/Selection/plotDataMC_TDR_miniTree.C | grep -v "//" | grep "string ${sample}" | cut -d \" -f 2 | sed -e "s/%s/LUMI/g"`
  echo "Deleting ${miniTree}"
  rm ${TMPDIR}/${miniTree}
done
rm ${TMPDIR}/PLOT.exe
rm ${TMPDIR}/DrawDataMC.h 
rm ${TMPDIR}/setTDRStyle.C



exit 0
