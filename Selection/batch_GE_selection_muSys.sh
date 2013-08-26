#! /usr/local/bin/bash -l
#$ -l ct=30000
#$ -P P_cmsf
#$ -l vmem=1G
#$ -l fsize=5G
#$ -q long
#$ -l sps=1
#$ -l dcache=1
###$ -l hpss=1
#$ -N Selection
### Merge the stdout et stderr in a single file
#$ -j y
### fichiers .e et .o copied to current working directory
#$ -cwd
#$ -m a
### set array job indices 'min-max:interval'
#$ -t 1-100

syntax="${0} {parameter}"
#if [[ -z ${6} ]]
if [[ -z ${2} ]]
then
	echo ${syntax}
	exit 1
fi

sample=${1}
output=${2}
ntotjob="100"
ijob=`echo "${SGE_TASK_ID} - 1" | bc -ql`
isZgammaMC=${3:-""}
lumi_set=${4:-""}
pu_set=${5:-""}
low_m_mumu_cut=${6:-""}
high_m_mumu_cut=${7:-""}
photon_energy_correction_scheme=${8:-""}
extra_photon_scale=${9:-""}
applyMuonScaleCorrection=${10:-""}
muon_correction_sys=${11:-""}
extra_resolution=${12:-""}

echo "Running job with the following parameters:"
echo -e "sample= \t${sample}"
echo -e "output= \t${output}"
echo -e "ntotjob= \t${ntotjob}"
echo -e "ijob= \t${ijob}"
echo -e "isZgammaMC= \t${isZgammaMC}"
echo -e "lumi_set= \t${lumi_set}"
echo -e "pu_set= \t${pu_set}"
echo -e "low_m_mumu_cut= \t${low_m_mumu_cut}"
echo -e "high_m_mumu_cut= \t${high_m_mumu_cut}"
echo -e "photon_energy_correction_scheme= \t${photon_energy_correction_scheme}"
echo -e "extra_photon_scale= \t${extra_photon_scale}"
echo -e "applyMuonScaleCorrection= \t${applyMuonScaleCorrection}"
echo -e "muon_correction_sys= \t${muon_correction_sys}"
echo -e "extra_resolution= \t${extra_resolution}"


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
cp ${SPSDIR}/Zmumugamma/Selection/Mmumu.txt ${TMPDIR}/
cp ${SPSDIR}/Zmumugamma/Selection/rochcor_v2.h ${TMPDIR}/
cp ${SPSDIR}/Zmumugamma/Selection/Apply_v31.exe ${TMPDIR}/
cp ${SPSDIR}/Zmumugamma/Selection/listFiles_${1} ${TMPDIR}/

echo "pwd; ls -als"
pwd; ls -als
echo ""

echo "USER=${USER}"

# EXECUTE JOB
echo "EXECUTE JOB"
cd ${TMPDIR}/
for itoy in `seq -w 76 99`
do
	echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
	echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% itoy = ${itoy} %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
	echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
	echo "ijob= ${ijob}"
	./Apply_v31.exe ${sample} ${output} ${ntotjob} ${ijob} ${isZgammaMC} ${lumi_set} ${pu_set} ${low_m_mumu_cut} ${high_m_mumu_cut} ${photon_energy_correction_scheme} ${extra_photon_scale} ${applyMuonScaleCorrection} ${muon_correction_sys} 2> ${output}_part${ijob}_toy${itoy}.err | tee ${output}_part${ijob}_toy${itoy}.out
	
#	./Apply_v28.exe ${1} ${2} 100 ${ijob} ${3} ${4} ${5} ${6} ${7} ${8} ${9} ${10} ${11} 2> ${2}_part${ijob}_toy${itoy}.err | tee ${2}_part${ijob}_toy${itoy}.out
	mv miniTree_${2}_part${ijob}.root miniTree_${2}_part${ijob}_toy${itoy}.root
	mv miniFriend_${2}_part${ijob}.root miniFriend_${2}_part${ijob}_toy${itoy}.root
	echo ""
	echo ""
done

echo "pwd; ls -als"
pwd; ls -als
echo ""

# GET BACK OUTPUT FILES TO SPS
echo "GET BACK OUTPUT FILES TO SPS AND REMOVE THEM FROM DISTANT DIR"
mv ${TMPDIR}/miniTree_${2}_part${ijob}*root ${SPSDIR}/Zmumugamma/Selection/
mv ${TMPDIR}/miniFriend_${2}_part${ijob}*root ${SPSDIR}/Zmumugamma/Selection/
mv ${TMPDIR}/${2}_part${ijob}*out ${SPSDIR}/Zmumugamma/Selection/
mv ${TMPDIR}/${2}_part${ijob}*err ${SPSDIR}/Zmumugamma/Selection/
rm ${TMPDIR}/Apply_v31.exe

#"cd ${SPSDIR}/Zmumugamma/Selection/"
#cd ${SPSDIR}/Zmumugamma/Selection/
#echo "pwd; ls -als"
#pwd; ls -als
#echo ""



exit 0
