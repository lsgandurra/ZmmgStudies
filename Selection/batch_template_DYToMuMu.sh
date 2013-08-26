#!/usr/local/bin/bash
#PBS -l platform=LINUX,u_sps_cmsf
#PBS -l T=4286000
#PBS -l scratch=4GB
#PBS -q T
#PBS -l M=2GB
#PBS -N DYToMuMu
##PBS -o nonFSR.out
##PBS -e nonFSR.err
#PBS -mb
#PBS -me

# LOAD CORRECT ENVIRONMENT VARIABLES FROM SPS
echo "LOAD CORRECT ENVIRONMENT VARIABLES FROM SPS"
export HOMEDIR=/afs/in2p3.fr/home/o/obondu
source ${HOMEDIR}/428v2.sh
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
cp ${SPSDIR}/Zmumugamma/Selection/Apply_v13.exe ${TMPBATCH}/

echo "pwd; ls -als"
pwd; ls -als
echo ""

# COPY INPUT FILES TO WORKER
#for ifile in `seq -w 1 10`
#do
#	cp /sps/cms/obondu/CMSSW_4_2_3_patch2/src/Zmumugamma/RecoSamples/DYToMuMu_0${ifile}.root ${TMPBATCH}/
#done

# EXECUTE JOB
echo "EXECUTE JOB"
cd ${TMPBATCH}/
#Apply_v13.exe DYToMuMu v13_FSR_DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_v2 1 2> v13_FSR_DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_v2.err | tee v13_FSR_DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_v2.out
#Apply_v13.exe DYToMuMu v13_nonFSR_DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_v2 2 2> v13_nonFSR_DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_v2.err | tee v13_nonFSR_DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_v2.out


#for i in `seq 292 -1 1`; do echo ${i}; file=`head -n ${i} ~/list_DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia___ | tail -n 1`; ./TEST DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia FSR_DY_powheg_${i} 1 1.0 0.0 ${file} 2> FSR_DY_powheg_${i}.err | tee FSR_DY_powheg_${i}.out; done
#for i in `seq 292 -1 1`; do echo ${i}; file=`head -n ${i} ~/list_DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia___ | tail -n 1`; ./TEST DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia nonFSR_DY_powheg_${i} 2 1.0 0.0 ${file} 2> nonFSR_DY_powheg_${i}.err | tee nonFSR_DY_powheg_${i}.out; done

#hadd miniTree_FSR_DY_powheg_ALL.root miniTree_FSR_DY_powheg_1.root miniTree_FSR_DY_powheg_2.root
#hadd miniTree_nonFSR_DY_powheg_ALL.root miniTree_nonFSR_DY_powheg_1.root miniTree_nonFSR_DY_powheg_2.root
#cp miniTree_FSR_DY_powheg_ALL.root temp.root
#for i in `seq 3 292`; do echo ${i}; rm temp.root; hadd temp.root miniTree_FSR_DY_powheg_ALL.root miniTree_FSR_DY_powheg_${i}.root; cp temp.root miniTree_FSR_DY_powheg_ALL.root; done
#for i in `seq 3 292`; do echo ${i}; rm temp.root; hadd temp.root miniTree_nonFSR_DY_powheg_ALL.root miniTree_nonFSR_DY_powheg_${i}.root; cp temp.root miniTree_nonFSR_DY_powheg_ALL.root; done


echo "pwd; ls -als"
pwd; ls -als
echo ""

# GET BACK OUTPUT FILES TO SPS
echo "GET BACK OUTPUT FILES TO SPS AND REMOVE THEM FROM DISTANT DIR"
mv ${TMPBATCH}/*v13_*_DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_v2* ${SPSDIR}/Zmumugamma/Selection/
#mv ${TMPBATCH}/miniTree*root /sps/cms/obondu/CMSSW_4_2_3_patch2/src/Zmumugamma/Selection/temp_DYToMuMu_batch/
#mv ${TMPBATCH}/*FSR_DY_powheg_*out /sps/cms/obondu/CMSSW_4_2_3_patch2/src/Zmumugamma/Selection/temp_DYToMuMu_batch/
#mv ${TMPBATCH}/*FSR_DY_powheg_*err /sps/cms/obondu/CMSSW_4_2_3_patch2/src/Zmumugamma/Selection/temp_DYToMuMu_batch/
#rm ${TMPBATCH}/TEST

#echo "cd /sps/cms/obondu/CMSSW_4_2_3_patch2/src/Zmumugamma/Selection/"
#cd /sps/cms/obondu/CMSSW_4_2_3_patch2/src/Zmumugamma/Selection/

echo "pwd; ls -als"
pwd; ls -als
echo ""

exit 0

