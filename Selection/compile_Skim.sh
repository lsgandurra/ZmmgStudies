#!/bin/bash
syntax="Syntax ${0} {executableName}"
if [[ -z ${1} ]]
then
  echo ${syntax}
#  exit 1
fi
exeFile=${1}

g++ SkimSelection_miniTree.C -L/afs/in2p3.fr/grid/toolkit/cms2/slc5_amd64_gcc434/lcg/root/5.27.06b-cms14/lib -lTMVA -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMLP -lTreePlayer -lMinuit -pthread -lm -ldl -rdynamic -pthread -m64 -I/afs/in2p3.fr/grid/toolkit/cms2/slc5_amd64_gcc434/lcg/root/5.27.06b-cms14/include -L`pwd` -lToto `root-config --libs --cflags` -o ${exeFile}

#g++ Selection_miniTree.C -L/afs/in2p3.fr/grid/toolkit/cms2/slc5_amd64_gcc434/lcg/root/5.27.06b-cms14/lib -lTMVA -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMLP -lTreePlayer -lMinuit -pthread -lm -ldl -rdynamic -pthread -m64 -I/afs/in2p3.fr/grid/toolkit/cms2/slc5_amd64_gcc434/lcg/root/5.27.06b-cms14/include -lProof -lProofPlayer -L`pwd` -lToto `root-config --libs --cflags` -o ${exeFile}
#g++ Selection_miniTree.C -L`pwd` -lToto `root-config --libs --cflags` -o ${exeFile}


