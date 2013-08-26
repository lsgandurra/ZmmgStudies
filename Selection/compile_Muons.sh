#!/bin/bash
syntax="Syntax ${0} {executableName}"
if [[ -z ${1} ]]
then
  echo ${syntax}
#  exit 1
fi
exeFile=${1}


#g++ Selection_miniTree.C -L${ROOTSYS}lib -lTMVA -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMLP -lTreePlayer -lMinuit -pthread -lm -ldl -rdynamic -pthread -m64 -I${ROOTSYS}include -L`pwd` -lToto `root-config --libs --cflags` -o ${exeFile}

#g++ -L`pwd` `root-config --libs --cflags` -c rochcor_v2.C
#g++ -L`pwd` -lToto `root-config --libs --cflags` -c Muons_miniTree.C
#g++ -L`pwd` -lToto `root-config --libs --cflags` Muons_miniTree.C rochcor_v2.o -o ${exeFile}
g++ Muons_miniTree.C -L`pwd` -lToto `root-config --libs --cflags` -o ${exeFile}



