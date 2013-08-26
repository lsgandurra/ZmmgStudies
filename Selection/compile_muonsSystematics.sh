#!/bin/bash
syntax="Syntax ${0} {executableName}"
if [[ -z ${1} ]]
then
  echo ${syntax}
#  exit 1
fi
exeFile=${1}

#g++ Muons_Systematics.C -L`pwd` -lToto `root-config --libs --cflags` -g -O0 -o ${exeFile}
#g++ Muons_Systematics.C `root-config --libs --cflags` -g -O0 -o ${exeFile}
#g++ Muons_Systematics.C -L`pwd` -lRooFit -lRooFitCore `root-config --libs --cflags` -g -O0 -o ${exeFile}
g++ Muons_Systematics.C -lRooFit -lRooFitCore `root-config --libs --cflags` -g -O0 -o ${exeFile}




