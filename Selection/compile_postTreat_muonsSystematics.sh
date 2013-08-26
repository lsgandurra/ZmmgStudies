#!/bin/bash
syntax="Syntax ${0} {executableName}"
if [[ -z ${1} ]]
then
  echo ${syntax}
#  exit 1
fi
exeFile=${1}

g++ postTreat_Muons_Systematics.C -lRooFit -lRooFitCore `root-config --libs --cflags` -g -O0 -o ${exeFile}



