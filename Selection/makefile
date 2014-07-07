# définition des cibles particulières
.PHONY: clean, mrproper
 
# désactivation des règles implicites
.SUFFIXES:
 
# définition des variables
CC = g++
#CFLAGS = -W -Wall -v -L`pwd` -lRooFitCore -lRooFit `root-config --libs --cflags`
#CFLAGS = -L`pwd` -lRooFitCore -lRooFit `root-config --libs --cflags` 
CFLAGS = -L${ROOTSYS}lib -lTMVA -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMLP -lTreePlayer -lMinuit -pthread -lm -ldl -rdynamic -pthread -m64 -I${ROOTSYS}include -L`pwd` -lToto `root-config --libs --cflags` -g 
  
# all
all: Selection_miniTree.o Muons_miniTree_test.o corrections.o rochcor2012jan22.o rochcor_wasym_v4.o pileUp_getNvtxRECO_test.o pileUp_scaleNvtxRECO.o pileUp_getNvtxRECO.o SumOutputMinitree.o
	$(CC) Selection_miniTree.o corrections.o rochcor2012jan22.o rochcor_wasym_v4.o -o Selection_miniTree.exe $(CFLAGS)
	$(CC) Muons_miniTree_test.o corrections.o rochcor2012jan22.o rochcor_wasym_v4.o -o Muons_miniTree_test.exe $(CFLAGS) 
	$(CC) pileUp_getNvtxRECO_test.o -o pileUp_getNvtxRECO_test.exe $(CFLAGS)
	$(CC) pileUp_scaleNvtxRECO.o -o pileUp_scaleNvtxRECO.exe $(CFLAGS)
	$(CC) pileUp_getNvtxRECO.o -o pileUp_getNvtxRECO.exe $(CFLAGS)
	$(CC) SumOutputMinitree.o -o SumOutputMinitree.exe $(CFLAGS)

Selection_miniTree.o: Selection_miniTree.C corrections.h rochcor2012jan22.h
	$(CC) -c Selection_miniTree.C -o Selection_miniTree.o $(CFLAGS)

Muons_miniTree_test.o: Muons_miniTree_test.C corrections.h rochcor2012jan22.h
	$(CC) -c Muons_miniTree_test.C -o Muons_miniTree_test.o $(CFLAGS)

corrections.o : corrections.C
	$(CC) -c corrections.C -o corrections.o $(CFLAGS)

rochcor2012jan22.o : rochcor2012jan22.C
	$(CC) -c rochcor2012jan22.C -o rochcor2012jan22.o $(CFLAGS)

rochcor_wasym_v4.o: rochcor_wasym_v4.C
	$(CC) -c rochcor_wasym_v4.C -o rochcor_wasym_v4.o $(CFLAGS)

pileUp_getNvtxRECO_test.o : pileUp_getNvtxRECO_test.cpp
	$(CC) -c pileUp_getNvtxRECO_test.cpp -o pileUp_getNvtxRECO_test.o $(CFLAGS)

pileUp_scaleNvtxRECO.o : pileUp_scaleNvtxRECO.C 
	$(CC) -c pileUp_scaleNvtxRECO.C -o pileUp_scaleNvtxRECO.o $(CFLAGS)  

pileUp_getNvtxRECO.o : pileUp_getNvtxRECO.C
	$(CC) -c pileUp_getNvtxRECO.C -o pileUp_getNvtxRECO.o $(CFLAGS)


SumOutputMinitree.o : SumOutputMinitree.cpp
	$(CC) -c SumOutputMinitree.cpp -o SumOutputMinitree.o $(CFLAGS)
 
# clean
clean:
	rm -rf *.bak rm -rf *.o
 
# mrproper
mrproper: clean
	rm -rf Programme


