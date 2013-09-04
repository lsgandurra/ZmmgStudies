# définition des cibles particulières
.PHONY: clean, mrproper
 
# désactivation des règles implicites
.SUFFIXES:
 
# définition des variables
CC = g++
#CFLAGS = -W -Wall -v -L`pwd` -lRooFitCore -lRooFit `root-config --libs --cflags`
#CFLAGS = -L`pwd` -lRooFitCore -lRooFit `root-config --libs --cflags` 
CFLAGS = `root-config --cflags --libs` -L $ROOTSYS/lib -lRooFit -lHtml -lMinuit -lRooFitCore -W 
  
# all
all: SFits.o fitFunctions.o functions.o RooCruijff.o Fit_Selection.o PvaluesPlotter.o Data_MC_Pt_Comparison.o Pvalues_vs_percent.o Mean_vs_percent.o Data_MC_r9_Comparison.o Data_MC_Mmumugamma_Comparison.o Data_MC_Var_Comparison.o GlobalPvalue.o fitFunctionSystematics.o SameStatBins.o ProfilePlotter.o Mean_vs_injectedResolution.o Data_MC_Var_Comparison_temp.o smearingPlotter.o muonsSystematics.o ProfilePlotter2.o Data_MC_PU_Comparison.o SurfaceGenerator.o Surface_fit.o ThreeBodyVsTwoBody.o CheckAllFits.o
	$(CC) SFits.o fitFunctions.o functions.o RooCruijff.o -o SFits.exe $(CFLAGS)
	$(CC) Fit_Selection.o functions.o -o Fit_Selection.exe $(CFLAGS)
	$(CC) PvaluesPlotter.o functions.o -o PvaluesPlotter.exe $(CFLAGS)
	$(CC) Data_MC_Pt_Comparison.o functions.o -o Data_MC_Pt_Comparison.exe $(CFLAGS)
	$(CC) Data_MC_r9_Comparison.o functions.o -o Data_MC_r9_Comparison.exe $(CFLAGS)
	$(CC) Data_MC_Mmumugamma_Comparison.o functions.o -o Data_MC_Mmumugamma_Comparison.exe $(CFLAGS)
	$(CC) Data_MC_Var_Comparison.o functions.o -o Data_MC_Var_Comparison.exe $(CFLAGS)
	$(CC) Pvalues_vs_percent.o functions.o -o Pvalues_vs_percent.exe $(CFLAGS)
	$(CC) Mean_vs_percent.o functions.o -o Mean_vs_percent.exe $(CFLAGS)
	$(CC) GlobalPvalue.o -o GlobalPvalue.exe $(CFLAGS)
	$(CC) fitFunctionSystematics.o RooCruijff.o functions.o -o fitFunctionSystematics.exe $(CFLAGS)
	$(CC) SameStatBins.o functions.o -o SameStatBins.exe $(CFLAGS)
	$(CC) ProfilePlotter.o functions.o -o ProfilePlotter.exe $(CFLAGS)
	$(CC) ProfilePlotter2.o functions.o -o ProfilePlotter2.exe $(CFLAGS)
	$(CC) Mean_vs_injectedResolution.o functions.o -o Mean_vs_injectedResolution.exe $(CFLAGS)
	$(CC) Data_MC_Var_Comparison_temp.o functions.o -o Data_MC_Var_Comparison_temp.exe $(CFLAGS)
	$(CC) smearingPlotter.o functions.o -o smearingPlotter.exe $(CFLAGS)
	$(CC) muonsSystematics.o functions.o -o muonsSystematics.exe $(CFLAGS)
	$(CC) Data_MC_PU_Comparison.o functions.o -o  Data_MC_PU_Comparison.exe $(CFLAGS) 
	$(CC) SurfaceGenerator.o functions.o -o SurfaceGenerator.exe $(CFLAGS)
	$(CC) Surface_fit.o fitFunctions.o functions.o RooCruijff.o -o Surface_fit.exe $(CFLAGS) 
	$(CC) ThreeBodyVsTwoBody.o functions.o -o ThreeBodyVsTwoBody.exe $(CFLAGS)
	$(CC) CheckAllFits.o functions.o -o CheckAllFits.exe $(CFLAGS) 

SFits.o: SFits.cpp fitFunctions.h functions.h RooCruijff.hh
	$(CC) -c SFits.cpp -o SFits.o $(CFLAGS)

fitFunctions.o : fitFunctions.cpp RooCruijff.hh
	$(CC) -c fitFunctions.cpp -o fitFunctions.o $(CFLAGS)

functions.o : functions.cpp
	$(CC) -c functions.cpp -o functions.o $(CFLAGS)

RooCruijff.o : RooCruijff.cc
	$(CC) -c RooCruijff.cc -o RooCruijff.o $(CFLAGS) 

Fit_Selection.o : Fit_Selection.cpp functions.h
	$(CC) -c Fit_Selection.cpp -o Fit_Selection.o $(CFLAGS)

PvaluesPlotter.o : PvaluesPlotter.cpp functions.h 
	$(CC) -c PvaluesPlotter.cpp -o PvaluesPlotter.o $(CFLAGS)
	
Data_MC_Pt_Comparison.o : Data_MC_Pt_Comparison.cpp functions.h 
	$(CC) -c Data_MC_Pt_Comparison.cpp -o Data_MC_Pt_Comparison.o $(CFLAGS)

Data_MC_r9_Comparison.o : Data_MC_r9_Comparison.cpp functions.h 
	$(CC) -c Data_MC_r9_Comparison.cpp -o Data_MC_r9_Comparison.o $(CFLAGS)

Data_MC_Mmumugamma_Comparison.o : Data_MC_Mmumugamma_Comparison.cpp functions.h
	$(CC) -c Data_MC_Mmumugamma_Comparison.cpp -o Data_MC_Mmumugamma_Comparison.o $(CFLAGS)

Data_MC_Var_Comparison.o : Data_MC_Var_Comparison.cpp functions.h
	$(CC) -c Data_MC_Var_Comparison.cpp -o Data_MC_Var_Comparison.o $(CFLAGS)

Pvalues_vs_percent.o : Pvalues_vs_percent.cpp functions.h
	$(CC) -c Pvalues_vs_percent.cpp -o Pvalues_vs_percent.o $(CFLAGS)

Mean_vs_percent.o : Mean_vs_percent.cpp functions.h
	$(CC) -c Mean_vs_percent.cpp -o Mean_vs_percent.o $(CFLAGS)

GlobalPvalue.o : GlobalPvalue.cpp 
	$(CC) -c GlobalPvalue.cpp -o GlobalPvalue.o $(CFLAGS)

fitFunctionSystematics.o : fitFunctionSystematics.cpp RooCruijff.hh functions.h
	$(CC) -c fitFunctionSystematics.cpp -o fitFunctionSystematics.o $(CFLAGS)

SameStatBins.o : SameStatBins.cpp functions.h
	$(CC) -c SameStatBins.cpp -o SameStatBins.o $(CFLAGS)

ProfilePlotter.o : ProfilePlotter.cpp functions.h
	$(CC) -c ProfilePlotter.cpp -o ProfilePlotter.o $(CFLAGS)

Mean_vs_injectedResolution.o : Mean_vs_injectedResolution.cpp functions.h
	$(CC) -c Mean_vs_injectedResolution.cpp -o Mean_vs_injectedResolution.o $(CFLAGS)

Data_MC_Var_Comparison_temp.o : Data_MC_Var_Comparison_temp.cpp functions.h
	$(CC) -c Data_MC_Var_Comparison_temp.cpp -o Data_MC_Var_Comparison_temp.o $(CFLAGS)

smearingPlotter.o : smearingPlotter.cpp functions.h
	$(CC) -c smearingPlotter.cpp -o smearingPlotter.o $(CFLAGS)

muonsSystematics.o : muonsSystematics.cpp functions.h
	$(CC) -c muonsSystematics.cpp -o muonsSystematics.o $(CFLAGS)

ProfilePlotter2.o : ProfilePlotter2.cpp functions.h
	$(CC) -c ProfilePlotter2.cpp -o ProfilePlotter2.o $(CFLAGS)

Data_MC_PU_Comparison.o : Data_MC_PU_Comparison.cpp functions.h
	$(CC) -c Data_MC_PU_Comparison.cpp -o Data_MC_PU_Comparison.o $(CFLAGS)

SurfaceGenerator.o : SurfaceGenerator.cpp functions.h
	$(CC) -c SurfaceGenerator.cpp -o SurfaceGenerator.o $(CFLAGS)
	
Surface_fit.o : Surface_fit.cpp fitFunctions.h RooCruijff.hh functions.h 
	$(CC) -c Surface_fit.cpp -o Surface_fit.o $(CFLAGS)

ThreeBodyVsTwoBody.o : ThreeBodyVsTwoBody.cpp functions.h
	$(CC) -c ThreeBodyVsTwoBody.cpp -o ThreeBodyVsTwoBody.o $(CFLAGS)

CheckAllFits.o : CheckAllFits.cpp functions.h
	$(CC) -c CheckAllFits.cpp -o CheckAllFits.o $(CFLAGS)

# clean
clean:
	rm -rf *.bak rm -rf *.o
 
# mrproper
mrproper: clean
	rm -rf Programme




