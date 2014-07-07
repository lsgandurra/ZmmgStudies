#ifndef FUNCTIONS
#define FUNCTIONS

// --- Root libraries --- //
#include "TF1.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TH2F.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TBranch.h"
#include "TChain.h"
#include "TFile.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TString.h"
#include "TBits.h"
#include "TMath.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TPaveStats.h"
#include "TPaveText.h"
#include "TMultiGraph.h"
#include "TGraphAsymmErrors.h"
#include "TLatex.h"
#include "TGraphErrors.h"
#include "TPad.h"
#include "TGraph.h"
#include "TText.h"
#include "TLine.h"
#include "TGaxis.h"
#include "TROOT.h"
#include "TRint.h"
#include "TDirectory.h"
#include "TProfile.h"
#include "TProfile3D.h"
#include "TPaletteAxis.h"

// --- Roofit libraries --- //
#include "RooAbsPdf.h"
#include "RooCBShape.h"
#include "RooLandau.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooHistPdf.h"
#include "RooHist.h"
#include "RooMsgService.h"
#include "RooNLLVar.h"
#include "RooPlot.h"
#include "RooRandom.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "RooFFTConvPdf.h"
#include "RooLognormal.h"
#include "RooGaussian.h"
#include "RooGamma.h"
#include "RooBifurGauss.h"
#include "RooAddPdf.h"
#include "RooGenericPdf.h"
#include "RooProdPdf.h"
#include "RooDerivative.h"
#include "RooFitResult.h"
#include "RooMCStudy.h"
#include "RooChi2MCSModule.h"
#include "RooConstVar.h"
#include "RooChebychev.h"
#include "RooKeysPdf.h"
#include "RooVoigtian.h"
#include "RooBreitWigner.h"
#include "RooBinning.h"
#include "RooStringVar.h"
//#include "RooCruijff.hh"

// --- c++ libraries --- //
#include <utility>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstring>
#include <iomanip>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

// --- others --- //
//#include "CrystalBall.C"
//#include "setTDRStyle.C"
//#include "CMSStyle.C"

using namespace RooFit;
using namespace std;

bool FileExists(const char* FileName);
string doubleToString(double x);
int rowsNumberInFile(string filename);
void plotsRecording(string directoryName, string fileName, TCanvas * c1);
//void chi2Chain(char * buffer, double chi2);

Double_t effSigma(TH1 * hist);
//double sigmaR(TF1* function, double xmin, double xmax);
//double sigmaL(TF1* function, double xmin, double xmax);

void rangeEstimator(double percentage, TChain * chain, TString cut, double &minRange, double &maxRange, string variableName, float variable, vector <double> &fitParameters);
//void rangeEstimator(double percentage, double centralValue, TChain * chain, int Endcaps, double * MinRange, double * MaxRange, double Varmin, double Varmax);
//void rangeEstimator2(double percentage, TChain * chain, int Endcaps, double * MinRange, double * MaxRange, double Varmin, double Varmax);
//void symetricRangeEstimator(TChain * chain, double centralValue, double * MinRange, double * MaxRange, double Entries, double percentage, TString temp);
//void symetricRangeEstimator2(TChain * chain, double lastX, double * MinRange, double * MaxRange, double Entries, double percentage, TString temp);
//void symetricRangeEstimator3(TChain * chain, double centralValue, double * MinRange, double * MaxRange, double Entries, double percentage, TString temp);


Double_t chiSquare(RooPlot* plot, char* pdfname, char* histname, vector <double> &fitParameters, int nbPar);
RooHist* residHist(RooPlot* plot, char *histname, char* curvename, bool normalize, string recordingDirectory, int iteration);
RooHist* pullHist(RooPlot* plot_, char* histname, char* pdfname, string dossierSauvegardePull, int iteration);

#endif
