#ifndef FITFUNCTIONS
#define FITFUNCTIONS

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
//#include "RooCruijff.hh"

// --- c++ libraries --- //
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstring>
#include <iomanip>
#include <cmath>
#include <sys/types.h>
#include <sys/stat.h>

// --- others --- //
#include "RooCruijff.hh"
//#include "CrystalBall.C"
//#include "setTDRStyle.C"
//#include "CMSStyle.C"

using namespace RooFit;
using namespace std;

void voigtian(RooDataSet *dataset, RooDataSet *dataset2, RooRealVar &variable, RooPlot *fitFrame, RooBinning b, double rangeMin, double rangeMax, vector <double> &fitParameters);
void voigtian_surface(RooDataSet *dataset, RooDataSet *dataset2, RooRealVar &variable, RooPlot *fitFrame, RooBinning b, double rangeMin, double rangeMax, double xMaxHisto, vector <double> &fitParameters);
void cruijff(RooDataSet *dataset, RooDataSet *dataset2, RooRealVar &variable, RooPlot *fitFrame, RooBinning b, double rangeMin, double rangeMax, vector <double> &fitParameters);
void voigtianXcb(RooDataSet *dataset, RooDataSet *dataset2, RooRealVar &variable, RooPlot *fitFrame, RooBinning b, double rangeMin, double rangeMax, vector <double> &fitParameters);
void bwXcb(RooDataSet *dataset, RooDataSet *dataset2, RooRealVar &variable, RooPlot *fitFrame, RooBinning b, double rangeMin, double rangeMax, vector <double> &fitParameters);
void voigtianXgauss(RooDataSet *dataset, RooDataSet *dataset2, RooRealVar &variable, RooPlot *fitFrame, RooBinning b, double rangeMin, double rangeMax, vector <double> &fitParameters);

#endif
