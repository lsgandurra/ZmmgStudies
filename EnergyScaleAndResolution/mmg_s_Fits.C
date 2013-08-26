//#include "TROOT.h "
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
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
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstring>
#include <iomanip>
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
#include "TROOT.h"
#include "TRint.h"
#include "TDirectory.h"
#include "RooFitResult.h"
#include "RooMCStudy.h"
#include "RooChi2MCSModule.h"
#include "RooConstVar.h"
#include "RooChebychev.h"
#include "RooKeysPdf.h"
#include "RooVoigtian.h"
#include "CrystalBall.C"
#include "RooBreitWigner.h"
#include "RooBinning.h"
#include "setTDRStyle.C"
#include "CMSStyle.C"
#include <cmath>
#include <sys/types.h>
#include <sys/stat.h>
#include "RooCruijff.hh"
//#include "PValue.h"
#pragma optimize 0

using namespace RooFit;
using namespace std;

string DoubleToString(double x);

Double_t myfunction(Double_t *x, Double_t *par);
Double_t myfunction2(Double_t *x, Double_t *par);
void myfunc(int EndCaps);
void myfit(TH1D *h1);

float PtCor(float brem, int EndCaps, int ParametersVersion); //ParametersVersion = 1,2,3,4 (=electrons old, new, photons old, new);
Double_t effSigma(TH1 * hist);
void ChaineChi2(char * buffer, double chi2);
double SigmaR(TF1* crystalBall, double Xmin, double Xmax);
double SigmaL(TF1* crystalBall, double Xmin, double Xmax);
float fEta(float eta);
void CrystalBallMethode(double * mean_value, double * mean_error, double * sigma_value, double * ChiSquare, double * DegreesOfFreedom, double param[5], TH1D Data);

void RooCrystalBall(double * mean_value, double * mean_error, Double_t * sigmaEff_value, double * sigma_value, double * sigma_value_error, double * sigmaR_value, double * sigmaL_value, double * ChiSquare, TH1D * hh, int j, int EndCaps, TString temp, TTree* Tree_Data, double mean, double rms, TCanvas * c2, RooPlot* Erawframe, TF1 * f, string nomDossier, string nomFichier, double RangeMin, double RangeMax, double * JanChi2, double * DegreesOfFreedom, double * pValue, int r9sup, double MinVar, double MaxVar, string variableX, int isMC);

void RooLogNormal(double * mean_value, double * mean_error, Double_t * sigmaEff_value, double * sigma_value, double * sigma_value_error, double * sigmaR_value, double * sigmaL_value, double * ChiSquare, TH1D * hh, int j, int EndCaps, TString temp, TTree* Tree_Data, double mean, double rms, TCanvas * c2, RooPlot* Erawframe, TF1 * f, string nomDossier, string nomFichier, double RangeMin, double RangeMax, double * JanChi2, double * DegreesOfFreedom, double * pValue, int r9sup, double MinVar, double MaxVar, string variableX, int isMC);

void RooGauss(double * mean_value, double * mean_error, Double_t * sigmaEff_value, double * sigma_value, double * sigma_value_error, double * sigmaR_value, double * sigmaL_value, double * ChiSquare, TH1D * hh, int j, int EndCaps, TString temp, TTree* Tree_Data, double mean, double rms, TCanvas * c2, RooPlot* Erawframe, TF1 * f, string nomDossier, string nomFichier, double RangeMin, double RangeMax, double * JanChi2, double * DegreesOfFreedom, double * pValue, int r9sup, double MinVar, double MaxVar, string variableX, int isMC);

void RooGamma2(double * mean_value, double * mean_error, Double_t * sigmaEff_value, double * sigma_value, double * sigma_value_error, double * sigmaR_value, double * sigmaL_value, double * ChiSquare, TH1D * hh, int j, int EndCaps, TString temp, TTree* Tree_Data, double mean, double rms, TCanvas * c2, RooPlot* Erawframe, TF1 * f, string nomDossier, string nomFichier, double RangeMin, double RangeMax, double * JanChi2, double * DegreesOfFreedom, double * pValue, int r9sup, double MinVar, double MaxVar, string variableX, int isMC);

void RooBifurcatedGauss(double * mean_value, double * mean_error, Double_t * sigmaEff_value, double * sigma_value, double * sigma_value_error, double * sigmaR_value, double * sigmaL_value, double * ChiSquare, double * minLogLikelihood, double * differenceLogLikelihood, TH1D * hh, int j, int EndCaps, TString temp, TTree* Tree_Data, double mean, double rms, TCanvas * c2, RooPlot* Erawframe, TF1 * f, string nomDossier, string nomFichier, double RangeMin, double RangeMax, double * JanChi2, double * DegreesOfFreedom, double * pValue, int r9sup, double MinVar, double MaxVar, string variableX, int isMC);

void RooSumGauss(double * mean_value, double * mean_error, Double_t * sigmaEff_value, double * sigma_value, double * sigma_value_error, double * sigmaR_value, double * sigmaL_value, double * ChiSquare, TH1D * hh, int j, int EndCaps, TString temp, TTree* Tree_Data, double mean, double rms, TCanvas * c2, RooPlot* Erawframe, TF1 * f, string nomDossier, string nomFichier, double RangeMin, double RangeMax, double * JanChi2, double * DegreesOfFreedom, double * pValue, int r9sup, double MinVar, double MaxVar, string variableX, int isMC);

void RooGenericPDF(double * mean_value, double * mean_error, Double_t * sigmaEff_value, double * sigma_value, double * sigma_value_error, double * sigmaR_value, double * sigmaL_value, double * ChiSquare, TH1D * hh, int j, int EndCaps, TString temp, TTree* Tree_Data, double mean, double rms, TCanvas * c2, RooPlot* Erawframe, TF1 * f, string nomDossier, string nomFichier);

void RooLandau2(double * mean_value, double * mean_error, Double_t * sigmaEff_value, double * sigma_value, double * sigma_value_error, double * sigmaR_value, double * sigmaL_value, double * ChiSquare, TH1D * hh, int j, int EndCaps, TString temp, TTree* Tree_Data, double mean, double rms, TCanvas * c2, RooPlot* Erawframe, TF1 * f, string nomDossier, string nomFichier, double RangeMin, double RangeMax, double * JanChi2, double * DegreesOfFreedom, double * pValue, int r9sup, double MinVar, double MaxVar, string variableX, int isMC);

void RooLandauConvGaussian(double * mean_value, double * mean_error, Double_t * sigmaEff_value, double * sigma_value, double * sigma_value_error, double * sigmaR_value, double * sigmaL_value, double * ChiSquare, TH1D * hh, int j, int EndCaps, TString temp, TTree* Tree_Data, double mean, double rms, TCanvas * c2, RooPlot* Erawframe, TF1 * f, string nomDossier, string nomFichier, double RangeMin, double RangeMax, double * JanChi2, double * DegreesOfFreedom, double * pValue, int r9sup, double MinVar, double MaxVar, string variableX, int isMC);

void RooKernel(double * mean_value, double * mean_error, Double_t * sigmaEff_value, double * sigma_value, double * sigma_value_error, double * sigmaR_value, double * sigmaL_value, double * ChiSquare, double * minLogLikelihood, double * differenceLogLikelihood, TH1D * hh, int j, int EndCaps, TString temp, TTree* Tree_Data, double mean, double rms, TCanvas * c2, RooPlot* Erawframe, TF1 * f, string nomDossier, string nomFichier, double RangeMin, double RangeMax, double * JanChi2, double * DegreesOfFreedom, double * pValue, int r9sup, double MinVar, double MaxVar, string variableX, int isMC);

void RooVoigtian2(double * mean_value, double * mean_error, Double_t * sigmaEff_value, double * sigma_value, double * sigma_value_error, double * sigmaR_value, double * sigmaL_value, double * ChiSquare, double * minLogLikelihood, double * differenceLogLikelihood, TH1D * hh, int j, int EndCaps, TString temp, TTree* Tree_Data, double mean, double rms, TCanvas * c2, RooPlot* Erawframe, TF1 * f, string nomDossier, string nomFichier, double RangeMin, double RangeMax, double * JanChi2, double * DegreesOfFreedom, double * pValue, int r9sup, double MinVar, double MaxVar, string variableX, int isMC);


void Cruijff(double * mean_value, double * mean_error, Double_t * sigmaEff_value, double * sigma_value, double * sigma_value_error, double * sigmaR_value, double * sigmaL_value, double * ChiSquare, double * minLogLikelihood, double * differenceLogLikelihood, TH1D * hh, int j, int EndCaps, TString temp, TTree* Tree_Data, double mean, double rms, TCanvas * c2, RooPlot* Erawframe, TF1 * f, string nomDossier, string nomFichier, double RangeMin, double RangeMax, double * JanChi2, double * DegreesOfFreedom, double * pValue, int r9sup, double MinVar, double MaxVar, string variableX, int isMC);

void roofit_bwcb_fft2(double * mean_value, double * mean_error, Double_t * sigmaEff_value, double * sigma_value, double * sigma_value_error, double * sigmaR_value, double * sigmaL_value, double * ChiSquare, double * minLogLikelihood, double * differenceLogLikelihood, TH1D * hh, int j, int EndCaps, TString temp, TTree* Tree_Data, double mean, double rms, TCanvas * c2, RooPlot* Erawframe, TF1 * f, string nomDossier, string nomFichier, double RangeMin, double RangeMax, double * JanChi2, double * DegreesOfFreedom, double * pValue, int r9sup, double MinVar, double MaxVar, string variableX, int isMC);

void enregistrementPlots(string nomDossier, string nomFichier, int EndCaps, int iteration, TCanvas * c1);

void RangeEstimator(double pourcentage, double centralValue, TChain * chain, int Endcaps, double * MinRange, double * MaxRange, double Varmin, double Varmax);

void RangeEstimator2(double pourcentage, TChain * chain, int Endcaps, double * MinRange, double * MaxRange, double Varmin, double Varmax);

void RangeEstimator3(double pourcentage, TChain * chain, TString temp, int Endcaps, double * MinRange, double * MaxRange);

void SymetricRangeEstimator(TChain * chain, double centralValue, double * MinRange, double * MaxRange, double Entries, double pourcentage, TString temp);

void SymetricRangeEstimator2(TChain * chain, double lastX, double * MinRange, double * MaxRange, double Entries, double pourcentage, TString temp);

void SymetricRangeEstimator3(TChain * chain, double centralValue, double * MinRange, double * MaxRange, double Entries, double pourcentage, TString temp);

Double_t chiSquare(RooPlot* plot_, char* pdfname, char* histname, int nFitParam, double* JanChi2, double* DegreesOfFreedom, double* pValue, int* fewBins);

RooHist* residHist(RooPlot* plot_, char *histname, char* curvename, bool normalize, string dossierSauvegardePull, int iteration);

RooHist* pullHist(RooPlot* plot_, char* histname, char* pdfname, string dossierSauvegardePull, int iteration) { return residHist(plot_, histname, pdfname, true, dossierSauvegardePull, iteration); }


int NbLignesFichier(string fichier);

string DoubleToString(double x)
{

        std::string s;
        {
                std::ostringstream oss; 
                oss << x;
                s = oss.str();
        }
        std::cout << "x = " << x << " s = " << s << std::endl;

        return s;
}


Double_t myfunction(Double_t *x, Double_t *par)
{
        float bremLowThr  = 1.10000002384185791015625;
        float bremHighThr = 8.0;
        if ( x[0] < bremLowThr  ) x[0] = bremLowThr;
        if ( x[0] > bremHighThr ) x[0] = bremHighThr;

        float threshold = par[4];

        float y = par[0]*threshold*threshold + par[1]*threshold + par[2];
        float yprime = 2*par[0]*threshold + par[1];
        float a = par[3];
        float b = yprime - 2*a*threshold;
        float c = y - a*threshold*threshold - b*threshold;

        float bremCorrection = 1.0;
        if ( x[0] < threshold ) bremCorrection = par[0]*x[0]*x[0] + par[1]*x[0] + par[2];
        else bremCorrection = a*x[0]*x[0] + b*x[0] + c;

        return bremCorrection;

}

Double_t myfunction2(Double_t *x, Double_t *par)
{
        float bremLowThr  = 0.89999997615814208984375;
  	float bremHighThr = 6.5;
	if ( x[0] < bremLowThr  ) x[0] = bremLowThr;
        if ( x[0] > bremHighThr ) x[0] = bremHighThr;

        float threshold = par[4];

        float y = par[0]*threshold*threshold + par[1]*threshold + par[2];
        float yprime = 2*par[0]*threshold + par[1];
        float a = par[3];
        float b = yprime - 2*a*threshold;
        float c = y - a*threshold*threshold - b*threshold;

        float bremCorrection = 1.0;
        if ( x[0] < threshold ) bremCorrection = par[0]*x[0]*x[0] + par[1]*x[0] + par[2];
        else bremCorrection = a*x[0]*x[0] + b*x[0] + c;

        return bremCorrection;

}


void myfunc(int Endcaps)
{

	TF1 *f1;
      if(Endcaps == 0) 
      {	
		TF1 *f1 = new TF1("myfunc",myfunction,0,100,5);
		//f1->Draw();
      }
      if(Endcaps == 1) 
      {		
		TF1 *f1 = new TF1("myfunc",myfunction2,0,100,5);
      		//f1->Draw();
      }
	f1->Draw();

	//f1->SetParameters(-0.002362,0.004314,1.001,0.0003413,3.124);
      //f1->SetParNames("constant","coefficient");
}

void myfit(TH1D *h1)
{


      //TH1F *h1=new TH1F("h1","test",100,0,12);
      //h1->FillRandom("myfunc",20000);
      TF1 *f1=(TF1*)gROOT->GetFunction("myfunc");
      //f1->SetParameters(800,1);
      f1->SetFillColor(19);
      f1->SetFillStyle(0);
      f1->SetLineWidth(3);


	TPaveStats *ptstats = new TPaveStats(0.6845638,0.7226277,0.9865772,0.9912587,"brNDC");
        ptstats->SetName("stats");
        ptstats->SetBorderSize(2);
        ptstats->SetFillColor(kWhite);
        ptstats->SetTextColor(19);
        ptstats->SetTextAlign(12);
        ptstats->SetTextFont(42);
        ptstats->SetOptStat(0);
        ptstats->SetOptFit(11111111);
        ptstats->Draw();


      h1->Fit("myfunc");
}




float PtCor(float brem, int EndCaps, int ParametersVersion) 
{
  // brem == phiWidth/etaWidth of the SuperCluster 
  // e  == energy of the SuperCluster 
  // first parabola (for br > threshold) 
  // p0 + p1*x + p2*x^2 
  // second parabola (for br <= threshold) 
  // ax^2 + bx + c, make y and y' the same in threshold 
  // y = p0 + p1*threshold + p2*threshold^2  
  // yprime = p1 + 2*p2*threshold 
  // a = p3 
  // b = yprime - 2*a*threshold 
  // c = y - a*threshold^2 - b*threshold 
  float p0 = 0;
  float p1 = 0;
  float p2 = 0;
  float p3 = 0;
  float p4 = 0;

/*
  int offset;
  if ( EndCaps == 0 ) offset = 0; //Barrel
  else if ( EndCaps == 1 ) offset = 20; //End caps

  else {
    // not supported, produce no correction
    return 1.0;
  }
*/

  //Make No Corrections if brem is invalid! 
	//cout<<endl<<setprecision( 10 )<<"brem in function= "<<brem<<endl;
  if ( brem == 0 ) return 1.0;

 if(EndCaps == 1) //End caps
 {
  float bremLowThr  = 0.89999997615814208984375;
  float bremHighThr = 6.5;
  if ( brem < bremLowThr  ) brem = bremLowThr;
  if ( brem > bremHighThr ) brem = bremHighThr;

  if(ParametersVersion == 1)
  {
  	p0 = -0.07945;
	p1 = 0.1298;
	p2 = 0.9147;
	p3 = -0.001565;
	p4 = 0.9;
  }

  if(ParametersVersion == 2)
  {
  	p0 = -0.12139999866485595703125;
	p1 = 0.2362000048160552978515625;
	p2 = 0.884700000286102294921875;
	p3 = -0.00193000002764165401458740234375;
	p4 = 1.05700004100799560546875;

  }

 if(ParametersVersion == 3)
  {
	p0 = -0.0728;
	p1 = 0.1323;
	p2 = 0.9201;
	p3 = 0.002501;
	p4 = 1.118;
  }

 if(ParametersVersion == 4)
  {
	p0 = -0.07667;
	p1 = 0.1407;
	p2 = 0.9157;
	p3 = 0.00251;
	p4 = 1.117;
  }



 }

if(EndCaps == 0) //Barrel
 {
  float bremLowThr  = 1.10000002384185791015625;
  float bremHighThr = 8.0;
  if ( brem < bremLowThr  ) brem = bremLowThr;
  if ( brem > bremHighThr ) brem = bremHighThr;

	
  if(ParametersVersion == 1)
  {
	p0 = -0.05185;
	p1 = 0.1354;
	p2 = 0.9165;
	p3 = -0.0005626;
	p4 = 1.385;    
  }

  if(ParametersVersion == 2)
  {
	p0 = -0.05289;
	p1 = 0.1374;
	p2 = 0.9141;
	p3 = -0.000669;
	p4 = 1.38;
  }

 if(ParametersVersion == 3)
  {
	p0 = -0.0625;
	p1 = 0.1331;
	p2 = 0.9329;
	p3 = -0.0007823;
	p4 = 1.1;
  }

 if(ParametersVersion == 4)
  {
	p0 = -0.002362;
	p1 = 0.004314;
	p2 = 1.001;
	p3 = 0.0003413;
	p4 = 3.124;
  }

if(ParametersVersion == 5)
  {
        p0 = -0.0004081;
        p1 = -0.005385;
        p2 = 1.008;
        p3 = -0.02381;
        p4 = 123.4;
  }

 }

  float threshold = p4;

  float y = p0*threshold*threshold + p1*threshold + p2;
  float yprime = 2*p0*threshold + p1;
  float a = p3;
  float b = yprime - 2*a*threshold;
  float c = y - a*threshold*threshold - b*threshold;

  float bremCorrection = 1.0;
  if ( brem < threshold ) bremCorrection = p0*brem*brem + p1*brem + p2;
  else bremCorrection = a*brem*brem + b*brem + c;

  return bremCorrection;
}




Double_t effSigma(TH1 * hist)
{

  TAxis *xaxis = hist->GetXaxis();
  Int_t nb = xaxis->GetNbins();
  if(nb < 10) {
    cout << "effsigma: Not a valid histo. nbins = " << nb << endl;
    return 0.;
  }
  
  Double_t bwid = xaxis->GetBinWidth(1);
  if(bwid == 0) {
    cout << "effsigma: Not a valid histo. bwid = " << bwid << endl;
    return 0.;
  }
  Double_t xmax = xaxis->GetXmax();
  Double_t xmin = xaxis->GetXmin();
  Double_t ave = hist->GetMean();
  Double_t rms = hist->GetRMS();

  Double_t total=0.;
  for(Int_t i=0; i<nb+2; i++) {
    total+=hist->GetBinContent(i);
  }
  if(total < 100.) {
    cout << "effsigma: Too few entries " << total << endl;
    return 0.;
  }
  Int_t ierr=0;
  Int_t ismin=999;
  
  Double_t rlim=0.683*total;
  Int_t nrms=rms/(bwid);    // Set scan size to +/- rms
  if(nrms > nb/10) nrms=nb/10; // Could be tuned...

  Double_t widmin=9999999.;
  for(Int_t iscan=-nrms;iscan<nrms+1;iscan++) { // Scan window centre
    Int_t ibm=(ave-xmin)/bwid+1+iscan;
    Double_t x=(ibm-0.5)*bwid+xmin;
    Double_t xj=x;
    Double_t xk=x;
    Int_t jbm=ibm;
    Int_t kbm=ibm;
    Double_t bin=hist->GetBinContent(ibm);
    total=bin;
    for(Int_t j=1;j<nb;j++){
      if(jbm < nb) {
        jbm++;
        xj+=bwid;
        bin=hist->GetBinContent(jbm);
        total+=bin;
        if(total > rlim) break;
      }
      else ierr=1;
      if(kbm > 0) {
        kbm--;
        xk-=bwid;
        bin=hist->GetBinContent(kbm);
        total+=bin;
        if(total > rlim) break;
      }
      else ierr=1;
    }
    Double_t dxf=(total-rlim)*bwid/bin;
    Double_t wid=(xj-xk+bwid-dxf)*0.5;
    if(wid < widmin) {
      widmin=wid;
      ismin=iscan;
    }   
  }
  if(ismin == nrms || ismin == -nrms) ierr=3;
  if(ierr != 0) cout << "effsigma: Error of type " << ierr << endl;
  
  return widmin;
  
}




void ChaineChi2(char * buffer, double chi2)
{
	
	string chaine1("#chi^{2} / ndf = ");
        string chi2Chaine;
        string chaineFinale;
	
	ostringstream DoubleToString;
        DoubleToString << chi2;
        chi2Chaine = DoubleToString.str();
        chaineFinale = chaine1 + chi2Chaine;

        //size_t size = chaineFinale.size() + 1; 
        //char * buffer = new char[ size ];
        //strncpy( buffer, chaineFinale.c_str(), size );
	strncpy( buffer, chaineFinale.c_str(), 25 );

}

double SigmaR(TF1* crystalBall, double Xmin, double Xmax)
{

	double y = crystalBall->GetMaximum(Xmin, Xmax) * 1.0/exp(1.0);
	double MaxX = crystalBall->GetMaximumX(Xmin, Xmax);
	double sigmaR = crystalBall->GetX(y, MaxX, Xmax) - MaxX;

	return sigmaR;
}

double SigmaL(TF1* crystalBall, double Xmin, double Xmax)
{

        double y = crystalBall->GetMaximum(Xmin, Xmax) * 1.0/exp(1.0);
        double MaxX = crystalBall->GetMaximumX(Xmin, Xmax);
	double sigmaL = MaxX - crystalBall->GetX(y, Xmin, MaxX);

        return sigmaL;
}

float fEta(float eta)
{

  
  float ieta = fabs(eta)*(5/0.087);
  float p0 = 40.2198;  // should be 40.2198
  float p1 = -3.03103e-6;  // should be -3.03103e-6
  float feta = 1;

  if ( ieta < p0 || fabs(eta) > 1.4442 ) feta = 1.0;
  else feta = 1.0/(1.0 + p1*(ieta-p0)*(ieta-p0));

  //correctedEnergy = energy/(1.0 + p1*(ieta-p0)*(ieta-p0));
  return feta;

}


void CrystalBallMethode(double * mean_value, double * mean_error, double * sigma_value, double * ChiSquare, double * DegreesOfFreedom, double param[5], TH1D Data)
{
	TCanvas *cf = new TCanvas("cf", "cf",0,0,600,600);

	int b0 = 0;
        int b1 = 0;
        int b3 = 0;
        int b4 = 0;
        double fit0 = 10;
        double fit1 = 0.6;
        double fit2 = 5;
        double fit3 = 1;
        double fit4 = 0.02;
	
	double ChiSquareBest = 10000000;
	double mean_errorBest = 10000000;
	*mean_error = 10000000;	

	double fit0Temp = 0;
	double fit1Temp = 0;
	double fit2Temp = 0;
	double fit3Temp = 0;
	double fit4Temp = 0;

	double fit0Best = 0;
	double fit1Best = 0;
	double fit2Best = 0;
	double fit3Best = 0;
	double fit4Best = 0;
	
	ChiSquareBest = 10000000;
        mean_errorBest = 10000000;

	for(int r = 0; r<400; r++)
       	{
               	cout<<endl<<"r = "<<r<<endl;
               	//fit0 = r;
		fit0Temp = r;
               	Data.Draw("");
		TF1* func = new TF1("func",CrystalBall,0.65, 1.1,5);
               	func->SetParameters(fit0Temp, fit1, fit2, fit3, fit4);
               	func->SetLineColor(3);
         	func->SetLineWidth(3);
               	Data.Fit(func);
		*mean_error  = func->GetParError(3);
               	*ChiSquare = func->GetChisquare();
		*DegreesOfFreedom = func->GetNDF();
               	if(*DegreesOfFreedom != 0)
                {
			cout<<endl<<"ChiSquare / DegreesOfFreedom = "<<*ChiSquare / *DegreesOfFreedom <<endl;
			if( (*ChiSquare / *DegreesOfFreedom) < ChiSquareBest)
                	{
				//if(mean_error <= mean_errorBest)
				if(*mean_error <0.01)
				{
                       			ChiSquareBest = *ChiSquare / *DegreesOfFreedom;
                       			fit0Best = r;
					mean_errorBest = *mean_error;
                       			cout<<endl<<"fit0Best = "<<fit0Best<<endl;
                       			cout<<endl<<"ChiSquareBest = "<<ChiSquareBest<<endl;
					b0 = 1;
				}
                	}
		}
               	func->Delete();
        }

        //cout<<endl<<"fit0Best = "<<fit0Best<<endl;
        //cout<<endl<<"ChiSquareBest = "<<ChiSquareBest<<endl;
	if(b0 == 1) fit0 = fit0Best;
        ChiSquareBest = 10000000;	
	mean_errorBest = 10000000;

	for(double s = -1.0; s<1; s+=0.05)
       	{
            	cout<<endl<<"s = "<<s<<endl;
               	//fit1 = s;
		fit1Temp = s;
               	Data.Draw("");
               	TF1* func = new TF1("func",CrystalBall,0.65, 1.1,5);
               	func->SetParameters(fit0, fit1Temp, fit2, fit3, fit4);
               	func->SetLineColor(3);
               	func->SetLineWidth(3);
               	Data.Fit(func);
		*mean_error  = func->GetParError(3);
                *ChiSquare = func->GetChisquare();
		*DegreesOfFreedom = func->GetNDF();
		if(*DegreesOfFreedom != 0)
                {
			cout<<endl<<"ChiSquare / DegreesOfFreedom = "<<*ChiSquare / *DegreesOfFreedom <<endl;
                        if( (*ChiSquare / *DegreesOfFreedom) < ChiSquareBest)
                        {
                                //if(mean_error <= mean_errorBest)
                                if(*mean_error <0.01)
				{
                                        ChiSquareBest = *ChiSquare / *DegreesOfFreedom;
                                        fit1Best = s;
					mean_errorBest = *mean_error;
                                        cout<<endl<<"fit1Best = "<<fit1Best<<endl;
                                        cout<<endl<<"ChiSquareBest = "<<ChiSquareBest<<endl;
                                      	b1 = 1;
				}
                        }
		}
               	func->Delete();
        }
        if(b1 == 1) fit1 = fit1Best;
        ChiSquareBest = 10000000;
	mean_errorBest = 10000000;

	for(double t = 0.900; t<1.010; t+=0.001)
        {
               	cout<<endl<<"t = "<<t<<endl;
               	//fit3 = t;
		fit3Temp = t;
               	Data.Draw("");
               	TF1* func = new TF1("func",CrystalBall,0.65, 1.1,5);
               	func->SetParameters(fit0, fit1, fit2, fit3Temp, fit4);
               	func->SetLineColor(3);
               	func->SetLineWidth(3);
               	Data.Fit(func);
		*mean_error  = func->GetParError(3);
                *ChiSquare = func->GetChisquare();
		*DegreesOfFreedom = func->GetNDF();
                if(*DegreesOfFreedom != 0)
                {
			cout<<endl<<"ChiSquare / DegreesOfFreedom = "<<*ChiSquare / *DegreesOfFreedom <<endl;
			if( (*ChiSquare / *DegreesOfFreedom) < ChiSquareBest)
                        {
                                //if(mean_error <= mean_errorBest)
                                if(*mean_error <0.01)
				{
                                        ChiSquareBest = *ChiSquare / *DegreesOfFreedom;
                                        fit3Best = t;
					mean_errorBest = *mean_error;
                                        cout<<endl<<"fit3Best = "<<fit3Best<<endl;
                                        cout<<endl<<"ChiSquareBest = "<<ChiSquareBest<<endl;
                                       	b3 = 1;
				}
                        }
		}				
                func->Delete();

        }
        if(b3 == 1) fit3 = fit3Best;
        ChiSquareBest = 10000000;
	mean_errorBest = 10000000;

	for(double u = 0.001; u<0.1; u+=0.001)
        {
               	cout<<endl<<"u = "<<u<<endl;
               	//fit4 = u;
		fit4Temp = u;
               	Data.Draw("");
               	TF1* func = new TF1("func",CrystalBall,0.65, 1.1,5);
               	func->SetParameters(fit0, fit1, fit2, fit3, fit4Temp);
               	func->SetLineColor(3);
               	func->SetLineWidth(3);
               	Data.Fit(func);
		*mean_error  = func->GetParError(3);
                *ChiSquare = func->GetChisquare();
		*DegreesOfFreedom = func->GetNDF();
		if(*DegreesOfFreedom != 0)
                {
			cout<<endl<<"ChiSquare / DegreesOfFreedom = "<<*ChiSquare / *DegreesOfFreedom <<endl;
                        if( (*ChiSquare / *DegreesOfFreedom) < ChiSquareBest)
                        {
                	        //if(mean_error <= mean_errorBest)
                                if(*mean_error <0.01)
				{
                                        ChiSquareBest = *ChiSquare / *DegreesOfFreedom;
                                        fit4Best = u;
					mean_errorBest = *mean_error;
                                        cout<<endl<<"fit4Best = "<<fit4Best<<endl;
                                        cout<<endl<<"ChiSquareBest = "<<ChiSquareBest<<endl;
					b4 = 1;
                                }
                        }
		}
                func->Delete();
        }
        if( b4 == 1 ) fit4 = fit4Best;
	cout<<endl<<"fit0Best = "<<fit0Best<<endl<<"fit1Best = "<<fit1Best<<endl<<"fit3Best = "<<fit3Best<<endl<<"fit4Best = "<<fit4Best<<endl;

	cf->Clear();
       	Data.Draw("");
       	TF1* func = new TF1("func",CrystalBall,0.65, 1.1,5);
       	func->SetParameters(fit0, fit1, fit2, fit3, fit4);
	func->SetLineColor(3);
        func->SetLineWidth(3);
       	Data.Fit(func);
       	*mean_value  = func->GetParameter(3);
       	*mean_error  = func->GetParError(3); 
       	*sigma_value = func->GetParameter(4);
	*ChiSquare = func->GetChisquare();
	*DegreesOfFreedom = func->GetNDF();	
	

	param[0] = fit0;
	param[1] = fit1;
	param[2] = fit2;
	param[3] = fit3;
	param[4] = fit4;

	
	cf->Clear();
	func->Delete();
	//cf->Delete();
	cf->Close();
}

void RooCrystalBall(double * mean_value, double * mean_error, Double_t * sigmaEff_value, double * sigma_value, double * sigma_value_error, double * sigmaR_value, double * sigmaL_value, double * ChiSquare, TH1D * hh, int j, int EndCaps, TString temp, TTree* Tree_Data, double mean, double rms, TCanvas * c2, RooPlot* Erawframe, TF1 * f, string nomDossier, string nomFichier, double RangeMin, double RangeMax, double * JanChi2, double * DegreesOfFreedom, double * pValue, int r9sup, double MinVar, double MaxVar, string variableX, int isMC)
{

	RooRealVar * mmg_s = new RooRealVar("mmg_s", "s", -0.5, 0.5);
        RooRealVar * Photon_SC_Eta = new RooRealVar("Photon_SC_Eta", "Photon_SC_Eta", -3.0, 3.0, "");
        RooRealVar * Photon_SC_Phi = new RooRealVar("Photon_SC_Phi", "Photon_SC_Phi", -3.5, 3.5, "");
	RooRealVar * Photon_r9 = new RooRealVar("Photon_r9", "Photon_r9", 0.0, 1.0, "");
        RooRealVar * isJanLooseMMG = new RooRealVar("isJanLooseMMG", "isJanLooseMMG", 0.0, 1.0, "");
        RooRealVar * isMultipleCandidate = new RooRealVar("isMultipleCandidate", "isMultipleCandidate", 0.0, 1.0, "");
	RooRealVar * Photon_Et = new RooRealVar("Photon_Et", "Photon_Et", 0.0, 250.0, "");	
	RooRealVar * Photon_SC_rawEt = new RooRealVar("Photon_SC_rawEt", "Photon_SC_rawEt", 0.0, 250.0, ""); 
        RooRealVar * Photon_E = new RooRealVar("Photon_E", "Photon_E", 0.0, 1000.0, ""); 
        RooRealVar * Photon_SC_brem = new RooRealVar("Photon_SC_brem", "Photon_SC_brem", 0.0, 15.0, ""); 
        RooRealVar * weight_pileUp = new RooRealVar("weight_pileUp", "weight_pileUp", 0.0, 10000000); 

	RooArgSet *ntplVars = new RooArgSet(*mmg_s, *Photon_SC_Eta, *Photon_SC_Phi, *Photon_r9, *Photon_Et, *Photon_SC_rawEt, *Photon_E, *Photon_SC_brem, *isJanLooseMMG);

        ntplVars->add(*isMultipleCandidate);
        ntplVars->add(*weight_pileUp);



        //RooArgSet *ntplVars = new RooArgSet(*mmg_s, *Photon_SC_Eta, *Photon_SC_Phi, *Photon_r9, *Photon_Et, *Photon_SC_rawEt, *Photon_E, *Photon_SC_brem, *isJanLooseMMG, *isMultipleCandidate, *weight_pileUp);



        //RooArgSet *ntplVars = new RooArgSet(*mmg_s, *Photon_SC_Eta, *Photon_SC_Phi, *Photon_r9, *Photon_Et, *isJanLooseMMG, *isMultipleCandidate, *weight_pileUp);
	//RooArgSet ntplVars(mmg_s, Photon_SC_Eta, Photon_r9, Photon_Et, isJanLooseMMG, Photon_Et);
        //RooDataSet *Data = new RooDataSet("Data", "Data", Tree_Data, *ntplVars); //non weighté
	RooDataSet *Data = new RooDataSet("Data", "Data", Tree_Data, *ntplVars, "", "weight_pileUp"); //weighté
	//RooDataSet Data("Data", "Data", Tree_Data, ntplVars);


        RooDataSet* Data_subset = (RooDataSet*)Data->reduce(*ntplVars, temp);
	//RooDataSet Data_subset = (RooDataSet)Data.reduce(ntplVars, temp);
        //RooPlot* Erawframe = mmg_s.frame();
        
/*	
	Erawframe = mmg_s->frame();
	Data_subset->plotOn(Erawframe,Name("myhist"),Binning(300));
        Erawframe->Draw();
*/

	RooRealVar * Eraw_CB_m0 = new RooRealVar("Eraw_CB_m0", "CB #Delta m_{0}", mean, mean - rms, mean + rms, "GeV");
        RooRealVar * Eraw_CB_sigma = new RooRealVar("Eraw_CB_sigma", "CB ", rms, 0.0, 2 * rms, "GeV");
        RooRealVar * Eraw_CB_alpha = new RooRealVar("Eraw_CB_alpha", "CB #alpha", 0.9, -5.0, 5.0);
        RooRealVar * Eraw_CB_n = new RooRealVar("Eraw_CB_n", "CB n", 10.0, 1.0, 500.0);


/*
	//Ne pas supprimer !!!!//
	RooRealVar Eraw_CB_m0("Eraw_CB_m0", "CB #Delta m_{0}", 1.0, 0.9, 1.05, "GeV");
        RooRealVar Eraw_CB_sigma("Eraw_CB_sigma", "CB ", 0.0215, 0.0, 0.5, "GeV");
        RooRealVar Eraw_CB_alpha("Eraw_CB_alpha", "CB #alpha", 0.9, 0.0, 3.0);
        RooRealVar Eraw_CB_n("Eraw_CB_n", "CB n", 10.0, 0roo.0, 60.0);
*/
        RooCBShape * Eraw_CrystalBall = new RooCBShape("Eraw_CrystalBall","mmg_s_CrystalBall", *mmg_s, *Eraw_CB_m0, *Eraw_CB_sigma, *Eraw_CB_alpha, *Eraw_CB_n);

	int fewBins = 1;
	Double_t Chi2J;

	for(int RightBinning = 600; RightBinning > 10; RightBinning -= 10)
	{
		cout<<endl<<endl<<endl<<"RightBinning = "<<RightBinning<<endl<<endl<<endl;
		Erawframe->Clear();
		Data_subset = (RooDataSet*)Data->reduce(*ntplVars, temp);		
		Eraw_CrystalBall = new RooCBShape("Eraw_CrystalBall","mmg_s_CrystalBall", *mmg_s, *Eraw_CB_m0, *Eraw_CB_sigma, *Eraw_CB_alpha, *Eraw_CB_n); 

		Erawframe = mmg_s->frame();
		Data_subset->plotOn(Erawframe,Name("myhist"),Binning(RightBinning));
		//Erawframe->Draw();
        	//Eraw_CrystalBall.fitTo(*Data_subset, Range(0.6,1.2));
		//Eraw_CrystalBall.fitTo(Data_subset, Range(0.5,1.15)); //ATTENTION NE PAS SUPPRIMER!!!!
		//Eraw_CrystalBall.fitTo(Data_subset, Range(mean - 5 * rms,mean + 5 * rms));
		//double maxDistri = hh->GetMaximumBin() * 0.01;
		//Eraw_CrystalBall.fitTo(*Data_subset, Range(mean - 0.35 * rms,mean + 0.35 * rms));
		//Eraw_CrystalBall.fitTo(*Data_subset, Range(maxDistri - 0.6 * rms,maxDistri + 0.6 * rms));
		Eraw_CrystalBall->fitTo(*Data_subset, Range(RangeMin, RangeMax));
		//RooArgSet* Eraw_CrystalBall_param = Eraw_CrystalBall.getVariables();
		//Eraw_CrystalBall_param->Print("v");
		//Eraw_CrystalBall.plotOn(Erawframe);
        	Eraw_CrystalBall->plotOn(Erawframe,Name("mycurve"));
		Erawframe->Draw();


		Chi2J = chiSquare(Erawframe,(char *)"mycurve",(char *)"myhist",4, JanChi2, DegreesOfFreedom, pValue, &fewBins);
		if(fewBins == 0) break; 


	}


	//RooDataHist data("data", "dataset with mmg_s",*mmg_s,hh);
	//TF1 * f = Eraw_CrystalBall.asTF( RooArgList(mmg_s) );
	f = Eraw_CrystalBall->asTF( RooArgList(*mmg_s) );
	//ftest = Eraw_CrystalBall.asTF( RooArgList(*mmg_s) );
	//cout<<endl<<"ftest->GetXmax() = "<<ftest->GetXmax()<<endl;
	//cout<<endl<<"ftest->GetMaximum = "<<ftest->GetMaximum(0.0, 2.0,1.E-10,100,false)<<endl;
	*mean_value  = Eraw_CB_m0->getVal();
        *mean_error  = Eraw_CB_m0->getError();
        *sigmaEff_value = effSigma(hh);
        *sigma_value = Eraw_CB_sigma->getVal();
	*sigma_value_error = Eraw_CB_sigma->getError();
        *sigmaR_value = SigmaR(f, 0.0, 2.0);
        *sigmaL_value = SigmaL(f, 0.0, 2.0);
        *ChiSquare = Erawframe->chiSquare();

 
        double entries = hh->GetEntries();
        
        int entriesInt = (int) entries;


        TLatex latexLabel;
        latexLabel.SetTextSize(0.026);
        latexLabel.SetNDC();
        latexLabel.DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
	if(isMC == 1) latexLabel.DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV");
        if(isMC == 0) latexLabel.DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV, #int L = 4.93 fb^{-1}");	
        //if(EndCaps == 0) latexLabel.DrawLatex(0.16, 0.85, "Barrel");
        //if(EndCaps == 1) latexLabel.DrawLatex(0.16, 0.85, "EndCaps");   
        string * EndCapsR9Chain(0);
        EndCapsR9Chain = new string;
        if(EndCaps == 0) *EndCapsR9Chain = "Barrel, ";
        if(EndCaps == 1) *EndCapsR9Chain = "Endcaps, ";
        if(r9sup == 0) *EndCapsR9Chain += "Low r9";
        if(r9sup == 1) *EndCapsR9Chain += "High r9";
	if(r9sup == 2) *EndCapsR9Chain += "All r9";
        string * tempLegChain(0);
        tempLegChain = new string;
        *tempLegChain = DoubleToString(MinVar) + " <" + variableX + "< "  + DoubleToString(MaxVar);     
        latexLabel.DrawLatex(0.16, 0.80, EndCapsR9Chain->c_str());
        latexLabel.DrawLatex(0.16, 0.75, tempLegChain->c_str());

        latexLabel.DrawLatex(0.16, 0.70, Form("Entries = %d",entriesInt));
        cout<<endl<<"entries = "<<entries<<endl;
        cout<<endl<<"Mean = "<<hh->GetMean()<<endl;
        latexLabel.DrawLatex(0.16, 0.65, Form("Mean = %f +/- %f",hh->GetMean(), hh->GetMeanError()));
        latexLabel.DrawLatex(0.16, 0.60, Form("RMS = %f +/- %f",hh->GetRMS(), hh->GetRMSError()));
	latexLabel.DrawLatex(0.61, 0.90, Form("#color[4]{m_{0} = %f +/- %f}",Eraw_CB_m0->getVal(), Eraw_CB_m0->getError()));
        latexLabel.DrawLatex(0.61, 0.85, Form("#color[4]{#sigma = %f +/- %f}",Eraw_CB_sigma->getVal(), Eraw_CB_sigma->getError()));
        latexLabel.DrawLatex(0.61, 0.80, Form("#color[4]{#alpha = %f +/- %f}",Eraw_CB_alpha->getVal(), Eraw_CB_alpha->getError()));
        latexLabel.DrawLatex(0.61, 0.75, Form("#color[4]{n = %f +/- %f}",Eraw_CB_n->getVal(), Eraw_CB_n->getError()));
        latexLabel.DrawLatex(0.61, 0.70, Form("#color[4]{Default #chi^{2}/ndf = %f}",Erawframe->chiSquare()));
	latexLabel.DrawLatex(0.61, 0.65, Form("#color[4]{Jan #chi^{2}/ndf = %f}",Chi2J));
	latexLabel.DrawLatex(0.61, 0.60, Form("#color[4]{p-value = %0.6e}",*pValue));

	enregistrementPlots(nomDossier, nomFichier, EndCaps, j, c2);




	c2->Clear();

	////////// chi2 Pulls & Residuals //////////

	string nomDossierPull = "";
        if(EndCaps == 0) nomDossierPull = nomDossier + "EB/";
        if(EndCaps == 1) nomDossierPull = nomDossier + "EE/";

        RooHist* residuals = residHist(Erawframe,(char *)"myhist",(char *)"mycurve", false, nomDossierPull, j);

        //RooHist* pulls = pullHist(Erawframe,"myhist", "mycurve");
        RooHist* pulls = residHist(Erawframe,(char *)"myhist", (char *)"mycurve", true, nomDossierPull, j);

	residuals->Draw("AP");
	nomFichier = "Chi2Residuals";
	
	residuals->GetYaxis()->SetTitle("#chi^{2} Residuals");
        residuals->GetXaxis()->SetLabelFont(42);
        residuals->GetXaxis()->SetTitleFont(42);
        residuals->GetXaxis()->SetLabelSize(0.03);
        residuals->GetXaxis()->SetTitle("s");
        residuals->GetYaxis()->SetLabelFont(42);
        residuals->GetYaxis()->SetTitleOffset(1.24);
        residuals->GetYaxis()->SetTitleFont(42);
        residuals->GetYaxis()->SetLabelSize(0.03);
        //residuals->SetMarkerColor(4);
        //residuals->SetMarkerStyle(21);
        //residuals->SetMarkerSize(0.6);

	TLine *lineResid = new TLine(-0.5,0,0.5,0);
        lineResid->SetLineStyle(3);
        lineResid->Draw();

        TLatex *textResid = new TLatex();

        textResid = new TLatex();
        textResid->SetNDC();
        textResid->SetTextAlign(11);
        textResid->SetTextFont(42);
        textResid->SetTextSizePixels(17);
        textResid->SetTextSize(0.028);
        //if(EndCaps == 0) textResid->DrawLatex(0.80, 0.88, "Barrel");
        //if(EndCaps == 1) textResid->DrawLatex(0.77, 0.88, "End Caps");
	textResid->DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
	if(isMC == 1) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV");
        if(isMC == 0) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV, #int L = 4.93 fb^{-1}");
        textResid->DrawLatex(0.16, 0.80, EndCapsR9Chain->c_str());
        textResid->DrawLatex(0.16, 0.75, tempLegChain->c_str());
	residuals->GetXaxis()->SetLimits(-0.5,0.5);
	
	enregistrementPlots(nomDossier, nomFichier, EndCaps, j, c2);

	c2->Clear();

	pulls->Draw("AP");
	nomFichier = "Chi2Pulls";
	
	pulls->GetYaxis()->SetTitle("#chi^{2} Pulls");
        pulls->GetXaxis()->SetLabelFont(42);
        pulls->GetXaxis()->SetTitleFont(42);
        pulls->GetXaxis()->SetLabelSize(0.03);
        pulls->GetXaxis()->SetTitle("s");
        pulls->GetYaxis()->SetLabelFont(42);
        pulls->GetYaxis()->SetTitleOffset(1.24);
        pulls->GetYaxis()->SetTitleFont(42);
        pulls->GetYaxis()->SetLabelSize(0.03);
        //pulls->SetMarkerColor(4);
        //pulls->SetMarkerStyle(21);
        //pulls->SetMarkerSize(0.6);

	lineResid = new TLine(-0.5,0,0.5,0);
        lineResid->SetLineStyle(3);
        lineResid->Draw();

        textResid = new TLatex();
        textResid->SetNDC();
        textResid->SetTextAlign(11);
        textResid->SetTextFont(42);
        textResid->SetTextSizePixels(17);
        textResid->SetTextSize(0.028);
        //if(EndCaps == 0) textResid->DrawLatex(0.80, 0.88, "Barrel");
        //if(EndCaps == 1) textResid->DrawLatex(0.77, 0.88, "End Caps");	
	textResid->DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
	if(isMC == 1) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV");
        if(isMC == 0) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV, #int L = 4.93 fb^{-1}");
        textResid->DrawLatex(0.16, 0.80, EndCapsR9Chain->c_str());
        textResid->DrawLatex(0.16, 0.75, tempLegChain->c_str());
	pulls->GetXaxis()->SetLimits(-0.5,0.5);	

	enregistrementPlots(nomDossier, nomFichier, EndCaps, j, c2);	

	c2->Clear();

	delete EndCapsR9Chain;
        delete tempLegChain;
	residuals->Delete();
	pulls->Delete();	
	textResid->Delete();
	lineResid->Delete();
	Data_subset->Delete();
	Data->Delete();
	ntplVars->Delete();
	mmg_s->Delete();
        Photon_SC_Eta->Delete();
        Photon_SC_Phi->Delete();
	Photon_r9->Delete();
        isJanLooseMMG->Delete();
        isMultipleCandidate->Delete();
	Photon_Et->Delete();
	Photon_SC_rawEt->Delete();
        Photon_E->Delete();
        Photon_SC_brem->Delete();
        weight_pileUp->Delete();
        Eraw_CB_m0->Delete();
        Eraw_CB_sigma->Delete();
        Eraw_CB_alpha->Delete();
        Eraw_CB_n->Delete();
	Eraw_CrystalBall->Delete();



}


void RooLogNormal(double * mean_value, double * mean_error, Double_t * sigmaEff_value, double * sigma_value, double * sigma_value_error, double * sigmaR_value, double * sigmaL_value, double * ChiSquare, TH1D * hh, int j, int EndCaps, TString temp, TTree* Tree_Data, double mean, double rms, TCanvas * c2, RooPlot* Erawframe, TF1 * f, string nomDossier, string nomFichier, double RangeMin, double RangeMax, double * JanChi2, double * DegreesOfFreedom, double * pValue, int r9sup, double MinVar, double MaxVar, string variableX, int isMC)
{

	RooRealVar * mmg_s = new RooRealVar("mmg_s", "s", -0.5, 0.5);
        RooRealVar * Photon_SC_Eta = new RooRealVar("Photon_SC_Eta", "Photon_SC_Eta", -3.0, 3.0, "");
        RooRealVar * Photon_SC_Phi = new RooRealVar("Photon_SC_Phi", "Photon_SC_Phi", -3.5, 3.5, "");
        RooRealVar * Photon_r9 = new RooRealVar("Photon_r9", "Photon_r9", 0.0, 1.0, "");
        RooRealVar * isJanLooseMMG = new RooRealVar("isJanLooseMMG", "isJanLooseMMG", 0.0, 1.0, "");
        RooRealVar * isMultipleCandidate = new RooRealVar("isMultipleCandidate", "isMultipleCandidate", 0.0, 1.0, "");
        RooRealVar * Photon_Et = new RooRealVar("Photon_Et", "Photon_Et", 0.0, 250.0, "");
        RooRealVar * Photon_SC_rawEt = new RooRealVar("Photon_SC_rawEt", "Photon_SC_rawEt", 0.0, 250.0, "");
        RooRealVar * Photon_E = new RooRealVar("Photon_E", "Photon_E", 0.0, 1000.0, "");
        RooRealVar * Photon_SC_brem = new RooRealVar("Photon_SC_brem", "Photon_SC_brem", 0.0, 15.0, "");
        RooRealVar * weight_pileUp = new RooRealVar("weight_pileUp", "weight_pileUp", 0.0, 10000000);

        RooArgSet *ntplVars = new RooArgSet(*mmg_s, *Photon_SC_Eta, *Photon_SC_Phi, *Photon_r9, *Photon_Et, *Photon_SC_rawEt, *Photon_E, *Photon_SC_brem, *isJanLooseMMG);

        ntplVars->add(*isMultipleCandidate);
        ntplVars->add(*weight_pileUp);

	RooDataSet *Data = new RooDataSet("Data", "Data", Tree_Data, *ntplVars, "", "weight_pileUp");
	RooDataSet* Data_subset = (RooDataSet*)Data->reduce(*ntplVars, temp);


	RooRealVar * Eraw_LN_m0 = new RooRealVar("mmg_s_LN_m0", "CB #Delta m_{0}", 1.0, -5.0, 5.0, "GeV");
        RooRealVar * Eraw_LN_k = new RooRealVar("mmg_s_LN_sigma", "CB ", 0.45, 0.01, 2.0, "GeV");

        RooLognormal * Eraw_LogNormal = new RooLognormal("Eraw_LogNormal", "mmg_s_LogNormal", *mmg_s, *Eraw_LN_m0, *Eraw_LN_k);


	int fewBins = 1;
        Double_t Chi2J;

	for(int RightBinning = 600; RightBinning > 10; RightBinning -= 10)
        {
                cout<<endl<<endl<<endl<<"RightBinning = "<<RightBinning<<endl<<endl<<endl;
                Erawframe->Clear();
                Data_subset = (RooDataSet*)Data->reduce(*ntplVars, temp);
		Eraw_LogNormal = new RooLognormal("Eraw_LogNormal", "mmg_s_LogNormal", *mmg_s, *Eraw_LN_m0, *Eraw_LN_k);

                Erawframe = mmg_s->frame();
                Data_subset->plotOn(Erawframe,Name("myhist"),Binning(RightBinning));
                Eraw_LogNormal->fitTo(*Data_subset, Range(RangeMin, RangeMax));
                Eraw_LogNormal->plotOn(Erawframe,Name("mycurve"));
                Erawframe->Draw();

                Chi2J = chiSquare(Erawframe,(char *)"mycurve",(char *)"myhist",2, JanChi2, DegreesOfFreedom, pValue, &fewBins);
                if(fewBins == 0) break;
        }

	//RooDataHist data("data", "dataset with mmg_s",mmg_s,hh);
        //TF1 * f = Eraw_LogNormal.asTF( RooArgList(mmg_s) );
        f = Eraw_LogNormal->asTF( RooArgList(*mmg_s) );
        *mean_value  = Eraw_LN_m0->getVal();
        *mean_error  = Eraw_LN_m0->getError();
        *sigmaEff_value = effSigma(hh);
        *sigma_value = Eraw_LN_k->getVal();
        *sigma_value_error = Eraw_LN_k->getError();
        *sigmaR_value = SigmaR(f, 0.0, 2.0);
        *sigmaL_value = SigmaL(f, 0.0, 2.0);
        *ChiSquare = Erawframe->chiSquare();

        double entries = hh->GetEntries();

        int entriesInt = (int) entries;

	TLatex latexLabel;
        latexLabel.SetTextSize(0.026);
        latexLabel.SetNDC();
        latexLabel.DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
        if(isMC == 1) latexLabel.DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV");
        if(isMC == 0) latexLabel.DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV, #int L = 4.93 fb^{-1}");
        //if(EndCaps == 0) latexLabel.DrawLatex(0.16, 0.85, "Barrel");
        //if(EndCaps == 1) latexLabel.DrawLatex(0.16, 0.85, "EndCaps");   
        string * EndCapsR9Chain(0);
        EndCapsR9Chain = new string;
        if(EndCaps == 0) *EndCapsR9Chain = "Barrel, ";
        if(EndCaps == 1) *EndCapsR9Chain = "Endcaps, ";
        if(r9sup == 0) *EndCapsR9Chain += "Low r9";
        if(r9sup == 1) *EndCapsR9Chain += "High r9";
        if(r9sup == 2) *EndCapsR9Chain += "All r9";
	string * tempLegChain(0);
        tempLegChain = new string;
        *tempLegChain = DoubleToString(MinVar) + " <" + variableX + "< "  + DoubleToString(MaxVar);
        latexLabel.DrawLatex(0.16, 0.80, EndCapsR9Chain->c_str());
        latexLabel.DrawLatex(0.16, 0.75, tempLegChain->c_str());

        latexLabel.DrawLatex(0.16, 0.70, Form("Entries = %d",entriesInt));
        cout<<endl<<"entries = "<<entries<<endl;
        cout<<endl<<"Mean = "<<hh->GetMean()<<endl;
        latexLabel.DrawLatex(0.16, 0.65, Form("Mean = %f +/- %f",hh->GetMean(), hh->GetMeanError()));
        latexLabel.DrawLatex(0.16, 0.60, Form("RMS = %f +/- %f",hh->GetRMS(), hh->GetRMSError()));
	latexLabel.DrawLatex(0.61, 0.90, Form("#color[4]{m_{0} = %f +/- %f}",Eraw_LN_m0->getVal(), Eraw_LN_m0->getError()));
        latexLabel.DrawLatex(0.61, 0.85, Form("#color[4]{k = %f +/- %f}",Eraw_LN_k->getVal(), Eraw_LN_k->getError()));
	latexLabel.DrawLatex(0.61, 0.80, Form("#color[4]{Default #chi^{2}/ndf = %f}",Erawframe->chiSquare()));
        latexLabel.DrawLatex(0.61, 0.75, Form("#color[4]{Jan #chi^{2}/ndf = %f}",Chi2J));
        latexLabel.DrawLatex(0.61, 0.70, Form("#color[4]{p-value = %0.6e}",*pValue));

        enregistrementPlots(nomDossier, nomFichier, EndCaps, j, c2);

	c2->Clear();

	////////// chi2 Pulls & Residuals //////////

        string nomDossierPull = "";
        if(EndCaps == 0) nomDossierPull = nomDossier + "EB/";
        if(EndCaps == 1) nomDossierPull = nomDossier + "EE/";

        RooHist* residuals = residHist(Erawframe,(char *)"myhist",(char *)"mycurve", false, nomDossierPull, j);

        //RooHist* pulls = pullHist(Erawframe,"myhist", "mycurve");
        RooHist* pulls = residHist(Erawframe,(char *)"myhist", (char *)"mycurve", true, nomDossierPull, j);

        residuals->Draw("AP");
        nomFichier = "Chi2Residuals";
      
        residuals->GetYaxis()->SetTitle("#chi^{2} Residuals");
        residuals->GetXaxis()->SetLabelFont(42);
        residuals->GetXaxis()->SetTitleFont(42);
        residuals->GetXaxis()->SetLabelSize(0.03);
        residuals->GetXaxis()->SetTitle("s");
        residuals->GetYaxis()->SetLabelFont(42);
        residuals->GetYaxis()->SetTitleOffset(1.24);
        residuals->GetYaxis()->SetTitleFont(42);
        residuals->GetYaxis()->SetLabelSize(0.03);
        //residuals->SetMarkerColor(4);
        //residuals->SetMarkerStyle(21);
        //residuals->SetMarkerSize(0.6);

        TLine *lineResid = new TLine(-0.5,0,0.5,0);
        lineResid->SetLineStyle(3); 
        lineResid->Draw();

        TLatex *textResid = new TLatex();
	textResid = new TLatex();
        textResid->SetNDC();
        textResid->SetTextAlign(11);
        textResid->SetTextFont(42);
        textResid->SetTextSizePixels(17);
        textResid->SetTextSize(0.028);
        //if(EndCaps == 0) textResid->DrawLatex(0.80, 0.88, "Barrel");
        //if(EndCaps == 1) textResid->DrawLatex(0.77, 0.88, "End Caps");
        textResid->DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
        if(isMC == 1) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV");
        if(isMC == 0) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV, #int L = 4.93 fb^{-1}");
        textResid->DrawLatex(0.16, 0.80, EndCapsR9Chain->c_str());
        textResid->DrawLatex(0.16, 0.75, tempLegChain->c_str());
        residuals->GetXaxis()->SetLimits(-0.5,0.5);

        enregistrementPlots(nomDossier, nomFichier, EndCaps, j, c2);

        c2->Clear();

        pulls->Draw("AP");
        nomFichier = "Chi2Pulls";

	pulls->GetYaxis()->SetTitle("#chi^{2} Pulls");
        pulls->GetXaxis()->SetLabelFont(42);
        pulls->GetXaxis()->SetTitleFont(42);
        pulls->GetXaxis()->SetLabelSize(0.03);
        pulls->GetXaxis()->SetTitle("s");
        pulls->GetYaxis()->SetLabelFont(42);
        pulls->GetYaxis()->SetTitleOffset(1.24);
        pulls->GetYaxis()->SetTitleFont(42);
        pulls->GetYaxis()->SetLabelSize(0.03);
        //pulls->SetMarkerColor(4);
        //pulls->SetMarkerStyle(21);
        //pulls->SetMarkerSize(0.6);

        lineResid = new TLine(-0.5,0,0.5,0);
        lineResid->SetLineStyle(3);
        lineResid->Draw();

        textResid = new TLatex();
        textResid->SetNDC();
        textResid->SetTextAlign(11);
        textResid->SetTextFont(42);
        textResid->SetTextSizePixels(17);
        textResid->SetTextSize(0.028);
        //if(EndCaps == 0) textResid->DrawLatex(0.80, 0.88, "Barrel");
        //if(EndCaps == 1) textResid->DrawLatex(0.77, 0.88, "End Caps");        
        textResid->DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
        if(isMC == 1) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV");
        if(isMC == 0) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV, #int L = 4.93 fb^{-1}");
        textResid->DrawLatex(0.16, 0.80, EndCapsR9Chain->c_str());
        textResid->DrawLatex(0.16, 0.75, tempLegChain->c_str());
        pulls->GetXaxis()->SetLimits(-0.5,0.5);

        enregistrementPlots(nomDossier, nomFichier, EndCaps, j, c2);

        c2->Clear();
	
	delete EndCapsR9Chain;
        delete tempLegChain;
        residuals->Delete();
        pulls->Delete();
        textResid->Delete();
        lineResid->Delete();
        Data_subset->Delete();
        Data->Delete();
        ntplVars->Delete();
        mmg_s->Delete();
        Photon_SC_Eta->Delete();
        Photon_SC_Phi->Delete();
        Photon_r9->Delete();
        isJanLooseMMG->Delete();
        isMultipleCandidate->Delete();
        Photon_Et->Delete();
        Photon_SC_rawEt->Delete();
        Photon_E->Delete();
        Photon_SC_brem->Delete();
        weight_pileUp->Delete();
        Eraw_LN_m0->Delete();
        Eraw_LN_k->Delete();
        Eraw_LogNormal->Delete();

}


void RooGauss(double * mean_value, double * mean_error, Double_t * sigmaEff_value, double * sigma_value, double * sigma_value_error, double * sigmaR_value, double * sigmaL_value, double * ChiSquare, TH1D * hh, int j, int EndCaps, TString temp, TTree* Tree_Data, double mean, double rms, TCanvas * c2, RooPlot* Erawframe, TF1 * f, string nomDossier, string nomFichier, double RangeMin, double RangeMax, double * JanChi2, double * DegreesOfFreedom, double * pValue, int r9sup, double MinVar, double MaxVar, string variableX, int isMC)
{

	RooRealVar * mmg_s = new RooRealVar("mmg_s", "s", -0.5, 0.5);
        RooRealVar * Photon_SC_Eta = new RooRealVar("Photon_SC_Eta", "Photon_SC_Eta", -3.0, 3.0, "");
        RooRealVar * Photon_SC_Phi = new RooRealVar("Photon_SC_Phi", "Photon_SC_Phi", -3.5, 3.5, "");
        RooRealVar * Photon_r9 = new RooRealVar("Photon_r9", "Photon_r9", 0.0, 1.0, "");
        RooRealVar * isJanLooseMMG = new RooRealVar("isJanLooseMMG", "isJanLooseMMG", 0.0, 1.0, "");
        RooRealVar * isMultipleCandidate = new RooRealVar("isMultipleCandidate", "isMultipleCandidate", 0.0, 1.0, "");
        RooRealVar * Photon_Et = new RooRealVar("Photon_Et", "Photon_Et", 0.0, 250.0, "");
        RooRealVar * Photon_SC_rawEt = new RooRealVar("Photon_SC_rawEt", "Photon_SC_rawEt", 0.0, 250.0, "");
        RooRealVar * Photon_E = new RooRealVar("Photon_E", "Photon_E", 0.0, 1000.0, "");
        RooRealVar * Photon_SC_brem = new RooRealVar("Photon_SC_brem", "Photon_SC_brem", 0.0, 15.0, "");
        RooRealVar * weight_pileUp = new RooRealVar("weight_pileUp", "weight_pileUp", 0.0, 10000000);

        RooArgSet *ntplVars = new RooArgSet(*mmg_s, *Photon_SC_Eta, *Photon_SC_Phi, *Photon_r9, *Photon_Et, *Photon_SC_rawEt, *Photon_E, *Photon_SC_brem, *isJanLooseMMG);

        ntplVars->add(*isMultipleCandidate);
        ntplVars->add(*weight_pileUp);

	RooDataSet *Data = new RooDataSet("Data", "Data", Tree_Data, *ntplVars, "", "weight_pileUp");
	RooDataSet* Data_subset = (RooDataSet*)Data->reduce(*ntplVars, temp);


	RooRealVar * sigmean = new RooRealVar("sigmean","B^{#pm} mass",1.0,0.7,1.2);
        RooRealVar * sigwidth = new RooRealVar("sigwidth","B^{#pm} width",0.5,0.0,1.0) ;

        // --- Build Gaussian PDF ---
        RooGaussian * Eraw_Gaussian = new RooGaussian("Eraw_Gaussian","mmg_s_Gaussian",*mmg_s,*sigmean,*sigwidth) ;



	int fewBins = 1;
        Double_t Chi2J;

	for(int RightBinning = 600; RightBinning > 10; RightBinning -= 10)
        {
                cout<<endl<<endl<<endl<<"RightBinning = "<<RightBinning<<endl<<endl<<endl;
                Erawframe->Clear();
                Data_subset = (RooDataSet*)Data->reduce(*ntplVars, temp);
		Eraw_Gaussian = new RooGaussian("Eraw_Gaussian","mmg_s_Gaussian",*mmg_s,*sigmean,*sigwidth);

                Erawframe = mmg_s->frame();
                Data_subset->plotOn(Erawframe,Name("myhist"),Binning(RightBinning));
                Eraw_Gaussian->fitTo(*Data_subset, Range(RangeMin, RangeMax));
                Eraw_Gaussian->plotOn(Erawframe,Name("mycurve"));
                Erawframe->Draw();

                Chi2J = chiSquare(Erawframe,(char *)"mycurve",(char *)"myhist",2, JanChi2, DegreesOfFreedom, pValue, &fewBins);
                if(fewBins == 0) break;
        }

	//RooDataHist data("data", "dataset with mmg_s",mmg_s,hh);
        //TF1 * f = Eraw_Gaussian.asTF( RooArgList(mmg_s) );
        f = Eraw_Gaussian->asTF( RooArgList(*mmg_s) );
        *mean_value  = sigmean->getVal();
        *mean_error  = sigmean->getError();
        *sigmaEff_value = effSigma(hh);
        *sigma_value = sigwidth->getVal();
        *sigma_value_error = sigwidth->getError();
        *sigmaR_value = SigmaR(f, 0.0, 2.0);
        *sigmaL_value = SigmaL(f, 0.0, 2.0);
        *ChiSquare = Erawframe->chiSquare();
	

        double entries = hh->GetEntries();

        int entriesInt = (int) entries;

	TLatex latexLabel;
        latexLabel.SetTextSize(0.026);
        latexLabel.SetNDC();
        latexLabel.DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
        if(isMC == 1) latexLabel.DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV");
        if(isMC == 0) latexLabel.DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV, #int L = 4.93 fb^{-1}");
        //if(EndCaps == 0) latexLabel.DrawLatex(0.16, 0.85, "Barrel");
        //if(EndCaps == 1) latexLabel.DrawLatex(0.16, 0.85, "EndCaps");   
        string * EndCapsR9Chain(0);
        EndCapsR9Chain = new string;
        if(EndCaps == 0) *EndCapsR9Chain = "Barrel, ";
        if(EndCaps == 1) *EndCapsR9Chain = "Endcaps, ";
        if(r9sup == 0) *EndCapsR9Chain += "Low r9";
        if(r9sup == 1) *EndCapsR9Chain += "High r9";
        if(r9sup == 2) *EndCapsR9Chain += "All r9";
	string * tempLegChain(0);
        tempLegChain = new string;
        *tempLegChain = DoubleToString(MinVar) + " <" + variableX + "< "  + DoubleToString(MaxVar);
        latexLabel.DrawLatex(0.16, 0.80, EndCapsR9Chain->c_str());
        latexLabel.DrawLatex(0.16, 0.75, tempLegChain->c_str());

        latexLabel.DrawLatex(0.16, 0.70, Form("Entries = %d",entriesInt));
        cout<<endl<<"entries = "<<entries<<endl;
        cout<<endl<<"Mean = "<<hh->GetMean()<<endl;
        latexLabel.DrawLatex(0.16, 0.65, Form("Mean = %f +/- %f",hh->GetMean(), hh->GetMeanError()));
        latexLabel.DrawLatex(0.16, 0.60, Form("RMS = %f +/- %f",hh->GetRMS(), hh->GetRMSError()));
	latexLabel.DrawLatex(0.61, 0.90, Form("#color[4]{m_{0} = %f +/- %f}",sigmean->getVal(), sigmean->getError()));
        latexLabel.DrawLatex(0.61, 0.85, Form("#color[4]{#sigma = %f +/- %f}",sigwidth->getVal(), sigwidth->getError()));
	latexLabel.DrawLatex(0.61, 0.80, Form("#color[4]{Default #chi^{2}/ndf = %f}",Erawframe->chiSquare()));
        latexLabel.DrawLatex(0.61, 0.75, Form("#color[4]{Jan #chi^{2}/ndf = %f}",Chi2J));
        latexLabel.DrawLatex(0.61, 0.70, Form("#color[4]{p-value = %0.6e}",*pValue));

        enregistrementPlots(nomDossier, nomFichier, EndCaps, j, c2);

	c2->Clear();

	////////// chi2 Pulls & Residuals //////////

        string nomDossierPull = "";
        if(EndCaps == 0) nomDossierPull = nomDossier + "EB/";
        if(EndCaps == 1) nomDossierPull = nomDossier + "EE/";

        RooHist* residuals = residHist(Erawframe,(char *)"myhist",(char *)"mycurve", false, nomDossierPull, j);

        //RooHist* pulls = pullHist(Erawframe,"myhist", "mycurve");
        RooHist* pulls = residHist(Erawframe,(char *)"myhist", (char *)"mycurve", true, nomDossierPull, j);

        residuals->Draw("AP");
        nomFichier = "Chi2Residuals";
      
        residuals->GetYaxis()->SetTitle("#chi^{2} Residuals");
        residuals->GetXaxis()->SetLabelFont(42);
        residuals->GetXaxis()->SetTitleFont(42);
        residuals->GetXaxis()->SetLabelSize(0.03);
        residuals->GetXaxis()->SetTitle("s");
        residuals->GetYaxis()->SetLabelFont(42);
        residuals->GetYaxis()->SetTitleOffset(1.24);
        residuals->GetYaxis()->SetTitleFont(42);
        residuals->GetYaxis()->SetLabelSize(0.03);
        //residuals->SetMarkerColor(4);
        //residuals->SetMarkerStyle(21);
        //residuals->SetMarkerSize(0.6);

        TLine *lineResid = new TLine(-0.5,0,0.5,0);
        lineResid->SetLineStyle(3); 
        lineResid->Draw();

        TLatex *textResid = new TLatex();
	textResid = new TLatex();
        textResid->SetNDC();
        textResid->SetTextAlign(11);
        textResid->SetTextFont(42);
        textResid->SetTextSizePixels(17);
        textResid->SetTextSize(0.028);
        //if(EndCaps == 0) textResid->DrawLatex(0.80, 0.88, "Barrel");
        //if(EndCaps == 1) textResid->DrawLatex(0.77, 0.88, "End Caps");
        textResid->DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
        if(isMC == 1) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV");
        if(isMC == 0) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV, #int L = 4.93 fb^{-1}");
        textResid->DrawLatex(0.16, 0.80, EndCapsR9Chain->c_str());
        textResid->DrawLatex(0.16, 0.75, tempLegChain->c_str());
        residuals->GetXaxis()->SetLimits(-0.5,0.5);

        enregistrementPlots(nomDossier, nomFichier, EndCaps, j, c2);

        c2->Clear();

        pulls->Draw("AP");
        nomFichier = "Chi2Pulls";

	pulls->GetYaxis()->SetTitle("#chi^{2} Pulls");
        pulls->GetXaxis()->SetLabelFont(42);
        pulls->GetXaxis()->SetTitleFont(42);
        pulls->GetXaxis()->SetLabelSize(0.03);
        pulls->GetXaxis()->SetTitle("s");
        pulls->GetYaxis()->SetLabelFont(42);
        pulls->GetYaxis()->SetTitleOffset(1.24);
        pulls->GetYaxis()->SetTitleFont(42);
        pulls->GetYaxis()->SetLabelSize(0.03);
        //pulls->SetMarkerColor(4);
        //pulls->SetMarkerStyle(21);
        //pulls->SetMarkerSize(0.6);

        lineResid = new TLine(-0.5,0,0.5,0);
        lineResid->SetLineStyle(3);
        lineResid->Draw();

        textResid = new TLatex();
        textResid->SetNDC();
        textResid->SetTextAlign(11);
        textResid->SetTextFont(42);
        textResid->SetTextSizePixels(17);
        textResid->SetTextSize(0.028);
        //if(EndCaps == 0) textResid->DrawLatex(0.80, 0.88, "Barrel");
        //if(EndCaps == 1) textResid->DrawLatex(0.77, 0.88, "End Caps");        
        textResid->DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
        if(isMC == 1) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV");
        if(isMC == 0) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV, #int L = 4.93 fb^{-1}");
        textResid->DrawLatex(0.16, 0.80, EndCapsR9Chain->c_str());
        textResid->DrawLatex(0.16, 0.75, tempLegChain->c_str());
        pulls->GetXaxis()->SetLimits(-0.5,0.5);

        enregistrementPlots(nomDossier, nomFichier, EndCaps, j, c2);

        c2->Clear();
	
	delete EndCapsR9Chain;
        delete tempLegChain;
        residuals->Delete();
        pulls->Delete();
        textResid->Delete();
        lineResid->Delete();
        Data_subset->Delete();
        Data->Delete();
        ntplVars->Delete();
        mmg_s->Delete();
        Photon_SC_Eta->Delete();
        Photon_SC_Phi->Delete();
        Photon_r9->Delete();
        isJanLooseMMG->Delete();
        isMultipleCandidate->Delete();
        Photon_Et->Delete();
        Photon_SC_rawEt->Delete();
        Photon_E->Delete();
        Photon_SC_brem->Delete();
        weight_pileUp->Delete();
        sigmean->Delete();
        sigwidth->Delete();
        Eraw_Gaussian->Delete();

}


void RooGamma2(double * mean_value, double * mean_error, Double_t * sigmaEff_value, double * sigma_value, double * sigma_value_error, double * sigmaR_value, double * sigmaL_value, double * ChiSquare, TH1D * hh, int j, int EndCaps, TString temp, TTree* Tree_Data, double mean, double rms, TCanvas * c2, RooPlot* Erawframe, TF1 * f, string nomDossier, string nomFichier, double RangeMin, double RangeMax, double * JanChi2, double * DegreesOfFreedom, double * pValue, int r9sup, double MinVar, double MaxVar, string variableX, int isMC)
{

	RooRealVar * mmg_s = new RooRealVar("mmg_s", "s", -0.5, 0.5);
        RooRealVar * Photon_SC_Eta = new RooRealVar("Photon_SC_Eta", "Photon_SC_Eta", -3.0, 3.0, "");
        RooRealVar * Photon_SC_Phi = new RooRealVar("Photon_SC_Phi", "Photon_SC_Phi", -3.5, 3.5, "");
        RooRealVar * Photon_r9 = new RooRealVar("Photon_r9", "Photon_r9", 0.0, 1.0, "");
        RooRealVar * isJanLooseMMG = new RooRealVar("isJanLooseMMG", "isJanLooseMMG", 0.0, 1.0, "");
        RooRealVar * isMultipleCandidate = new RooRealVar("isMultipleCandidate", "isMultipleCandidate", 0.0, 1.0, "");
        RooRealVar * Photon_Et = new RooRealVar("Photon_Et", "Photon_Et", 0.0, 250.0, "");
        RooRealVar * Photon_SC_rawEt = new RooRealVar("Photon_SC_rawEt", "Photon_SC_rawEt", 0.0, 250.0, "");
        RooRealVar * Photon_E = new RooRealVar("Photon_E", "Photon_E", 0.0, 1000.0, "");
        RooRealVar * Photon_SC_brem = new RooRealVar("Photon_SC_brem", "Photon_SC_brem", 0.0, 15.0, "");
        RooRealVar * weight_pileUp = new RooRealVar("weight_pileUp", "weight_pileUp", 0.0, 10000000);

        RooArgSet *ntplVars = new RooArgSet(*mmg_s, *Photon_SC_Eta, *Photon_SC_Phi, *Photon_r9, *Photon_Et, *Photon_SC_rawEt, *Photon_E, *Photon_SC_brem, *isJanLooseMMG);

        ntplVars->add(*isMultipleCandidate);
        ntplVars->add(*weight_pileUp);

	RooDataSet *Data = new RooDataSet("Data", "Data", Tree_Data, *ntplVars, "", "weight_pileUp");
	RooDataSet* Data_subset = (RooDataSet*)Data->reduce(*ntplVars, temp);



 
	// --- Parameters ---
/*	RooRealVar beta("beta","beta",20.0,0.0,200.0);
	RooRealVar gamma("gamma","gamma",1.0,0.0,5.0) ;
	RooRealVar mu("mu","mu",0.5,0.0,1.0) ; 
*/
	RooRealVar * beta = new RooRealVar("beta","beta",20.0,0.0,200.0);
        RooRealVar * gamma = new RooRealVar("gamma","gamma",1.0,0.0,5.0) ;
        RooRealVar * mu = new RooRealVar("mu","mu",0.5,0.0,1.0) ;


	// --- Build Gaussian PDF ---
	RooGamma * Eraw_Gamma = new RooGamma("Eraw_Gamma","mmg_s_Gamma",*mmg_s,*beta,*gamma,*mu);

	int fewBins = 1;
        Double_t Chi2J;

	for(int RightBinning = 600; RightBinning > 10; RightBinning -= 10)
        {
                cout<<endl<<endl<<endl<<"RightBinning = "<<RightBinning<<endl<<endl<<endl;
                Erawframe->Clear();
                Data_subset = (RooDataSet*)Data->reduce(*ntplVars, temp);
		Eraw_Gamma = new RooGamma("Eraw_Gamma","mmg_s_Gamma",*mmg_s,*beta,*gamma,*mu);

                Erawframe = mmg_s->frame();
                Data_subset->plotOn(Erawframe,Name("myhist"),Binning(RightBinning));
                Eraw_Gamma->fitTo(*Data_subset, Range(RangeMin, RangeMax));
                Eraw_Gamma->plotOn(Erawframe,Name("mycurve"));
                Erawframe->Draw();

                Chi2J = chiSquare(Erawframe,(char *)"mycurve",(char *)"myhist",3, JanChi2, DegreesOfFreedom, pValue, &fewBins);
                if(fewBins == 0) break;
        }



        //RooDataHist data("data", "dataset with mmg_s",mmg_s,hh);
        //TF1 * f = Eraw_Gamma.asTF( RooArgList(mmg_s) );
       	f = Eraw_Gamma->asTF( RooArgList(*mmg_s) ); 
	//*mean_value  = sigmean.getVal();
        //*mean_error  = sigmean.getError();
        *sigmaEff_value = effSigma(hh);
        //*sigma_value = sigwidth.getVal();
        *sigmaR_value = SigmaR(f, 0.0, 2.0);
        *sigmaL_value = SigmaL(f, 0.0, 2.0);
        *ChiSquare = Erawframe->chiSquare();


        double entries = hh->GetEntries();

        int entriesInt = (int) entries;

	TLatex latexLabel;
        latexLabel.SetTextSize(0.026);
        latexLabel.SetNDC();
        latexLabel.DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
        if(isMC == 1) latexLabel.DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV");
        if(isMC == 0) latexLabel.DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV, #int L = 4.93 fb^{-1}");
        //if(EndCaps == 0) latexLabel.DrawLatex(0.16, 0.85, "Barrel");
        //if(EndCaps == 1) latexLabel.DrawLatex(0.16, 0.85, "EndCaps");   
        string * EndCapsR9Chain(0);
        EndCapsR9Chain = new string;
        if(EndCaps == 0) *EndCapsR9Chain = "Barrel, ";
        if(EndCaps == 1) *EndCapsR9Chain = "Endcaps, ";
        if(r9sup == 0) *EndCapsR9Chain += "Low r9";
        if(r9sup == 1) *EndCapsR9Chain += "High r9";
        if(r9sup == 2) *EndCapsR9Chain += "All r9";
	string * tempLegChain(0);
        tempLegChain = new string;
        *tempLegChain = DoubleToString(MinVar) + " <" + variableX + "< "  + DoubleToString(MaxVar);
        latexLabel.DrawLatex(0.16, 0.80, EndCapsR9Chain->c_str());
        latexLabel.DrawLatex(0.16, 0.75, tempLegChain->c_str());

        latexLabel.DrawLatex(0.16, 0.70, Form("Entries = %d",entriesInt));
        cout<<endl<<"entries = "<<entries<<endl;
        cout<<endl<<"Mean = "<<hh->GetMean()<<endl;
        latexLabel.DrawLatex(0.16, 0.65, Form("Mean = %f +/- %f",hh->GetMean(), hh->GetMeanError()));
        latexLabel.DrawLatex(0.16, 0.60, Form("RMS = %f +/- %f",hh->GetRMS(), hh->GetRMSError()));
        latexLabel.DrawLatex(0.61, 0.90, Form("#color[4]{#beta = %f +/- %f}",beta->getVal(), beta->getError()));
        latexLabel.DrawLatex(0.61, 0.85, Form("#color[4]{#gamma = %f +/- %f}",gamma->getVal(), gamma->getError()));
        latexLabel.DrawLatex(0.61, 0.80, Form("#color[4]{#mu = %f +/- %f}",mu->getVal(), mu->getError()));
	latexLabel.DrawLatex(0.61, 0.75, Form("#color[4]{Default #chi^{2}/ndf = %f}",Erawframe->chiSquare()));
        latexLabel.DrawLatex(0.61, 0.70, Form("#color[4]{Jan #chi^{2}/ndf = %f}",Chi2J));
        latexLabel.DrawLatex(0.61, 0.65, Form("#color[4]{p-value = %0.6e}",*pValue));

        enregistrementPlots(nomDossier, nomFichier, EndCaps, j, c2);

	c2->Clear();

	////////// chi2 Pulls & Residuals //////////

        string nomDossierPull = "";
        if(EndCaps == 0) nomDossierPull = nomDossier + "EB/";
        if(EndCaps == 1) nomDossierPull = nomDossier + "EE/";

        RooHist* residuals = residHist(Erawframe,(char *)"myhist",(char *)"mycurve", false, nomDossierPull, j);

        //RooHist* pulls = pullHist(Erawframe,"myhist", "mycurve");
        RooHist* pulls = residHist(Erawframe,(char *)"myhist", (char *)"mycurve", true, nomDossierPull, j);

        residuals->Draw("AP");
        nomFichier = "Chi2Residuals";
      
        residuals->GetYaxis()->SetTitle("#chi^{2} Residuals");
        residuals->GetXaxis()->SetLabelFont(42);
        residuals->GetXaxis()->SetTitleFont(42);
        residuals->GetXaxis()->SetLabelSize(0.03);
        residuals->GetXaxis()->SetTitle("s");
        residuals->GetYaxis()->SetLabelFont(42);
        residuals->GetYaxis()->SetTitleOffset(1.24);
        residuals->GetYaxis()->SetTitleFont(42);
        residuals->GetYaxis()->SetLabelSize(0.03);
        //residuals->SetMarkerColor(4);
        //residuals->SetMarkerStyle(21);
        //residuals->SetMarkerSize(0.6);

        TLine *lineResid = new TLine(-0.5,0,0.5,0);
        lineResid->SetLineStyle(3); 
        lineResid->Draw();

        TLatex *textResid = new TLatex();
	textResid = new TLatex();
        textResid->SetNDC();
        textResid->SetTextAlign(11);
        textResid->SetTextFont(42);
        textResid->SetTextSizePixels(17);
        textResid->SetTextSize(0.028);
        //if(EndCaps == 0) textResid->DrawLatex(0.80, 0.88, "Barrel");
        //if(EndCaps == 1) textResid->DrawLatex(0.77, 0.88, "End Caps");
        textResid->DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
        if(isMC == 1) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV");
        if(isMC == 0) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV, #int L = 4.93 fb^{-1}");
        textResid->DrawLatex(0.16, 0.80, EndCapsR9Chain->c_str());
        textResid->DrawLatex(0.16, 0.75, tempLegChain->c_str());
        residuals->GetXaxis()->SetLimits(-0.5,0.5);

        enregistrementPlots(nomDossier, nomFichier, EndCaps, j, c2);

        c2->Clear();

        pulls->Draw("AP");
        nomFichier = "Chi2Pulls";

	pulls->GetYaxis()->SetTitle("#chi^{2} Pulls");
        pulls->GetXaxis()->SetLabelFont(42);
        pulls->GetXaxis()->SetTitleFont(42);
        pulls->GetXaxis()->SetLabelSize(0.03);
        pulls->GetXaxis()->SetTitle("s");
        pulls->GetYaxis()->SetLabelFont(42);
        pulls->GetYaxis()->SetTitleOffset(1.24);
        pulls->GetYaxis()->SetTitleFont(42);
        pulls->GetYaxis()->SetLabelSize(0.03);
        //pulls->SetMarkerColor(4);
        //pulls->SetMarkerStyle(21);
        //pulls->SetMarkerSize(0.6);

        lineResid = new TLine(-0.5,0,0.5,0);
        lineResid->SetLineStyle(3);
        lineResid->Draw();

        textResid = new TLatex();
        textResid->SetNDC();
        textResid->SetTextAlign(11);
        textResid->SetTextFont(42);
        textResid->SetTextSizePixels(17);
        textResid->SetTextSize(0.028);
        //if(EndCaps == 0) textResid->DrawLatex(0.80, 0.88, "Barrel");
        //if(EndCaps == 1) textResid->DrawLatex(0.77, 0.88, "End Caps");        
        textResid->DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
        if(isMC == 1) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV");
        if(isMC == 0) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV, #int L = 4.93 fb^{-1}");
        textResid->DrawLatex(0.16, 0.80, EndCapsR9Chain->c_str());
        textResid->DrawLatex(0.16, 0.75, tempLegChain->c_str());
        pulls->GetXaxis()->SetLimits(-0.5,0.5);

        enregistrementPlots(nomDossier, nomFichier, EndCaps, j, c2);

        c2->Clear();
	
	delete EndCapsR9Chain;
        delete tempLegChain;
        residuals->Delete();
        pulls->Delete();
        textResid->Delete();
        lineResid->Delete();
        Data_subset->Delete();
        Data->Delete();
        ntplVars->Delete();
        mmg_s->Delete();
        Photon_SC_Eta->Delete();
        Photon_SC_Phi->Delete();
        Photon_r9->Delete();
        isJanLooseMMG->Delete();
        isMultipleCandidate->Delete();
        Photon_Et->Delete();
        Photon_SC_rawEt->Delete();
        Photon_E->Delete();
        Photon_SC_brem->Delete();
        weight_pileUp->Delete();
        beta->Delete();
        gamma->Delete();
        mu->Delete();
        Eraw_Gamma->Delete();

}

void RooBifurcatedGauss(double * mean_value, double * mean_error, Double_t * sigmaEff_value, double * sigma_value, double * sigma_value_error, double * sigmaR_value, double * sigmaL_value, double * ChiSquare, double * minLogLikelihood, double * differenceLogLikelihood, TH1D * hh, int j, int EndCaps, TString temp, TTree* Tree_Data, double mean, double rms, TCanvas * c2, RooPlot* Erawframe, TF1 * f, string nomDossier, string nomFichier, double RangeMin, double RangeMax, double * JanChi2, double * DegreesOfFreedom, double * pValue, int r9sup, double MinVar, double MaxVar, string variableX, int isMC)
{
	RooRealVar * mmg_s = new RooRealVar("mmg_s", "s", -0.5, 0.5);
        RooRealVar * Photon_SC_Eta = new RooRealVar("Photon_SC_Eta", "Photon_SC_Eta", -3.0, 3.0, "");
        RooRealVar * Photon_SC_Phi = new RooRealVar("Photon_SC_Phi", "Photon_SC_Phi", -3.5, 3.5, "");
        RooRealVar * Photon_r9 = new RooRealVar("Photon_r9", "Photon_r9", 0.0, 1.0, "");
        RooRealVar * isJanLooseMMG = new RooRealVar("isJanLooseMMG", "isJanLooseMMG", 0.0, 1.0, "");
        RooRealVar * isMultipleCandidate = new RooRealVar("isMultipleCandidate", "isMultipleCandidate", 0.0, 1.0, "");
        RooRealVar * Photon_Et = new RooRealVar("Photon_Et", "Photon_Et", 0.0, 250.0, "");
        RooRealVar * Photon_SC_rawEt = new RooRealVar("Photon_SC_rawEt", "Photon_SC_rawEt", 0.0, 250.0, "");
        RooRealVar * Photon_E = new RooRealVar("Photon_E", "Photon_E", 0.0, 1000.0, "");
        RooRealVar * Photon_SC_brem = new RooRealVar("Photon_SC_brem", "Photon_SC_brem", 0.0, 15.0, "");
        RooRealVar * weight_pileUp = new RooRealVar("weight_pileUp", "weight_pileUp", 0.0, 100);

        RooArgSet *ntplVars = new RooArgSet(*mmg_s, *Photon_SC_Eta, *Photon_SC_Phi, *Photon_r9, *Photon_Et, *Photon_SC_rawEt, *Photon_E, *Photon_SC_brem, *isJanLooseMMG);



        ntplVars->add(*isMultipleCandidate);
        ntplVars->add(*weight_pileUp);


        //RooDataSet *Data = new RooDataSet("Data", "Data", Tree_Data, *ntplVars, "", "weight_pileUp"); //weighté
	//RooDataSet *Data = new RooDataSet("Data", "Data", Tree_Data, *ntplVars, "", "");
	RooDataSet *Data_subset = new RooDataSet("Data_subset", "Data_subset", Tree_Data, *ntplVars, temp, "weight_pileUp");

	RooRealVar * BifurGauss_mean = new RooRealVar("BifurGauss_mean","BifurGauss_mean",0.0,-0.3,0.2);
        RooRealVar * BifurGauss_sigmaR = new RooRealVar("BifurGauss_sigmaR","BifurGauss_sigmaR",0.5,0.0,1.0);
        RooRealVar * BifurGauss_sigmaL = new RooRealVar("BifurGauss_sigmaL","BifurGauss_sigmaL",0.5,0.0,1.0);


        // --- Build Bifurcated Gaussian PDF ---
        RooBifurGauss * BFG = new RooBifurGauss("Eraw_BifurGauss","mmg_s_BifurGauss",*mmg_s,*BifurGauss_mean,*BifurGauss_sigmaR,*BifurGauss_sigmaL);


	int fewBins = 1;
        Double_t Chi2J;

	RooFitResult *res;
	Double_t minNll;
		

        //for(int RightBinning = 600; RightBinning > 10; RightBinning -= 10)
        int RightBinning = 25;
	for(int i = 0; i<4; i++)
	{
                Erawframe->Clear();
                //Data_subset = (RooDataSet*)Data->reduce(*ntplVars, temp);
		BFG = new RooBifurGauss("BFG","BFG",*mmg_s,*BifurGauss_mean,*BifurGauss_sigmaR,*BifurGauss_sigmaR);

                Erawframe = mmg_s->frame();
                Data_subset->plotOn(Erawframe,Name("myhist"),Binning(RightBinning),DataError(RooAbsData::SumW2));
                //BFG->fitTo(*Data_subset, Range(RangeMin, RangeMax),Minos(true));
                //res = BFG->fitTo(*Data_subset, Range(RangeMin, RangeMax),Minos(true),Save(),SumW2Error(kFALSE));
		res = BFG->fitTo(*Data_subset, Range(RangeMin, RangeMax),Save(),SumW2Error(kTRUE));
		res->Print();
		minNll = res->minNll();
		BFG->plotOn(Erawframe,Name("mycurve"));
		Erawframe->Draw();


                Chi2J = chiSquare(Erawframe,(char *)"mycurve",(char *)"myhist",3, JanChi2, DegreesOfFreedom, pValue, &fewBins);
                if(fewBins == 0) break; 


        }
	cout << endl << "ICI !!!!!!!!!!!!!!!!!! "<<Data_subset->weight() << endl ;

        f = BFG->asTF( RooArgList(*mmg_s) );
	*mean_value  = BifurGauss_mean->getVal();
        *mean_error  = BifurGauss_mean->getError();
        *sigmaEff_value = effSigma(hh);
        //*sigma_value = sigwidth.getVal();
        *sigmaR_value = SigmaR(f, 0.0, 2.0);
        *sigmaL_value = SigmaL(f, 0.0, 2.0);
        *ChiSquare = Erawframe->chiSquare();
	*minLogLikelihood = minNll;


        double entries = hh->GetEntries();

        int entriesInt = (int) entries;

	
	double fxmax = f->GetMaximumX(0.8,1.2,1.E-10,100,false);
        //cout<<"////////// ---- FXMAX BG = "<<f->GetMaximumX(0.8,1.2,1.E-10,100,false)<<" ---- //////////"<<endl;

	TLatex latexLabel;
        latexLabel.SetTextSize(0.026);
        latexLabel.SetNDC();
	latexLabel.DrawLatex(0.13, 0.96, "CMS Preliminary 2011, #sqrt{s} = 7 TeV");
        if(isMC == 1) latexLabel.DrawLatex(0.17, 0.88, "Simulation");
        if(isMC == 0) latexLabel.DrawLatex(0.17, 0.88, "Data, #int L = 4.93 fb^{-1}"); //lumi a changer !!!

	if(EndCaps == 0) latexLabel.DrawLatex(0.17, 0.83,"ECAL Barrel");
        if(EndCaps == 1) latexLabel.DrawLatex(0.17, 0.83,"ECAL Endcaps");
        if(r9sup == 0 && EndCaps == 0) latexLabel.DrawLatex(0.17, 0.78,"E_{T}^{#gamma} > 25 GeV, r9 < 0,94");
        if(r9sup == 0 && EndCaps == 1) latexLabel.DrawLatex(0.17, 0.78,"E_{T}^{#gamma} > 25 GeV, r9 < 0,95");
        if(r9sup == 1 && EndCaps == 0) latexLabel.DrawLatex(0.17, 0.78,"E_{T}^{#gamma} > 25 GeV, r9 > 0,94");
        if(r9sup == 1 && EndCaps == 1) latexLabel.DrawLatex(0.17, 0.78,"E_{T}^{#gamma} > 25 GeV, r9 > 0,95");
        if(r9sup == 2) latexLabel.DrawLatex(0.17, 0.78,"E_{T}^{#gamma} > 25 GeV, All r9");
	
	latexLabel.DrawLatex(0.17, 0.73, Form("Entries = %d",entriesInt));


	latexLabel.DrawLatex(0.61, 0.90, Form("#color[4]{mean = %f #pm %f}",BifurGauss_mean->getVal(), BifurGauss_mean->getError()));
	latexLabel.DrawLatex(0.61, 0.85, Form("#color[4]{#sigma_{L} = %f #pm %f}",BifurGauss_sigmaL->getVal(), BifurGauss_sigmaL->getError()));
	latexLabel.DrawLatex(0.61, 0.80, Form("#color[4]{#sigma_{R} = %f #pm %f}",BifurGauss_sigmaL->getVal(), BifurGauss_sigmaL->getError()));	


	enregistrementPlots(nomDossier, nomFichier, EndCaps, j, c2);

	c2->Clear();

	////////// RooMCStudy //////////

	
	RooMCStudy* mcs = new RooMCStudy(*BFG,*mmg_s,FitModel(*BFG),Silence(),Binned()) ;
        RooChi2MCSModule chi2mod ;
        mcs->addModule(chi2mod) ;

        mcs->generateAndFit(2000,entriesInt) ;

        //RooPlot* frame = mcs->plotNLL(Bins(50));
        RooPlot* frame = mcs->plotNLL(Name("myhistTest"));
	frame->Draw();

	//RooRealVar * VArTest = (RooRealVar*) frame->getPlotVar(); //artefact

	
	RooHist* histTest = frame->getHist("myhistTest");
	cout<<endl<<"histTest->GetMean() = "<<histTest->GetMean()<<endl;

	*differenceLogLikelihood = fabs(histTest->GetMean() - minNll);

	TLatex *textResid = new TLatex();

        textResid = new TLatex();
        textResid->SetNDC();
        textResid->SetTextAlign(11);
        textResid->SetTextFont(42);
        textResid->SetTextSizePixels(17);
        textResid->SetTextSize(0.028);
        //if(EndCaps == 0) textResid->DrawLatex(0.80, 0.88, "Barrel");
        //if(EndCaps == 1) textResid->DrawLatex(0.77, 0.88, "End Caps");
        textResid->DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
        if(isMC == 1) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV");
        if(isMC == 0) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV, #int L = 4.93 fb^{-1}");
        //textResid->DrawLatex(0.16, 0.80, EndCapsR9Chain->c_str());
        //textResid->DrawLatex(0.16, 0.75, tempLegChain->c_str());

	textResid->DrawLatex(0.61, 0.90, Form("#color[4]{Fit_minll = %f}",minNll));
	textResid->DrawLatex(0.61, 0.85, Form("#color[4]{Distrib_mean = %f}",histTest->GetMean()));
	textResid->DrawLatex(0.61, 0.80, Form("#color[3]{Difference = %f}",*differenceLogLikelihood));
	
	//RooAbsRealLValue VArTest = frame->getPlotVar();


	enregistrementPlots(nomDossier, "LogLikelihood", EndCaps, j, c2);

	c2->Clear();


	////////// chi2 Pulls & Residuals //////////

	string nomDossierPull = "";
        if(EndCaps == 0) nomDossierPull = nomDossier + "EB/";
        if(EndCaps == 1) nomDossierPull = nomDossier + "EE/";


        RooHist* residuals = residHist(Erawframe,(char *)"myhist",(char *)"mycurve", false, nomDossierPull, j);
	

        //RooHist* pulls = pullHist(Erawframe,"myhist", "mycurve");
        RooHist* pulls = residHist(Erawframe,(char *)"myhist", (char *)"mycurve", true, nomDossierPull, j);


	residuals->Draw("AP");
	nomFichier = "Chi2Residuals";
	
	residuals->GetYaxis()->SetTitle("#chi^{2} Residuals");
        residuals->GetXaxis()->SetLabelFont(42);
        residuals->GetXaxis()->SetTitleFont(42);
        residuals->GetXaxis()->SetLabelSize(0.03);
        residuals->GetXaxis()->SetTitle("s");
        residuals->GetYaxis()->SetLabelFont(42);
        residuals->GetYaxis()->SetTitleOffset(1.24);
        residuals->GetYaxis()->SetTitleFont(42);
        residuals->GetYaxis()->SetLabelSize(0.03);
        //residuals->SetMarkerColor(4);
        //residuals->SetMarkerStyle(21);
        //residuals->SetMarkerSize(0.6);

	TLine *lineResid = new TLine(-0.5,0,0.5,0);
        lineResid->SetLineStyle(3);
        lineResid->Draw();


        //TLatex *textResid = new TLatex();

        textResid = new TLatex();
        textResid->SetNDC();
        textResid->SetTextAlign(11);
        textResid->SetTextFont(42);
        textResid->SetTextSizePixels(17);
        textResid->SetTextSize(0.028);
        //if(EndCaps == 0) textResid->DrawLatex(0.80, 0.88, "Barrel");
        //if(EndCaps == 1) textResid->DrawLatex(0.77, 0.88, "End Caps");
	textResid->DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
	if(isMC == 1) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV");
        if(isMC == 0) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV, #int L = 4.93 fb^{-1}");
        //textResid->DrawLatex(0.16, 0.80, EndCapsR9Chain->c_str());
        //textResid->DrawLatex(0.16, 0.75, tempLegChain->c_str());
	residuals->GetXaxis()->SetLimits(-0.5,0.5);
	
	enregistrementPlots(nomDossier, nomFichier, EndCaps, j, c2);

	c2->Clear();

	pulls->Draw("AP");
	nomFichier = "Chi2Pulls";
	
	pulls->GetYaxis()->SetTitle("#chi^{2} Pulls");
        pulls->GetXaxis()->SetLabelFont(42);
        pulls->GetXaxis()->SetTitleFont(42);
        pulls->GetXaxis()->SetLabelSize(0.03);
        pulls->GetXaxis()->SetTitle("s");
        pulls->GetYaxis()->SetLabelFont(42);
        pulls->GetYaxis()->SetTitleOffset(1.24);
        pulls->GetYaxis()->SetTitleFont(42);
        pulls->GetYaxis()->SetLabelSize(0.03);
        //pulls->SetMarkerColor(4);
        //pulls->SetMarkerStyle(21);
        //pulls->SetMarkerSize(0.6);

	lineResid = new TLine(-0.5,0,0.5,0);
        lineResid->SetLineStyle(3);
        lineResid->Draw();

        textResid = new TLatex();
        textResid->SetNDC();
        textResid->SetTextAlign(11);
        textResid->SetTextFont(42);
        textResid->SetTextSizePixels(17);
        textResid->SetTextSize(0.028);
        //if(EndCaps == 0) textResid->DrawLatex(0.80, 0.88, "Barrel");
        //if(EndCaps == 1) textResid->DrawLatex(0.77, 0.88, "End Caps");	
	textResid->DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
	if(isMC == 1) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV");
        if(isMC == 0) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV, #int L = 4.93 fb^{-1}");
        //textResid->DrawLatex(0.16, 0.80, EndCapsR9Chain->c_str());
        //textResid->DrawLatex(0.16, 0.75, tempLegChain->c_str());
	pulls->GetXaxis()->SetLimits(-0.5,0.5);	



	c2->Clear();

        //delete EndCapsR9Chain;
        //delete tempLegChain;
        res->Delete();
        mcs->Delete();
        frame->Delete();
        residuals->Delete();
        pulls->Delete();
        textResid->Delete();
        lineResid->Delete();
        Data_subset->Delete();
        //Data->Delete();
        ntplVars->Delete();
        mmg_s->Delete();
        Photon_SC_Eta->Delete();
        Photon_SC_Phi->Delete();
        Photon_r9->Delete();
        isJanLooseMMG->Delete();
        isMultipleCandidate->Delete();
        Photon_Et->Delete();
        Photon_SC_rawEt->Delete();
        Photon_E->Delete();
        Photon_SC_brem->Delete();
        weight_pileUp->Delete();
        BifurGauss_mean->Delete();
        BifurGauss_sigmaL->Delete();
        BifurGauss_sigmaR->Delete();
        BFG->Delete();

}


void RooSumGauss(double * mean_value, double * mean_error, Double_t * sigmaEff_value, double * sigma_value, double * sigma_value_error, double * sigmaR_value, double * sigmaL_value, double * ChiSquare, TH1D * hh, int j, int EndCaps, TString temp, TTree* Tree_Data, double mean, double rms, TCanvas * c2, RooPlot* Erawframe, TF1 * f, string nomDossier, string nomFichier, double RangeMin, double RangeMax, double * JanChi2, double * DegreesOfFreedom, double * pValue, int r9sup, double MinVar, double MaxVar, string variableX, int isMC)
{

	RooRealVar * mmg_s = new RooRealVar("mmg_s", "s", -0.5, 0.5);
        RooRealVar * Photon_SC_Eta = new RooRealVar("Photon_SC_Eta", "Photon_SC_Eta", -3.0, 3.0, ""); 
        RooRealVar * Photon_SC_Phi = new RooRealVar("Photon_SC_Phi", "Photon_SC_Phi", -3.5, 3.5, ""); 
        RooRealVar * Photon_r9 = new RooRealVar("Photon_r9", "Photon_r9", 0.0, 1.0, ""); 
        RooRealVar * isJanLooseMMG = new RooRealVar("isJanLooseMMG", "isJanLooseMMG", 0.0, 1.0, ""); 
        RooRealVar * isMultipleCandidate = new RooRealVar("isMultipleCandidate", "isMultipleCandidate", 0.0, 1.0, ""); 
        RooRealVar * Photon_Et = new RooRealVar("Photon_Et", "Photon_Et", 0.0, 250.0, "");     
        RooRealVar * Photon_SC_rawEt = new RooRealVar("Photon_SC_rawEt", "Photon_SC_rawEt", 0.0, 250.0, ""); 
        RooRealVar * Photon_E = new RooRealVar("Photon_E", "Photon_E", 0.0, 1000.0, "");
        RooRealVar * Photon_SC_brem = new RooRealVar("Photon_SC_brem", "Photon_SC_brem", 0.0, 15.0, "");
        RooRealVar * weight_pileUp = new RooRealVar("weight_pileUp", "weight_pileUp", 0.0, 10000000); 

        RooArgSet *ntplVars = new RooArgSet(*mmg_s, *Photon_SC_Eta, *Photon_SC_Phi, *Photon_r9, *Photon_Et, *Photon_SC_rawEt, *Photon_E, *Photon_SC_brem, *isJanLooseMMG);

        ntplVars->add(*isMultipleCandidate);
        ntplVars->add(*weight_pileUp);
	
	RooDataSet *Data = new RooDataSet("Data", "Data", Tree_Data, *ntplVars, "", "weight_pileUp");
	RooDataSet* Data_subset = (RooDataSet*)Data->reduce(*ntplVars, temp);

 
	// --- Parameters ---
	//RooRealVar sigmean1("sigmean1","sigmean1",1.0,0.6,1.4);
	//RooRealVar sigwidth1("sigwidth1","sigwidth1",0.5,0.0,1.0) ;
 
	RooRealVar * sigmean1 = new RooRealVar("sigmean1","sigmean1",1.0,0.6,1.4);
        RooRealVar * sigwidth1 = new RooRealVar("sigwidth1","sigwidth1",0.5,0.0,1.0) ;

	// --- Build Gaussian PDF ---
	RooGaussian Eraw_Gaussian1("Eraw_Gaussian1","mmg_s_Gaussian1",*mmg_s,*sigmean1,*sigwidth1) ;
	
	// --- Parameters ---
        //RooRealVar sigmean2("sigmean2","sigmean2",0.9,0.6,1.4);
	//RooRealVar sigwidth2("sigwidth2","sigwidth2",0.5,0.0,1.0) ;
        RooRealVar * sigmean2 = new RooRealVar("sigmean2","sigmean2",0.9,0.6,1.4);
	RooRealVar * sigwidth2 = new RooRealVar("sigwidth2","sigwidth2",0.5,0.0,1.0) ;


	// --- Build Gaussian PDF ---
        RooGaussian Eraw_Gaussian2("Eraw_Gaussian2","mmg_s_Gaussian2",*mmg_s,*sigmean2,*sigwidth2) ;


	// --- Parameters ---
        //RooRealVar sigmean3("sigmean3","sigmean3",1.1,0.6,1.4);
	//RooRealVar sigwidth3("sigwidth3","sigwidth3",0.5,0.0,1.0) ;
	RooRealVar * sigmean3 = new RooRealVar("sigmean3","sigmean3",1.1,0.6,1.4);
	RooRealVar * sigwidth3 = new RooRealVar("sigwidth3","sigwidth3",0.5,0.0,1.0) ;


        // --- Build Gaussian PDF ---
        RooGaussian Eraw_Gaussian3("Eraw_Gaussian3","mmg_s_Gaussian3",*mmg_s,*sigmean3,*sigwidth3) ;


	// --- Parameters ---
        //RooRealVar sigmean4("sigmean4","sigmean4",0.8,0.6,1.4);
        //RooRealVar sigwidth4("sigwidth4","sigwidth4",0.5,0.0,1.0);
	RooRealVar * sigmean4 = new RooRealVar("sigmean4","sigmean4",0.8,0.6,1.4);
	RooRealVar * sigwidth4 = new RooRealVar("sigwidth4","sigwidth4",0.5,0.0,1.0);


        // --- Build Gaussian PDF ---
        RooGaussian Eraw_Gaussian4("Eraw_Gaussian4","mmg_s_Gaussian4",*mmg_s,*sigmean4,*sigwidth4) ;	


	// --- Parameters ---
        //RooRealVar sigmean5("sigmean5","sigmean5",1.2,0.6,1.4);
        //RooRealVar sigwidth5("sigwidth5","sigwidth5",0.5,0.0,1.0);
	RooRealVar * sigmean5 = new RooRealVar("sigmean5","sigmean5",1.2,0.6,1.4);
	RooRealVar * sigwidth5 = new RooRealVar("sigwidth5","sigwidth5",0.5,0.0,1.0);

        // --- Build Gaussian PDF ---
        RooGaussian Eraw_Gaussian5("Eraw_Gaussian5","mmg_s_Gaussian5",*mmg_s,*sigmean5,*sigwidth5) ;


	// --- Parameters ---
        //RooRealVar sigmean6("sigmean6","sigmean6",0.7,0.6,1.4);
        //RooRealVar sigwidth6("sigwidth6","sigwidth6",0.5,0.0,1.0);
	RooRealVar * sigmean6 = new RooRealVar("sigmean6","sigmean6",0.7,0.6,1.4);
	RooRealVar * sigwidth6 = new RooRealVar("sigwidth6","sigwidth6",0.5,0.0,1.0);

        // --- Build Gaussian PDF ---
        RooGaussian Eraw_Gaussian6("Eraw_Gaussian6","mmg_s_Gaussian6",*mmg_s,*sigmean6,*sigwidth6) ;


	// --- Parameters ---
        //RooRealVar sigmean7("sigmean7","sigmean7",1.3,0.6,1.4);
        //RooRealVar sigwidth7("sigwidth7","sigwidth7",0.5,0.0,1.0);
	RooRealVar * sigmean7 = new RooRealVar("sigmean7","sigmean7",1.3,0.6,1.4);
	RooRealVar * sigwidth7 = new RooRealVar("sigwidth7","sigwidth7",0.5,0.0,1.0);

        // --- Build Gaussian PDF ---
        RooGaussian Eraw_Gaussian7("Eraw_Gaussian7","mmg_s_Gaussian7",*mmg_s,*sigmean7,*sigwidth7) ;




	// --- Build Sum of Gaussian PDF ---
/*
	RooRealVar * g1Frac = new RooRealVar("g1Frac","g1Frac",0.1428);
	RooRealVar * g2Frac = new RooRealVar("g2Frac","g2Frac",0.1428);
	RooRealVar * g3Frac = new RooRealVar("g3Frac","g3Frac",0.1428);
	RooRealVar * g4Frac = new RooRealVar("g4Frac","g4Frac",0.1428);
	RooRealVar * g5Frac = new RooRealVar("g5Frac","g5Frac",0.1428);
        RooRealVar * g6Frac = new RooRealVar("g6Frac","g6Frac",0.1428);
        RooRealVar * g7Frac = new RooRealVar("g7Frac","g7Frac",0.1428);	
*/
	double princ = 0.99;

	double frac = (1.0 - princ) / 6.0;

	RooRealVar * g1Frac = new RooRealVar("g1Frac","g1Frac",princ);
        RooRealVar * g2Frac = new RooRealVar("g2Frac","g2Frac",frac);
        RooRealVar * g3Frac = new RooRealVar("g3Frac","g3Frac",frac);
        RooRealVar * g4Frac = new RooRealVar("g4Frac","g4Frac",frac);
        RooRealVar * g5Frac = new RooRealVar("g5Frac","g5Frac",frac);
        RooRealVar * g6Frac = new RooRealVar("g6Frac","g6Frac",frac);
        RooRealVar * g7Frac = new RooRealVar("g7Frac","g7Frac",frac);	



	RooAddPdf SumGaussians("SumGaussians","SumGaussians", RooArgList(Eraw_Gaussian1,Eraw_Gaussian2,Eraw_Gaussian3,Eraw_Gaussian4,Eraw_Gaussian5,Eraw_Gaussian6,Eraw_Gaussian7), RooArgList(*g1Frac,*g2Frac,*g3Frac,*g4Frac,*g5Frac,*g6Frac,*g7Frac));


        //Eraw_CrystalBall.fitTo(*Data_subset, Range(0.6,1.2));
        //Eraw_CrystalBall.fitTo(Data_subset, Range(0.5,1.15)); //ATTENTION NE PAS SUPPRIMER!!!!
        double maxDistri = hh->GetMaximumBin() * 0.01;
        //SumGaussians.fitTo(*Data_subset, Range(mean - 0.35 * rms,mean + 0.35 * rms));
        SumGaussians.fitTo(*Data_subset, Range(maxDistri - 0.6 * rms,maxDistri + 0.6 * rms));
        //SumGaussians.fitTo(Data_subset, Range(0.4,1.6));
	SumGaussians.plotOn(Erawframe);
        Erawframe->Draw();



        //RooDataHist data("data", "dataset with mmg_s",mmg_s,hh);
        //TF1 * f = SumGaussians.asTF( RooArgList(mmg_s) );
	f = SumGaussians.asTF( RooArgList(*mmg_s) );        

	//*mean_value  = sigmean.getVal();
        //*mean_error  = sigmean.getError();
        *sigmaEff_value = effSigma(hh);
        //*sigma_value = sigwidth.getVal();
        *sigmaR_value = SigmaR(f, 0.0, 2.0);
        *sigmaL_value = SigmaL(f, 0.0, 2.0);
        *ChiSquare = Erawframe->chiSquare();


	cout<<"SigmaR(f, 0.0, 2.0) = "<<SigmaR(f, 0.0, 2.0)<<endl;
	cout<<"SigmaL(f, 0.0, 2.0) = "<<SigmaL(f, 0.0, 2.0)<<endl;

	double fxmax = f->GetMaximumX(0.4,1.6,1.E-10,100,false);

        double entries = hh->GetEntries();

        int entriesInt = (int) entries;

	
	TLatex latexLabel;
        latexLabel.SetTextSize(0.028);
        latexLabel.SetNDC();
        latexLabel.DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
        if(EndCaps == 0) latexLabel.DrawLatex(0.16, 0.85, "Barrel");
        if(EndCaps == 1) latexLabel.DrawLatex(0.16, 0.85, "EndCaps");
        latexLabel.DrawLatex(0.16, 0.80, Form("Entries = %d",entriesInt));
        cout<<endl<<"entries = "<<entries<<endl;
        cout<<endl<<"Mean = "<<hh->GetMean()<<endl;
        latexLabel.DrawLatex(0.16, 0.75, Form("Mean = %f +/- %f",hh->GetMean(), hh->GetMeanError()));
        latexLabel.DrawLatex(0.16, 0.70, Form("RMS = %f +/- %f",hh->GetRMS(), hh->GetRMSError()));
        latexLabel.DrawLatex(0.61, 0.90, Form("#color[4]{m_{0} = %f}",fxmax));
        //latexLabel.DrawLatex(0.61, 0.85, Form("#color[4]{ = %f +/- %f}",sigwidth.getVal(), sigwidth.getError()));
        latexLabel.DrawLatex(0.61, 0.80, Form("#color[4]{#chi^{2} = %f}",Erawframe->chiSquare()));


	enregistrementPlots(nomDossier, nomFichier, EndCaps, j, c2);

	Data_subset->Delete();
        Data->Delete();
        ntplVars->Delete();
        mmg_s->Delete();
        Photon_SC_Eta->Delete();
        Photon_r9->Delete();
        isJanLooseMMG->Delete();
	isMultipleCandidate->Delete();
        Photon_Et->Delete();
	sigmean1->Delete();
        sigwidth1->Delete();
	sigmean2->Delete();
        sigwidth2->Delete();
	sigmean3->Delete();
        sigwidth3->Delete();
	sigmean4->Delete();
        sigwidth4->Delete();
	sigmean5->Delete();
        sigwidth5->Delete();
	sigmean6->Delete();
        sigwidth6->Delete();
	sigmean7->Delete();
	sigwidth7->Delete();
	g1Frac->Delete();
	g2Frac->Delete();
	g3Frac->Delete();
	g4Frac->Delete();
	g5Frac->Delete();
	g6Frac->Delete();
	g7Frac->Delete();



}




void RooGenericPDF(double * mean_value, double * mean_error, Double_t * sigmaEff_value, double * sigma_value, double * sigma_value_error, double * sigmaR_value, double * sigmaL_value, double * ChiSquare, TH1D * hh, int j, int EndCaps, TString temp, TTree* Tree_Data, double mean, double rms, TCanvas * c2, RooPlot* Erawframe, TF1 * f, string nomDossier, string nomFichier)
{

	RooRealVar * mmg_s = new RooRealVar("mmg_s", "s", -0.5, 0.5);
        RooRealVar * Photon_SC_Eta = new RooRealVar("Photon_SC_Eta", "Photon_SC_Eta", -3.0, 3.0, "");
        RooRealVar * Photon_r9 = new RooRealVar("Photon_r9", "Photon_r9", 0.0, 1.0, "");
        RooRealVar * isJanLooseMMG = new RooRealVar("isJanLooseMMG", "isJanLooseMMG", 0.0, 1.0, "");
	RooRealVar * isMultipleCandidate = new RooRealVar("isMultipleCandidate", "isMultipleCandidate", 0.0, 1.0, "");
        RooRealVar * Photon_Et = new RooRealVar("Photon_Et", "Photon_Et", 0.0, 250.0, "");


        RooArgSet *ntplVars = new RooArgSet(*mmg_s, *Photon_SC_Eta, *Photon_r9, *Photon_Et, *isJanLooseMMG, *isMultipleCandidate);
        //RooArgSet ntplVars(mmg_s, Photon_SC_Eta, Photon_r9, Photon_Et, isJanLooseMMG, Photon_Et);
        RooDataSet *Data = new RooDataSet("Data", "Data", Tree_Data, *ntplVars);
        //RooDataSet Data("Data", "Data", Tree_Data, ntplVars);


        RooDataSet* Data_subset = (RooDataSet*)Data->reduce(*ntplVars, temp);
        //RooDataSet Data_subset = (RooDataSet)Data.reduce(ntplVars, temp);
        //RooPlot* Erawframe = mmg_s.frame();
        Erawframe = mmg_s->frame();
        Data_subset->plotOn(Erawframe);
        //Data_subset.plotOn(Erawframe);

        Erawframe->Draw();



 
	// --- Parameters 1/Gaussian ---
	//RooRealVar GenericPDF_mean("GenericPDF_mean","GenericPDF_mean",1.0,0.6,1.4);
	//RooRealVar GenericPDF_mean("GenericPDF_mean","GenericPDF_mean", mean, mean - rms, mean + rms, "GeV");
	//RooRealVar GenericPDF_width("GenericPDF_width","GenericPDF_width",0.5,0.0,1.0);

	RooRealVar * GenericPDF_mean = new RooRealVar("GenericPDF_mean","GenericPDF_mean", mean, mean - rms, mean + rms, "GeV");
        RooRealVar * GenericPDF_width = new RooRealVar("GenericPDF_width","GenericPDF_width",0.5,0.0,1.0);

 
	// --- Build 1/Gaussian PDF ---
	RooGenericPdf Eraw_GenericPdf("Eraw_GenericPdf","mmg_s_GenericPdf","(GenericPDF_width * sqrt(2*3.14159265358979323846))*exp((mmg_s - GenericPDF_mean)*(mmg_s - GenericPDF_mean)/(2*GenericPDF_width*GenericPDF_width))",RooArgSet(*mmg_s, *GenericPDF_mean,*GenericPDF_width));



	// --- Parameters CB ---
/*	RooRealVar Eraw_CB_m0("Eraw_CB_m0", "CB #Delta m_{0}", mean, mean - rms, mean + rms, "GeV");
        RooRealVar Eraw_CB_sigma("Eraw_CB_sigma", "CB ", rms, 0.0, 2 * rms, "GeV");
        RooRealVar Eraw_CB_alpha("Eraw_CB_alpha", "CB #alpha", 0.9, -5.0, 5.0);
        RooRealVar Eraw_CB_n("Eraw_CB_n", "CB n", 10.0, 0.0, 500.0);    
*/
	RooRealVar * Eraw_CB_m0 = new RooRealVar("Eraw_CB_m0", "CB #Delta m_{0}", mean, mean - rms, mean + rms, "GeV");
        RooRealVar * Eraw_CB_sigma = new RooRealVar("Eraw_CB_sigma", "CB ", rms, 0.0, 2 * rms, "GeV");
        RooRealVar * Eraw_CB_alpha = new RooRealVar("Eraw_CB_alpha", "CB #alpha", 0.9, -5.0, 5.0);
        RooRealVar * Eraw_CB_n = new RooRealVar("Eraw_CB_n", "CB n", 10.0, 0.0, 500.0); 

/*
        //Ne pas supprimer !!!!//
        RooRealVar Eraw_CB_m0("Eraw_CB_m0", "CB #Delta m_{0}", 1.0, 0.9, 1.05, "GeV");
        RooRealVar Eraw_CB_sigma("Eraw_CB_sigma", "CB ", 0.0215, 0.0, 0.5, "GeV");
        RooRealVar Eraw_CB_alpha("Eraw_CB_alpha", "CB #alpha", 0.9, 0.0, 3.0);
        RooRealVar Eraw_CB_n("Eraw_CB_n", "CB n", 10.0, 0.0, 60.0);
*/

        //RooRealVar Eraw_CB_m0("mmg_s_CB_m0", "CB #Delta m_{0}", 1.0, -5.0, 5.0, "GeV");
        //RooRealVar Eraw_CB_sigma("mmg_s_CB_sigma", "CB ", 0.45, 0.01, 0.5, "GeV");
        //RooRealVar Eraw_CB_alpha("mmg_s_CB_alpha", "CB #alpha", -1.0, -10.01, 10.0);
        //RooRealVar Eraw_CB_n("mmg_s_CB_n", "CB n", 2.0, 0.5, 500.0);


	// --- Build CB PDF ---
        RooCBShape Eraw_CrystalBall("Eraw_CrystalBall","mmg_s_CrystalBall", *mmg_s, *Eraw_CB_m0, *Eraw_CB_sigma, *Eraw_CB_alpha, *Eraw_CB_n);

	RooProdPdf prod("CB/Gaussian","CB/Gaussian",RooArgList(Eraw_CrystalBall,Eraw_GenericPdf));



        //Eraw_CrystalBall.fitTo(*Data_subset, Range(0.6,1.2));
        //prod.fitTo(Data_subset, Range(0.6,1.2)); //ATTENTION NE PAS SUPPRIMER!!!!
	double maxDistri = hh->GetMaximumBin() * 0.01;
        //prod.fitTo(*Data_subset, Range(mean - 0.35 * rms,mean + 0.35 * rms));
        prod.fitTo(*Data_subset, Range(maxDistri - 0.6 * rms,maxDistri + 0.6 * rms));
        prod.plotOn(Erawframe);
        Erawframe->Draw();



        //RooDataHist data("data", "dataset with mmg_s",mmg_s,hh);
        //TF1 * f = prod.asTF( RooArgList(mmg_s) );
        f = prod.asTF( RooArgList(*mmg_s) );
	*mean_value  = Eraw_CB_m0->getVal();
        *mean_error  = Eraw_CB_m0->getError();
        *sigmaEff_value = effSigma(hh);
        *sigma_value = Eraw_CB_sigma->getVal();
        *sigmaR_value = SigmaR(f, 0.0, 2.0);
        *sigmaL_value = SigmaL(f, 0.0, 2.0);
        *ChiSquare = Erawframe->chiSquare();


        double entries = hh->GetEntries();

        int entriesInt = (int) entries;

	
	TLatex latexLabel;
        latexLabel.SetTextSize(0.028);
        latexLabel.SetNDC();
        latexLabel.DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
        if(EndCaps == 0) latexLabel.DrawLatex(0.16, 0.85, "Barrel");
        if(EndCaps == 1) latexLabel.DrawLatex(0.16, 0.85, "EndCaps");
        latexLabel.DrawLatex(0.16, 0.80, Form("Entries = %d",entriesInt));
        cout<<endl<<"entries = "<<entries<<endl;
        cout<<endl<<"Mean = "<<hh->GetMean()<<endl;
        
	latexLabel.DrawLatex(0.16, 0.75, Form("Mean = %f +/- %f",hh->GetMean(), hh->GetMeanError()));
        latexLabel.DrawLatex(0.16, 0.70, Form("RMS = %f +/- %f",hh->GetRMS(), hh->GetRMSError()));
        latexLabel.SetTextSize(0.023);
        latexLabel.DrawLatex(0.61, 0.90, Form("#color[4]{mean_{1/Gauss} = %f +/- %f}",GenericPDF_mean->getVal(), GenericPDF_mean->getError()));
	latexLabel.DrawLatex(0.61, 0.85, Form("#color[4]{width_{1/Gauss} = %f +/- %f}",GenericPDF_width->getVal(), GenericPDF_width->getError()));
	latexLabel.DrawLatex(0.61, 0.80, Form("#color[4]{m_{0}_{CB} = %f +/- %f}",Eraw_CB_m0->getVal(), Eraw_CB_m0->getError()));
	latexLabel.DrawLatex(0.61, 0.75, Form("#color[4]{#sigma_{CB} = %f +/- %f}",Eraw_CB_sigma->getVal(), Eraw_CB_sigma->getError()));
	latexLabel.DrawLatex(0.61, 0.70, Form("#color[4]{#alpha_{CB} = %f +/- %f}",Eraw_CB_alpha->getVal(), Eraw_CB_alpha->getError()));
	latexLabel.DrawLatex(0.61, 0.65, Form("#color[4]{n_{CB} = %f +/- %f}",Eraw_CB_n->getVal(), Eraw_CB_n->getError()));
        latexLabel.DrawLatex(0.61, 0.60, Form("#color[4]{#chi^{2} = %f}",Erawframe->chiSquare()));


	enregistrementPlots(nomDossier, nomFichier, EndCaps, j, c2);

	Data_subset->Delete();
        Data->Delete();
        ntplVars->Delete();
        mmg_s->Delete();
        Photon_SC_Eta->Delete();
        Photon_r9->Delete();
        isJanLooseMMG->Delete();
	isMultipleCandidate->Delete();
        Photon_Et->Delete();
	GenericPDF_mean->Delete();
	GenericPDF_width->Delete();
	Eraw_CB_m0->Delete();
        Eraw_CB_sigma->Delete();
        Eraw_CB_alpha->Delete();
        Eraw_CB_n->Delete();


}

void RooLandau2(double * mean_value, double * mean_error, Double_t * sigmaEff_value, double * sigma_value, double * sigma_value_error, double * sigmaR_value, double * sigmaL_value, double * ChiSquare, TH1D * hh, int j, int EndCaps, TString temp, TTree* Tree_Data, double mean, double rms, TCanvas * c2, RooPlot* Erawframe, TF1 * f, string nomDossier, string nomFichier, double RangeMin, double RangeMax, double * JanChi2, double * DegreesOfFreedom, double * pValue, int r9sup, double MinVar, double MaxVar, string variableX, int isMC)
{

	RooRealVar * mmg_s = new RooRealVar("mmg_s", "s", -0.5, 0.5);
        RooRealVar * Photon_SC_Eta = new RooRealVar("Photon_SC_Eta", "Photon_SC_Eta", -3.0, 3.0, "");
        RooRealVar * Photon_SC_Phi = new RooRealVar("Photon_SC_Phi", "Photon_SC_Phi", -3.5, 3.5, "");
	RooRealVar * Photon_r9 = new RooRealVar("Photon_r9", "Photon_r9", 0.0, 1.0, "");
        RooRealVar * isJanLooseMMG = new RooRealVar("isJanLooseMMG", "isJanLooseMMG", 0.0, 1.0, "");
        RooRealVar * isMultipleCandidate = new RooRealVar("isMultipleCandidate", "isMultipleCandidate", 0.0, 1.0, "");
	RooRealVar * Photon_Et = new RooRealVar("Photon_Et", "Photon_Et", 0.0, 250.0, "");	
	RooRealVar * Photon_SC_rawEt = new RooRealVar("Photon_SC_rawEt", "Photon_SC_rawEt", 0.0, 250.0, ""); 
        RooRealVar * Photon_E = new RooRealVar("Photon_E", "Photon_E", 0.0, 1000.0, ""); 
        RooRealVar * Photon_SC_brem = new RooRealVar("Photon_SC_brem", "Photon_SC_brem", 0.0, 15.0, ""); 
        RooRealVar * weight_pileUp = new RooRealVar("weight_pileUp", "weight_pileUp", 0.0, 10000000); 

	RooArgSet *ntplVars = new RooArgSet(*mmg_s, *Photon_SC_Eta, *Photon_SC_Phi, *Photon_r9, *Photon_Et, *Photon_SC_rawEt, *Photon_E, *Photon_SC_brem, *isJanLooseMMG);

        ntplVars->add(*isMultipleCandidate);
        ntplVars->add(*weight_pileUp);



        //RooArgSet *ntplVars = new RooArgSet(*mmg_s, *Photon_SC_Eta, *Photon_SC_Phi, *Photon_r9, *Photon_Et, *Photon_SC_rawEt, *Photon_E, *Photon_SC_brem, *isJanLooseMMG, *isMultipleCandidate, *weight_pileUp);



        //RooArgSet *ntplVars = new RooArgSet(*mmg_s, *Photon_SC_Eta, *Photon_SC_Phi, *Photon_r9, *Photon_Et, *isJanLooseMMG, *isMultipleCandidate, *weight_pileUp);
	//RooArgSet ntplVars(mmg_s, Photon_SC_Eta, Photon_r9, Photon_Et, isJanLooseMMG, Photon_Et);
        //RooDataSet *Data = new RooDataSet("Data", "Data", Tree_Data, *ntplVars); //non weighté
	RooDataSet *Data = new RooDataSet("Data", "Data", Tree_Data, *ntplVars, "", "weight_pileUp"); //weighté
	//RooDataSet Data("Data", "Data", Tree_Data, ntplVars);


        RooDataSet* Data_subset = (RooDataSet*)Data->reduce(*ntplVars, temp);
	//RooDataSet Data_subset = (RooDataSet)Data.reduce(ntplVars, temp);
        //RooPlot* Erawframe = mmg_s.frame();
        
/*	
	Erawframe = mmg_s->frame();
	Data_subset->plotOn(Erawframe,Name("myhist"),Binning(300));
        Erawframe->Draw();
*/

	RooRealVar * Eraw_Landau_m0 = new RooRealVar("Eraw_Landau_m0", "Landau #Delta m_{0}", mean, mean - rms, mean + rms, "GeV");
        RooRealVar * Eraw_Landau_sigma = new RooRealVar("Eraw_Landau_sigma", "Landau ", rms, 0.0, 2 * rms, "GeV");


/*
	//Ne pas supprimer !!!!//
	RooRealVar Eraw_Landau_m0("Eraw_Landau_m0", "Landau #Delta m_{0}", 1.0, 0.9, 1.05, "GeV");
        RooRealVar Eraw_Landau_sigma("Eraw_Landau_sigma", "Landau ", 0.0215, 0.0, 0.5, "GeV");
        RooRealVar Eraw_Landau_alpha("Eraw_Landau_alpha", "Landau #alpha", 0.9, 0.0, 3.0);
        RooRealVar Eraw_Landau_n("Eraw_Landau_n", "Landau n", 10.0, 0roo.0, 60.0);
*/
        RooLandau * Eraw_Landau = new RooLandau("Eraw_Landau","mmg_s_Landau", *mmg_s, *Eraw_Landau_m0, *Eraw_Landau_sigma);

	int fewBins = 1;
	Double_t Chi2J;

	for(int RightBinning = 600; RightBinning > 10; RightBinning -= 10)
	{
		cout<<endl<<endl<<endl<<"RightBinning = "<<RightBinning<<endl<<endl<<endl;
		Erawframe->Clear();
		Data_subset = (RooDataSet*)Data->reduce(*ntplVars, temp);		
		Eraw_Landau = new RooLandau("Eraw_Landau","mmg_s_Landau", *mmg_s, *Eraw_Landau_m0, *Eraw_Landau_sigma); 

		Erawframe = mmg_s->frame();
		Data_subset->plotOn(Erawframe,Name("myhist"),Binning(RightBinning));
		//Erawframe->Draw();
        	//Eraw_Landau.fitTo(*Data_subset, Range(0.6,1.2));
		//Eraw_Landau.fitTo(Data_subset, Range(0.5,1.15)); //ATTENTION NE PAS SUPPRIMER!!!!
		//Eraw_Landau.fitTo(Data_subset, Range(mean - 5 * rms,mean + 5 * rms));
		//double maxDistri = hh->GetMaximumBin() * 0.01;
		//Eraw_Landau.fitTo(*Data_subset, Range(mean - 0.35 * rms,mean + 0.35 * rms));
		//Eraw_Landau.fitTo(*Data_subset, Range(maxDistri - 0.6 * rms,maxDistri + 0.6 * rms));
		Eraw_Landau->fitTo(*Data_subset, Range(RangeMin, RangeMax));
		//RooArgSet* Eraw_Landau_param = Eraw_Landau.getVariables();
		//Eraw_Landau_param->Print("v");
		//Eraw_Landau.plotOn(Erawframe);
        	Eraw_Landau->plotOn(Erawframe,Name("mycurve"));
		Erawframe->Draw();


		Chi2J = chiSquare(Erawframe,(char *)"mycurve",(char *)"myhist",2, JanChi2, DegreesOfFreedom, pValue, &fewBins);
		if(fewBins == 0) break; 


	}


	//RooDataHist data("data", "dataset with mmg_s",*mmg_s,hh);
	//TF1 * f = Eraw_Landau.asTF( RooArgList(mmg_s) );
	f = Eraw_Landau->asTF( RooArgList(*mmg_s) );
	//ftest = Eraw_Landau.asTF( RooArgList(*mmg_s) );
	//cout<<endl<<"ftest->GetXmax() = "<<ftest->GetXmax()<<endl;
	//cout<<endl<<"ftest->GetMaximum = "<<ftest->GetMaximum(0.0, 2.0,1.E-10,100,false)<<endl;
	*mean_value  = Eraw_Landau_m0->getVal();
        *mean_error  = Eraw_Landau_m0->getError();
        *sigmaEff_value = effSigma(hh);
        *sigma_value = Eraw_Landau_sigma->getVal();
	*sigma_value_error = Eraw_Landau_sigma->getError();
        *sigmaR_value = SigmaR(f, 0.0, 2.0);
        *sigmaL_value = SigmaL(f, 0.0, 2.0);
        *ChiSquare = Erawframe->chiSquare();

 
        double entries = hh->GetEntries();
        
        int entriesInt = (int) entries;


        TLatex latexLabel;
        latexLabel.SetTextSize(0.026);
        latexLabel.SetNDC();
        latexLabel.DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
	if(isMC == 1) latexLabel.DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV");
        if(isMC == 0) latexLabel.DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV, #int L = 4.93 fb^{-1}");
        //if(EndCaps == 0) latexLabel.DrawLatex(0.16, 0.85, "Barrel");
        //if(EndCaps == 1) latexLabel.DrawLatex(0.16, 0.85, "EndCaps");   
        string * EndCapsR9Chain(0);
        EndCapsR9Chain = new string;
        if(EndCaps == 0) *EndCapsR9Chain = "Barrel, ";
        if(EndCaps == 1) *EndCapsR9Chain = "Endcaps, ";
        if(r9sup == 0) *EndCapsR9Chain += "Low r9";
        if(r9sup == 1) *EndCapsR9Chain += "High r9";
        if(r9sup == 2) *EndCapsR9Chain += "All r9";
	string * tempLegChain(0);
        tempLegChain = new string;
        *tempLegChain = DoubleToString(MinVar) + " <" + variableX + "< "  + DoubleToString(MaxVar);     
        latexLabel.DrawLatex(0.16, 0.80, EndCapsR9Chain->c_str());
        latexLabel.DrawLatex(0.16, 0.75, tempLegChain->c_str());

        latexLabel.DrawLatex(0.16, 0.70, Form("Entries = %d",entriesInt));
        cout<<endl<<"entries = "<<entries<<endl;
        cout<<endl<<"Mean = "<<hh->GetMean()<<endl;
        latexLabel.DrawLatex(0.16, 0.65, Form("Mean = %f +/- %f",hh->GetMean(), hh->GetMeanError()));
        latexLabel.DrawLatex(0.16, 0.60, Form("RMS = %f +/- %f",hh->GetRMS(), hh->GetRMSError()));
	latexLabel.DrawLatex(0.61, 0.90, Form("#color[4]{m_{0} = %f +/- %f}",Eraw_Landau_m0->getVal(), Eraw_Landau_m0->getError()));
        latexLabel.DrawLatex(0.61, 0.85, Form("#color[4]{#sigma_{landau} = %f +/- %f}",Eraw_Landau_sigma->getVal(), Eraw_Landau_sigma->getError()));
        latexLabel.DrawLatex(0.61, 0.80, Form("#color[4]{Default #chi^{2}/ndf = %f}",Erawframe->chiSquare()));
	latexLabel.DrawLatex(0.61, 0.75, Form("#color[4]{Jan #chi^{2}/ndf = %f}",Chi2J));
	latexLabel.DrawLatex(0.61, 0.70, Form("#color[4]{p-value = %0.6e}",*pValue));

	enregistrementPlots(nomDossier, nomFichier, EndCaps, j, c2);




	c2->Clear();

	////////// chi2 Pulls & Residuals //////////

	string nomDossierPull = "";
        if(EndCaps == 0) nomDossierPull = nomDossier + "EB/";
        if(EndCaps == 1) nomDossierPull = nomDossier + "EE/";

        RooHist* residuals = residHist(Erawframe,(char *)"myhist",(char *)"mycurve", false, nomDossierPull, j);

        //RooHist* pulls = pullHist(Erawframe,"myhist", "mycurve");
        RooHist* pulls = residHist(Erawframe,(char *)"myhist", (char *)"mycurve", true, nomDossierPull, j);

	residuals->Draw("AP");
	nomFichier = "Chi2Residuals";
	
	residuals->GetYaxis()->SetTitle("#chi^{2} Residuals");
        residuals->GetXaxis()->SetLabelFont(42);
        residuals->GetXaxis()->SetTitleFont(42);
        residuals->GetXaxis()->SetLabelSize(0.03);
        residuals->GetXaxis()->SetTitle("s");
        residuals->GetYaxis()->SetLabelFont(42);
        residuals->GetYaxis()->SetTitleOffset(1.24);
        residuals->GetYaxis()->SetTitleFont(42);
        residuals->GetYaxis()->SetLabelSize(0.03);
        //residuals->SetMarkerColor(4);
        //residuals->SetMarkerStyle(21);
        //residuals->SetMarkerSize(0.6);

	TLine *lineResid = new TLine(-0.5,0,0.5,0);
        lineResid->SetLineStyle(3);
        lineResid->Draw();

        TLatex *textResid = new TLatex();

        textResid = new TLatex();
        textResid->SetNDC();
        textResid->SetTextAlign(11);
        textResid->SetTextFont(42);
        textResid->SetTextSizePixels(17);
        textResid->SetTextSize(0.028);
        //if(EndCaps == 0) textResid->DrawLatex(0.80, 0.88, "Barrel");
        //if(EndCaps == 1) textResid->DrawLatex(0.77, 0.88, "End Caps");
	textResid->DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
	if(isMC == 1) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV");
        if(isMC == 0) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV, #int L = 4.93 fb^{-1}");
        textResid->DrawLatex(0.16, 0.80, EndCapsR9Chain->c_str());
        textResid->DrawLatex(0.16, 0.75, tempLegChain->c_str());
	residuals->GetXaxis()->SetLimits(-0.5,0.5);
	
	enregistrementPlots(nomDossier, nomFichier, EndCaps, j, c2);

	c2->Clear();

	pulls->Draw("AP");
	nomFichier = "Chi2Pulls";
	
	pulls->GetYaxis()->SetTitle("#chi^{2} Pulls");
        pulls->GetXaxis()->SetLabelFont(42);
        pulls->GetXaxis()->SetTitleFont(42);
        pulls->GetXaxis()->SetLabelSize(0.03);
        pulls->GetXaxis()->SetTitle("s");
        pulls->GetYaxis()->SetLabelFont(42);
        pulls->GetYaxis()->SetTitleOffset(1.24);
        pulls->GetYaxis()->SetTitleFont(42);
        pulls->GetYaxis()->SetLabelSize(0.03);
        //pulls->SetMarkerColor(4);
        //pulls->SetMarkerStyle(21);
        //pulls->SetMarkerSize(0.6);

	lineResid = new TLine(-0.5,0,0.5,0);
        lineResid->SetLineStyle(3);
        lineResid->Draw();

        textResid = new TLatex();
        textResid->SetNDC();
        textResid->SetTextAlign(11);
        textResid->SetTextFont(42);
        textResid->SetTextSizePixels(17);
        textResid->SetTextSize(0.028);
        //if(EndCaps == 0) textResid->DrawLatex(0.80, 0.88, "Barrel");
        //if(EndCaps == 1) textResid->DrawLatex(0.77, 0.88, "End Caps");	
	textResid->DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
	if(isMC == 1) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV");
        if(isMC == 0) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV, #int L = 4.93 fb^{-1}");
        textResid->DrawLatex(0.16, 0.80, EndCapsR9Chain->c_str());
        textResid->DrawLatex(0.16, 0.75, tempLegChain->c_str());
	pulls->GetXaxis()->SetLimits(-0.5,0.5);	

	enregistrementPlots(nomDossier, nomFichier, EndCaps, j, c2);	

	c2->Clear();

	delete EndCapsR9Chain;
        delete tempLegChain;
	residuals->Delete();
	pulls->Delete();	
	textResid->Delete();
	lineResid->Delete();
	Data_subset->Delete();
	Data->Delete();
	ntplVars->Delete();
	mmg_s->Delete();
        Photon_SC_Eta->Delete();
        Photon_SC_Phi->Delete();
	Photon_r9->Delete();
        isJanLooseMMG->Delete();
        isMultipleCandidate->Delete();
	Photon_Et->Delete();
	Photon_SC_rawEt->Delete();
        Photon_E->Delete();
        Photon_SC_brem->Delete();
        weight_pileUp->Delete();
        Eraw_Landau_m0->Delete();
        Eraw_Landau_sigma->Delete();
	Eraw_Landau->Delete();



}

void RooLandauConvGaussian(double * mean_value, double * mean_error, Double_t * sigmaEff_value, double * sigma_value, double * sigma_value_error, double * sigmaR_value, double * sigmaL_value, double * ChiSquare, TH1D * hh, int j, int EndCaps, TString temp, TTree* Tree_Data, double mean, double rms, TCanvas * c2, RooPlot* Erawframe, TF1 * f, string nomDossier, string nomFichier, double RangeMin, double RangeMax, double * JanChi2, double * DegreesOfFreedom, double * pValue, int r9sup, double MinVar, double MaxVar, string variableX, int isMC)
{

	RooRealVar * mmg_s = new RooRealVar("mmg_s", "s", -0.5, 0.5);
        RooRealVar * Photon_SC_Eta = new RooRealVar("Photon_SC_Eta", "Photon_SC_Eta", -3.0, 3.0, "");
        RooRealVar * Photon_SC_Phi = new RooRealVar("Photon_SC_Phi", "Photon_SC_Phi", -3.5, 3.5, "");
	RooRealVar * Photon_r9 = new RooRealVar("Photon_r9", "Photon_r9", 0.0, 1.0, "");
        RooRealVar * isJanLooseMMG = new RooRealVar("isJanLooseMMG", "isJanLooseMMG", 0.0, 1.0, "");
        RooRealVar * isMultipleCandidate = new RooRealVar("isMultipleCandidate", "isMultipleCandidate", 0.0, 1.0, "");
	RooRealVar * Photon_Et = new RooRealVar("Photon_Et", "Photon_Et", 0.0, 250.0, "");	
	RooRealVar * Photon_SC_rawEt = new RooRealVar("Photon_SC_rawEt", "Photon_SC_rawEt", 0.0, 250.0, ""); 
        RooRealVar * Photon_E = new RooRealVar("Photon_E", "Photon_E", 0.0, 1000.0, ""); 
        RooRealVar * Photon_SC_brem = new RooRealVar("Photon_SC_brem", "Photon_SC_brem", 0.0, 15.0, ""); 
        RooRealVar * weight_pileUp = new RooRealVar("weight_pileUp", "weight_pileUp", 0.0, 10000000); 

	RooArgSet *ntplVars = new RooArgSet(*mmg_s, *Photon_SC_Eta, *Photon_SC_Phi, *Photon_r9, *Photon_Et, *Photon_SC_rawEt, *Photon_E, *Photon_SC_brem, *isJanLooseMMG);

        ntplVars->add(*isMultipleCandidate);
        ntplVars->add(*weight_pileUp);



        //RooArgSet *ntplVars = new RooArgSet(*mmg_s, *Photon_SC_Eta, *Photon_SC_Phi, *Photon_r9, *Photon_Et, *Photon_SC_rawEt, *Photon_E, *Photon_SC_brem, *isJanLooseMMG, *isMultipleCandidate, *weight_pileUp);



        //RooArgSet *ntplVars = new RooArgSet(*mmg_s, *Photon_SC_Eta, *Photon_SC_Phi, *Photon_r9, *Photon_Et, *isJanLooseMMG, *isMultipleCandidate, *weight_pileUp);
	//RooArgSet ntplVars(mmg_s, Photon_SC_Eta, Photon_r9, Photon_Et, isJanLooseMMG, Photon_Et);
        //RooDataSet *Data = new RooDataSet("Data", "Data", Tree_Data, *ntplVars); //non weighté
	RooDataSet *Data = new RooDataSet("Data", "Data", Tree_Data, *ntplVars, "", "weight_pileUp"); //weighté
	//RooDataSet Data("Data", "Data", Tree_Data, ntplVars);


        RooDataSet* Data_subset = (RooDataSet*)Data->reduce(*ntplVars, temp);
	//RooDataSet Data_subset = (RooDataSet)Data.reduce(ntplVars, temp);
        //RooPlot* Erawframe = mmg_s.frame();
        
/*	
	Erawframe = mmg_s->frame();
	Data_subset->plotOn(Erawframe,Name("myhist"),Binning(300));
        Erawframe->Draw();
*/

	RooRealVar * Eraw_Landau_m0 = new RooRealVar("Eraw_Landau_m0", "Landau #Delta m_{0}", mean, mean - rms, mean + rms, "GeV");
        RooRealVar * Eraw_Landau_sigma = new RooRealVar("Eraw_Landau_sigma", "Landau ", rms , 0.0,1.0, "GeV");

	// --- Build Landau PDF ---
        RooLandau * Eraw_Landau = new RooLandau("Eraw_Landau","mmg_s_Landau", *mmg_s, *Eraw_Landau_m0, *Eraw_Landau_sigma);


	/*
        RooRealVar * sigmean = new RooRealVar("sigmean","B^{#pm} mass",0.0,0.0,0.0);
        RooRealVar * sigwidth = new RooRealVar("sigwidth","B^{#pm} width",rms, 0.0,1.0);


        // --- Build Gaussian PDF ---
        RooGaussian * Eraw_Gaussian = new RooGaussian("Eraw_Gaussian","mmg_s_Gaussian",*mmg_s,*sigmean,*sigwidth) ;


        // CONVOLUTION
        RooFFTConvPdf * model = new RooFFTConvPdf("model", "model", *mmg_s, *Eraw_Landau, *Eraw_Gaussian);
*/

        int fewBins = 1; 
        Double_t Chi2J;

        Erawframe = mmg_s->frame();
        Data_subset->plotOn(Erawframe,Name("myhist"));
        Eraw_Landau->fitTo(*Data_subset, Range(RangeMin, RangeMax));
        Eraw_Landau->plotOn(Erawframe,Name("mycurve"));
        Erawframe->Draw();

        RooRealVar * sigmean = new RooRealVar("sigmean","B^{#pm} mass",0.0,0.0,0.0);
        RooRealVar * sigwidth = new RooRealVar("sigwidth","B^{#pm} width",(Eraw_Landau_sigma->getVal() / 4.0), 0.0,(Eraw_Landau_sigma->getVal() / 1.5));

        // --- Build Gaussian PDF ---
        RooGaussian * Eraw_Gaussian = new RooGaussian("Eraw_Gaussian","mmg_s_Gaussian",*mmg_s,*sigmean,*sigwidth) ;


        // CONVOLUTION
        RooFFTConvPdf * model = new RooFFTConvPdf("model", "model", *mmg_s, *Eraw_Landau, *Eraw_Gaussian);




	for(int RightBinning = 600; RightBinning > 10; RightBinning -= 10)
	{
		cout<<endl<<endl<<endl<<"RightBinning = "<<RightBinning<<endl<<endl<<endl;
		Erawframe->Clear();
		Data_subset = (RooDataSet*)Data->reduce(*ntplVars, temp);		
		Eraw_Landau = new RooLandau("Eraw_Landau","mmg_s_Landau", *mmg_s, *Eraw_Landau_m0, *Eraw_Landau_sigma);
		Eraw_Gaussian = new RooGaussian("Eraw_Gaussian","mmg_s_Gaussian",*mmg_s,*sigmean,*sigwidth);
		model = new RooFFTConvPdf("model", "model", *mmg_s, *Eraw_Landau, *Eraw_Gaussian); 

		Erawframe = mmg_s->frame();
		Data_subset->plotOn(Erawframe,Name("myhist"),Binning(RightBinning));
		//Erawframe->Draw();
        	//model.fitTo(*Data_subset, Range(0.6,1.2));
		//model.fitTo(Data_subset, Range(0.5,1.15)); //ATTENTION NE PAS SUPPRIMER!!!!
		//model.fitTo(Data_subset, Range(mean - 5 * rms,mean + 5 * rms));
		//double maxDistri = hh->GetMaximumBin() * 0.01;
		//model.fitTo(*Data_subset, Range(mean - 0.35 * rms,mean + 0.35 * rms));
		//model.fitTo(*Data_subset, Range(maxDistri - 0.6 * rms,maxDistri + 0.6 * rms));
		model->fitTo(*Data_subset, Range(RangeMin, RangeMax));
		//RooArgSet* model_param = model.getVariables();
		//model_param->Print("v");
		//model.plotOn(Erawframe);
        	model->plotOn(Erawframe,Name("mycurve"));
		Erawframe->Draw();


		Chi2J = chiSquare(Erawframe,(char *)"mycurve",(char *)"myhist",4, JanChi2, DegreesOfFreedom, pValue, &fewBins);
		if(fewBins == 0) break; 


	}

	// Create function that represents derivative
        //RooAbsReal * deriv = (RooAbsReal*) model->derivative(*mmg_s,1) ;
        RooDerivative * deriv = model->derivative(*mmg_s,1) ;

        // Find point where derivative is zero
        Double_t x_max = deriv->findRoot(*mmg_s,RangeMin, RangeMax,0) ;


	//RooDataHist data("data", "dataset with mmg_s",*mmg_s,hh);
	//TF1 * f = model.asTF( RooArgList(mmg_s) );
	f = model->asTF( RooArgList(*mmg_s) );
	//ftest = model.asTF( RooArgList(*mmg_s) );
	//cout<<endl<<"ftest->GetXmax() = "<<ftest->GetXmax()<<endl;
	//cout<<endl<<"ftest->GetMaximum = "<<ftest->GetMaximum(0.0, 2.0,1.E-10,100,false)<<endl;
	*mean_value  = x_max;
        *mean_error  = Eraw_Landau_m0->getError();
        *sigmaEff_value = effSigma(hh);
        *sigma_value = Eraw_Landau_sigma->getVal();
	*sigma_value_error = Eraw_Landau_sigma->getError();
        *sigmaR_value = SigmaR(f, 0.0, 2.0);
        *sigmaL_value = SigmaL(f, 0.0, 2.0);
        *ChiSquare = Erawframe->chiSquare();

 
        double entries = hh->GetEntries();
        
        int entriesInt = (int) entries;


        TLatex latexLabel;
        latexLabel.SetTextSize(0.026);
        latexLabel.SetNDC();
        latexLabel.DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
	if(isMC == 1) latexLabel.DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV");
        if(isMC == 0) latexLabel.DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV, #int L = 4.93 fb^{-1}");
        //if(EndCaps == 0) latexLabel.DrawLatex(0.16, 0.85, "Barrel");
        //if(EndCaps == 1) latexLabel.DrawLatex(0.16, 0.85, "EndCaps");   
        string * EndCapsR9Chain(0);
        EndCapsR9Chain = new string;
        if(EndCaps == 0) *EndCapsR9Chain = "Barrel, ";
        if(EndCaps == 1) *EndCapsR9Chain = "Endcaps, ";
        if(r9sup == 0) *EndCapsR9Chain += "Low r9";
        if(r9sup == 1) *EndCapsR9Chain += "High r9";
        if(r9sup == 2) *EndCapsR9Chain += "All r9";
	string * tempLegChain(0);
        tempLegChain = new string;
        *tempLegChain = DoubleToString(MinVar) + " <" + variableX + "< "  + DoubleToString(MaxVar);     
        latexLabel.DrawLatex(0.16, 0.80, EndCapsR9Chain->c_str());
        latexLabel.DrawLatex(0.16, 0.75, tempLegChain->c_str());

        latexLabel.DrawLatex(0.16, 0.70, Form("Entries = %d",entriesInt));
        cout<<endl<<"entries = "<<entries<<endl;
        cout<<endl<<"Mean = "<<hh->GetMean()<<endl;
        latexLabel.DrawLatex(0.16, 0.65, Form("Mean = %f +/- %f",hh->GetMean(), hh->GetMeanError()));
        latexLabel.DrawLatex(0.16, 0.60, Form("RMS = %f +/- %f",hh->GetRMS(), hh->GetRMSError()));
	latexLabel.DrawLatex(0.61, 0.90, Form("#color[4]{m_{0 Landau} = %f +/- %f}",Eraw_Landau_m0->getVal(), Eraw_Landau_m0->getError()));
        latexLabel.DrawLatex(0.61, 0.85, Form("#color[4]{#sigma_{Landau} = %f +/- %f}",Eraw_Landau_sigma->getVal(), Eraw_Landau_sigma->getError()));
	latexLabel.DrawLatex(0.61, 0.80, Form("#color[4]{m_{0 Gaus} = %f +/- %f}",sigmean->getVal(), sigmean->getError()));
        latexLabel.DrawLatex(0.61, 0.75, Form("#color[4]{#sigma_{Gaus} = %f +/- %f}",sigwidth->getVal(), sigwidth->getError()));
	latexLabel.DrawLatex(0.61, 0.70, Form("#color[4]{m_{0 Conv} = %f}",*mean_value));
        latexLabel.DrawLatex(0.61, 0.65, Form("#color[4]{Default #chi^{2}/ndf = %f}",Erawframe->chiSquare()));
        latexLabel.DrawLatex(0.61, 0.60, Form("#color[4]{Jan #chi^{2}/ndf = %f}",Chi2J));
        latexLabel.DrawLatex(0.61, 0.55, Form("#color[4]{p-value = %0.6e}",*pValue));

	enregistrementPlots(nomDossier, nomFichier, EndCaps, j, c2);




	c2->Clear();

	////////// chi2 Pulls & Residuals //////////

	string nomDossierPull = "";
        if(EndCaps == 0) nomDossierPull = nomDossier + "EB/";
        if(EndCaps == 1) nomDossierPull = nomDossier + "EE/";

        RooHist* residuals = residHist(Erawframe,(char *)"myhist",(char *)"mycurve", false, nomDossierPull, j);

        //RooHist* pulls = pullHist(Erawframe,"myhist", "mycurve");
        RooHist* pulls = residHist(Erawframe,(char *)"myhist", (char *)"mycurve", true, nomDossierPull, j);

	residuals->Draw("AP");
	nomFichier = "Chi2Residuals";
	
	residuals->GetYaxis()->SetTitle("#chi^{2} Residuals");
        residuals->GetXaxis()->SetLabelFont(42);
        residuals->GetXaxis()->SetTitleFont(42);
        residuals->GetXaxis()->SetLabelSize(0.03);
        residuals->GetXaxis()->SetTitle("s");
        residuals->GetYaxis()->SetLabelFont(42);
        residuals->GetYaxis()->SetTitleOffset(1.24);
        residuals->GetYaxis()->SetTitleFont(42);
        residuals->GetYaxis()->SetLabelSize(0.03);
        //residuals->SetMarkerColor(4);
        //residuals->SetMarkerStyle(21);
        //residuals->SetMarkerSize(0.6);

	TLine *lineResid = new TLine(-0.5,0,0.5,0);
        lineResid->SetLineStyle(3);
        lineResid->Draw();

        TLatex *textResid = new TLatex();

        textResid = new TLatex();
        textResid->SetNDC();
        textResid->SetTextAlign(11);
        textResid->SetTextFont(42);
        textResid->SetTextSizePixels(17);
        textResid->SetTextSize(0.028);
        //if(EndCaps == 0) textResid->DrawLatex(0.80, 0.88, "Barrel");
        //if(EndCaps == 1) textResid->DrawLatex(0.77, 0.88, "End Caps");
	textResid->DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
	if(isMC == 1) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV");
        if(isMC == 0) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV, #int L = 4.93 fb^{-1}");
        textResid->DrawLatex(0.16, 0.80, EndCapsR9Chain->c_str());
        textResid->DrawLatex(0.16, 0.75, tempLegChain->c_str());
	residuals->GetXaxis()->SetLimits(-0.5,0.5);
	
	enregistrementPlots(nomDossier, nomFichier, EndCaps, j, c2);

	c2->Clear();

	pulls->Draw("AP");
	nomFichier = "Chi2Pulls";
	
	pulls->GetYaxis()->SetTitle("#chi^{2} Pulls");
        pulls->GetXaxis()->SetLabelFont(42);
        pulls->GetXaxis()->SetTitleFont(42);
        pulls->GetXaxis()->SetLabelSize(0.03);
        pulls->GetXaxis()->SetTitle("s");
        pulls->GetYaxis()->SetLabelFont(42);
        pulls->GetYaxis()->SetTitleOffset(1.24);
        pulls->GetYaxis()->SetTitleFont(42);
        pulls->GetYaxis()->SetLabelSize(0.03);
        //pulls->SetMarkerColor(4);
        //pulls->SetMarkerStyle(21);
        //pulls->SetMarkerSize(0.6);

	lineResid = new TLine(-0.5,0,0.5,0);
        lineResid->SetLineStyle(3);
        lineResid->Draw();

        textResid = new TLatex();
        textResid->SetNDC();
        textResid->SetTextAlign(11);
        textResid->SetTextFont(42);
        textResid->SetTextSizePixels(17);
        textResid->SetTextSize(0.028);
        //if(EndCaps == 0) textResid->DrawLatex(0.80, 0.88, "Barrel");
        //if(EndCaps == 1) textResid->DrawLatex(0.77, 0.88, "End Caps");	
	textResid->DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
	if(isMC == 1) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV");
        if(isMC == 0) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV, #int L = 4.93 fb^{-1}");
        textResid->DrawLatex(0.16, 0.80, EndCapsR9Chain->c_str());
        textResid->DrawLatex(0.16, 0.75, tempLegChain->c_str());
	pulls->GetXaxis()->SetLimits(-0.5,0.5);	

	enregistrementPlots(nomDossier, nomFichier, EndCaps, j, c2);	

	c2->Clear();

	delete EndCapsR9Chain;
        delete tempLegChain;
	residuals->Delete();
	pulls->Delete();	
	textResid->Delete();
	lineResid->Delete();
	Data_subset->Delete();
	Data->Delete();
	ntplVars->Delete();
	mmg_s->Delete();
        Photon_SC_Eta->Delete();
        Photon_SC_Phi->Delete();
	Photon_r9->Delete();
        isJanLooseMMG->Delete();
        isMultipleCandidate->Delete();
	Photon_Et->Delete();
	Photon_SC_rawEt->Delete();
        Photon_E->Delete();
        Photon_SC_brem->Delete();
        weight_pileUp->Delete();
        Eraw_Landau_m0->Delete();
        Eraw_Landau_sigma->Delete();
	Eraw_Landau->Delete();
	sigmean->Delete();
	sigwidth->Delete();
	Eraw_Gaussian->Delete();
	model->Delete();



}


void RooKernel(double * mean_value, double * mean_error, Double_t * sigmaEff_value, double * sigma_value, double * sigma_value_error, double * sigmaR_value, double * sigmaL_value, double * ChiSquare, double * minLogLikelihood, double * differenceLogLikelihood, TH1D * hh, int j, int EndCaps, TString temp, TTree* Tree_Data, double mean, double rms, TCanvas * c2, RooPlot* Erawframe, TF1 * f, string nomDossier, string nomFichier, double RangeMin, double RangeMax, double * JanChi2, double * DegreesOfFreedom, double * pValue, int r9sup, double MinVar, double MaxVar, string variableX, int isMC)
{
	RooRealVar * mmg_s = new RooRealVar("mmg_s", "s", -0.5, 0.5);
        RooRealVar * Photon_SC_Eta = new RooRealVar("Photon_SC_Eta", "Photon_SC_Eta", -3.0, 3.0, "");
        RooRealVar * Photon_SC_Phi = new RooRealVar("Photon_SC_Phi", "Photon_SC_Phi", -3.5, 3.5, "");
        RooRealVar * Photon_r9 = new RooRealVar("Photon_r9", "Photon_r9", 0.0, 1.0, "");
        RooRealVar * isJanLooseMMG = new RooRealVar("isJanLooseMMG", "isJanLooseMMG", 0.0, 1.0, "");
	RooRealVar * isMultipleCandidate = new RooRealVar("isMultipleCandidate", "isMultipleCandidate", 0.0, 1.0, "");
        RooRealVar * Photon_Et = new RooRealVar("Photon_Et", "Photon_Et", 0.0, 250.0, "");
        RooRealVar * Photon_SC_rawEt = new RooRealVar("Photon_SC_rawEt", "Photon_SC_rawEt", 0.0, 250.0, "");
        RooRealVar * Photon_E = new RooRealVar("Photon_E", "Photon_E", 0.0, 1000.0, "");
        RooRealVar * Photon_SC_brem = new RooRealVar("Photon_SC_brem", "Photon_SC_brem", 0.0, 15.0, "");
        RooRealVar * weight_pileUp = new RooRealVar("weight_pileUp", "weight_pileUp", 0.0, 10000000);

        RooArgSet *ntplVars = new RooArgSet(*mmg_s, *Photon_SC_Eta, *Photon_SC_Phi, *Photon_r9, *Photon_Et, *Photon_SC_rawEt, *Photon_E, *Photon_SC_brem, *isJanLooseMMG);



        ntplVars->add(*isMultipleCandidate);
        ntplVars->add(*weight_pileUp);


        RooDataSet *Data = new RooDataSet("Data", "Data", Tree_Data, *ntplVars, "", "weight_pileUp"); //weighté


        RooDataSet* Data_subset = (RooDataSet*)Data->reduce(*ntplVars, temp);

	RooKeysPdf * K = new RooKeysPdf("K","K",*mmg_s,*Data_subset);
	Erawframe = mmg_s->frame();
	
	//Data_subset->plotOn(Erawframe,Name("myhist"),Binning(60));
	Data_subset->plotOn(Erawframe,Name("myhist"));
	RooFitResult* r = K->fitTo(*Data_subset,Minos(true),Save());
        K->plotOn(Erawframe,Name("mycurve"));
        Erawframe->Draw();


	double entries = hh->GetEntries();

        int entriesInt = (int) entries;


        TLatex latexLabel;
        latexLabel.SetTextSize(0.026);
        latexLabel.SetNDC();
        latexLabel.DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
        if(isMC == 1) latexLabel.DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV");
        if(isMC == 0) latexLabel.DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV, #int L = 4.93 fb^{-1}");
        //if(EndCaps == 0) latexLabel.DrawLatex(0.16, 0.85, "Barrel");
        //if(EndCaps == 1) latexLabel.DrawLatex(0.16, 0.85, "EndCaps");   
        string * EndCapsR9Chain(0);
        EndCapsR9Chain = new string;
        if(EndCaps == 0) *EndCapsR9Chain = "Barrel, ";
        if(EndCaps == 1) *EndCapsR9Chain = "Endcaps, ";
        if(r9sup == 0) *EndCapsR9Chain += "Low r9";
        if(r9sup == 1) *EndCapsR9Chain += "High r9";
        if(r9sup == 2) *EndCapsR9Chain += "All r9";
        string * tempLegChain(0);
        tempLegChain = new string;
        *tempLegChain = DoubleToString(MinVar) + " <" + variableX + "< "  + DoubleToString(MaxVar);
        latexLabel.DrawLatex(0.16, 0.80, EndCapsR9Chain->c_str());
        latexLabel.DrawLatex(0.16, 0.75, tempLegChain->c_str());

        latexLabel.DrawLatex(0.16, 0.70, Form("Entries = %d",entriesInt));
        cout<<endl<<"entries = "<<entries<<endl;
        cout<<endl<<"Mean = "<<hh->GetMean()<<endl;
        latexLabel.DrawLatex(0.16, 0.65, Form("Mean = %f +/- %f",hh->GetMean(), hh->GetMeanError()));
        latexLabel.DrawLatex(0.16, 0.60, Form("RMS = %f +/- %f",hh->GetRMS(), hh->GetRMSError()));
        /*latexLabel.DrawLatex(0.61, 0.90, Form("#color[4]{m_{0 Landau} = %f +/- %f}",Eraw_Landau_m0->getVal(), Eraw_Landau_m0->getError()));
        latexLabel.DrawLatex(0.61, 0.85, Form("#color[4]{#sigma_{Landau} = %f +/- %f}",Eraw_Landau_sigma->getVal(), Eraw_Landau_sigma->getError()));
        latexLabel.DrawLatex(0.61, 0.80, Form("#color[4]{m_{0 Gaus} = %f +/- %f}",sigmean->getVal(), sigmean->getError()));
        latexLabel.DrawLatex(0.61, 0.75, Form("#color[4]{#sigma_{Gaus} = %f +/- %f}",sigwidth->getVal(), sigwidth->getError()));
        latexLabel.DrawLatex(0.61, 0.70, Form("#color[4]{m_{0 Conv} = %f}",*mean_value));
        latexLabel.DrawLatex(0.61, 0.65, Form("#color[4]{Default #chi^{2}/ndf = %f}",Erawframe->chiSquare()));
        latexLabel.DrawLatex(0.61, 0.60, Form("#color[4]{Jan #chi^{2}/ndf = %f}",Chi2J));
        latexLabel.DrawLatex(0.61, 0.55, Form("#color[4]{p-value = %0.6e}",*pValue));
	*/

	enregistrementPlots(nomDossier, nomFichier, EndCaps, j, c2);

	c2->Clear();

	delete EndCapsR9Chain;
        delete tempLegChain;
        //residuals->Delete();
        //pulls->Delete();     
        //textResid->Delete();
        //lineResid->Delete();
        Data_subset->Delete();
        Data->Delete();
        ntplVars->Delete();
        mmg_s->Delete();
        Photon_SC_Eta->Delete();
        Photon_SC_Phi->Delete();
        Photon_r9->Delete();
        isJanLooseMMG->Delete();
	isMultipleCandidate->Delete();
        Photon_Et->Delete();
        Photon_SC_rawEt->Delete();
        Photon_E->Delete();
        Photon_SC_brem->Delete();
        weight_pileUp->Delete();
        K->Delete();

}

void RooVoigtian2(double * mean_value, double * mean_error, Double_t * sigmaEff_value, double * sigma_value, double * sigma_value_error, double * sigmaR_value, double * sigmaL_value, double * ChiSquare, double * minLogLikelihood, double * differenceLogLikelihood, TH1D * hh, int j, int EndCaps, TString temp, TTree* Tree_Data, double mean, double rms, TCanvas * c2, RooPlot* Erawframe, TF1 * f, string nomDossier, string nomFichier, double RangeMin, double RangeMax, double * JanChi2, double * DegreesOfFreedom, double * pValue, int r9sup, double MinVar, double MaxVar, string variableX, int isMC)
{
	//RooRealVar * mmg_s = new RooRealVar("mmg_s", "s", -0.5, 0.5);
        RooRealVar * mmg_s = new RooRealVar("mmg_s", "s", RangeMin, RangeMax);
	RooRealVar * Photon_SC_Eta = new RooRealVar("Photon_SC_Eta", "Photon_SC_Eta", -3.0, 3.0, "");
        RooRealVar * Photon_SC_Phi = new RooRealVar("Photon_SC_Phi", "Photon_SC_Phi", -3.5, 3.5, "");
        RooRealVar * Photon_r9 = new RooRealVar("Photon_r9", "Photon_r9", 0.0, 1.0, "");
        RooRealVar * isJanLooseMMG = new RooRealVar("isJanLooseMMG", "isJanLooseMMG", 0.0, 1.0, "");
        RooRealVar * isMultipleCandidate = new RooRealVar("isMultipleCandidate", "isMultipleCandidate", 0.0, 1.0, "");
        RooRealVar * Photon_Et = new RooRealVar("Photon_Et", "Photon_Et", 0.0, 250.0, "");
        RooRealVar * Photon_SC_rawEt = new RooRealVar("Photon_SC_rawEt", "Photon_SC_rawEt", 0.0, 250.0, "");
        RooRealVar * Photon_E = new RooRealVar("Photon_E", "Photon_E", 0.0, 1000.0, "");
        RooRealVar * Photon_SC_brem = new RooRealVar("Photon_SC_brem", "Photon_SC_brem", 0.0, 15.0, "");
        RooRealVar * weight_pileUp = new RooRealVar("weight_pileUp", "weight_pileUp", 0.0, 100);

        RooArgSet *ntplVars = new RooArgSet(*mmg_s, *Photon_SC_Eta, *Photon_SC_Phi, *Photon_r9, *Photon_Et, *Photon_SC_rawEt, *Photon_E, *Photon_SC_brem, *isJanLooseMMG);



        ntplVars->add(*isMultipleCandidate);
        ntplVars->add(*weight_pileUp);


        //RooDataSet *Data = new RooDataSet("Data", "Data", Tree_Data, *ntplVars, "", "weight_pileUp"); //weighté
	//RooDataSet *Data = new RooDataSet("Data", "Data", Tree_Data, *ntplVars, "", "");
	RooDataSet *Data_subset = new RooDataSet("Data_subset", "Data_subset", Tree_Data, *ntplVars, temp, "weight_pileUp");


	RooRealVar * meanV = new RooRealVar("meanV","meanV",0.0,-0.1,0.1);
        RooRealVar * sigma = new RooRealVar("sigma","sigma",0.5,0.0,1.0);
        RooRealVar * width = new RooRealVar("width","width",0.5,0.0,1.0);


	// --- Build Bifurcated Gaussian PDF ---
	RooVoigtian * Voigtian = new RooVoigtian("Voigtian","Voigtian",*mmg_s,*meanV,*sigma,*width);


	int fewBins = 1;
        Double_t Chi2J;

	RooFitResult *res;
	Double_t minNll;
		

        //for(int RightBinning = 600; RightBinning > 10; RightBinning -= 10)
	//int RightBinning = 50;
        //RightBinning = (int) ((RangeMax - RangeMin) / 0.02);
	RooBinning b(-17.0,17.0) ;
        b.addUniform(1700,-17.0,17.0) ;
	for(int i = 0; i<4; i++)
	{
                Erawframe->Clear();
                //Data_subset = (RooDataSet*)Data->reduce(*ntplVars, temp);
		Voigtian = new RooVoigtian("Voigtian","Voigtian",*mmg_s,*meanV,*sigma,*width);

                //Erawframe = mmg_s->frame();
                Erawframe = mmg_s->frame(-0.5,0.5);
		Data_subset->plotOn(Erawframe,Name("myhist"),Binning(b),DataError(RooAbsData::SumW2));
                //Voigtian->fitTo(*Data_subset, Range(RangeMin, RangeMax),Minos(true));
                //res = Voigtian->fitTo(*Data_subset, Range(RangeMin, RangeMax),Minos(true),Save(),SumW2Error(kFALSE));
		res = Voigtian->fitTo(*Data_subset, Range(RangeMin, RangeMax),Save(),SumW2Error(kTRUE));
		//res = Voigtian->fitTo(*Data_subset,Save(),SumW2Error(kTRUE));
		res->Print();
		minNll = res->minNll();
		Voigtian->plotOn(Erawframe,Name("mycurve"));
		Erawframe->Draw();


                Chi2J = chiSquare(Erawframe,(char *)"mycurve",(char *)"myhist",3, JanChi2, DegreesOfFreedom, pValue, &fewBins);
                if(fewBins == 0) break; 


        }

	c2->Clear();
        Erawframe->Clear();

        ///// Draw a nice plot /////

        RooRealVar * mmg_s2 = new RooRealVar("mmg_s", "s", -0.5, 0.5);
        RooArgSet *ntplVars2 = new RooArgSet(*mmg_s2, *Photon_SC_Eta, *Photon_SC_Phi, *Photon_r9, *Photon_Et, *Photon_SC_rawEt, *Photon_E, *Photon_SC_brem, *isJanLooseMMG);
        ntplVars2->add(*isMultipleCandidate);
        ntplVars2->add(*weight_pileUp); 
        RooDataSet *Data_subset2 = new RooDataSet("Data_subset2", "Data_subset2", Tree_Data, *ntplVars2, temp, "weight_pileUp");
        Erawframe = mmg_s2->frame();
        Data_subset2->plotOn(Erawframe,Name("myhist"),Binning(b));
        Voigtian->plotOn(Erawframe,Name("mycurve"),Range(RangeMin, RangeMax));
        Erawframe->Draw();

        ///// /////


	cout << endl << "ICI !!!!!!!!!!!!!!!!!! "<<Data_subset->weight() << endl ;

        f = Voigtian->asTF( RooArgList(*mmg_s) );
	*mean_value  = meanV->getVal();
        *mean_error  = meanV->getError();
        *sigmaEff_value = effSigma(hh);
        //*sigma_value = sigwidth.getVal();
        *sigmaR_value = SigmaR(f, 0.0, 2.0);
        *sigmaL_value = SigmaL(f, 0.0, 2.0);
        *ChiSquare = Erawframe->chiSquare();
	*minLogLikelihood = minNll;


        double entries = hh->GetEntries();

        int entriesInt = (int) entries;

	
	double fxmax = f->GetMaximumX(0.8,1.2,1.E-10,100,false);
        //cout<<"////////// ---- FXMAX BG = "<<f->GetMaximumX(0.8,1.2,1.E-10,100,false)<<" ---- //////////"<<endl;

	TLatex latexLabel;
        latexLabel.SetTextSize(0.030);
        latexLabel.SetNDC();
	latexLabel.DrawLatex(0.13, 0.96, "CMS Preliminary 2011, #sqrt{s} = 7 TeV");
        if(isMC == 1) latexLabel.DrawLatex(0.17, 0.88, "Simulation");
        if(isMC == 0) latexLabel.DrawLatex(0.17, 0.88, "Data, #int L = 4.93 fb^{-1}"); //lumi a changer !!!

	if(EndCaps == 0) latexLabel.DrawLatex(0.17, 0.83,"ECAL Barrel");
        if(EndCaps == 1) latexLabel.DrawLatex(0.17, 0.83,"ECAL Endcaps");
        if(r9sup == 0 && EndCaps == 0) latexLabel.DrawLatex(0.17, 0.78,"E_{T}^{#gamma} > 25 GeV, r9 < 0,94");
        if(r9sup == 0 && EndCaps == 1) latexLabel.DrawLatex(0.17, 0.78,"E_{T}^{#gamma} > 25 GeV, r9 < 0,95");
        if(r9sup == 1 && EndCaps == 0) latexLabel.DrawLatex(0.17, 0.78,"E_{T}^{#gamma} > 25 GeV, r9 > 0,94");
        if(r9sup == 1 && EndCaps == 1) latexLabel.DrawLatex(0.17, 0.78,"E_{T}^{#gamma} > 25 GeV, r9 > 0,95");
        if(r9sup == 2) latexLabel.DrawLatex(0.17, 0.78,"E_{T}^{#gamma} > 25 GeV, All r9");
	
	//latexLabel.DrawLatex(0.17, 0.73, Form("Entries = %d",entriesInt));


	latexLabel.DrawLatex(0.65, 0.88, Form("#color[1]{#mu(E_{#gamma}) = %f #pm %f}",meanV->getVal(), meanV->getError()));
	latexLabel.DrawLatex(0.65, 0.83, Form("#color[1]{#sigma = %f #pm %f}",sigma->getVal(), sigma->getError()));
	latexLabel.DrawLatex(0.65, 0.78, Form("#color[1]{width = %f #pm %f}",width->getVal(), width->getError()));	


	enregistrementPlots(nomDossier, nomFichier, EndCaps, j, c2);

	c2->Clear();
/*
	////////// RooMCStudy //////////

	
	RooMCStudy* mcs = new RooMCStudy(*Voigtian,*mmg_s,FitModel(*Voigtian),Silence(),Binned()) ;
        RooChi2MCSModule chi2mod ;
        mcs->addModule(chi2mod) ;

        mcs->generateAndFit(2000,entriesInt) ;

        //RooPlot* frame = mcs->plotNLL(Bins(50));
        RooPlot* frame = mcs->plotNLL(Name("myhistTest"));
	frame->Draw();

	//RooRealVar * VArTest = (RooRealVar*) frame->getPlotVar(); //artefact

	
	RooHist* histTest = frame->getHist("myhistTest");
	cout<<endl<<"histTest->GetMean() = "<<histTest->GetMean()<<endl;

	*differenceLogLikelihood = fabs(histTest->GetMean() - minNll);
*/
	TLatex *textResid = new TLatex();
/*
        textResid = new TLatex();
        textResid->SetNDC();
        textResid->SetTextAlign(11);
        textResid->SetTextFont(42);
        textResid->SetTextSizePixels(17);
        textResid->SetTextSize(0.028);
        //if(EndCaps == 0) textResid->DrawLatex(0.80, 0.88, "Barrel");
        //if(EndCaps == 1) textResid->DrawLatex(0.77, 0.88, "End Caps");
        textResid->DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
        if(isMC == 1) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV");
        if(isMC == 0) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV, #int L = 4.93 fb^{-1}");
        //textResid->DrawLatex(0.16, 0.80, EndCapsR9Chain->c_str());
        //textResid->DrawLatex(0.16, 0.75, tempLegChain->c_str());

	textResid->DrawLatex(0.61, 0.90, Form("#color[4]{Fit_minll = %f}",minNll));
	textResid->DrawLatex(0.61, 0.85, Form("#color[4]{Distrib_mean = %f}",histTest->GetMean()));
	textResid->DrawLatex(0.61, 0.80, Form("#color[3]{Difference = %f}",*differenceLogLikelihood));
	
	//RooAbsRealLValue VArTest = frame->getPlotVar();


	enregistrementPlots(nomDossier, "LogLikelihood", EndCaps, j, c2);

	c2->Clear();

*/
	////////// chi2 Pulls & Residuals //////////

	string nomDossierPull = "";
        if(EndCaps == 0) nomDossierPull = nomDossier + "EB/";
        if(EndCaps == 1) nomDossierPull = nomDossier + "EE/";


        RooHist* residuals = residHist(Erawframe,(char *)"myhist",(char *)"mycurve", false, nomDossierPull, j);
	

        //RooHist* pulls = pullHist(Erawframe,"myhist", "mycurve");
        RooHist* pulls = residHist(Erawframe,(char *)"myhist", (char *)"mycurve", true, nomDossierPull, j);


	residuals->Draw("AP");
	nomFichier = "Chi2Residuals";
	
	residuals->GetYaxis()->SetTitle("#chi^{2} Residuals");
        residuals->GetXaxis()->SetLabelFont(42);
        residuals->GetXaxis()->SetTitleFont(42);
        residuals->GetXaxis()->SetLabelSize(0.03);
        residuals->GetXaxis()->SetTitle("s");
        residuals->GetYaxis()->SetLabelFont(42);
        residuals->GetYaxis()->SetTitleOffset(1.24);
        residuals->GetYaxis()->SetTitleFont(42);
        residuals->GetYaxis()->SetLabelSize(0.03);
        //residuals->SetMarkerColor(4);
        //residuals->SetMarkerStyle(21);
        //residuals->SetMarkerSize(0.6);

	TLine *lineResid = new TLine(-0.5,0,0.5,0);
        lineResid->SetLineStyle(3);
        lineResid->Draw();


        //TLatex *textResid = new TLatex();

        textResid = new TLatex();
        textResid->SetNDC();
        textResid->SetTextAlign(11);
        textResid->SetTextFont(42);
        textResid->SetTextSizePixels(17);
        textResid->SetTextSize(0.028);
        //if(EndCaps == 0) textResid->DrawLatex(0.80, 0.88, "Barrel");
        //if(EndCaps == 1) textResid->DrawLatex(0.77, 0.88, "End Caps");
	textResid->DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
	if(isMC == 1) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV");
        if(isMC == 0) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV, #int L = 4.93 fb^{-1}");
        //textResid->DrawLatex(0.16, 0.80, EndCapsR9Chain->c_str());
        //textResid->DrawLatex(0.16, 0.75, tempLegChain->c_str());
	residuals->GetXaxis()->SetLimits(-0.5,0.5);
	
	enregistrementPlots(nomDossier, nomFichier, EndCaps, j, c2);

	c2->Clear();

	pulls->Draw("AP");
	nomFichier = "Chi2Pulls";
	
	pulls->GetYaxis()->SetTitle("#chi^{2} Pulls");
        pulls->GetXaxis()->SetLabelFont(42);
        pulls->GetXaxis()->SetTitleFont(42);
        pulls->GetXaxis()->SetLabelSize(0.03);
        pulls->GetXaxis()->SetTitle("s");
        pulls->GetYaxis()->SetLabelFont(42);
        pulls->GetYaxis()->SetTitleOffset(1.24);
        pulls->GetYaxis()->SetTitleFont(42);
        pulls->GetYaxis()->SetLabelSize(0.03);
        //pulls->SetMarkerColor(4);
        //pulls->SetMarkerStyle(21);
        //pulls->SetMarkerSize(0.6);

	lineResid = new TLine(-0.5,0,0.5,0);
        lineResid->SetLineStyle(3);
        lineResid->Draw();

        textResid = new TLatex();
        textResid->SetNDC();
        textResid->SetTextAlign(11);
        textResid->SetTextFont(42);
        textResid->SetTextSizePixels(17);
        textResid->SetTextSize(0.028);
        //if(EndCaps == 0) textResid->DrawLatex(0.80, 0.88, "Barrel");
        //if(EndCaps == 1) textResid->DrawLatex(0.77, 0.88, "End Caps");	
	textResid->DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
	if(isMC == 1) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV");
        if(isMC == 0) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV, #int L = 4.93 fb^{-1}");
        //textResid->DrawLatex(0.16, 0.80, EndCapsR9Chain->c_str());
        //textResid->DrawLatex(0.16, 0.75, tempLegChain->c_str());
	pulls->GetXaxis()->SetLimits(-0.5,0.5);	



	c2->Clear();

        //delete EndCapsR9Chain;
        //delete tempLegChain;
        res->Delete();
	res = 0;
        //mcs->Delete();
        //frame->Delete();
        residuals->Delete();
	residuals = 0;
        pulls->Delete();
	pulls = 0;
        textResid->Delete();
	textResid = 0;
        lineResid->Delete();
	lineResid = 0;
        Data_subset->Delete();
	Data_subset = 0;
        //Data->Delete();
        ntplVars->Delete();
	ntplVars = 0;
        mmg_s->Delete();
	mmg_s = 0;
        Photon_SC_Eta->Delete();
	Photon_SC_Eta = 0;
        Photon_SC_Phi->Delete();
	Photon_SC_Phi = 0;
        Photon_r9->Delete();
	Photon_r9 = 0;
        isJanLooseMMG->Delete();
	isJanLooseMMG = 0;
        isMultipleCandidate->Delete();
	isMultipleCandidate = 0;
        Photon_Et->Delete();
	Photon_Et = 0;
        Photon_SC_rawEt->Delete();
	Photon_SC_rawEt = 0;
        Photon_E->Delete();
	Photon_E = 0;
        Photon_SC_brem->Delete();
	Photon_SC_brem = 0;
        weight_pileUp->Delete();
	weight_pileUp = 0;
        meanV->Delete();
	meanV = 0;
        sigma->Delete();
	sigma = 0;
        width->Delete();
	width = 0;
        Voigtian->Delete();
	Voigtian = 0;

	Data_subset2->Delete();
        Data_subset2 = 0;
        ntplVars2->Delete();
        ntplVars2 = 0;
        mmg_s2->Delete();
        mmg_s2 = 0;

}

void Cruijff(double * mean_value, double * mean_error, Double_t * sigmaEff_value, double * sigma_value, double * sigma_value_error, double * sigmaR_value, double * sigmaL_value, double * ChiSquare, double * minLogLikelihood, double * differenceLogLikelihood, TH1D * hh, int j, int EndCaps, TString temp, TTree* Tree_Data, double mean, double rms, TCanvas * c2, RooPlot* Erawframe, TF1 * f, string nomDossier, string nomFichier, double RangeMin, double RangeMax, double * JanChi2, double * DegreesOfFreedom, double * pValue, int r9sup, double MinVar, double MaxVar, string variableX, int isMC)
{

	//RooRealVar * mmg_s = new RooRealVar("mmg_s", "s", -0.5, 0.5);
        RooRealVar * mmg_s = new RooRealVar("mmg_s", "s", RangeMin, RangeMax);
	RooRealVar * Photon_SC_Eta = new RooRealVar("Photon_SC_Eta", "Photon_SC_Eta", -3.0, 3.0, "");
        RooRealVar * Photon_SC_Phi = new RooRealVar("Photon_SC_Phi", "Photon_SC_Phi", -3.5, 3.5, "");
        RooRealVar * Photon_r9 = new RooRealVar("Photon_r9", "Photon_r9", 0.0, 1.0, "");
        RooRealVar * isJanLooseMMG = new RooRealVar("isJanLooseMMG", "isJanLooseMMG", 0.0, 1.0, "");
        RooRealVar * isMultipleCandidate = new RooRealVar("isMultipleCandidate", "isMultipleCandidate", 0.0, 1.0, "");
        RooRealVar * Photon_Et = new RooRealVar("Photon_Et", "Photon_Et", 0.0, 250.0, "");
        RooRealVar * Photon_SC_rawEt = new RooRealVar("Photon_SC_rawEt", "Photon_SC_rawEt", 0.0, 250.0, "");
        RooRealVar * Photon_E = new RooRealVar("Photon_E", "Photon_E", 0.0, 1000.0, "");
        RooRealVar * Photon_SC_brem = new RooRealVar("Photon_SC_brem", "Photon_SC_brem", 0.0, 15.0, "");
        RooRealVar * weight_pileUp = new RooRealVar("weight_pileUp", "weight_pileUp", 0.0, 100);

        RooArgSet *ntplVars = new RooArgSet(*mmg_s, *Photon_SC_Eta, *Photon_SC_Phi, *Photon_r9, *Photon_Et, *Photon_SC_rawEt, *Photon_E, *Photon_SC_brem, *isJanLooseMMG);


        ntplVars->add(*isMultipleCandidate);
        ntplVars->add(*weight_pileUp);


	RooDataSet *Data_subset = new RooDataSet("Data_subset", "Data_subset", Tree_Data, *ntplVars, temp, "weight_pileUp");

	RooRealVar * m0 = new RooRealVar("m0","m0",0.0,-0.1,0.1);
        RooRealVar * sigmaL = new RooRealVar("sigmaL","sigmaL",0.5,0.0,1.0);
	RooRealVar * sigmaR = new RooRealVar("sigmaR","sigmaR",0.5,0.0,1.0);
	RooRealVar * alphaL = new RooRealVar("alphaL","alphaL",0.5,0.0,1.0);
	RooRealVar * alphaR = new RooRealVar("alphaR","alphaR",0.5,0.0,1.0);

	// --- Build Cruijff PDF ---
	RooCruijff * Cruijff = new RooCruijff("Cruijff","Cruijff",*mmg_s,*m0,*sigmaL,*sigmaR,*alphaL,*alphaR);


	int fewBins = 1;
        Double_t Chi2J;

	RooFitResult *res;
	Double_t minNll;
		
        //int RightBinning = 50;
	//RightBinning = (int) ((RangeMax - RangeMin) / 0.02);
	RooBinning b(-17.0,17.0) ;
        b.addUniform(1700,-17.0,17.0) ;
	for(int i = 0; i<4; i++)
	{
                Erawframe->Clear();
		Cruijff = new RooCruijff("Cruijff","Cruijff",*mmg_s,*m0,*sigmaL,*sigmaR,*alphaL,*alphaR);

                //Erawframe = mmg_s->frame();
		Erawframe = mmg_s->frame(-0.5,0.5);
                Data_subset->plotOn(Erawframe,Name("myhist"),Binning(b),DataError(RooAbsData::SumW2));
		res = Cruijff->fitTo(*Data_subset, Range(RangeMin, RangeMax),Save(),SumW2Error(kTRUE));
		//res = Cruijff->fitTo(*Data_subset,Save(),SumW2Error(kTRUE));
		res->Print();
		minNll = res->minNll();
		Cruijff->plotOn(Erawframe,Name("mycurve"));
		Erawframe->Draw();


                Chi2J = chiSquare(Erawframe,(char *)"mycurve",(char *)"myhist",5, JanChi2, DegreesOfFreedom, pValue, &fewBins);
                if(fewBins == 0) break; 


        }

	c2->Clear();
        Erawframe->Clear();

        ///// Draw a nice plot /////

        RooRealVar * mmg_s2 = new RooRealVar("mmg_s", "s", -0.5, 0.5);
        RooArgSet *ntplVars2 = new RooArgSet(*mmg_s2, *Photon_SC_Eta, *Photon_SC_Phi, *Photon_r9, *Photon_Et, *Photon_SC_rawEt, *Photon_E, *Photon_SC_brem, *isJanLooseMMG);
        ntplVars2->add(*isMultipleCandidate);
        ntplVars2->add(*weight_pileUp); 
        RooDataSet *Data_subset2 = new RooDataSet("Data_subset2", "Data_subset2", Tree_Data, *ntplVars2, temp, "weight_pileUp");
        Erawframe = mmg_s2->frame();
        Data_subset2->plotOn(Erawframe,Name("myhist"),Binning(b));
        Cruijff->plotOn(Erawframe,Name("mycurve"),Range(RangeMin, RangeMax));
        Erawframe->Draw();

        ///// /////

	cout << endl << "ICI !!!!!!!!!!!!!!!!!! "<<Data_subset->weight() << endl ;

        f = Cruijff->asTF( RooArgList(*mmg_s) );
	*mean_value  = m0->getVal();
        *mean_error  = m0->getError();
        *sigmaEff_value = effSigma(hh);
        //*sigma_value = sigwidth.getVal();
        *sigmaR_value = SigmaR(f, 0.0, 2.0);
        *sigmaL_value = SigmaL(f, 0.0, 2.0);
        *ChiSquare = Erawframe->chiSquare();
	*minLogLikelihood = minNll;


        double entries = hh->GetEntries();

        int entriesInt = (int) entries;

	
	double fxmax = f->GetMaximumX(0.8,1.2,1.E-10,100,false);
        //cout<<"////////// ---- FXMAX BG = "<<f->GetMaximumX(0.8,1.2,1.E-10,100,false)<<" ---- //////////"<<endl;

	TLatex latexLabel;
        latexLabel.SetTextSize(0.030);
        latexLabel.SetNDC();
	latexLabel.DrawLatex(0.13, 0.96, "CMS Preliminary 2011, #sqrt{s} = 7 TeV");
        if(isMC == 1) latexLabel.DrawLatex(0.17, 0.88, "Simulation");
        if(isMC == 0) latexLabel.DrawLatex(0.17, 0.88, "Data, #int L = 4.93 fb^{-1}"); //lumi a changer !!!

	if(EndCaps == 0) latexLabel.DrawLatex(0.17, 0.83,"ECAL Barrel");
        if(EndCaps == 1) latexLabel.DrawLatex(0.17, 0.83,"ECAL Endcaps");
        if(r9sup == 0 && EndCaps == 0) latexLabel.DrawLatex(0.17, 0.78,"E_{T}^{#gamma} > 25 GeV, r9 < 0,94");
        if(r9sup == 0 && EndCaps == 1) latexLabel.DrawLatex(0.17, 0.78,"E_{T}^{#gamma} > 25 GeV, r9 < 0,95");
        if(r9sup == 1 && EndCaps == 0) latexLabel.DrawLatex(0.17, 0.78,"E_{T}^{#gamma} > 25 GeV, r9 > 0,94");
        if(r9sup == 1 && EndCaps == 1) latexLabel.DrawLatex(0.17, 0.78,"E_{T}^{#gamma} > 25 GeV, r9 > 0,95");
        if(r9sup == 2) latexLabel.DrawLatex(0.17, 0.78,"E_{T}^{#gamma} > 25 GeV, All r9");
	
	//latexLabel.DrawLatex(0.17, 0.73, Form("Entries = %d",entriesInt));


	latexLabel.DrawLatex(0.65, 0.88, Form("#color[1]{#mu(E_{#gamma}) = %f #pm %f}",m0->getVal(), m0->getError()));
	latexLabel.DrawLatex(0.65, 0.83, Form("#color[1]{#sigma_{L} = %f #pm %f}",sigmaL->getVal(), sigmaL->getError()));
	latexLabel.DrawLatex(0.65, 0.78, Form("#color[1]{#sigma_{R} = %f #pm %f}",sigmaR->getVal(), sigmaR->getError()));
	latexLabel.DrawLatex(0.65, 0.73, Form("#color[1]{#alpha_{L} = %f #pm %f}",alphaL->getVal(), alphaL->getError()));
	latexLabel.DrawLatex(0.65, 0.68, Form("#color[1]{#alpha_{R} = %f #pm %f}",alphaR->getVal(), alphaR->getError()));

	enregistrementPlots(nomDossier, nomFichier, EndCaps, j, c2);

	c2->Clear();
/*
	////////// RooMCStudy //////////

	
	RooMCStudy* mcs = new RooMCStudy(*Cruijff,*mmg_s,FitModel(*Cruijff),Silence(),Binned()) ;
        RooChi2MCSModule chi2mod ;
        mcs->addModule(chi2mod) ;

        mcs->generateAndFit(2000,entriesInt) ;

        //RooPlot* frame = mcs->plotNLL(Bins(50));
        RooPlot* frame = mcs->plotNLL(Name("myhistTest"));
	frame->Draw();

	//RooRealVar * VArTest = (RooRealVar*) frame->getPlotVar(); //artefact

	
	RooHist* histTest = frame->getHist("myhistTest");
	cout<<endl<<"histTest->GetMean() = "<<histTest->GetMean()<<endl;

	*differenceLogLikelihood = fabs(histTest->GetMean() - minNll);
*/
	TLatex *textResid = new TLatex();
/*
        textResid = new TLatex();
        textResid->SetNDC();
        textResid->SetTextAlign(11);
        textResid->SetTextFont(42);
        textResid->SetTextSizePixels(17);
        textResid->SetTextSize(0.028);
        //if(EndCaps == 0) textResid->DrawLatex(0.80, 0.88, "Barrel");
        //if(EndCaps == 1) textResid->DrawLatex(0.77, 0.88, "End Caps");
        textResid->DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
        if(isMC == 1) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV");
        if(isMC == 0) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV, #int L = 4.93 fb^{-1}");
        //textResid->DrawLatex(0.16, 0.80, EndCapsR9Chain->c_str());
        //textResid->DrawLatex(0.16, 0.75, tempLegChain->c_str());

	textResid->DrawLatex(0.61, 0.90, Form("#color[4]{Fit_minll = %f}",minNll));
	textResid->DrawLatex(0.61, 0.85, Form("#color[4]{Distrib_mean = %f}",histTest->GetMean()));
	textResid->DrawLatex(0.61, 0.80, Form("#color[3]{Difference = %f}",*differenceLogLikelihood));
	
	//RooAbsRealLValue VArTest = frame->getPlotVar();


	enregistrementPlots(nomDossier, "LogLikelihood", EndCaps, j, c2);

	c2->Clear();

*/
	////////// chi2 Pulls & Residuals //////////

	string nomDossierPull = "";
        if(EndCaps == 0) nomDossierPull = nomDossier + "EB/";
        if(EndCaps == 1) nomDossierPull = nomDossier + "EE/";


        RooHist* residuals = residHist(Erawframe,(char *)"myhist",(char *)"mycurve", false, nomDossierPull, j);
	

        //RooHist* pulls = pullHist(Erawframe,"myhist", "mycurve");
        RooHist* pulls = residHist(Erawframe,(char *)"myhist", (char *)"mycurve", true, nomDossierPull, j);


	residuals->Draw("AP");
	nomFichier = "Chi2Residuals";
	
	residuals->GetYaxis()->SetTitle("#chi^{2} Residuals");
        residuals->GetXaxis()->SetLabelFont(42);
        residuals->GetXaxis()->SetTitleFont(42);
        residuals->GetXaxis()->SetLabelSize(0.03);
        residuals->GetXaxis()->SetTitle("s");
        residuals->GetYaxis()->SetLabelFont(42);
        residuals->GetYaxis()->SetTitleOffset(1.24);
        residuals->GetYaxis()->SetTitleFont(42);
        residuals->GetYaxis()->SetLabelSize(0.03);
        //residuals->SetMarkerColor(4);
        //residuals->SetMarkerStyle(21);
        //residuals->SetMarkerSize(0.6);

	TLine *lineResid = new TLine(-0.5,0,0.5,0);
        lineResid->SetLineStyle(3);
        lineResid->Draw();


        //TLatex *textResid = new TLatex();

        textResid = new TLatex();
        textResid->SetNDC();
        textResid->SetTextAlign(11);
        textResid->SetTextFont(42);
        textResid->SetTextSizePixels(17);
        textResid->SetTextSize(0.028);
        //if(EndCaps == 0) textResid->DrawLatex(0.80, 0.88, "Barrel");
        //if(EndCaps == 1) textResid->DrawLatex(0.77, 0.88, "End Caps");
	textResid->DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
	if(isMC == 1) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV");
        if(isMC == 0) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV, #int L = 4.93 fb^{-1}");
        //textResid->DrawLatex(0.16, 0.80, EndCapsR9Chain->c_str());
        //textResid->DrawLatex(0.16, 0.75, tempLegChain->c_str());
	residuals->GetXaxis()->SetLimits(-0.5,0.5);
	
	enregistrementPlots(nomDossier, nomFichier, EndCaps, j, c2);

	c2->Clear();

	pulls->Draw("AP");
	nomFichier = "Chi2Pulls";
	
	pulls->GetYaxis()->SetTitle("#chi^{2} Pulls");
        pulls->GetXaxis()->SetLabelFont(42);
        pulls->GetXaxis()->SetTitleFont(42);
        pulls->GetXaxis()->SetLabelSize(0.03);
        pulls->GetXaxis()->SetTitle("s");
        pulls->GetYaxis()->SetLabelFont(42);
        pulls->GetYaxis()->SetTitleOffset(1.24);
        pulls->GetYaxis()->SetTitleFont(42);
        pulls->GetYaxis()->SetLabelSize(0.03);
        //pulls->SetMarkerColor(4);
        //pulls->SetMarkerStyle(21);
        //pulls->SetMarkerSize(0.6);

	lineResid = new TLine(-0.5,0,0.5,0);
        lineResid->SetLineStyle(3);
        lineResid->Draw();

        textResid = new TLatex();
        textResid->SetNDC();
        textResid->SetTextAlign(11);
        textResid->SetTextFont(42);
        textResid->SetTextSizePixels(17);
        textResid->SetTextSize(0.028);
        //if(EndCaps == 0) textResid->DrawLatex(0.80, 0.88, "Barrel");
        //if(EndCaps == 1) textResid->DrawLatex(0.77, 0.88, "End Caps");	
	textResid->DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
	if(isMC == 1) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV");
        if(isMC == 0) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV, #int L = 4.93 fb^{-1}");
        //textResid->DrawLatex(0.16, 0.80, EndCapsR9Chain->c_str());
        //textResid->DrawLatex(0.16, 0.75, tempLegChain->c_str());
	pulls->GetXaxis()->SetLimits(-0.5,0.5);	



	c2->Clear();

        //delete EndCapsR9Chain;
        //delete tempLegChain;
        res->Delete();
	res = 0;
        //mcs->Delete();
        //frame->Delete();
        residuals->Delete();
	residuals = 0;
        pulls->Delete();
	pulls = 0;
        textResid->Delete();
	textResid = 0;
        lineResid->Delete();
	lineResid = 0;
        Data_subset->Delete();
	Data_subset = 0;
        //Data->Delete();
        ntplVars->Delete();
	ntplVars = 0;
        mmg_s->Delete();
	mmg_s = 0;
        Photon_SC_Eta->Delete();
	Photon_SC_Eta = 0;
        Photon_SC_Phi->Delete();
	Photon_SC_Phi = 0;
        Photon_r9->Delete();
	Photon_r9 = 0;
        isJanLooseMMG->Delete();
	isJanLooseMMG = 0;
        isMultipleCandidate->Delete();
	isMultipleCandidate = 0;
        Photon_Et->Delete();
	Photon_Et = 0;
        Photon_SC_rawEt->Delete();
	Photon_SC_rawEt = 0;
        Photon_E->Delete();
	Photon_E = 0;
        Photon_SC_brem->Delete();
	Photon_SC_brem = 0;
        weight_pileUp->Delete();
	weight_pileUp = 0;
        m0->Delete();
	m0 = 0;
        sigmaL->Delete();
	sigmaL = 0;
        sigmaR->Delete();
        sigmaR = 0;
	alphaL->Delete();
        alphaL = 0;
	alphaR->Delete();
        alphaR = 0;
	Cruijff->Delete();
	Cruijff = 0;

	Data_subset2->Delete();
        Data_subset2 = 0;
        ntplVars2->Delete();
        ntplVars2 = 0;
        mmg_s2->Delete();
        mmg_s2 = 0;
}

void roofit_bwcb_fft2(double * mean_value, double * mean_error, Double_t * sigmaEff_value, double * sigma_value, double * sigma_value_error, double * sigmaR_value, double * sigmaL_value, double * ChiSquare, double * minLogLikelihood, double * differenceLogLikelihood, TH1D * hh, int j, int EndCaps, TString temp, TTree* Tree_Data, double mean, double rms, TCanvas * c2, RooPlot* Erawframe, TF1 * f, string nomDossier, string nomFichier, double RangeMin, double RangeMax, double * JanChi2, double * DegreesOfFreedom, double * pValue, int r9sup, double MinVar, double MaxVar, string variableX, int isMC)
{
	//RooRealVar * mmg_s = new RooRealVar("mmg_s", "s", -0.5, 0.5);
        RooRealVar * mmg_s = new RooRealVar("mmg_s", "s", RangeMin, RangeMax);
	RooRealVar * Photon_SC_Eta = new RooRealVar("Photon_SC_Eta", "Photon_SC_Eta", -3.0, 3.0, "");
        RooRealVar * Photon_SC_Phi = new RooRealVar("Photon_SC_Phi", "Photon_SC_Phi", -3.5, 3.5, "");
        RooRealVar * Photon_r9 = new RooRealVar("Photon_r9", "Photon_r9", 0.0, 1.0, "");
        RooRealVar * isJanLooseMMG = new RooRealVar("isJanLooseMMG", "isJanLooseMMG", 0.0, 1.0, "");
        RooRealVar * isMultipleCandidate = new RooRealVar("isMultipleCandidate", "isMultipleCandidate", 0.0, 1.0, "");
        RooRealVar * Photon_Et = new RooRealVar("Photon_Et", "Photon_Et", 0.0, 250.0, "");
        RooRealVar * Photon_SC_rawEt = new RooRealVar("Photon_SC_rawEt", "Photon_SC_rawEt", 0.0, 250.0, "");
        RooRealVar * Photon_E = new RooRealVar("Photon_E", "Photon_E", 0.0, 1000.0, "");
        RooRealVar * Photon_SC_brem = new RooRealVar("Photon_SC_brem", "Photon_SC_brem", 0.0, 15.0, "");
        RooRealVar * weight_pileUp = new RooRealVar("weight_pileUp", "weight_pileUp", 0.0, 100);

        RooArgSet *ntplVars = new RooArgSet(*mmg_s, *Photon_SC_Eta, *Photon_SC_Phi, *Photon_r9, *Photon_Et, *Photon_SC_rawEt, *Photon_E, *Photon_SC_brem, *isJanLooseMMG);



        ntplVars->add(*isMultipleCandidate);
        ntplVars->add(*weight_pileUp);


	RooDataSet *Data_subset = new RooDataSet("Data_subset", "Data_subset", Tree_Data, *ntplVars, temp, "weight_pileUp");


	RooRealVar * BWmean = new RooRealVar("BWmean","BWmean",mean, mean - rms, mean + rms, "GeV");
        RooRealVar * BWsigma = new RooRealVar("BWsigma","BWsigma",rms, 0.0, 2 * rms, "GeV");

	RooRealVar * CBmean = new RooRealVar("CBmean","CBmean",mean, mean - rms, mean + rms, "GeV");
        RooRealVar * CBsigma = new RooRealVar("CBsigma","CBsigma",rms, 0.0, 2 * rms, "GeV");
	RooRealVar * CBalpha = new RooRealVar("CBalpha","CBalpha",0.9,-5.0,5.0);
	RooRealVar * CBn = new RooRealVar("CBn","CBn",10.0,1.0,500.0);

  	RooBreitWigner * BW = new RooBreitWigner("BW","Breit Wigner theory",*mmg_s,*BWmean,*BWsigma);
	RooCBShape * cryBall = new RooCBShape("cryBall","Crystal Ball resolution model", *mmg_s, *CBmean, *CBsigma, *CBalpha, *CBn) ;
	RooFFTConvPdf * bwxCryBall  = new RooFFTConvPdf("bwxCryBall", "FFT Conv CryBall and BW", *mmg_s, *BW, *cryBall); 

	int fewBins = 1;
        Double_t Chi2J;

	RooFitResult *res;
	Double_t minNll;
		

        //for(int RightBinning = 600; RightBinning > 10; RightBinning -= 10)
	//int RightBinning = 50;
        //RightBinning = (int) ((RangeMax - RangeMin) / 0.02);
	RooBinning b(-17.0,17.0) ;
	b.addUniform(1700,-17.0,17.0) ;
	for(int i = 0; i<4; i++)
	{
                Erawframe->Clear();
                //Data_subset = (RooDataSet*)Data->reduce(*ntplVars, temp);
		bwxCryBall = new RooFFTConvPdf("bwxCryBall","bwxCryBall",*mmg_s, *BW, *cryBall);

                Erawframe = mmg_s->frame(-0.5,0.5);
                //mmg_s->frame(-0.5,0.5) ;
		//Data_subset->plotOn(Erawframe,Name("myhist"),Binning(RightBinning),DataError(RooAbsData::SumW2));
		//Data_subset->plotOn(Erawframe,Name("myhist"),DataError(RooAbsData::SumW2));
		Data_subset->plotOn(Erawframe,Name("myhist"),Binning(b),DataError(RooAbsData::SumW2));
		//bwxCryBall->fitTo(*Data_subset, Range(RangeMin, RangeMax),Minos(true));
                //res = bwxCryBall->fitTo(*Data_subset, Range(RangeMin, RangeMax),Minos(true),Save(),SumW2Error(kFALSE));
		res = bwxCryBall->fitTo(*Data_subset, Range(RangeMin, RangeMax),Save(),SumW2Error(kTRUE));
		//res = bwxCryBall->fitTo(*Data_subset,Save(),SumW2Error(kTRUE));
		res->Print();
		minNll = res->minNll();
		bwxCryBall->plotOn(Erawframe,Name("mycurve"));
		//cryBall->plotOn(Erawframe,LineColor(kRed),LineStyle(kDashed));
		//BW->plotOn(Erawframe,LineColor(kGreen),LineStyle(kDashed));
		Erawframe->Draw();


                Chi2J = chiSquare(Erawframe,(char *)"mycurve",(char *)"myhist",6, JanChi2, DegreesOfFreedom, pValue, &fewBins);
                if(fewBins == 0) break; 


        }

	c2->Clear();
        Erawframe->Clear();

        ///// Draw a nice plot /////

        RooRealVar * mmg_s2 = new RooRealVar("mmg_s", "s", -0.5, 0.5);
        RooArgSet *ntplVars2 = new RooArgSet(*mmg_s2, *Photon_SC_Eta, *Photon_SC_Phi, *Photon_r9, *Photon_Et, *Photon_SC_rawEt, *Photon_E, *Photon_SC_brem, *isJanLooseMMG);
        ntplVars2->add(*isMultipleCandidate);
        ntplVars2->add(*weight_pileUp); 
        RooDataSet *Data_subset2 = new RooDataSet("Data_subset2", "Data_subset2", Tree_Data, *ntplVars2, temp, "weight_pileUp");
        Erawframe = mmg_s2->frame();
        Data_subset2->plotOn(Erawframe,Name("myhist"),Binning(b));
        bwxCryBall->plotOn(Erawframe,Name("mycurve"),Range(RangeMin, RangeMax));
        Erawframe->Draw();

        ///// /////

	cout << endl << "ICI !!!!!!!!!!!!!!!!!! "<<Data_subset->weight() << endl ;

	RooDerivative * deriv = bwxCryBall->derivative(*mmg_s,1) ;

        // Find point where derivative is zero
        //Double_t x_max = deriv->findRoot(*mmg_s,RangeMin, RangeMax,0) ;
	Double_t x_max = deriv->findRoot(*mmg_s,-0.1, 0.1,0) ;

        f = bwxCryBall->asTF( RooArgList(*mmg_s) );
	*mean_value  = x_max;//CBmean->getVal();
        *mean_error  = CBmean->getError();
        *sigmaEff_value = effSigma(hh);
        //*sigma_value = sigwidth.getVal();
        *sigmaR_value = SigmaR(f, 0.0, 2.0);
        *sigmaL_value = SigmaL(f, 0.0, 2.0);
        *ChiSquare = Erawframe->chiSquare();
	*minLogLikelihood = minNll;


        double entries = hh->GetEntries();

        int entriesInt = (int) entries;

	
	double fxmax = f->GetMaximumX(0.8,1.2,1.E-10,100,false);
        //cout<<"////////// ---- FXMAX BG = "<<f->GetMaximumX(0.8,1.2,1.E-10,100,false)<<" ---- //////////"<<endl;

	TLatex latexLabel;
        latexLabel.SetTextSize(0.030);
        latexLabel.SetNDC();
	latexLabel.DrawLatex(0.13, 0.96, "CMS Preliminary 2011, #sqrt{s} = 7 TeV");
        if(isMC == 1) latexLabel.DrawLatex(0.17, 0.88, "Simulation");
        if(isMC == 0) latexLabel.DrawLatex(0.17, 0.88, "Data, #int L = 4.93 fb^{-1}"); //lumi a changer !!!

	if(EndCaps == 0) latexLabel.DrawLatex(0.17, 0.83,"ECAL Barrel");
        if(EndCaps == 1) latexLabel.DrawLatex(0.17, 0.83,"ECAL Endcaps");
        if(r9sup == 0 && EndCaps == 0) latexLabel.DrawLatex(0.17, 0.78,"E_{T}^{#gamma} > 25 GeV, r9 < 0,94");
        if(r9sup == 0 && EndCaps == 1) latexLabel.DrawLatex(0.17, 0.78,"E_{T}^{#gamma} > 25 GeV, r9 < 0,95");
        if(r9sup == 1 && EndCaps == 0) latexLabel.DrawLatex(0.17, 0.78,"E_{T}^{#gamma} > 25 GeV, r9 > 0,94");
        if(r9sup == 1 && EndCaps == 1) latexLabel.DrawLatex(0.17, 0.78,"E_{T}^{#gamma} > 25 GeV, r9 > 0,95");
        if(r9sup == 2) latexLabel.DrawLatex(0.17, 0.78,"E_{T}^{#gamma} > 25 GeV, All r9");
	
	//latexLabel.DrawLatex(0.17, 0.73, Form("Entries = %d",entriesInt));


	latexLabel.DrawLatex(0.65, 0.88, Form("#color[1]{mean_{CB} = %f #pm %f}",CBmean->getVal(), CBmean->getError()));
	latexLabel.DrawLatex(0.65, 0.83, Form("#color[1]{#sigma_{CB} = %f #pm %f}",CBsigma->getVal(), CBsigma->getError()));
	latexLabel.DrawLatex(0.65, 0.78, Form("#color[1]{#alpha_{CB} = %f #pm %f}",CBalpha->getVal(), CBalpha->getError()));
	latexLabel.DrawLatex(0.65, 0.73, Form("#color[1]{n_{CB} = %f #pm %f}",CBn->getVal(), CBn->getError()));	
	latexLabel.DrawLatex(0.65, 0.68, Form("#color[1]{mean_{BW} = %f #pm %f}",BWmean->getVal(), BWmean->getError()));
	latexLabel.DrawLatex(0.65, 0.63, Form("#color[1]{#sigma_{BW} = %f #pm %f}",BWsigma->getVal(),BWsigma->getError()));
	latexLabel.DrawLatex(0.65, 0.58, Form("#color[1]{x_{max} = %f}", x_max));

	enregistrementPlots(nomDossier, nomFichier, EndCaps, j, c2);

	c2->Clear();
/*
	////////// RooMCStudy //////////

	
	RooMCStudy* mcs = new RooMCStudy(*bwxCryBall,*mmg_s,FitModel(*bwxCryBall),Silence(),Binned()) ;
        RooChi2MCSModule chi2mod ;
        mcs->addModule(chi2mod) ;

        mcs->generateAndFit(2000,entriesInt) ;

        //RooPlot* frame = mcs->plotNLL(Bins(50));
        RooPlot* frame = mcs->plotNLL(Name("myhistTest"));
	frame->Draw();

	//RooRealVar * VArTest = (RooRealVar*) frame->getPlotVar(); //artefact

	
	RooHist* histTest = frame->getHist("myhistTest");
	cout<<endl<<"histTest->GetMean() = "<<histTest->GetMean()<<endl;

	*differenceLogLikelihood = fabs(histTest->GetMean() - minNll);
*/
	TLatex *textResid = new TLatex();
/*
        textResid = new TLatex();
        textResid->SetNDC();
        textResid->SetTextAlign(11);
        textResid->SetTextFont(42);
        textResid->SetTextSizePixels(17);
        textResid->SetTextSize(0.028);
        //if(EndCaps == 0) textResid->DrawLatex(0.80, 0.88, "Barrel");
        //if(EndCaps == 1) textResid->DrawLatex(0.77, 0.88, "End Caps");
        textResid->DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
        if(isMC == 1) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV");
        if(isMC == 0) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV, #int L = 4.93 fb^{-1}");
        //textResid->DrawLatex(0.16, 0.80, EndCapsR9Chain->c_str());
        //textResid->DrawLatex(0.16, 0.75, tempLegChain->c_str());

	textResid->DrawLatex(0.61, 0.90, Form("#color[4]{Fit_minll = %f}",minNll));
	textResid->DrawLatex(0.61, 0.85, Form("#color[4]{Distrib_mean = %f}",histTest->GetMean()));
	textResid->DrawLatex(0.61, 0.80, Form("#color[3]{Difference = %f}",*differenceLogLikelihood));
	
	//RooAbsRealLValue VArTest = frame->getPlotVar();


	enregistrementPlots(nomDossier, "LogLikelihood", EndCaps, j, c2);

	c2->Clear();

*/
	////////// chi2 Pulls & Residuals //////////

	string nomDossierPull = "";
        if(EndCaps == 0) nomDossierPull = nomDossier + "EB/";
        if(EndCaps == 1) nomDossierPull = nomDossier + "EE/";


        RooHist* residuals = residHist(Erawframe,(char *)"myhist",(char *)"mycurve", false, nomDossierPull, j);
	

        //RooHist* pulls = pullHist(Erawframe,"myhist", "mycurve");
        RooHist* pulls = residHist(Erawframe,(char *)"myhist", (char *)"mycurve", true, nomDossierPull, j);


	residuals->Draw("AP");
	nomFichier = "Chi2Residuals";
	
	residuals->GetYaxis()->SetTitle("#chi^{2} Residuals");
        residuals->GetXaxis()->SetLabelFont(42);
        residuals->GetXaxis()->SetTitleFont(42);
        residuals->GetXaxis()->SetLabelSize(0.03);
        residuals->GetXaxis()->SetTitle("s");
        residuals->GetYaxis()->SetLabelFont(42);
        residuals->GetYaxis()->SetTitleOffset(1.24);
        residuals->GetYaxis()->SetTitleFont(42);
        residuals->GetYaxis()->SetLabelSize(0.03);
        //residuals->SetMarkerColor(4);
        //residuals->SetMarkerStyle(21);
        //residuals->SetMarkerSize(0.6);

	TLine *lineResid = new TLine(-0.5,0,0.5,0);
        lineResid->SetLineStyle(3);
        lineResid->Draw();


        //TLatex *textResid = new TLatex();

        textResid = new TLatex();
        textResid->SetNDC();
        textResid->SetTextAlign(11);
        textResid->SetTextFont(42);
        textResid->SetTextSizePixels(17);
        textResid->SetTextSize(0.028);
        //if(EndCaps == 0) textResid->DrawLatex(0.80, 0.88, "Barrel");
        //if(EndCaps == 1) textResid->DrawLatex(0.77, 0.88, "End Caps");
	textResid->DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
	if(isMC == 1) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV");
        if(isMC == 0) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV, #int L = 4.93 fb^{-1}");
        //textResid->DrawLatex(0.16, 0.80, EndCapsR9Chain->c_str());
        //textResid->DrawLatex(0.16, 0.75, tempLegChain->c_str());
	residuals->GetXaxis()->SetLimits(-0.5,0.5);
	
	enregistrementPlots(nomDossier, nomFichier, EndCaps, j, c2);

	c2->Clear();

	pulls->Draw("AP");
	nomFichier = "Chi2Pulls";
	
	pulls->GetYaxis()->SetTitle("#chi^{2} Pulls");
        pulls->GetXaxis()->SetLabelFont(42);
        pulls->GetXaxis()->SetTitleFont(42);
        pulls->GetXaxis()->SetLabelSize(0.03);
        pulls->GetXaxis()->SetTitle("s");
        pulls->GetYaxis()->SetLabelFont(42);
        pulls->GetYaxis()->SetTitleOffset(1.24);
        pulls->GetYaxis()->SetTitleFont(42);
        pulls->GetYaxis()->SetLabelSize(0.03);
        //pulls->SetMarkerColor(4);
        //pulls->SetMarkerStyle(21);
        //pulls->SetMarkerSize(0.6);

	lineResid = new TLine(-0.5,0,0.5,0);
        lineResid->SetLineStyle(3);
        lineResid->Draw();

        textResid = new TLatex();
        textResid->SetNDC();
        textResid->SetTextAlign(11);
        textResid->SetTextFont(42);
        textResid->SetTextSizePixels(17);
        textResid->SetTextSize(0.028);
        //if(EndCaps == 0) textResid->DrawLatex(0.80, 0.88, "Barrel");
        //if(EndCaps == 1) textResid->DrawLatex(0.77, 0.88, "End Caps");	
	textResid->DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
	if(isMC == 1) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV");
        if(isMC == 0) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV, #int L = 4.93 fb^{-1}");
        //textResid->DrawLatex(0.16, 0.80, EndCapsR9Chain->c_str());
        //textResid->DrawLatex(0.16, 0.75, tempLegChain->c_str());
	pulls->GetXaxis()->SetLimits(-0.5,0.5);	



	c2->Clear();

        //delete EndCapsR9Chain;
        //delete tempLegChain;
        res->Delete();
	res = 0;
        //mcs->Delete();
        //frame->Delete();
        residuals->Delete();
	residuals = 0;
        pulls->Delete();
	pulls = 0;
        textResid->Delete();
	textResid = 0;
        lineResid->Delete();
	lineResid = 0;
        Data_subset->Delete();
	Data_subset = 0;
        //Data->Delete();
        ntplVars->Delete();
	ntplVars = 0;
        mmg_s->Delete();
	mmg_s = 0;
        Photon_SC_Eta->Delete();
	Photon_SC_Eta = 0;
        Photon_SC_Phi->Delete();
	Photon_SC_Phi = 0;
        Photon_r9->Delete();
	Photon_r9 = 0;
        isJanLooseMMG->Delete();
	isJanLooseMMG = 0;
        isMultipleCandidate->Delete();
	isMultipleCandidate = 0;
        Photon_Et->Delete();
	Photon_Et = 0;
        Photon_SC_rawEt->Delete();
	Photon_SC_rawEt = 0;
        Photon_E->Delete();
	Photon_E = 0;
        Photon_SC_brem->Delete();
	Photon_SC_brem = 0;
        weight_pileUp->Delete();
	weight_pileUp = 0;
        CBmean->Delete();
	CBmean = 0;
        CBsigma->Delete();
	CBsigma = 0;
	CBalpha->Delete();
	CBalpha = 0;
	CBn->Delete();
	CBn = 0;
	BWmean->Delete();
	BWmean = 0;
	BWsigma->Delete();
	BWsigma = 0;
	cryBall->Delete();
	cryBall = 0;
	BW->Delete();
	BW = 0;
        bwxCryBall->Delete();
	bwxCryBall = 0;

	Data_subset2->Delete();
        Data_subset2 = 0;
        ntplVars2->Delete();
        ntplVars2 = 0;
        mmg_s2->Delete();
        mmg_s2 = 0;


}

void enregistrementPlots(string nomDossier, string nomFichier, int EndCaps, int iteration, TCanvas * c1)
{    
     
	if(EndCaps == 0) nomDossier += "EB/";
	if(EndCaps == 1) nomDossier += "EE/";
	//mkdir(nomDossier.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	system(Form("mkdir -p %s", nomDossier.c_str()));
	if(iteration != 10000) nomFichier += Form("%d",iteration);
        c1->Print(Form("%s%s.root",nomDossier.c_str(),nomFichier.c_str()));
        c1->Print(Form("%s%s.C",nomDossier.c_str(),nomFichier.c_str()));
        c1->Print(Form("%s%s.pdf",nomDossier.c_str(),nomFichier.c_str()));
        c1->Print(Form("%s%s.ps",nomDossier.c_str(),nomFichier.c_str()));
        c1->Print(Form("%s%s.png",nomDossier.c_str(),nomFichier.c_str()));
        return;
}

void RangeEstimator(double pourcentage, double centralValue, TChain * chain, int Endcaps, double * MinRange, double * MaxRange, double Varmin, double Varmax)
{

	float mmg_s;
        float Photon_SC_Eta;
        float Photon_r9;
        float isJanLooseMMG;
        float isMultipleCandidate;
        float Photon_Et;
        chain->SetBranchAddress("mmg_s",&mmg_s);
        chain->SetBranchAddress("Photon_SC_Eta",&Photon_SC_Eta);
        chain->SetBranchAddress("Photon_r9",&Photon_r9);
        chain->SetBranchAddress("isJanLooseMMG",&isJanLooseMMG);
        chain->SetBranchAddress("isMultipleCandidate",&isMultipleCandidate);
        chain->SetBranchAddress("Photon_Et",&Photon_Et);
        vector <float> mmg_sVector;

	for (int ievt = 0 ; ievt < chain->GetEntries() ; ievt++)
        {
                chain->GetEntry(ievt);

                if(Endcaps == 0)
                {

                        if(( (fabs(Photon_SC_Eta)) > 0.018 || ( (fabs(Photon_SC_Eta)) < 0.423 && (fabs(Photon_SC_Eta)) > 0.461 ) || ( (fabs(Photon_SC_Eta)) < 0.770 && (fabs(Photon_SC_Eta)) > 0.806 ) || ( (fabs(Photon_SC_Eta)) < 1.127 && (fabs(Photon_SC_Eta)) > 1.163 ) ) && (fabs(Photon_SC_Eta)) < 1.4442 && isMultipleCandidate == 0 && isJanLooseMMG != 0 && Photon_r9 > 0.94 && Photon_Et > Varmin && Photon_Et <= Varmax)
                        {

                                mmg_sVector.push_back(mmg_s);

                        }


                }

                if(Endcaps == 1)
                {

        //cout<<endl<<BremVector[1000]<<endl;
                        if((fabs(Photon_SC_Eta)) > 1.566 && isJanLooseMMG != 0 && isMultipleCandidate == 0 && Photon_r9 > 0.95 && Photon_Et > Varmin && Photon_Et <= Varmax)  mmg_sVector.push_back(mmg_s);

                }


        }
        sort(mmg_sVector.begin(), mmg_sVector.end());

	cout<<endl<<"mmg_sVector.size() = "<<mmg_sVector.size()<<endl;

	int meanVector = 0;
	for(int h = 0; h < mmg_sVector.size(); h++)
	{
		if(mmg_sVector[h] >= centralValue)
		{
			meanVector = h;
			h = mmg_sVector.size();
			break;
		}


	}

        
	cout<<endl<<"meanVector = "<<meanVector<<endl;

        double stopVectorMoins = meanVector - mmg_sVector.size() * pourcentage / 2.0;
        double stopVectorPlus = meanVector + mmg_sVector.size() * pourcentage / 2.0;
        cout<<endl<<"stopVectorPlus = "<<stopVectorPlus<<endl;
        int intStopVectorMoins = (int) stopVectorMoins;
        int intStopVectorPlus = (int) stopVectorPlus;
        cout<<endl<<"intStopVectorMoins = "<<intStopVectorMoins<<endl;
        cout<<endl<<"intStopVectorPlus = "<<intStopVectorPlus<<endl;

        cout<<endl<<"mmg_sVector[meanVector] = "<<mmg_sVector[meanVector]<<endl;

        cout<<endl<<"Le range est compris entre : "<<mmg_sVector[intStopVectorMoins]<<" et "<<mmg_sVector[intStopVectorPlus]<<endl;

	*MinRange = mmg_sVector[intStopVectorMoins];
	*MaxRange = mmg_sVector[intStopVectorPlus]; 


}


void RangeEstimator2(double pourcentage, TChain * chain, int Endcaps, double * MinRange, double * MaxRange, double Varmin, double Varmax)
{

	float mmg_s;
        float Photon_SC_Eta;
        float Photon_r9;
        float isJanLooseMMG;
        float isMultipleCandidate;
        float Photon_Et;
        chain->SetBranchAddress("mmg_s",&mmg_s);
        chain->SetBranchAddress("Photon_SC_Eta",&Photon_SC_Eta);
        chain->SetBranchAddress("Photon_r9",&Photon_r9);
        chain->SetBranchAddress("isJanLooseMMG",&isJanLooseMMG);
        chain->SetBranchAddress("isMultipleCandidate",&isMultipleCandidate);
        chain->SetBranchAddress("Photon_Et",&Photon_Et);
        vector <float> mmg_sVector;

        for (int ievt = 0 ; ievt < chain->GetEntries() ; ievt++)
        {
                chain->GetEntry(ievt);

                if(Endcaps == 0)
                {

                        if(( (fabs(Photon_SC_Eta)) > 0.018 || ( (fabs(Photon_SC_Eta)) < 0.423 && (fabs(Photon_SC_Eta)) > 0.461 ) || ( (fabs(Photon_SC_Eta)) < 0.770 && (fabs(Photon_SC_Eta)) > 0.806 ) || ( (fabs(Photon_SC_Eta)) < 1.127 && (fabs(Photon_SC_Eta)) > 1.163 ) ) && (fabs(Photon_SC_Eta)) < 1.4442 && isMultipleCandidate == 0 && isJanLooseMMG != 0 && Photon_r9 > 0.94 && Photon_Et > Varmin && Photon_Et <= Varmax)
                        {

                                mmg_sVector.push_back(mmg_s);

                        }


                }

                if(Endcaps == 1)
                {

        //cout<<endl<<BremVector[1000]<<endl;
                        if((fabs(Photon_SC_Eta)) > 1.566 && isJanLooseMMG != 0 && isMultipleCandidate == 0 && Photon_r9 > 0.95 && Photon_Et > Varmin && Photon_Et <= Varmax)  mmg_sVector.push_back(mmg_s);

                }


        }
        sort(mmg_sVector.begin(), mmg_sVector.end());

        double Min = mmg_sVector.size() * (pourcentage / 4.0);
        int MinInt = (int) Min;
        double Max = mmg_sVector.size() * (3.0 * pourcentage / 4.0);
        int MaxInt = (int) Max;

        *MinRange = mmg_sVector[MinInt];
        *MaxRange = mmg_sVector[MaxInt];

}

void RangeEstimator3(double pourcentage, TChain * chain, TString temp, int Endcaps, double * MinRange, double * MaxRange)
{
      

        TChain * ReducedChain = (TChain *) chain->CopyTree(temp);
      
        float mmg_s;
        ReducedChain->SetBranchAddress("mmg_s",&mmg_s);

        vector <float> mmg_sVector;


        for (int ievt = 0 ; ievt < ReducedChain->GetEntries() ; ievt++)
        {
                ReducedChain->GetEntry(ievt);
                mmg_sVector.push_back(mmg_s);

        }


        sort(mmg_sVector.begin(), mmg_sVector.end());

        size_t interval_entries = TMath::Ceil(pourcentage * mmg_sVector.size());

        vector<float>::iterator lower = mmg_sVector.begin();
        vector<float>::iterator upper = mmg_sVector.begin() + interval_entries - 1; 

        double dx = *upper - *lower;

        for(vector<float>::iterator first = lower, last = upper; last < mmg_sVector.end(); first++, last++)
        {
                if((*last - *first) < dx)
                {

                        lower = first;
                        upper = last;
                        dx = *upper - *lower;
                }
                      
        }

        *MinRange = *lower;
        *MaxRange = *upper;

}



void SymetricRangeEstimator(TChain * chain, double centralValue, double * MinRange, double * MaxRange, double Entries, double pourcentage, TString temp)
{

	double min = 0.0;
	double max = 0.0;
	string minresult;
	string maxresult;
	std::ostringstream oss;
	std::ostringstream oss2;	
	TString temp2;
	

	for(double sigma = 0.001; sigma < 1.0; sigma += 0.001)
	{
		temp2.Clear();
		oss.str("");
		oss2.str("");
		minresult.clear();
		maxresult.clear();

		min = centralValue - sigma; 
		max = centralValue + sigma;

        	oss << min;
        	minresult = oss.str();
		
                oss2 << max; 
                maxresult = oss2.str();	
		
		TH1D *histo = new TH1D("histo","histo", 200, 0.0, 2.0);
		
		temp2 = temp;
		temp2 += " && mmg_s > ";
		temp2 += minresult;
		temp2 += " && mmg_s < ";
		temp2 += maxresult;

		chain->Draw("mmg_s>>histo", temp2);
		
		if(histo->GetEntries() >= (Entries * pourcentage))
		{
			*MinRange = centralValue - sigma;
			*MaxRange = centralValue + sigma;
			sigma = 1.0;
		}

		histo->Delete();	
	}



}

void SymetricRangeEstimator2(TChain * chain, double lastX, double * MinRange, double * MaxRange, double Entries, double pourcentage, TString temp)
{

        double min = 0.0;
        double max = 0.0;
        int itermin = 0;
        int itermax = 0;
	int loop1 = 0;
	int loop2 = 0;
        string minresult;
        string maxresult;
        std::ostringstream oss;
        std::ostringstream oss2;
        TString temp2;

	double limite = ((1.0 - pourcentage) / 2.0) * Entries;

        for(double sigma = 0.0; sigma < 2.0; sigma += 0.01)
        {
                temp2.Clear();
                oss.str("");
                oss2.str("");
                minresult.clear();
                maxresult.clear();

                min = sigma;
                max = lastX - sigma;

                oss << min;
                minresult = oss.str();

                oss2 << max;
                maxresult = oss2.str();

                TH1D *histo = new TH1D("histo","histo", 2000, 0.0, 2.0);
                TH1D *histo2 = new TH1D("histo2","histo2", 2000, 0.0, 2.0);

		temp2 = temp;
                temp2 += " && mmg_s < ";
                temp2 += minresult;

                chain->Draw("mmg_s>>histo", temp2);


                if(histo->GetEntries() >= limite && itermin == 0)
		{
                        *MinRange = sigma;
                        itermin = 1;
		}

                temp2.Clear();

                temp2 = temp;
                temp2 += " && mmg_s > ";
                temp2 += maxresult;

		chain->Draw("mmg_s>>histo2", temp2);


                if(histo2->GetEntries() >= limite && itermax == 0)
		{
		        *MaxRange =  lastX - sigma;
                        itermax = 1;
                }

		if (itermin == 1 && itermax == 1) sigma = 2.0;

		

		if((histo->GetEntries() < (limite * 0.1)) && (histo2->GetEntries() < (limite * 0.1)) && loop1 == 0)
		{
			 sigma +=0.05;
		}
		else loop1 = 1;
/*	
                if(loop1 == 1 && (histo->GetEntries() < (limite * 0.9)) && (histo2->GetEntries() < (limite * 0.9)))
		{
                        sigma +=0.01;
                }
		else if(loop1 == 1) loop2 = 1;

		if(loop1 == 1 && loop2 == 1 && (histo->GetEntries() > limite) || (histo2->GetEntries() > limite) && (histo->GetEntries() > (limite * 0.9)) || (histo2->GetEntries() > (limite * 0.9)))
                {
                         sigma +=0.01;
                }

 */
		

                histo->Delete();
                histo2->Delete();
        }


}

void SymetricRangeEstimator3(TChain * chain, double centralValue, double * MinRange, double * MaxRange, double Entries, double pourcentage, TString temp)
{
	

        double min = 0.0;
        double max = 0.0;
	int iterMin = 0;
	int iterMax = 0;
        string minresult;
        string maxresult;
        std::ostringstream oss;
        std::ostringstream oss2;
        TString temp2;


        for(double sigma = 0.001; sigma < 1.0; sigma += 0.001)
        //for(double sigma = 0.01; sigma < 1.0; sigma += 0.01)
	{
		temp2.Clear();
                oss.str("");
                oss2.str("");
                minresult.clear();
                maxresult.clear();

                min = centralValue - sigma; 
                max = centralValue + sigma;

                oss << min; 
                minresult = oss.str();
      
                oss2 << max; 
                maxresult = oss2.str(); 
      
		if(iterMin == 0)
		{
                	TH1D *histo = new TH1D("histo","histo", 200, 0.0, 2.0);
      
                	temp2 = temp;
                	temp2 += " && mmg_s > ";
                	temp2 += minresult;
                	temp2 += " && mmg_s < ";
                	temp2 += centralValue;

                	chain->Draw("mmg_s>>histo", temp2);
			if(histo->GetEntries() >= (Entries * (pourcentage / 2.0)))
                	{
                        	*MinRange = centralValue - sigma;
                        	iterMin = 1;
                	}
			histo->Delete();
		}

		if(iterMax == 0)
                {		

			TH1D *histo2 = new TH1D("histo2","histo2", 200, 0.0, 2.0);
                	
			temp2 = temp;
                	temp2 += " && mmg_s > ";
                	temp2 += centralValue;
                	temp2 += " && mmg_s < ";
                	temp2 += maxresult;
                
			chain->Draw("mmg_s>>histo2", temp2);            
	
			if(histo2->GetEntries() >= (Entries * (pourcentage / 2.0)))
                	{
                        	*MaxRange = centralValue + sigma;
                        	iterMax = 1;
                	}

                	histo2->Delete();
		}
		
		if(iterMin == 1 && iterMax == 1) sigma = 1.0;

        }
		if(iterMin == 0) *MinRange = 0.0;
		if(iterMax == 0) *MaxRange = 1.5;		

}


Double_t chiSquare(RooPlot* plot_, char* pdfname, char* histname, int nFitParam, double* JanChi2, double* DegreesOfFreedom, double* pValue, int* fewBins)
{
  // Calculate the chi^2/NDOF of this curve with respect to the histogram
  // 'hist' accounting nFitParam floating parameters in case the curve
  // was the result of a fit

  // Find curve object
  RooCurve* curve = (RooCurve*) plot_->findObject(pdfname, RooCurve::Class());
  //RooCurve* curve = plot_->getCurve(pdfname);  
  //curve->Print();

  if (!curve) {
    cout<<endl << "cit::RooChi2Calculator(plotname=" << plot_->GetName()
         << ")::chiSquare(..) cannot find curve" << endl ;
    return 0 ;
  }

  // Find histogram object
  RooHist* hist = (RooHist*) plot_->findObject(histname, RooHist::Class()) ;
  //RooHist* hist = plot_->getHist(histname);

  if (!hist) {
    cout<<endl << "cit::RooChi2Calculator(plotname=" << plot_->GetName()
         << ")::chiSquare(..) cannot find histogram" << endl ;
    return 0 ;
  }


  Int_t i,np = hist->GetN() ;
  Double_t x,y,/*eyl,eyh,*/ xl,xh ;

  // Find starting and ending bin of histogram based on range of RooCurve
  Double_t xstart,xstop ;

#if ROOT_VERSION_CODE >= ROOT_VERSION(4,0,1)
  curve->GetPoint(0,xstart,y) ;
  curve->GetPoint(curve->GetN()-1,xstop,y) ;
#else
  const_cast<RooCurve*>(curve)->GetPoint(0,xstart,y) ;
  const_cast<RooCurve*>(curve)->GetPoint(curve->GetN() - 1,xstop,y) ;
#endif

  Int_t nbin(0) ;

  Double_t chisq(0) ;
  for (i=0 ; i<np ; i++) {   

    // Retrieve histogram contents
    hist->GetPoint(i,x,y) ;
    xl = x - hist->GetEXlow()[i] ;
    xh = x + hist->GetEXhigh()[i] ;
    // eyl = hist->GetEYlow()[i] ;
    // eyh = hist->GetEYhigh()[i] ;

    // Check if the whole bin is in range of curve
    if (xl < xstart || xstop < xh) continue ;

    //if(y != 0 && y < 35.0)
    if(y == 0)
    {
    	cout<<endl<<"Trop peu d'entree : "<<y<<" dans le bin : "<<i<<"  >>>Need to reduce the binning for the p-value calculation!"<<endl;
	*fewBins = 1;
	break;
	
    }
    else *fewBins = 0;

    nbin++ ;

    // Integrate function over this bin.
    // Start a hack to work around a bug in RooCurve::interpolate
    // that sometimes gives a wrong result.
    Double_t avg = curve->average(xl, xh);
    Double_t avg2 = 0.5 * (curve->average(xl, x) + curve->average(x, xh));
    if (avg + avg2 > 0 &&
	(avg2 - avg) / (avg2 + avg) > 0.1) {
      avg = curve->interpolate(x);
    }
    // End of hack around the bug in RooCurve::interpolate

    // JV: Adjust observed and expected number of events for bin width to represent
    // number of events.
    Double_t norm = (xh - xl) / plot_->getFitRangeBinW();
    y *= norm;
    avg *= norm;

    if (avg < 5.) {
      cout << "cit::RooChi2Calculator(plotname=" << plot_->GetName()
			    << ")::chiSquare(..) expectation in bin "
			    << i << " is " << avg << " < 5!" << endl ;
    }

    // JV: Use the expected number of events for the y uncertainty,
    // See (33.34) of http://pdg.lbl.gov/2011/reviews/rpp2011-rev-statistics.pdf

    // Add pull^2 to chisq
    if (avg != 0) {      
      Double_t resid = y - avg;
      chisq += (resid * resid / avg) ;
    }
  }

  // Return chisq/nDOF 
  *JanChi2 = chisq / (nbin - nFitParam);
  *DegreesOfFreedom = (nbin - nFitParam);
  *pValue =  TMath::Prob(chisq, nbin - nFitParam);

  return chisq / (nbin - nFitParam) ;
}

  
RooHist* residHist(RooPlot* plot_, char *histname, char* curvename, bool normalize, string dossierSauvegardePull, int iteration)
{
  if(normalize == true)
  {
        FILE *fPullsX = fopen(Form("%sPullsX%d.txt",dossierSauvegardePull.c_str(), iteration), "w" );
        delete fPullsX;
        FILE *fPullsErrorX = fopen(Form("%sPullsErrorX%d.txt",dossierSauvegardePull.c_str(), iteration), "w" );
        delete fPullsErrorX;
        FILE *fPullsY = fopen(Form("%sPullsY%d.txt",dossierSauvegardePull.c_str(), iteration), "w" );
        delete fPullsY;
        FILE *fPullsErrorY = fopen(Form("%sPullsErrorY%d.txt",dossierSauvegardePull.c_str(), iteration), "w" );
        delete fPullsErrorY;
  }

  // Create and return RooHist containing  residuals w.r.t to given curve->
  // If normalize is true, the residuals are normalized by the histogram
  // errors creating a RooHist with pull values


  // Find curve object
  RooCurve* curve = (RooCurve*) plot_->findObject(curvename, RooCurve::Class());
  if (!curve) {
    cout << "cit::RooChi2Calculator(plotname=" << plot_->GetName()
         << ")::residHist(..) cannot find curve" << endl ;
    return 0 ;
  }

  // Find histogram object
  RooHist* hist = (RooHist*) plot_->findObject(histname, RooHist::Class()) ;
  if (!hist) {
    cout << "cit::RooChi2Calculator(plotname=" << plot_->GetName()
         << ")::residHist(..) cannot find histogram" << endl ;
    return 0 ;
  }


  // Copy all non-content properties from hist
  RooHist* ret = new RooHist(plot_->getFitRangeBinW()) ;
  if (normalize) {
    ret->SetName(Form("pull_%s_%s", hist->GetName(), curve->GetName())) ;
    ret->SetTitle(Form("Pull of %s and %s", hist->GetTitle(), curve->GetTitle())) ;
  } else {
    ret->SetName(Form("resid_%s_%s", hist->GetName(), curve->GetName())) ;
    ret->SetTitle(Form("Residual of %s and %s", hist->GetTitle(), curve->GetTitle())) ;
  }

  // Determine range of curve
  Double_t xstart, xstop, y ;
#if ROOT_VERSION_CODE >= ROOT_VERSION(4,0,1)
  curve->GetPoint(0, xstart, y) ;
  curve->GetPoint(curve->GetN()-1,xstop,y) ;
#else
  const_cast<RooCurve&>(curve)->GetPoint(0, xstart, y) ;
  const_cast<RooCurve&>(curve)->GetPoint(curve->GetN()-1, xstop, y) ;
#endif
  // cout << "cit::RooChi2Calculator::residHist dumping curve:\n";
  // for (int i=0; i<curve->GetN(); ++i){
  //   Double_t xi, yi;
  //   curve->GetPoint(i, xi, yi);
  //   printf("i=%d x,y: %.3g, %.3g\n", i, xi, yi);
  // }

  // cout << "cit::RooChi2Calculator::residHist  adding bins with error:\n";

  // Add histograms, calculate Poisson confidence interval on sum value
  for(Int_t i=0 ; i < hist->GetN() ; i++) {
    Double_t x, point;
#if ROOT_VERSION_CODE >= ROOT_VERSION(4,0,1)
    hist->GetPoint(i,x,point) ;
#else
    const_cast<RooHist&>(hist)->GetPoint(i,x,point) ;
#endif
    Double_t xl = x - hist->GetErrorXlow(i);
    Double_t xh = x + hist->GetErrorXhigh(i);

    // Only calculate pull for bins inside curve range
    if (xl < xstart || xstop < xh) continue ;

    Double_t norm = (xh - xl) / plot_->getFitRangeBinW();
    point *= norm;

    // Start a hack to work around a bug in RooCurve::interpolate
    // that sometimes gives a wrong result.
    Double_t avg = curve->average(xl, xh);
    Double_t avg2 = 0.5 * (curve->average(xl, x) + curve->average(x, xh));
    Double_t yexpected;
    if (avg + avg2 > 0 && (avg2 - avg) / (avg2 + avg) > 0.1) {
      yexpected = curve->interpolate(x);
    } else {
      yexpected = avg;
    }
    // End of hack around the bug in RooCurve::interpolate

    // Correct the expected number of events in this bin for the non-uniform
    // bin width.
    yexpected *= norm;

    Double_t yy = point - yexpected;
    // Normalize to the number of events per bin taking into account
    // variable bin width.
    Double_t dy = TMath::Sqrt(yexpected);
    if (normalize) {
	if (dy==0.) {
	  cout << "cit::RooChi2Calculator::residHist(histname ="
               << hist->GetName() << ", ...) WARNING: point "
               << i << " has zero error, setting residual to zero"
               << endl ;
	  yy=0 ;
	  dy=0 ;
	} else {
	  yy /= dy;
	  dy = 1.;
	}
    }
    // printf("bin=%3d n=%5.3g nu=%5.3g x=%5.3g .. %5.3g y=%5.3g +/- %5.3g "
    //	   "norm=%5.3g\n", i, point, yexpected, xl, xh, yy, dy, norm);
    ret->addBinWithError(x,yy,dy,dy);
  
    if(normalize == true)
    {
        ofstream monFluxPullsX(Form("%sPullsX%d.txt",dossierSauvegardePull.c_str(), iteration), ios::app);
        monFluxPullsX << yy <<endl;
        monFluxPullsX.close();

        ofstream monFluxPullsErrorX(Form("%sPullsErrorX%d.txt",dossierSauvegardePull.c_str(), iteration), ios::app);
        monFluxPullsErrorX << dy <<endl;
        monFluxPullsErrorX.close();

        ofstream monFluxPullsY(Form("%sPullsY%d.txt",dossierSauvegardePull.c_str(), iteration), ios::app);
        monFluxPullsY << x <<endl;
        monFluxPullsY.close();

        ofstream monFluxPullsErrorY(Form("%sPullsErrorY%d.txt",dossierSauvegardePull.c_str(), iteration), ios::app);
        monFluxPullsErrorY << dy <<endl;
        monFluxPullsErrorY.close();

    }

  }
  return ret;
}


int NbLignesFichier(string fichier)
{
        ifstream in(fichier.c_str()); //Ouverture en mode lecture de fichier

        string ligne; //Création d'une chaine de caractere
        int nbLignes = 0;

        while(std::getline(in, ligne)) nbLignes++;

        in.close(); //On ferme le fichier

        return nbLignes;
}




int main(int argc, char *argv[])
{

	cout << "argc= " << argc << endl;
        for(int iarg = 0 ; iarg < argc; iarg++)
        {
                cout << "argv[" << iarg << "]= " << argv[iarg] << endl;
        }

	if( argc == 1 )
	{
		cerr << "arguments should be passed !! EndCaps, r9sup, SetOfCorrections, nomFitMethode, PileUpVersion, LowMmumuLim, HightMmumuLim, FitPourcentage, phiCracks, etaCracks, variableX, isMC" <<endl;
		return 1;

	}
	
	int EndCaps = 0;//0 1 (0 = Barrel, 1 = Endcaps)
	int r9sup = 2;//0 1 2 (0 = low r9, 1 = high r9, 2 = no r9 cuts)
	string SetOfCorrections = "Regression";//"ElectronTunedCorrections" "ETHZCorrections" "Regression" "DefaultCorrections"
	string nomFitMethode = "roofit_bwcb_fft2"; //"RooBifurcatedGauss"; //RooVoigtian2 //Cruijff
	string PileUpVersion = "2011";
	string LowMmumuLim = "40";
	string HightMmumuLim = "80";
	int FitPourcentage = 90;
	bool phiCracks = true;
	bool etaCracks = true;
	string variableX = "Photon_Et"; 
	int isMC = 0; //0 1
	string MZbinning = "5GeV"; //"05GeV" "1GeV" "2GeV" "5GeV" 
	string MuonCorrection = "Rochester"; //"Rochester" "NoMuonCorrection"
	string SurfaceMethod = "Profile_Surface"; //"Profile_Surface" "Fitted_Surface"
	string Category = "OneBin"; //"Vgamma24" "Vgamma8" "OneBin"
	string date = "rereco"; //"16Jan" "30Nov" "Vgamma" "Fall11"

	if( argc > 1 ) 
	{
		std::stringstream ss ( argv[1] );
                ss >> EndCaps;
	}
	if( argc > 2 ) 
	{
		std::stringstream ss ( argv[2] );
                ss >> r9sup;
	}
	if( argc > 3 ) 
	{
		SetOfCorrections = argv[3];
	}
	if( argc > 4 ) 
	{
		nomFitMethode = argv[4];
	}
	if( argc > 5 ) 
	{
		PileUpVersion = argv[5];
	}
	if( argc > 6 ) 
	{
		LowMmumuLim = argv[6];
	}
	if( argc > 7 ) 
	{
		HightMmumuLim = argv[7];
	}
	if( argc > 8 ) 
	{
		std::stringstream ss ( argv[8] );
                ss >> FitPourcentage;
	}
	if( argc > 9 ) 
	{
		phiCracks = argv[9];	
	}
	if( argc > 10 ) 
	{
		etaCracks = argv[10];
	}
	if( argc > 11 ) 
	{
		variableX = argv[11];
	}
	if( argc > 12 ) 
	{
		std::stringstream ss ( argv[12] );
                ss >> isMC;
	}
	if( argc > 13 ) 
        {
        	MZbinning = argv[13];   
        }
	if( argc > 14 ) 
        {
                MuonCorrection = argv[14];   
        }
	if( argc > 15 )
        {
                SurfaceMethod = argv[15];
        }	
	if( argc > 16 ) 
        {
                Category = argv[16];
        }	
	if( argc > 17 ) 
        {
                date = argv[17];
        }		


	//gSystem->Load("../../src/libToto.so");

	//gROOT->SetStyle("Plain");	
	gROOT->Reset();
	TGaxis::SetMaxDigits(3);
	setTDRStyle();
	//CMSstyle();
	//TCanvas *c1 = new TCanvas("c1", "c1",0,0,600,600);
	string isMCChain = "";
	if(isMC == 0) isMCChain = "Data";
	if(isMC == 1) isMCChain = "MC";	

	string compressSetOfCorrectionsChain = "";
        if(SetOfCorrections == "ETHZCorrections") compressSetOfCorrectionsChain = "ETHZ";
        //if(SetOfCorrections == "ElectronTunedCorrections") compressSetOfCorrectionsChain = "START42_V11";
	if(SetOfCorrections == "ElectronTunedCorrections") compressSetOfCorrectionsChain = "1.0";
	if(SetOfCorrections == "Regression") compressSetOfCorrectionsChain = "Regression";

	double FitPourcentageD = FitPourcentage / 100.0; 
/*
	string DossierRangeVersion = "";
        if(RangeVersion == "Louis") DossierRangeVersion = "Olivier21";
        if(RangeVersion == "Jan") DossierRangeVersion = "Olivier20";
*/
        string DossierCracks = "";
        if(phiCracks == true && etaCracks == true) DossierCracks = "WithCracks";
        if(phiCracks == false && etaCracks == false) DossierCracks = "WithoutCracks";
        if(phiCracks == true && etaCracks == false) DossierCracks = "WithoutEtaCracks";
        if(phiCracks == false && etaCracks == true) DossierCracks = "WithoutPhiCracks";


	double MeanTab[100] = {0.0};
	double MeanErrorTab[100] = {0.0};
	double SigmaRTab[100] = {0.0};
        double SigmaLTab[100] = {0.0};
	double SigmaTab[100] = {0.0};
	double SigmaErrorTab[100] = {0.0};
	double SigmaEffTab[100] = {0.0};
	double MinVar[100] = {0.0};
	double MaxVar[100] = {0.0};
	double xValue[100] = {0.0};
        double xErrorR[100] = {0.0};
	double xErrorL[100] = {0.0};
	//Chi2
	double yValue[100] = {0};
        double xErrorRight[100] = {0};
        double xErrorLeft[100] = {0};
	
	double ChiSquareTab[100] = {0.0};

	double JanChiSquareTab[100] = {0.0};
        double PValueTab[100] = {0.0};
        double PValueErrorTab[100] = {0.0};
	double DegreesOfFreedomTab[100] = {0.0};
	
	double minLogLikelihoodTab[100] = {0.0};
	double differenceLogLikelihoodTab[100] = {0.0};

	int entriesTot = 0;
	int anti2 = 0;
	int n = 0;
        double Mean = 0;
        double Sigma = 0;

	int fitRoo = 1;
/*
	////////// CONFIG //////////

	
	//int EndCaps = 0;  //On est dans le Barrel
	int EndCaps = 1; //On est dans End Caps

	//int r9sup = 1; // r9>0.94(0.95)
	int r9sup = 0; // r9<0.94(0.95)

	//int fitRoo = 0; // Use minuit fit
        int fitRoo = 1; // Use Roofit
	
	if(EndCaps == 0) n = 6; // Number of bins  
	if(EndCaps == 1) n = 6; // Number of bins

	// Roofit Fit Fonctions //

	//string nomFitMethode = "RooCrystalBall";
	//string nomFitMethode = "RooLogNormal";
	//string nomFitMethode = "RooGauss";
	//string nomFitMethode = "RooGamma";
	string nomFitMethode = "RooBifurcatedGauss";
	//string nomFitMethode = "RooSumGauss";
	//string nomFitMethode = "RooGenericPDF";
 	
	//

	/// PileUpVersion ///

	//string PileUpVersion = "May10";
	//string PileUpVersion = "Promptv4";
	//string PileUpVersion = "July05";
	//string PileUpVersion = "Aug05";	
	//string PileUpVersion = "Oct03";
	//string PileUpVersion = "2011A";
	//string PileUpVersion = "2011B";
	string PileUpVersion = "2011"; 

	///
*/
	double pourcentage = 0.78;


	// Plots limits //
	
	double xminChi2 , xminSigmaLSigmaR , xminSigma , xminSigmaTg , xminSigmaEff , xminSigmaEffTg , xminmmg_s , xminJanChi2 , xminPValue;
        xminChi2 = xminSigmaLSigmaR = xminSigma = xminSigmaTg = xminSigmaEff = xminSigmaEffTg = xminmmg_s = xminJanChi2 = xminPValue = 0.0; 

	if(variableX == "Photon_SC_Eta") xminSigmaLSigmaR = xminSigma = xminSigmaTg = xminSigmaEff = xminSigmaEffTg = xminmmg_s = xminJanChi2 = xminPValue = -3.0;

        double xmaxChi2 , xmaxSigmaLSigmaR , xmaxSigma , xmaxSigmaTg , xmaxSigmaEff , xmaxSigmaEffTg , xmaxmmg_s , xmaxJanChi2 , xmaxPValue;

        if(variableX == "Photon_SC_rawEt") xmaxChi2 = xmaxSigmaLSigmaR = xmaxSigma = xmaxSigmaTg = xmaxSigmaEff = xmaxSigmaEffTg = xmaxmmg_s = xmaxJanChi2 = xmaxPValue = 250.0;
        if(variableX == "Photon_Et") xmaxChi2 = xmaxSigmaLSigmaR = xmaxSigma = xmaxSigmaTg = xmaxSigmaEff = xmaxSigmaEffTg = xmaxmmg_s = xmaxJanChi2 = xmaxPValue = 250.0;
        if(variableX == "Photon_E") xmaxChi2 = xmaxSigmaLSigmaR = xmaxSigma = xmaxSigmaTg = xmaxSigmaEff = xmaxSigmaEffTg = xmaxmmg_s = xmaxJanChi2 = xmaxPValue = 1000.0;
        if(variableX == "Photon_SC_Eta") xmaxChi2 = xmaxSigmaLSigmaR = xmaxSigma = xmaxSigmaTg = xmaxSigmaEff = xmaxSigmaEffTg = xmaxmmg_s = xmaxJanChi2 = xmaxPValue = 3.0;
        if(variableX == "Photon_SC_brem") xmaxChi2 = xmaxSigmaLSigmaR = xmaxSigma = xmaxSigmaTg = xmaxSigmaEff = xmaxSigmaEffTg = xmaxmmg_s = xmaxJanChi2 = xmaxPValue = 15.0;


        double yminChi2 , yminSigmaLSigmaR , yminSigma , yminSigmaTg , yminSigmaEff , yminSigmaEffTg , yminmmg_s, yminJanChi2, yminPValue;
        yminChi2 = 0.0; 
        yminSigmaLSigmaR = -0.15;
        yminSigma = -15; 
        yminSigmaTg = 0.0; 
        yminSigmaEff = -0.15;
        yminSigmaEffTg = 0.0; 
        yminJanChi2 = 0.0; 
        yminPValue = 0.0;     
        if(r9sup == 2) yminmmg_s = -15.0;
	if(r9sup == 1) yminmmg_s = -15.0;
        if(r9sup == 0) yminmmg_s = -15.0;
	

        double ymaxChi2 , ymaxSigmaLSigmaR , ymaxSigma , ymaxSigmaTg , ymaxSigmaEff , ymaxSigmaEffTg , ymaxmmg_s, ymaxJanChi2, ymaxPValue;
        ymaxChi2 = 15.0;
        ymaxSigmaLSigmaR = 0.15;
        ymaxSigma = 15;
        ymaxSigmaTg = 10.0;
        ymaxSigmaEff = 0.15;
        ymaxSigmaEffTg = 0.1; 
        ymaxJanChi2 = 15.0;
        ymaxPValue = 1.0;
	if(r9sup == 2) ymaxmmg_s = 15.0;     
        if(r9sup == 1) ymaxmmg_s = 15.0; 
        if(r9sup == 0) ymaxmmg_s = 15.0;
	

	//

	////////// END CONFIG //////////


	double MinRange = 0; 
	double MaxRange = 0;
	double maxDistri = 0;
	double lastX = 0;
	double Entries = 0;
	double mean_value = 0;
        double mean_error = 0;
        double sigma_value = 0;
	double sigma_value_error = 0;
	Double_t sigmaEff_value = 0;
	double sigmaR_value = 0;
	double sigmaL_value = 0;
        double ChiSquare = 0;
	double params[5] = {0};
	double JanChi2 = 0;
        double DegreesOfFreedom = 0;
        double pValue = 0;
	double minLogLikelihood = 0;
	double differenceLogLikelihood = 0;

	double mean = 0;
        double rms = 0;
	double meanMC = 0.0;

	double phiCrackSize = (double)(21.5) / (double) (1290.0);
        double phiCrackPosition = (double)(TMath::Pi()) / (double)(9.0);
        double phiOffset = -(double)(10.0 * TMath::Pi()) / (double)(180.0);

	string EndcapsChain;

	cout<<endl<<"variableX = "<<variableX<<", Category = "<<Category<<endl;

	if(variableX == "Photon_SC_rawEt") EndcapsChain = "LimitesAllEtRaw.txt";
        if(variableX == "Photon_Et" && Category == "Vgamma24") EndcapsChain = "LimitesAllPtVgamma24.txt";
	if(variableX == "Photon_Et" && Category == "Vgamma8") EndcapsChain = "LimitesAllPtVgamma8.txt";
	if(variableX == "Photon_Et" && Category == "OneBin") EndcapsChain = "LimitesAllPtOneBin.txt";
        if(variableX == "Photon_E" && EndCaps == 0) EndcapsChain = "LimitesEnergyBARREL.txt";
	if(variableX == "Photon_E" && EndCaps == 1) EndcapsChain = "LimitesEnergyENDCAPS.txt";
        if(variableX == "Photon_SC_Eta" && EndCaps == 0) EndcapsChain = "LimitesAllEtaBARREL.txt";
        if(variableX == "Photon_SC_Eta" && EndCaps == 1) EndcapsChain = "LimitesAllEtaENDCAPS.txt";
        if(variableX == "Photon_SC_brem") EndcapsChain = "LimitesAllBrem.txt";

	

	cout<<"EndcapsChain = "<<EndcapsChain<<endl;
	n = NbLignesFichier(EndcapsChain.c_str()) - 1;
	ifstream monFlux(EndcapsChain.c_str());

	double nombre;
	
	string TestChain_v1 = "";
	//TestChain_v1 = Form("miniTree_DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_PU_S6_v2_%s_%s_%s_%s_%s_1_v3_partALL.root",compressSetOfCorrectionsChain.c_str(),MuonCorrection.c_str(),LowMmumuLim.c_str(),HightMmumuLim.c_str(),MZbinning.c_str()); 

	//cout<<endl<<"TestChain_v1 = "<<TestChain_v1<<endl<<endl;

	TChain * DataChain = new TChain("miniTree");
	
	if(isMC == 0)
	{
        	DataChain->Add("miniTree_2011A_03Oct2011V1ReReco_toto_v2_Triangle_0_Regression_Rochester_v2_partALL.root");
		DataChain->Add("miniTree_2011A_PromptSkimV5ReReco_toto_v2_Triangle_0_Regression_Rochester_v2_partALL.root");
		DataChain->Add("miniTree_2011A_05Jul2011ReReco_toto_v2_Triangle_0_Regression_Rochester_v2_partALL.root");
		DataChain->Add("miniTree_2011B_PromptSkimV1ReReco_toto_v2_Triangle_0_Regression_Rochester_v2_partALL.root");
	}

	if(isMC == 1)
        {
                DataChain->Add("miniTree_DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Fall11-PU_S6_START42_V14B-v1_September12_Triangle_1_Regression_Rochester_v2_partALL.root");
                DataChain->Add("miniTree_DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Fall11-PU_S6_START42_V14B-v1_September12_Triangle_2_Regression_Rochester_v2_partALL.root");
                DataChain->Add("miniTree_TTJets_TuneZ2_7TeV-madgraph-tauola_Triangle_3_Regression_Rochester_v2_partALL.root");
                DataChain->Add("miniTree_WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_Triangle_3_Regression_Rochester_v2_partALL.root");
        }

	TString temp = "";
	TString tempVarChain = "";

	string nomDossier = Form("Results_November_2012_v1_Rochv4_NewSelection_v3/%s/MmumuSelection%s_%s/%s/%s/%s/%s/MZbinning_%s/%s/%s/%dPourcents/",variableX.c_str(),LowMmumuLim.c_str(),HightMmumuLim.c_str(),DossierCracks.c_str(),SetOfCorrections.c_str(),MuonCorrection.c_str(),isMCChain.c_str(),MZbinning.c_str(),SurfaceMethod.c_str(),Category.c_str(),FitPourcentage);

	cout<<endl<<"nomDossier = "<<nomDossier<<endl;
	

	string nomDossierFits = "";
	string nomDossierTGraph = ""; 
	string nomFichier = "";	

	for(int j = 0; j < n; j++) //mettre j<n
        {
		TCanvas *c2 = new TCanvas("c2", "c2",0,0,600,600);

/*		//TCanvas *c1 = new TCanvas("c1", "c1",0,0,600,600);
                cout<<endl<<"j + 1 = "<<j+1<<endl;
		
		TFile file1("HistoBin/Histogram_172EB_162EE_r9_4.root","UPDATE");


                //c1->Clear();
		//c1->Close();


		if(EndCaps == 0) TH1D *MC = (TH1D*)file1.Get(Form("ErawCEtaEB_%d",j));
		if(EndCaps == 1) TH1D *MC = (TH1D*)file1.Get(Form("ErawCEtaEE_%d",j));

		if(EndCaps == 0) TH1D *EtaTH1 = (TH1D*)file1.Get(Form("EtaEB_%d",j));
		if(EndCaps == 1) TH1D *EtaTH1 = (TH1D*)file1.Get(Form("EtaEE_%d",j));

		xValue[j] = VarTH1->GetMean(1);
*/

                cout<<endl<<"////////// ----- Bin "<<j+1<<" sur "<<n<<" ----- //////////"<<endl<<endl;
	/*	
		TChain * DataChain = new TChain("photons");
        	//DataChain->Add("/tmp/lsgandur/theMiniTree.root");
		
		//DataChain->Add("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/DeuxGammas/theMiniTree.root");

		DataChain->Add("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/DeuxGammas/theMiniTree.root");
*/
		TH1D *MC = new TH1D("MC","Data Z#rightarrow#mu#mu#gamma", 1000, -0.5, 0.5);

                temp.Clear();
		tempVarChain.Clear();

		temp += "isJanLooseMMG == 1";

		if(EndCaps == 0 && r9sup == 1 && variableX != "Photon_SC_Eta") temp += " && (abs(Photon_SC_Eta)) < 1.4442 && Photon_r9 > 0.94";
                if(EndCaps == 0 && r9sup == 0 && variableX != "Photon_SC_Eta") temp += " && (abs(Photon_SC_Eta)) < 1.4442 && Photon_r9 < 0.94";
                if(EndCaps == 1 && r9sup == 1 && variableX != "Photon_SC_Eta") temp += " && (abs(Photon_SC_Eta)) > 1.566 && Photon_r9 > 0.95";
                if(EndCaps == 1 && r9sup == 0 && variableX != "Photon_SC_Eta") temp += " && (abs(Photon_SC_Eta)) > 1.566 && Photon_r9 < 0.95";

		if(EndCaps == 0 && r9sup == 2 && variableX != "Photon_SC_Eta") temp += " && (abs(Photon_SC_Eta)) < 1.4442";
                if(EndCaps == 1 && r9sup == 2 && variableX != "Photon_SC_Eta") temp += " && (abs(Photon_SC_Eta)) > 1.566";

                if(EndCaps == 0 && r9sup == 1 && variableX == "Photon_SC_Eta") temp += " && Photon_r9 > 0.94";
                if(EndCaps == 0 && r9sup == 0 && variableX == "Photon_SC_Eta") temp += " && Photon_r9 < 0.94";
                if(EndCaps == 1 && r9sup == 1 && variableX == "Photon_SC_Eta") temp += " && Photon_r9 > 0.95";
                if(EndCaps == 1 && r9sup == 0 && variableX == "Photon_SC_Eta") temp += " && Photon_r9 < 0.95";

                if(EndCaps == 0 && etaCracks == false) temp += " && ( (abs(Photon_SC_Eta)) > 0.018 || ( (abs(Photon_SC_Eta)) < 0.423 && (abs(Photon_SC_Eta)) > 0.461 ) || ( (abs(Photon_SC_Eta)) < 0.770 && (abs(Photon_SC_Eta)) > 0.806 ) || ( (abs(Photon_SC_Eta)) < 1.127 && (abs(Photon_SC_Eta)) > 1.163 ) )";

                if(EndCaps == 0 && phiCracks == false)
                {
                        temp += Form(" && ((abs(Photon_SC_Phi) + (%f)) < (%f) ||",phiOffset,phiCrackSize);
                        temp += Form(" ( (1.0*(%f) - (%f)) < (abs(Photon_SC_Phi) + (%f)) && (abs(Photon_SC_Phi) + (%f)) < (1.0*(%f) + (%f)) ) ||",phiCrackPosition,phiCrackSize,phiOffset,phiOffset,phiCrackPosition,phiCrackSize);
                        temp += Form(" ( (2.0*(%f) - (%f)) < (abs(Photon_SC_Phi) + (%f)) && (abs(Photon_SC_Phi) + (%f)) < (2.0*(%f) + (%f)) ) ||",phiCrackPosition,phiCrackSize,phiOffset,phiOffset,phiCrackPosition,phiCrackSize);
                        temp += Form(" ( (3.0*(%f) - (%f)) < (abs(Photon_SC_Phi) + (%f)) && (abs(Photon_SC_Phi) + (%f)) < (3.0*(%f) + (%f)) ) ||",phiCrackPosition,phiCrackSize,phiOffset,phiOffset,phiCrackPosition,phiCrackSize);
                        temp += Form(" ( (4.0*(%f) - (%f)) < (abs(Photon_SC_Phi) + (%f)) && (abs(Photon_SC_Phi) + (%f)) < (4.0*(%f) + (%f)) ) ||",phiCrackPosition,phiCrackSize,phiOffset,phiOffset,phiCrackPosition,phiCrackSize);
                        temp += Form(" ( (5.0*(%f) - (%f)) < (abs(Photon_SC_Phi) + (%f)) && (abs(Photon_SC_Phi) + (%f)) < (5.0*(%f) + (%f)) ) ||",phiCrackPosition,phiCrackSize,phiOffset,phiOffset,phiCrackPosition,phiCrackSize);
                        temp += Form(" ( (6.0*(%f) - (%f)) < (abs(Photon_SC_Phi) + (%f)) && (abs(Photon_SC_Phi) + (%f)) < (6.0*(%f) + (%f)) ) ||",phiCrackPosition,phiCrackSize,phiOffset,phiOffset,phiCrackPosition,phiCrackSize);
                        temp += Form(" ( (7.0*(%f) - (%f)) < (abs(Photon_SC_Phi) + (%f)) && (abs(Photon_SC_Phi) + (%f)) < (7.0*(%f) + (%f)) ) ||",phiCrackPosition,phiCrackSize,phiOffset,phiOffset,phiCrackPosition,phiCrackSize);
                        temp += Form(" ( (8.0*(%f) - (%f)) < (abs(Photon_SC_Phi) + (%f)) && (abs(Photon_SC_Phi) + (%f)) < (8.0*(%f) + (%f)) ) ||",phiCrackPosition,phiCrackSize,phiOffset,phiOffset,phiCrackPosition,phiCrackSize);
                        temp += Form(" ( (9.0*(%f) - (%f)) < (abs(Photon_SC_Phi) + (%f)) ) )",phiCrackPosition,phiCrackSize,phiOffset);

                }

		temp +=" && " + variableX + " > ";
                //temp += " && Photon_Et > ";

		if(j == 0)
		{
			monFlux >> nombre;
			MinVar[j] = nombre;
			monFlux >> nombre;
                	MaxVar[j] = nombre;
		}
	
		if(j > 0)
		{
			MinVar[j] = MaxVar[j-1];
			monFlux >> nombre;
			MaxVar[j] = nombre;

		}

		temp += MinVar[j];
		temp += " && " + variableX + " <= ";
                //temp += " && Photon_Et < ";
		temp += MaxVar[j];


		cout<<endl<<"temp = "<<temp<<endl;

		DataChain->Draw("mmg_s>>MC",temp);
		c2->Print("alors.png");
		cout<<"MC->GetMaximum= "<<MC->GetMaximumBin() * 0.001 + 0.5<<endl;	
		maxDistri = MC->GetMaximumBin() * 0.001 - 0.5;
		lastX = 0.5;
		Entries = MC->GetEntries();		
		meanMC = MC->GetMean();
		cout<<"MC->GetMean() = "<<MC->GetMean()<<endl;
		//maxDistri = MC->GetMean();

		TH1D *VarXvalue = 0; //= new TH1D();
		if(variableX == "Photon_SC_rawEt") VarXvalue = new TH1D("VarXvalue","Var Z#rightarrow#mu#mu#gamma", 1000, 0, 250.0);
                if(variableX == "Photon_Et") VarXvalue = new TH1D("VarXvalue","Var Z#rightarrow#mu#mu#gamma", 1000, 0, 250.0);
                if(variableX == "Photon_E") VarXvalue = new TH1D("VarXvalue","Var Z#rightarrow#mu#mu#gamma", 1000, 0, 1000.0);
                if(variableX == "Photon_SC_Eta") VarXvalue = new TH1D("VarXvalue","Var Z#rightarrow#mu#mu#gamma", 1000, -3.0, 3.0);
                if(variableX == "Photon_SC_brem") VarXvalue = new TH1D("VarXvalue","Var Z#rightarrow#mu#mu#gamma", 1000, 0, 15.0);

		tempVarChain = variableX + ">>VarXvalue";
                DataChain->Draw(tempVarChain,temp);
                //DataChain->Draw("Photon_Et>>VarXvalue",temp);

		xValue[j] = VarXvalue->GetMean(1);

		cout<<endl<<"xValue = "<<xValue[j]<<endl;

		c2->Clear();
		
		pourcentage = FitPourcentageD;


		cout<<endl<<"POURCENTAGE = "<<pourcentage<<endl;

		//RangeEstimator(pourcentage, maxDistri, DataChain, EndCaps, &MinRange, &MaxRange, MinVar[j], MaxVar[j]);
	
		//RangeEstimator2(pourcentage, DataChain, EndCaps, &MinRange, &MaxRange, MinVar[j], MaxVar[j]);	
		
		RangeEstimator3(pourcentage, DataChain, temp, EndCaps, &MinRange, &MaxRange);

		//SymetricRangeEstimator(DataChain, maxDistri, &MinRange, &MaxRange, Entries, pourcentage, temp);	
		//SymetricRangeEstimator2(DataChain, lastX, &MinRange, &MaxRange, Entries, pourcentage, temp);
		
		//SymetricRangeEstimator3(DataChain, meanMC, &MinRange, &MaxRange, Entries, pourcentage, temp);	

		//SymetricRangeEstimator3(DataChain, maxDistri, &MinRange, &MaxRange, Entries, pourcentage, temp);


		//MinRange = -0.5;
		//MaxRange = 0.5;

		cout<<endl<<"MinRange = "<<MinRange<<endl;
		cout<<endl<<"MaxRange = "<<MaxRange<<endl;

			

		if(fitRoo == 1)
		{
			mean = MC->GetMean();
                	rms = MC->GetRMS();
			TF1 * f = new TF1();
        		RooPlot* Erawframe = new RooPlot(-1.0,1.0);
			nomFichier = "mmg_s";
			nomDossierFits = nomDossier;

			if(nomFitMethode == "RooCrystalBall") 
			{
				if(r9sup == 2) nomDossierFits += "Fits/mmg_sCBr9All";
				if(r9sup == 1) nomDossierFits += "Fits/mmg_sCBr9sup";
				if(r9sup == 0) nomDossierFits += "Fits/mmg_sCBr9inf";
				RooCrystalBall(&mean_value, &mean_error, &sigmaEff_value, &sigma_value, &sigma_value_error, &sigmaR_value, &sigmaL_value, &ChiSquare, MC, j, EndCaps, temp, DataChain, mean, rms, c2, Erawframe, f, nomDossierFits, nomFichier, MinRange, MaxRange, &JanChi2, &DegreesOfFreedom, &pValue, r9sup, MinVar[j], MaxVar[j], variableX, isMC);
			}
			if(nomFitMethode == "RooLogNormal") 
                        {
				if(r9sup == 2) nomDossierFits += "Fits/mmg_sLNr9All";
                                if(r9sup == 1) nomDossierFits += "Fits/mmg_sLNr9sup";
				if(r9sup == 0) nomDossierFits += "Fits/mmg_sLNr9inf";
				RooLogNormal(&mean_value, &mean_error, &sigmaEff_value, &sigma_value, &sigma_value_error, &sigmaR_value, &sigmaL_value, &ChiSquare, MC, j, EndCaps, temp, DataChain, mean, rms, c2, Erawframe, f, nomDossierFits, nomFichier, MinRange, MaxRange, &JanChi2, &DegreesOfFreedom, &pValue, r9sup, MinVar[j], MaxVar[j], variableX, isMC);
			}
			if(nomFitMethode == "RooGauss")   
                        {
				if(r9sup == 2) nomDossierFits += "Fits/mmg_sGr9All";
                                if(r9sup == 1) nomDossierFits += "Fits/mmg_sGr9sup";
				if(r9sup == 0) nomDossierFits += "Fits/mmg_sGr9inf";
                                RooGauss(&mean_value, &mean_error, &sigmaEff_value, &sigma_value, &sigma_value_error, &sigmaR_value, &sigmaL_value, &ChiSquare, MC, j, EndCaps, temp, DataChain, mean, rms, c2, Erawframe, f, nomDossierFits, nomFichier, MinRange, MaxRange, &JanChi2, &DegreesOfFreedom, &pValue, r9sup, MinVar[j], MaxVar[j], variableX, isMC);
                        }
			if(nomFitMethode == "RooGamma")   
                        {
				if(r9sup == 2) nomDossierFits += "Fits/mmg_sGammar9All";
                                if(r9sup == 1) nomDossierFits += "Fits/mmg_sGammar9sup";
				if(r9sup == 0) nomDossierFits += "Fits/mmg_sGammar9inf";
                                RooGamma2(&mean_value, &mean_error, &sigmaEff_value, &sigma_value, &sigma_value_error, &sigmaR_value, &sigmaL_value, &ChiSquare, MC, j, EndCaps, temp, DataChain, mean, rms, c2, Erawframe, f, nomDossierFits, nomFichier, MinRange, MaxRange, &JanChi2, &DegreesOfFreedom, &pValue, r9sup, MinVar[j], MaxVar[j], variableX, isMC);
                        }
			if(nomFitMethode == "RooBifurcatedGauss")   
                        {
				if(r9sup == 2) nomDossierFits += "Fits/mmg_sBFGr9All";
                                if(r9sup == 1) nomDossierFits += "Fits/mmg_sBFGr9sup";
				if(r9sup == 0) nomDossierFits += "Fits/mmg_sBFGr9inf";
                                RooBifurcatedGauss(&mean_value, &mean_error, &sigmaEff_value, &sigma_value_error, &sigma_value, &sigmaR_value, &sigmaL_value, &ChiSquare, &minLogLikelihood, &differenceLogLikelihood, MC, j, EndCaps, temp, DataChain, mean, rms, c2, Erawframe, f, nomDossierFits, nomFichier, MinRange, MaxRange, &JanChi2, &DegreesOfFreedom, &pValue, r9sup, MinVar[j], MaxVar[j], variableX, isMC);
                        }
			if(nomFitMethode == "RooSumGauss")   
                        {
				if(r9sup == 2) nomDossierFits += "Fits/mmg_sSumGr9All";
                                if(r9sup == 1) nomDossierFits += "Fits/mmg_sSumGr9sup";
				if(r9sup == 0) nomDossierFits += "Fits/mmg_sSumGr9inf";
                                RooSumGauss(&mean_value, &mean_error, &sigmaEff_value, &sigma_value, &sigma_value_error, &sigmaR_value, &sigmaL_value, &ChiSquare, MC, j, EndCaps, temp, DataChain, mean, rms, c2, Erawframe, f, nomDossierFits, nomFichier, MinRange, MaxRange, &JanChi2, &DegreesOfFreedom, &pValue, r9sup, MinVar[j], MaxVar[j], variableX, isMC);
                        }
			if(nomFitMethode == "RooGenericPDF")   
                        {
				if(r9sup == 2) nomDossierFits += "Fits/mmg_sCBGr9All";
                                if(r9sup == 1) nomDossierFits += "Fits/mmg_sCBGr9sup";
				if(r9sup == 0) nomDossierFits += "Fits/mmg_sCBGr9inf";
                                RooGenericPDF(&mean_value, &mean_error, &sigmaEff_value, &sigma_value, &sigma_value_error, &sigmaR_value, &sigmaL_value, &ChiSquare, MC, j, EndCaps, temp, DataChain, mean, rms, c2, Erawframe, f, nomDossierFits, nomFichier);
                        }
			if(nomFitMethode == "RooLandau2")   
                        {
				if(r9sup == 2) nomDossierFits += "Fits/mmg_sLandaur9All";
                                if(r9sup == 1) nomDossierFits += "Fits/mmg_sLandaur9sup";
                                if(r9sup == 0) nomDossierFits += "Fits/mmg_sLandaur9inf";
                                RooLandau2(&mean_value, &mean_error, &sigmaEff_value, &sigma_value_error, &sigma_value, &sigmaR_value, &sigmaL_value, &ChiSquare, MC, j, EndCaps, temp, DataChain, mean, rms, c2, Erawframe, f, nomDossierFits, nomFichier, MinRange, MaxRange, &JanChi2, &DegreesOfFreedom, &pValue, r9sup, MinVar[j], MaxVar[j], variableX, isMC);
			}
			if(nomFitMethode == "RooLandauConvGaussian")   
                        {
				if(r9sup == 2) nomDossierFits += "Fits/mmg_sLandauConvGausr9All";
                                if(r9sup == 1) nomDossierFits += "Fits/mmg_sLandauConvGausr9sup";
                                if(r9sup == 0) nomDossierFits += "Fits/mmg_sLandauConvGausr9inf";
                                RooLandauConvGaussian(&mean_value, &mean_error, &sigmaEff_value, &sigma_value_error, &sigma_value, &sigmaR_value, &sigmaL_value, &ChiSquare, MC, j, EndCaps, temp, DataChain, mean, rms, c2, Erawframe, f, nomDossierFits, nomFichier, MinRange, MaxRange, &JanChi2, &DegreesOfFreedom, &pValue, r9sup, MinVar[j], MaxVar[j], variableX, isMC);
			}			
			if(nomFitMethode == "RooKernel")
			{
				if(r9sup == 2) nomDossierFits += "Fits/mmg_sKernelr9All";
                                if(r9sup == 1) nomDossierFits += "Fits/mmg_sKernelr9sup";
                                if(r9sup == 0) nomDossierFits += "Fits/mmg_sKernelr9inf";
                                RooKernel(&mean_value, &mean_error, &sigmaEff_value, &sigma_value_error, &sigma_value, &sigmaR_value, &sigmaL_value, &ChiSquare, &minLogLikelihood, &differenceLogLikelihood, MC, j, EndCaps, temp, DataChain, mean, rms, c2, Erawframe, f, nomDossierFits, nomFichier, MinRange, MaxRange, &JanChi2, &DegreesOfFreedom, &pValue, r9sup, MinVar[j], MaxVar[j], variableX, isMC);
			}
			if(nomFitMethode == "RooVoigtian2")
			{
				if(r9sup == 2) nomDossierFits += "Fits/mmg_sVoigtianr9All";
                                if(r9sup == 1) nomDossierFits += "Fits/mmg_sVoigtianr9sup";
                                if(r9sup == 0) nomDossierFits += "Fits/mmg_sVoigtianr9inf";
                                RooVoigtian2(&mean_value, &mean_error, &sigmaEff_value, &sigma_value_error, &sigma_value, &sigmaR_value, &sigmaL_value, &ChiSquare, &minLogLikelihood, &differenceLogLikelihood, MC, j, EndCaps, temp, DataChain, mean, rms, c2, Erawframe, f, nomDossierFits, nomFichier, MinRange, MaxRange, &JanChi2, &DegreesOfFreedom, &pValue, r9sup, MinVar[j], MaxVar[j], variableX, isMC);
			}
			if(nomFitMethode == "Cruijff")
                        {
                                if(r9sup == 2) nomDossierFits += "Fits/mmg_sCruijffr9All";
                                if(r9sup == 1) nomDossierFits += "Fits/mmg_sCruijffr9sup";
                                if(r9sup == 0) nomDossierFits += "Fits/mmg_sCruijffr9inf";
                                Cruijff(&mean_value, &mean_error, &sigmaEff_value, &sigma_value_error, &sigma_value, &sigmaR_value, &sigmaL_value, &ChiSquare, &minLogLikelihood, &differenceLogLikelihood, MC, j, EndCaps, temp, DataChain, mean, rms, c2, Erawframe, f, nomDossierFits, nomFichier, MinRange, MaxRange, &JanChi2, &DegreesOfFreedom, &pValue, r9sup, MinVar[j], MaxVar[j], variableX, isMC);
                        }	
			if(nomFitMethode == "roofit_bwcb_fft2")
                        {
                                if(r9sup == 2) nomDossierFits += "Fits/mmg_sBWxCBr9All";
                                if(r9sup == 1) nomDossierFits += "Fits/mmg_sBWxCBr9sup";
                                if(r9sup == 0) nomDossierFits += "Fits/mmg_sBWxCBr9inf";
                                roofit_bwcb_fft2(&mean_value, &mean_error, &sigmaEff_value, &sigma_value_error, &sigma_value, &sigmaR_value, &sigmaL_value, &ChiSquare, &minLogLikelihood, &differenceLogLikelihood, MC, j, EndCaps, temp, DataChain, mean, rms, c2, Erawframe, f, nomDossierFits, nomFichier, MinRange, MaxRange, &JanChi2, &DegreesOfFreedom, &pValue, r9sup, MinVar[j], MaxVar[j], variableX, isMC);
			}

			MeanTab[j] = mean_value; //Mean
                        MeanErrorTab[j] = mean_error; //MeanError 
                        ChiSquareTab[j] = ChiSquare; //Chi2
			SigmaEffTab[j] = sigmaEff_value;
			SigmaTab[j] = sigma_value;
			SigmaErrorTab[j] = sigma_value_error;
                        SigmaRTab[j] = sigmaR_value;
                        SigmaLTab[j] = sigmaL_value;
			JanChiSquareTab[j] = JanChi2;
                        PValueTab[j] = pValue;
                        DegreesOfFreedomTab[j] = DegreesOfFreedom;
			minLogLikelihoodTab[j] = minLogLikelihood;
			differenceLogLikelihoodTab[j] = differenceLogLikelihood;		
	
/*

			c2->Clear();
			TH1D *htest = f->GetHistogram(); //ICI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        		//MC->Draw("");
        		//htest->Draw("SAMES");
			//htest->Draw("");             
			//htest->SetLineColor(4);
			f->Draw("");
			f->SetLineColor(4);	
			cout<<endl<<"f->GetXmax() = "<<f->GetXmax()<<endl;			
	

			Double_t res[200]; 
        		double pvalue = MC->Chi2Test(htest,"WW P",res);
        		cout<<endl<<"P-value = "<<pvalue<<endl;
			c2->Print("testPvalue.png");
*/
			f->Delete();
			f = 0;
			Erawframe->Delete();
			Erawframe = 0;

		}


		else
		{

			CrystalBallMethode(&mean_value, &mean_error, &sigma_value, &ChiSquare, &DegreesOfFreedom, params, *MC);//faire des if et des form()
		

			//TCanvas *c1 = new TCanvas("c1", "c1",0,0,600,600);
			c2->ToggleEventStatus();

			TPaveStats *ptstats = new TPaveStats(0.6845638,0.7226277,0.9865772,0.9912587,"brNDC");
        		ptstats->SetName("stats");
        		ptstats->SetBorderSize(2);
        		ptstats->SetFillColor(kWhite);
        		ptstats->SetTextColor(4);
        		ptstats->SetTextAlign(12);
        		ptstats->SetTextFont(42);
        		ptstats->SetOptStat(1110);
        		ptstats->SetOptFit(11111111);
        		ptstats->Draw();
        		MC->GetListOfFunctions()->Add(ptstats);
        		ptstats->SetParent(MC->GetListOfFunctions());
        		MC->SetLineColor(1);
        		MC->SetMarkerStyle(6);
        		MC->GetXaxis()->SetTitle("E_{SCraw*Ceta*PtCor}/E_{TRUE}");
       			MC->GetXaxis()->SetLabelFont(42);
        		MC->GetXaxis()->SetTitleFont(42);
        		MC->GetYaxis()->SetTitle("Number of events");
        		MC->GetYaxis()->SetLabelFont(42);
        		MC->GetYaxis()->SetTitleOffset(1.27);
        		MC->GetYaxis()->SetTitleFont(42);

        		MC->Draw("");
        		TF1* func = new TF1("func",CrystalBall,0.65, 1.1,5);
        		func->SetParameters(params[0], params[1], params[2], params[3], params[4]);
        		func->SetLineColor(4);
        		func->SetLineWidth(3);
        		MC->Fit(func);
			
			TPaveText *pt;
        		if(MC->GetMaximum(100000) < 1000) pt = new TPaveText(0.1057047,0.9124088,0.6694631,0.9912587,"blNDC");
        		if(MC->GetMaximum(100000) >= 1000) pt = new TPaveText(0.1812081,0.9124088,0.6694631,0.9912587,"blNDC");
        		pt->SetName("title");
        		pt->SetBorderSize(2);
        		pt->SetFillColor(kWhite);
        		pt->SetTextFont(42);
			//TLatex *text = new TLatex();	
			TText * text;
			if(EndCaps == 0) text = pt->AddText("MC #gamma#gamma, Barrel");
                	if(EndCaps == 1) text = pt->AddText("MC #gamma#gamma, End Caps");

        		pt->Draw();
	
			c2->SetTickx(1);
        		c2->SetTicky(1);
			c2->SetGridx(1);
        		c2->SetGridy(1);
        		c2->Modified();
        		c2->cd();
        		c2->SetSelected(c2);
        		c2->ToggleToolBar();


			MeanTab[j] = mean_value; //Mean
                	MeanErrorTab[j] = mean_error; //MeanError 
                	if(DegreesOfFreedom != 0) ChiSquareTab[j] = ChiSquare / DegreesOfFreedom; //Chi2
                	else ChiSquareTab[j] = 0;
                	SigmaRTab[j] = SigmaR(func, 0, 3);
                	SigmaLTab[j] = SigmaL(func, 0, 3);


			if(EndCaps == 0) c2->Print(Form("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/Fits/mmg_sEB/mmg_s%d.png",j));
                	if(EndCaps == 0) c2->Print(Form("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/Fits/mmg_sEB/mmg_s%d.ps",j));
                	if(EndCaps == 0) c2->Print(Form("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/Fits/mmg_sEB/mmg_s%d.pdf",j));
                	if(EndCaps == 0) c2->Print(Form("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/Fits/mmg_sEB/mmg_s%d.gif",j));
                	if(EndCaps == 0) c2->Print(Form("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/Fits/mmg_sEB/mmg_s%d.C",j));

               		if(EndCaps == 1) c2->Print(Form("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/Fits/mmg_sEE/mmg_s%d.png",j));
			if(EndCaps == 1) c2->Print(Form("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/Fits/mmg_sEE/mmg_s%d.ps",j));
			if(EndCaps == 1) c2->Print(Form("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/Fits/mmg_sEE/mmg_s%d.pdf",j));
			if(EndCaps == 1) c2->Print(Form("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/Fits/mmg_sEE/mmg_s%d.gif",j));
               		if(EndCaps == 1) c2->Print(Form("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/Fits/mmg_sEE/mmg_s%d.C",j));
	


			func->Delete();
			c2->Clear();
               	 	//c1->Close();
			//delete c1;
		}

		//DataChain->Delete();
                MC->Delete();
		MC = 0;
		//VarTH1->Delete();
       		VarXvalue->Delete();
		VarXvalue = 0;
		//c2->Delete();
		delete c2;
		c2 = 0;
	}
	//delete VarXvalue;
        //delete MC;
	//delete f;
	//delete Erawframe;
	DataChain->Delete();
	DataChain = 0;
	//f->Delete();
	//Erawframe->Delete();

	//file1.Close();
	//monFlux.close();


	cout<<endl<<"////////// ----- Generation TGraphs ----- //////////"<<endl<<endl;


	//gROOT->SetStyle("Plain");


	TCanvas* c1 = new TCanvas("c1", "c1",0,0,600,600);
	

//Generation des bins

	for(int k = 0; k < n ; k++) 
	{
		xErrorR[k] = MaxVar[k] - xValue[k];
		xErrorL[k] = xValue[k] - MinVar[k];
	
		////////// Modifs mean et sigma plots de Jan a enlever si non comparaison  //////////
		MeanTab[k] = MeanTab[k] * 100.0;
		MeanErrorTab[k] = MeanErrorTab[k] * 100.0; 
		SigmaTab[k] = SigmaTab[k] * 100.0;
		SigmaErrorTab[k] = SigmaErrorTab[k] * 100.0;

	}

/*
//Generation du Chi2/ndf moyen

	cout<<endl<<endl;
        double chi2mean = 0; 
	char * buffer = new char[ 25 ];
 
        for(int v = 0; v < n; v++) 
        {
                chi2mean += ChiSquareTab[v];
        }

        cout<<endl<<"chi2mean = "<<chi2mean/5;
        cout<<endl;     
	
	chi2mean = chi2mean / 5;

*/


//Enregistrement parametres fichiers
	
	if(EndCaps == 0) nomDossierFits += "EB/";
	if(EndCaps == 1) nomDossierFits += "EE/";
	
	FILE *f1 = fopen(Form("%sMeanTab.txt",nomDossierFits.c_str()), "w" );
	delete f1;
	FILE *f2 = fopen(Form("%sMeanErrorTab.txt",nomDossierFits.c_str()), "w" );
        delete f2;
	FILE *f3 = fopen(Form("%sSigmaRTab.txt",nomDossierFits.c_str()), "w" );
        delete f3;
	FILE *f4 = fopen(Form("%sSigmaLTab.txt",nomDossierFits.c_str()), "w" );
        delete f4;
	FILE *f5 = fopen(Form("%sSigmaTab.txt",nomDossierFits.c_str()), "w" );
        delete f5;
	FILE *f6 = fopen(Form("%sMinVar.txt",nomDossierFits.c_str()), "w" );
        delete f6;
	FILE *f7 = fopen(Form("%sMaxVar.txt",nomDossierFits.c_str()), "w" );
        delete f7;
	FILE *f8 = fopen(Form("%sxValue.txt",nomDossierFits.c_str()), "w" );
        delete f8;
	FILE *f9 = fopen(Form("%sxErrorR.txt",nomDossierFits.c_str()), "w" );
        delete f9;
	FILE *f10 = fopen(Form("%sxErrorL.txt",nomDossierFits.c_str()), "w" );
        delete f10;
	FILE *f11 = fopen(Form("%syValue.txt",nomDossierFits.c_str()), "w" );
        delete f11;
	FILE *f12 = fopen(Form("%sxErrorRight.txt",nomDossierFits.c_str()), "w" );
        delete f12;
	FILE *f13 = fopen(Form("%sxErrorLeft.txt",nomDossierFits.c_str()), "w" );
        delete f13;
	FILE *f14 = fopen(Form("%sChiSquareTab.txt",nomDossierFits.c_str()), "w" );
        delete f14;
	FILE *f15 = fopen(Form("%sSigmaEffTab.txt",nomDossierFits.c_str()), "w");
	delete f15;
	FILE *f16 = fopen(Form("%sSigmaErrorTab.txt",nomDossierFits.c_str()), "w" );
        delete f16;
	FILE *f17 = fopen(Form("%sJanChiSquareTab.txt",nomDossierFits.c_str()), "w" );
        delete f17;
        FILE *f18 = fopen(Form("%sPValueTab.txt",nomDossierFits.c_str()), "w" );
        delete f18;
        FILE *f19 = fopen(Form("%sDegreesOfFreedomTab.txt",nomDossierFits.c_str()), "w" );
        delete f19;
	FILE *f20 = fopen(Form("%sminLogLikelihoodTab.txt",nomDossierFits.c_str()), "w" );
        delete f20;
	FILE *f21 = fopen(Form("%sdifferenceLogLikelihoodTab.txt",nomDossierFits.c_str()), "w" );
        delete f21;
	

	for(int m = 0; m <n ; m++)
	{

                ofstream monFlux2(Form("%sMeanTab.txt",nomDossierFits.c_str()), ios::app);

                monFlux2 << MeanTab[m] <<endl;

                monFlux2.close();

		ofstream monFlux3(Form("%sMeanErrorTab.txt",nomDossierFits.c_str()), ios::app);

                monFlux3 << MeanErrorTab[m] <<endl;

                monFlux3.close();

		ofstream monFlux4(Form("%sSigmaRTab.txt",nomDossierFits.c_str()), ios::app);

                monFlux4 << SigmaRTab[m] <<endl;

                monFlux4.close();

		ofstream monFlux5(Form("%sSigmaLTab.txt",nomDossierFits.c_str()), ios::app);

                monFlux5 << SigmaLTab[m] <<endl;

                monFlux5.close();

		ofstream monFlux6(Form("%sSigmaTab.txt",nomDossierFits.c_str()), ios::app);

                monFlux6 << SigmaTab[m] <<endl;

                monFlux6.close();

		ofstream monFlux7(Form("%sMinVar.txt",nomDossierFits.c_str()), ios::app);

                monFlux7 << MinVar[m] <<endl;

                monFlux7.close();

                ofstream monFlux8(Form("%sMaxVar.txt",nomDossierFits.c_str()), ios::app);

                monFlux8 << MaxVar[m] <<endl;

                monFlux8.close();

                ofstream monFlux9(Form("%sxValue.txt",nomDossierFits.c_str()), ios::app);

                monFlux9 << xValue[m] <<endl;

                monFlux9.close();

                ofstream monFlux10(Form("%sxErrorR.txt",nomDossierFits.c_str()), ios::app);

                monFlux10 << xErrorR[m] <<endl;

                monFlux10.close();

                ofstream monFlux11(Form("%sxErrorL.txt",nomDossierFits.c_str()), ios::app);

                monFlux11 << xErrorL[m] <<endl;

                monFlux11.close();

		ofstream monFlux12(Form("%syValue.txt",nomDossierFits.c_str()), ios::app);

                monFlux12 << yValue[m] <<endl;

                monFlux12.close();

                ofstream monFlux13(Form("%sxErrorRight.txt",nomDossierFits.c_str()), ios::app);

                monFlux13 << xErrorRight[m] <<endl;

                monFlux13.close();
		
		ofstream monFlux14(Form("%sxErrorLeft.txt",nomDossierFits.c_str()), ios::app);

                monFlux14 << xErrorLeft[m] <<endl;

                monFlux14.close();

		ofstream monFlux15(Form("%sChiSquareTab.txt",nomDossierFits.c_str()), ios::app);

                monFlux15 << ChiSquareTab[m] <<endl;

                monFlux15.close();

		ofstream monFlux16(Form("%sSigmaEffTab.txt",nomDossierFits.c_str()), ios::app);

                monFlux16 << SigmaEffTab[m] <<endl;

                monFlux16.close();

		ofstream monFlux17(Form("%sSigmaErrorTab.txt",nomDossierFits.c_str()), ios::app);

                monFlux17 << SigmaErrorTab[m] <<endl;

                monFlux17.close();

		ofstream monFlux18(Form("%sJanChiSquareTab.txt",nomDossierFits.c_str()), ios::app);

                monFlux18 << JanChiSquareTab[m] <<endl;

                monFlux18.close();

                ofstream monFlux19(Form("%sPValueTab.txt",nomDossierFits.c_str()), ios::app);

                monFlux19 << PValueTab[m] <<endl;

                monFlux19.close();

                ofstream monFlux20(Form("%sDegreesOfFreedomTab.txt",nomDossierFits.c_str()), ios::app);

                monFlux20 << DegreesOfFreedomTab[m] <<endl;

                monFlux20.close();

		ofstream monFlux21(Form("%sminLogLikelihoodTab.txt",nomDossierFits.c_str()), ios::app);

                monFlux21 << minLogLikelihoodTab[m] <<endl;

                monFlux21.close();
	
		ofstream monFlux22(Form("%sdifferenceLogLikelihoodTab.txt",nomDossierFits.c_str()), ios::app);

                monFlux22 << differenceLogLikelihoodTab[m] <<endl;

                monFlux22.close();


	}
	


//Graphe E_{RAW*Ceta*PtCor}/E_{TRUE} Converted (MC) vs Var MC 2gammas

	nomDossierTGraph = nomDossier;
	if(r9sup == 2)
        {
                if(nomFitMethode == "RooCrystalBall") nomDossierTGraph += "TGraph/mmg_sCBr9All";
                if(nomFitMethode == "RooLogNormal") nomDossierTGraph += "TGraph/mmg_sLNr9All";
                if(nomFitMethode == "RooGauss") nomDossierTGraph += "TGraph/mmg_sGr9All";
                if(nomFitMethode == "RooGamma") nomDossierTGraph += "TGraph/mmg_sGammar9All";
                if(nomFitMethode == "RooBifurcatedGauss") nomDossierTGraph += "TGraph/mmg_sBFGr9All";
                if(nomFitMethode == "RooSumGauss") nomDossierTGraph += "TGraph/mmg_sSumGr9All";
                if(nomFitMethode == "RooGenericPDF") nomDossierTGraph += "TGraph/mmg_sCBGr9All";
                if(nomFitMethode == "RooLandau2") nomDossierTGraph += "TGraph/mmg_sRooLandaur9All";
                if(nomFitMethode == "RooLandauConvGaussian") nomDossierTGraph += "TGraph/mmg_sLandauConvGausr9All";
        	if(nomFitMethode == "RooKernel") nomDossierTGraph += "TGraph/mmg_sKernelr9All";
		if(nomFitMethode == "RooVoigtian2") nomDossierTGraph += "TGraph/mmg_sVoigtianr9All";
		if(nomFitMethode == "Cruijff") nomDossierTGraph += "TGraph/mmg_sCruijffr9All";
		if(nomFitMethode == "roofit_bwcb_fft2") nomDossierTGraph += "TGraph/mmg_sBWxCBr9All";
	}

	if(r9sup == 1)
	{
		if(nomFitMethode == "RooCrystalBall") nomDossierTGraph += "TGraph/mmg_sCBr9sup";
		if(nomFitMethode == "RooLogNormal") nomDossierTGraph += "TGraph/mmg_sLNr9sup";
		if(nomFitMethode == "RooGauss") nomDossierTGraph += "TGraph/mmg_sGr9sup";
		if(nomFitMethode == "RooGamma") nomDossierTGraph += "TGraph/mmg_sGammar9sup";
		if(nomFitMethode == "RooBifurcatedGauss") nomDossierTGraph += "TGraph/mmg_sBFGr9sup";
		if(nomFitMethode == "RooSumGauss") nomDossierTGraph += "TGraph/mmg_sSumGr9sup";
		if(nomFitMethode == "RooGenericPDF") nomDossierTGraph += "TGraph/mmg_sCBGr9sup";
		if(nomFitMethode == "RooLandau2") nomDossierTGraph += "TGraph/mmg_sRooLandaur9sup";
		if(nomFitMethode == "RooLandauConvGaussian") nomDossierTGraph += "TGraph/mmg_sLandauConvGausr9sup";	
		if(nomFitMethode == "RooKernel") nomDossierTGraph += "TGraph/mmg_sKernelr9sup";
		if(nomFitMethode == "RooVoigtian2") nomDossierTGraph += "TGraph/mmg_sVoigtianr9sup";
		if(nomFitMethode == "Cruijff") nomDossierTGraph += "TGraph/mmg_sCruijffr9sup";
		if(nomFitMethode == "roofit_bwcb_fft2") nomDossierTGraph += "TGraph/mmg_sBWxCBr9sup";
	}

	if(r9sup == 0)
        {
                if(nomFitMethode == "RooCrystalBall") nomDossierTGraph += "TGraph/mmg_sCBr9inf";
                if(nomFitMethode == "RooLogNormal") nomDossierTGraph += "TGraph/mmg_sLNr9inf";
                if(nomFitMethode == "RooGauss") nomDossierTGraph += "TGraph/mmg_sGr9inf";
                if(nomFitMethode == "RooGamma") nomDossierTGraph += "TGraph/mmg_sGammar9inf";
                if(nomFitMethode == "RooBifurcatedGauss") nomDossierTGraph += "TGraph/mmg_sBFGr9inf";
                if(nomFitMethode == "RooSumGauss") nomDossierTGraph += "TGraph/mmg_sSumGr9inf";
                if(nomFitMethode == "RooGenericPDF") nomDossierTGraph += "TGraph/mmg_sCBGr9inf";  
		if(nomFitMethode == "RooLandau2") nomDossierTGraph += "TGraph/mmg_sRooLandaur9inf";
                if(nomFitMethode == "RooLandauConvGaussian") nomDossierTGraph += "TGraph/mmg_sLandauConvGausr9inf";   
        	if(nomFitMethode == "RooKernel") nomDossierTGraph += "TGraph/mmg_sKernelr9inf";
		if(nomFitMethode == "RooVoigtian2") nomDossierTGraph += "TGraph/mmg_sVoigtianr9inf";
		if(nomFitMethode == "Cruijff") nomDossierTGraph += "TGraph/mmg_sCruijffr9inf";
		if(nomFitMethode == "roofit_bwcb_fft2") nomDossierTGraph += "TGraph/mmg_sBWxCBr9inf";
	}	

	//TGraphErrors* MCtg = new TGraphErrors(n,xValue,MeanTab,xError,MeanErrorTab);
	TGraphAsymmErrors * MCtg = new TGraphAsymmErrors(n,xValue, MeanTab, xErrorL, xErrorR, MeanErrorTab, MeanErrorTab);

        c1->ToggleEventStatus();
	MCtg->SetTitle("");
        MCtg->SetLineColor(4);
	if(variableX == "Photon_SC_rawEt") MCtg->GetXaxis()->SetTitle("E_{T RAW}");
        if(variableX == "Photon_Et") MCtg->GetXaxis()->SetTitle("P_{T}");
        if(variableX == "Photon_E") MCtg->GetXaxis()->SetTitle("E");
        if(variableX == "Photon_SC_Eta") MCtg->GetXaxis()->SetTitle("#eta");
        if(variableX == "Photon_SC_brem") MCtg->GetXaxis()->SetTitle("#sigma_{#phi}/#sigma_{#eta}");
	//MCtg->GetXaxis()->SetTitle("Var");
   	MCtg->GetXaxis()->SetLabelFont(42);
   	MCtg->GetXaxis()->SetTitleFont(42);
	MCtg->GetXaxis()->SetLabelSize(0.03);
	MCtg->GetYaxis()->SetTitle("s (%)");
	MCtg->GetYaxis()->SetLabelFont(42);
   	MCtg->GetYaxis()->SetTitleOffset(1.24);
   	MCtg->GetYaxis()->SetTitleFont(42);
	MCtg->GetYaxis()->SetLabelSize(0.03);
	//MCtg->SetMarkerColor(4);
        //MCtg->SetMarkerStyle(21);
        //MCtg->SetMarkerSize(0.6);
	MCtg->Draw("AP");  

 
        TPaveText *pt = new TPaveText(0.1073826,0.9125874,0.897651,0.9912587,"blNDC");
   	pt->SetName("title");
   	pt->SetBorderSize(2);
   	pt->SetFillColor(kWhite);
   	pt->SetTextFont(42);
	TText * text2;
	text2 = pt->AddText("");
        pt->Draw();

/*
	ChaineChi2(buffer, chi2meanData);	

	pt = new TPaveText(16.13571,1.197034,20.84785,1.244703,"br");
	pt->SetFillColor(19);
   	pt->SetTextColor(4);
   	text = pt->AddText("Data s raw");
   	TLine *line = pt->AddLine(0,0.5238107,0,0.5238107);
   	text = pt->AddText(buffer);
	pt->Draw();
   
	ChaineChi2(buffer, chi2mean);

  	pt = new TPaveText(16.13571,1.149365,20.84785,1.197034,"br");
   	pt->SetFillColor(19);
   	pt->SetTextColor(2);
   	text = pt->AddText("MC s raw");
   	line = pt->AddLine(0,0.5238107,0,0.5238107);
   	text = pt->AddText(buffer);
	pt->Draw();

	ChaineChi2(buffer, chi2mean2);
   
   	pt = new TPaveText(16.13571,1.101696,20.84785,1.149365,"br");
   	pt->SetFillColor(19);
   	pt->SetTextColor(3);
   	text = pt->AddText("MC E_{RAW*Ceta*VarCor}/E_{TRUE} Converted (MC)");
   	line = pt->AddLine(0,0.5238107,0,0.5238107);
   	text = pt->AddText(buffer);
	pt->Draw();
*/
	TLatex *textF = new TLatex();
        textF->SetNDC();
        textF->SetTextAlign(11);
        //textF->SetTextFont(42);
        textF->SetTextSizePixels(17);
	textF->SetTextSize(0.038);

        //textF->DrawLatex(0.135, 0.93, "s vs Var, MC Z#rightarrow#mu#mu#gamma");
	
	gStyle->SetPadBorderMode(0);



	TLatex *textL = new TLatex();

/*        textL->SetNDC();
        textL->SetTextAlign(11);
        textL->SetTextFont(42);
        textL->SetTextSizePixels(17);
	textL->SetTextSize(0.030);
        textL->DrawLatex(0.16, 0.85, "#color[4]{Data}");
        gStyle->SetPadBorderMode(0);

	textL = new TLatex();
	textL->SetNDC();
        textL->SetTextAlign(11);
        textL->SetTextFont(42);
        textL->SetTextSizePixels(17);
	textL->SetTextSize(0.030);
        textL->DrawLatex(0.30, 0.85, "#color[2]{MC (s)}");
        gStyle->SetPadBorderMode(0);

	textL = new TLatex();
        textL->SetNDC();
        textL->SetTextAlign(11);
        textL->SetTextFont(42);
        textL->SetTextSizePixels(17);
	textL->SetTextSize(0.030);
        textL->DrawLatex(0.48, 0.85, "#color[3]{MC (E_{RAW*Ceta*VarCor}/E_{TRUE} Converted (MC))}");
        gStyle->SetPadBorderMode(0);
*/
	//textL = new TLatex();
        textL->SetNDC();
        textL->SetTextAlign(11);
        textL->SetTextFont(42);
        textL->SetTextSizePixels(17);
        textL->SetTextSize(0.028);
        //if(EndCaps == 0) textL->DrawLatex(0.75, 0.83, "Barrel");
	//if(EndCaps == 1) textL->DrawLatex(0.72, 0.83, "End Caps");
	//if(EndCaps == 2) textL->DrawLatex(0.75, 0.83, "");
        textL->DrawLatex(0.66, 0.88, "CMS Preliminary 2011");
        if(EndCaps == 0 && r9sup == 0) textL->DrawLatex(0.66, 0.83, "Barrel, Low r9");
        if(EndCaps == 0 && r9sup == 1) textL->DrawLatex(0.66, 0.83, "Barrel, High r9");
        if(EndCaps == 1 && r9sup == 0) textL->DrawLatex(0.66, 0.83, "Endcaps, Low r9");
        if(EndCaps == 1 && r9sup == 1) textL->DrawLatex(0.66, 0.83, "Endcaps, High r9");
	if(EndCaps == 0 && r9sup == 2) textL->DrawLatex(0.66, 0.83, "Barrel, All r9");
	if(EndCaps == 1 && r9sup == 2) textL->DrawLatex(0.66, 0.83, "Endcaps, All r9");
	gStyle->SetPadBorderMode(0);

/*

	TGraphErrors* Datad = new TGraphErrors(n,xValueData,MeanTabData,xErrorData,MeanErrorTabData);
	Datad->SetLineColor(4);	
	Datad->Draw("SAMES");

	TGraphErrors* Datae = new TGraphErrors(n,xValue,MeanTab2,xError,MeanErrorTab2);
        Datae->SetLineColor(3);
        Datae->Draw("SAMES");

*/	
	if(EndCaps == 0) MCtg->GetYaxis()->SetRangeUser(yminmmg_s,ymaxmmg_s);//a changer
	if(EndCaps == 1) MCtg->GetYaxis()->SetRangeUser(yminmmg_s,ymaxmmg_s);//a changer
	if(EndCaps == 2) MCtg->GetYaxis()->SetRangeUser(yminmmg_s,ymaxmmg_s);
	
	MCtg->GetXaxis()->SetLimits(xminmmg_s,xmaxmmg_s);

	c1->SetTickx(1);
   	c1->SetTicky(1);
   	c1->SetGridx(1);
	c1->SetGridy(1);
	c1->Modified();
   	c1->cd();
   	c1->SetSelected(c1);
   	c1->ToggleToolBar();

	nomFichier = "mmg_sVsVar";
	enregistrementPlots(nomDossierTGraph, nomFichier, EndCaps, 10000, c1);

	MCtg->Delete();
	MCtg = 0;

	c1->Clear();


////////// PValue VS PT TGraph //////////


	TGraphAsymmErrors * PValueTg = new TGraphAsymmErrors(n,xValue, PValueTab, xErrorL, xErrorR, PValueErrorTab, PValueErrorTab);

        c1->ToggleEventStatus();
	PValueTg->SetFillColor(1);
        PValueTg->SetLineColor(4);
        PValueTg->SetMarkerColor(4);
        PValueTg->SetMarkerStyle(21);
        PValueTg->SetMarkerSize(0.6);
        PValueTg->SetLineColor(4);
	PValueTg->SetTitle("");
        PValueTg->SetLineColor(4);
	if(variableX == "Photon_SC_rawEt") PValueTg->GetXaxis()->SetTitle("E_{T RAW}");
        if(variableX == "Photon_Et") PValueTg->GetXaxis()->SetTitle("P_{T}");
        if(variableX == "Photon_E") PValueTg->GetXaxis()->SetTitle("E");
        if(variableX == "Photon_SC_Eta") PValueTg->GetXaxis()->SetTitle("#eta");
        if(variableX == "Photon_SC_brem") PValueTg->GetXaxis()->SetTitle("#sigma_{#phi}/#sigma_{#eta}");
	//PValueTg->GetXaxis()->SetTitle("Var");
   	PValueTg->GetXaxis()->SetLabelFont(42);
   	PValueTg->GetXaxis()->SetTitleFont(42);
	PValueTg->GetXaxis()->SetLabelSize(0.03);
	PValueTg->GetYaxis()->SetTitle("p-value");
	PValueTg->GetYaxis()->SetLabelFont(42);
   	PValueTg->GetYaxis()->SetTitleOffset(1.24);
   	PValueTg->GetYaxis()->SetTitleFont(42);
	PValueTg->GetYaxis()->SetLabelSize(0.03);
	//PValueTg->SetMarkerColor(4);
        //PValueTg->SetMarkerStyle(21);
        //PValueTg->SetMarkerSize(0.6);
	PValueTg->Draw("AP");  

	//pt->Delete();
	//delete pt;
	//pt = 0;
 
        pt = new TPaveText(0.1073826,0.9125874,0.897651,0.9912587,"blNDC");
   	pt->SetName("title");
   	pt->SetBorderSize(2);
   	pt->SetFillColor(kWhite);
   	pt->SetTextFont(42);
	//TText * text2;
	text2 = pt->AddText("");
        pt->Draw();

	//TLatex *textF = new TLatex();
        textF->SetNDC();
        textF->SetTextAlign(11);
        //textF->SetTextFont(42);
        textF->SetTextSizePixels(17);
	textF->SetTextSize(0.038);

	
	gStyle->SetPadBorderMode(0);


	//textL->Delete();
	//textL = 0;

	textL = new TLatex();
        textL->SetNDC();
        textL->SetTextAlign(11);
        textL->SetTextFont(42);
        textL->SetTextSizePixels(17);
        textL->SetTextSize(0.028);
        //if(EndCaps == 0) textL->DrawLatex(0.75, 0.83, "Barrel");
	//if(EndCaps == 1) textL->DrawLatex(0.72, 0.83, "End Caps");
	//if(EndCaps == 2) textL->DrawLatex(0.75, 0.83, "");
        textL->DrawLatex(0.66, 0.88, "CMS Preliminary 2011");
        if(EndCaps == 0 && r9sup == 0) textL->DrawLatex(0.66, 0.83, "Barrel, Low r9");
        if(EndCaps == 0 && r9sup == 1) textL->DrawLatex(0.66, 0.83, "Barrel, High r9");
        if(EndCaps == 1 && r9sup == 0) textL->DrawLatex(0.66, 0.83, "Endcaps, Low r9");
        if(EndCaps == 1 && r9sup == 1) textL->DrawLatex(0.66, 0.83, "Endcaps, High r9");
	if(EndCaps == 0 && r9sup == 2) textL->DrawLatex(0.66, 0.83, "Barrel, All r9");
        if(EndCaps == 1 && r9sup == 2) textL->DrawLatex(0.66, 0.83, "Endcaps, All r9");
	gStyle->SetPadBorderMode(0);

	PValueTg->GetYaxis()->SetRangeUser(yminPValue,ymaxPValue);//a changer
	
	PValueTg->GetXaxis()->SetLimits(xminPValue,xmaxPValue);

	c1->SetTickx(1);
   	c1->SetTicky(1);
   	c1->SetGridx(1);
	c1->SetGridy(1);
	c1->Modified();
   	c1->cd();
   	c1->SetSelected(c1);
   	c1->ToggleToolBar();

	nomFichier = "PValueVsVar";
	enregistrementPlots(nomDossierTGraph, nomFichier, EndCaps, 10000, c1);

	PValueTg->Delete();
	PValueTg = 0;

	c1->Clear();
/*	
        if(EndCaps == 0) 
	{	
		c1->Print("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/TGraph/mmg_sEB/mmg_sVsVar.png");
        	c1->Print("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/TGraph/mmg_sEB/mmg_sVsVar.ps");
        	c1->Print("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/TGraph/mmg_sEB/mmg_sVsVar.pdf");
        	c1->Print("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/TGraph/mmg_sEB/mmg_sVsVar.gif");
		c1->Print("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/TGraph/mmg_sEB/mmg_sVsVar.C");
	}

	if(EndCaps == 1) 
        {       
                c1->Print("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/TGraph/mmg_sEE/mmg_sVsVar.png");
                c1->Print("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/TGraph/mmg_sEE/mmg_sVsVar.ps");
                c1->Print("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/TGraph/mmg_sEE/mmg_sVsVar.pdf");
                c1->Print("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/TGraph/mmg_sEE/mmg_sVsVar.gif");
                c1->Print("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/TGraph/mmg_sEE/mmg_sVsVar.C");
        }
*/	
/*		

	//test fonction Anne-Fleur
	double xValuefVar[100] = {0};
        double fVarTab1[100] = {0};
	double fVarTab2[100] = {0};
	double fVarTab3[100] = {0};
	double fVarTab4[100] = {0};
        double fVarTab5[100] = {0};
	double  xErrorFVar[100] = {0};
        double  fVarError[100] = {0};
        for(int bin = 0; bin<100; bin ++)
        {
                xValuefVar[bin] = (21.0 * bin) / 100.0;

                if(EndCaps == 0) fVarTab1[bin] = VarCor(xValuefVar[bin], 0, 1);
                if(EndCaps == 1) fVarTab1[bin] = VarCor(xValuefVar[bin], 1, 1);
		if(EndCaps == 0) fVarTab2[bin] = VarCor(xValuefVar[bin], 0, 2);
                if(EndCaps == 1) fVarTab2[bin] = VarCor(xValuefVar[bin], 1, 2);
		if(EndCaps == 0) fVarTab3[bin] = VarCor(xValuefVar[bin], 0, 3);
                if(EndCaps == 1) fVarTab3[bin] = VarCor(xValuefVar[bin], 1, 3);
		if(EndCaps == 0) fVarTab4[bin] = VarCor(xValuefVar[bin], 0, 4);
                if(EndCaps == 1) fVarTab4[bin] = VarCor(xValuefVar[bin], 1, 4);
		if(EndCaps == 0) fVarTab5[bin] = VarCor(xValuefVar[bin], 0, 5);
                if(EndCaps == 1) fVarTab5[bin] = VarCor(xValuefVar[bin], 1, 5);	
		
        }

	fVarTab1[0] = fVarTab1[1];
	fVarTab2[0] = fVarTab2[1];
	fVarTab3[0] = fVarTab3[1];
	fVarTab4[0] = fVarTab4[1];

        TGraphErrors* VarCorGraph1 = new TGraphErrors(100,xValuefVar,fVarTab1,xErrorFVar,fVarError);
        VarCorGraph1->SetLineColor(417);
        VarCorGraph1->SetLineWidth(2);
        VarCorGraph1->Draw("CSAME");

	TGraphErrors* VarCorGraph2 = new TGraphErrors(100,xValuefVar,fVarTab2,xErrorFVar,fVarError);
        VarCorGraph2->SetLineColor(433);
        VarCorGraph2->SetLineWidth(2);
        VarCorGraph2->Draw("CSAME");

	TGraphErrors* VarCorGraph3 = new TGraphErrors(100,xValuefVar,fVarTab3,xErrorFVar,fVarError);
        VarCorGraph3->SetLineColor(617);
        VarCorGraph3->SetLineWidth(2);
        VarCorGraph3->Draw("CSAME");

	TGraphErrors* VarCorGraph4 = new TGraphErrors(100,xValuefVar,fVarTab4,xErrorFVar,fVarError);
        VarCorGraph4->SetLineColor(633);
        VarCorGraph4->SetLineWidth(2);
        VarCorGraph4->Draw("CSAME");

	TGraphErrors* VarCorGraph5 = new TGraphErrors(100,xValuefVar,fVarTab5,xErrorFVar,fVarError);
        VarCorGraph5->SetLineColor(418);
        VarCorGraph5->SetLineWidth(2);
        VarCorGraph5->Draw("CSAME");


	TLatex *tex = new TLatex();

	tex = new TLatex(0.20,0.83,"#color[417]{pre CMSSW_3_10_0 Electron VarCor}");
	tex->SetNDC();
	tex->SetTextFont(42);
	tex->SetTextSize(0.030);
	tex->SetLineWidth(2);
	tex->Draw();

	tex = new TLatex(0.20,0.79,"#color[433]{Electron VarCor}");
        tex->SetNDC();
        tex->SetTextFont(42);
        tex->SetTextSize(0.030);
        tex->SetLineWidth(2);
        tex->Draw();

	tex = new TLatex(0.20,0.75,"#color[617]{CMSSW_4_1_3 Photon VarCor}");
        tex->SetNDC();
        tex->SetTextFont(42);
        tex->SetTextSize(0.030);
        tex->SetLineWidth(2);
        tex->Draw();

	tex = new TLatex(0.20,0.71,"#color[633]{CMSSW_4_1_3 New Photon VarCor}");
        tex->SetNDC();
        tex->SetTextFont(42);
        tex->SetTextSize(0.030);
        tex->SetLineWidth(2);
        tex->Draw();

	tex = new TLatex(0.20,0.67,"#color[418]{New Photon VarCor}");
        tex->SetNDC();
        tex->SetTextFont(42);
        tex->SetTextSize(0.030);
        tex->SetLineWidth(2);
        tex->Draw();


        gStyle->SetPadBorderMode(0);

	c1->Modified();
        c1->cd();

	if(EndCaps == 0)
        {
                c1->Print("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/TGraph/mmg_sEB/mmg_sVsVarVarCor.png");
                c1->Print("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/TGraph/mmg_sEB/mmg_sVsVarVarCor.ps");
                c1->Print("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/TGraph/mmg_sEB/mmg_sVsVarVarCor.pdf");
                c1->Print("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/TGraph/mmg_sEB/mmg_sVsVarVarCor.gif");
                c1->Print("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/TGraph/mmg_sEB/mmg_sVsVarVarCor.C");
        }

        if(EndCaps == 1)
        {
                c1->Print("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/TGraph/mmg_sEE/mmg_sVsVarVarCor.png");
                c1->Print("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/TGraph/mmg_sEE/mmg_sVsVarVarCor.ps");
                c1->Print("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/TGraph/mmg_sEE/mmg_sVsVarVarCor.pdf");
                c1->Print("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/TGraph/mmg_sEE/mmg_sVsVarVarCor.gif");
                c1->Print("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/TGraph/mmg_sEE/mmg_sVsVarVarCor.C");
        }	

	c1->Clear();
	MCtg->Delete();
*/

//*****************TH1D*****************//

/*
	gStyle->SetOptTitle(0);
        gROOT->ProcessLine(".x setTDRStyle.C");

	TH1D * MCth = new TH1D("MCth","test", 10000, xminmmg_s, xmaxmmg_s);
        int numBin = 0;
        //double i2 = 0;
        int numBinTrue = 0;
*/
/*
	ifstream monFlux20("xValue.txt");
	double xValueTabTemp[1000] = NULL;

	for(int i3 = 0; i3 < 172; i3++)
        {
		cout<<endl<<"i3 = "<<i3;
		monFlux20 >> xValueTabTemp[i3];
                cout<<endl<<"xValueTabTemp = "<<xValueTabTemp[i3];
        }

	monFlux20.close();
*/

/*
        for(double i2 = 0.0; i2 < 100.0 ; i2 +=0.0100)
        {

		if(fabs(xValue[numBinTrue]-i2)<=0.0006)
		{
		
			//cout<<endl<<"xValue ="<<xValue[numBinTrue];
                	//cout<<endl<<"i2 = "<<i2;
                	//cout<<endl<<"difference = "<<xValue[numBinTrue]-i2;
			//cout<<endl<<"difference = "<<setprecision(6)<<abs(xValue[numBinTrue]-i2);
			//cout<<endl<<"numBinTrue = "<<numBinTrue;
                	//cout<<endl<<"i2 = "<<i2;
		        MCth->SetBinContent(numBin,MeanTab[numBinTrue]);
                        MCth->SetBinError(numBin,MeanErrorTab[numBinTrue]);
                        numBinTrue++;
                }

                numBin++;

        }


        MCth->SetMarkerColor(2);
        MCth->SetMarkerStyle(7);
        //Data->SetMarkerSize(0.5);

	c1->ToggleEventStatus();
        MCth->SetTitle("");
        MCth->SetLineColor(4);
        MCth->GetXaxis()->SetTitle("Var");
        MCth->GetXaxis()->SetLabelFont(42);
        MCth->GetXaxis()->SetTitleFont(42);
        MCth->GetXaxis()->SetLabelSize(0.03);
        MCth->GetYaxis()->SetTitle("s");
        MCth->GetYaxis()->SetLabelFont(42);
        MCth->GetYaxis()->SetTitleOffset(1.24);
        MCth->GetYaxis()->SetTitleFont(42);
        MCth->GetYaxis()->SetLabelSize(0.03);
        MCth->Draw("");


        pt = new TPaveText(0.1073826,0.9125874,0.897651,0.9912587,"blNDC");
        pt->SetName("title");
        pt->SetBorderSize(2);
        pt->SetFillColor(kWhite);
        pt->SetTextFont(42);
        TText * text3;
	text3 = pt->AddText("");
        pt->Draw();

	c1->Print("TH1D.png");

	myfunc(EndCaps);
        myfit(MCth);

	TLatex *textH = new TLatex();
        textH->SetNDC();
        textH->SetTextAlign(11);
        //textF->SetTextFont(42);
        textH->SetTextSizePixels(17);
        textH->SetTextSize(0.038);

        //textH->DrawLatex(0.100, 0.93, "s, MC Z#rightarrow#mu#mu#gamma");

        gStyle->SetPadBorderMode(0);
        TLatex *textI = new TLatex();

	textI = new TLatex();
        textI->SetNDC();
        textI->SetTextAlign(11);
        textI->SetTextFont(42);
        textI->SetTextSizePixels(17);
        textI->SetTextSize(0.038);
        if(EndCaps == 0) textI->DrawLatex(0.14, 0.83, "Barrel");
        if(EndCaps == 1) textI->DrawLatex(0.14, 0.83, "End Caps");
        if(EndCaps == 2) textI->DrawLatex(0.14, 0.83, "");
        gStyle->SetPadBorderMode(0);

	if(EndCaps == 0) MCth->GetYaxis()->SetRangeUser(yminmmg_s,ymaxmmg_s);//a changer
        if(EndCaps == 1) MCth->GetYaxis()->SetRangeUser(yminmmg_s,ymaxmmg_s);//a changer
        if(EndCaps == 2) MCth->GetYaxis()->SetRangeUser(yminmmg_s,ymaxmmg_s);

        MCth->GetXaxis()->SetLimits(xminmmg_s,xmaxmmg_s);

        c1->SetTickx(1);
        c1->SetTicky(1);
        c1->SetGridx(1);
        c1->SetGridy(1);
	c1->Modified();
        c1->cd();
        c1->SetSelected(c1);
        c1->ToggleToolBar();

	nomFichier = "mmg_sVsVarTH1";
        enregistrementPlots(nomDossierTGraph, nomFichier, EndCaps, 10000, c1);
*/
/*
	if(EndCaps == 0)
        {
                c1->Print("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/TGraph/mmg_sEB/mmg_sVsVarTH1.png");
                c1->Print("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/TGraph/mmg_sEB/mmg_sVsVarTH1.ps");
                c1->Print("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/TGraph/mmg_sEB/mmg_sVsVarTH1.pdf");
                c1->Print("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/TGraph/mmg_sEB/mmg_sVsVarTH1.gif");
                c1->Print("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/TGraph/mmg_sEB/mmg_sVsVarTH1.C");
        }

        if(EndCaps == 1)
        {
                c1->Print("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/TGraph/mmg_sEE/mmg_sVsVarTH1.png");
                c1->Print("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/TGraph/mmg_sEE/mmg_sVsVarTH1.ps");
                c1->Print("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/TGraph/mmg_sEE/mmg_sVsVarTH1.pdf");
                c1->Print("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/TGraph/mmg_sEE/mmg_sVsVarTH1.gif");
                c1->Print("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/TGraph/mmg_sEE/mmg_sVsVarTH1.C");
        }
*/
//	c1->Clear();
        //c1->Delete();
//	MCth->Delete();


//*****************Graph chi2 et sigmaRL*****************//
	
	gROOT->SetStyle("Plain");
       	gStyle->SetOptTitle(0);

	
	//c1 = new TCanvas("c1", "c1",0,0,600,600);
	c1->Range(0,0,1,1);
   	c1->SetBorderSize(2);
   	c1->SetFrameFillColor(0);	
	
	TPad *c1_1 = new TPad("c1_1", "c1_1",0.01,0.51,0.99,0.99);
   	c1_1->Draw();
   	c1_1->cd();
   	c1_1->Range(-12.5,-0.4455269,12.5,1.229671);
   	c1_1->SetBorderSize(2);
   	c1_1->SetFrameFillColor(0);
	//c1_1->SetLogy();
        c1_1->SetTickx(1);
        c1_1->SetTicky(1);	
	c1_1->SetGridx(1);
        c1_1->SetGridy(1);

        TGraph* chi2 = new TGraph(n,xValue,ChiSquareTab);
	chi2->SetFillColor(1);
   	chi2->SetLineColor(4);
   	chi2->SetMarkerColor(4);
   	chi2->SetMarkerStyle(21);
   	chi2->SetMarkerSize(0.6);

        chi2->GetYaxis()->SetTitle("#chi^{2} / ndf");
	chi2->GetYaxis()->CenterTitle(true);
	chi2->GetYaxis()->SetLabelSize(0.06);
	chi2->GetXaxis()->SetLabelSize(0.06);
	chi2->GetYaxis()->SetTitleSize(0.07);
	chi2->GetYaxis()->SetTitleOffset(0.5);
	chi2->GetXaxis()->SetLabelFont(42);
	chi2->GetYaxis()->SetTitleFont(42);
	chi2->GetYaxis()->SetLabelFont(42);
	chi2->Draw("AP");
	chi2->SetTitle("");
/*
        TGraph* chi2MC = new TGraph(n,xValue,ChiSquareTab);
	chi2MC->SetFillColor(1);
        chi2MC->SetLineColor(2);
        chi2MC->SetMarkerColor(2);
        chi2MC->SetMarkerStyle(21);
        chi2MC->SetMarkerSize(0.6);
	chi2MC->Draw("P");

        TGraph* chi2MC2 = new TGraph(n,xValue,ChiSquareTab2);   
	chi2MC2->SetFillColor(1);
        chi2MC2->SetLineColor(3);
        chi2MC2->SetMarkerColor(3);
        chi2MC2->SetMarkerStyle(21);
        chi2MC2->SetMarkerSize(0.6);

	chi2MC2->Draw("P");
*/	
	chi2->GetYaxis()->SetRangeUser(yminChi2,ymaxChi2);
	chi2->GetXaxis()->SetLimits(xminChi2,xmaxChi2);
	//chi2->GetXaxis()->SetRangeUser(0,17);

	//pt->Delete();
	//pt = 0;

	pt = new TPaveText(0,0,0,0,"blNDC");
        pt->SetName("title");
        pt->SetBorderSize(2);
        pt->SetFillColor(kWhite);
        TText * text4 = pt->AddText("");
        pt->Draw();


	TLine *line = new TLine(0,1,xmaxChi2,1);
        line->SetLineStyle(3);
        line->Draw();

	c1_1->Modified();
   	c1->cd();


	TPad *c1_2 = new TPad("c1_2", "c1_2",0.01,0.01,0.99,0.49);
  	c1_2->Draw();
   	c1_2->cd();
   	c1_2->Range(-12.5,-1.375519,12.5,1.380519);
   	c1_2->SetBorderSize(2);
   	c1_2->SetFrameFillColor(0);
	c1_2->SetTickx(1);
        c1_2->SetTicky(1);
	c1_2->SetGridx(1);
        c1_2->SetGridy(1);


	TGraphAsymmErrors * sigmaRLGraph = new TGraphAsymmErrors(n,xValue, yValue, xErrorLeft, xErrorRight, SigmaLTab, SigmaRTab);
	sigmaRLGraph->SetFillColor(4);
        sigmaRLGraph->SetLineColor(4);
        sigmaRLGraph->GetYaxis()->SetTitle("_{L} - #sigma_{R}");
        sigmaRLGraph->GetYaxis()->CenterTitle(true);
        sigmaRLGraph->GetYaxis()->SetLabelSize(0.06);
        sigmaRLGraph->GetXaxis()->SetLabelSize(0.06);
        sigmaRLGraph->GetYaxis()->SetTitleSize(0.07);
        sigmaRLGraph->GetYaxis()->SetTitleOffset(0.5);
        sigmaRLGraph->GetXaxis()->SetLabelFont(42);
        sigmaRLGraph->GetYaxis()->SetTitleFont(42);
        sigmaRLGraph->GetYaxis()->SetLabelFont(42);
	sigmaRLGraph->SetLineWidth(1);
        sigmaRLGraph->Draw("AP");	

	sigmaRLGraph->SetTitle("");

/*
	TGraphAsymmErrors * sigmaRLGraph = new TGraphAsymmErrors(n,xValue, yValue, xErrorLeft, xErrorRight, SigmaLTab, SigmaRTab);
	sigmaRLGraph->SetLineColor(2);
        sigmaRLGraph->SetLineWidth(1);
        sigmaRLGraph->Draw("SAMES");

	TGraphAsymmErrors * sigmaRLGraph2 = new TGraphAsymmErrors(n,xValue, yValue, xErrorLeft, xErrorRight, SigmaLTab2, SigmaRTab2);
	sigmaRLGraph2->SetLineColor(3);
        sigmaRLGraph2->SetLineWidth(1);
        sigmaRLGraph2->Draw("SAMES");
*/
	//pt->Delete();
	//pt = 0;
	
	pt = new TPaveText(0,0,0,0,"blNDC");
        pt->SetName("title");
        pt->SetBorderSize(2);
        pt->SetFillColor(kWhite);
        TText *text5 = pt->AddText("");
        pt->Draw();
		
	sigmaRLGraph->GetYaxis()->SetRangeUser(yminSigmaLSigmaR,ymaxSigmaLSigmaR);


	sigmaRLGraph->GetXaxis()->SetLimits(xminSigmaLSigmaR,xmaxSigmaLSigmaR);
        TLine *line2 = new TLine(0,0,xmaxSigmaLSigmaR,0);
        line2->SetLineStyle(3);
        line2->Draw();
	
	
	c1_2->Modified();
   	c1->cd();
   	c1->Modified();
   	c1->cd();
   	c1->SetSelected(c1);

	nomFichier = "Chi2SigmaLRVsVar";
        enregistrementPlots(nomDossierTGraph, nomFichier, EndCaps, 10000, c1);

/*
	if(EndCaps == 0)
        {
                c1->Print("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/TGraph/mmg_sEB/Chi2SigmaLRmmg_sVsVar.png");
                c1->Print("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/TGraph/mmg_sEB/Chi2SigmaLRmmg_sVsVar.ps");
                c1->Print("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/TGraph/mmg_sEB/Chi2SigmaLRmmg_sVsVar.pdf");
                c1->Print("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/TGraph/mmg_sEB/Chi2SigmaLRmmg_sVsVar.gif");
        	c1->Print("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/TGraph/mmg_sEB/Chi2SigmaLRmmg_sVsVar.C");
	}

	if(EndCaps == 1)
        {
                c1->Print("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/TGraph/mmg_sEE/Chi2SigmaLRmmg_sVsVar.png");
                c1->Print("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/TGraph/mmg_sEE/Chi2SigmaLRmmg_sVsVar.ps");
                c1->Print("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/TGraph/mmg_sEE/Chi2SigmaLRmmg_sVsVar.pdf");
                c1->Print("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/TGraph/mmg_sEE/Chi2SigmaLRmmg_sVsVar.gif");
                c1->Print("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/TGraph/mmg_sEE/Chi2SigmaLRmmg_sVsVar.C");
        }
*/
	c1->Clear();
        //c1->Delete();
	chi2->Delete();
	chi2 = 0;
	sigmaRLGraph->Delete();
	sigmaRLGraph = 0;

//*****************Graph sigma et sigmaEff*****************//


       	gStyle->SetOptTitle(0);

	
	//c1 = new TCanvas("c1", "c1",0,0,600,600);
	c1->Range(0,0,1,1);
   	c1->SetBorderSize(2);
   	c1->SetFrameFillColor(0);	
	
	//c1_1->Delete();
	//c1_1 = 0;

	c1_1 = new TPad("c1_1", "c1_1",0.01,0.51,0.99,0.99);
   	c1_1->Draw();
   	c1_1->cd();
   	c1_1->Range(-12.5,-0.4455269,12.5,1.229671);
   	c1_1->SetBorderSize(2);
   	c1_1->SetFrameFillColor(0);
	//c1_1->SetLogy();
        c1_1->SetTickx(1);
        c1_1->SetTicky(1);
	c1_1->SetGridx(1);
        c1_1->SetGridy(1);	
	
	TGraphAsymmErrors * sigmaGraph = new TGraphAsymmErrors(n,xValue, yValue, xErrorLeft, xErrorRight, SigmaTab, SigmaTab);
        sigmaGraph->SetFillColor(4);
        sigmaGraph->SetLineColor(4);
        sigmaGraph->GetYaxis()->SetTitle(" (%)");
        sigmaGraph->GetYaxis()->CenterTitle(true);
        sigmaGraph->GetYaxis()->SetLabelSize(0.06);
        sigmaGraph->GetXaxis()->SetLabelSize(0.06);
        sigmaGraph->GetYaxis()->SetTitleSize(0.07);
        sigmaGraph->GetYaxis()->SetTitleOffset(0.5);
        sigmaGraph->GetXaxis()->SetLabelFont(42);
        sigmaGraph->GetYaxis()->SetTitleFont(42);
        sigmaGraph->GetYaxis()->SetLabelFont(42);
        sigmaGraph->SetLineWidth(1);
        sigmaGraph->Draw("AP");

        sigmaGraph->SetTitle("");

	//pt->Delete();
	//pt = 0;

	pt = new TPaveText(0,0,0,0,"blNDC");
        pt->SetName("title");
        pt->SetBorderSize(2);
        pt->SetFillColor(kWhite);
        TText *text6 = pt->AddText("");
        pt->Draw();

        sigmaGraph->GetYaxis()->SetRangeUser(yminSigma,ymaxSigma);


        sigmaGraph->GetXaxis()->SetLimits(xminSigma,xmaxSigma);
        TLine *line3 = new TLine(0,0,xmaxSigma,0);
        line3->SetLineStyle(3);
        line3->Draw();


	c1_1->Modified();
   	c1->cd();

	//c1_2->Delete();
	//c1_2 = 0;

	c1_2 = new TPad("c1_2", "c1_2",0.01,0.01,0.99,0.49);
  	c1_2->Draw();
   	c1_2->cd();
   	c1_2->Range(-12.5,-1.375519,12.5,1.380519);
   	c1_2->SetBorderSize(2);
   	c1_2->SetFrameFillColor(0);
	c1_2->SetTickx(1);
        c1_2->SetTicky(1);
	c1_2->SetGridx(1);
        c1_2->SetGridy(1);


	TGraphAsymmErrors * sigmaEffGraph = new TGraphAsymmErrors(n,xValue, yValue, xErrorLeft, xErrorRight, SigmaEffTab, SigmaEffTab);
	sigmaEffGraph->SetFillColor(4);
        sigmaEffGraph->SetLineColor(4);
        sigmaEffGraph->GetYaxis()->SetTitle("_{eff}");
        sigmaEffGraph->GetYaxis()->CenterTitle(true);
        sigmaEffGraph->GetYaxis()->SetLabelSize(0.06);
        sigmaEffGraph->GetXaxis()->SetLabelSize(0.06);
        sigmaEffGraph->GetYaxis()->SetTitleSize(0.07);
        sigmaEffGraph->GetYaxis()->SetTitleOffset(0.5);
        sigmaEffGraph->GetXaxis()->SetLabelFont(42);
        sigmaEffGraph->GetYaxis()->SetTitleFont(42);
        sigmaEffGraph->GetYaxis()->SetLabelFont(42);
	sigmaEffGraph->SetLineWidth(1);
        sigmaEffGraph->Draw("AP");	

	sigmaEffGraph->SetTitle("");

	//pt->Delete();
        //pt = 0;

	pt = new TPaveText(0,0,0,0,"blNDC");
        pt->SetName("title");
        pt->SetBorderSize(2);
        pt->SetFillColor(kWhite);
        TText *text7 = pt->AddText("");
        pt->Draw();
		
	sigmaEffGraph->GetYaxis()->SetRangeUser(yminSigmaEff,ymaxSigmaEff);


	sigmaEffGraph->GetXaxis()->SetLimits(xminSigmaEff,xmaxSigmaEff);
        TLine *line4 = new TLine(0,0,xmaxSigmaEff,0);
        line4->SetLineStyle(3);
        line4->Draw();
	
	
	c1_2->Modified();
   	c1->cd();
   	c1->Modified();
   	c1->cd();
   	c1->SetSelected(c1);
	
	nomFichier = "SigmaSigmaEffVsVar";
        enregistrementPlots(nomDossierTGraph, nomFichier, EndCaps, 10000, c1);

/*
	if(EndCaps == 0)
        {
                c1->Print("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/TGraph/mmg_sEB/SigmaSigmaEffrawCVarOverEtrueVsVar.png");
                c1->Print("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/TGraph/mmg_sEB/SigmaSigmaEffrawCVarOverEtrueVsVar.ps");
                c1->Print("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/TGraph/mmg_sEB/SigmaSigmaEffrawCVarOverEtrueVsVar.pdf");
                c1->Print("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/TGraph/mmg_sEB/SigmaSigmaEffrawCVarOverEtrueVsVar.gif");
        	c1->Print("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/TGraph/mmg_sEB/SigmaSigmaEffrawCVarOverEtrueVsVar.C");
	}

	if(EndCaps == 1)
        {
                c1->Print("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/TGraph/mmg_sEE/SigmaSigmaEffrawCVarOverEtrueVsVar.png");
                c1->Print("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/TGraph/mmg_sEE/SigmaSigmaEffrawCVarOverEtrueVsVar.ps");
                c1->Print("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/TGraph/mmg_sEE/SigmaSigmaEffrawCVarOverEtrueVsVar.pdf");
                c1->Print("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/TGraph/mmg_sEE/SigmaSigmaEffrawCVarOverEtrueVsVar.gif");
                c1->Print("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/TGraph/mmg_sEE/SigmaSigmaEffrawCVarOverEtrueVsVar.C");
        }
*/
	c1->Clear();
        //c1->Delete();
	sigmaGraph->Delete();
	sigmaGraph = 0;
	sigmaEffGraph->Delete();
	sigmaEffGraph = 0;

//*******TH1D sigmaEff******//


	//c1 = new TCanvas("c1", "c1",0,0,600,600);

	TGraph* sigmaEffTg = new TGraph(n,xValue,SigmaEffTab);

	//pt->Delete();
	//pt = 0;

        pt = new TPaveText(0.1073826,0.9125874,0.897651,0.9912587,"blNDC");
        pt->SetName("title");
        pt->SetBorderSize(2);
        pt->SetFillColor(kWhite);
        pt->SetTextFont(42);
        TText * text = pt->AddText("");
        pt->Draw();



	sigmaEffTg->SetFillColor(1);
        sigmaEffTg->SetLineColor(4);
        sigmaEffTg->SetMarkerColor(4);
        sigmaEffTg->SetMarkerStyle(21);
        sigmaEffTg->SetMarkerSize(0.6);

        sigmaEffTg->GetYaxis()->SetTitle("_{eff}");
        sigmaEffTg->GetYaxis()->SetLabelSize(0.03);
        sigmaEffTg->GetXaxis()->SetLabelSize(0.03);
        //sigmaEffTg->GetYaxis()->SetTitleSize(0.07);
        sigmaEffTg->GetYaxis()->SetTitleOffset(1.24);
	if(variableX == "Photon_SC_rawEt") sigmaEffTg->GetXaxis()->SetTitle("E_{T RAW}");
        if(variableX == "Photon_Et") sigmaEffTg->GetXaxis()->SetTitle("P_{T}");
        if(variableX == "Photon_E") sigmaEffTg->GetXaxis()->SetTitle("E");
        if(variableX == "Photon_SC_Eta") sigmaEffTg->GetXaxis()->SetTitle("#eta");
        if(variableX == "Photon_SC_brem") sigmaEffTg->GetXaxis()->SetTitle("#sigma_{#phi}/#sigma_{#eta}");
	//sigmaEffTg->GetXaxis()->SetTitle("Var");
        sigmaEffTg->GetXaxis()->SetLabelFont(42);
        sigmaEffTg->GetYaxis()->SetTitleFont(42);
        sigmaEffTg->GetYaxis()->SetLabelFont(42);
        sigmaEffTg->Draw("AP");
        sigmaEffTg->SetTitle("");


	sigmaEffTg->GetYaxis()->SetRangeUser(yminSigmaEffTg,ymaxSigmaEffTg);
        sigmaEffTg->GetXaxis()->SetLimits(xminSigmaEffTg,xmaxSigmaEffTg);
        //sigmaEffTg->GetXaxis()->SetRangeUser(0,17);

	//pt->Delete();
        //pt = 0;

        pt = new TPaveText(0,0,0,0,"blNDC");
        pt->SetName("title");
        pt->SetBorderSize(2);
        pt->SetFillColor(kWhite);
        TText *text8 = pt->AddText("");
        pt->Draw();
	//if(EndCaps == 0) textL->DrawLatex(0.75, 0.83, "Barrel");
        //if(EndCaps == 1) textL->DrawLatex(0.72, 0.83, "End Caps");
	textL->DrawLatex(0.66, 0.88, "CMS Preliminary 2011");
        if(EndCaps == 0 && r9sup == 0) textL->DrawLatex(0.66, 0.83, "Barrel, Low r9");
        if(EndCaps == 0 && r9sup == 1) textL->DrawLatex(0.66, 0.83, "Barrel, High r9");
        if(EndCaps == 1 && r9sup == 0) textL->DrawLatex(0.66, 0.83, "Endcaps, Low r9");
        if(EndCaps == 1 && r9sup == 1) textL->DrawLatex(0.66, 0.83, "Endcaps, High r9");
	if(EndCaps == 0 && r9sup == 2) textL->DrawLatex(0.66, 0.83, "Barrel, All r9");
        if(EndCaps == 1 && r9sup == 2) textL->DrawLatex(0.66, 0.83, "Endcaps, All r9");

	c1->SetTickx(1);
        c1->SetTicky(1);
	c1->SetGridx(1);
        c1->SetGridy(1);
        c1->Modified();
        c1->cd();
        c1->SetSelected(c1);
        c1->ToggleToolBar();

	nomFichier = "SigmaEffVsVarTg";
        enregistrementPlots(nomDossierTGraph, nomFichier, EndCaps, 10000, c1);
	
	c1->Clear();
	
	sigmaEffTg->Delete();
	sigmaEffTg = 0;

/*

	if(EndCaps == 0)
        {
                c1->Print("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/TGraph/mmg_sEB/SigmaEffSCrawCetaOverEtrueVsVarTg.png");
                c1->Print("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/TGraph/mmg_sEB/SigmaEffSCrawCetaOverEtrueVsVarTg.ps");
                c1->Print("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/TGraph/mmg_sEB/SigmaEffSCrawCetaOverEtrueVsVarTg.pdf");
                c1->Print("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/TGraph/mmg_sEB/SigmaEffSCrawCetaOverEtrueVsVarTg.gif");
                c1->Print("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/TGraph/mmg_sEB/SigmaEffSCrawCetaOverEtrueVsVarTg.C");
        }

        if(EndCaps == 1)
        {
                c1->Print("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/TGraph/mmg_sEE/SigmaEffSCrawCetaOverEtrueVsVarTg.png");
                c1->Print("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/TGraph/mmg_sEE/SigmaEffSCrawCetaOverEtrueVsVarTg.ps");
                c1->Print("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/TGraph/mmg_sEE/SigmaEffSCrawCetaOverEtrueVsVarTg.pdf");
                c1->Print("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/TGraph/mmg_sEE/SigmaEffSCrawCetaOverEtrueVsVarTg.gif");
                c1->Print("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsmmg_sJan/Var/TGraph/mmg_sEE/SigmaEffSCrawCetaOverEtrueVsVarTg.C");
        }

*/

////////// TGraphAsym sigma CB /////////

	TGraphAsymmErrors * SigmaCBtg = new TGraphAsymmErrors(n,xValue, SigmaTab, xErrorL, xErrorR, SigmaErrorTab, SigmaErrorTab);

        SigmaCBtg->SetFillColor(1);
        SigmaCBtg->SetLineColor(4);
        //SigmaCBtg->SetMarkerColor(4);
        //SigmaCBtg->SetMarkerStyle(21);
        //SigmaCBtg->SetMarkerSize(0.6);

        SigmaCBtg->GetYaxis()->SetTitle("_{CB} (%)");
        SigmaCBtg->GetYaxis()->SetLabelSize(0.03);
        SigmaCBtg->GetXaxis()->SetLabelSize(0.03);
        //SigmaCBtg->GetYaxis()->SetTitleSize(0.07);
        SigmaCBtg->GetYaxis()->SetTitleOffset(1.24);
        if(variableX == "Photon_SC_rawEt") SigmaCBtg->GetXaxis()->SetTitle("E_{T RAW}");
        if(variableX == "Photon_Et") SigmaCBtg->GetXaxis()->SetTitle("P_{T}");
        if(variableX == "Photon_E") SigmaCBtg->GetXaxis()->SetTitle("E");
        if(variableX == "Photon_SC_Eta") SigmaCBtg->GetXaxis()->SetTitle("#eta");
        if(variableX == "Photon_SC_brem") SigmaCBtg->GetXaxis()->SetTitle("#sigma_{#phi}/#sigma_{#eta}");
	//SigmaCBtg->GetXaxis()->SetTitle("Var");
        SigmaCBtg->GetXaxis()->SetLabelFont(42);
        SigmaCBtg->GetYaxis()->SetTitleFont(42);
        SigmaCBtg->GetYaxis()->SetLabelFont(42);
        SigmaCBtg->Draw("AP");


        SigmaCBtg->GetYaxis()->SetRangeUser(yminSigmaTg,ymaxSigmaTg);
        SigmaCBtg->GetXaxis()->SetLimits(xminSigmaTg,xmaxSigmaTg);
        //SigmaCBtg->GetXaxis()->SetRangeUser(0,17);
	//if(EndCaps == 0) textL->DrawLatex(0.75, 0.83, "Barrel");
        //if(EndCaps == 1) textL->DrawLatex(0.72, 0.83, "End Caps");
	textL->DrawLatex(0.66, 0.88, "CMS Preliminary 2011");
        if(EndCaps == 0 && r9sup == 0) textL->DrawLatex(0.66, 0.83, "Barrel, Low r9");
        if(EndCaps == 0 && r9sup == 1) textL->DrawLatex(0.66, 0.83, "Barrel, High r9");
        if(EndCaps == 1 && r9sup == 0) textL->DrawLatex(0.66, 0.83, "Endcaps, Low r9");
        if(EndCaps == 1 && r9sup == 1) textL->DrawLatex(0.66, 0.83, "Endcaps, High r9");
	if(EndCaps == 0 && r9sup == 2) textL->DrawLatex(0.66, 0.83, "Barrel, All r9");
        if(EndCaps == 1 && r9sup == 2) textL->DrawLatex(0.66, 0.83, "Endcaps, All r9");

	c1->SetTickx(1);
        c1->SetTicky(1);
	c1->SetGridx(1);
        c1->SetGridy(1);
        c1->Modified();
        c1->cd();
        c1->SetSelected(c1);
        c1->ToggleToolBar();

        nomFichier = "SigmaCBVsVarTg";
        enregistrementPlots(nomDossierTGraph, nomFichier, EndCaps, 10000, c1);

	c1->Clear();

	SigmaCBtg->Delete();
	SigmaCBtg = 0;

////////// Graph Janchi2 et p-value ///////////

	
	//c1 = new TCanvas("c1", "c1",0,0,600,600);
	c1->Range(0,0,1,1);
   	c1->SetBorderSize(2);
   	c1->SetFrameFillColor(0);	

	//c1_1->Delete();
	//c1_1 = 0;
	
	c1_1 = new TPad("c1_1", "c1_1",0.01,0.51,0.99,0.99);
   	c1_1->Draw();
   	c1_1->cd();
   	c1_1->Range(-12.5,-0.4455269,12.5,1.229671);
   	c1_1->SetBorderSize(2);
   	c1_1->SetFrameFillColor(0);
	//c1_1->SetLogy();
        c1_1->SetTickx(1);
        c1_1->SetTicky(1);	
	c1_1->SetGridx(1);
        c1_1->SetGridy(1);


        TGraph* Janchi2 = new TGraph(n,xValue,JanChiSquareTab);
	Janchi2->SetFillColor(1);
   	Janchi2->SetLineColor(4);
   	Janchi2->SetMarkerColor(4);
   	Janchi2->SetMarkerStyle(21);
   	Janchi2->SetMarkerSize(0.6);

        Janchi2->GetYaxis()->SetTitle("Jan #chi^{2} / ndf");
	Janchi2->GetYaxis()->CenterTitle(true);
	Janchi2->GetYaxis()->SetLabelSize(0.06);
	Janchi2->GetXaxis()->SetLabelSize(0.06);
	Janchi2->GetYaxis()->SetTitleSize(0.07);
	Janchi2->GetYaxis()->SetTitleOffset(0.5);
	Janchi2->GetXaxis()->SetLabelFont(42);
	Janchi2->GetYaxis()->SetTitleFont(42);
	Janchi2->GetYaxis()->SetLabelFont(42);
	Janchi2->Draw("AP");
	Janchi2->SetTitle("");
	
	Janchi2->GetYaxis()->SetRangeUser(yminJanChi2,ymaxJanChi2);
	Janchi2->GetXaxis()->SetLimits(xminJanChi2,xmaxJanChi2);
	//Janchi2->GetXaxis()->SetRangeUser(0,17);

	//pt->Delete();
	//pt = 0;

	pt = new TPaveText(0,0,0,0,"blNDC");
        pt->SetName("title");
        pt->SetBorderSize(2);
        pt->SetFillColor(kWhite);
        //TText * text11 = pt->AddText("");
        pt->Draw();

	//line->Delete();
	//line = 0;

	line = new TLine(0,1,xmaxJanChi2,1);
        line->SetLineStyle(3);
        line->Draw();

	c1_1->Modified();
   	c1->cd();

	//c1_2->Delete();
	//c1_2 = 0; 	

	c1_2 = new TPad("c1_2", "c1_2",0.01,0.01,0.99,0.49);
  	c1_2->Draw();
   	c1_2->cd();
   	c1_2->Range(-12.5,-1.375519,12.5,1.380519);
   	c1_2->SetBorderSize(2);
   	c1_2->SetFrameFillColor(0);
	c1_2->SetTickx(1);
        c1_2->SetTicky(1);
	c1_2->SetGridx(1);
        c1_2->SetGridy(1);


	TGraph* PValueGraph = new TGraph(n,xValue,PValueTab);
	PValueGraph->SetFillColor(1);
        PValueGraph->SetLineColor(4);
        PValueGraph->SetMarkerColor(4);
        PValueGraph->SetMarkerStyle(21);
        PValueGraph->SetMarkerSize(0.6);
        PValueGraph->SetLineColor(4);
        PValueGraph->GetYaxis()->SetTitle("p-value");
        PValueGraph->GetYaxis()->CenterTitle(true);
        PValueGraph->GetYaxis()->SetLabelSize(0.06);
        PValueGraph->GetXaxis()->SetLabelSize(0.06);
        PValueGraph->GetYaxis()->SetTitleSize(0.07);
        PValueGraph->GetYaxis()->SetTitleOffset(0.7);
        PValueGraph->GetXaxis()->SetLabelFont(42);
        PValueGraph->GetYaxis()->SetTitleFont(42);
        PValueGraph->GetYaxis()->SetLabelFont(42);
	PValueGraph->SetLineWidth(1);
        PValueGraph->Draw("AP");	


	PValueGraph->SetTitle("");
	//pt->Delete();
        //pt = 0;

	pt = new TPaveText(0,0,0,0,"blNDC");
        pt->SetName("title");
        pt->SetBorderSize(2);
        pt->SetFillColor(kWhite);
        //TText *text12 = pt->AddText("");
        pt->Draw();
		
	PValueGraph->GetYaxis()->SetRangeUser(yminPValue,ymaxPValue);


	PValueGraph->GetXaxis()->SetLimits(xminPValue,xmaxPValue);
	
	c1_2->Modified();
   	c1->cd();
   	c1->Modified();
   	c1->cd();
   	c1->SetSelected(c1);

	nomFichier = "JanChi2PValueVsVar";
        enregistrementPlots(nomDossierTGraph, nomFichier, EndCaps, 10000, c1);

	c1->Clear();
	

	Janchi2->Delete();
	Janchi2 = 0;
	PValueGraph->Delete();
	PValueGraph = 0;
	
/*
	line4->Delete();
	line4 = 0;
	line3->Delete();
        line3 = 0;
	line2->Delete();
        line2 = 0;
	line->Delete();
        line = 0;
	c1_1->Delete();
	c1_1 = 0;
	pt->Delete();
        pt = 0;
	textF->Delete();
	textF = 0;
	textL->Delete();
        textL = 0;
*/
	//delete c1;
	//c1 = 0;

	

	return 0;

}



