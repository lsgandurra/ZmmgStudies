// Root headers
#include "TGaxis.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TChain.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TLatex.h"
#include "TClonesArray.h"
#include "TGraphAsymmErrors.h"
#include "TLegendEntry.h"
#include "TF1.h"
#include "TMath.h"
#include "TTree.h"
#include "TFile.h"
//#include "TObjectTable.h"
// C++ headers
#include <sstream>
#include <algorithm>
#include <math.h>
#include <iostream>
#include <string>
// RooFit headers
//#include "RooGlobalFunc.h"
//#include "RooChi2MCSModule.h"
#include "RooHist.h"
#include "RooArgList.h"
#include "RooFitResult.h"
#include "RooVoigtian.h"
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooDataSet.h"
#include "RooMsgService.h"
// Dipslay headers
#include "CMSStyle.C"
// namespaces
using namespace std;
using namespace RooFit;

void RangeEstimator3(double pourcentage, TChain * chain, TString temp, int Endcaps, double * MinRange, double * MaxRange);
void enregistrementPlots(string nomDossier, string nomFichier, int EndCaps, int iteration, TCanvas * c1);
Double_t effSigma(TH1 * hist);
Double_t chiSquare(RooPlot* plot_, char* pdfname, char* histname, int nFitParam, double* JanChi2, double* DegreesOfFreedom, double* pValue, int* fewBins);
double SigmaR(TF1* crystalBall, double Xmin, double Xmax);
double SigmaL(TF1* crystalBall, double Xmin, double Xmax);
void RooVoigtian2(double * mean_value, double * mean_error, double * sigma_value, double * sigma_error, double * width_value, double * width_error, Double_t * sigmaEff_value, double * sigmaR_value, double * sigmaL_value, double * ChiSquare, double * minLogLikelihood, TH1F * hh, int j, int EndCaps, TString temp, TChain* Tree_Data, double mean, double rms, TCanvas * c1, RooPlot* Erawframe, TF1 * f, string nomDossier, string nomFichier, double RangeMin, double RangeMax, double * JanChi2, double * DegreesOfFreedom, double * pValue, int r9sup, double MinVar, double MaxVar, string variableX, int isMC);

int main(int argc, char *argv[])
{
//	gObjectTable->Print();
  cout << "argc= " << argc << endl;
  for(int iarg = 0 ; iarg < argc; iarg++) cout << "argv[" << iarg << "]= " << argv[iarg] << endl;
  if( argc == 1 ){cerr << "arguments should be passed !! sample icatmin icatmax itoymin itoymax data/mc" << endl; return 1;}
  string sample = argv[1];
  int icatmin = 0;
  int icatmax = -99;
  if( argc > 2 )
  {
    std::stringstream ss ( argv[2] );
    ss >> icatmin;
  }
  if( argc > 3 )
  {
    std::stringstream ss ( argv[3] );
    ss >> icatmax;
  }
  int itoyMin = 1;
  int itoyMax = 100;
	if( argc > 4 )
	{
		std::stringstream ss ( argv[4] );
		ss >> itoyMin;
	}
	if( argc > 5 )
	{
		std::stringstream ss ( argv[5] );
		ss >> itoyMax;
	}
	string data="mc";
	if( argc > 6 )
	{
		data = argv[6];
	}
//  **************************************************************
//  **************************************************************
  cout << "########################" << endl;
  cout << "########################" << endl;
  cout << "### Muon Systematics ###" << endl;
  cout << "########################" << endl;
  cout << "########################" << endl;
  cout << "### initialize stuff" << endl;
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  CMSstyle();
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);
  int ntoy = itoyMax - itoyMin;
  string selectionBase = "isLooseMMG && Photon_Et > 25.";
  TTree* outtree = new TTree("miniTree", "miniTree");
  Int_t Itoy, Icat;
  Float_t Mean, Rms;
  Float_t Toys_Min, Toys_Max, Toys_Mean;
  Float_t Voigtian_mean, Voigtian_mean_error;
  Float_t Voigtian_sigma, Voigtian_sigma_error;
  Float_t Voigtian_width, Voigtian_width_error;
  Float_t Voigtian_sigmaEff, Voigtian_sigmaLeft, Voigtian_sigmaRight;
  Float_t Voigtian_chiSquare, Voigtian_minLL, Voigtian_chiSquareJan, Voigtian_dof, Voigtian_pvalue;
  outtree->Branch("Itoy", &Itoy, "Itoy/I");
  outtree->Branch("Icat", &Icat, "Icat/I");
  outtree->Branch("Mean", &Mean, "Mean/F");
  outtree->Branch("Rms", &Rms, "Rms/F");
  outtree->Branch("Toys_Min", &Toys_Min, "Toys_Min/F");
  outtree->Branch("Toys_Max", &Toys_Max, "Toys_Max/F");
  outtree->Branch("Toys_Mean", &Toys_Mean, "Toys_Mean/F");
  outtree->Branch("Voigtian_mean", &Voigtian_mean, "Voigtian_mean/F");
  outtree->Branch("Voigtian_mean_error", &Voigtian_mean_error, "Voigtian_mean_error/F");
  outtree->Branch("Voigtian_sigma", &Voigtian_sigma, "Voigtian_sigma/F");
  outtree->Branch("Voigtian_sigma_error", &Voigtian_sigma_error, "Voigtian_sigma_error/F");
  outtree->Branch("Voigtian_width", &Voigtian_width, "Voigtian_width/F");
  outtree->Branch("Voigtian_width_error", &Voigtian_width_error, "Voigtian_width_error/F");
  outtree->Branch("Voigtian_sigmaEff", &Voigtian_sigmaEff, "Voigtian_sigmaEff/F");
  outtree->Branch("Voigtian_sigmaLeft", &Voigtian_sigmaLeft, "Voigtian_sigmaLeft/F");
  outtree->Branch("Voigtian_sigmaRight", &Voigtian_sigmaRight, "Voigtian_sigmaRight/F");
  outtree->Branch("Voigtian_chiSquare", &Voigtian_chiSquare, "Voigtian_chiSquare/F");
  outtree->Branch("Voigtian_minLL", &Voigtian_minLL, "Voigtian_minLL/F");
  outtree->Branch("Voigtian_chiSquareJan", &Voigtian_chiSquareJan, "Voigtian_chiSquareJan/F");
  outtree->Branch("Voigtian_dof", &Voigtian_dof, "Voigtian_dof/F");
  outtree->Branch("Voigtian_pvalue", &Voigtian_pvalue, "Voigtian_pvalue/F");

  Itoy = Icat = 0;
  Mean = Rms = -99.;
  Toys_Min = Toys_Max = Toys_Mean = -99.;
  Voigtian_mean = Voigtian_mean_error = -99.;
  Voigtian_sigma = Voigtian_sigma_error = -99.;
  Voigtian_width = Voigtian_width_error = -99.;
  Voigtian_sigmaEff = Voigtian_sigmaLeft = Voigtian_sigmaRight = -99.;
  Voigtian_chiSquare = Voigtian_minLL = Voigtian_chiSquareJan = Voigtian_dof = Voigtian_pvalue = -99.;


//  **************************************************************
//  **************************************************************
  cout << "### define categories" << endl;
  vector<string> categoryCut;
  vector<string> categoryName;
  vector<pair<double, double > > pourcentage;
  vector<int> isEE;
  vector<int> ishighR9;
  categoryCut.push_back("Photon_r9 < .95 && Photon_isEE");
  categoryName.push_back("EE_lowR9");
//  pourcentage.push_back(make_pair(0.82, 0.88));
  pourcentage.push_back(make_pair(0.82, 0.88));
  isEE.push_back(1);
  ishighR9.push_back(0);
  categoryCut.push_back("Photon_r9 > .95 && Photon_isEE");
  categoryName.push_back("EE_highR9");
//  pourcentage.push_back(make_pair(0.7, 0.92));
  pourcentage.push_back(make_pair(0.77, 0.92));
  isEE.push_back(1);
  ishighR9.push_back(1);
  categoryCut.push_back("Photon_r9 < .94 && Photon_isEB");
  categoryName.push_back("EB_lowR9");
//  pourcentage.push_back(make_pair(0.9, 0.92));
  pourcentage.push_back(make_pair(0.87, 0.80));
  isEE.push_back(0);
  ishighR9.push_back(0);
  categoryCut.push_back("Photon_r9 > .94 && Photon_isEB");
  categoryName.push_back("EB_highR9");
//  pourcentage.push_back(make_pair(0.9, 0.76));
  pourcentage.push_back(make_pair(0.90, 0.82));
  isEE.push_back(0);
  ishighR9.push_back(1);
  categoryCut.push_back("Photon_isEB");
  categoryName.push_back("EB_allR9");
//  pourcentage.push_back(make_pair(0.9, 0.7));
  pourcentage.push_back(make_pair(0.81, 0.73));
  isEE.push_back(0);
  ishighR9.push_back(2);
  categoryCut.push_back("Photon_isEE");
  categoryName.push_back("EE_allR9");
//  pourcentage.push_back(make_pair(0.9, 0.93));
  pourcentage.push_back(make_pair(0.86, 0.96));
  isEE.push_back(1);
  ishighR9.push_back(2);
//  categoryCut.push_back("1");
//  categoryName.push_back("EBEE_allR9");
//  **************************************************************
//  **************************************************************
  cout << "### Loading files" << endl;
	string protocol = "dcap://ccdcapcms.in2p3.fr:22125";
  string path= Form("/pnfs/in2p3.fr/data/cms/t2data/store/user/obondu/2011/%s/miniTrees/", data.c_str());
  if( icatmax == -99 ) icatmax = categoryCut.size();
  cout << "### Loop over " << icatmax - icatmin << " categories" << endl;
  for( int icat=icatmin ; icat < icatmax ; icat++)
  {
    // init
    string selection = selectionBase + " && " + categoryCut[icat];
    cout << "###\t### icat= " << icat << "\tselection = " << selection << endl;
    string limits = "(25,-0.5,0.5)";
    string limits_fine = "(1000,-0.5,0.5)";
    vector<unsigned long> n_itoy;
    n_itoy.push_back(1.0);
    //  **************************************************************
    //  **************************************************************
    // actual loop
    cout << "###\t### Loop over " << ntoy << " toys to find min/max/mean" << endl;
    for(int itoy = itoyMin; itoy < itoyMax ; itoy++)
    {
      TChain *chain_itoy = new TChain("miniTree");
			if( data == "mc")
			{
//	      chain_itoy->Add(Form("%s%sminiTree_v31_FSR_DYToMuMu_rochcor_ethz_toy%s%i_partALL.root", protocol.c_str(), path.c_str(), (itoy < 10)?"0":"", itoy));
//	      chain_itoy->Add(Form("%s%sminiTree_v31_nonFSR_DYToMuMu_rochcor_ethz_toy%s%i_partALL.root", protocol.c_str(), path.c_str(), (itoy < 10)?"0":"", itoy));
	      chain_itoy->Add(Form("%s%sminiTree_v31_FSR_DYToMuMu_rochcor_MITregression_toy%s%i_partALL.root", protocol.c_str(), path.c_str(), (itoy < 10)?"0":"", itoy));
	      chain_itoy->Add(Form("%s%sminiTree_v31_nonFSR_DYToMuMu_rochcor_MITregression_toy%s%i_partALL.root", protocol.c_str(), path.c_str(), (itoy < 10)?"0":"", itoy));
			} else {
	      chain_itoy->Add(Form("%s%sminiTree_v31_Run2011A-16Jan2012-v1_rochcor_MITregression_toy%s%i_partALL.root", protocol.c_str(), path.c_str(), (itoy < 10)?"0":"", itoy));
	      chain_itoy->Add(Form("%s%sminiTree_v31_Run2011B-16Jan2012-v1_rochcor_MITregression_toy%s%i_partALL.root", protocol.c_str(), path.c_str(), (itoy < 10)?"0":"", itoy));
			}
      unsigned long n_itoy_temp = chain_itoy->GetEntries(selection.c_str());
			cout << "###\t### itoy= " << itoy << "\tentries= " << n_itoy_temp << endl;
      TH1F* h_itoy = new TH1F(Form("hist_itoy_%i", itoy), Form("hist_itoy_%i", itoy), 25, -0.5, 0.5);
      TH1F* h_itoy_fine = new TH1F(Form("hist_itoy_fine_%i", itoy), Form("hist_itoy_fine_%i", itoy), 1000, -0.5, 0.5);
			string temp_histName = Form("hist_itoy_%i", itoy);
      string temp_histName_fine = Form("hist_itoy_fine_%i", itoy);
      string var_itoy = "mmg_s >> " + temp_histName;
      string var_itoy_fine = "mmg_s >> " + temp_histName_fine;
      chain_itoy->Draw(var_itoy.c_str(), Form("(%s)", selection.c_str()));
      chain_itoy->Draw(var_itoy_fine.c_str(), Form("(%s)", selection.c_str()));
      n_itoy.push_back(n_itoy_temp);
      h_itoy->Scale((double)1./(double)n_itoy[itoy]);
      double mean_value = 0.;
      double mean_error = 0.;
      double sigma_value = 0.;
      double sigma_error = 0.;
      double width_value = 0.;
      double width_error = 0.;
      Double_t sigmaEff_value = 0.;
      double sigmaR_value = 0.;
      double sigmaL_value = 0.;
      double ChiSquare = 0.;
      double minLogLikelihood = 0.;
      int j = itoy;
      int EndCaps = isEE[icat];
      TString temp = selection.c_str();
      double mean, rms; // not used
//      TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);;
      RooPlot* Erawframe = new RooPlot();;
      TF1 * f = new TF1();
      string nomDossier = "plots_" + sample + "_";
      string temp_ = Form("category_%i", icat);
      string nomFichier = temp_ + "_itoy_";
      double RangeMin = 0.;
      double RangeMax = 0.;
      double JanChi2 = 0.;
      double DegreesOfFreedom = 0.;
      double pValue = 0.;
      int r9sup = ishighR9[icat];
      string variableX; // not used
      int isMC = (data=="mc")?1:0;
     double pourcentage_ = (isMC==1) ? pourcentage[icat].second : pourcentage[icat].first;

      RangeEstimator3(pourcentage_, chain_itoy, temp, EndCaps, &RangeMin, &RangeMax);
      RooVoigtian2( &mean_value, &mean_error, &sigma_value, &sigma_error, &width_value, &width_error, &sigmaEff_value, &sigmaR_value, &sigmaL_value, &ChiSquare, &minLogLikelihood, h_itoy_fine, j, EndCaps, temp, chain_itoy, mean, rms, c1, Erawframe, f, nomDossier, nomFichier, RangeMin, RangeMax, &JanChi2, &DegreesOfFreedom, &pValue, r9sup, 1.0, 1.0, variableX, isMC);
      Itoy = itoy;
      Icat = icat;
      Mean = h_itoy_fine->GetMean();
      Rms =  h_itoy_fine->GetRMS();
      Voigtian_mean = mean_value;
      Voigtian_mean_error = mean_error;
      Voigtian_sigma = sigma_value;
      Voigtian_sigma_error = sigma_error;
      Voigtian_width = width_value;
      Voigtian_width_error = width_error;
      Voigtian_sigmaEff = sigmaEff_value;
      Voigtian_sigmaLeft = sigmaL_value;
      Voigtian_sigmaRight = sigmaR_value;
      Voigtian_chiSquare = ChiSquare;
      Voigtian_minLL = minLogLikelihood;
      Voigtian_chiSquareJan = JanChi2;
      Voigtian_dof = DegreesOfFreedom;
      Voigtian_pvalue = pValue;
      outtree->Fill();


			Erawframe->Delete();
			Erawframe = 0;
			f->Delete();
			f  = 0;
			h_itoy->Delete();
			h_itoy = 0;
			h_itoy_fine->Delete();
			h_itoy_fine = 0;
			chain_itoy->Delete();
			chain_itoy = 0;
//			gObjectTable->Print();
		}
	}

  cout << "### Writing stuff out " << endl;
  TFile* outfile = new TFile(Form("outfile_%s.root", sample.c_str()), "recreate");
  TTree *tree = (TTree*)outtree->CopyTree("");
  outfile->Write();
  tree->SetDirectory(0);
  outtree->SetDirectory(0);
  delete tree;
  delete outtree;
  tree = 0;
  outtree = 0;
  outfile->Close();
  delete outfile;
  outfile = 0;

	
	delete c1;
	c1 = 0;
//	gObjectTable->Print();
	return 0;
}

// Adapted from Louis' code, http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/LouisSgandurra/ZtoMuMuGamma/mmg_s_Fits_ApprovalVersion.C?revision=1.1&view=markup
void RooVoigtian2(double * mean_value, double * mean_error, double * sigma_value, double * sigma_error, double * width_value, double * width_error, Double_t * sigmaEff_value, double * sigmaR_value, double * sigmaL_value, double * ChiSquare, double * minLogLikelihood, TH1F * hh, int j, int EndCaps, TString temp, TChain* Tree_Data, double mean, double rms, TCanvas * c1, RooPlot* Erawframe, TF1 * f, string nomDossier, string nomFichier, double RangeMin, double RangeMax, double * JanChi2, double * DegreesOfFreedom, double * pValue, int r9sup, double MinVar, double MaxVar, string variableX, int isMC)
{
	RooRealVar * mmg_s = new RooRealVar("mmg_s", "s", -0.5, 0.5);
	RooRealVar * Photon_SC_Eta = new RooRealVar("Photon_SC_Eta", "Photon_SC_Eta", -3.0, 3.0, "");
	RooRealVar * Photon_SC_Phi = new RooRealVar("Photon_SC_Phi", "Photon_SC_Phi", -3.5, 3.5, "");
	RooRealVar * Photon_r9 = new RooRealVar("Photon_r9", "Photon_r9", 0.0, 1.0, "");
	RooRealVar * isLooseMMG = new RooRealVar("isLooseMMG", "isLooseMMG", 0.0, 1.0, "");
	RooRealVar * isMultipleCandidate = new RooRealVar("isMultipleCandidate", "isMultipleCandidate", 0.0, 1.0, "");
	RooRealVar * Photon_Et = new RooRealVar("Photon_Et", "Photon_Et", 0.0, 250.0, "");
	RooRealVar * Photon_SC_rawEt = new RooRealVar("Photon_SC_rawEt", "Photon_SC_rawEt", 0.0, 250.0, "");
	RooRealVar * Photon_E = new RooRealVar("Photon_E", "Photon_E", 0.0, 1000.0, "");
	RooRealVar * Photon_SC_brem = new RooRealVar("Photon_SC_brem", "Photon_SC_brem", 0.0, 15.0, "");
	RooRealVar * weight_pileUp = new RooRealVar("weight_pileUp", "weight_pileUp", 0.0, 100);
	RooRealVar * Photon_isEB = new RooRealVar("Photon_isEB", "Photon_isEB", 0.0, 1.0, "");
	RooRealVar * Photon_isEE = new RooRealVar("Photon_isEE", "Photon_isEE", 0.0, 1.0, "");

	RooArgSet *ntplVars = new RooArgSet(*mmg_s, *Photon_SC_Eta, *Photon_SC_Phi, *Photon_r9, *Photon_Et, *Photon_SC_rawEt, *Photon_E, *Photon_SC_brem, *isLooseMMG);
	ntplVars->add(*isMultipleCandidate);
	ntplVars->add(*weight_pileUp);
	ntplVars->add(*Photon_isEB);
	ntplVars->add(*Photon_isEE);

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
	int RightBinning = 25;
	for(int i = 0; i<4; i++)
	{
		Erawframe->Clear();
		//Data_subset = (RooDataSet*)Data->reduce(*ntplVars, temp);
		Voigtian = new RooVoigtian("Voigtian","Voigtian",*mmg_s,*meanV,*sigma,*width);


		Erawframe = mmg_s->frame();
		Data_subset->plotOn(Erawframe,Name("myhist"),Binning(RightBinning),DataError(RooAbsData::SumW2));
		//Voigtian->fitTo(*Data_subset, Range(RangeMin, RangeMax),Minos(true));
		//res = Voigtian->fitTo(*Data_subset, Range(RangeMin, RangeMax),Minos(true),Save(),SumW2Error(kFALSE));

		res = Voigtian->fitTo(*Data_subset, Range(RangeMin, RangeMax),Save(),SumW2Error(kTRUE));
//		res = Voigtian->fitTo(*Data_subset, Range(RangeMin, RangeMax));
//		res->Print();
		minNll = res->minNll();

		Voigtian->plotOn(Erawframe,Name("mycurve"));
		Erawframe->Draw();


		Chi2J = chiSquare(Erawframe,(char *)"mycurve",(char *)"myhist",3, JanChi2, DegreesOfFreedom, pValue, &fewBins);
		if(fewBins == 0) break; 


	}
//	cout << endl << "ICI !!!!!!!!!!!!!!!!!! "<<Data_subset->weight() << endl ;

	f = Voigtian->asTF( RooArgList(*mmg_s) );
	*mean_value  = meanV->getVal();
	*mean_error  = meanV->getError();
	*sigma_value = sigma->getVal();
	*sigma_error = sigma->getError();
	*width_value = width->getVal();
	*width_error = width->getError();
	*sigmaEff_value = effSigma(hh);
	*sigmaR_value = SigmaR(f, 0.0, 2.0);
	*sigmaL_value = SigmaL(f, 0.0, 2.0);
	*ChiSquare = Erawframe->chiSquare();
	*minLogLikelihood = minNll;
/*
	cout << "mean_value= " << *mean_value << endl;
	cout << "mean_error= " << *mean_error << endl;
	cout << "sigmaEff_value= " << *sigmaEff_value << endl;
	cout << "sigmaR_value= " << *sigmaR_value << endl;
	cout << "sigmaL_value= " << *sigmaL_value << endl;
	cout << "ChiSquare= " << *ChiSquare << endl;
	cout << "minLogLikelihood= " << *minLogLikelihood << endl;
*/
	double entries = hh->GetEntries();

	int entriesInt = (int) entries;

	
	double fxmax = f->GetMaximumX(0.8,1.2,1.E-10,100,false);
	//cout<<"////////// ---- FXMAX BG = "<<f->GetMaximumX(0.8,1.2,1.E-10,100,false)<<" ---- //////////"<<endl;


	TLatex latexLabel;
	latexLabel.SetTextSize(0.026);
	latexLabel.SetNDC();
	latexLabel.DrawLatex(0.13, 0.96, "CMS Preliminary 2011, #sqrt{s} = 7 TeV");
	if(isMC == 1) latexLabel.DrawLatex(0.17, 0.88, "Simulation");
	if(isMC == 0) latexLabel.DrawLatex(0.17, 0.88, "Data, #int L = 4,89 fb^{-1}"); //lumi a changer !!!
	if(EndCaps == 0) latexLabel.DrawLatex(0.17, 0.83,"ECAL Barrel");
	if(EndCaps == 1) latexLabel.DrawLatex(0.17, 0.83,"ECAL Endcaps");
	if(r9sup == 0 && EndCaps == 0) latexLabel.DrawLatex(0.17, 0.78,"E_{T}^{#gamma} > 25 GeV, r9 < 0,94");
	if(r9sup == 0 && EndCaps == 1) latexLabel.DrawLatex(0.17, 0.78,"E_{T}^{#gamma} > 25 GeV, r9 < 0,95");
	if(r9sup == 1 && EndCaps == 0) latexLabel.DrawLatex(0.17, 0.78,"E_{T}^{#gamma} > 25 GeV, r9 > 0,94");
	if(r9sup == 1 && EndCaps == 1) latexLabel.DrawLatex(0.17, 0.78,"E_{T}^{#gamma} > 25 GeV, r9 > 0,95");
	if(r9sup == 2) latexLabel.DrawLatex(0.17, 0.78,"E_{T}^{#gamma} > 25 GeV, All r9");
	latexLabel.DrawLatex(0.17, 0.73, Form("Entries = %d",entriesInt));
	latexLabel.DrawLatex(0.61, 0.90, Form("#color[4]{mean = %f^{+ %f}_{ %f}}",meanV->getVal(), meanV->getErrorHi(),meanV->getErrorLo()));
	latexLabel.DrawLatex(0.61, 0.85, Form("#color[4]{#sigma = %f^{+ %f}_{ %f}}",sigma->getVal(), sigma->getErrorHi(),sigma->getErrorLo()));
	latexLabel.DrawLatex(0.61, 0.80, Form("#color[4]{width = %f^{+ %f}_{ %f}}",width->getVal(), width->getErrorHi(),width->getErrorLo()));	


	enregistrementPlots(nomDossier, nomFichier, EndCaps, j, c1);
	c1->Clear();



	////////// chi2 Pulls & Residuals //////////
/*
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
	residuals->GetXaxis()->SetTitle("s_{RECO}");
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
	if(isMC == 0) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV, #int L = 4,89 fb^{-1}");
	//textResid->DrawLatex(0.16, 0.80, EndCapsR9Chain->c_str());
	//textResid->DrawLatex(0.16, 0.75, tempLegChain->c_str());
	residuals->GetXaxis()->SetLimits(-0.5,0.5);
	
	enregistrementPlots(nomDossier, nomFichier, EndCaps, j, c1);

	c1->Clear();

	pulls->Draw("AP");
	nomFichier = "Chi2Pulls";
	
	pulls->GetYaxis()->SetTitle("#chi^{2} Pulls");
	pulls->GetXaxis()->SetLabelFont(42);
	pulls->GetXaxis()->SetTitleFont(42);
	pulls->GetXaxis()->SetLabelSize(0.03);
	pulls->GetXaxis()->SetTitle("s_{RECO}");
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
	if(isMC == 0) textResid->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV, #int L = 4,89 fb^{-1}");
	//textResid->DrawLatex(0.16, 0.80, EndCapsR9Chain->c_str());
	//textResid->DrawLatex(0.16, 0.75, tempLegChain->c_str());
	pulls->GetXaxis()->SetLimits(-0.5,0.5);	
	c1->Clear();
*/
	//delete EndCapsR9Chain;
	//delete tempLegChain;
	res->Delete();
	res = 0;
//	mcs->Delete();
//	frame->Delete();
//	residuals->Delete();
//	pulls->Delete();
//	textResid->Delete();
//	lineResid->Delete();
	Data_subset->Delete();
//	delete Data_subset;
	Data_subset = 0;
	//Data->Delete();
	ntplVars->Delete();
//	delete ntplVars;
	ntplVars = 0;
	mmg_s->Delete();
////	delete mmg_s;
	mmg_s = 0;
	Photon_SC_Eta->Delete();
//	delete Photon_SC_Eta;
	Photon_SC_Eta = 0;
	Photon_SC_Phi->Delete();
//	delete Photon_SC_Phi;
	Photon_SC_Phi = 0;
	Photon_r9->Delete();
//	delete Photon_r9;
	Photon_r9 = 0;
	isLooseMMG->Delete();
//	delete isLooseMMG;
	isLooseMMG = 0;
	isMultipleCandidate->Delete();
//	delete isMultipleCandidate;
	isMultipleCandidate = 0;
	Photon_Et->Delete();
//	delete Photon_Et;
	Photon_Et = 0;
	Photon_SC_rawEt->Delete();
//	delete Photon_SC_rawEt;
	Photon_SC_rawEt = 0;
	Photon_E->Delete();
//	delete Photon_E;
	Photon_E = 0;
	Photon_SC_brem->Delete();
//	delete Photon_SC_brem;
	Photon_SC_brem = 0;
	weight_pileUp->Delete();
//	delete weight_pileUp;
	weight_pileUp = 0;
	Photon_isEB->Delete();
//	delete Photon_isEB;
	Photon_isEB = 0;
	Photon_isEE->Delete();
//	delete Photon_isEE;
	Photon_isEE = 0;
	meanV->Delete();
//	delete meanV;
	meanV = 0;
	sigma->Delete();
//	delete sigma;
	sigma = 0;
	width->Delete();
//	delete width;
	width = 0;
	Voigtian->Delete();
//	delete Voigtian;
	Voigtian = 0;
	return;
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

//#if ROOT_VERSION_CODE >= ROOT_VERSION(4,0,1)
  curve->GetPoint(0,xstart,y) ;
  curve->GetPoint(curve->GetN()-1,xstop,y) ;
/*
#else
  const_cast<RooCurve*>(curve)->GetPoint(0,xstart,y) ;
  const_cast<RooCurve*>(curve)->GetPoint(curve->GetN() - 1,xstop,y) ;
#endif
*/
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
    if(y != 0 && y < 35.0)
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
	hist = 0;
	curve = 0;

  return chisq / (nbin - nFitParam) ;
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
//	cout << "nb= " << nb << "\txmax= " << xmax << "\txmin= " << xmin << "\tave= " << ave << "\trms= " << rms << "\tnentries= " << hist->GetEntries() <<  endl;
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

	xaxis = 0;  
  return widmin;
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
        c1->Print(Form("%s%s.png",nomDossier.c_str(),nomFichier.c_str()));
        return;
}

void RangeEstimator3(double pourcentage, TChain * chain, TString temp, int Endcaps, double * MinRange, double * MaxRange)
{
	TChain * ReducedChain = (TChain *) chain->CopyTree(temp);
	float mmg_s;
	ReducedChain->SetBranchAddress("mmg_s",&mmg_s);
	vector <float> ErecoOverEtrueVector;
	for (int ievt = 0 ; ievt < ReducedChain->GetEntries() ; ievt++)
	{
		ReducedChain->GetEntry(ievt);
		ErecoOverEtrueVector.push_back(mmg_s);
	}
	sort(ErecoOverEtrueVector.begin(), ErecoOverEtrueVector.end());
	size_t interval_entries = TMath::Ceil(pourcentage * ErecoOverEtrueVector.size());
	vector<float>::iterator lower = ErecoOverEtrueVector.begin();
	vector<float>::iterator upper = ErecoOverEtrueVector.begin() + interval_entries - 1; 
	double dx = *upper - *lower;
	for(vector<float>::iterator first = lower, last = upper; last < ErecoOverEtrueVector.end(); first++, last++)
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
	ReducedChain = 0;
	return;
}




