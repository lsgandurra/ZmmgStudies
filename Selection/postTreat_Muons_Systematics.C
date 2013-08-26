// Root headers
#include "TSystem.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TLatex.h"
// C++ headers
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
// RooFit headers
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
// Dipslay headers
#include "CMSStyle.C"
// namespaces
using namespace std;
using namespace RooFit;

int main(int argc, char *argv[])
{
  cout << "argc= " << argc << endl;
  for(int iarg = 0 ; iarg < argc; iarg++) cout << "argv[" << iarg << "]= " << argv[iarg] << endl;
  if( argc == 1 ){cerr << "arguments should be passed !! sample" << endl; return 1;}
  string sample = argv[1];
  int icatmin = 0;
	int icatmax = 1;
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
//  **************************************************************
//  **************************************************************
	cout << "###################################" << endl;
	cout << "###################################" << endl;
	cout << "### Post-treat Muon Systematics ###" << endl;
	cout << "###################################" << endl;
	cout << "###################################" << endl;
	CMSstyle();
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);
	TGaxis::SetMaxDigits(3);
	
	bool isMC = true;
	vector<string> categoryCut;
  vector<string> categoryName;
	vector<int> isEE;
	vector<int> ishighR9;
	// cat 0
	categoryCut.push_back("Photon_r9 < .95 && Photon_isEE");
  categoryName.push_back("EE_lowR9");
  isEE.push_back(1);
  ishighR9.push_back(0);
	// cat 1
  categoryCut.push_back("Photon_r9 > .95 && Photon_isEE");
  categoryName.push_back("EE_highR9");
  isEE.push_back(1);
  ishighR9.push_back(1);
	// cat 2
  categoryCut.push_back("Photon_r9 < .94 && Photon_isEB");
  categoryName.push_back("EB_lowR9");
  isEE.push_back(0);
  ishighR9.push_back(0);
	// cat 3
  categoryCut.push_back("Photon_r9 > .94 && Photon_isEB");
  categoryName.push_back("EB_highR9");
  isEE.push_back(0);
  ishighR9.push_back(1);
	// cat 4
  categoryCut.push_back("Photon_isEB");
  categoryName.push_back("EB_allR9");
  isEE.push_back(0);
  ishighR9.push_back(2);
	// cat 5
  categoryCut.push_back("Photon_isEE");
  categoryName.push_back("EE_allR9");
  isEE.push_back(1);
  ishighR9.push_back(2);
	int EndCaps = isEE[icatmin];
	int r9sup = ishighR9[icatmin];

	cout << "### Selection= " << categoryCut[icatmin] << endl;
	TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);
//	TFile *f = new TFile(Form("outfile_ethz_%i_%i_ALL.root", icatmin, icatmax));
//	TFile *f = new TFile(Form("outfile_data_MITregression_v04_%i_%i_ALL.root", icatmin, icatmax));
	TFile *f = new TFile(Form("outfile_mc_MITregression_v04_%i_%i_ALL.root", icatmin, icatmax));
	TTree *t = (TTree*)f->Get("miniTree");
	cout << "n= " << t->GetEntries() << endl;
	cout << "min= " << t->GetMinimum("Voigtian_mean") << endl;
	cout << "max= " << t->GetMaximum("Voigtian_mean") << endl;
	cout << "emin= " << (int)(floor(log10(fabs(t->GetMinimum("Voigtian_mean"))))) << endl;
	cout << "emax= " << (int)(floor(log10(fabs(t->GetMaximum("Voigtian_mean"))))) << endl;
	int exponent = min((int)(floor(log10(fabs(t->GetMinimum("Voigtian_mean"))))), (int)(floor(log10(fabs(t->GetMaximum("Voigtian_mean"))))));
	cout << "exponent= " << exponent << endl;
	int emin = floor(t->GetMinimum("Voigtian_mean") * pow(10., -exponent));
	int emax = ceil(t->GetMaximum("Voigtian_mean") * pow(10., -exponent));
	cout << "emin= " << emin << endl;
	cout << "emax= " << emax << endl;
	cout << "expo= " << exponent << endl;
	double max = t->GetMaximum("Voigtian_mean");
	double min = t->GetMinimum("Voigtian_mean");
	double range = (max - min)*1.5;
	double mean = (double)(max + min)/(double)(2.0);
//	double xmin = mean - (double)(range)/(double)(2.0);
//	double xmax = mean + (double)(range)/(double)(2.0);
	int nbins = 20;
	double xmin = emin * pow(10., exponent);
	double xmax = emax * pow(10., exponent);
	cout << "xmax= " << xmax << endl;
	cout << "xmin= " << xmin << endl;
//	if( icatmin == 1 ) xmax = -0.003;
	if( icatmin == 2 ) xmin = 0.015;
	if( icatmin == 3 ) xmax = 0.007;
	if( icatmin == 4 ) xmax = 0.01;
	double binwidth = (double)(xmax - xmin)/(double)(nbins);
	int binwidth_exponent = floor(log10(binwidth));
	double ebinwidth = binwidth * pow(10., -binwidth_exponent);
	std::ostringstream tempString;
	tempString << setprecision (2) << fixed << ebinwidth;
	string ebinwidth_string = tempString.str();
	cout << "binwidth= " << binwidth << endl;
	cout << "binwidth_exponent= " << binwidth_exponent << endl;
	cout << "ebinwidth= " << ebinwidth << endl;
	cout << "xmin= " << xmin << endl;
	cout << "xmax= " << xmax << endl;
	RooRealVar *Icat = new RooRealVar("Icat", "Icat", 0., 5.);
	RooRealVar *Voigtian_mean = new RooRealVar("Voigtian_mean", "fitted mean", xmin, xmax);
	RooRealVar *g_mean = new RooRealVar("g_mean", "g_mean", xmin, xmax);
	RooRealVar *g_sigma = new RooRealVar("g_sigam", "g_sigma", 0.0000001, range);
	RooGaussian *g = new RooGaussian("g", "g", *Voigtian_mean, *g_mean, *g_sigma);
	RooDataSet *d = new RooDataSet("d", "d", t, RooArgSet(*Voigtian_mean) );
	g->fitTo(*d);
	RooPlot *frame = (RooPlot*)Voigtian_mean->frame();
	frame->GetYaxis()->SetTitle(Form("Toys / (%s #times 10^{%i})", ebinwidth_string.c_str(), binwidth_exponent));
	d->plotOn(frame, Binning(nbins));
	g->plotOn(frame, LineColor(kRed));
	frame->Draw();


	TLatex latexLabel;
	latexLabel.SetTextSize(0.03);
	latexLabel.SetNDC();
	latexLabel.DrawLatex(0.25, 0.96, "CMS Private 2011, #sqrt{s} = 7 TeV");
  if(isMC == 1) latexLabel.DrawLatex(0.20, 0.88, "Simulation");
  if(isMC == 0) latexLabel.DrawLatex(0.20, 0.88, "Data, #int L = 4.89 fb^{-1}"); //lumi a changer !!!
  if(EndCaps == 0) latexLabel.DrawLatex(0.20, 0.83,"ECAL Barrel");
  if(EndCaps == 1) latexLabel.DrawLatex(0.20, 0.83,"ECAL Endcaps");
  if(r9sup == 0 && EndCaps == 0) latexLabel.DrawLatex(0.20, 0.78,"E_{T}^{#gamma} > 25 GeV, R_{9} < 0.94");
  if(r9sup == 0 && EndCaps == 1) latexLabel.DrawLatex(0.20, 0.78,"E_{T}^{#gamma} > 25 GeV, R_{9} < 0.95");
  if(r9sup == 1 && EndCaps == 0) latexLabel.DrawLatex(0.20, 0.78,"E_{T}^{#gamma} > 25 GeV, R_{9} > 0.94");
  if(r9sup == 1 && EndCaps == 1) latexLabel.DrawLatex(0.20, 0.78,"E_{T}^{#gamma} > 25 GeV, R_{9} > 0.95");
  if(r9sup == 2) latexLabel.DrawLatex(0.20, 0.78,"E_{T}^{#gamma} > 25 GeV, All R_{9}");
//  latexLabel.DrawLatex(0.17, 0.73, Form("Entries = %d",entriesInt));
  latexLabel.DrawLatex(0.55, 0.90, Form("#mu = %f #pm %f",g_mean->getVal(), g_mean->getError()));
  latexLabel.DrawLatex(0.55, 0.85, Form("#sigma = %f #pm %f",g_sigma->getVal(), g_sigma->getError()));
//  latexLabel.DrawLatex(0.61, 0.80, Form("#color[4]{width = %f^{+ %f}_{ %f}}",width->getVal(), width->getErrorHi(),width->getErrorLo()));

	c1->Print(Form("muSys_%s_%s.gif", sample.c_str(), categoryName[icatmin].c_str()));
	c1->Print(Form("muSys_%s_%s.C", sample.c_str(), categoryName[icatmin].c_str()));
	c1->Print(Form("muSys_%s_%s.root", sample.c_str(), categoryName[icatmin].c_str()));
	c1->Print(Form("muSys_%s_%s.pdf", sample.c_str(), categoryName[icatmin].c_str()));

	cout << "### Cleaning ###" << endl;
	t->SetDirectory(0);
	delete t;
	t = 0;
	f->Close();
	delete f;
	f = 0;

	g_mean->Delete();
	g_mean = 0;
	g_sigma->Delete();
	g_sigma = 0;

	g->Delete();
	g = 0;
	d->Delete();
	d = 0;

	Icat->Delete();
	Icat = 0;
	Voigtian_mean->Delete();
	Voigtian_mean = 0;

	delete c1;
	c1 = 0;

	return 0;
}



