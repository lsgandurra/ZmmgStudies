#include <TLorentzVector.h>
#include <TVector3.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TFormula.h>
#include <TF1.h>
#include <TF2.h>
#include <TSystem.h>
#include <TClonesArray.h>
#include <TLeaf.h>
#include <TChain.h>
#include <TObject.h>
#include <string.h>
#include <algorithm>
#include <TBranch.h>
#include <TString.h>
#include <TBits.h>
#include <TMath.h>
#include "TROOT.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <TLatex.h>
#include <THStack.h>
#include <TLegendEntry.h>
#include <TMinuit.h>
#include <TPaveStats.h>
#include <TCut.h>

//#include "DrawDataMC.h"

//#include "/sps/cms/obondu/CMSSW_4_1_2/src/UserCode/IpnTreeProducer/interface/TRootBardak.h"
//#include "/sps/cms/obondu/CMSSW_4_1_2/src/UserCode/IpnTreeProducer/interface/TRootBeamSpot.h"
//#include "/sps/cms/obondu/CMSSW_4_1_2/src/UserCode/IpnTreeProducer/interface/TRootCluster.h"
//#include "/sps/cms/obondu/CMSSW_4_1_2/src/UserCode/IpnTreeProducer/interface/TRootEcalRecHit.h"
//#include "/sps/cms/obondu/CMSSW_4_1_2/src/UserCode/IpnTreeProducer/interface/TRootElectron.h"
//#include "/sps/cms/obondu/CMSSW_4_1_2/src/UserCode/IpnTreeProducer/interface/TRootEvent.h"
//#include "/sps/cms/obondu/CMSSW_4_1_2/src/UserCode/IpnTreeProducer/interface/TRootJet.h"
//#include "/sps/cms/obondu/CMSSW_4_1_2/src/UserCode/IpnTreeProducer/interface/TRootMCParticle.h"
//#include "/sps/cms/obondu/CMSSW_4_1_2/src/UserCode/IpnTreeProducer/interface/TRootMCPhoton.h"
//#include "/sps/cms/obondu/CMSSW_4_1_2/src/UserCode/IpnTreeProducer/interface/TRootMET.h"
//#include "/sps/cms/obondu/CMSSW_4_1_2/src/UserCode/IpnTreeProducer/interface/TRootMuon.h"
//#include "/sps/cms/obondu/CMSSW_4_1_2/src/UserCode/IpnTreeProducer/interface/TRootParticle.h"
//#include "/sps/cms/obondu/CMSSW_4_1_2/src/UserCode/IpnTreeProducer/interface/TRootPhoton.h"
//#include "/sps/cms/obondu/CMSSW_4_1_2/src/UserCode/IpnTreeProducer/interface/TRootRun.h"
//#include "/sps/cms/obondu/CMSSW_4_1_2/src/UserCode/IpnTreeProducer/interface/TRootSignalEvent.h"
//#include "/sps/cms/obondu/CMSSW_4_1_2/src/UserCode/IpnTreeProducer/interface/TRootSuperCluster.h"
//#include "/sps/cms/obondu/CMSSW_4_1_2/src/UserCode/IpnTreeProducer/interface/TRootTopTop.h"
//#include "/sps/cms/obondu/CMSSW_4_1_2/src/UserCode/IpnTreeProducer/interface/TRootTrack.h"
//#include "/sps/cms/obondu/CMSSW_4_1_2/src/UserCode/IpnTreeProducer/interface/TRootVertex.h"

#include "CMSStyle.C"

using namespace std;


///sps/cms/hbrun/miniTree41X/do_vertexWeight.C

//void DrawDataMCplot(TTree *Data_miniTree, TTree *FSR_DYToMuMu_miniTree, TTree *nonFSR_DYToMuMu_miniTree, TTree *QCDMu_miniTree, TTree *TTJets_miniTree, TTree *WJetsToLNu_miniTree, string var, string pic, string limits, string cut, string name, string Title, bool inlog, bool drawUnderOverFsubleading, TCanvas *c1, bool doFit){
void DrawDataMCplot(TTree *Data_miniTree, TTree *FSR_DYToMuMu_miniTree, TTree *nonFSR_DYToMuMu_miniTree, TTree *TTJets_miniTree, string var, string pic, string limits, string cut, string name, string Title, bool inlog, bool drawUnderOverFsubleading, TCanvas *c1, double integratedLuminosity, bool doFit){

	CMSstyle();
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);

  // Get Histo_Data from eventTree
  TH1F *Histo_Data_temp = new TH1F();
  string variable_Data = var + ">>Histo_Data_temp" + limits;
  Data_miniTree->Draw(variable_Data.c_str(), cut.c_str());
  TH1F *Histo_Data = (TH1F*)gDirectory->Get("Histo_Data_temp");
  c1->Clear();

  // Get Histo_FSR_DYToMuMu from eventTree
  TH1F *Histo_FSR_DYToMuMu_temp = new TH1F();
  string variable_FSR_DYToMuMu = var + ">>Histo_FSR_DYToMuMu_temp" + limits;
//	string cut_FSR_DYToMuMu = "(" + cut + ") * weight_DYToMuMu(nGenVertices+1)";
	string cut_FSR_DYToMuMu = "(" + cut + ") * weight_pileUp * weight_Xsection";
//	string cut_FSR_DYToMuMu = "(" + cut + ") * weight_Xsection";
  FSR_DYToMuMu_miniTree->Draw(variable_FSR_DYToMuMu.c_str(), cut_FSR_DYToMuMu.c_str());
  TH1F *Histo_FSR_DYToMuMu = (TH1F*)gDirectory->Get("Histo_FSR_DYToMuMu_temp");
  c1->Clear();


  // Get Histo_nonFSR_DYToMuMu from eventTree
  TH1F *Histo_nonFSR_DYToMuMu_temp = new TH1F();
  string variable_nonFSR_DYToMuMu = var + ">>Histo_nonFSR_DYToMuMu_temp" + limits;
	string cut_nonFSR_DYToMuMu = "(" + cut + ") * weight_pileUp * weight_Xsection";
//	string cut_nonFSR_DYToMuMu = "(" + cut + ") * weight_Xsection";
//	string cut_nonFSR_DYToMuMu = "(" + cut + ") * weight_DYToMuMu(nGenVertices+1)";
//	string cut_nonFSR_DYToMuMu = cut;
  nonFSR_DYToMuMu_miniTree->Draw(variable_nonFSR_DYToMuMu.c_str(), cut_nonFSR_DYToMuMu.c_str());
  TH1F *Histo_nonFSR_DYToMuMu = (TH1F*)gDirectory->Get("Histo_nonFSR_DYToMuMu_temp");
  c1->Clear();


  // Get Histo_TTJets from eventTree
  TH1F *Histo_TTJets_temp = new TH1F();
  string variable_TTJets = var + ">>Histo_TTJets_temp" + limits;
	string cut_TTJets = "(" + cut + ") * weight_pileUp * weight_Xsection";
//	string cut_TTJets = "(" + cut + ") * weight_TTJets(nGenVertices+1)";
//	string cut_TTJets = cut;
  TTJets_miniTree->Draw(variable_TTJets.c_str(), cut_TTJets.c_str());
  TH1F *Histo_TTJets = (TH1F*)gDirectory->Get("Histo_TTJets_temp");
  c1->Clear();

/*
  // Get Histo_WJetsToLNu from eventTree
  TH1F *Histo_WJetsToLNu_temp = new TH1F();
  string variable_WJetsToLNu = var + ">>Histo_WJetsToLNu_temp" + limits;
	string cut_WJetsToLNu = "(" + cut + ") * weight_pileUp * weight_Xsection";
//	string cut_WJetsToLNu = "(" + cut + ") * weight_WJetsToLNu(nGenVertices+1)";
//	string cut_WJetsToLNu = cut;
  WJetsToLNu_miniTree->Draw(variable_WJetsToLNu.c_str(), cut_WJetsToLNu.c_str());
  TH1F *Histo_WJetsToLNu = (TH1F*)gDirectory->Get("Histo_WJetsToLNu_temp");
  c1->Clear();
*/
/*
  // Get Histo_QCDMu from eventTree
  TH1F *Histo_QCDMu_temp = new TH1F();
  string variable_QCDMu = var + ">>Histo_QCDMu_temp" + limits;
	string cut_QCDMu = "(" + cut + ") * weight_pileUp * weight_Xsection";
//	string cut_QCDMu = "(" + cut + ") * weight_QCDMu(nGenVertices+1)";
//	string cut_QCDMu = cut;
  QCDMu_miniTree->Draw(variable_QCDMu.c_str(), cut_QCDMu.c_str());
  TH1F *Histo_QCDMu = (TH1F*)gDirectory->Get("Histo_QCDMu_temp");
  c1->Clear();
*/
  // Get the number of entries for further normalization
//  double a = Histo_Data->Integral();
/*
  double b_FSR_DYToMuMu = Histo_FSR_DYToMuMu->Integral();
  if( (a==0.0) || (b_FSR_DYToMuMu==0.0) ){
    cout << "no entries to plots" <<endl;
    return; 
  }*/
  // Normalize
  Histo_Data->Sumw2(); // In order to have the correct error bars on data after renormalization
  // // Normalize MC and Data to 1
  //Histo_Data->Scale((double)((double)1.0/(double)a));
  //Histo_MC->Scale((double)((double)1.0/(double)b));
  // // Normalize MC to Data number of entries
//  double integratedLuminosity = 510.8801;
//	double integratedLuminosity = 203.9036;
//  double integratedLuminosity = 714.78373;
//  double integratedLuminosity = (214.439 + 663.204);
//  double integratedLuminosity = 216.122 + 924.829 + 368.037 + 658.886 + 1.842*1000;
//  double integratedLuminosity = 975.94723;
//	double integratedLuminosity = 1091.14253;
//	double integratedLuminosity = 1078.19387;

//  double XSectionFSR_DYToMuMu = 1300.0 * 1.2416;
  double XSectionFSR_DYToMuMu = 1626.0;
//  double XSectionFSR_DYToMuMu = 1300.0 * 1.2416;
  double XSectionnonFSR_DYToMuMu = 1626.0;
//  double XSectionTTJets = 94.0;
  double XSectionTTJets = 94.76;
//	double XSectionWJetsToLNu = 7899.0;
//	double XSectionQCDMu = 349988.0;

//  double InitialNumberFSR_DYToMuMu = 2148325.0;
//  double InitialNumbernonFSR_DYToMuMu = 2148325.0;
  double InitialNumberFSR_DYToMuMu = 29743564.0;
  double InitialNumbernonFSR_DYToMuMu = 29743564.0;
//  double InitialNumberTTJets = 1089625.0;
  double InitialNumberTTJets = 3701947.0;
//	double InitialNumberWJetsToLNu = 5413258.0;
//	double InitialNumberQCDMu = 8797418.0;

//  Histo_FSR_DYToMuMu->Scale((double)(  (double)((double)(XSectionFSR_DYToMuMu) / (double)(InitialNumberFSR_DYToMuMu)) * (double)integratedLuminosity));
//  Histo_nonFSR_DYToMuMu->Scale((double)(  (double)((double)(XSectionnonFSR_DYToMuMu) / (double)(InitialNumbernonFSR_DYToMuMu)) * (double)integratedLuminosity));
//  Histo_TTJets->Scale((double)(  (double)((double)(XSectionTTJets) / (double)(InitialNumberTTJets)) * (double)integratedLuminosity));
//  Histo_WJetsToLNu->Scale((double)(  (double)((double)(XSectionWJetsToLNu) / (double)(InitialNumberWJetsToLNu)) * (double)integratedLuminosity));
//  Histo_QCDMu->Scale((double)(  (double)((double)(XSectionQCDMu) / (double)(InitialNumberQCDMu)) * (double)integratedLuminosity));
  // Adding histograms for binned samples
//  Histo_QCD_Pt15->Add(Histo_QCD_Pt30);
//  Histo_QCD_Pt15->Add(Histo_QCD_Pt80);
//  Histo_QCD_Pt15->Add(Histo_QCD_Pt170);
//  Histo_QCD_Pt15->Add(Histo_QCD_Pt300);
//  Histo_QCD_Pt15->Add(Histo_QCD_Pt470);

	double integral_data = Data_miniTree->GetEntries(cut.c_str());
	FSR_DYToMuMu_miniTree->Draw("weight_pileUp >> temp_pileup(100,0,100)", cut.c_str());
  FSR_DYToMuMu_miniTree->Draw("weight_Xsection >> temp_weight(100,0,10)", cut.c_str());
	double integral_signal = FSR_DYToMuMu_miniTree->GetEntries(cut.c_str()) * (((TH1F*)gDirectory->Get("temp_weight"))->GetMean()) * (((TH1F*)gDirectory->Get("temp_pileup"))->GetMean());
	double integral_bg = 0;
	nonFSR_DYToMuMu_miniTree->Draw("weight_pileUp >> temp_pileup(100,0,100)", cut.c_str());
  nonFSR_DYToMuMu_miniTree->Draw("weight_Xsection >> temp_weight(100,0,10)", cut.c_str());
	integral_bg += nonFSR_DYToMuMu_miniTree->GetEntries(cut.c_str()) * (((TH1F*)gDirectory->Get("temp_weight"))->GetMean()) * (((TH1F*)gDirectory->Get("temp_pileup"))->GetMean());
	TTJets_miniTree->Draw("weight_pileUp >> temp_pileup(100,0,100)", cut.c_str());
  TTJets_miniTree->Draw("weight_Xsection >> temp_weight(100,0,10)", cut.c_str());
	integral_bg += TTJets_miniTree->GetEntries(cut.c_str()) * (((TH1F*)gDirectory->Get("temp_weight"))->GetMean()) * (((TH1F*)gDirectory->Get("temp_pileup"))->GetMean());
/*	WJetsToLNu_miniTree->Draw("weight_pileUp >> temp_pileup(100,0,100)", cut.c_str());
  WJetsToLNu_miniTree->Draw("weight_Xsection >> temp_weight(100,0,10)", cut.c_str());
	integral_bg += WJetsToLNu_miniTree->GetEntries(cut.c_str()) * (((TH1F*)gDirectory->Get("temp_weight"))->GetMean()) * (((TH1F*)gDirectory->Get("temp_pileup"))->GetMean());
	QCDMu_miniTree->Draw("weight_pileUp >> temp_pileup(100,0,100)", cut.c_str());
  QCDMu_miniTree->Draw("weight_Xsection >> temp_weight(100,0,10)", cut.c_str());
	integral_bg += QCDMu_miniTree->GetEntries(cut.c_str()) * (((TH1F*)gDirectory->Get("temp_weight"))->GetMean()) * (((TH1F*)gDirectory->Get("temp_pileup"))->GetMean());
*/	double integral_mc = integral_bg + integral_signal; 
	c1->Clear();

//	cout << "data= " << integral_data << endl;
//	cout << "mc= " << integral_mc << " = " << integral_signal << " (signal) + " << integral_bg << " (bg)" << endl;
//	double integral_data = Data_miniTree->GetEntries(cut.c_str());
//	cout << "Data_miniTree->GetEntries(cut.c_str())= " << Data_miniTree->GetEntries(cut.c_str()) << "\t\tHisto_Data->Integral()= " << Histo_Data->Integral(0, Histo_Data->GetNbinsX()+1) << endl;
//	double integral_signal = FSR_DYToMuMu_miniTree->GetEntries(cut.c_str());
//	cout << setprecision(5) << fixed;
//	cout << "FSR_DYToMuMu_miniTree->GetEntries(cut.c_str())=" << FSR_DYToMuMu_miniTree->GetEntries(cut.c_str()) << endl;
//	cout << "Histo_FSR_DYToMuMu->Integral()= " << Histo_FSR_DYToMuMu->Integral() << endl;
//	cout << "Histo_FSR_DYToMuMu->Integral(0, Histo_FSR_DYToMuMu->GetNbinsX()+1)= " << Histo_FSR_DYToMuMu->Integral(0, Histo_FSR_DYToMuMu->GetNbinsX()+1) << endl;
//	cout << "FSR_DYToMuMu_miniTree->GetEntries(cut_FSR_DYToMuMu)= " << FSR_DYToMuMu_miniTree->GetEntries(cut_FSR_DYToMuMu.c_str()) << endl;
//	 * (double)(  (double)((double)(XSectionFSR_DYToMuMu) / (double)(InitialNumberFSR_DYToMuMu)) * (double)integratedLuminosity) << "\t\tHisto_FSR_DYToMuMu->Integral()= " << Histo_FSR_DYToMuMu->Integral(0, Histo_FSR_DYToMuMu->GetNbinsX()+1) << "\t\tFSR_DYToMuMu_miniTree->GetEntries(cut_FSR_DYToMuMu)= " << FSR_DYToMuMu_miniTree->GetEntries(cut_FSR_DYToMuMu.c_str()) << endl;
//	FSR_DYToMuMu_miniTree->Draw("weight_pileUp >> temp_pileup(100,0,100)", cut.c_str());
//	cout << "test0= " << (((TH1F*)gDirectory->Get("temp_pileup"))->Integral()) << endl;
//	cout << "test1= " << (((TH1F*)gDirectory->Get("temp_pileup"))->GetMean()) << endl;
//	FSR_DYToMuMu_miniTree->Draw("weight_Xsection >> temp_weight(100,0,10)", cut.c_str());
//	cout << "test2= " << (((TH1F*)gDirectory->Get("temp_weight"))->Integral()) << endl;
//	cout << "test3= " << (((TH1F*)gDirectory->Get("temp_weight"))->GetMean()) << endl;
//	cout << "test4= " << FSR_DYToMuMu_miniTree->GetEntries(cut.c_str()) * (((TH1F*)gDirectory->Get("temp_weight"))->GetMean()) * (((TH1F*)gDirectory->Get("temp_pileup"))->GetMean()) << endl;
//	double integral_bg = nonFSR_DYToMuMu_miniTree->GetEntries(cut.c_str()) + TTJets_miniTree->GetEntries(cut.c_str()) + WJetsToLNu_miniTree->GetEntries(cut.c_str()) + QCDMu_miniTree->GetEntries(cut.c_str());
//	double integral_mc = FSR_DYToMuMu_miniTree->GetEntries(cut.c_str()) + nonFSR_DYToMuMu_miniTree->GetEntries(cut.c_str()) + TTJets_miniTree->GetEntries(cut.c_str()) + WJetsToLNu_miniTree->GetEntries(cut.c_str()) + QCDMu_miniTree->GetEntries(cut.c_str());
/*
	Histo_WJetsToLNu->Add(Histo_QCDMu);
	Histo_TTJets->Add(Histo_WJetsToLNu);
*/	Histo_nonFSR_DYToMuMu->Add(Histo_TTJets);
	Histo_FSR_DYToMuMu->Add(Histo_nonFSR_DYToMuMu);


//	double integral_data = Histo_Data->Integral();
//	double integral_signal = (double)(Histo_FSR_DYToMuMu->Integral()) - (double)(Histo_nonFSR_DYToMuMu->Integral());
//	double integral_bg = (double)(Histo_nonFSR_DYToMuMu->Integral());
//	double integral_mc = (double)(Histo_FSR_DYToMuMu->Integral());

//	double integral_data = Histo_Data->Integral(0, Histo_Data->GetNbinsX()+1);
//	double integral_signal = (double)(Histo_FSR_DYToMuMu->Integral(0, Histo_FSR_DYToMuMu->GetNbinsX()+1)) - (double)(Histo_nonFSR_DYToMuMu->Integral(0, Histo_nonFSR_DYToMuMu->GetNbinsX()+1));
//	double integral_bg = (double)(Histo_nonFSR_DYToMuMu->Integral(0, Histo_nonFSR_DYToMuMu->GetNbinsX()+1));
//	double integral_mc = (double)(Histo_FSR_DYToMuMu->Integral(0, Histo_FSR_DYToMuMu->GetNbinsX()+1));

	// Total MC histo for comupting min/max
//	TH1F *Histo_allMC = new TH1F(*Histo_QCD_Mu_Pt20to30);
//	Histo_allMC->Add(Histo_QCD_Pt15);
//	Histo_allMC->Add(Histo_InclusiveMu15);
//	Histo_allMC->Add(Histo_FSR_ZmumuJet_Pt0to15);
//	Histo_allMC->Add(Histo_ZJets_7TeV);
//	Histo_allMC->Add(Histo_WJets_7TeV);
//	Histo_allMC->Add(Histo_TTbarJets_Tauola);
//	Histo_allMC->Add(Histo_FSR_DYToMuMu);


  // Get the maxs and the mins to further correct the Y-axis
  double dataMax = Histo_Data->GetMaximum();
  double YMax = dataMax;

  double FSR_DYToMuMuMax = Histo_FSR_DYToMuMu->GetMaximum();
  YMax = max(YMax, FSR_DYToMuMuMax);

//	double allMCMax = Histo_allMC->GetMaximum();
//	YMax = max(YMax, allMCMax);

  double dataMin = YMax;
  double nonFSR_DYToMuMuMin = YMax;
  double FSR_DYToMuMuMin = YMax;
  double TTJetsMin = YMax;
//  double WJetsToLNuMin = YMax;
//  double QCDMuMin = YMax;

	double allMCMin = YMax;

  double YMin = YMax;

  // Gets the actual minimum for each histogram, and not the unfilled bin if any

  for( int ibin=1 ; ibin<Histo_Data->GetNbinsX() ; ibin++ ){
		if( ((Histo_Data->GetBinContent(ibin))!=0) ){
			YMax = max(YMax, (Histo_Data->GetBinContent(ibin) + Histo_Data->GetBinError(ibin)));
//			cout << "YMax= " << YMax << endl;
		}
//		cout << "ibin= " << ibin << "\tcontent= " << Histo_Data->GetBinContent(ibin) << "\terror= " << Histo_Data->GetBinError(ibin) << endl;
    if( ((Histo_Data->GetBinContent(ibin))!=0) && ((Histo_Data->GetBinContent(ibin))<dataMin) ){
      dataMin = Histo_Data->GetBinContent(ibin);
    }
  }
  YMin = min(YMin, dataMin);
  for( int ibin=1 ; ibin<Histo_FSR_DYToMuMu->GetNbinsX() ; ibin++ ){
    if( ((Histo_FSR_DYToMuMu->GetBinContent(ibin))!=0) && ((Histo_FSR_DYToMuMu->GetBinContent(ibin))<FSR_DYToMuMuMin) ){
      FSR_DYToMuMuMin = Histo_FSR_DYToMuMu->GetBinContent(ibin);
    }
  }
  YMin = min(YMin, FSR_DYToMuMuMin);

 for( int ibin=1 ; ibin<Histo_nonFSR_DYToMuMu->GetNbinsX() ; ibin++ ){
    if( ((Histo_nonFSR_DYToMuMu->GetBinContent(ibin))!=0) && ((Histo_nonFSR_DYToMuMu->GetBinContent(ibin))<nonFSR_DYToMuMuMin) ){
      nonFSR_DYToMuMuMin = Histo_nonFSR_DYToMuMu->GetBinContent(ibin);
    }
  }
  YMin = min(YMin, nonFSR_DYToMuMuMin);

  for( int ibin=1 ; ibin<Histo_TTJets->GetNbinsX() ; ibin++ ){
    if( ((Histo_TTJets->GetBinContent(ibin))!=0) && ((Histo_TTJets->GetBinContent(ibin))<TTJetsMin) ){
      TTJetsMin = Histo_TTJets->GetBinContent(ibin);
    }
  }
  YMin = min(YMin, TTJetsMin);
/*
  for( int ibin=1 ; ibin<Histo_WJetsToLNu->GetNbinsX() ; ibin++ ){
    if( ((Histo_WJetsToLNu->GetBinContent(ibin))!=0) && ((Histo_WJetsToLNu->GetBinContent(ibin))<WJetsToLNuMin) ){
      WJetsToLNuMin = Histo_WJetsToLNu->GetBinContent(ibin);
    }
  }
  YMin = min(YMin, WJetsToLNuMin);

  for( int ibin=1 ; ibin<Histo_QCDMu->GetNbinsX() ; ibin++ ){
    if( ((Histo_QCDMu->GetBinContent(ibin))!=0) && ((Histo_QCDMu->GetBinContent(ibin))<QCDMuMin) ){
      QCDMuMin = Histo_QCDMu->GetBinContent(ibin);
    }
  }
  YMin = min(YMin, QCDMuMin);
*/

//  cout << "YMax= "<< YMax << "\t\tYMin= " << YMin << endl;
  double YMin_lin = (double)YMin / (double)10.0;
//  double Range_lin = ((double)(YMax - YMin_lin)) / ((double)(0.8));
  double Range_lin = ((double)(YMax - YMin_lin)) / ((double)(1.0));
  double YMax_lin = 0.2*Range_lin + YMax;
/*
  double Range_lin = ((double)(YMax - YMin)) / ((double)(0.77));
  double YMax_lin = 0.2*Range_lin + YMax;
  double YMin_lin = max(YMin - 0.03*Range_lin, (double)YMin / (double)10.0);
*/
  double Range_log = ((double)(log10(YMax) - log10(YMin))) / ((double)(0.77));
//  cout << "Range_lin= " << Range_lin << "\t\tRange_log= " << Range_log << endl;
  double YMax_log = pow(10.0, 0.2*Range_log + log10(YMax));
  double YMin_log = pow(10.0, log10(YMin) - 0.03*Range_log);
//  cout << "YMin_lin= " << YMin_lin << "\t\tYMax_lin= " << YMax_lin << endl;
//  cout << "YMin_log= " << YMin_log << "\t\tYMax_log= " << YMax_log << endl;


  // Setup the histo and canvas names and title
  string data_name = "Data_" + pic + "_" + name;
  string mc_name = "MC_" + pic + "_" + name;
  string canvas_name = "DataMC_" + pic + "_" + name;
  std::ostringstream binWidthOSS;
  binWidthOSS << (double)Histo_Data->GetBinWidth(1);
  string binWidth = binWidthOSS.str();
  string YaxisTitle = "";
  if((Title.rfind("[") < Title.size()) && (Title.rfind("]") < Title.size())){
//    string unit = Title.substr(Title.rfind("[")+1, Title.size()-Title.rfind("]")-2);
    string unit = Title.substr(Title.rfind("[")+1, Title.rfind("]")-Title.rfind("[")-1);
    YaxisTitle = "Events / " + binWidth + " " + unit;
  } else {
    YaxisTitle = "Events / " + binWidth;
  }
  Histo_Data->SetName(data_name.c_str());
//  Histo_QCDMu->SetName(mc_name.c_str());
	Histo_FSR_DYToMuMu->SetName(mc_name.c_str());
	
  c1->SetName(canvas_name.c_str());
  c1->SetTitle(canvas_name.c_str());

  // Draw the comparison plots
//	// Template empty histo
//	TH1F *Histo_template = new TH1F("Histo_template", "Histo_template", Histo_Data->GetNbinsX(), Histo_Data->GetXaxis()->GetXmin(),  Histo_Data->GetXaxis()->GetXmax());
//	Histo_template->SetAxisRange(Histo_Data->GetXaxis()->GetXmin(),  Histo_Data->GetXaxis()->GetXmax(), "X");
//	Histo_template->SetAxisRange(YMin_lin, YMax_lin,"Y");
//	Histo_template->SetMaximum(YMax_lin);
//	Histo_template->SetMinimum(YMin_lin);
//	Histo_template->Draw();
//	c1->Update();
//	c1->Draw();

  TLegend *legend = new TLegend(0.65, 0.82, 0.90, 0.93, "");
	if( doFit ){
/*
		legend->SetX1(0.47);
		legend->SetX2(0.65);
*/
		legend->SetX1(0.15);
		legend->SetX2(0.40);
		legend->SetY1(0.65);
		legend->SetY2(0.78);

	}
	legend->SetTextSize(0.025);
  legend->SetFillColor(kWhite);
  legend->SetLineColor(kWhite);
  legend->SetShadowColor(kWhite);

  // // First: draw the data to get correct Y-axis scale
	gPad->Update();
  Histo_Data->GetXaxis()->SetTitle(Title.c_str());
  Histo_Data->GetYaxis()->SetTitle(YaxisTitle.c_str());
  Histo_Data->SetLineColor(kBlack);
  Histo_Data->SetMarkerColor(kBlack);
  Histo_Data->SetMarkerSize(0.7);
  Histo_Data->SetMarkerStyle(20);
  Histo_Data->SetMaximum(YMax_lin);
  Histo_Data->SetMinimum(YMin_lin);
//  Histo_Data->GetYaxis()->SetRangeUser(YMin_lin, YMax_lin);
//  Histo_Data->Draw("E1sames");
  Histo_Data->Draw("E1");
  legend->AddEntry(Histo_Data->GetName(), "Data", "lp");

  // // Second: draw MC on the same canvas
//  Histo_InclusiveMu15->SetLineColor(kBlack);
//  Histo_InclusiveMu15->SetFillColor(kGreen-6);
//  Histo_InclusiveMu15->SetFillStyle(3001);
//  Histo_InclusiveMu15->SetMaximum(YMax_lin);
//  Histo_InclusiveMu15->SetMinimum(YMin_lin);
//  Histo_InclusiveMu15->Draw("same");  
//  TTJetsAddEntry(Histo_InclusiveMu15->GetName(), "InclusiveMu15", "f");
/*
  Histo_QCDMu->SetLineColor(kBlack);
  Histo_QCDMu->SetFillColor(kGreen-6);
  Histo_QCDMu->SetFillStyle(3001);
  Histo_QCDMu->SetMaximum(YMax_lin);
  Histo_QCDMu->SetMinimum(YMin_lin);

  Histo_WJetsToLNu->SetLineColor(kBlack);
  Histo_WJetsToLNu->SetFillColor(kMagenta+3);
  Histo_WJetsToLNu->SetFillStyle(3001);
  Histo_WJetsToLNu->SetMaximum(YMax_lin);
  Histo_WJetsToLNu->SetMinimum(YMin_lin);
*/
  Histo_TTJets->SetLineColor(kBlack);
  Histo_TTJets->SetFillColor(kBlue);
  Histo_TTJets->SetFillStyle(3001);
  Histo_TTJets->SetMaximum(YMax_lin);
  Histo_TTJets->SetMinimum(YMin_lin);

//  Histo_FSR_DYToMuMu->SetLineColor(kBlack);
//  Histo_FSR_DYToMuMu->SetFillColor(kMagenta);
//  Histo_FSR_DYToMuMu->SetFillStyle(3001);
//  Histo_FSR_DYToMuMu->SetMaximum(YMax_lin);
//  Histo_FSR_DYToMuMu->SetMinimum(YMin_lin);

  Histo_FSR_DYToMuMu->SetLineColor(kBlack);
  Histo_FSR_DYToMuMu->SetFillColor(kRed);
  Histo_FSR_DYToMuMu->SetFillStyle(3001);
  Histo_FSR_DYToMuMu->SetMaximum(YMax_lin);
  Histo_FSR_DYToMuMu->SetMinimum(YMin_lin);
  Histo_FSR_DYToMuMu->GetXaxis()->SetTitle(Title.c_str());
  Histo_FSR_DYToMuMu->GetYaxis()->SetTitle(YaxisTitle.c_str());

	Histo_nonFSR_DYToMuMu->SetLineColor(kBlack);
  Histo_nonFSR_DYToMuMu->SetFillColor(kOrange);
  Histo_nonFSR_DYToMuMu->SetFillStyle(3001);
  Histo_nonFSR_DYToMuMu->SetMaximum(YMax_lin);
  Histo_nonFSR_DYToMuMu->SetMinimum(YMin_lin);

//	Histo_QCD_Mu_Pt20to30->Draw("same");
  Histo_FSR_DYToMuMu->Draw("same");
  Histo_nonFSR_DYToMuMu->Draw("same");
  Histo_TTJets->Draw("same");
//  Histo_WJetsToLNu->Draw("same");
//  Histo_QCDMu->Draw("same");
  legend->AddEntry(Histo_FSR_DYToMuMu->GetName(), "Z#mu#muJet FSR", "f");
  legend->AddEntry(Histo_nonFSR_DYToMuMu->GetName(), "Z#mu#muJets (no FSR)", "f");
  legend->AddEntry(Histo_TTJets->GetName(), "t#bar{t}Jets", "f");
//  legend->AddEntry(Histo_WJetsToLNu->GetName(), "WJets", "f");
//  legend->AddEntry(Histo_QCDMu->GetName(), "QCD #mu", "f");
//  legend->AddEntry(Histo_FSR_DYToMuMu->GetName(), "PhotonJet", "f");
//  legend->AddEntry(Histo_QCD_Mu_Pt20to30->GetName(), "QCD Mu", "f");
//  legend->AddEntry(Histo_ZJets_7TeV->GetName(), "ZJets madgraph", "f");

  // // Third: re-draw Data so that data appears in front of MC
  Histo_Data->Draw("E1same");

/*
if( doFit ){
	Histo_Data->SetName("DATA");
	Histo_FSR_ZmumuJet_Pt0to15->SetName("MC");
	gStyle->SetOptStat(2201);
	string fitOpts = "MQ+";

	TF1* f_Data = new TF1("f_Data", "gaus");
	TF1* f_MC = new TF1("f_MC", "gaus");
	f_Data->SetLineColor(kBlue);
	f_MC->SetLineColor(kGreen+2);
	f_Data->SetLineWidth(2);
	f_MC->SetLineWidth(2);
	Histo_Data->Fit("f_Data", fitOpts.c_str());
	Histo_FSR_ZmumuJet_Pt0to15->Fit("f_MC", fitOpts.c_str());

	Histo_FSR_ZmumuJet_Pt0to15->SetMaximum(YMax_lin);
	Histo_FSR_ZmumuJet_Pt0to15->SetMinimum(YMin_lin);	
	Histo_QCD_Mu_Pt20to30->SetMaximum(YMax_lin);
	Histo_QCD_Mu_Pt20to30->SetMinimum(YMin_lin);	
	Histo_WJets_7TeV->SetMaximum(YMax_lin);
	Histo_WJets_7TeV->SetMinimum(YMin_lin);	
	Histo_FSR_DYToMuMu->SetMaximum(YMax_lin);
	Histo_FSR_DYToMuMu->SetMinimum(YMin_lin);	
	Histo_nonFSR_ZmumuJet_Pt0to15->SetMaximum(YMax_lin);
	Histo_nonFSR_ZmumuJet_Pt0to15->SetMinimum(YMin_lin);	
	Histo_TTbarJets_Tauola->SetMaximum(YMax_lin);
	Histo_TTbarJets_Tauola->SetMinimum(YMin_lin);	
	Histo_Data->SetMaximum(YMax_lin);
	Histo_Data->SetMinimum(YMin_lin);	
  Histo_FSR_ZmumuJet_Pt0to15->Draw("sames");
	//  Histo_ZJets_7TeV->Draw("sames");
  Histo_QCD_Mu_Pt20to30->Draw("sames");
  Histo_WJets_7TeV->Draw("sames");
  Histo_FSR_DYToMuMu->Draw("sames");
  Histo_nonFSR_ZmumuJet_Pt0to15->Draw("sames");
  Histo_TTbarJets_Tauola->Draw("sames");
  Histo_Data->GetXaxis()->SetTitle(Title.c_str());
  Histo_Data->GetYaxis()->SetTitle(YaxisTitle.c_str());
  Histo_Data->Draw("E1sames");
//	Histo_Data->SetMaximum(YMax_lin);
//	Histo_Data->SetMinimum(YMin_lin);	
//	c1->Update();	


//	TPaveStats* statsHisto_MC = (TPaveStats*) Histo_ZJets_7TeV->GetListOfFunctions()->FindObject("stats");
//  statsHisto_ZJets_7TeV->SetLineColor(kWhite);
//  statsHisto_ZJets_7TeV->SetTextColor(kWhite);
//	statsHisto_ZJets_7TeV->Draw();


	gPad->Update();
 	TPaveStats* statsHisto_QCD_Mu_Pt20to30 = (TPaveStats*) Histo_QCD_Mu_Pt20to30->GetListOfFunctions()->FindObject("stats");
  statsHisto_QCD_Mu_Pt20to30->SetLineColor(kWhite);
  statsHisto_QCD_Mu_Pt20to30->SetTextColor(kWhite);
	statsHisto_QCD_Mu_Pt20to30->Draw();
 
	gPad->Update();
 	TPaveStats* statsHisto_WJets_7TeV = (TPaveStats*) Histo_WJets_7TeV->GetListOfFunctions()->FindObject("stats");
  statsHisto_WJets_7TeV->SetLineColor(kWhite);
  statsHisto_WJets_7TeV->SetTextColor(kWhite);
	statsHisto_WJets_7TeV->Draw();
  
	gPad->Update();
	TPaveStats* statsHisto_FSR_DYToMuMu = (TPaveStats*) Histo_FSR_DYToMuMu->GetListOfFunctions()->FindObject("stats");
  statsHisto_FSR_DYToMuMu->SetLineColor(kWhite);
  statsHisto_FSR_DYToMuMu->SetTextColor(kWhite);
	statsHisto_FSR_DYToMuMu->Draw();
  
	gPad->Update();
	TPaveStats* statsHisto_nonFSR_ZmumuJet_Pt0to15 = (TPaveStats*) Histo_nonFSR_ZmumuJet_Pt0to15->GetListOfFunctions()->FindObject("stats");
  statsHisto_nonFSR_ZmumuJet_Pt0to15->SetLineColor(kWhite);
  statsHisto_nonFSR_ZmumuJet_Pt0to15->SetTextColor(kWhite);
	statsHisto_nonFSR_ZmumuJet_Pt0to15->Draw();

	gPad->Update();
	TPaveStats* statsHisto_TTbarJets_Tauola = (TPaveStats*) Histo_TTbarJets_Tauola->GetListOfFunctions()->FindObject("stats");
  statsHisto_TTbarJets_Tauola->SetLineColor(kWhite);
  statsHisto_TTbarJets_Tauola->SetTextColor(kWhite);
	statsHisto_TTbarJets_Tauola->Draw();
 

	gPad->Update();
	TPaveStats* statsHisto_Data = (TPaveStats*) Histo_Data->GetListOfFunctions()->FindObject("stats");
	statsHisto_Data->SetLineColor(kBlue);
	statsHisto_Data->SetTextColor(kBlue);
	statsHisto_Data->SetTextSize(0.020);
	statsHisto_Data->SetX1NDC(0.68);
	statsHisto_Data->SetX2NDC(0.93);
	statsHisto_Data->SetY1NDC(0.80);
	statsHisto_Data->SetY2NDC(0.93);
	statsHisto_Data->Draw();


	gPad->Update();
	TPaveStats* statsHisto_FSR_ZmumuJet_Pt0to15 = (TPaveStats*) Histo_FSR_ZmumuJet_Pt0to15->GetListOfFunctions()->FindObject("stats");
	statsHisto_FSR_ZmumuJet_Pt0to15->SetLineColor(kGreen+2);
	statsHisto_FSR_ZmumuJet_Pt0to15->SetTextColor(kGreen+2);
	statsHisto_FSR_ZmumuJet_Pt0to15->SetTextSize(0.020);
	statsHisto_FSR_ZmumuJet_Pt0to15->SetX1NDC(0.68);
	statsHisto_FSR_ZmumuJet_Pt0to15->SetX2NDC(0.93);
	statsHisto_FSR_ZmumuJet_Pt0to15->SetY1NDC(0.67);
	statsHisto_FSR_ZmumuJet_Pt0to15->SetY2NDC(0.80);
	statsHisto_FSR_ZmumuJet_Pt0to15->Draw();

//	c1->Update();	

}
*/
  // // Fourth: redraw axis so that axis appears in front of everything
  gPad->RedrawAxis();

  // // Fifth: draw legend
  legend->Draw();
	c1->Update();

  TLatex latexLabel;
  std::ostringstream intLumiString;
  intLumiString << setprecision (2) << fixed << (double)integratedLuminosity/(double)1000.0;
  string intLumiText = "#intL= " + intLumiString.str() + " fb^{-1}";
  latexLabel.SetTextSize(0.03);
  latexLabel.SetNDC();
  latexLabel.DrawLatex(0.13, 0.96, "CMS Preliminary 2011");
  latexLabel.DrawLatex(0.42, 0.96, "#sqrt{s} = 7 TeV");
  latexLabel.DrawLatex(0.57, 0.96, intLumiText.c_str());

	TLatex latexYields;
	latexYields.SetTextSize(0.03);
	latexYields.SetNDC();
	std::ostringstream tempString;
	tempString << setprecision (0) << fixed << integral_data;
	string tempText = "N_{data}= " + tempString.str();
	latexYields.DrawLatex(0.18, 0.90, tempText.c_str());
	std::ostringstream tempString2;
	tempString2 << setprecision (2) << fixed << integral_mc;
	tempText = "N_{MC}= " + tempString2.str();
	latexYields.DrawLatex(0.18, 0.86, tempText.c_str());
	std::ostringstream tempString3;
  tempString3 << setprecision (1) << fixed << 100.0*(double)(integral_signal)/(double)(integral_mc);
	tempText = "purity= " + tempString3.str() + "\%";
	latexYields.DrawLatex(0.16, 0.82, tempText.c_str());


  // // Sixth: update canvas
  c1->Update();
  c1->Draw();

  // Print the canvas
  string PicName="gif/DataMC_" + pic + "_" + name + ".gif";
  c1->Print(PicName.c_str());
  PicName="eps/DataMC_" + pic + "_" + name + ".eps";
  c1->Print(PicName.c_str());
  PicName="C/DataMC_" + pic + "_" + name + ".C";
  c1->Print(PicName.c_str());
  PicName="png/DataMC_" + pic + "_" + name + ".png";
  c1->Print(PicName.c_str());
  string convert = "convert eps/DataMC_" + pic + "_" + name + ".eps" + " pdf/DataMC_" + pic + "_" + name + ".pdf";
  system(convert.c_str());
  if (inlog==true) {
    c1->cd(1);
    Histo_Data->SetMaximum(YMax_log);
    Histo_Data->SetMinimum(YMin_log);
    Histo_Data->GetYaxis()->SetRangeUser(YMin_log, YMax_log);
    Histo_FSR_DYToMuMu->SetMaximum(YMax_log);
    Histo_FSR_DYToMuMu->SetMinimum(YMin_log);
    Histo_FSR_DYToMuMu->GetYaxis()->SetRangeUser(YMin_log, YMax_log);
//    Histo_QCD_Pt15->SetMaximum(YMax_log);
//    Histo_QCD_Pt15->SetMinimum(YMin_log);
//    Histo_QCD_Pt15->GetYaxis()->SetRangeUser(YMin_log, YMax_log);

    Histo_TTJets->SetMaximum(YMax_log);
    Histo_TTJets->SetMinimum(YMin_log);
    Histo_TTJets->GetYaxis()->SetRangeUser(YMin_log, YMax_log);
/*
    Histo_WJetsToLNu->SetMaximum(YMax_log);
    Histo_WJetsToLNu->SetMinimum(YMin_log);
    Histo_WJetsToLNu->GetYaxis()->SetRangeUser(YMin_log, YMax_log);

    Histo_QCDMu->SetMaximum(YMax_log);
    Histo_QCDMu->SetMinimum(YMin_log);
    Histo_QCDMu->GetYaxis()->SetRangeUser(YMin_log, YMax_log);
*/
//    Histo_InclusiveMu15->SetMaximum(YMax_log);
//    Histo_InclusiveMu15->SetMinimum(YMin_log);
//    Histo_InclusiveMu15->GetYaxis()->SetRangeUser(YMin_log, YMax_log);

    Histo_nonFSR_DYToMuMu->SetMaximum(YMax_log);
    Histo_nonFSR_DYToMuMu->SetMinimum(YMin_log);
    Histo_nonFSR_DYToMuMu->GetYaxis()->SetRangeUser(YMin_log, YMax_log);

    c1->SetLogy(1);
    c1->Update();
    c1->Draw();
    string PicName_log="gif/DataMC_" + pic + "_" + name + "_log.gif";
    c1->Print(PicName_log.c_str());
    PicName="eps/DataMC_" + pic + "_" + name + "_log.eps";
    c1->Print(PicName.c_str());
    PicName="C/DataMC_" + pic + "_" + name + "_log.C";
    c1->Print(PicName.c_str());
    PicName="png/DataMC_" + pic + "_" + name + "_log.png";
    c1->Print(PicName.c_str());
    string convert = "convert eps/DataMC_" + pic + "_" + name + "_log.eps" + " pdf/DataMC_" + pic + "_" + name + "_log.pdf";
    system(convert.c_str());
    c1->SetLogy(0);
    c1->Update();
  }

  // Clean the memory
  c1->Clear();
  legend->Clear();
//	Histo_template->Delete();
  Histo_Data_temp->Delete();
  Histo_Data->Delete();
  Histo_FSR_DYToMuMu_temp->Delete();
  Histo_FSR_DYToMuMu->Delete();

  Histo_nonFSR_DYToMuMu_temp->Delete();
  Histo_nonFSR_DYToMuMu->Delete();

  Histo_TTJets_temp->Delete();
  Histo_TTJets->Delete();
/*
  Histo_WJetsToLNu_temp->Delete();
  Histo_WJetsToLNu->Delete();

  Histo_QCDMu_temp->Delete();
  Histo_QCDMu->Delete();
*/
}


