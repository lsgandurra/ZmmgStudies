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

#include "DrawDataData.h"

#include "interface/TRootBardak.h"
#include "interface/TRootBeamSpot.h"
#include "interface/TRootCluster.h"
#include "interface/TRootEcalRecHit.h"
#include "interface/TRootElectron.h"
#include "interface/TRootEvent.h"
#include "interface/TRootJet.h"
#include "interface/TRootMCParticle.h"
#include "interface/TRootMCPhoton.h"
#include "interface/TRootMET.h"
#include "interface/TRootMuon.h"
#include "interface/TRootParticle.h"
#include "interface/TRootPhoton.h"
#include "interface/TRootRun.h"
#include "interface/TRootSignalEvent.h"
#include "interface/TRootSuperCluster.h"
#include "interface/TRootTopTop.h"
#include "interface/TRootTrack.h"
#include "interface/TRootVertex.h"

#include "CMSStyle.C"

using namespace std;

double fonction_affine(double *x, double *par){
	return x[0]*par[0] - x[1]*par[1];
}

void DrawDataDataplot(TTree *Data_miniTree, TTree *OLD_Data_miniTree, string var, string pic, string limits, string cut, string name, string Title, bool inlog, bool drawUnderOverFsubleading, TCanvas *c1){

	CMSstyle();
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);

  // Get Histo_Data from eventTree
  TH1F *Histo_Data_temp = new TH1F();
  string variable_Data = var + ">>Histo_Data_temp" + limits;
  Data_miniTree->Draw(variable_Data.c_str(), cut.c_str());
  TH1F *Histo_Data = (TH1F*)gDirectory->Get("Histo_Data_temp");
  c1->Clear();
/*
  // Get Histo_DYToMuMu from eventTree
  TH1F *Histo_DYToMuMu_temp = new TH1F();
  string variable_DYToMuMu = var + ">>Histo_DYToMuMu_temp" + limits;
  DYToMuMu_miniTree->Draw(variable_DYToMuMu.c_str(), cut.c_str());
  TH1F *Histo_DYToMuMu = (TH1F*)gDirectory->Get("Histo_DYToMuMu_temp");
  c1->Clear();


  // Get Histo_OLD_DYToMuMu from eventTree
  TH1F *Histo_OLD_DYToMuMu_temp = new TH1F();
  string variable_OLD_DYToMuMu = var + ">>Histo_OLD_DYToMuMu_temp" + limits;
  OLD_DYToMuMu_miniTree->Draw(variable_OLD_DYToMuMu.c_str(), cut.c_str());
  TH1F *Histo_OLD_DYToMuMu = (TH1F*)gDirectory->Get("Histo_OLD_DYToMuMu_temp");
  c1->Clear();
*/

  // Get Histo_OLD_Data from eventTree
  TH1F *Histo_OLD_Data_temp = new TH1F();
  string variable_OLD_Data = var + ">>Histo_OLD_Data_temp" + limits;
  OLD_Data_miniTree->Draw(variable_OLD_Data.c_str(), cut.c_str());
  TH1F *Histo_OLD_Data = (TH1F*)gDirectory->Get("Histo_OLD_Data_temp");
  c1->Clear();

  // Get Histo_WJetsToLNu from eventTree
//  TH1F *Histo_WJetsToLNu_temp = new TH1F();
//  string variable_WJetsToLNu = var + ">>Histo_WJetsToLNu_temp" + limits;
//  WJetsToLNu_miniTree->Draw(variable_WJetsToLNu.c_str(), cut.c_str());
//  TH1F *Histo_WJetsToLNu = (TH1F*)gDirectory->Get("Histo_WJetsToLNu_temp");
//  c1->Clear();

  // Get Histo_QCDMu from eventTree
//  TH1F *Histo_QCDMu_temp = new TH1F();
//  string variable_QCDMu = var + ">>Histo_QCDMu_temp" + limits;
//  QCDMu_miniTree->Draw(variable_QCDMu.c_str(), cut.c_str());
//  TH1F *Histo_QCDMu = (TH1F*)gDirectory->Get("Histo_QCDMu_temp");
//  c1->Clear();

  // Get the number of entries for further normalization
//  double a = Histo_Data->Integral();
/*
  double b_DYToMuMu = Histo_DYToMuMu->Integral();
  if( (a==0.0) || (b_DYToMuMu==0.0) ){
    cout << "no entries to plots" <<endl;
    return; 
  }*/
  // Normalize
  Histo_Data->Sumw2(); // In order to have the correct error bars on data after renormalization
  Histo_OLD_Data->Sumw2(); // In order to have the correct error bars on data after renormalization
  // // Normalize MC and Data to 1
  //Histo_Data->Scale((double)((double)1.0/(double)a));
  //Histo_MC->Scale((double)((double)1.0/(double)b));
  // // Normalize MC to Data number of entries
//  double integratedLuminosity = 191.09326;

//  double XSectionDYToMuMu = 1614.0;
//  double XSectionOLD_DYToMuMu = 1614.0;
//  double XSectionOLD_Data = 121.0;
//	double XSectionWJetsToLNu = 24640.0;
//	double XSectionQCDMu = 84679.30;

//  double InitialNumberDYToMuMu = 1995236.0;
//  double InitialNumberOLD_DYToMuMu = 1995236.0;
//  double InitialNumberOLD_Data = 1164208.0;
//	double InitialNumberWJetsToLNu = 15110974.0;
//	double InitialNumberQCDMu = 29434562.0;

// Normalize everything to 1
	double N_Data = Histo_Data->Integral();
	double N_OLD_Data = Histo_OLD_Data->Integral();
//	double N_DYToMuMu = Histo_DYToMuMu->Integral();
//	double N_OLD_DYToMuMu = Histo_OLD_DYToMuMu->Integral();

	Histo_Data->Scale((double)((double)1.0/(double)N_Data));
	Histo_OLD_Data->Scale((double)((double)1.0/(double)N_OLD_Data));
//	Histo_DYToMuMu->Scale((double)((double)1.0/(double)N_DYToMuMu));
//	Histo_OLD_DYToMuMu->Scale((double)((double)1.0/(double)N_OLD_DYToMuMu));


//  Histo_DYToMuMu->Scale((double)(  (double)((double)(XSectionDYToMuMu) / (double)(InitialNumberDYToMuMu)) * (double)integratedLuminosity));
//  Histo_OLD_DYToMuMu->Scale((double)(  (double)((double)(XSectionOLD_DYToMuMu) / (double)(InitialNumberOLD_DYToMuMu)) * (double)integratedLuminosity));
//  Histo_OLD_Data->Scale((double)(  (double)((double)(XSectionOLD_Data) / (double)(InitialNumberOLD_Data)) * (double)integratedLuminosity));
//  Histo_WJetsToLNu->Scale((double)(  (double)((double)(XSectionWJetsToLNu) / (double)(InitialNumberWJetsToLNu)) * (double)integratedLuminosity));
//  Histo_QCDMu->Scale((double)(  (double)((double)(XSectionQCDMu) / (double)(InitialNumberQCDMu)) * (double)integratedLuminosity));
  // Adding histograms for binned samples
//  Histo_QCD_Pt15->Add(Histo_QCD_Pt30);
//  Histo_QCD_Pt15->Add(Histo_QCD_Pt80);
//  Histo_QCD_Pt15->Add(Histo_QCD_Pt170);
//  Histo_QCD_Pt15->Add(Histo_QCD_Pt300);
//  Histo_QCD_Pt15->Add(Histo_QCD_Pt470);

//	Histo_WJetsToLNu->Add(Histo_QCDMu);
//	Histo_OLD_Data->Add(Histo_WJetsToLNu);
//	Histo_OLD_DYToMuMu->Add(Histo_OLD_Data);
//	Histo_DYToMuMu->Add(Histo_OLD_DYToMuMu);

	// Total MC histo for comupting min/max
//	TH1F *Histo_allMC = new TH1F(*Histo_QCD_Mu_Pt20to30);
//	Histo_allMC->Add(Histo_QCD_Pt15);
//	Histo_allMC->Add(Histo_InclusiveMu15);
//	Histo_allMC->Add(Histo_ZmumuJet_Pt0to15);
//	Histo_allMC->Add(Histo_ZJets_7TeV);
//	Histo_allMC->Add(Histo_WJets_7TeV);
//	Histo_allMC->Add(Histo_TTbarJets_Tauola);
//	Histo_allMC->Add(Histo_DYToMuMu);


  // Get the maxs and the mins to further correct the Y-axis
  double dataMax = Histo_Data->GetMaximum();
  double YMax = dataMax;

//  double DYToMuMuMax = Histo_DYToMuMu->GetMaximum();
//  YMax = max(YMax, DYToMuMuMax);

	double OLD_dataMax = Histo_OLD_Data->GetMaximum();
	YMax = max(YMax, OLD_dataMax);

//	double OLD_DYToMuMuMax = Histo_OLD_DYToMuMu->GetMaximum();
//  YMax = max(YMax, OLD_DYToMuMuMax);

//	double allMCMax = Histo_allMC->GetMaximum();
//	YMax = max(YMax, allMCMax);

  double dataMin = YMax;
//  double OLD_DYToMuMuMin = YMax;
//  double DYToMuMuMin = YMax;
  double OLD_DataMin = YMax;
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
/*
  for( int ibin=1 ; ibin<Histo_DYToMuMu->GetNbinsX() ; ibin++ ){
    if( ((Histo_DYToMuMu->GetBinContent(ibin))!=0) && ((Histo_DYToMuMu->GetBinContent(ibin))<DYToMuMuMin) ){
      DYToMuMuMin = Histo_DYToMuMu->GetBinContent(ibin);
    }
  }
  YMin = min(YMin, DYToMuMuMin);

 for( int ibin=1 ; ibin<Histo_OLD_DYToMuMu->GetNbinsX() ; ibin++ ){
    if( ((Histo_OLD_DYToMuMu->GetBinContent(ibin))!=0) && ((Histo_OLD_DYToMuMu->GetBinContent(ibin))<OLD_DYToMuMuMin) ){
      OLD_DYToMuMuMin = Histo_OLD_DYToMuMu->GetBinContent(ibin);
    }
  }
  YMin = min(YMin, OLD_DYToMuMuMin);
*/
  for( int ibin=1 ; ibin<Histo_OLD_Data->GetNbinsX() ; ibin++ ){
    if( ((Histo_OLD_Data->GetBinContent(ibin))!=0) && ((Histo_OLD_Data->GetBinContent(ibin))<OLD_DataMin) ){
      OLD_DataMin = Histo_OLD_Data->GetBinContent(ibin);
    }
  }
  YMin = min(YMin, OLD_DataMin);

//  for( int ibin=1 ; ibin<Histo_WJetsToLNu->GetNbinsX() ; ibin++ ){
//    if( ((Histo_WJetsToLNu->GetBinContent(ibin))!=0) && ((Histo_WJetsToLNu->GetBinContent(ibin))<WJetsToLNuMin) ){
//      WJetsToLNuMin = Histo_WJetsToLNu->GetBinContent(ibin);
//    }
//  }
//  YMin = min(YMin, WJetsToLNuMin);

//  for( int ibin=1 ; ibin<Histo_QCDMu->GetNbinsX() ; ibin++ ){
//    if( ((Histo_QCDMu->GetBinContent(ibin))!=0) && ((Histo_QCDMu->GetBinContent(ibin))<QCDMuMin) ){
//      QCDMuMin = Histo_QCDMu->GetBinContent(ibin);
//    }
//  }
//  YMin = min(YMin, QCDMuMin);


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



	c1->Divide(1,2);
	c1->cd(1);
//	gPad->SetNumber(1);
	gPad->SetPad(0,0.2,1,1);
//	gPad->SetBottomMargin(0);
	gPad->Draw();

	c1->cd(2);
//	gPad->SetNumber(2);
	gPad->SetPad(0,0.,1,0.2);
//	gPad->SetTopMargin(0);
//	gPad->SetBottomMargin(0.3);
	gPad->Draw();


	

/*
TPad *pad =new TPad("haut","haut",0,0.4,1,1);
	    pad->SetNumber(1);
//	    cout << pad->GetBottomMargin() << endl;
//	    pad->SetBottomMargin(0);
	    pad->Draw();
	    
	    TPad *pad2 =new TPad("milieu","milieu",0,0.2,1,0.4);
	    pad2->SetNumber(2);
//	    pad2->SetTopMargin(0);
//	   pad2->SetBottomMargin(0.3);
	    pad2->Draw();

		TPad *pad3=new TPad("bas", "bas", 0, 0, 1, 0.2);
		pad3->SetNumber(3);
		pad3->Draw();
*/

	c1->cd(1);



  // Setup the histo and canvas names and title
  string data_name = "Data_" + pic + "_" + name;
  string mc_name = "MC_" + pic + "_" + name;
  string canvas_name = "DataData_" + pic + "_" + name;
  std::ostringstream binWidthOSS;
  binWidthOSS << (double)Histo_Data->GetBinWidth(1);
  string binWidth = binWidthOSS.str();
  string YaxisTitle = "";
  if((Title.rfind("[") < Title.size()) && (Title.rfind("]") < Title.size())){
//    string unit = Title.substr(Title.rfind("[")+1, Title.size()-Title.rfind("]")-2);
    string unit = Title.substr(Title.rfind("[")+1, Title.rfind("]")-Title.rfind("[")-1);
    YaxisTitle = "a.u. / " + binWidth + " " + unit;
  } else {
    YaxisTitle = "a.u. / " + binWidth;
  }
  Histo_Data->SetName(data_name.c_str());
//  Histo_QCDMu->SetName(mc_name.c_str());
//	Histo_DYToMuMu->SetName(mc_name.c_str());
	
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

  // // Second: draw MC on the same canvas
//  Histo_InclusiveMu15->SetLineColor(kBlack);
//  Histo_InclusiveMu15->SetFillColor(kGreen-6);
//  Histo_InclusiveMu15->SetFillStyle(3001);
//  Histo_InclusiveMu15->SetMaximum(YMax_lin);
//  Histo_InclusiveMu15->SetMinimum(YMin_lin);
//  Histo_InclusiveMu15->Draw("same");  
//  OLD_DataAddEntry(Histo_InclusiveMu15->GetName(), "InclusiveMu15", "f");

//  Histo_QCDMu->SetLineColor(kBlack);
//  Histo_QCDMu->SetFillColor(kGreen-6);
//  Histo_QCDMu->SetFillStyle(3001);
//  Histo_QCDMu->SetMaximum(YMax_lin);
//  Histo_QCDMu->SetMinimum(YMin_lin);

//  Histo_WJetsToLNu->SetLineColor(kBlack);
//  Histo_WJetsToLNu->SetFillColor(kMagenta+3);
//  Histo_WJetsToLNu->SetFillStyle(3001);
//  Histo_WJetsToLNu->SetMaximum(YMax_lin);
//  Histo_WJetsToLNu->SetMinimum(YMin_lin);
  Histo_Data->Draw("E1");

	Histo_OLD_Data->SetLineColor(kRed);
	Histo_OLD_Data->SetMarkerColor(kRed);
	Histo_OLD_Data->SetMarkerSize(0.7);
	Histo_OLD_Data->SetMarkerStyle(20);
//  Histo_OLD_Data->SetFillColor(kBlue);
//  Histo_OLD_Data->SetFillStyle(3001);
  Histo_OLD_Data->SetMaximum(YMax_lin);
  Histo_OLD_Data->SetMinimum(YMin_lin);
  Histo_OLD_Data->Draw("E1same");

//  Histo_DYToMuMu->SetLineColor(kBlack);
//  Histo_DYToMuMu->SetFillColor(kMagenta);
//  Histo_DYToMuMu->SetFillStyle(3001);
//  Histo_DYToMuMu->SetMaximum(YMax_lin);
//  Histo_DYToMuMu->SetMinimum(YMin_lin);

/*
  Histo_DYToMuMu->SetLineColor(kBlack);
  Histo_DYToMuMu->SetFillColor(kRed);
  Histo_DYToMuMu->SetFillStyle(3001);
  Histo_DYToMuMu->SetMaximum(YMax_lin);
  Histo_DYToMuMu->SetMinimum(YMin_lin);
  Histo_DYToMuMu->GetXaxis()->SetTitle(Title.c_str());
  Histo_DYToMuMu->GetYaxis()->SetTitle(YaxisTitle.c_str());
  Histo_DYToMuMu->Draw("HISTsame");

	Histo_OLD_DYToMuMu->SetLineColor(kBlack);
  Histo_OLD_DYToMuMu->SetFillColor(kOrange);
  Histo_OLD_DYToMuMu->SetFillStyle(3001);
  Histo_OLD_DYToMuMu->SetMaximum(YMax_lin);
  Histo_OLD_DYToMuMu->SetMinimum(YMin_lin);
  Histo_OLD_DYToMuMu->Draw("HISTsame");
*/
//	Histo_QCD_Mu_Pt20to30->Draw("same");
//  Histo_WJetsToLNu->Draw("same");
//  Histo_QCDMu->Draw("same");
  legend->AddEntry(Histo_Data->GetName(), "Data Run2011A", "lp");
  legend->AddEntry(Histo_OLD_Data->GetName(), "Data Run2011B", "lp");
//  legend->AddEntry(Histo_DYToMuMu->GetName(), "Z#mu#muJet 41X", "f");
//  legend->AddEntry(Histo_OLD_DYToMuMu->GetName(), "Z#mu#muJets 39X", "f");
//  legend->AddEntry(Histo_WJetsToLNu->GetName(), "WJets", "f");
//  legend->AddEntry(Histo_QCDMu->GetName(), "QCD #mu", "f");
//  legend->AddEntry(Histo_DYToMuMu->GetName(), "PhotonJet", "f");
//  legend->AddEntry(Histo_QCD_Mu_Pt20to30->GetName(), "QCD Mu", "f");
//  legend->AddEntry(Histo_ZJets_7TeV->GetName(), "ZJets madgraph", "f");

  // // Third: re-draw Data so that data appears in front of MC
  Histo_Data->Draw("E1same");
  Histo_OLD_Data->Draw("E1same");

  // // Fourth: redraw axis so that axis appears in front of everything
  gPad->RedrawAxis();

  // // Fifth: draw legend
  legend->Draw();
	c1->Update();

	c1->cd(2);
	TH1F *Histo_DataRatio = (TH1F*) Histo_Data->Clone("Dataratio");
	Histo_DataRatio->Reset();
//	Histo_DataRatio->Sumw2();
	Histo_DataRatio->Divide(Histo_Data, Histo_OLD_Data, 1.0, 1.0);
	Histo_DataRatio->SetMinimum(0.0);
	Histo_DataRatio->SetMaximum(2);
	Histo_DataRatio->GetYaxis()->SetTitle("DATA / OLD DATA");
	Histo_DataRatio->SetTitleOffset(0.35,"Y");
	Histo_DataRatio->SetTitleSize(0.11,"Y");
	Histo_DataRatio->SetLabelSize(0.1,"Y"); 
	Histo_DataRatio->GetYaxis()->CenterTitle();
	Histo_DataRatio->SetNdivisions(509 ,"Y");
	Histo_DataRatio->Draw("EP");
	TLine *l = new TLine(Histo_DataRatio->GetXaxis()->GetXmin(),1.,Histo_DataRatio->GetXaxis()->GetXmax(),1.);
	l->SetLineWidth(1.5); 
	l->Draw("same");
	c1->Update();
/*
	c1->cd(3);
 	TH1F *Histo_DYToMuMuRatio = (TH1F*) Histo_DYToMuMu->Clone("DYToMuMuratio");
	Histo_DYToMuMuRatio->Reset();
//	Histo_DYToMuMuRatio->Sumw2();
	Histo_DYToMuMuRatio->Divide(Histo_DYToMuMu, Histo_OLD_DYToMuMu, 1.0, 1.0);
	Histo_DYToMuMuRatio->SetMinimum(0.0);
	Histo_DYToMuMuRatio->SetMaximum(2);
	Histo_DYToMuMuRatio->GetYaxis()->SetTitle("MC / OLD MC");
	Histo_DYToMuMuRatio->SetTitleOffset(0.35,"Y");
	Histo_DYToMuMuRatio->SetTitleSize(0.11,"Y");
	Histo_DYToMuMuRatio->SetLabelSize(0.1,"Y"); 
	Histo_DYToMuMuRatio->GetYaxis()->CenterTitle();
	Histo_DYToMuMuRatio->SetNdivisions(509 ,"Y");
	Histo_DYToMuMuRatio->Draw("EP");
	TLine *l1 = new TLine(Histo_DYToMuMuRatio->GetXaxis()->GetXmin(),1.,Histo_DYToMuMuRatio->GetXaxis()->GetXmax(),1.);
	l1->SetLineWidth(1.5); 
	l1->Draw("same");
	c1->Update();
*/
	c1->cd(1);
  TLatex latexLabel;
//  std::ostringstream intLumiString;
//  intLumiString << setprecision (2) << fixed << integratedLuminosity;
//  string intLumiText = "#intL= " + intLumiString.str() + " pb^{-1}";
  latexLabel.SetTextSize(0.03);
  latexLabel.SetNDC();
  latexLabel.DrawLatex(0.13, 0.96, "CMS Preliminary 2011");
  latexLabel.DrawLatex(0.42, 0.96, "#sqrt{s} = 7 TeV");
//  latexLabel.DrawLatex(0.57, 0.96, intLumiText.c_str());

  // // Sixth: update canvas
  c1->Update();
  c1->Draw();

  // Print the canvas
  string PicName="gif/DataData_" + pic + "_" + name + ".gif";
  c1->Print(PicName.c_str());
  PicName="eps/DataData_" + pic + "_" + name + ".eps";
  c1->Print(PicName.c_str());
  string convert = "convert eps/DataData_" + pic + "_" + name + ".eps" + " pdf/DataData_" + pic + "_" + name + ".pdf";
  system(convert.c_str());


  if (inlog==true) {
    c1->cd(1);
    Histo_Data->SetMaximum(YMax_log);
    Histo_Data->SetMinimum(YMin_log);
    Histo_Data->GetYaxis()->SetRangeUser(YMin_log, YMax_log);
//    Histo_DYToMuMu->SetMaximum(YMax_log);
//    Histo_DYToMuMu->SetMinimum(YMin_log);
//    Histo_DYToMuMu->GetYaxis()->SetRangeUser(YMin_log, YMax_log);
//    Histo_QCD_Pt15->SetMaximum(YMax_log);
//    Histo_QCD_Pt15->SetMinimum(YMin_log);
//    Histo_QCD_Pt15->GetYaxis()->SetRangeUser(YMin_log, YMax_log);

    Histo_OLD_Data->SetMaximum(YMax_log);
    Histo_OLD_Data->SetMinimum(YMin_log);
    Histo_OLD_Data->GetYaxis()->SetRangeUser(YMin_log, YMax_log);

//    Histo_WJetsToLNu->SetMaximum(YMax_log);
//    Histo_WJetsToLNu->SetMinimum(YMin_log);
//    Histo_WJetsToLNu->GetYaxis()->SetRangeUser(YMin_log, YMax_log);

//    Histo_QCDMu->SetMaximum(YMax_log);
//    Histo_QCDMu->SetMinimum(YMin_log);
//    Histo_QCDMu->GetYaxis()->SetRangeUser(YMin_log, YMax_log);

//    Histo_InclusiveMu15->SetMaximum(YMax_log);
//    Histo_InclusiveMu15->SetMinimum(YMin_log);
//    Histo_InclusiveMu15->GetYaxis()->SetRangeUser(YMin_log, YMax_log);

//    Histo_OLD_DYToMuMu->SetMaximum(YMax_log);
//    Histo_OLD_DYToMuMu->SetMinimum(YMin_log);
//    Histo_OLD_DYToMuMu->GetYaxis()->SetRangeUser(YMin_log, YMax_log);

    c1->cd(1)->SetLogy(1);
    c1->Update();
    c1->Draw();
    string PicName_log="gif/DataData_" + pic + "_" + name + "_log.gif";
    c1->Print(PicName_log.c_str());
    PicName="eps/DataData_" + pic + "_" + name + "_log.eps";
    c1->Print(PicName.c_str());
    string convert = "convert eps/DataData_" + pic + "_" + name + "_log.eps" + " pdf/DataData_" + pic + "_" + name + "_log.pdf";
    system(convert.c_str());
    c1->cd(1)->SetLogy(0);
    c1->Update();
  }


  // Clean the memory
  c1->Clear();
  legend->Clear();
//	Histo_template->Delete();
  Histo_Data_temp->Delete();
  Histo_Data->Delete();
//  Histo_DYToMuMu_temp->Delete();
//  Histo_DYToMuMu->Delete();

//  Histo_OLD_DYToMuMu_temp->Delete();
//  Histo_OLD_DYToMuMu->Delete();

  Histo_OLD_Data_temp->Delete();
  Histo_OLD_Data->Delete();

//  Histo_WJetsToLNu_temp->Delete();
//  Histo_WJetsToLNu->Delete();

//  Histo_QCDMu_temp->Delete();
//  Histo_QCDMu->Delete();

}

