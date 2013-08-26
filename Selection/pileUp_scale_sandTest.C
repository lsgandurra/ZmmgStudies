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

#include "setTDRStyle.C"

#pragma optimize 0

using namespace std;

//int main()
int pileUp_scale_sandTest(string output = "puTest.root")
{

  	setTDRStyle();
  	gStyle->SetOptTitle(0);
  	gStyle->SetOptStat(0);

	TCanvas *c1 = new TCanvas("c1", "c1",0,0,600,600);	
/*
	TFile* data_file = new TFile(Data.c_str());
	TH1D* data_histo = (TH1D*)data_file->Get("pileup");
	double data_integral = data_histo->Integral();
	cout << "data_integral= " << data_integral << endl;
	TH1D* data_pdf = new TH1D("data_pdf", "data_pdf", 60, 0, 60);
	data_pdf->Add(data_histo);
	data_pdf->Scale((double)(1.0)/(double)(data_integral) );
	cout << data_pdf->GetBinContent(0) << "\t" << data_pdf->GetBinContent(1) << endl;
	data_pdf->Draw();
*/

	
	TChain * data_tree = new TChain("miniTree");
	data_tree->Add("miniTree_Run2012A_13Jul2012_v1_NewMuonID_NewSelection_0_injRe0_v4_partALL.root");
	data_tree->Add("miniTree_Run2012A_recover_06Aug2012_v1_NewMuonID_NewSelection_0_injRe0_v4_partALL.root");
        data_tree->Add("miniTree_Run2012B_13Jul2012_v4_NewMuonID_NewSelection_0_injRe0_v4_partALL.root");
        data_tree->Add("miniTree_Run2012C-24Aug2012-v1_NewMuonID_NewSelection_0_injRe0_v4_partALL.root");
        data_tree->Add("miniTree_Run2012C-EcalRecover_11Dec2012-v1_NewMuonID_NewSelection_0_injRe0_v4_partALL.root");
        data_tree->Add("miniTree_Run2012C_PromptReco_v2_NewMuonID_NewSelection_0_injRe0_v4_partALL.root");
        data_tree->Add("miniTree_Run2012D_PromptReco_v1_NewMuonID_NewSelection_0_injRe0_v4_partALL.root");
	TH1D* data_histo = new TH1D("data_histo","data_histo",60, 0, 60);
        data_tree->Draw("nVertices>>data_histo");	
	double data_integral = data_histo->Integral();
	TH1D* data_pdf = new TH1D("data_pdf", "data_pdf", 60, 0, 60);
        data_pdf->Add(data_histo);
        data_pdf->Scale((double)(1.0)/(double)(data_integral) );
        cout << data_pdf->GetBinContent(0) << "\t" << data_pdf->GetBinContent(1) << endl;
        data_pdf->Draw();	
	


	//TFile* MC_file = new TFile(MC.c_str());
	//TTree* MC_tree = (TTree*) MC_file->Get("miniTree");
	TChain * MC_tree = new TChain("miniTree");
	MC_tree->Add("miniTree_DYToMuMu_Summer12_NewMuonID_NewSelection_1_injRe0_v4_partALL.root");
	MC_tree->Add("miniTree_DYToMuMu_Summer12_NewMuonID_NewSelection_2_injRe0_v4_partALL.root");
	TH1D* MC_histo_temp = new TH1D();
	MC_tree->Draw("nVertices>>MC_histo_temp(60, 0, 60)");
	TH1D* MC_histo = (TH1D*)gDirectory->Get("MC_histo_temp");

	double MC_integral = MC_histo->Integral();
	cout << "MC_integral= " << MC_integral << endl;
	TH1D* MC_pdf = new TH1D("MC_pdf", "MC_pdf", 60, 0, 60);
	MC_pdf->Add(MC_histo);
	MC_pdf->Scale((double)(1.0)/(double)(MC_integral) );
	TH1D* MC_weights = new TH1D("MC_weights", "MC_weights", 60, 0, 60);
	MC_weights->Divide(data_pdf, MC_pdf);

	TH1D* MC_reweighted = new TH1D("MC_reweighted", "MC_reweighted", 60, 0, 60);
	MC_reweighted->Add(MC_histo);
	MC_reweighted->Multiply(MC_weights);
	double MC_reweighted_integral = MC_reweighted->Integral();
	MC_reweighted->Scale((double)(1.0)/(double)(MC_reweighted_integral));

	cout << "MC" << endl;
	for(int i=0; i < MC_weights->GetNbinsX(); i++){ cout << MC_weights->GetBinContent(i) << ", ";  }
	cout << endl;	

	TFile* weights_file = new TFile(Form("%s",output.c_str()), "recreate");
	weights_file->cd();
	MC_weights->Write();
	MC_reweighted->Write();
	data_pdf->Write();
	weights_file->Close();

		


	return 0;
}
