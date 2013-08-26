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

int main(int argc, char *argv[])
//int pileUp_scaleNvtxRECO(string Data = "MyDataPileupHistogram_69400_ABCDNew_final.root", string MC = "miniTree_pileUp_DYToMuMu_Summer12_NewMuonID_partALL.root", string output = "pileUpweights_DYToMuMu_Summer12.root")
{
	cout << "argc= " << argc << endl;
        for(int iarg = 0 ; iarg < argc; iarg++)
        {
                cout << "argv[" << iarg << "]= " << argv[iarg] << endl;
        }
	
	if( argc == 1 )
        {
                cerr << "arguments should be passed !! (Data) (MC) (output) " << endl;
                return 1;
        }
		
	string Data = "";
	string MC = "";
	string output = "";
	string MCOfficials = "true";

	if( argc > 1 ) Data = argv[1];
        if( argc > 2 ) MC = argv[2];
        if( argc > 3 ) output = argv[3];
	if( argc > 4 ) output = argv[4];

  	setTDRStyle();
  	gStyle->SetOptTitle(0);
  	gStyle->SetOptStat(0);

	TCanvas *c1 = new TCanvas("c1", "c1",0,0,600,600);	

	TFile* data_file = new TFile(Data.c_str());
	TH1D* data_histo = (TH1D*)data_file->Get("pileup");
	double data_integral = data_histo->Integral();
	cout << "data_integral= " << data_integral << endl;
	TH1D* data_pdf = new TH1D("data_pdf", "data_pdf", 60, 0, 60);
	data_pdf->Add(data_histo);
	data_pdf->Scale((double)(1.0)/(double)(data_integral) );
	cout << data_pdf->GetBinContent(0) << "\t" << data_pdf->GetBinContent(1) << endl;
	data_pdf->Draw();

	/*
	TChain * MC_tree = new TChain("miniTree");
        MC_tree->Add("miniTree_DYToMuMu_Summer12_NewMuonID_NewSelection_1_injRe0_v4_partALL.root");
        MC_tree->Add("miniTree_DYToMuMu_Summer12_NewMuonID_NewSelection_2_injRe0_v4_partALL.root");
        //MC_tree->Add("miniTree_TTJets_Summer12_S7_NewMuonID_NewSelection_3_injRe0_v5_partALL.root");
        //MC_tree->Add("miniTree_WJetsToLNu_Summer12_S10_NewMuonID_NewSelection_3_injRe0_v5_partALL.root");*/
	

	TFile* MC_file = 0;
	TTree* MC_tree = 0;
	TH1D* MC_histo_temp = new TH1D();
	TH1D* MC_histo = 0;


	if(MCOfficials == "false")
	{
		MC_file = new TFile(MC.c_str());
		MC_tree = (TTree*) MC_file->Get("miniTree");
		MC_tree->Draw("nVertices>>MC_histo_temp(60, 0, 60)");
		MC_histo = (TH1D*)gDirectory->Get("MC_histo_temp");
	}

	if(MCOfficials == "true")
        {
		MC_histo = new TH1D("MC_histo","MC_histo",60, 0, 60);
		Double_t Summer2012_S10[60] = {
                         2.560E-06,
                         5.239E-06,
                         1.420E-05,
                         5.005E-05,
                         1.001E-04,
                         2.705E-04,
                         1.999E-03,
                         6.097E-03,
                         1.046E-02,
                         1.383E-02,
                         1.685E-02,
                         2.055E-02,
                         2.572E-02,
                         3.262E-02,
                         4.121E-02,
                         4.977E-02,
                         5.539E-02,
                         5.725E-02,
                         5.607E-02,
                         5.312E-02,
                         5.008E-02,
                         4.763E-02,
                         4.558E-02,
                         4.363E-02,
                         4.159E-02,
                         3.933E-02,
                         3.681E-02,
                         3.406E-02,
                         3.116E-02,
                         2.818E-02,
                         2.519E-02,
                         2.226E-02,
                         1.946E-02,
                         1.682E-02,
                         1.437E-02,
                         1.215E-02,
                         1.016E-02,
                         8.400E-03,
                         6.873E-03,
                         5.564E-03,
                         4.457E-03,
                         3.533E-03,
                         2.772E-03,
                         2.154E-03,
                         1.656E-03,
                         1.261E-03,
                         9.513E-04,
                         7.107E-04,
                         5.259E-04,
                         3.856E-04,
                         2.801E-04,
                         2.017E-04,
                         1.439E-04,
                         1.017E-04,
                         7.126E-05,
                         4.948E-05,
                         3.405E-05,
                         2.322E-05,
                         1.570E-05,
                         5.005E-06};	


			Double_t Summer2012_S7[60] = {
			    2.344E-05,
			    2.344E-05,
			    2.344E-05,
			    2.344E-05,
			    4.687E-04,
			    4.687E-04,
			    7.032E-04,
			    9.414E-04,
			    1.234E-03,
			    1.603E-03,
			    2.464E-03,
			    3.250E-03,
			    5.021E-03,
			    6.644E-03,
			    8.502E-03,
			    1.121E-02,
			    1.518E-02,
			    2.033E-02,
			    2.608E-02,
			    3.171E-02,
			    3.667E-02,
			    4.060E-02,
			    4.338E-02,
			    4.520E-02,
			    4.641E-02,
			    4.735E-02,
			    4.816E-02,
			    4.881E-02,
			    4.917E-02,
			    4.909E-02,
			    4.842E-02,
			    4.707E-02,
			    4.501E-02,
			    4.228E-02,
			    3.896E-02,
			    3.521E-02,
			    3.118E-02,
			    2.702E-02,
			    2.287E-02,
			    1.885E-02,
			    1.508E-02,
			    1.166E-02,
			    8.673E-03,
			    6.190E-03,
			    4.222E-03,
			    2.746E-03,
			    1.698E-03,
			    9.971E-04,
			    5.549E-04,
			    2.924E-04,
			    1.457E-04,
			    6.864E-05,
			    3.054E-05,
			    1.282E-05,
			    5.081E-06,
			    1.898E-06,
			    6.688E-07,
			    2.221E-07,
			    6.947E-08,
			    2.047E-08
			   };

			for (int i = 1; i <= 60; i++) 
			{
    				//MC_histo->SetBinContent(i,Summer2012_S10[i-1]); 
  				MC_histo->SetBinContent(i,Summer2012_S7[i-1]);
			}	

	}

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
