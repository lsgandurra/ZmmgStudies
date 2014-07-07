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
		//MC_histo = new TH1D("MC_histo","MC_histo",35, 0, 35);
		Double_t Summer2011_S13[35] = {
		               1.30976e-05,
                               0.000148266,
                               0.00226073,
                               0.030543,
                               0.0868303,
                               0.120295,
                               0.124687,
                               0.110419,
                               0.0945742,
                               0.0837875,
                               0.0774277,
                               0.0740595,
                               0.0676844,
                               0.0551203,
                               0.0378357,
                               0.0210203,
                               0.00918262,
                               0.00309786,
                               0.000808509,
                               0.000168568,
                               3.02344e-05,
                               5.16455e-06,
                               8.83185e-07,
                               1.43975e-07,
                               2.07228e-08,
                               2.51393e-09,
                               2.52072e-10,
                               2.07328e-11,
                               1.39369e-12,
                               7.63843e-14,
                               3.4069e-15,
                               1.23492e-16,
                               3.63059e-18,
                               8.53277e-20,
                               1.33668e-22};

		

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


			Double_t Summer2012_RD1_AB[60] = {
			    5.62384e-12,
			    5.37563e-11,
			    3.24138e-08,
			    2.27084e-05,
			    3.38197e-05,
			    0.000334948,
			    0.00228691,
			    0.00603839,
			    0.0123791,
			    0.020876,
			    0.0302324,
			    0.0401893,
			    0.0486242,
			    0.0545462,
			    0.0591351,
			    0.0624621,
			    0.0628103,
			    0.0609193,
			    0.058327,
			    0.0558273,
			    0.0540694,
			    0.053264,
			    0.0529594,
			    0.0523153,
			    0.050163,
			    0.0456089,
			    0.0386456,
			    0.0300892,
			    0.021204,
			    0.0133216,
			    0.00737675,
			    0.00357903,
			    0.00152009,
			    0.000566906,
			    0.000186918,
			    5.5118e-05,
			    1.47867e-05,
			    3.6868e-06,
			    8.71796e-07,
			    1.97817e-07,
			    4.30634e-08	
			};

			Double_t Summer2012_RD1_C[60] = {
	                1.64489e-08,
	                1.39807e-07,
	                1.101e-06,
	                7.51692e-06,
	                4.23423e-05,
	                0.000191918,
	                0.000692231,
	                0.0019829,
	                0.00458417,
	                0.00963676,
	                0.0230679,
	                0.0500937,
	                0.0733726,
	                0.0815892,
	                0.0828382,
	                0.0810634,
	                0.075564,
	                0.0687302,
	                0.0631491,
	                0.0583158,
	                0.0540233,
	                0.0502324,
	                0.0465909,
	                0.0426951,
	                0.0379946,
	                0.031998,
	                0.0247511,
	                0.0171058,
	                0.010353,
	                0.00542871,
	                0.00246064,
	                0.000968498,
	                0.000334249,
	                0.000102514,
	                2.84199e-05,
	                7.29052e-06,
	                1.79136e-06,
	                4.4065e-07,
	                1.12483e-07,
	                2.99082e-08,
	                8.06412e-09	
			};

			Double_t Summer2012_RD1_D[60] = {
	                2.54656e-11,
	                6.70051e-11,
	                2.74201e-06,
	                6.9111e-06,
	                5.00919e-06,
	                6.24538e-05,
	                0.000338679,
	                0.000892795,
	                0.00237358,
	                0.00686023,
	                0.0144954,
	                0.026012,
	                0.0360377,
	                0.0420151,
	                0.0457901,
	                0.0482319,
	                0.0503176,
	                0.052569,
	                0.0546253,
	                0.0561205,
	                0.0568903,
	                0.0570889,
	                0.0566598,
	                0.0553747,
	                0.0531916,
	                0.0501454,
	                0.0463101,
	                0.0417466,
	                0.0364842,
	                0.0306443,
	                0.0245417,
	                0.0186276,
	                0.0133446,
	                0.00900314,
	                0.00571947,
	                0.00342706,
	                0.00194292,
	                0.00104671,
	                0.000538823,
	                0.000266973,
	                0.000128572,
	                6.09778e-05,
	                2.89549e-05,
	                1.40233e-05,
	                7.04619e-06,
	                3.71289e-06,
	                2.055e-06,
	                1.18713e-06,
	                7.08603e-07,
	                4.32721e-07,
	                2.6817e-07,
	                1.67619e-07,
	                1.05157e-07,
	                6.59446e-08,
	                4.11915e-08,
	                2.55494e-08	
                        };


			//for (int i = 1; i <= 35; i++)
			for (int i = 1; i <= 60; i++) 
			{
				
    				//MC_histo->SetBinContent(i,Summer2012_S10[i-1]); 
  				//MC_histo->SetBinContent(i,Summer2012_S7[i-1]);
				//MC_histo->SetBinContent(i,Summer2012_RD1_AB[i-1]);
				//MC_histo->SetBinContent(i,Summer2012_RD1_A[i-1]);
				//MC_histo->SetBinContent(i,Summer2012_RD1_B[i-1]);
				//MC_histo->SetBinContent(i,Summer2012_RD1_C[i-1]);
				//MC_histo->SetBinContent(i,Summer2012_RD1_D[i-1]);
				MC_histo->SetBinContent(i,Summer2012_S10[i-1]);
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
	for(int i=1; i <= MC_weights->GetNbinsX(); i++)
	{ cout << MC_weights->GetBinContent(i) << ", ";  }
	cout << endl;	

	TFile* weights_file = new TFile(Form("%s",output.c_str()), "recreate");
	weights_file->cd();
	MC_weights->Write();
	MC_reweighted->Write();
	data_pdf->Write();
	weights_file->Close();

		


	return 0;
}
