// ROOT HEADERS
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
#include <TBranch.h>
#include <TString.h>
#include <TBits.h>
#include <TMath.h>
#include "TROOT.h"
#include <TLatex.h>
#include "TProof.h"

// C++ HEADERS
#include <string.h>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>

// DISPLAY HEADERS
#include "CMSStyle.C"

#pragma optimize 0

using namespace std;


//int plot2D()
int main()
{
//	gROOT->ProcessLine(".x setTDRStyle.C");
	CMSstyle();
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);
	gStyle->SetMarkerSize(0.1);
	gStyle->SetPalette(1);

//  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadRightMargin(0.15);
  gStyle->SetPadBottomMargin(0.17);
  gStyle->SetPadLeftMargin(0.17);



	TCanvas *c1 = new TCanvas("c1", "c1");

	vector<string> x_file;
	vector<string> x_variable;
	vector<string> x_display;
	vector<int> x_nbins;
	vector< vector<double> > x_bins;
	vector<double> y_min;
	vector<double> y_max;
	vector<string> y_file;
	vector<string> y_variable;
	vector<string> y_display;
/*
	x_file.push_back("MuonM_Phi");
	x_variable.push_back("MuonM_Phi");
	x_display.push_back("#phi^{#mu^{-}}");
	x_nbins.push_back(12);
	double t_MuonM_Phi[] = {-3.15, -2.5, -2., -1.5, -1., -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.15};
	vector<double> vt_MuonM_Phi(t_MuonM_Phi+0, t_MuonM_Phi+13);
	x_bins.push_back(vt_MuonM_Phi);
	y_file.push_back("mmg_sX100");
	y_variable.push_back("mmg_s * 100");
  y_display.push_back("s (%)");
	y_min.push_back(-100.);
	y_max.push_back(100.);

	x_file.push_back("MuonM_Eta");
	x_variable.push_back("MuonM_Eta");
	x_display.push_back("#eta^{#mu^{-}}");
	x_nbins.push_back(12);
	double t_MuonM_Eta[] = {-3., -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.};
	vector<double> vt_MuonM_Eta(t_MuonM_Eta+0, t_MuonM_Eta+13);
	x_bins.push_back(vt_MuonM_Eta);
	y_file.push_back("mmg_sX100");
	y_variable.push_back("mmg_s * 100");
  y_display.push_back("s (%)");
	y_min.push_back(-100.);
	y_max.push_back(100.);

	x_file.push_back("MuonM_Pt");
	x_variable.push_back("MuonM_Pt");
	x_display.push_back("E_{T}^{#mu^{-}}");
	x_nbins.push_back(6);
	double t_MuonM_Pt[] = {10., 12., 15., 20., 25., 30., 100.};
	vector<double> vt_MuonM_Pt(t_MuonM_Pt+0, t_MuonM_Pt+7);
	x_bins.push_back(vt_MuonM_Pt);
	y_file.push_back("mmg_sX100");
	y_variable.push_back("mmg_s * 100");
  y_display.push_back("s (%)");
	y_min.push_back(-100.);
	y_max.push_back(100.);


	x_file.push_back("MuonP_Phi");
	x_variable.push_back("MuonP_Phi");
	x_display.push_back("#phi^{#mu^{+}}");
	x_nbins.push_back(12);
	double t_MuonP_Phi[] = {-3.15, -2.5, -2., -1.5, -1., -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.15};
	vector<double> vt_MuonP_Phi(t_MuonP_Phi+0, t_MuonP_Phi+13);
	x_bins.push_back(vt_MuonP_Phi);
	y_file.push_back("mmg_sX100");
	y_variable.push_back("mmg_s * 100");
  y_display.push_back("s (%)");
	y_min.push_back(-100.);
	y_max.push_back(100.);

	x_file.push_back("MuonP_Eta");
	x_variable.push_back("MuonP_Eta");
	x_display.push_back("#eta^{#mu^{+}}");
	x_nbins.push_back(12);
	double t_MuonP_Eta[] = {-3., -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.};
	vector<double> vt_MuonP_Eta(t_MuonP_Eta+0, t_MuonP_Eta+13);
	x_bins.push_back(vt_MuonP_Eta);
	y_file.push_back("mmg_sX100");
	y_variable.push_back("mmg_s * 100");
  y_display.push_back("s (%)");
	y_min.push_back(-100.);
	y_max.push_back(100.);

	x_file.push_back("MuonP_Pt");
	x_variable.push_back("MuonP_Pt");
	x_display.push_back("E_{T}^{#mu^{+}}");
	x_nbins.push_back(6);
	double t_MuonP_Pt[] = {10., 12., 15., 20., 25., 30., 100.};
	vector<double> vt_MuonP_Pt(t_MuonP_Pt+0, t_MuonP_Pt+7);
	x_bins.push_back(vt_MuonP_Pt);
	y_file.push_back("mmg_sX100");
	y_variable.push_back("mmg_s * 100");
  y_display.push_back("s (%)");
	y_min.push_back(-100.);
	y_max.push_back(100.);
*/
/*
	x_variable.push_back("MuonF_Phi");
	x_display.push_back("#phi^{#mu^{far}}");
	x_variable.push_back("MuonF_Eta");
	x_display.push_back("#eta^{#mu^{far}}");
	x_variable.push_back("MuonF_Pt");
	x_display.push_back("E_{T}^{#mu^{far}}");

	x_variable.push_back("MuonN_Phi");
	x_display.push_back("#phi^{#mu^{near}}");
	x_variable.push_back("MuonN_Eta");
	x_display.push_back("#eta^{#mu^{near}}");
	x_variable.push_back("MuonN_Pt");
	x_display.push_back("E_{T}^{#mu^{near}}");

	x_variable.push_back("MuonL_Phi");
	x_display.push_back("#phi^{#mu^{lead}}");
	x_variable.push_back("MuonL_Eta");
	x_display.push_back("#eta^{#mu^{lead}}");
	x_variable.push_back("MuonL_Pt");
	x_display.push_back("E_{T}^{#mu^{lead}}");

	x_variable.push_back("MuonS_Phi");
	x_display.push_back("#phi^{#mu^{trail}}");
	x_variable.push_back("MuonS_Eta");
	x_display.push_back("#eta^{#mu^{trail}}");
	x_variable.push_back("MuonS_Pt");
	x_display.push_back("E_{T}^{#mu^{trail}}");
*/

	x_file.push_back("deltaRNear");
	x_variable.push_back("deltaRNear");
	x_display.push_back("#Delta R (#gamma, #mu^{near})");
	x_nbins.push_back(10);
	double t_deltaRNear[] = {0., 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0};
	vector<double> vt_deltaRNear(t_deltaRNear+0, t_deltaRNear+11);
	x_bins.push_back(vt_deltaRNear);
	y_file.push_back("mmg_sX100");
	y_variable.push_back("mmg_s * 100");
  y_display.push_back("s (%)");
	y_min.push_back(-100.);
	y_max.push_back(100.);

	x_file.push_back("deltaRFar");
	x_variable.push_back("deltaRFar");
	x_display.push_back("#Delta R (#gamma, #mu^{far})");
	x_nbins.push_back(10);
	double t_deltaRFar[] = {0., 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0};
	vector<double> vt_deltaRFar(t_deltaRFar+0, t_deltaRFar+11);
	x_bins.push_back(vt_deltaRFar);
	y_file.push_back("mmg_sX100");
	y_variable.push_back("mmg_s * 100");
  y_display.push_back("s (%)");
	y_min.push_back(-100.);
	y_max.push_back(100.);

	x_file.push_back("deltaRMinus");
	x_variable.push_back("deltaRMinus");
	x_display.push_back("#Delta R (#gamma, #mu^{-})");
	x_nbins.push_back(10);
	double t_deltaRMinus[] = {0., 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0};
	vector<double> vt_deltaRMinus(t_deltaRMinus+0, t_deltaRMinus+11);
	x_bins.push_back(vt_deltaRMinus);
	y_file.push_back("mmg_sX100");
	y_variable.push_back("mmg_s * 100");
  y_display.push_back("s (%)");
	y_min.push_back(-100.);
	y_max.push_back(100.);

	x_file.push_back("deltaRPlus");
	x_variable.push_back("deltaRPlus");
	x_display.push_back("#Delta R (#gamma, #mu^{+})");
	x_nbins.push_back(10);
	double t_deltaRPlus[] = {0., 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0};
	vector<double> vt_deltaRPlus(t_deltaRPlus+0, t_deltaRPlus+11);
	x_bins.push_back(vt_deltaRPlus);
	y_file.push_back("mmg_sX100");
	y_variable.push_back("mmg_s * 100");
  y_display.push_back("s (%)");
	y_min.push_back(-100.);
	y_max.push_back(100.);

	x_file.push_back("deltaRLeading");
	x_variable.push_back("deltaRLeading");
	x_display.push_back("#Delta R (#gamma, #mu^{leading})");
	x_nbins.push_back(10);
	double t_deltaRLeading[] = {0., 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0};
	vector<double> vt_deltaRLeading(t_deltaRLeading+0, t_deltaRLeading+11);
	x_bins.push_back(vt_deltaRLeading);
	y_file.push_back("mmg_sX100");
	y_variable.push_back("mmg_s * 100");
  y_display.push_back("s (%)");
	y_min.push_back(-100.);
	y_max.push_back(100.);

	x_file.push_back("deltaRSubleading");
	x_variable.push_back("deltaRSubleading");
	x_display.push_back("#Delta R (#gamma, #mu^{subleading})");
	x_nbins.push_back(10);
	double t_deltaRSubleading[] = {0., 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0};
	vector<double> vt_deltaRSubleading(t_deltaRSubleading+0, t_deltaRSubleading+11);
	x_bins.push_back(vt_deltaRSubleading);
	y_file.push_back("mmg_sX100");
	y_variable.push_back("mmg_s * 100");
  y_display.push_back("s (%)");
	y_min.push_back(-100.);
	y_max.push_back(100.);

	x_file.push_back("Photon_Et");
	x_variable.push_back("Photon_Et");
	x_display.push_back("E_{T}^{#gamma}");
	x_nbins.push_back(6);
	double t_Photon_Et[] = {10., 12., 15., 20., 25., 30., 100.};
	vector<double> vt_Photon_Et(t_Photon_Et+0, t_Photon_Et+7);
	x_bins.push_back(vt_Photon_Et);
	y_file.push_back("mmg_sX100");
	y_variable.push_back("mmg_s * 100");
  y_display.push_back("s (%)");
	y_min.push_back(-100.);
	y_max.push_back(100.);
/*
	x_file.push_back("Photon_Eta");
	x_variable.push_back("Photon_Eta");
	x_display.push_back("#eta^{#gamma}");
	x_nbins.push_back(8);
	double t_Photon_Eta[] = {-3., -2.5, -1.566, -1.4442, 0.0, 1.4442, 1.566, 2.5, 3.};
	vector<double> vt_Photon_Eta(t_Photon_Eta+0, t_Photon_Eta+9);
	x_bins.push_back(vt_Photon_Eta);
	y_file.push_back("mmg_sX100");
	y_variable.push_back("mmg_s * 100");
  y_display.push_back("s (%)");
	y_min.push_back(-100.);
	y_max.push_back(100.);
	

	x_file.push_back("Photon_Phi");
	x_variable.push_back("Photon_Phi");
	x_display.push_back("#phi^{#gamma}");
	x_nbins.push_back(12);
	double t_Photon_Phi[] = {-3.15, -2.5, -2., -1.5, -1., -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.15};
	vector<double> vt_Photon_Phi(t_Photon_Phi+0, t_Photon_Phi+13);
	x_bins.push_back(vt_Photon_Phi);
	y_file.push_back("mmg_sX100");
	y_variable.push_back("mmg_s * 100");
  y_display.push_back("s (%)");
	y_min.push_back(-100.);
	y_max.push_back(100.);

	x_file.push_back("Photon_Et");
	x_variable.push_back("Photon_Et");
	x_display.push_back("E_{T}^{#gamma}");
	x_nbins.push_back(6);
	x_bins.push_back(vt_Photon_Et);
	y_file.push_back("MuonF_Pt_RECO_o_true");
	y_variable.push_back("sqrt((MuonF_Pt * MuonF_Pt)/(MuonF_MC_Px*MuonF_MC_Px + MuonF_MC_Py*MuonF_MC_Py))");
	y_display.push_back("p_{T}^{RECO #mu far} / p_{T}^{GEN #mu far}");
	y_min.push_back(.5);
	y_max.push_back(1.5);

	x_file.push_back("Photon_Et");
	x_variable.push_back("Photon_Et");
	x_display.push_back("E_{T}^{#gamma}");
	x_nbins.push_back(6);
	x_bins.push_back(vt_Photon_Et);
	y_file.push_back("MuonN_Pt_RECO_o_true");
	y_variable.push_back("sqrt((MuonN_Pt * MuonN_Pt)/(MuonN_MC_Px*MuonN_MC_Px + MuonN_MC_Py*MuonN_MC_Py))");
	y_display.push_back("p_{T}^{RECO #mu near} / p_{T}^{GEN #mu near}");
	y_min.push_back(.5);
	y_max.push_back(1.5);

	x_file.push_back("Photon_Et");
	x_variable.push_back("Photon_Et");
	x_display.push_back("E_{T}^{#gamma}");
	x_nbins.push_back(6);
	x_bins.push_back(vt_Photon_Et);
	y_file.push_back("Mmumu");
	y_variable.push_back("Mmumu");
  y_display.push_back("m_{#mu#mu} (GeV)");
	y_min.push_back(30.);
	y_max.push_back(100.);
*/


	x_file.push_back("Photon_Et");
	x_variable.push_back("Photon_Et");
	x_display.push_back("E_{T}^{#gamma}");
	x_nbins.push_back(6);
	x_bins.push_back(vt_Photon_Et);
	y_file.push_back("mmg_s_Photon_MC_o_mmg_s");
	y_variable.push_back("mmg_s_Photon_MC / mmg_s");
	y_display.push_back("s'_{RECO}, #gamma^{RECO} #rightarrow #gamma^{GEN} / s_{RECO}");
	y_min.push_back(.5);
	y_max.push_back(1.5);

	x_file.push_back("Photon_Et");
	x_variable.push_back("Photon_Et");
	x_display.push_back("E_{T}^{#gamma}");
	x_nbins.push_back(6);
	x_bins.push_back(vt_Photon_Et);
	y_file.push_back("mmg_s_Muons_MC_o_mmg_s");
	y_variable.push_back("mmg_s_Muons_MC / mmg_s");
	y_display.push_back("s'_{RECO}, #mu^{RECO} #rightarrow #mu^{GEN} / s_{RECO}");
	y_min.push_back(.5);
	y_max.push_back(1.5);

	x_file.push_back("Photon_Et");
	x_variable.push_back("Photon_Et");
	x_display.push_back("E_{T}^{#gamma}");
	x_nbins.push_back(6);
	x_bins.push_back(vt_Photon_Et);
	y_file.push_back("mmg_s_MMG_MC_o_mmg_s");
	y_variable.push_back("mmg_s_MMG_MC / mmg_s");
	y_display.push_back("s'_{RECO}, (#gamma^{RECO},#mu^{RECO}) #rightarrow (#gamma^{GEN},#mu^{GEN}) / s_{RECO}");
	y_min.push_back(.5);
	y_max.push_back(1.5);

	x_file.push_back("Photon_Et");
	x_variable.push_back("Photon_Et");
	x_display.push_back("E_{T}^{#gamma}");
	x_nbins.push_back(6);
	x_bins.push_back(vt_Photon_Et);
	y_file.push_back("mmg_s_MZ_o_mmg_s");
	y_variable.push_back("mmg_s_MZ / mmg_s");
	y_display.push_back("s'_{RECO}, m_{Z^{0}} #rightarrow m_{#mu#mu#gamma}^{GEN} / s_{RECO}");
	y_min.push_back(.5);
	y_max.push_back(1.5);

	x_file.push_back("Photon_Et");
	x_variable.push_back("Photon_Et");
	x_display.push_back("E_{T}^{#gamma}");
	x_nbins.push_back(6);
	x_bins.push_back(vt_Photon_Et);
	y_file.push_back("mmg_s_MZ_Photon_MC_o_mmg_s");
	y_variable.push_back("mmg_s_MZ_Photon_MC / mmg_s");
	y_display.push_back("s'_{RECO}, (m_{Z^{0}}, #gamma^{RECO}) #rightarrow (m_{#mu#mu#gamma}^{GEN}, #gamma^{GEN}) / s_{RECO}");
	y_min.push_back(.5);
	y_max.push_back(1.5);

	x_file.push_back("Photon_Et");
	x_variable.push_back("Photon_Et");
	x_display.push_back("E_{T}^{#gamma}");
	x_nbins.push_back(6);
	x_bins.push_back(vt_Photon_Et);
	y_file.push_back("mmg_s_MZ_Muons_MC_o_mmg_s");
	y_variable.push_back("mmg_s_MZ_Muons_MC / mmg_s");
	y_display.push_back("s'_{RECO}, (m_{Z^{0}}, #mu^{RECO}) #rightarrow (m_{#mu#mu#gamma}^{GEN}, #mu^{GEN} / s_{RECO}");
	y_min.push_back(.5);
	y_max.push_back(1.5);



	vector<string> cat_cut;
	vector<string> cat_title;

	cat_cut.push_back("1");
	cat_title.push_back("tight");

//	cat_cut.push_back("Photon_isEB");
//	cat_title.push_back("EB-tight");
	cat_cut.push_back("Photon_isEB && Photon_r9 > .94");
	cat_title.push_back("EB-tight-highR9");
	cat_cut.push_back("Photon_isEB && Photon_r9 < .94");
	cat_title.push_back("EB-tight-lowR9");
//	cat_cut.push_back("Photon_isEE");
//	cat_title.push_back("EE-tight");
	cat_cut.push_back("Photon_isEE && Photon_r9 > .95");
	cat_title.push_back("EE-tight-highR9");
	cat_cut.push_back("Photon_isEE && Photon_r9 < .95");
	cat_title.push_back("EE-tight-lowR9");

	string tight_cut = "isTightMMG && isMultipleCandidate == 0";

//	for(int idata = 0 ; idata < 2 ; idata++)
	for(int idata = 1 ; idata < 2 ; idata++)
	{
		TChain *data = new TChain("miniTree");
//		data->Add( idata == 0 ? "miniTree_v16_DATA_875pb-1.root" : "miniTree_v16_DYToMuMu.root");
//		data->Add( idata == 0 ? "miniTree_v17_Run2011A-ZMu-Jul05-Aug05-PromptSkim-v6.root" : "miniTree_v17_DYToMuMu.root");
//			data->Add( idata == 0 ? "miniTree_v12_2011_40_80_partALL.root" : "miniTree_v10_DYToMuMu_S6_2011_40_80_partALL.root");
			data->Add( idata == 0 ? "miniTree_v12_2011_40_80_partALL.root" : "miniTree_v11_ALL_DYToMuMu_S6_2011_40_80_partALL.root");
		for(int icat = 0 ; icat < cat_cut.size() ; icat++)
		{
			for(int ivar = 0 ; ivar < x_variable.size() ; ivar++)
			{
				TH2F* hist = new TH2F("hist", "hist", x_nbins[ivar], &x_bins[ivar][0], 5000, y_min[ivar], y_max[ivar]);
//				hist->SetBins(x_nbins[ivar], &x_bins[ivar][0]);
//				data->Draw(Form("mmg_s:%s", x_variable[ivar].c_str()), Form("(%s && %s )* weight_pileUp", tight_cut.c_str(), cat_cut[icat].c_str()), "COLZ");
				data->Draw(Form("(%s):(%s)>>hist", y_variable[ivar].c_str(), x_variable[ivar].c_str()), Form("(%s && %s )* weight_pileUp", tight_cut.c_str(), cat_cut[icat].c_str()), "goff");
//				data->Draw(Form("(%s):(%s)>>hist", y_variable[ivar].c_str(), x_variable[ivar].c_str()), Form("(%s && %s )* weight_pileUp", tight_cut.c_str(), cat_cut[icat].c_str()), "COLZ");
				//hist = (TH2F*)(gROOT->FindObject("htemp"));
//				double *arr = &x_bins[0];
//				std::copy(x_bins[ivar].begin(), x_bins[ivar].end(), arr);
//cout << "oups" << endl;
//				hist->SetBit(TH1::kCanRebin);
//				hist->SetBins(x_nbins[ivar], arr);
//				hist->SetBins(x_nbins[ivar], &x_bins[ivar][0]);
//				hist->ProfileX()->Draw();
//				TProfile* hist_pfx = (TProfile*)(gROOT->FindObject("hist_pfx"));
				TProfile* hist_pfx = (TProfile*)(hist->ProfileX());
				hist_pfx->SetMarkerStyle(25);
				hist_pfx->SetMarkerSize(1);
				hist_pfx->Draw();
				hist_pfx->GetYaxis()->SetTitle(y_display[ivar].c_str());
				hist_pfx->GetXaxis()->SetTitle(x_display[ivar].c_str());
//				c1->Clear();
//				TProfile *profile = (TProfile*) hist->ProfileX();
//				profile->Draw();

				TLatex latexLabel;
			  latexLabel.SetNDC();
			  latexLabel.SetTextSize(0.04);
			  latexLabel.DrawLatex(0.13, 0.96, "CMS Preliminary 2011");
			  latexLabel.DrawLatex(0.40, 0.96, "#sqrt{s} = 7 TeV");
				latexLabel.SetTextSize(0.05);
			  latexLabel.DrawLatex(0.55, 0.96, idata == 0 ? "DATA" : "Simulation");
			  latexLabel.SetTextSize(0.04);
				latexLabel.DrawLatex(0.75, 0.96, cat_title[icat].c_str());
		
				c1->Print(Form("eps/%s_%s_VS_%s_%s.eps", idata == 0 ? "DATA" : "MC", y_file[ivar].c_str(), x_file[ivar].c_str(), cat_title[icat].c_str()));
				c1->Print(Form("gif/%s_%s_VS_%s_%s.gif", idata == 0 ? "DATA" : "MC", y_file[ivar].c_str(), x_file[ivar].c_str(), cat_title[icat].c_str()));
				c1->Print(Form("png/%s_%s_VS_%s_%s.png", idata == 0 ? "DATA" : "MC", y_file[ivar].c_str(), x_file[ivar].c_str(), cat_title[icat].c_str()));
				c1->Print(Form("C/%s_%s_VS_%s_%s.C", idata == 0 ? "DATA" : "MC", y_file[ivar].c_str(), x_file[ivar].c_str(), cat_title[icat].c_str()));
				system(Form("convert eps/%s_%s_VS_%s_%s.eps pdf/%s_%s_VS_%s_%s.pdf", idata == 0 ? "DATA" : "MC", y_file[ivar].c_str(), x_file[ivar].c_str(), cat_title[icat].c_str(), idata == 0 ? "DATA" : "MC", y_file[ivar].c_str(), x_file[ivar].c_str(), cat_title[icat].c_str()));
	
				hist->Clear();
				hist_pfx->Clear();
				c1->Clear();
//return 1;
			} // END OF LOOP OVER VARIABLES
		}	// END OF LOOP OVER CATEGORIES
	}	// END OF LOOP OVER DATA/MC

/*
	cat_cut.push_back("");
	cat_title.push_back("");
	cat_cut.push_back("");
	cat_title.push_back("");
*/
return 0;
}
