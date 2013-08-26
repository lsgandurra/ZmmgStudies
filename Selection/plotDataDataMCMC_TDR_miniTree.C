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

#include "DrawDataDataMCMC.h"

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



#pragma optimize 0

using namespace std;
	
//int plotDataMC_TDR_miniTree()
int main()
{
	string selection = "loose";
	gSystem->Load("libToto.so");
	gROOT->ProcessLine(".x setTDRStyle.C");
//	string Data = "miniTree_v2_DATA_2011_ALL.root"; 
//	string Data = "miniTree_v2_DoubleMu-Run2011A_160404-163369_vALL_v02.root"; 
//	string Data = "miniTree_v2_DoubleMu-Run2011A_160404-163757_vALL.root";
	string Data = "/sps/cms/obondu/CMSSW_4_1_2/src/Zmumugamma/Selection/miniTree_v4_DYToMuMu_ALL.root";
	string DYToMuMu = "/sps/cms/obondu/CMSSW_4_1_2/src/Zmumugamma/Selection/miniTree_v4_DoubleMu-Run2011A_160404-163869_vALL.root";
	string OLD_Data = "../Selection/miniTree_Run2010_ALL_v3.root";
	string OLD_DYToMuMu = "../Selection/miniTree_DYToMuMu_v3.root";
//	string OLD_DYToMuMu = "miniTree_v3_DYToMuMu_ALL.root";
//	string OLD_Data = "miniTree_v3_Run2010_ALL.root";
//	string OLD_Data = "miniTree_test.root";
//	string WJetsToLNu = "miniTree_v3_WJets.root";
//	string QCDMu = "miniTree_v3_QCD.root"; 

	TFile *Data_File = new TFile(Data.c_str());
	TTree* Data_miniTree = (TTree*) Data_File->Get("miniTree");
	TTree* Data_miniTree_allmuons = (TTree*) Data_File->Get("miniTree_allmuons");
	TTree* Data_miniTree_allphotons = (TTree*) Data_File->Get("miniTree_allphotons");
	TFile *DYToMuMu_File = new TFile(DYToMuMu.c_str());
	TTree* DYToMuMu_miniTree = (TTree*) DYToMuMu_File->Get("miniTree");
	TTree* DYToMuMu_miniTree_allmuons = (TTree*) DYToMuMu_File->Get("miniTree_allmuons");
	TTree* DYToMuMu_miniTree_allphotons = (TTree*) DYToMuMu_File->Get("miniTree_allphotons");
	TFile *OLD_DYToMuMu_File = new TFile(OLD_DYToMuMu.c_str());
	TTree* OLD_DYToMuMu_miniTree = (TTree*) OLD_DYToMuMu_File->Get("miniTree");
	TTree* OLD_DYToMuMu_miniTree_allmuons = (TTree*) OLD_DYToMuMu_File->Get("miniTree_allmuons");
	TTree* OLD_DYToMuMu_miniTree_allphotons = (TTree*) OLD_DYToMuMu_File->Get("miniTree_allphotons");
	TFile *OLD_Data_File = new TFile(OLD_Data.c_str());
	TTree* OLD_Data_miniTree = (TTree*) OLD_Data_File->Get("miniTree");
	TTree* OLD_Data_miniTree_allmuons = (TTree*) OLD_Data_File->Get("miniTree_allmuons");
	TTree* OLD_Data_miniTree_allphotons = (TTree*) OLD_Data_File->Get("miniTree_allphotons");
//	TFile *WJetsToLNu_File = new TFile(WJetsToLNu.c_str());
//	TTree* WJetsToLNu_miniTree = (TTree*) WJetsToLNu_File->Get("miniTree");
//	TTree* WJetsToLNu_miniTree_allmuons = (TTree*) WJetsToLNu_File->Get("miniTree_allmuons");
//	TTree* WJetsToLNu_miniTree_allphotons = (TTree*) WJetsToLNu_File->Get("miniTree_allphotons");

//	TFile *QCDMu_File = new TFile(QCDMu.c_str());
//	TTree* QCDMu_miniTree = (TTree*) QCDMu_File->Get("miniTree");
//	TTree* QCDMu_miniTree_allmuons = (TTree*) QCDMu_File->Get("miniTree_allmuons");
//	TTree* QCDMu_miniTree_allphotons = (TTree*) QCDMu_File->Get("miniTree_allphotons");


	TCanvas *c1 = new TCanvas("Default", "Default");

//////	DrawDataDataMCMCplot(Data_miniTree_allmuons, DYToMuMu_miniTree_allmuons, OLD_DYToMuMu_miniTree_allmuons, QCDMu_miniTree_allmuons, OLD_Data_miniTree_allmuons, WJetsToLNu_miniTree_allmuons, "Ptmumu", "Ptmumu", "(100,0,200)", "isMM", "dimuon", "p_{T}^{#mu#mu} [GeV]", true, false, c1);

	vector<string> set_of_cuts;
	vector<string> name;

//	set_of_cuts.push_back("isMMGCandidate");
//	name.push_back("selected-00-beforeFSRcuts");
/*
	set_of_cuts.push_back("isAfterFSRCut1");
	name.push_back("selected-01");
	set_of_cuts.push_back("isAfterFSRCut2");
	name.push_back("selected-02");
	set_of_cuts.push_back("isAfterFSRCut3");
	name.push_back("selected-03");
	set_of_cuts.push_back("isAfterFSRCut4");
	name.push_back("selected-04");
	set_of_cuts.push_back("isVeryLooseMMG");
	name.push_back("selected-veryloose");
*/
	set_of_cuts.push_back("isLooseMMG");
  name.push_back("selected-loose");
/*
	set_of_cuts.push_back("isLooseMMG && Photon_isEB");
  name.push_back("selected-loose-EB");
	set_of_cuts.push_back("isLooseMMG && Photon_isEE");
  name.push_back("selected-loose-EE");
	set_of_cuts.push_back("isTightMMG");
  name.push_back("selected-tight");
	set_of_cuts.push_back("isTightMMG && Photon_isEB");
  name.push_back("selected-tight-EB");
	set_of_cuts.push_back("isTightMMG && Photon_isEE");
  name.push_back("selected-tight-EE");

	set_of_cuts.push_back("isLooseMMG && isMultipleCandidate==0");
  name.push_back("selected-loose-nomultiple");
	set_of_cuts.push_back("isLooseMMG && Photon_isEB && isMultipleCandidate==0");
  name.push_back("selected-loose-EB-nomultiple");
	set_of_cuts.push_back("isLooseMMG && Photon_isEE && isMultipleCandidate==0");
  name.push_back("selected-loose-EE-nomultiple");
	set_of_cuts.push_back("isTightMMG && isMultipleCandidate==0");
  name.push_back("selected-tight-nomultiple");
	set_of_cuts.push_back("isTightMMG && Photon_isEB && isMultipleCandidate==0");
  name.push_back("selected-tight-EB-nomultiple");
	set_of_cuts.push_back("isTightMMG && Photon_isEE && isMultipleCandidate==0");
  name.push_back("selected-tight-EE-nomultiple");
*/

	for(int i=0; i<set_of_cuts.size() ; i++)
	{
//	DrawDataDataMCMCplot(Data_miniTree, DYToMuMu_miniTree, OLD_Data_miniTree, OLD_DYToMuMu_miniTree, "Ptmumu", "Ptmumu", "(100,0,200)", set_of_cuts[i], name[i], "p_{T}^{#mu#mu} [GeV]", true, false, c1);
/*
//		DrawDataDataMCMCplot_TH1I(Data_miniTree, DYToMuMu_miniTree, OLD_Data_miniTree, OLD_DYToMuMu_miniTree, OLD_Data_miniTree, WJetsToLNu_miniTree, "NbMuons", "NbMuons", "(10,0,10)", set_of_cuts[i], name[i], "# of muons", true, false, c1);
//		DrawDataDataMCMCplot_TH1I(Data_miniTree, DYToMuMu_miniTree, OLD_Data_miniTree, OLD_DYToMuMu_miniTree, OLD_Data_miniTree, WJetsToLNu_miniTree, "NbPhotons", "NbPhotons", "(10,0,10)", set_of_cuts[i], name[i], "# of photons", true, false, c1);

*/
/*

//		DrawDataDataMCMCplot(Data_miniTree, DYToMuMu_miniTree, OLD_Data_miniTree, OLD_DYToMuMu_miniTree, OLD_Data_miniTree, WJetsToLNu_miniTree, "MuonL_Pt", "MuonL_Pt", "(50,0,100)", set_of_cuts[i], name[i], "p_{T} #mu_{leading} [GeV]", true, false, c1);
//		DrawDataDataMCMCplot(Data_miniTree, DYToMuMu_miniTree, OLD_Data_miniTree, OLD_DYToMuMu_miniTree, OLD_Data_miniTree, WJetsToLNu_miniTree, "MuonS_Pt", "MuonS_Pt", "(50,0,100)", set_of_cuts[i], name[i], "p_{T} #mu_{trailing} [GeV]", true, false, c1);
//		DrawDataDataMCMCplot(Data_miniTree, DYToMuMu_miniTree, OLD_Data_miniTree, OLD_DYToMuMu_miniTree, OLD_Data_miniTree, WJetsToLNu_miniTree, "MuonL_Eta", "MuonL_Eta", "(50,-3,3)", set_of_cuts[i], name[i], "#eta #mu_{leading}", true, false, c1);
//		DrawDataDataMCMCplot(Data_miniTree, DYToMuMu_miniTree, OLD_Data_miniTree, OLD_DYToMuMu_miniTree, OLD_Data_miniTree, WJetsToLNu_miniTree, "MuonS_Eta", "MuonS_Eta", "(50,-3,3)", set_of_cuts[i], name[i], "#eta #mu_{trailing}", true, false, c1);

//		DrawDataDataMCMCplot(Data_miniTree, DYToMuMu_miniTree, OLD_Data_miniTree, OLD_DYToMuMu_miniTree, OLD_Data_miniTree, WJetsToLNu_miniTree, "MuonL_Phi", "MuonL_Phi", "(50,-3.15,3.15)", set_of_cuts[i], name[i], "#phi #mu_{leading}", true, false, c1);
//		DrawDataDataMCMCplot(Data_miniTree, DYToMuMu_miniTree, OLD_Data_miniTree, OLD_DYToMuMu_miniTree, OLD_Data_miniTree, WJetsToLNu_miniTree, "MuonS_Phi", "MuonS_Phi", "(50,-3.15,3.15)", set_of_cuts[i], name[i], "#phi #mu_{trailing}", true, false, c1);
//		DrawDataDataMCMCplot(Data_miniTree, DYToMuMu_miniTree, OLD_Data_miniTree, OLD_DYToMuMu_miniTree, OLD_Data_miniTree, WJetsToLNu_miniTree, "MuonL_Charge", "MuonL_Charge", "(20,-2,2)", set_of_cuts[i], name[i], "charge #mu_{leading}", true, false, c1);
//		DrawDataDataMCMCplot(Data_miniTree, DYToMuMu_miniTree, OLD_Data_miniTree, OLD_DYToMuMu_miniTree, OLD_Data_miniTree, WJetsToLNu_miniTree, "MuonS_Charge", "MuonS_Charge", "(20,-2,2)", set_of_cuts[i], name[i], "charge #mu_{trailing}", true, false, c1);

//		DrawDataDataMCMCplot(Data_miniTree, DYToMuMu_miniTree, OLD_Data_miniTree, OLD_DYToMuMu_miniTree, OLD_Data_miniTree, WJetsToLNu_miniTree, "MuonF_Pt", "MuonF_Pt", "(50,0,100)", set_of_cuts[i], name[i], "p_{T} #mu_{far} [GeV]", true, false, c1);
//		DrawDataDataMCMCplot(Data_miniTree, DYToMuMu_miniTree, OLD_Data_miniTree, OLD_DYToMuMu_miniTree, OLD_Data_miniTree, WJetsToLNu_miniTree, "MuonN_Pt", "MuonN_Pt", "(50,0,100)", set_of_cuts[i], name[i], "p_{T} #mu_{near} [GeV]", true, false, c1);
//		DrawDataDataMCMCplot(Data_miniTree, DYToMuMu_miniTree, OLD_Data_miniTree, OLD_DYToMuMu_miniTree, OLD_Data_miniTree, WJetsToLNu_miniTree, "MuonF_Eta", "MuonF_Eta", "(50,-3,3)", set_of_cuts[i], name[i], "#eta #mu_{far}", true, false, c1);
//		DrawDataDataMCMCplot(Data_miniTree, DYToMuMu_miniTree, OLD_Data_miniTree, OLD_DYToMuMu_miniTree, OLD_Data_miniTree, WJetsToLNu_miniTree, "MuonN_Eta", "MuonN_Eta", "(50,-3,3)", set_of_cuts[i], name[i], "#eta #mu_{near}", true, false, c1);

		DrawDataDataMCMCplot(Data_miniTree, DYToMuMu_miniTree, OLD_Data_miniTree, OLD_DYToMuMu_miniTree, "MuonF_Phi", "MuonF_Phi", "(50,-3.15,3.15)", set_of_cuts[i], name[i], "#phi #mu_{far}", true, false, c1);
		DrawDataDataMCMCplot(Data_miniTree, DYToMuMu_miniTree, OLD_Data_miniTree, OLD_DYToMuMu_miniTree, "MuonN_Phi", "MuonN_Phi", "(50,-3.15,3.15)", set_of_cuts[i], name[i], "#phi #mu_{near}", true, false, c1);
		DrawDataDataMCMCplot(Data_miniTree, DYToMuMu_miniTree, OLD_Data_miniTree, OLD_DYToMuMu_miniTree, "MuonF_Charge", "MuonF_Charge", "(20,-2,2)", set_of_cuts[i], name[i], "charge #mu_{far}", true, false, c1);
		DrawDataDataMCMCplot(Data_miniTree, DYToMuMu_miniTree, OLD_Data_miniTree, OLD_DYToMuMu_miniTree, "MuonN_Charge", "MuonN_Charge", "(20,-2,2)", set_of_cuts[i], name[i], "charge #mu_{near}", true, false, c1);
*/
//		DrawDataDataMCMCplot(Data_miniTree, DYToMuMu_miniTree, OLD_Data_miniTree, OLD_DYToMuMu_miniTree, "Mmumu", "Mmumu", "(150,0,300)", set_of_cuts[i], name[i], "m_{#mu#mu} [GeV]", true, false, c1);
//		DrawDataDataMCMCplot(Data_miniTree, DYToMuMu_miniTree, OLD_Data_miniTree, OLD_DYToMuMu_miniTree, "Mmumu", "Mmumu_extended", "(30,30,90)", set_of_cuts[i], name[i], "m_{#mu#mu} [GeV]", true, false, c1);
		DrawDataDataMCMCplot(Data_miniTree, DYToMuMu_miniTree, OLD_Data_miniTree, OLD_DYToMuMu_miniTree, "Mmumugamma", "Mmumugamma", "(300,0,300)", set_of_cuts[i], name[i], "m_{#mu#mu#gamma} [GeV]", true, false, c1);
		DrawDataDataMCMCplot(Data_miniTree, DYToMuMu_miniTree, OLD_Data_miniTree, OLD_DYToMuMu_miniTree, "Mmumugamma", "Mmumugamma_extended", "(60,60,120)", set_of_cuts[i], name[i], "m_{#mu#mu#gamma} [GeV]", true, false, c1);
		DrawDataDataMCMCplot(Data_miniTree, DYToMuMu_miniTree, OLD_Data_miniTree, OLD_DYToMuMu_miniTree, "Mmumugamma", "Mmumugamma_zoom", "(48,85,97)", set_of_cuts[i], name[i], "m_{#mu#mu#gamma} [GeV]", true, false, c1);

		DrawDataDataMCMCplot(Data_miniTree, DYToMuMu_miniTree, OLD_Data_miniTree, OLD_DYToMuMu_miniTree, "Mmumugamma_SCraw", "Mmumugamma_SCraw", "(50,0,200)", set_of_cuts[i], name[i], "m_{#mu#mu#gamma}^{SCraw} [GeV]", true, false, c1);
//		DrawDataDataMCMCplot(Data_miniTree, DYToMuMu_miniTree, OLD_Data_miniTree, OLD_DYToMuMu_miniTree, "Mmumugamma_SCraw_fEta", "Mmumugamma_SCraw_fEta", "(50,0,400)", set_of_cuts[i], name[i], "m_{#mu#mu#gamma}^{SCraw x fEta} [GeV]", true, false, c1);
	

		DrawDataDataMCMCplot(Data_miniTree, DYToMuMu_miniTree, OLD_Data_miniTree, OLD_DYToMuMu_miniTree, "deltaRNear", "deltaRNear", "(100,0,1)", set_of_cuts[i], name[i], "#Delta R(#gamma, #mu_{near})", true, false, c1);
		DrawDataDataMCMCplot(Data_miniTree, DYToMuMu_miniTree, OLD_Data_miniTree, OLD_DYToMuMu_miniTree, "deltaRFar", "deltaRFar", "(100,0,5)", set_of_cuts[i], name[i], "#Delta R(#gamma, #mu_{far})", true, false, c1);
//		DrawDataDataMCMCplot(Data_miniTree, DYToMuMu_miniTree, OLD_Data_miniTree, OLD_DYToMuMu_miniTree, "deltaRPlus", "deltaRPlus", "(100,0,10)", set_of_cuts[i], name[i], "#Delta R(#gamma, #mu_{plus})", true, false, c1);
//		DrawDataDataMCMCplot(Data_miniTree, DYToMuMu_miniTree, OLD_Data_miniTree, OLD_DYToMuMu_miniTree, "deltaRMinus", "deltaRMinus", "(100,0,10)", set_of_cuts[i], name[i], "#Delta R(#gamma, #mu_{minus})", true, false, c1);
//		DrawDataDataMCMCplot(Data_miniTree, DYToMuMu_miniTree, OLD_Data_miniTree, OLD_DYToMuMu_miniTree, "deltaRLeading", "deltaRLeading", "(100,0,10)", set_of_cuts[i], name[i], "#Delta R(#gamma, #mu_{leading})", true, false, c1);
//		DrawDataDataMCMCplot(Data_miniTree, DYToMuMu_miniTree, OLD_Data_miniTree, OLD_DYToMuMu_miniTree, "deltaRSubleading", "deltaRSubleading", "(100,0,10)", set_of_cuts[i], name[i], "#Delta R(#gamma, #mu_{trailing})", true, false, c1);


//		DrawDataDataMCMCplot(Data_miniTree, DYToMuMu_miniTree, OLD_Data_miniTree, OLD_DYToMuMu_miniTree, "Photon_Eta", "Photon_Eta", "(50,-3,3)", set_of_cuts[i], name[i], "#eta^{#gamma}", true, false, c1);
		DrawDataDataMCMCplot(Data_miniTree, DYToMuMu_miniTree, OLD_Data_miniTree, OLD_DYToMuMu_miniTree, "Photon_Eta", "Photon_Eta", "(16,-3.2,3.2)", set_of_cuts[i], name[i], "#eta^{#gamma}", true, false, c1);

		DrawDataDataMCMCplot(Data_miniTree, DYToMuMu_miniTree, OLD_Data_miniTree, OLD_DYToMuMu_miniTree, "Photon_Phi", "Photon_Phi", "(21,-3.15,3.15)", set_of_cuts[i], name[i], "#phi^{#gamma}", true, false, c1);

/*
*/
//		DrawDataDataMCMCplot(Data_miniTree, DYToMuMu_miniTree, OLD_Data_miniTree, OLD_DYToMuMu_miniTree, "Photon_E", "Photon_E", "(100, 0, 100)", set_of_cuts[i], name[i], "E^{#gamma} [GeV]", true, false, c1);
//		DrawDataDataMCMCplot(Data_miniTree, DYToMuMu_miniTree, OLD_Data_miniTree, OLD_DYToMuMu_miniTree, "Photon_E", "Photon_E_extended", "(100, 0, 150)", set_of_cuts[i], name[i], "E^{#gamma} [GeV]", true, false, c1);
		DrawDataDataMCMCplot(Data_miniTree, DYToMuMu_miniTree, OLD_Data_miniTree, OLD_DYToMuMu_miniTree, "Photon_SC_rawE", "Photon_SC_rawE", "(50, 0, 400)", set_of_cuts[i], name[i], "E^{SC raw} [GeV]", true, false, c1);
		DrawDataDataMCMCplot(Data_miniTree, DYToMuMu_miniTree, OLD_Data_miniTree, OLD_DYToMuMu_miniTree, "Photon_SC_rawEt", "Photon_SC_rawEt", "(50, 0, 400)", set_of_cuts[i], name[i], "E^{SC raw}_{T} [GeV]", true, false, c1);

		DrawDataDataMCMCplot(Data_miniTree, DYToMuMu_miniTree, OLD_Data_miniTree, OLD_DYToMuMu_miniTree, "Photon_E", "Photon_E", "(50, 0, 250)", set_of_cuts[i], name[i], "E^{#gamma} [GeV]", true, false, c1);

//		DrawDataDataMCMCplot(Data_miniTree, DYToMuMu_miniTree, OLD_Data_miniTree, OLD_DYToMuMu_miniTree, "Photon_Et", "Photon_Et", "(100, 0, 100)", set_of_cuts[i], name[i], "E_{T}^{#gamma} [GeV]", true, false, c1);
		DrawDataDataMCMCplot(Data_miniTree, DYToMuMu_miniTree, OLD_Data_miniTree, OLD_DYToMuMu_miniTree, "Photon_Et", "Photon_Et", "(50, 0, 100)", set_of_cuts[i], name[i], "E_{T}^{#gamma} [GeV]", true, false, c1);
/*
		DrawDataDataMCMCplot_TH2F(Data_miniTree, DYToMuMu_miniTree, OLD_Data_miniTree, OLD_DYToMuMu_miniTree, "Mmumu", "Mmumugamma", "(900,0,300,900,0,300)", set_of_cuts[i], name[i], "m_{#mu#mu} [GeV]", "m_{#mu#mu#gamma} [GeV]", "Mmumu_VS_Mmumugamma", false, false, c1);
		DrawDataDataMCMCplot_TH2F(Data_miniTree, DYToMuMu_miniTree, OLD_Data_miniTree, OLD_DYToMuMu_miniTree, "Mmumu", "Mmumugamma", "(450,20,100,450,80,110)", set_of_cuts[i], name[i], "m_{#mu#mu} [GeV]", "m_{#mu#mu#gamma} [GeV]", "Mmumu_VS_Mmumugamma_extended", false, false, c1);

DrawDataDataMCMCplot_TH2F(Data_miniTree, DYToMuMu_miniTree, OLD_Data_miniTree, OLD_DYToMuMu_miniTree, "((91.1876**2 - Mmumu**2 ) / (Mmumugamma**2 - Mmumu**2))*Photon_E", "Photon_E", "(400,0,100,400,0,100)", set_of_cuts[i], name[i], "E_{true} = k*E_{reco} [GeV]", "E_{reco} [GeV]", "Etrue_VS_Ereco", false, false, c1, true);
DrawDataDataMCMCplot_TH2F(Data_miniTree, DYToMuMu_miniTree, OLD_Data_miniTree, OLD_DYToMuMu_miniTree, "((91.1876**2 - Mmumu**2 ) / (Mmumugamma**2 - Mmumu**2))*Photon_E", "Photon_E", "(600,0,150,600,0,150)", set_of_cuts[i], name[i], "E_{true} = k*E_{reco} [GeV]", "E_{reco} [GeV]", "Etrue_VS_Ereco_extended", false, false, c1, true);
*/


	DrawDataDataMCMCplot(Data_miniTree, DYToMuMu_miniTree, OLD_Data_miniTree, OLD_DYToMuMu_miniTree, "mmg_k", "mmg_k", "(40,0,2.0)", set_of_cuts[i], name[i], "k = E_{muons} / E_{reco}", true, false, c1);
	DrawDataDataMCMCplot(Data_miniTree, DYToMuMu_miniTree, OLD_Data_miniTree, OLD_DYToMuMu_miniTree, "mmg_s", "mmg_s", "(40,-1.0,1.0)", set_of_cuts[i], name[i], "s = E_{reco} / E_{muons} - 1", true, false, c1);

	DrawDataDataMCMCplot(Data_miniTree, DYToMuMu_miniTree, OLD_Data_miniTree, OLD_DYToMuMu_miniTree, "mmg_k", "mmg_k_extended", "(20,0,2.0)", set_of_cuts[i], name[i], "k = E_{muons} / E_{reco}", true, false, c1);
	DrawDataDataMCMCplot(Data_miniTree, DYToMuMu_miniTree, OLD_Data_miniTree, OLD_DYToMuMu_miniTree, "mmg_s", "mmg_s_extended", "(20,-1.0,1.0)", set_of_cuts[i], name[i], "s = E_{reco} / E_{muons} - 1", true, false, c1);
	DrawDataDataMCMCplot(Data_miniTree, DYToMuMu_miniTree, OLD_Data_miniTree, OLD_DYToMuMu_miniTree, "mmg_k", "mmg_k_zoom", "(20,0.6,1.6)", set_of_cuts[i], name[i], "k = E_{muons} / E_{reco}", true, false, c1);
	DrawDataDataMCMCplot(Data_miniTree, DYToMuMu_miniTree, OLD_Data_miniTree, OLD_DYToMuMu_miniTree, "mmg_s", "mmg_s_zoom", "(20,-0.5,0.5)", set_of_cuts[i], name[i], "s = E_{reco} / E_{muons} - 1", true, false, c1);

	DrawDataDataMCMCplot(Data_miniTree, DYToMuMu_miniTree, OLD_Data_miniTree, OLD_DYToMuMu_miniTree, "mmg_k_SCraw", "mmg_k_SCraw", "(60,-1.0,3.0)", set_of_cuts[i], name[i], "k_{SCraw} = E_{muons} / E_{SCraw}", true, false, c1);
	DrawDataDataMCMCplot(Data_miniTree, DYToMuMu_miniTree, OLD_Data_miniTree, OLD_DYToMuMu_miniTree, "mmg_s_SCraw", "mmg_s_SCraw", "(60,-2.0,2.0)", set_of_cuts[i], name[i], "s_{SCraw} = E_{SCraw} / E_{muons} - 1", true, false, c1);

//	DrawDataDataMCMCplot(Data_miniTree, DYToMuMu_miniTree, OLD_Data_miniTree, OLD_DYToMuMu_miniTree, "mmg_k_SCraw_fEta", "mmg_k_SCraw_fEta", "(60,-1.0,3.0)", set_of_cuts[i], name[i], "k_{SCraw x fEta} = E_{muons} / E_{SCraw x fEta}", true, false, c1);
//	DrawDataDataMCMCplot(Data_miniTree, DYToMuMu_miniTree, OLD_Data_miniTree, OLD_DYToMuMu_miniTree, "mmg_s_SCraw_fEta", "mmg_s_SCraw_fEta", "(60,-2.0,2.0)", set_of_cuts[i], name[i], "s_{SCraw x fEta} = E_{SCraw x fEta} /  E_{muons} -1", true, false, c1);

/*
//DrawDataDataMCMCplot_TH2F(Data_miniTree, DYToMuMu_miniTree, OLD_Data_miniTree, OLD_DYToMuMu_miniTree, "Photon_convNTracks", "MuonN_isoR03_sumPt", "(16,0,4,200,0,10)", set_of_cuts[i], name[i], "Photon_convNTracks", "MuonN_isoR03_sumPt", "Photon_convNTracks_VS_MuonN_isoR03_sumPt", false, false, c1, false);
//DrawDataDataMCMCplot_TH2F(Data_miniTree, DYToMuMu_miniTree, OLD_Data_miniTree, OLD_DYToMuMu_miniTree, "Photon_convNTracks", "MuonN_isoR03_sumPt", "(12,0,3,200,0,3)", set_of_cuts[i], name[i], "Photon_convNTracks", "MuonN_isoR03_sumPt", "Photon_convNTracks_VS_MuonN_isoR03_sumPt_extended", false, false, c1, false);
	//	DrawDataDataMCMCplot(Data_miniTree, DYToMuMu_miniTree, OLD_Data_miniTree, OLD_DYToMuMu_miniTree, "MuonN_isoR03_sumPt", "MuonN_isoR03_sumPt", "(100,0,10)", set_of_cuts[i], name[i], "MuonN_isoR03_sumPt", true, false, c1);
*/
	}
	return 0;	
}

