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
#include "TGaxis"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstring>
#include <iomanip>
#pragma optimize 0
using namespace std;

int EventsList(string inputMiniTree = "/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/cvs_developpment/Selection_NewMuID/miniTree_partALL")
{
	TChain* plaf = new TChain("miniTree");
	//plaf->Add(Form("%s.root",inputMiniTree.c_str()));


	//plaf->Add("/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/cvs_developpment/Selection_NewMuID/miniTree_Run2012A_13Jul2012_v1_NewMuonID_NewSelection_0_injRe0.00_v1_partALL.root");
	//plaf->Add("/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/cvs_developpment/Selection_NewMuID/miniTree_Run2012A_recover_06Aug2012_v1_NewMuonID_NewSelection_0_injRe0.00_v1_partALL.root");
	//plaf->Add("/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/cvs_developpment/Selection_NewMuID/miniTree_Run2012B_13Jul2012_v4_NewMuonID_NewSelection_0_injRe0.00_v1_partALL.root");
	//plaf->Add("/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/cvs_developpment/Selection_NewMuID/miniTree_Run2012C-24Aug2012-v1_NewMuonID_NewSelection_0_injRe0.00_v1_partALL.root");
	//plaf->Add("/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/cvs_developpment/Selection_NewMuID/miniTree_Run2012C-EcalRecover_11Dec2012-v1_NewMuonID_NewSelection_0_injRe0.00_v1_partALL.root");
	//plaf->Add("/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/cvs_developpment/Selection_NewMuID/miniTree_Run2012C_PromptReco_v2_NewMuonID_NewSelection_0_injRe0.00_v1_partALL.root");
	//plaf->Add("/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/cvs_developpment/Selection_NewMuID/miniTree_Run2012D_PromptReco_v1_NewMuonID_NewSelection_0_injRe0.00_v1_partALL.root");

	//plaf->Add("/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/cvs_developpment/Selection_NewMuID/miniTree_DYToMuMu_Summer12_NewMuonID_NewSelection_1_injRe0.00_v1_partALL.root");
	//plaf->Add("/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/cvs_developpment/Selection_NewMuID/miniTree_DYToMuMu_Summer12_NewMuonID_NewSelection_2_injRe0.00_v1_partALL.root");
	//plaf->Add("/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/cvs_developpment/Selection_NewMuID/miniTree_TTJets_Summer12_S7_NewMuonID_NewSelection_3_injRe0.00_v1_partALL.root");
	//plaf->Add("/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/cvs_developpment/Selection_NewMuID/miniTree_WJetsToLNu_Summer12_S10_NewMuonID_NewSelection_3_injRe0.00_v1_partALL.root");
	plaf->Add("miniTree_parked_Run2012A_22Jan2013_v1_3_NewSelection_0_injRe0_v1_22Jan2013_partALL.root");

	plaf->SetScanField(0);
	//plaf->Scan("iRunID:iLumiID:iEventID:Mmumugamma:Mmumu:Photon_E:Photon_SC_Eta:Photon_r9", "isJanLooseMMG == 1 && (Photon_isEB == 1 || Photon_isEE == 1) && Photon_Et > 25", "colsize=50"); 
	TString cut = "isJanLooseMMG == 1";
	//cut += " && ( hltnames == \"HLT_Mu17_TkMu8_v9\" || hltnames == \"HLT_Mu17_TkMu8_v10\" || hltnames == \"HLT_Mu17_TkMu8_v11\" || hltnames == \"HLT_Mu17_TkMu8_v12\" || hltnames == \"HLT_Mu17_TkMu8_v13\" || hltnames == \"HLT_Mu17_TkMu8_v14\" || hltnames == \"HLT_Mu17_TkMu8_v15\" || hltnames == \"HLT_Mu17_TkMu8_v16\" )";
        cut += " && Photon_Et > 25";

	cout<<endl<<"cut = "<<cut<<endl;

	plaf->Scan("iRunID:iLumiID:iEventID", cut, "colsize=50");	

	return 0;
}

