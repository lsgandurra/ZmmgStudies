#include "TTree.h"
#include "TClonesArray.h"
#include "TBranch.h"
#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "TBits.h"
#include "TMath.h"
#include "TSystem.h"
#include "TCanvas.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstring>
#include <iomanip>
#pragma optimize 0
using namespace std;

int EventsList(string inputMiniTree = "miniTree_partALL")
{
	TString cut = "isJanLooseMMG == 1";
	cut += " && ( hltnames == \"HLT_Mu17_TkMu8_v9\" || hltnames == \"HLT_Mu17_TkMu8_v10\" || hltnames == \"HLT_Mu17_TkMu8_v11\" || hltnames == \"HLT_Mu17_TkMu8_v12\" || hltnames == \"HLT_Mu17_TkMu8_v13\" || hltnames == \"HLT_Mu17_TkMu8_v14\" || hltnames == \"HLT_Mu17_TkMu8_v15\" || hltnames == \"HLT_Mu17_TkMu8_v16\" )";
	//cut += " && Photon_Et > 25"; //FIXME

	//cout<<endl<<"cut = "<<cut<<endl;

	TChain* plaf = new TChain("miniTree");
	plaf->Add(Form("%s.root",inputMiniTree.c_str()));
	plaf->SetScanField(0);
	plaf->Scan("iRunID:iLumiID:iEventID", cut, "colsize=50"); 

	return 0;
}
