#include "TChain.h"
#include "TFile.h"
#include <iostream>
using namespace std;


int hltnamesTest()
{
	TChain * chain = new TChain("miniTree");

	chain->Add("miniTree_totouples_Run2012A_22Jan2013_v1_November2013_0_injRe0_v1_partALL.root");
        chain->Add("miniTree_totouples_parked_Run2012B_22Jan2013_v1_November2013_0_injRe0_v1_partALL.root");
        chain->Add("miniTree_totouples_parked_Run2012C_22Jan2013_v1_November2013_0_injRe0_v1_partALL.root");
        chain->Add("miniTree_totouples_parked_Run2012D_22Jan2013_v1_November2013_0_injRe0_v1_partALL.root");
        chain->Add("miniTree_totouples_DYToMuMu_Summer12_November2013_1_injRe0_v1_partALL.root");
        chain->Add("miniTree_totouples_DYToMuMu_Summer12_November2013_2_injRe0_v1_partALL.root");
        chain->Add("miniTree_totouples_TTJets_SemiLeptMGDecays_Summer12_November2013_3_injRe0_v1_partALL.root");
        chain->Add("miniTree_totouples_TTJets_FullLeptMGDecays_Summer12_November2013_v2_3_injRe0_v1_partALL.root");
        chain->Add("miniTree_totouples_WJetsToLNu_Summer12_S10_November2013_3_injRe0_v1_partALL.root");
	

	string cut = "";

	for(int i = 0; i < 25; i++)
	{
		cut = Form(" hltnames == \"HLT_Mu17_TkMu8_v%d\" ",i);

		cout << cut << " : " << chain->GetEntries(cut.c_str()) << endl;
	}

	return 0;
}

