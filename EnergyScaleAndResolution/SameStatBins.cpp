// -------------------------------------------------------
// program implemented by Louis Sgandurra (October 2012)
// -------------------------------------------------------

#include "functions.h"

int main(int argc, char *argv[])
{
        for(int iarg = 0 ; iarg < argc; iarg++)
        {
                cout << "argv[" << iarg << "]= " << argv[iarg] << endl;
        }

        if( argc == 1 ) 
        {
                cerr << "arguments should be passed : dataType, eta, r9, variable, nBins" <<endl; 
                return 1;

        }

	string dataType = "data";
        string eta = "all"; 
        string r9 = "all";
	string variable = "nVertices";
	int nBins = 5;

	if( argc > 1 ) dataType = argv[1];
        if( argc > 2 ) eta = argv[2];
        if( argc > 3 ) r9 = argv[3];
        if( argc > 4 ) variable = argv[4];
        if( argc > 5 ) 
        {
                std::stringstream ss ( argv[5] );
                ss >> nBins;
        }

	int injectedResolutionTemp = 0;
	TChain * chain = new TChain("miniTree");
	
	if(dataType == "data")
        {
                //chain->Add("../cvs_developpment/Selection_NewMuID/miniTree_Run2012A_13Jul2012_v1_NewMuonID_NewSelection_0_v2_partALL.root");
                //chain->Add("../cvs_developpment/Selection_NewMuID/miniTree_Run2012A_recover_06Aug2012_v1_NewMuonID_NewSelection_0_v2_partALL.root");
                //chain->Add("../cvs_developpment/Selection_NewMuID/miniTree_Run2012B_13Jul2012_v4_NewMuonID_NewSelection_0_v2_partALL.root");
                //chain->Add("../cvs_developpment/Selection_NewMuID/miniTree_Run2012C-24Aug2012-v1_NewMuonID_NewSelection_0_v2_partALL.root");
                //chain->Add("../cvs_developpment/Selection_NewMuID/miniTree_Run2012C-EcalRecover_11Dec2012-v1_NewMuonID_NewSelection_0_v2_partALL.root");                
		//chain->Add("../cvs_developpment/Selection_NewMuID/miniTree_Run2012C_PromptReco_v2_NewMuonID_NewSelection_0_v2_partALL.root");
                //chain->Add("../cvs_developpment/Selection_NewMuID/miniTree_Run2012D_PromptReco_v1_NewMuonID_NewSelection_0_v2_partALL.root");
                //chain->Add("../cvs_developpment/Selection_NewMuID/miniTree_2011A_03Oct2011V1ReReco_toto_v2_NewSelection_0_v6_partALL.root");
                //chain->Add("../cvs_developpment/Selection_NewMuID/miniTree_2011A_05Jul2011ReReco_toto_v2_NewSelection_0_v6_partALL.root");
                //chain->Add("../cvs_developpment/Selection_NewMuID/miniTree_2011A_PromptSkimV5ReReco_toto_v2_NewSelection_0_v6_partALL.root");
                //chain->Add("../cvs_developpment/Selection_NewMuID/miniTree_2011B_PromptSkimV1ReReco_toto_v2_NewSelection_0_v6_partALL.root");
        
		chain->Add(Form("../cvs_developpment/Selection_NewMuID/miniTree_Run2012A_13Jul2012_v1_NewMuonID_NewSelection_0_injRe%d_v3_NewPU_partALL.root",injectedResolutionTemp));
		chain->Add(Form("../cvs_developpment/Selection_NewMuID/miniTree_Run2012A_recover_06Aug2012_v1_NewMuonID_NewSelection_0_injRe%d_v3_NewPU_partALL.root",injectedResolutionTemp));
                chain->Add(Form("../cvs_developpment/Selection_NewMuID/miniTree_Run2012B_13Jul2012_v4_NewMuonID_NewSelection_0_injRe%d_v3_NewPU_partALL.root",injectedResolutionTemp));
                chain->Add(Form("../cvs_developpment/Selection_NewMuID/miniTree_Run2012C-24Aug2012-v1_NewMuonID_NewSelection_0_injRe%d_v3_NewPU_partALL.root",injectedResolutionTemp));
                chain->Add(Form("../cvs_developpment/Selection_NewMuID/miniTree_Run2012C-EcalRecover_11Dec2012-v1_NewMuonID_NewSelection_0_injRe%d_v3_NewPU_partALL.root",injectedResolutionTemp));      
                chain->Add(Form("../cvs_developpment/Selection_NewMuID/miniTree_Run2012C_PromptReco_v2_NewMuonID_NewSelection_0_injRe%d_v3_NewPU_partALL.root",injectedResolutionTemp));
                chain->Add(Form("../cvs_developpment/Selection_NewMuID/miniTree_Run2012D_PromptReco_v1_NewMuonID_NewSelection_0_injRe%d_v3_NewPU_partALL.root",injectedResolutionTemp));

	}
        if(dataType == "MC")
        {
                chain->Add("../cvs_developpment/Selection_NewMuID/miniTree_DYToMuMu_Summer12_NewMuonID_NewSelection_1_v2_partALL.root");
                chain->Add("../cvs_developpment/Selection_NewMuID/miniTree_DYToMuMu_Summer12_NewMuonID_NewSelection_2_v2_partALL.root");
                //chain->Add("../cvs_developpment/Selection_NewMuID/miniTree_TTJets_Summer12_S7_NewMuonID_NewSelection_3_v2_partALL.root");
                //chain->Add("../cvs_developpment/Selection_NewMuID/miniTree_WJetsToLNu_Summer12_S10_NewMuonID_NewSelection_3_v2_partALL.root");
                //chain->Add("../cvs_developpment/Selection_NewMuID/miniTree_DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Fall11-PU_S6_START42_V14B-v1_September12_NewSelection_1_v6_partALL.root");
                //chain->Add("../cvs_developpment/Selection_NewMuID/miniTree_DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Fall11-PU_S6_START42_V14B-v1_September12_NewSelection_2_v6_partALL.root");
                //chain->Add("../cvs_developpment/Selection_NewMuID/miniTree_TTJets_TuneZ2_7TeV-madgraph-tauola_NewSelection_3_v1_partALL.root");
                //chain->Add("../cvs_developpment/Selection_NewMuID/miniTree_WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_NewSelection_3_v1_partALL.root");
        }

	TString cut = "isJanLooseMMG == 1 && Photon_Et > 25";	

       	if(r9 == "low" && eta == "Barrel") cut += " && Photon_isEB == 1 && Photon_r9 < 0.94";
        if(r9 == "high" && eta == "Barrel") cut += " && Photon_isEB == 1 && Photon_r9 > 0.94";
        if(r9 == "low" && eta == "Endcaps") cut += " && Photon_isEE == 1 && Photon_r9 < 0.95";
        if(r9 == "high" && eta == "Endcaps") cut += " && Photon_isEE == 1 && Photon_r9 > 0.95";
        if(r9 == "all" && eta == "Barrel") cut += " && Photon_isEB == 1";
        if(r9 == "all" && eta == "Endcaps") cut += " && Photon_isEE == 1";
	if(r9 == "all" && eta == "all") cut += " && (Photon_isEB == 1 || Photon_isEE == 1)";	

	cout<<endl<<"cut = "<<cut<<endl;
	cout<<endl<<"chain->GetEntries(cut) = "<<chain->GetEntries(cut)<<endl;

	TChain * ReducedChain = (TChain *) chain->CopyTree(cut);

	int nVertices;
        float Photon_Et;
	float Photon_E;
	float Photon_SC_brem;
	float Photon_SC_Eta;
	float Photon_SC_rawEt;
	
	ReducedChain->SetBranchAddress("nVertices",&nVertices);
	ReducedChain->SetBranchAddress("Photon_Et",&Photon_Et);
	ReducedChain->SetBranchAddress("Photon_E",&Photon_E);
	ReducedChain->SetBranchAddress("Photon_SC_brem",&Photon_SC_brem);
	ReducedChain->SetBranchAddress("Photon_SC_Eta",&Photon_SC_Eta);
	ReducedChain->SetBranchAddress("Photon_SC_rawEt",&Photon_SC_rawEt);

	vector <double> variableTab;

	cout<<endl<<"ReducedChain->GetEntries() = "<<ReducedChain->GetEntries()<<endl;
	if(variable == "nVertices")
	{
		for (int ievt = 0 ; ievt < ReducedChain->GetEntries() ; ievt++)
		{
			ReducedChain->GetEntry(ievt);
			variableTab.push_back(nVertices);
		}
	}
	if(variable == "Photon_Et")
        {
                for (int ievt = 0 ; ievt < ReducedChain->GetEntries() ; ievt++)
                {
			ReducedChain->GetEntry(ievt);
                        variableTab.push_back(Photon_Et);
                }
        }
	if(variable == "Photon_E")
        {
                for (int ievt = 0 ; ievt < ReducedChain->GetEntries() ; ievt++)
                {
			ReducedChain->GetEntry(ievt);
                        variableTab.push_back(Photon_E);
                }
        }
	if(variable == "Photon_SC_brem")
        {
                for (int ievt = 0 ; ievt < ReducedChain->GetEntries() ; ievt++)
                {
			ReducedChain->GetEntry(ievt);
                        variableTab.push_back(Photon_SC_brem);
                }
        }
	if(variable == "Photon_SC_Eta")
        {
                for (int ievt = 0 ; ievt < ReducedChain->GetEntries() ; ievt++)
                {
			ReducedChain->GetEntry(ievt);
                        variableTab.push_back(Photon_SC_Eta);
                }
        }
	if(variable == "Photon_SC_rawEt")
        {
                for (int ievt = 0 ; ievt < ReducedChain->GetEntries() ; ievt++)
                {
			ReducedChain->GetEntry(ievt);
                        variableTab.push_back(Photon_SC_rawEt);
                }
        }	


        sort(variableTab.begin(), variableTab.end());

	int entriesPerBins = variableTab.size() / nBins;

	cout<<endl<<"entriesPerBins = "<<entriesPerBins<<endl;

	ofstream file(Form("Limites_%s_%sR9_%s_%s_%dbins.txt",variable.c_str(),r9.c_str(),eta.c_str(),dataType.c_str(),nBins));

	file<<0.0<<endl;	

	for(int i = 1; i <= nBins; i++)
	{
		cout<<endl<<"i = "<<i;
		if(i != nBins) file<<variableTab[i*entriesPerBins -1]<<endl;
		if(i == nBins) file<<variableTab[variableTab.size() - 1];
	}

	
	file.close();	


	return 0;
}	













