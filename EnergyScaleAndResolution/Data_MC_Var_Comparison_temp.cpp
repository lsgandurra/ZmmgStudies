#include "functions.h"
#include "setTDRStyle.C"
//#include "style-Egamma.C"

//int data_MC_Pt_Comparison_entries(int EndCaps = 0, int r9sup = 1, int log = 0)
int main(int argc, char *argv[])
{

	for(int iarg = 0 ; iarg < argc; iarg++)
        {
                cout << "argv[" << iarg << "]= " << argv[iarg] << endl;
        }

        if( argc == 1 ) 
        {
                cerr << "arguments should be passed : directoryName, eta, r9, xVariable, log" <<endl; 
                return 1;

        }    

	string directoryName = "testLog";
	string eta = "all";
	string r9 = "all";
	string xVariable = "Photon_Et";
	int log = 1;	
	string normalization = "integral"; // "integral"	
	string cutVariable = "Photon_Et";
	string cutVariableValue = "0";
	string isAfterFSRCut = "isAfterFSRCut1";

	if( argc > 1 ) directoryName = argv[1];
	if( argc > 2 ) eta = argv[2];
        if( argc > 3 ) r9 = argv[3];	
	if( argc > 4 ) xVariable = argv[4];
	if( argc > 5 ) 
        {
                std::stringstream ss ( argv[5] );
                ss >> log;
        }
	if( argc > 6 ) normalization = argv[6];
	if( argc > 7 ) cutVariable = argv[7];
	if( argc > 8 ) cutVariableValue = argv[8];	
	if( argc > 9 ) isAfterFSRCut = argv[9];


	int nBins = 30;
	double xMin, xMax;
	double xMinLeg, xMaxLeg;

	xMinLeg = 0.67;
	xMaxLeg = 0.91;	
	string xVariableName, yVariableName;

	if(xVariable == "Photon_Et")
	{
		nBins = 30;
		xMin = 0;
		xMax = 150;
		xVariableName = "E_{T}^{#gamma} [GeV]";
		yVariableName = "Events / 5 GeV";
	}

	if(xVariable == "Mmumugamma")
        {
                nBins = 60;
                xMin = 60;
                xMax = 120;
                xVariableName = "M_{#mu#mu#gamma} [GeV]";
                yVariableName = "Events / 1 GeV";
        }
	
	if(xVariable == "Photon_r9")
        {
                nBins = 50;
                xMin = 0;
                xMax = 1;
                xVariableName = "r9";
                yVariableName = "Events / 0.02";
		xMinLeg = 0.20;
		xMaxLeg = 0.44;
        }
	
	if(xVariable == "nVertices")
        {
                nBins = 50; 
                xMin = 0;
                xMax = 50;
                xVariableName = "nVertices";
                yVariableName = "Events / 1";
                //xMinLeg = 0.20;
                //xMaxLeg = 0.44;
        }	
	
	if(xVariable == "nGenVertices")
	{
		nBins = 50; 
                xMin = 0;
                xMax = 50; 
                xVariableName = "nGenVertices";
                yVariableName = "Events / 1";
                //xMinLeg = 0.20;
                //xMaxLeg = 0.44;	
	}
	
	if(xVariable == "Photon_SC_Eta")
        {
                nBins = 100; 
                xMin = -2.6;
                xMax = 2.6; 
                xVariableName = "#eta_{#gamma}";
                yVariableName = "Events / 0.052";
                //xMinLeg = 0.20;
                //xMaxLeg = 0.44;
        }	
	
	if(xVariable == "MuonM_Eta")
        {
                nBins = 100; 
                xMin = -2.6;
                xMax = 2.6; 
                xVariableName = "#eta_{#mu^{-}}";
                yVariableName = "Events / 0.052";
                //xMinLeg = 0.20;
                //xMaxLeg = 0.44;
        }
	
	if(xVariable == "MuonP_Eta")
        {
                nBins = 100;
                xMin = -2.6;
                xMax = 2.6;
                xVariableName = "#eta_{#mu^{+}}";
                yVariableName = "Events / 0.052";
                //xMinLeg = 0.20;
                //xMaxLeg = 0.44;
        }


	gROOT->Reset();
        TGaxis::SetMaxDigits(3);
        setTDRStyle();
	//setEgammaStyle();

	string fileName = directoryName;
	directoryName += Form("/%s_cut_%s/%s/%s/Normalization_%s",cutVariable.c_str(), cutVariableValue.c_str(), isAfterFSRCut.c_str(), xVariable.c_str(),normalization.c_str());
	string directoryName_2 = directoryName;

	if(eta == "Barrel") directoryName += "/Barrel_";
        if(eta == "Endcaps") directoryName += "/Endcaps_";
        if(eta == "all") directoryName += "/BarrelAndEncaps_";

        if(r9 == "low") directoryName += "lowR9/";
        if(r9 == "high") directoryName += "highR9/";
        if(r9 == "all") directoryName += "AllR9/";

	if(log == 1) fileName += "_log";

	//TString cut = "(Photon_Et > 25 && isJanLooseMMG == 1";
	//TString cut = "(isJanLooseMMG == 1";
	TString cut = Form("(%s > %s && %s == 1",cutVariable.c_str(), cutVariableValue.c_str(), isAfterFSRCut.c_str());

	if(r9 == "low" && eta == "Barrel") cut += " && Photon_isEB == 1 && Photon_r9 < 0.94";
        if(r9 == "high" && eta == "Barrel") cut += " && Photon_isEB == 1 && Photon_r9 > 0.94";
        if(r9 == "low" && eta == "Endcaps") cut += " && Photon_isEE == 1 && Photon_r9 < 0.95";
        if(r9 == "high" && eta == "Endcaps") cut += " && Photon_isEE == 1 && Photon_r9 > 0.95";
        if(r9 == "all" && eta == "Barrel") cut += " && Photon_isEB == 1";
        if(r9 == "all" && eta == "Endcaps") cut += " && Photon_isEE == 1";
	if(r9 == "all" && eta == "all") cut += " && (Photon_isEE == 1 || Photon_isEB == 1)";

	//cut += " && ( hltnames == \"HLT_Mu17_TkMu8_v9\" || hltnames == \"HLT_Mu17_TkMu8_v10\" || hltnames == \"HLT_Mu17_TkMu8_v11\" || hltnames == \"HLT_Mu17_TkMu8_v12\" || hltnames == \"HLT_Mu17_TkMu8_v13\" || hltnames == \"HLT_Mu17_TkMu8_v14\" || hltnames == \"HLT_Mu17_TkMu8_v15\" || hltnames == \"HLT_Mu17_TkMu8_v16\" )"; //FIXME	

	//cut += " && hltnames == \"HLT_Mu17_TkMu8_v9\"";

	if(normalization == "lumi") cut += ")*weight_pileUp*weight_Xsection"; 
	if(normalization == "integral") cut += ")*weight_pileUp"; 
	if(normalization == "lumi2") cut += ")*weight_Xsection";
	if(normalization == "integral2") cut += ")";
	//if(normalization == "lumi") cut += ")*weight_Xsection";

	cout<<endl<<"cut = "<<cut<<endl;

        TChain * dataChain = new TChain("miniTree");
	TChain * dYToMuMuFSRChain = new TChain("miniTree");
	TChain * dYToMuMuNonFSRChain = new TChain("miniTree");
	TChain * ttJetsChain = new TChain("miniTree");
	TChain * wJetsChain = new TChain("miniTree");
/*
	dataChain->Add("/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/cvs_developpment/Selection_NewMuID/miniTree_Run2012A_13Jul2012_v1_NewMuonID_NewSelection_0_injRe0.00_v1_partALL.root");
	dataChain->Add("/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/cvs_developpment/Selection_NewMuID/miniTree_Run2012A_recover_06Aug2012_v1_NewMuonID_NewSelection_0_injRe0.00_v1_partALL.root");
	dataChain->Add("/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/cvs_developpment/Selection_NewMuID/miniTree_Run2012B_13Jul2012_v4_NewMuonID_NewSelection_0_injRe0.00_v1_partALL.root");
	dataChain->Add("/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/cvs_developpment/Selection_NewMuID/miniTree_Run2012C-24Aug2012-v1_NewMuonID_NewSelection_0_injRe0.00_v1_partALL.root");
	dataChain->Add("/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/cvs_developpment/Selection_NewMuID/miniTree_Run2012C-EcalRecover_11Dec2012-v1_NewMuonID_NewSelection_0_injRe0.00_v1_partALL.root");
	dataChain->Add("/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/cvs_developpment/Selection_NewMuID/miniTree_Run2012C_PromptReco_v2_NewMuonID_NewSelection_0_injRe0.00_v1_partALL.root");
	dataChain->Add("/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/cvs_developpment/Selection_NewMuID/miniTree_Run2012D_PromptReco_v1_NewMuonID_NewSelection_0_injRe0.00_v1_partALL.root");
	//dYToMuMuFSRChain->Add("/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/cvs_developpment/Selection_NewMuID/miniTree_DYToMuMu_Summer12_NewMuonID_NewSelection_1_injRe0_PUTest_v5_partALL.root");
	//dYToMuMuNonFSRChain->Add("/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/cvs_developpment/Selection_NewMuID/miniTree_DYToMuMu_Summer12_NewMuonID_NewSelection_2_injRe0_PUTest_v5_partALL.root");
	dYToMuMuFSRChain->Add("/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/cvs_developpment/Selection_NewMuID/miniTree_DYToMuMu_Summer12_NewMuonID_NewSelection_1_injRe0.00_v1_partALL.root");
	dYToMuMuNonFSRChain->Add("/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/cvs_developpment/Selection_NewMuID/miniTree_DYToMuMu_Summer12_NewMuonID_NewSelection_2_injRe0.00_v1_partALL.root");
	ttJetsChain->Add("/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/cvs_developpment/Selection_NewMuID/miniTree_TTJets_Summer12_S7_NewMuonID_NewSelection_3_injRe0.00_v1_partALL.root");
	wJetsChain->Add("/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/cvs_developpment/Selection_NewMuID/miniTree_WJetsToLNu_Summer12_S10_NewMuonID_NewSelection_3_injRe0.00_v1_partALL.root");
	//dYToMuMuFSRChain->Add("/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/cvs_developpment/Selection_NewMuID/miniTree_DYToMuMu_Summer12_NewMuonID_NewSelection_1_injRe0.00_v2_NewPU_partALL.root");
	//dYToMuMuNonFSRChain->Add("/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/cvs_developpment/Selection_NewMuID/miniTree_DYToMuMu_Summer12_NewMuonID_NewSelection_2_injRe0.00_v2_NewPU_partALL.root");
	//ttJetsChain->Add("/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/cvs_developpment/Selection_NewMuID/miniTree_TTJets_Summer12_S7_NewMuonID_NewSelection_3_injRe0.00_v2_NewPU_partALL.root");
	//wJetsChain->Add("/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/cvs_developpment/Selection_NewMuID/miniTree_WJetsToLNu_Summer12_S10_NewMuonID_NewSelection_3_injRe0.00_v2_NewPU_partALL.root");

*/

	dataChain->Add("/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/cvs_developpment/Selection_NewMuID/miniTree_Run2012A_13Jul2012_v1_NewMuonID_NewSelection_0_injRe0_v5_partALL.root");
        dataChain->Add("/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/cvs_developpment/Selection_NewMuID/miniTree_Run2012A_recover_06Aug2012_v1_NewMuonID_NewSelection_0_injRe0_v5_partALL.root");
        dataChain->Add("/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/cvs_developpment/Selection_NewMuID/miniTree_Run2012B_13Jul2012_v4_NewMuonID_NewSelection_0_injRe0_v5_partALL.root");
        dataChain->Add("/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/cvs_developpment/Selection_NewMuID/miniTree_Run2012C-24Aug2012-v1_NewMuonID_NewSelection_0_injRe0_v5_partALL.root");
        dataChain->Add("/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/cvs_developpment/Selection_NewMuID/miniTree_Run2012C-EcalRecover_11Dec2012-v1_NewMuonID_NewSelection_0_injRe0_v5_partALL.root");
        dataChain->Add("/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/cvs_developpment/Selection_NewMuID/miniTree_Run2012C_PromptReco_v2_NewMuonID_NewSelection_0_injRe0_v5_partALL.root");
        dataChain->Add("/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/cvs_developpment/Selection_NewMuID/miniTree_Run2012D_PromptReco_v1_NewMuonID_NewSelection_0_injRe0_v5_partALL.root");
        //dYToMuMuFSRChain->Add("/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/cvs_developpment/Selection_NewMuID/miniTree_DYToMuMu_Summer12_NewMuonID_NewSelection_1_injRe0_v5_partALL.root");
        //dYToMuMuNonFSRChain->Add("/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/cvs_developpment/Selection_NewMuID/miniTree_DYToMuMu_Summer12_NewMuonID_NewSelection_2_injRe0_v5_partALL.root");
        //ttJetsChain->Add("/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/cvs_developpment/Selection_NewMuID/miniTree_TTJets_Summer12_S7_NewMuonID_NewSelection_3_injRe0_v5_partALL.root");
        //wJetsChain->Add("/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/cvs_developpment/Selection_NewMuID/miniTree_WJetsToLNu_Summer12_S10_NewMuonID_NewSelection_3_injRe0_v5_partALL.root");
	//dYToMuMuFSRChain->Add("/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/cvs_developpment/Selection_NewMuID/miniTree_DYToMuMu_Summer12_NewMuonID_NewSelection_1_injRe0_PUAfterSelection_partALL.root");
	//dYToMuMuNonFSRChain->Add("/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/cvs_developpment/Selection_NewMuID/miniTree_DYToMuMu_Summer12_NewMuonID_NewSelection_2_injRe0_PUAfterSelection_partALL.root");
	//ttJetsChain->Add("/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/cvs_developpment/Selection_NewMuID/miniTree_TTJets_Summer12_S7_NewMuonID_NewSelection_3_injRe0_PUAfterSelection_partALL.root");
	//wJetsChain->Add("/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/cvs_developpment/Selection_NewMuID/miniTree_WJetsToLNu_Summer12_S10_NewMuonID_NewSelection_3_injRe0_PUAfterSelection_partALL.root");	

	dYToMuMuFSRChain->Add("/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/cvs_developpment/Selection_NewMuID/miniTree_DYToMuMu_Summer12_NewMuonID_NewSelection_1_injRe0_PU_true_May2013_MCOfficials_partALL.root");
        dYToMuMuNonFSRChain->Add("/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/cvs_developpment/Selection_NewMuID/miniTree_DYToMuMu_Summer12_NewMuonID_NewSelection_2_injRe0_PU_true_May2013_MCOfficials_partALL.root");
        ttJetsChain->Add("/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/cvs_developpment/Selection_NewMuID/miniTree_TTJets_Summer12_S7_NewMuonID_NewSelection_3_injRe0_PU_true_May2013_MCOfficials_partALL.root");
        wJetsChain->Add("/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/cvs_developpment/Selection_NewMuID/miniTree_WJetsToLNu_Summer12_S10_NewMuonID_NewSelection_3_injRe0_PU_true_May2013_MCOfficials_partALL.root");


	TCanvas *c1 = new TCanvas("c1", "c1",0,0,600,600);
	if(log == 1) 
	{	
		c1->SetLogy();
		cout<<endl<<"c1->SetLogy() !"<<endl;
	}

	TH1D *data = new TH1D("data","data", nBins, xMin, xMax);	
	TH1D *dYToMuMuFSR = new TH1D("dYToMuMuFSR","dYToMuMuFSR", nBins, xMin, xMax);
	TH1D *dYToMuMuNonFSR = new TH1D("dYToMuMuNonFSR","dYToMuMuNonFSR", nBins, xMin, xMax);
	TH1D *ttJets = new TH1D("ttJets","ttJets", nBins, xMin, xMax);
	TH1D *wJets = new TH1D("wJets","wJets", nBins, xMin, xMax);				
	
	cout<<endl<<"coucou";
	
	if(xVariable != "nGenVertices") dataChain->Draw(Form("%s>>data",xVariable.c_str()),cut);
	cout<<endl<<"coucou1";
	dYToMuMuFSRChain->Draw(Form("%s>>dYToMuMuFSR",xVariable.c_str()),cut);
	cout<<endl<<"coucou1b";
	dYToMuMuNonFSRChain->Draw(Form("%s>>dYToMuMuNonFSR",xVariable.c_str()),cut);
	cout<<endl<<"coucou1c";
	ttJetsChain->Draw(Form("%s>>ttJets",xVariable.c_str()),cut);
	cout<<endl<<"coucou1d";
	wJetsChain->Draw(Form("%s>>wJets",xVariable.c_str()),cut);	

	cout<<endl<<"coucou2";

	// --- 2012 Lumi --- //
        double lumidata = 808.472 + 82.136 + 4429.0 + 495.003 + 134.242 + 6397.0 + 7274.0;
        //double lumidata = 808.472 + 82.136 + 4429.0 + 495.003 + 134.242 + 6397.0;
	double lumiDY = 48819386.0 / 1914.894;
        double lumiTtJets = 6736135.0 / 234.0;
        double lumiwJets = 57709905.0 / 37509.25;	

	// --- 2011 Lumi --- //
	//double lumidata = (0.706370 + 0.385819 + 2.741 + 1.099) * 1000;
	//double lumiDY = 29743564.0 / 1665.835;
	//double lumiTtJets = 3701947.0 / 165.0;
	//double lumiwJets = 81345381.0 / 31314.0;	

	// --- 2011 Lumi old --- //
	//double lumiDY = 29743564.0 / 1626.0;
        //double lumiTtJets = 3701947.0 / 94.76;
        //double lumiwJets = 81345381.0 / 27770.0;

	int nData = dataChain->GetEntries(cut);
	
	double nMCw = dYToMuMuFSRChain->GetEntries(cut) * lumidata / lumiDY + dYToMuMuNonFSRChain->GetEntries(cut) * lumidata / lumiDY + ttJetsChain->GetEntries(cut) * lumidata / lumiTtJets + wJetsChain->GetEntries(cut) * lumidata / lumiwJets;	
	
	int nDYFSR = dYToMuMuFSRChain->GetEntries(cut);
	int nDYNonFSR = dYToMuMuNonFSRChain->GetEntries(cut);
	int nTtJets = ttJetsChain->GetEntries(cut);
	int nWJets = wJetsChain->GetEntries(cut);

	double nDYFSRw = dYToMuMuFSRChain->GetEntries(cut) * lumidata / lumiDY;
        double nDYNonFSRw = dYToMuMuNonFSRChain->GetEntries(cut) * lumidata / lumiDY;
        double nTtJetsw = ttJetsChain->GetEntries(cut) * lumidata / lumiTtJets;
        double nWJetsw = wJetsChain->GetEntries(cut) * lumidata / lumiwJets;


	double purity = 100.0 * (dYToMuMuFSRChain->GetEntries(cut) * (lumiDY / lumiDY)) / (dYToMuMuFSRChain->GetEntries(cut) * (lumiDY / lumiDY) + dYToMuMuNonFSRChain->GetEntries(cut) * (lumiDY / lumiDY) + ttJetsChain->GetEntries(cut) * (lumiDY / lumiTtJets) + wJetsChain->GetEntries(cut) * (lumiDY / lumiwJets) );

	cout<<endl<<"purity = "<<purity<<" %"<<endl;

	cout<<endl<<"dYToMuMuFSRChain->GetEntries(cut) = "<<dYToMuMuFSRChain->GetEntries(cut);
	cout<<endl<<"dYToMuMuNonFSRChain->GetEntries(cut) = "<<dYToMuMuNonFSRChain->GetEntries(cut);
	cout<<endl<<"ttJetsChain->GetEntries(cut) = "<<ttJetsChain->GetEntries(cut);
	cout<<endl<<"wJetsChain->GetEntries(cut) = "<<wJetsChain->GetEntries(cut);
	cout<<endl<<endl<<"nData = "<<nData<<", nMCw = "<<nMCw<<", r = "<<nData/nMCw<<endl<<endl;
		
	double eventPerPico = nData / lumidata;
	cout<<endl<<"eventPerPico = "<<eventPerPico<<endl;
	
	

	// --- Weight by integral --- //
	
	double weight = 1.0;
        if(normalization == "integral" || normalization == "integral2")
        {
                ttJets->Scale(lumiDY / lumiTtJets);
                wJets->Scale(lumiDY / lumiwJets);

                dYToMuMuFSR->Add(dYToMuMuNonFSR);
                dYToMuMuFSR->Add(ttJets);
                dYToMuMuFSR->Add(wJets);

                dYToMuMuNonFSR->Add(ttJets);
                dYToMuMuNonFSR->Add(wJets);

                ttJets->Add(wJets);

                //weight = data->GetEntries() / dYToMuMuFSR->GetEntries();      
                if(xVariable != "nGenVertices") weight = data->Integral() / dYToMuMuFSR->Integral();
		else weight = 1.0;

                cout<<endl<<"weight = "<<weight<<endl;

                dYToMuMuFSR->Scale(weight);
                dYToMuMuNonFSR->Scale(weight);
                ttJets->Scale(weight);
                wJets->Scale(weight);

	}
	
	// --- Weight by luminosity --- //
	if(normalization == "lumi" || normalization == "lumi2")
        {
		dYToMuMuFSR->Add(dYToMuMuNonFSR);
        	dYToMuMuFSR->Add(ttJets);
        	dYToMuMuFSR->Add(wJets);

        	dYToMuMuNonFSR->Add(ttJets);
        	dYToMuMuNonFSR->Add(wJets);

        	ttJets->Add(wJets);
	}
	

	nDYFSRw = dYToMuMuFSR->Integral() - dYToMuMuNonFSR->Integral();
        nDYNonFSRw = dYToMuMuNonFSR->Integral() - ttJets->Integral();
        nTtJetsw = ttJets->Integral() - wJets->Integral(); 
        nWJetsw = wJets->Integral(); 
	nMCw = nDYFSRw + nDYNonFSRw + nTtJetsw + nWJetsw;

	cout<<endl<<"dYToMuMuFSR->Integral() = "<<dYToMuMuFSR->Integral()<<endl<<"nMCw = "<<nMCw<<endl;
	cout<<endl<<"data->Integral() = "<<data->Integral()<<endl;

	c1->Clear();

	double meanEtData = data->GetMean();
        double meanEtMC = dYToMuMuFSR->GetMean();
	
	data->SetMarkerStyle(20);
	data->SetMarkerSize(0.5);

	dYToMuMuFSR->GetXaxis()->SetTitle(Form("%s",xVariableName.c_str()));
	dYToMuMuFSR->GetYaxis()->SetTitle(Form("%s",yVariableName.c_str()));
	/*
	dYToMuMuFSR->GetXaxis()->SetLabelFont(42);
	dYToMuMuFSR->GetXaxis()->SetTitleFont(42);
	dYToMuMuFSR->GetYaxis()->SetLabelFont(42);
	dYToMuMuFSR->GetYaxis()->SetTitleFont(42);
	dYToMuMuFSR->GetYaxis()->SetTitleOffset(1.65);
	*/
	dYToMuMuFSR->Draw("");
	dYToMuMuNonFSR->Draw("SAMES");
	ttJets->Draw("SAMES");
	wJets->Draw("SAMES");
	if(xVariable != "nGenVertices") data->Draw("E1SAMES");


	//dYToMuMuFSR->SetFillColor(2);
   	dYToMuMuFSR->SetFillColor(kGreen-7);
	//dYToMuMuFSR->SetFillStyle(3001);
	//dYToMuMuNonFSR->SetFillColor(3);
	dYToMuMuNonFSR->SetFillColor(kAzure-5);
        //dYToMuMuNonFSR->SetFillStyle(3001);
	//ttJets->SetFillColor(4);
	ttJets->SetFillColor(kBlue-1);
        //ttJets->SetFillStyle(3001);
	//wJets->SetFillColor(5);
	wJets->SetFillColor(kCyan+2);
        //wJets->SetFillStyle(3001);

	if(xVariable != "nGenVertices") data->SetName("data");
	dYToMuMuFSR->SetName("dYToMuMuFSR");
	dYToMuMuNonFSR->SetName("dYToMuMuNonFSR");
	ttJets->SetName("ttJets");
	wJets->SetName("wJets");

	TLegend leg(xMinLeg,0.74,xMaxLeg,0.94,"","brNDC");
	leg.SetTextFont(42);
        leg.SetTextSize(0.028);
        leg.SetFillColor(kWhite);
        leg.SetLineColor(kWhite);
        leg.SetShadowColor(kWhite);
        if(xVariable != "nGenVertices") leg.AddEntry(data->GetName(),"data","lep");
        leg.AddEntry(dYToMuMuFSR->GetName(),"Z#mu#mu + #gamma FSR","f");
        leg.AddEntry(dYToMuMuNonFSR->GetName(),"Z#mu#mu + #gamma non FSR","f");
        leg.AddEntry(ttJets->GetName(),"t#bar{t} + jets","f");
        leg.AddEntry(wJets->GetName(),"W + jets","f");
	leg.Draw();

	TLatex latexLabel;
	latexLabel.SetTextFont(42);
	latexLabel.SetTextSize(0.028);
	latexLabel.SetNDC();
	latexLabel.DrawLatex(0.25, 0.96, "CMS Preliminary 2012               #sqrt{s} = 8 TeV               L = 19.6 fb^{-1}");
	
	std::ostringstream cutString2;
        cutString2 << setprecision (2) << fixed << nMCw;
	string cutText = "N_{MC} = " + cutString2.str();	


	std::ostringstream cutString3;
        cutString3 << setprecision (2) << fixed << meanEtData;
        string cutText3 = "mean_{data} = " + cutString3.str();

	std::ostringstream cutString4;
        cutString4 << setprecision (2) << fixed << meanEtMC;
        string cutText4 = "mean_{MC} = " + cutString4.str();	

	std::ostringstream cutString5;
        cutString5 << setprecision (2) << fixed << eventPerPico;
        string cutText5 = "event rate = " + cutString5.str() + " / pb^{-1}";	

	std::ostringstream cutString6;
        cutString6 << setprecision (1) << fixed << purity;
        //cutString6 << setprecision (1) << purity;
	string cutText6 = "purity = " + cutString6.str() + " %";	

	/*
	latexLabel.DrawLatex(0.47, 0.88,Form("N_{data} = %d",nData));
	latexLabel.DrawLatex(0.47, 0.83,Form("%s",cutText.c_str()));
	latexLabel.DrawLatex(0.47, 0.78,Form("%s",cutText3.c_str()));
	latexLabel.DrawLatex(0.47, 0.73,Form("%s",cutText4.c_str()));		
	latexLabel.DrawLatex(0.47, 0.68,Form("%s",cutText6.c_str()));
	if(r9sup == 2 && EndCaps == 2) latexLabel.DrawLatex(0.47, 0.63,Form("%s",cutText5.c_str()));
	*/
		
	double Ymin = 0;
	double Ymax; 
	if(xVariable != "nGenVertices") Ymax = max(dYToMuMuFSR->GetMaximum(),data->GetMaximum()) + max(dYToMuMuFSR->GetMaximum(),data->GetMaximum()) * 0.1;	
	else Ymax = dYToMuMuFSR->GetMaximum() * 0.1;

	
	if(log == 0 && xVariable != "nGenVertices") dYToMuMuFSR->GetYaxis()->SetRangeUser(Ymin,Ymax);
        if(log == 1) dYToMuMuFSR->SetMinimum(pow(10.0,-2));

	plotsRecording(directoryName, fileName, c1);


	c1->Clear();

	string spsDirectory = "/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/Energy_scale_extraction";
	ofstream summaryFile(Form("%s/%s/Summary.txt",spsDirectory.c_str(),directoryName_2.c_str()), ios::app);

	if(log == 0 && xVariable == "Photon_Et")
	{
		summaryFile << eta << " " << r9 << " r9 & " << nData << " & " << nMCw << " & " << nDYFSRw << " & " << nDYNonFSRw << " & " << nTtJetsw << " & " << nWJetsw << " & " << nDYFSR << " & " << nDYNonFSR << " & " << nTtJets << " & " << nWJets << " & " << purity << " \\\\" << endl; 
	}

	summaryFile.close();


	data->Delete();
        data = 0;
	dYToMuMuFSR->Delete();
        dYToMuMuFSR = 0;
	dYToMuMuNonFSR->Delete();
        dYToMuMuNonFSR = 0;
	ttJets->Delete();
        ttJets = 0;
	wJets->Delete();
        wJets = 0;
	dataChain->Delete();
	dataChain = 0;
	dYToMuMuFSRChain->Delete();
        dYToMuMuFSRChain = 0;
        dYToMuMuNonFSRChain->Delete();
        dYToMuMuNonFSRChain = 0;
        ttJetsChain->Delete();
        ttJetsChain = 0;
        wJetsChain->Delete();
        wJetsChain = 0;
	
	delete c1;
	c1 = 0;	
	
	return 0;


}






 
