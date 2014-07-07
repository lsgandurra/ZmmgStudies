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
	string xVariable = "Photon_SC_Phi";
	int log = 1;	
	string normalization = "lumi"; // "integral"	
	string cutVariable = "Photon_Et";
	string cutVariableValue = "25";
	string analysisVersion = "2011";

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
	if( argc > 9 ) analysisVersion = argv[9];


	int nBins = 30;
	double xMin, xMax;
	double xMinLeg, xMaxLeg, yMinLeg, legTextSize;

	xMinLeg = 0.67;
	xMaxLeg = 0.91;	
	yMinLeg = 0.74;
	legTextSize = 0.028;
	string xVariableName, yVariableName;

	if(xVariable == "Photon_Et")
	{
		nBins = 30;
		xMin = 0;
		xMax = 150;
		//xMax = 1000;
		xVariableName = "E_{T}^{#gamma} [GeV]";
		yVariableName = "Events / 5 GeV";
	}

	if(xVariable == "Photon_E")
        {
                nBins = 40; 
                xMin = 0;
                xMax = 200;
                xVariableName = "E^{#gamma} [GeV]";
                yVariableName = "Events / 5 GeV";
        }	

	if(xVariable == "Photon_SC_brem")
        {
                nBins = 30;
                xMin = 0;
                xMax = 15;
                xVariableName = "#sigma_{#phi} / #sigma_{#eta}";
                yVariableName = "Events / 0.5";
        }

	if(xVariable == "Mmumugamma")
        {
                nBins = 60;
                xMin = 60;
                xMax = 120;
                xVariableName = "M_{#mu#mu#gamma} [GeV]";
                yVariableName = "Events / 1 GeV";
        }

	if(xVariable == "Mmumu")
        {
                nBins = 60; 
                xMin = 30; 
                xMax = 90;
                xVariableName = "M_{#mu#mu} [GeV]";
                yVariableName = "Events / 1 GeV";
        	//xMinLeg = 0.20;
                //xMaxLeg = 0.44;
		//yMinLeg = 0.78;		
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
	
	if(xVariable == "deltaRNear")
        {
                nBins = 50; 
                xMin = 0;
                xMax = 1;
                xVariableName = "#Delta_{R,Near}";
                yVariableName = "Events / 0.02";
                //xMinLeg = 0.20;
                //xMaxLeg = 0.44;
        }	

	if(xVariable == "deltaRFar")
        {
                nBins = 50;
                xMin = 0;
                xMax = 10;
                xVariableName = "#Delta_{R,Far}";
                yVariableName = "Events / 0.2";
                //xMinLeg = 0.20;
                //xMaxLeg = 0.44;
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
	
	if(xVariable == "Photon_SC_Eta")
        {
                nBins = 52; 
                xMin = -2.6;
                xMax = 2.6; 
                xVariableName = "#eta^{#gamma}";
                yVariableName = "Events / 0.1";
                //xMinLeg = 0.20;
                //xMaxLeg = 0.44;
		yMinLeg = 0.78;
		legTextSize = 0.026;
        }	

	if(xVariable == "Photon_SC_Phi")
        {
                nBins = 52;
                xMin = -3.3;
                xMax = 3.3;
                xVariableName = "#phi^{#gamma}";
                yVariableName = "Events / 0.1";
		yMinLeg = 0.78;
                legTextSize = 0.026;
        }

	
	if(xVariable == "MuonM_Eta")
        {
                nBins = 52; 
                xMin = -2.6;
                xMax = 2.6; 
                xVariableName = "#eta^{#mu^{-}}";
                yVariableName = "Events / 0.1";
                
		if(eta == "Endcaps")
		{
			xMinLeg = 0.43;
                	xMaxLeg = 0.67;
		}
		yMinLeg = 0.78;
		legTextSize = 0.026;
        }
	
	if(xVariable == "MuonP_Eta")
        {
                nBins = 52;
                xMin = -2.6;
                xMax = 2.6;
                xVariableName = "#eta^{#mu^{+}}";
                yVariableName = "Events / 0.1";

        	if(eta == "Endcaps")
                {   
                        xMinLeg = 0.43;
                        xMaxLeg = 0.67;
                }
		yMinLeg = 0.78;
		legTextSize = 0.026;
	}


	gROOT->Reset();
        TGaxis::SetMaxDigits(3);
        setTDRStyle();
	//setEgammaStyle();

	string fileName = directoryName;
	directoryName += Form("/%s_cut_%s/%s/Normalization_%s",cutVariable.c_str(), cutVariableValue.c_str(), xVariable.c_str(),normalization.c_str());
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
	TString cut = Form("(%s > %s && isJanLooseMMG == 1",cutVariable.c_str(), cutVariableValue.c_str());

	if(r9 == "low" && eta == "Barrel") cut += " && Photon_isEB == 1 && Photon_r9 < 0.94";
        if(r9 == "high" && eta == "Barrel") cut += " && Photon_isEB == 1 && Photon_r9 > 0.94";
        if(r9 == "low" && eta == "Endcaps") cut += " && Photon_isEE == 1 && Photon_r9 < 0.95";
        if(r9 == "high" && eta == "Endcaps") cut += " && Photon_isEE == 1 && Photon_r9 > 0.95";
        if(r9 == "all" && eta == "Barrel") cut += " && Photon_isEB == 1";
        if(r9 == "all" && eta == "Endcaps") cut += " && Photon_isEE == 1";
	if(r9 == "all" && eta == "all") cut += " && (Photon_isEE == 1 || Photon_isEB == 1)";

	//if(analysisVersion == "2012" || analysisVersion == "2011") cut += " && ( hltnames == \"HLT_Mu17_TkMu8_v9\" || hltnames == \"HLT_Mu17_TkMu8_v10\" || hltnames == \"HLT_Mu17_TkMu8_v11\" || hltnames == \"HLT_Mu17_TkMu8_v12\" || hltnames == \"HLT_Mu17_TkMu8_v13\" || hltnames == \"HLT_Mu17_TkMu8_v14\" || hltnames == \"HLT_Mu17_TkMu8_v15\" || hltnames == \"HLT_Mu17_TkMu8_v16\" )"; //FIXME	

	//if(analysisVersion == "2012") cut += " && ( hltnames == \"HLT_Mu17_Mu8_v1\" || hltnames == \"HLT_Mu17_Mu8_v2\" || hltnames == \"HLT_Mu17_Mu8_v3\" || hltnames == \"HLT_Mu17_Mu8_v4\" || hltnames == \"HLT_Mu17_Mu8_v5\" || hltnames == \"HLT_Mu17_Mu8_v6\" || hltnames == \"HLT_Mu17_Mu8_v7\" || hltnames == \"HLT_Mu17_Mu8_v8\" || hltnames == \"HLT_Mu17_Mu8_v9\" || hltnames == \"HLT_Mu17_Mu8_v10\" || hltnames == \"HLT_Mu17_Mu8_v11\" || hltnames == \"HLT_Mu17_Mu8_v12\" || hltnames == \"HLT_Mu17_Mu8_v13\" || hltnames == \"HLT_Mu17_Mu8_v14\" || hltnames == \"HLT_Mu17_Mu8_v15\" || hltnames == \"HLT_Mu17_Mu8_v16\" || hltnames == \"HLT_Mu17_Mu8_v17\" || hltnames == \"HLT_Mu17_Mu8_v18\" || hltnames == \"HLT_Mu17_Mu8_v19\" || hltnames == \"HLT_Mu17_Mu8_v20\" || hltnames == \"HLT_Mu17_Mu8_v21\" || hltnames == \"HLT_Mu17_Mu8_v22\" || hltnames == \"HLT_Mu17_Mu8_v23\") ";	
	if(analysisVersion == "2012") cut += " && isGoodHLT_Mu17_Mu8 == 1 && Mmumugamma < 80"; //FIXME !!!

	//if(analysisVersion == "2012") cut += " && ( hltnames == \"HLT_Mu17_TkMu8_v1\" || hltnames == \"HLT_Mu17_TkMu8_v2\" || hltnames == \"HLT_Mu17_TkMu8_v3\" || hltnames == \"HLT_Mu17_TkMu8_v4\" || hltnames == \"HLT_Mu17_TkMu8_v5\" || hltnames == \"HLT_Mu17_TkMu8_v6\" || hltnames == \"HLT_Mu17_TkMu8_v7\" || hltnames == \"HLT_Mu17_TkMu8_v8\" || hltnames == \"HLT_Mu17_TkMu8_v9\" || hltnames == \"HLT_Mu17_TkMu8_v10\" || hltnames == \"HLT_Mu17_TkMu8_v11\" || hltnames == \"HLT_Mu17_TkMu8_v12\" || hltnames == \"HLT_Mu17_TkMu8_v13\" || hltnames == \"HLT_Mu17_TkMu8_v14\" || hltnames == \"HLT_Mu17_TkMu8_v15\" || hltnames == \"HLT_Mu17_TkMu8_v16\" || hltnames == \"HLT_Mu17_TkMu8_v17\" || hltnames == \"HLT_Mu17_TkMu8_v18\" || hltnames == \"HLT_Mu17_TkMu8_v19\" || hltnames == \"HLT_Mu17_TkMu8_v20\" || hltnames == \"HLT_Mu17_TkMu8_v21\" || hltnames == \"HLT_Mu17_TkMu8_v22\" || hltnames == \"HLT_Mu17_TkMu8_v23\")";

	//cut += " && hltnames == \"HLT_Mu17_TkMu8_v9\"";

	//if(normalization == "lumi" && analysisVersion == "2012") cut += ")*weight_pileUp*weight_Xsection*weight_hlt_scaleFactors*weight_tightMuId_scaleFactors"; 
	if(normalization == "lumi" && analysisVersion == "2012") cut += ")*weight_pileUp*weight_Xsection"; //FIXME ^^ or vv 	
	//if(normalization == "lumi" && analysisVersion == "2012") cut += ")*weight_pileUp*weight_Xsection*weight_hlt_scaleFactors2_01";
	//if(normalization == "lumi" && analysisVersion == "2012") cut += ")*weight_Xsection";
	if(normalization == "integral" && analysisVersion == "2012") cut += ")*weight_pileUp*weight_hlt_scaleFactors*weight_tightMuId_scaleFactors"; 
	if(normalization == "lumi" && analysisVersion == "2011") cut += ")*weight_pileUp*weight_Xsection";
        if(normalization == "integral" && analysisVersion == "2011") cut += ")*weight_pileUp";

	//if(normalization == "lumi") cut += ")*weight_pileUp*weight_Xsection*weight_hlt_scaleFactors";
	//if(normalization == "integral") cut += ")*weight_pileUp*weight_hlt_scaleFactors";	
	if(normalization == "lumi2") cut += ")*weight_Xsection";
	if(normalization == "integral2") cut += ")";
	//if(normalization == "lumi") cut += ")*weight_Xsection";

	cout<<endl<<"cut = "<<cut<<endl;

        TChain * dataChain = new TChain("miniTree");
	TChain * dYToMuMuFSRChain = new TChain("miniTree");
	TChain * dYToMuMuNonFSRChain = new TChain("miniTree");
	TChain * ttJetsChain = new TChain("miniTree");
	TChain * wJetsChain = new TChain("miniTree");

	// thesis v1 27/03/14
	if(analysisVersion == "2011") //FIXME
	{
		dataChain->Add("miniTree_totouples_Run2011A_12Oct2013_v1_0_thesis_v1_partALL.root");
        	dataChain->Add("miniTree_totouples_Run2011B_12Oct2013_v1_0_thesis_v1_partALL.root");
		dYToMuMuFSRChain->Add("miniTree_totouple_DYToMuMu_Summer11_S13_v2_1_thesis_v3_noR9rescaling_partALL.root");
		dYToMuMuNonFSRChain->Add("miniTree_totouple_DYToMuMu_Summer11_S13_v2_2_thesis_v3_noR9rescaling_partALL.root");
		ttJetsChain->Add("miniTree_totouple_TTJets_Summer11_S13_3_thesis_v3_noR9rescaling_partALL.root");
		wJetsChain->Add("miniTree_totouple_WJetsToLNu_Summer11_S13_3_thesis_v3_noR9rescaling_partALL.root");
	}	

	if(analysisVersion == "2012") //FIXME
        {
                dataChain->Add("miniTree_totouples_Run2012A_22Jan2013_v1_November2013_0_thesis_v4f_partALL.root");
                dataChain->Add("miniTree_totouples_parked_Run2012B_22Jan2013_v1_November2013_0_thesis_v4f_partALL.root");
		dataChain->Add("miniTree_totouples_parked_Run2012C_22Jan2013_v1_November2013_0_thesis_v4f_partALL.root");
		dataChain->Add("miniTree_totouples_parked_Run2012D_22Jan2013_v1_November2013_0_thesis_v4f_partALL.root");
                dYToMuMuFSRChain->Add("miniTree_totouple_DYToMuMu_Summer12_S10_reg5_1_thesis_v4f_partALL.root");
		dYToMuMuNonFSRChain->Add("miniTree_totouple_DYToMuMu_Summer12_S10_reg5_2_thesis_v4f_partALL.root");
                ttJetsChain->Add("miniTree_totouple_TTJets_Summer12_S10_reg5_3_thesis_v4f_partALL.root");
                wJetsChain->Add("miniTree_totouples_WJetsToLNu_Summer12_S10_November2013_3_thesis_v4f_partALL.root");
        	/*dYToMuMuFSRChain->Add("../Selection/miniTree_totouples_DYToMuMu_Summer12_November2013_1_RD1_thesis_v1_partALL.root");
		dYToMuMuNonFSRChain->Add("../Selection/miniTree_totouples_DYToMuMu_Summer12_November2013_2_RD1_thesis_v1_partALL.root");
		ttJetsChain->Add("../Selection/miniTree_totouples_TTJets_FullLeptMGDecays_Summer12_November2013_v2_3_RD1_thesis_v1_partALL.root");
		ttJetsChain->Add("../Selection/miniTree_totouples_TTJets_SemiLeptMGDecays_Summer12_November2013_3_RD1_thesis_v1_partALL.root");
		ttJetsChain->Add("../Selection/miniTree_totouples_TTJets_HadronicMGDecays_Summer12_November2013_v2_temp218_3_RD1_thesis_v1_partALL.root");
		wJetsChain->Add("../Selection/miniTree_totouples_WJetsToLNu_Summer12_S10_November2013_3_thesis_v1_partALL.root");	
		*/
	}


/*
	dataChain->Add("miniTree_totouples_Run2012A_22Jan2013_v1_November2013_0_injRe0_v1_partALL.root");
	dataChain->Add("miniTree_totouples_parked_Run2012B_22Jan2013_v1_November2013_0_injRe0_v1_partALL.root");
	dataChain->Add("miniTree_totouples_parked_Run2012C_22Jan2013_v1_November2013_0_injRe0_v1_partALL.root");
	dataChain->Add("miniTree_totouples_parked_Run2012D_22Jan2013_v1_November2013_0_injRe0_v1_partALL.root");
*/	
	//First Pileup try
	/*
	dYToMuMuFSRChain->Add("miniTree_totouples_DYToMuMu_Summer12_November2013_1_injRe0_v1_partALL.root");
	dYToMuMuNonFSRChain->Add("miniTree_totouples_DYToMuMu_Summer12_November2013_2_injRe0_v1_partALL.root");
	ttJetsChain->Add("miniTree_totouples_TTJets_SemiLeptMGDecays_Summer12_November2013_3_injRe0_v1_partALL.root");
	ttJetsChain->Add("miniTree_totouples_TTJets_FullLeptMGDecays_Summer12_November2013_v2_3_injRe0_v1_partALL.root");
	wJetsChain->Add("miniTree_totouples_WJetsToLNu_Summer12_S10_November2013_3_injRe0_v1_partALL.root");
	*/

	//corr Pileup
/*
	dYToMuMuFSRChain->Add("miniTree_totouples_DYToMuMu_Summer12_November2013_1_injRe0_v2_partALL.root");
        dYToMuMuNonFSRChain->Add("miniTree_totouples_DYToMuMu_Summer12_November2013_2_injRe0_v2_partALL.root");
        ttJetsChain->Add("miniTree_totouples_TTJets_SemiLeptMGDecays_Summer12_November2013_3_injRe0_v2_partALL.root");
        ttJetsChain->Add("miniTree_totouples_TTJets_FullLeptMGDecays_Summer12_November2013_v2_3_injRe0_v2_partALL.root");
        wJetsChain->Add("miniTree_totouples_WJetsToLNu_Summer12_S10_November2013_3_injRe0_v2_partALL.root");
*/
	//pixelcorr Pileup
/*	
	dYToMuMuFSRChain->Add("miniTree_totouples_DYToMuMu_Summer12_November2013_1_injRe0_v3_partALL.root");
        dYToMuMuNonFSRChain->Add("miniTree_totouples_DYToMuMu_Summer12_November2013_2_injRe0_v3_partALL.root");
        ttJetsChain->Add("miniTree_totouples_TTJets_SemiLeptMGDecays_Summer12_November2013_3_injRe0_v3_partALL.root");
        ttJetsChain->Add("miniTree_totouples_TTJets_FullLeptMGDecays_Summer12_November2013_v2_3_injRe0_v3_partALL.root");
        wJetsChain->Add("miniTree_totouples_WJetsToLNu_Summer12_S10_November2013_3_injRe0_v3_partALL.root");
*/
/*
	//pixelcorr Pileup corrected
	dYToMuMuFSRChain->Add("miniTree_totouples_DYToMuMu_Summer12_November2013_1_injRe0_v4_partALL.root");
        dYToMuMuNonFSRChain->Add("miniTree_totouples_DYToMuMu_Summer12_November2013_2_injRe0_v4_partALL.root");
        ttJetsChain->Add("miniTree_totouples_TTJets_SemiLeptMGDecays_Summer12_November2013_3_injRe0_v4_partALL.root");
        ttJetsChain->Add("miniTree_totouples_TTJets_FullLeptMGDecays_Summer12_November2013_v2_3_injRe0_v4_partALL.root");
        wJetsChain->Add("miniTree_totouples_WJetsToLNu_Summer12_S10_November2013_3_injRe0_v4_partALL.root");	
*/
	
/*
	//pixelcorr Pileup corrected + minBias 72000mb
	dYToMuMuFSRChain->Add("miniTree_totouples_DYToMuMu_Summer12_November2013_1_injRe0_v6_partALL.root");
        dYToMuMuNonFSRChain->Add("miniTree_totouples_DYToMuMu_Summer12_November2013_2_injRe0_v6_partALL.root");
        ttJetsChain->Add("miniTree_totouples_TTJets_SemiLeptMGDecays_Summer12_November2013_3_injRe0_v6_partALL.root");
        ttJetsChain->Add("miniTree_totouples_TTJets_FullLeptMGDecays_Summer12_November2013_v2_3_injRe0_v6_partALL.root");
        wJetsChain->Add("miniTree_totouples_WJetsToLNu_Summer12_S10_November2013_3_injRe0_v6_partALL.root");
*/
/*
	//pixelcorr Pileup corrected + lumiOld
	dYToMuMuFSRChain->Add("miniTree_totouples_DYToMuMu_Summer12_November2013_1_injRe0_v7_partALL.root");
        dYToMuMuNonFSRChain->Add("miniTree_totouples_DYToMuMu_Summer12_November2013_2_injRe0_v7_partALL.root");
        ttJetsChain->Add("miniTree_totouples_TTJets_SemiLeptMGDecays_Summer12_November2013_3_injRe0_v7_partALL.root");
        ttJetsChain->Add("miniTree_totouples_TTJets_FullLeptMGDecays_Summer12_November2013_v2_3_injRe0_v7_partALL.root");
        wJetsChain->Add("miniTree_totouples_WJetsToLNu_Summer12_S10_November2013_3_injRe0_v7_partALL.root");
*/
/*
	//pixelcorr Pileup corrected + DY without skim  Data_MC_Comparison_22Jan_v9
	dYToMuMuFSRChain->Add("miniTree_totouples_DYToMuMu_Summer12_February2014_noskim_1_injRe0_v8_partALL.root");
	dYToMuMuNonFSRChain->Add("miniTree_totouples_DYToMuMu_Summer12_February2014_noskim_2_injRe0_v8_partALL.root");
	ttJetsChain->Add("miniTree_totouples_TTJets_SemiLeptMGDecays_Summer12_November2013_3_injRe0_v4_partALL.root");
	ttJetsChain->Add("miniTree_totouples_TTJets_FullLeptMGDecays_Summer12_November2013_v2_3_injRe0_v4_partALL.root");
	ttJetsChain->Add("miniTree_totouples_TTJets_HadronicMGDecays_Summer12_November2013_v2_temp218_3_injRe0_v4_partALL.root");
	wJetsChain->Add("miniTree_totouples_WJetsToLNu_Summer12_S10_November2013_3_injRe0_v4_partALL.root");
*/

/*
	//RD1 ok
	dataChain->Add("miniTree_Run2012A_13Jul2012_v1_NewMuonID_NewSelection_0_injRe0_v6_partALL.root");
        dataChain->Add("miniTree_Run2012A_recover_06Aug2012_v1_NewMuonID_NewSelection_0_injRe0_v6_partALL.root");
        dataChain->Add("miniTree_Run2012B_13Jul2012_v4_NewMuonID_NewSelection_0_injRe0_v6_partALL.root");
        dataChain->Add("miniTree_Run2012C-24Aug2012-v1_NewMuonID_NewSelection_0_injRe0_v6_partALL.root");
        dataChain->Add("miniTree_Run2012C-EcalRecover_11Dec2012-v1_NewMuonID_NewSelection_0_injRe0_v6_partALL.root");
        dataChain->Add("miniTree_Run2012C_PromptReco_v2_NewMuonID_NewSelection_0_injRe0_v6_partALL.root");
        dataChain->Add("miniTree_Run2012D_PromptReco_v1_NewMuonID_NewSelection_0_injRe0_v6_partALL.root");

        dYToMuMuFSRChain->Add("miniTree_DYToMuMu_Summer12_NewMuonID_NewSelection_1_injRe0_v6_partALL.root");
        dYToMuMuNonFSRChain->Add("miniTree_DYToMuMu_Summer12_NewMuonID_NewSelection_2_injRe0_v6_partALL.root");
        ttJetsChain->Add("miniTree_TTJets_Summer12_S7_NewMuonID_NewSelection_3_injRe0_v6_partALL.root");
        wJetsChain->Add("miniTree_WJetsToLNu_Summer12_S10_NewMuonID_NewSelection_3_injRe0_v6_partALL.root");
*/



/*
	dataChain->Add("miniTree_Run2012A_13Jul2012_v1_NewMuonID_NewSelection_0_injRe0_v6_partALL.root");
        dataChain->Add("miniTree_Run2012A_recover_06Aug2012_v1_NewMuonID_NewSelection_0_injRe0_v6_partALL.root");
        dataChain->Add("miniTree_Run2012B_13Jul2012_v4_NewMuonID_NewSelection_0_injRe0_v6_partALL.root");
        dataChain->Add("miniTree_Run2012C-24Aug2012-v1_NewMuonID_NewSelection_0_injRe0_v6_partALL.root");
        dataChain->Add("miniTree_Run2012C-EcalRecover_11Dec2012-v1_NewMuonID_NewSelection_0_injRe0_v6_partALL.root");
        dataChain->Add("miniTree_Run2012C_PromptReco_v2_NewMuonID_NewSelection_0_injRe0_v6_partALL.root");
        dataChain->Add("miniTree_Run2012D_PromptReco_v1_NewMuonID_NewSelection_0_injRe0_v6_partALL.root");
        //dYToMuMuFSRChain->Add("miniTree_DYToMuMu_Summer12_NewMuonID_NewSelection_1_injRe0_v5_partALL.root");
        //dYToMuMuNonFSRChain->Add("miniTree_DYToMuMu_Summer12_NewMuonID_NewSelection_2_injRe0_v5_partALL.root");
        //ttJetsChain->Add("miniTree_TTJets_Summer12_S7_NewMuonID_NewSelection_3_injRe0_v5_partALL.root");
        //wJetsChain->Add("miniTree_WJetsToLNu_Summer12_S10_NewMuonID_NewSelection_3_injRe0_v5_partALL.root");
	//dYToMuMuFSRChain->Add("miniTree_DYToMuMu_Summer12_NewMuonID_NewSelection_1_injRe0_PUAfterSelection_partALL.root");
	//dYToMuMuNonFSRChain->Add("miniTree_DYToMuMu_Summer12_NewMuonID_NewSelection_2_injRe0_PUAfterSelection_partALL.root");
	//ttJetsChain->Add("miniTree_TTJets_Summer12_S7_NewMuonID_NewSelection_3_injRe0_PUAfterSelection_partALL.root");
	//wJetsChain->Add("miniTree_WJetsToLNu_Summer12_S10_NewMuonID_NewSelection_3_injRe0_PUAfterSelection_partALL.root");
	
	//dYToMuMuFSRChain->Add("miniTree_DYToMuMu_Summer12_NewMuonID_NewSelection_1_injRe0_PU_true_May2013_partALL.root");
        //dYToMuMuNonFSRChain->Add("miniTree_DYToMuMu_Summer12_NewMuonID_NewSelection_2_injRe0_PU_true_May2013_partALL.root");
        //ttJetsChain->Add("miniTree_TTJets_Summer12_S7_NewMuonID_NewSelection_3_injRe0_PU_true_May2013_partALL.root");
        //wJetsChain->Add("miniTree_WJetsToLNu_Summer12_S10_NewMuonID_NewSelection_3_injRe0_PU_true_May2013_partALL.root");

	//dYToMuMuFSRChain->Add("miniTree_DYToMuMu_Summer12_NewMuonID_NewSelection_1_injRe0_PU_true_May2013_MCOfficials_partALL.root");
        //dYToMuMuNonFSRChain->Add("miniTree_DYToMuMu_Summer12_NewMuonID_NewSelection_2_injRe0_PU_true_May2013_MCOfficials_partALL.root");
        //ttJetsChain->Add("miniTree_TTJets_Summer12_S7_NewMuonID_NewSelection_3_injRe0_PU_true_May2013_MCOfficials_partALL.root");
        //wJetsChain->Add("miniTree_WJetsToLNu_Summer12_S10_NewMuonID_NewSelection_3_injRe0_PU_true_May2013_MCOfficials_partALL.root");	

	dYToMuMuFSRChain->Add("miniTree_DYToMuMu_Summer12_NewMuonID_NewSelection_1_injRe0_v6_partALL.root");
	dYToMuMuNonFSRChain->Add("miniTree_DYToMuMu_Summer12_NewMuonID_NewSelection_2_injRe0_v6_partALL.root");
	ttJetsChain->Add("miniTree_TTJets_Summer12_S7_NewMuonID_NewSelection_3_injRe0_v6_partALL.root");
	wJetsChain->Add("miniTree_WJetsToLNu_Summer12_S10_NewMuonID_NewSelection_3_injRe0_v6_partALL.root");	

*/	
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
	
	
	dataChain->Draw(Form("%s>>data",xVariable.c_str()),cut);
	dYToMuMuFSRChain->Draw(Form("%s>>dYToMuMuFSR",xVariable.c_str()),cut);
	dYToMuMuNonFSRChain->Draw(Form("%s>>dYToMuMuNonFSR",xVariable.c_str()),cut);
	ttJetsChain->Draw(Form("%s>>ttJets",xVariable.c_str()),cut);
	wJetsChain->Draw(Form("%s>>wJets",xVariable.c_str()),cut);	

	cout<<endl<<"coucou2";

	// --- 2012 Lumi 22Jan--- //
	//double lumidata = 889.362 + 4429.0 + 7114 + 7318;
	double lumidata = 808.472 + 82.136 + 4429.0 + 495.003 + 134.242 + 6397.0 + 7274.0; //FIXME ERASE this line !!!
	double lumiDY = 48851166.0 / 1914.894;
	double lumiTtJets = 24949110.0 / 107.66 + 12066717.0 / 25.81;
	double lumiwJets = 57709905.0 / 37509.25;

	// --- 2012 Lumi --- //
        //double lumidata = 808.472 + 82.136 + 4429.0 + 495.003 + 134.242 + 6397.0 + 7274.0;
        //double lumidata = 808.472 + 82.136 + 4429.0 + 495.003 + 134.242 + 6397.0;
	//double lumiDY = 48819386.0 / 1914.894;
        //double lumiTtJets = 6736135.0 / 234.0;
        //double lumiwJets = 57709905.0 / 37509.25;	

	// --- 2011 Lumi --- //
	//double lumidata = (0.706370 + 0.385819 + 2.741 + 1.099) * 1000;
	//double lumiDY = 29743564.0 / 1665.835;
	//double lumiTtJets = 3701947.0 / 165.0;
	//double lumiwJets = 81345381.0 / 31314.0;	

	// --- 2011 Lumi old --- //
	//double lumiDY = 29743564.0 / 1626.0;
        //double lumiTtJets = 3701947.0 / 94.76;
        //double lumiwJets = 81345381.0 / 27770.0;


	if(analysisVersion == "2011")
	{
		lumidata = 2333 + 2766; 
		lumiDY = 23760702.0 / 1665.835; 
		lumiTtJets = 33476020.0 / 165;
		lumiwJets = 81064825.0 / 31314;
	}

	if(analysisVersion == "2012")
        {
                lumidata = 889.362 + 4429.0 + 7114 + 7318;
                lumiDY = 48851166.0 / 1914.894;
                lumiTtJets = 6923652.0 / 245.8;
                lumiwJets = 57709905.0 / 37509.25;
        }

	cout << endl << "lumidata = "<<lumidata << endl;

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

	cout << endl << "nDYFSRw first = " << nDYFSRw <<endl; 

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
                weight = data->Integral() / dYToMuMuFSR->Integral();

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

	cout << endl << "nDYFSRw end = " << nDYFSRw << endl;

	cout<<endl<<"dYToMuMuFSR->Integral() = "<<dYToMuMuFSR->Integral()<<endl<<"nMCw = "<<nMCw<<endl;
	cout<<endl<<"data->Integral() = "<<data->Integral()<<endl;

	c1->Clear();

	double meanEtData = data->GetMean();
        double meanEtMC = dYToMuMuFSR->GetMean();
	double rmsEtData = data->GetRMS();
	double rmsEtMC = dYToMuMuFSR->GetRMS();
	
	data->SetMarkerStyle(20);
	data->SetMarkerSize(0.5);

	dYToMuMuFSR->GetXaxis()->SetTitle(Form("%s",xVariableName.c_str()));
	dYToMuMuFSR->GetYaxis()->SetTitle(Form("%s",yVariableName.c_str()));
	dYToMuMuFSR->GetYaxis()->SetTitleOffset(1.55);
        dYToMuMuFSR->GetXaxis()->SetTitleOffset(1.40);
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
	data->Draw("E1SAMES");


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

	data->SetName("data");
	dYToMuMuFSR->SetName("dYToMuMuFSR");
	dYToMuMuNonFSR->SetName("dYToMuMuNonFSR");
	ttJets->SetName("ttJets");
	wJets->SetName("wJets");

	TLegend leg(xMinLeg,yMinLeg,xMaxLeg,0.94,"","brNDC");
	leg.SetTextFont(42);
        leg.SetTextSize(legTextSize);
        leg.SetFillColor(kWhite);
        leg.SetLineColor(kWhite);
        leg.SetShadowColor(kWhite);
        leg.AddEntry(data->GetName(),"data","lep");
        leg.AddEntry(dYToMuMuFSR->GetName(),"Z#mu#mu + #gamma FSR","f");
        leg.AddEntry(dYToMuMuNonFSR->GetName(),"Z#mu#mu + #gamma non FSR","f");
        leg.AddEntry(ttJets->GetName(),"t#bar{t} + jets","f");
        leg.AddEntry(wJets->GetName(),"W + jets","f");
	leg.Draw();

	TLatex latexLabel;
	latexLabel.SetTextFont(42);
	latexLabel.SetTextSize(0.028);
	latexLabel.SetNDC();
	//latexLabel.DrawLatex(0.25, 0.96, "CMS Preliminary 2012               #sqrt{s} = 8 TeV               L = 19.6 fb^{-1}");
	//if(analysisVersion == "2011") latexLabel.DrawLatex(0.25, 0.96, "CMS Preliminary 2011               #sqrt{s} = 7 TeV                L = 5.1 fb^{-1}");
	//if(analysisVersion == "2012") latexLabel.DrawLatex(0.25, 0.96, "CMS Preliminary 2012              #sqrt{s} = 8 TeV               L = 19.75 fb^{-1}");		
	if(analysisVersion == "2011") latexLabel.DrawLatex(0.25, 0.96, "CMS Private 2011                 #sqrt{s} = 7 TeV                  L = 5.1 fb^{-1}");
	if(analysisVersion == "2012") latexLabel.DrawLatex(0.25, 0.96, "CMS Private 2012                #sqrt{s} = 8 TeV                 L = 19.75 fb^{-1}");

	std::ostringstream cutString2;
        cutString2 << setprecision (2) << fixed << nMCw;
	string cutText = "N_{MC} = " + cutString2.str();	


	std::ostringstream cutString3;
        cutString3 << setprecision (2) << fixed << meanEtData;
        string cutText3 = "mean_{data} = " + cutString3.str();

	std::ostringstream cutString4;
        cutString4 << setprecision (2) << fixed << meanEtMC;
        string cutText4 = "mean_{MC} = " + cutString4.str();	

	std::ostringstream cutString3b;
        cutString3b << setprecision (2) << fixed << rmsEtData;
        string cutText3b = "rms_{data} = " + cutString3b.str();

        std::ostringstream cutString4b;
        cutString4b << setprecision (2) << fixed << rmsEtMC;
        string cutText4b = "rms_{MC} = " + cutString4b.str();	


	std::ostringstream cutString5;
        cutString5 << setprecision (2) << fixed << eventPerPico;
        string cutText5 = "event rate = " + cutString5.str() + " / pb^{-1}";	

	std::ostringstream cutString6;
        cutString6 << setprecision (1) << fixed << purity;
        //cutString6 << setprecision (1) << purity;
	string cutText6 = "purity = " + cutString6.str() + " %";	

	latexLabel.SetTextSize(0.032);	
	//latexLabel.DrawLatex(0.47, 0.88,Form("N_{data} = %d",nData));
	/*
	if(xVariable == "Mmumugamma")
	{
		latexLabel.DrawLatex(0.20, 0.83,Form("%s",cutText3.c_str()));
		latexLabel.DrawLatex(0.20, 0.78,Form("%s",cutText4.c_str()));
		latexLabel.DrawLatex(0.20, 0.73,Form("%s",cutText3b.c_str()));		
		latexLabel.DrawLatex(0.20, 0.68,Form("%s",cutText4b.c_str()));
	}
	if(xVariable == "Mmumu")
        {
                latexLabel.DrawLatex(0.67, 0.68,Form("%s",cutText3.c_str()));
                latexLabel.DrawLatex(0.67, 0.63,Form("%s",cutText4.c_str()));
                latexLabel.DrawLatex(0.67, 0.58,Form("%s",cutText3b.c_str()));    
                latexLabel.DrawLatex(0.67, 0.53,Form("%s",cutText4b.c_str()));
        }
	if(xVariable == "Photon_Et")
        {
                latexLabel.DrawLatex(0.40, 0.83,Form("%s",cutText3.c_str()));
                latexLabel.DrawLatex(0.40, 0.78,Form("%s",cutText4.c_str()));
                latexLabel.DrawLatex(0.40, 0.73,Form("%s",cutText3b.c_str()));    
                latexLabel.DrawLatex(0.40, 0.68,Form("%s",cutText4b.c_str()));
        }
	*/
	//if(r9sup == 2 && EndCaps == 2) latexLabel.DrawLatex(0.47, 0.63,Form("%s",cutText5.c_str()));
	
		
	double Ymin = 0;
	double Ymax = max(dYToMuMuFSR->GetMaximum(),data->GetMaximum()) + max(dYToMuMuFSR->GetMaximum(),data->GetMaximum()) * 0.1;	
	
	if(log == 0) dYToMuMuFSR->GetYaxis()->SetRangeUser(Ymin,Ymax);
	if(log == 0 && (xVariable == "Photon_SC_Phi" || xVariable == "MuonM_Eta" || xVariable == "MuonP_Eta" || xVariable == "Photon_SC_Eta")) dYToMuMuFSR->GetYaxis()->SetRangeUser(Ymin,Ymax + Ymax * 0.2); 
        if(log == 1) dYToMuMuFSR->SetMinimum(pow(10.0,-2));
	if(log == 1 && (xVariable == "Photon_SC_Phi" || xVariable == "MuonM_Eta" || xVariable == "MuonP_Eta" || xVariable == "Photon_SC_Eta" || xVariable == "Mmumugamma" || xVariable == "Mmumu")) dYToMuMuFSR->SetMaximum(30 * Ymax);

	plotsRecording(directoryName, fileName, c1);


	c1->Clear();

	//string spsDirectory = "/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/Energy_scale_extraction";
	string spsDirectory = "/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/EnergyScaleAndResolution";
	system(Form("mkdir -p %s/%s",spsDirectory.c_str(),directoryName_2.c_str()));

	ofstream summaryFile(Form("%s/%s/Summary.txt",spsDirectory.c_str(),directoryName_2.c_str()), ios::app);
	
	string testChain = Form("%s/%s/Summary.txt",spsDirectory.c_str(),directoryName_2.c_str());
	cout << endl << "testChain = " << testChain << endl;


	//if(log == 0 && xVariable == "Photon_Et")
	if(log == 0 && xVariable == "Photon_SC_Eta")
	{
		
		cout << "coucou dans le if du summary file !!" << endl;
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






 
