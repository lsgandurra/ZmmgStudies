#include "fitFunctions.h"
#include "functions.h"
#include "setTDRStyle.C"

int main(int argc, char *argv[])
{
	
	// --- Initialization --- //

	for(int iarg = 0 ; iarg < argc; iarg++)
        {
                cout << "argv[" << iarg << "]= " << argv[iarg] << endl;
        }

	if( argc == 1 )
        {
                cerr << "arguments should be passed : inputFilesDirectory, directoryName, dataType, cutVariable, fitVariable, binFileName, eta, r9, fitFunction, fitPercentage" <<endl;
                return 1;

        }	

	string inputFilesDirectory = "../Zmumugamma_miniTrees_rereco_2012/";
	string directoryName = "Results_ResolutionTest_v1";
	string dataType = "MC";
	string fitVariable = "Mmumugamma";
	string cutVariable = "Photon_Et";
	string binFileName = "LimitesAllPtOneBin.txt";
	string eta = "Barrel"; 
	string r9 = "all";
	string fitFunction = "voigtianXcb"; //"bwXcb"; 
	string fitPercentageS = "100";
	double fitPercentage = 100;
	string injectedResolution = "0";
	int muSys = 0;
	string extraScale = "0";

	if( argc > 1 ) inputFilesDirectory = argv[1];
	if( argc > 2 ) directoryName = argv[2];
	if( argc > 3 ) dataType = argv[3];
	if( argc > 4 ) fitVariable = argv[4];
	if( argc > 5 ) cutVariable = argv[5];
	if( argc > 6 ) binFileName = argv[6];
	if( argc > 7 ) eta = argv[7];
	if( argc > 8 ) r9 = argv[8];
	if( argc > 9 ) fitFunction = argv[9];
	if( argc > 10 )
	{
		fitPercentageS = argv[10];
		std::stringstream ss ( argv[10] );
		ss >> fitPercentage;
	}
	if( argc > 11 ) injectedResolution = argv[11];
	if( argc > 12 ) 
	{	
		std::stringstream ss ( argv[12] );
		ss >> muSys;
	}
	if( argc > 13 ) extraScale = argv[13];

	int goodLine = 0;
        if(muSys > 0)
        {
		if(eta == "Barrel" && r9 == "low") goodLine = 1;
		if(eta == "Barrel" && r9 == "high") goodLine = 2;
		if(eta == "Barrel" && r9 == "all") goodLine = 3;
		if(eta == "Endcaps" && r9 == "low") goodLine = 4;
		if(eta == "Endcaps" && r9 == "high") goodLine = 5;
		if(eta == "Endcaps" && r9 == "all") goodLine = 6;
		ifstream goodRangesFile(Form("SelectedFits_%s.txt",dataType.c_str()));
		for(int i = 0; i < goodLine; i++)
		{
			goodRangesFile >> fitPercentage;
		}
		cout<<endl<<"fitPercentage = "<<fitPercentage;
		goodRangesFile.close();
        }

	if(extraScale != "0" && extraScale != "Zee_Regression" && extraScale != "Zee_RecoEnergy") //Results_v6
        {
		if(dataType == "data")
		{
			if(eta == "Barrel" && r9 == "low") extraScale = "0.999909";
			if(eta == "Barrel" && r9 == "high") extraScale = "0.988684";
			if(eta == "Barrel" && r9 == "all") extraScale = "0.994057";
			if(eta == "Endcaps" && r9 == "low") extraScale = "0.994712";
                        if(eta == "Endcaps" && r9 == "high") extraScale = "0.978923";
                        if(eta == "Endcaps" && r9 == "all") extraScale = "0.987336";
		}

		if(dataType == "MC")
                {
                        if(eta == "Barrel" && r9 == "low") extraScale = "0.997083";
                        if(eta == "Barrel" && r9 == "high") extraScale = "0.99926";
                        if(eta == "Barrel" && r9 == "all") extraScale = "0.998738";
                        if(eta == "Endcaps" && r9 == "low") extraScale = "0.989083";
                        if(eta == "Endcaps" && r9 == "high") extraScale = "0.995002";
                        if(eta == "Endcaps" && r9 == "all") extraScale = "0.991969";
                }

	}



	if(muSys == 0) directoryName += Form("/InjectedResolution_%sPercent/%s/%sPercents/%s_%s_%sR9_%s/",injectedResolution.c_str(),dataType.c_str(),fitPercentageS.c_str(),fitVariable.c_str(),eta.c_str(),r9.c_str(),fitFunction.c_str());
	if(muSys > 0) directoryName += Form("/%s/toy_%d/%s_%s_%sR9_%s/",dataType.c_str(),muSys,fitVariable.c_str(),eta.c_str(),r9.c_str(),fitFunction.c_str());
	system(Form("mkdir -p %s",directoryName.c_str()));
	if(muSys == 0) fitPercentage = fitPercentage / 100.0;
	string fileName = "fit";
	string fitVariableName = "";	
	string lumi = "19.6";
	if(fitVariable == "mmg_s") fitVariableName = "s";
	if(fitVariable == "mmg_s_true") fitVariableName = "s_{TRUE}";		
	if(fitVariable == "Mmumugamma") fitVariableName = "M_{#mu#mu#gamma} (GeV)";
	string cutVariableName = "";
	if(cutVariable == "Photon_Et") cutVariableName = "P_{T}^{#gamma} (GeV)";
	if(cutVariable == "Photon_E") cutVariableName = "E^{#gamma} (GeV)";
        if(cutVariable == "Photon_SC_rawEt") cutVariableName = "E_{T RAW}^{#gamma} (GeV)";
	if(cutVariable == "Photon_SC_Eta") cutVariableName = "#eta";
	if(cutVariable == "Photon_SC_brem") cutVariableName = "#sigma_{#phi}/#sigma_{#eta}";
	TString cut = "isJanLooseMMG == 1";
        int nBins = rowsNumberInFile(binFileName.c_str()) - 1;
	int accessNumber = 0;
	float variableF = 0;
        double rangeMin = -0.3;
        double rangeMax = 0.3;
	double tempNumber = 0.0;
        vector <double> lowCutVariable;
	vector <double> highCutVariable;
	vector <double> fitParameters; // numBin, nb fit parameters, mean, mean error, sigma, sigma error, ..., chi2, degrees of freedom, p-value, sigmaEff, x-value.
        TH1D * fitVariableTH1D;
	TH1D * cutVariableTH1D;
	TLatex latexLabel;
	RooHist* residuals = NULL;
	RooHist* pulls = NULL;

	double xMinCutVariable, xMaxCutVariable, xMinFitVariable, xMaxFitVariable, yMinFitVariable, yMaxFitVariable, yMaxChi2, yMaxSigmaEff;
        if(cutVariable == "Photon_Et") 
	{
		xMinCutVariable = 0.0; 
		xMaxCutVariable = 250.0;
	}
	if(fitVariable == "mmg_s") 
	{
		xMinFitVariable = -0.5; 
		xMaxFitVariable = 0.5;
	}
	if(fitVariable == "mmg_s_true") 
	{
		xMinFitVariable = -0.2; 
		xMaxFitVariable = 0.2;
	}
	if(fitVariable == "Mmumugamma")
	{
                xMinFitVariable = 70.0;
                xMaxFitVariable = 110.0;
        }
        yMinFitVariable = -3.0;
        yMaxFitVariable = 3.0;
        yMaxChi2 = 0;
        yMaxSigmaEff = 0;


	// --- style of plots --- //
	gROOT->Reset();
        TGaxis::SetMaxDigits(3);
        setTDRStyle();
	TCanvas * c1 = new TCanvas("c1", "c1",0,0,600,600);
	TLine line(xMinFitVariable,0,xMaxFitVariable,0);
        line.SetLineStyle(3);
	
	// --- Input root files --- //
	TChain * chain = new TChain("miniTree");
	TChain * reducedChain = 0;

	//chain->Add(Form("%s/*partALL.root",inputFilesDirectory.c_str()));
	if(muSys == 0 && extraScale == "0")
        {
                if(dataType == "data")
        	{
                	chain->Add(Form("miniTree_Run2012A_13Jul2012_v1_NewMuonID_NewSelection_0_injRe%s_v6_RecoEnergy_partALL.root",injectedResolution.c_str()));
                	chain->Add(Form("miniTree_Run2012A_recover_06Aug2012_v1_NewMuonID_NewSelection_0_injRe%s_v6_RecoEnergy_partALL.root",injectedResolution.c_str()));
                	chain->Add(Form("miniTree_Run2012B_13Jul2012_v4_NewMuonID_NewSelection_0_injRe%s_v6_RecoEnergy_partALL.root",injectedResolution.c_str()));
                	chain->Add(Form("miniTree_Run2012C-24Aug2012-v1_NewMuonID_NewSelection_0_injRe%s_v6_RecoEnergy_partALL.root",injectedResolution.c_str()));
                	chain->Add(Form("miniTree_Run2012C-EcalRecover_11Dec2012-v1_NewMuonID_NewSelection_0_injRe%s_v6_RecoEnergy_partALL.root",injectedResolution.c_str()));      
                	chain->Add(Form("miniTree_Run2012C_PromptReco_v2_NewMuonID_NewSelection_0_injRe%s_v6_RecoEnergy_partALL.root",injectedResolution.c_str()));
                	chain->Add(Form("miniTree_Run2012D_PromptReco_v1_NewMuonID_NewSelection_0_injRe%s_v6_RecoEnergy_partALL.root",injectedResolution.c_str()));
                	//chain->Add(Form("miniTree_2011A_03Oct2011V1ReReco_toto_v2_NewSelection_0_v6_RecoEnergy_partALL.root"));
                	//chain->Add(Form("miniTree_2011A_05Jul2011ReReco_toto_v2_NewSelection_0_v6_RecoEnergy_partALL.root"));
                	//chain->Add(Form("miniTree_2011A_PromptSkimV5ReReco_toto_v2_NewSelection_0_v6_RecoEnergy_partALL.root"));
                	//chain->Add(Form("miniTree_2011B_PromptSkimV1ReReco_toto_v2_NewSelection_0_v6_RecoEnergy_partALL.root"));
        	}
        	if(dataType == "MC")
        	{
                	//chain->Add(Form("miniTree_DYToMuMu_Summer12_NewMuonID_NewSelection_1_injRe%s_v6_RecoEnergy_partALL.root",injectedResolution.c_str()));
                	//chain->Add(Form("miniTree_DYToMuMu_Summer12_NewMuonID_NewSelection_2_injRe%s_v6_RecoEnergy_partALL.root",injectedResolution.c_str()));	
			chain->Add(Form("/sps/cms/sgandurr/CMSSW_5_3_7_RECO_5_3_3_v4/src/Selection_July2013/miniTree_DYToMuMu_Summer12_NewMuonID_NewSelection_1_injRe%s_v6_RecoEnergy_partALL.root",injectedResolution.c_str()));
                        chain->Add(Form("/sps/cms/sgandurr/CMSSW_5_3_7_RECO_5_3_3_v4/src/Selection_July2013/miniTree_DYToMuMu_Summer12_NewMuonID_NewSelection_2_injRe%s_v6_RecoEnergy_partALL.root",injectedResolution.c_str()));

			//chain->Add(Form("miniTree_TTJets_Summer12_S7_NewMuonID_NewSelection_3_injRe%s_v6_RecoEnergy_partALL.root",injectedResolution.c_str()));
                	//chain->Add(Form("miniTree_WJetsToLNu_Summer12_S10_NewMuonID_NewSelection_3_injRe%s_v6_RecoEnergy_partALL.root",injectedResolution.c_str()));
                	//chain->Add(Form("miniTree_DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Fall11-PU_S6_START42_V14B-v1_September12_NewSelection_1_v6_RecoEnergy_partALL.root"));
                	//chain->Add(Form("miniTree_DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Fall11-PU_S6_START42_V14B-v1_September12_NewSelection_2_v6_RecoEnergy_partALL.root"));
                	//chain->Add(Form("miniTree_TTJets_TuneZ2_7TeV-madgraph-tauola_NewSelection_3_v3_NewPU_partALL.root"));
		}

        }


	if(muSys > 0 && extraScale == "0")
	{
		if(dataType == "data")
                {
                        chain->Add(Form("miniTree_Run2012A_13Jul2012_v1_NewMuonID_NewSelection_0_MuSys_itoy%d_v6_RecoEnergy_partALL.root",muSys));
                        chain->Add(Form("miniTree_Run2012A_recover_06Aug2012_v1_NewMuonID_NewSelection_0_MuSys_itoy%d_v6_RecoEnergy_partALL.root",muSys));
                        chain->Add(Form("miniTree_Run2012B_13Jul2012_v4_NewMuonID_NewSelection_0_MuSys_itoy%d_v6_RecoEnergy_partALL.root",muSys));
                        chain->Add(Form("miniTree_Run2012C-24Aug2012-v1_NewMuonID_NewSelection_0_MuSys_itoy%d_v6_RecoEnergy_partALL.root",muSys));
                        chain->Add(Form("miniTree_Run2012C-EcalRecover_11Dec2012-v1_NewMuonID_NewSelection_0_MuSys_itoy%d_v6_RecoEnergy_partALL.root",muSys));      
                        chain->Add(Form("miniTree_Run2012C_PromptReco_v2_NewMuonID_NewSelection_0_MuSys_itoy%d_v6_RecoEnergy_partALL.root",muSys));
                        chain->Add(Form("miniTree_Run2012D_PromptReco_v1_NewMuonID_NewSelection_0_MuSys_itoy%d_v6_RecoEnergy_partALL.root",muSys));
                        //chain->Add(Form("miniTree_2011A_03Oct2011V1ReReco_toto_v2_NewSelection_0_v6_RecoEnergy_partALL.root"));
                        //chain->Add(Form("miniTree_2011A_05Jul2011ReReco_toto_v2_NewSelection_0_v6_RecoEnergy_partALL.root"));
                        //chain->Add(Form("miniTree_2011A_PromptSkimV5ReReco_toto_v2_NewSelection_0_v6_RecoEnergy_partALL.root"));
                        //chain->Add(Form("miniTree_2011B_PromptSkimV1ReReco_toto_v2_NewSelection_0_v6_RecoEnergy_partALL.root"));
                }           
                if(dataType == "MC")
                {
                        chain->Add(Form("miniTree_DYToMuMu_Summer12_NewMuonID_NewSelection_1_MuSys_itoy%d_v6_RecoEnergy_partALL.root",muSys));
                        chain->Add(Form("miniTree_DYToMuMu_Summer12_NewMuonID_NewSelection_2_MuSys_itoy%d_v6_RecoEnergy_partALL.root",muSys));
                        //chain->Add(Form("miniTree_TTJets_Summer12_S7_NewMuonID_NewSelection_3_MuSys_itoy%d_v6_RecoEnergy_partALL.root",muSys));
                        //chain->Add(Form("miniTree_WJetsToLNu_Summer12_S10_NewMuonID_NewSelection_3_MuSys_itoy%d_v6_RecoEnergy_partALL.root",muSys));
                        //chain->Add(Form("miniTree_DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Fall11-PU_S6_START42_V14B-v1_September12_NewSelection_1_v6_RecoEnergy_partALL.root"));        
                        //chain->Add(Form("miniTree_DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Fall11-PU_S6_START42_V14B-v1_September12_NewSelection_2_v6_RecoEnergy_partALL.root"));        
                        //chain->Add(Form("miniTree_TTJets_TuneZ2_7TeV-madgraph-tauola_NewSelection_3_v6_RecoEnergy_partALL.root"));
                }
	}

	if(extraScale != "0")
        {	
		if(dataType == "data")
                {
                        chain->Add(Form("miniTree_Run2012A_13Jul2012_v1_NewMuonID_NewSelection_0_injRe%s_%s_v6_RecoEnergy_partALL.root",injectedResolution.c_str(),extraScale.c_str()));
                        chain->Add(Form("miniTree_Run2012A_recover_06Aug2012_v1_NewMuonID_NewSelection_0_injRe%s_%s_v6_RecoEnergy_partALL.root",injectedResolution.c_str(),extraScale.c_str()));
                        chain->Add(Form("miniTree_Run2012B_13Jul2012_v4_NewMuonID_NewSelection_0_injRe%s_%s_v6_RecoEnergy_partALL.root",injectedResolution.c_str(),extraScale.c_str()));
                        chain->Add(Form("miniTree_Run2012C-24Aug2012-v1_NewMuonID_NewSelection_0_injRe%s_%s_v6_RecoEnergy_partALL.root",injectedResolution.c_str(),extraScale.c_str()));
                        chain->Add(Form("miniTree_Run2012C-EcalRecover_11Dec2012-v1_NewMuonID_NewSelection_0_injRe%s_%s_v6_RecoEnergy_partALL.root",injectedResolution.c_str(),extraScale.c_str()));
                        chain->Add(Form("miniTree_Run2012C_PromptReco_v2_NewMuonID_NewSelection_0_injRe%s_%s_v6_RecoEnergy_partALL.root",injectedResolution.c_str(),extraScale.c_str()));
                        chain->Add(Form("miniTree_Run2012D_PromptReco_v1_NewMuonID_NewSelection_0_injRe%s_%s_v6_RecoEnergy_partALL.root",injectedResolution.c_str(),extraScale.c_str()));
                }
                if(dataType == "MC")
                {
                        chain->Add(Form("miniTree_DYToMuMu_Summer12_NewMuonID_NewSelection_1_injRe%s_%s_v6_RecoEnergy_partALL.root",injectedResolution.c_str(),extraScale.c_str()));
                        chain->Add(Form("miniTree_DYToMuMu_Summer12_NewMuonID_NewSelection_2_injRe%s_%s_v6_RecoEnergy_partALL.root",injectedResolution.c_str(),extraScale.c_str()));
                        //chain->Add(Form("miniTree_TTJets_Summer12_S7_NewMuonID_NewSelection_3_injRe%s_%s_v6_RecoEnergy_partALL.root",injectedResolution.c_str(),extraScale.c_str()));
                        //chain->Add(Form("miniTree_WJetsToLNu_Summer12_S10_NewMuonID_NewSelection_3_injRe%s_%s_v6_RecoEnergy_partALL.root",injectedResolution.c_str(),extraScale.c_str()));	
		}


	}

	// --- openning of the file containing the bins of the x-variable --- //
        ifstream binFile(binFileName.c_str());
	
	// --- Fit loop --- //
	for(int i = 1; i <= nBins; i++)
	{
		fileName = Form("fit_%d",i);

		fitParameters.push_back(i);	
	
		// -- creation of the cut -- //
		cut.Clear();
		cut += "isJanLooseMMG == 1";
		if(r9 == "low" && eta == "Barrel") cut += " && Photon_isEB == 1 && Photon_r9 < 0.94";
		if(r9 == "high" && eta == "Barrel") cut += " && Photon_isEB == 1 && Photon_r9 > 0.94";
		if(r9 == "low" && eta == "Endcaps") cut += " && Photon_isEE == 1 && Photon_r9 < 0.95";
                if(r9 == "high" && eta == "Endcaps") cut += " && Photon_isEE == 1 && Photon_r9 > 0.95";
		if(r9 == "all" && eta == "Barrel") cut += " && Photon_isEB == 1";
		if(r9 == "all" && eta == "Endcaps") cut += " && Photon_isEE == 1";

		//cut += " && ( hltnames == \"HLT_Mu17_TkMu8_v9\" || hltnames == \"HLT_Mu17_TkMu8_v10\" || hltnames == \"HLT_Mu17_TkMu8_v11\" || hltnames == \"HLT_Mu17_TkMu8_v12\" || hltnames == \"HLT_Mu17_TkMu8_v13\" || hltnames == \"HLT_Mu17_TkMu8_v14\" || hltnames == \"HLT_Mu17_TkMu8_v15\" || hltnames == \"HLT_Mu17_TkMu8_v16\" )";//FIXME	

		cut += " && " + cutVariable + " > ";
		if(i == 1)
                {
                        binFile >> tempNumber;
			lowCutVariable.push_back(tempNumber);
			binFile >> tempNumber;
                        highCutVariable.push_back(tempNumber);
		} 
                if(i > 1) 
                {
                        lowCutVariable.push_back(highCutVariable[i-2]);
                        binFile >> tempNumber;
                        highCutVariable.push_back(tempNumber);

                }
		cut += lowCutVariable[i-1];
		cut += " && " + cutVariable + " <= ";
		cut += highCutVariable[i-1];

		cout<<endl<<"cut = "<<cut<<endl;

		reducedChain = (TChain *) chain->CopyTree(cut);

		// --- fit range estimation  --- //
		rangeEstimator(fitPercentage, chain, cut, rangeMin, rangeMax, fitVariable, variableF, fitParameters);

		if(fitVariable == "Mmumugamma")
		{
			rangeMin = 70.0; // FIXME
			rangeMax = 110.0; // FIXME
		}
		cout<<endl<<"rangeMin = "<<rangeMin<<", rangeMax = "<<rangeMax<<endl;
		cout<<endl<<"xMinFitVariable = "<<xMinFitVariable<<", xMaxFitVariable = "<<xMaxFitVariable<<endl;

		// --- Creation of a RooDataSet to fit --- //
	        RooRealVar variable(fitVariable.c_str(), fitVariableName.c_str(), rangeMin, rangeMax);	
		/*RooRealVar Photon_SC_Eta("Photon_SC_Eta", "Photon_SC_Eta", -3.0, 3.0, "");
	        RooRealVar Photon_SC_Phi("Photon_SC_Phi", "Photon_SC_Phi", -3.5, 3.5, "");
	        RooRealVar Photon_r9("Photon_r9", "Photon_r9", 0.0, 1.0, "");
	        RooRealVar isJanLooseMMG("isJanLooseMMG", "isJanLooseMMG", 0.0, 1.0, "");
	        RooRealVar Photon_Et("Photon_Et", "Photon_Et", 0.0, 250.0, "");
	        RooRealVar Photon_SC_rawEt("Photon_SC_rawEt", "Photon_SC_rawEt", 0.0, 250.0, "");
	        RooRealVar Photon_E("Photon_E", "Photon_E", 0.0, 1000.0, "");
	        RooRealVar Photon_SC_brem("Photon_SC_brem", "Photon_SC_brem", 0.0, 15.0, "");
	        RooRealVar weight_pileUp("weight_pileUp", "weight_pileUp", 0.0, 100);
		RooRealVar weight_hlt_scaleFactors("weight_hlt_scaleFactors", "weight_hlt_scaleFactors", 0.0, 100);
		RooRealVar Photon_isEB("Photon_isEB", "Photon_isEB", 0.0, 1.0, "");
		RooRealVar Photon_isEE("Photon_isEE", "Photon_isEE", 0.0, 1.0, "");
	        RooStringVar hltnames("hltnames","hltnames","HLT_Mu17_TkMu8");
		RooArgSet ntplVars(variable, Photon_SC_Eta, Photon_SC_Phi, Photon_r9, Photon_Et, Photon_SC_rawEt, Photon_E, Photon_SC_brem, isJanLooseMMG);
	        ntplVars.add(weight_pileUp);
		ntplVars.add(weight_hlt_scaleFactors);
		ntplVars.add(Photon_isEB);
		ntplVars.add(Photon_isEE);
		ntplVars.add(hltnames);
		*/
		RooRealVar weight_pileUp("weight_pileUp", "weight_pileUp", 0.0, 100);
		RooRealVar weight_hlt_scaleFactors("weight_hlt_scaleFactors", "weight_hlt_scaleFactors", 0.0, 100);
		RooArgSet ntplVars(variable, weight_pileUp, weight_hlt_scaleFactors);
	        //RooDataSet * dataset = new RooDataSet("dataset", "dataset", reducedChain, ntplVars, "", "weight_pileUp");
		RooDataSet * dataset = new RooDataSet("dataset", "dataset", reducedChain, ntplVars, "", "weight_pileUp * weight_hlt_scaleFactors");
		//RooDataSet * dataset = new RooDataSet("dataset", "dataset", chain, ntplVars, cut, "weight_pileUp * weight_hlt_scaleFactors");
		//RooDataSet * dataset = new RooDataSet("dataset", "dataset", chain, ntplVars, cut, "");//FIXME	

		// --- Second RooDataSet (only for cosmetic considerations) --- // 
	        RooRealVar variable2(fitVariable.c_str(), fitVariableName.c_str(), xMinFitVariable, xMaxFitVariable);
		/*RooArgSet ntplVars2(variable2, Photon_SC_Eta, Photon_SC_Phi, Photon_r9, Photon_Et, Photon_SC_rawEt, Photon_E, Photon_SC_brem, isJanLooseMMG);
	        ntplVars2.add(weight_pileUp);
		ntplVars2.add(weight_hlt_scaleFactors);
		ntplVars2.add(Photon_isEB);
                ntplVars2.add(Photon_isEE);
		ntplVars.add(hltnames);*/
		RooArgSet ntplVars2(variable2, weight_pileUp, weight_hlt_scaleFactors);
		//RooDataSet * dataset2 = new RooDataSet("dataset2", "dataset2", reducedChain, ntplVars2, "", "weight_pileUp");
		RooDataSet * dataset2 = new RooDataSet("dataset2", "dataset2", reducedChain, ntplVars2, "", "weight_pileUp * weight_hlt_scaleFactors");
	        //RooDataSet * dataset2 = new RooDataSet("dataset2", "dataset2", chain, ntplVars2, cut, "weight_pileUp * weight_hlt_scaleFactors");	
		//RooDataSet * dataset2 = new RooDataSet("dataset2", "dataset2", chain, ntplVars2, cut, "");//FIXME
		
	        RooPlot * fitFrame; 
	
	        fitFrame = variable.frame(xMinFitVariable,xMaxFitVariable);
	
		// --- Binning of 0.02 between -17.0 and 17.0 (necessary for high fit ranges) --- //
	        RooBinning b;
		if(fitVariable == "mmg_s") 
		{	
			b.setRange(-17.0,17.0);
			b.addUniform(1700,-17.0,17.0);
		}
	        if(fitVariable == "mmg_s_true") 
		{
			b.setRange(-17.0,17.0);
			b.addUniform(17000,-17.0,17.0);
		}
		if(fitVariable == "Mmumugamma") 
		{
			b.setRange(70.0,110.0);
			b.addUniform(80,70.0,110.0);	
		}
	
		// --- Call of the fit function --- // 
		if(fitFunction == "voigtian") voigtian(dataset, dataset2, variable, fitFrame, b, rangeMin, rangeMax, fitParameters);
		if(fitFunction == "cruijff") cruijff(dataset, dataset2, variable, fitFrame, b, rangeMin, rangeMax, fitParameters);
		if(fitFunction == "voigtianXcb") voigtianXcb(dataset, dataset2, variable, fitFrame, b, rangeMin, rangeMax, fitParameters);
		if(fitFunction == "bwXcb") bwXcb(dataset, dataset2, variable, fitFrame, b, rangeMin, rangeMax, fitParameters);
		if(fitFunction == "voigtianXgauss") voigtianXgauss(dataset, dataset2, variable, fitFrame, b, rangeMin, rangeMax, fitParameters);
	
		// --- Compute the chi2 and the associated p-value and save informations in fitParameters --- //	
		chiSquare(fitFrame,(char *)"mycurve",(char *)"myhist",fitParameters);

		// --- histogram of chi2 pulls and chi2 residuals --- //
		residuals = residHist(fitFrame,(char *)"myhist",(char *)"mycurve", false, directoryName, i);
		pulls = residHist(fitFrame,(char *)"myhist",(char *)"mycurve", true, directoryName, i);

		// --- Draw a nice plot --- //
		RooCurve * pdf = fitFrame->getCurve("mycurve");
		fitFrame = variable2.frame();
		dataset2->plotOn(fitFrame,Name("myhist"),Binning(b));
		fitFrame->Draw();
		pdf->Draw("SAME");
	
		latexLabel.SetTextSize(0.030);
		latexLabel.SetNDC();
		latexLabel.DrawLatex(0.13, 0.96, "CMS Preliminary 2012, #sqrt{s} = 8 TeV");
		
		if(dataType == "MC") latexLabel.DrawLatex(0.17, 0.88, "Simulation");
		if(dataType == "data") latexLabel.DrawLatex(0.17, 0.88, Form("Data, #int L = %s fb^{-1}",lumi.c_str()));
		latexLabel.DrawLatex(0.17, 0.83,Form("ECAL %s",eta.c_str()));	
		if(r9 == "low" && eta == "Barrel") latexLabel.DrawLatex(0.17, 0.78,"E_{T}^{#gamma} > 25 GeV, r9 < 0,94");
		if(r9 == "low" && eta == "Endcaps") latexLabel.DrawLatex(0.17, 0.78,"E_{T}^{#gamma} > 25 GeV, r9 < 0,95");
		if(r9 == "high" && eta == "Barrel") latexLabel.DrawLatex(0.17, 0.78,"E_{T}^{#gamma} > 25 GeV, r9 > 0,94");
		if(r9 == "high" && eta == "Endcaps") latexLabel.DrawLatex(0.17, 0.78,"E_{T}^{#gamma} > 25 GeV, r9 > 0,95");
		if(r9 == "all")	latexLabel.DrawLatex(0.17, 0.78,"E_{T}^{#gamma} > 25 GeV, All r9");

		accessNumber = (i-1) * (3+2*fitParameters[2]+5) + 3;	
		latexLabel.DrawLatex(0.65, 0.88, Form("#color[1]{#mu(E_{#gamma}) = %f #pm %f}",fitParameters[accessNumber], fitParameters[accessNumber + 1]));
	
		plotsRecording(directoryName, fileName, c1);		
	
		c1->Clear();

		// --- chi2 residuals graph --- //

		fileName = Form("Chi2Residuals_%d",i);

		residuals->GetYaxis()->SetTitle("#chi^{2} Residuals");
		residuals->GetXaxis()->SetTitle(fitVariableName.c_str());

		residuals->Draw("AP");
		line.Draw();
		
		latexLabel.DrawLatex(0.13, 0.96, "CMS Preliminary 2012, #sqrt{s} = 7 TeV");
		if(dataType == "MC") latexLabel.DrawLatex(0.17, 0.88, "Simulation");
                if(dataType == "data") latexLabel.DrawLatex(0.17, 0.88, Form("Data, #int L = %s fb^{-1}",lumi.c_str()));
                latexLabel.DrawLatex(0.17, 0.83,Form("ECAL %s",eta.c_str()));   
                if(r9 == "low") latexLabel.DrawLatex(0.17, 0.78,"E_{T}^{#gamma} > 25 GeV, r9 < 0,94");
                if(r9 == "high") latexLabel.DrawLatex(0.17, 0.78,"E_{T}^{#gamma} > 25 GeV, r9 < 0,95");
                if(r9 == "all") latexLabel.DrawLatex(0.17, 0.78,"E_{T}^{#gamma} > 25 GeV, All r9");	
		
		residuals->GetXaxis()->SetLimits(xMinFitVariable,xMaxFitVariable);
		
		plotsRecording(directoryName, fileName, c1);		

		c1->Clear();	

		residuals->Delete();
		residuals = 0;

		// --- chi2 pulls graph --- //

                fileName = Form("Chi2Pulls_%d",i);

                pulls->GetYaxis()->SetTitle("#chi^{2} Pulls");
                pulls->GetXaxis()->SetTitle(fitVariableName.c_str());

                pulls->Draw("AP");
                line.Draw();
     
                latexLabel.DrawLatex(0.13, 0.96, "CMS Preliminary 2012, #sqrt{s} = 7 TeV");
                if(dataType == "MC") latexLabel.DrawLatex(0.17, 0.88, "Simulation");
                if(dataType == "data") latexLabel.DrawLatex(0.17, 0.88, Form("Data, #int L = %s fb^{-1}",lumi.c_str()));
                latexLabel.DrawLatex(0.17, 0.83,Form("ECAL %s",eta.c_str()));
                if(r9 == "low") latexLabel.DrawLatex(0.17, 0.78,"E_{T}^{#gamma} > 25 GeV, r9 < 0,94");
                if(r9 == "high") latexLabel.DrawLatex(0.17, 0.78,"E_{T}^{#gamma} > 25 GeV, r9 < 0,95");
                if(r9 == "all") latexLabel.DrawLatex(0.17, 0.78,"E_{T}^{#gamma} > 25 GeV, All r9");

                pulls->GetXaxis()->SetLimits(xMinFitVariable,xMaxFitVariable);

                plotsRecording(directoryName, fileName, c1);

                c1->Clear();

		pulls->Delete();
		pulls = 0;


		// --- sigma effectif --- //
                fitVariableTH1D = new TH1D("fitVariableTH1D","fitVariableTH1D",1700,-17.0,17.0);
		chain->Draw(Form("%s>>fitVariableTH1D",fitVariable.c_str()),cut);	
		fitParameters.push_back(effSigma(fitVariableTH1D));	

		c1->Clear();

		// --- x-value of cutVariable --- //
	
		if(cutVariable == "Photon_Et") cutVariableTH1D  = new TH1D("cutVariableTH1D","cutVariableTH1D",1000,0.0,250.0);	
		if(cutVariable == "Photon_E") cutVariableTH1D  = new TH1D("cutVariableTH1D","cutVariableTH1D",1000,0.0,1000.0);
		if(cutVariable == "Photon_SC_rawEt") cutVariableTH1D  = new TH1D("cutVariableTH1D","cutVariableTH1D",1000,0.0,250.0);
		if(cutVariable == "Photon_SC_Eta") cutVariableTH1D  = new TH1D("cutVariableTH1D","cutVariableTH1D",1000,-3.0,3.0);
		if(cutVariable == "Photon_SC_brem") cutVariableTH1D  = new TH1D("cutVariableTH1D","cutVariableTH1D",1000,0.0,15.0);
		chain->Draw(Form("%s>>cutVariableTH1D",cutVariable.c_str()),cut);
		fitParameters.push_back(cutVariableTH1D->GetMean(1));

		c1->Clear();
	
		// --- Memory dump --- //

		cutVariableTH1D->Delete();
		cutVariableTH1D = 0;
	
		fitVariableTH1D->Delete();
                fitVariableTH1D = 0;
	
		dataset->Delete();
		dataset = 0;
	
		dataset2->Delete();
	        dataset2 = 0;	
	
		fitFrame->Delete();
		fitFrame = 0;
		
	}
	
	chain->Delete();
	chain = 0;


	// --- Fit parameters recording --- //


	ofstream fitParametersFile(Form("%sfitsInformations.txt",directoryName.c_str()));	

	int nbInfPerFit = 8 + fitParameters[2] * 2;

	for(int i = 0; i < (int) fitParameters.size(); i++)
	{

		if((i % nbInfPerFit) == 0) fitParametersFile << "/// --- Fit number "<< fitParameters[i] << " --- //" <<endl;
		if((i % nbInfPerFit) == 1) fitParametersFile << "entries = " << fitParameters[i] <<endl;
		if((i % nbInfPerFit) == 2) fitParametersFile << "nb of fit parameters : " << fitParameters[i] <<endl; 		
		if(fitFunction == "voigtian")
		{
			if((i % nbInfPerFit) == 3) fitParametersFile << "mean = " << fitParameters[i] <<endl;
			if((i % nbInfPerFit) == 4) fitParametersFile << "mean error = " << fitParameters[i] <<endl;
			if((i % nbInfPerFit) == 5) fitParametersFile << "sigma = " << fitParameters[i] <<endl;
			if((i % nbInfPerFit) == 6) fitParametersFile << "sigma error = " << fitParameters[i] <<endl;
			if((i % nbInfPerFit) == 7) fitParametersFile << "width = " << fitParameters[i] <<endl;
			if((i % nbInfPerFit) == 8) fitParametersFile << "width error = " << fitParameters[i] <<endl;

		}
		if(fitFunction == "cruijff")
                {
                        if((i % nbInfPerFit) == 3) fitParametersFile << "m0 = " << fitParameters[i] <<endl;
                        if((i % nbInfPerFit) == 4) fitParametersFile << "m0 error = " << fitParameters[i] <<endl;
                        if((i % nbInfPerFit) == 5) fitParametersFile << "sigmaL = " << fitParameters[i] <<endl;
                        if((i % nbInfPerFit) == 6) fitParametersFile << "sigmaL error = " << fitParameters[i] <<endl;
                        if((i % nbInfPerFit) == 7) fitParametersFile << "sigmaR = " << fitParameters[i] <<endl;
                        if((i % nbInfPerFit) == 8) fitParametersFile << "sigmaR error = " << fitParameters[i] <<endl;
			if((i % nbInfPerFit) == 9) fitParametersFile << "alphaL = " << fitParameters[i] <<endl;
                        if((i % nbInfPerFit) == 10) fitParametersFile << "alphaL error = " << fitParameters[i] <<endl;
			if((i % nbInfPerFit) == 11) fitParametersFile << "alphaR = " << fitParameters[i] <<endl;
                        if((i % nbInfPerFit) == 12) fitParametersFile << "alphaR error = " << fitParameters[i] <<endl;

                }
		if(fitFunction == "voigtianXcb")
                {
                        if((i % nbInfPerFit) == 3) fitParametersFile << "meanV = " << fitParameters[i] <<endl;
                        if((i % nbInfPerFit) == 4) fitParametersFile << "meanV error = " << fitParameters[i] <<endl;
                        if((i % nbInfPerFit) == 5) fitParametersFile << "sigmaV = " << fitParameters[i] <<endl;
                        if((i % nbInfPerFit) == 6) fitParametersFile << "sigmaV error = " << fitParameters[i] <<endl;
                        if((i % nbInfPerFit) == 7) fitParametersFile << "widthV = " << fitParameters[i] <<endl;
                        if((i % nbInfPerFit) == 8) fitParametersFile << "widthV error = " << fitParameters[i] <<endl;
			if((i % nbInfPerFit) == 9) fitParametersFile << "meanCB = " << fitParameters[i] <<endl;
                        if((i % nbInfPerFit) == 10) fitParametersFile << "meanCB error = " << fitParameters[i] <<endl;
			if((i % nbInfPerFit) == 11) fitParametersFile << "sigmaCB = " << fitParameters[i] <<endl;
                        if((i % nbInfPerFit) == 12) fitParametersFile << "sigmaCB error = " << fitParameters[i] <<endl;
			if((i % nbInfPerFit) == 13) fitParametersFile << "alphaCB = " << fitParameters[i] <<endl;
                        if((i % nbInfPerFit) == 14) fitParametersFile << "alphaCB error = " << fitParameters[i] <<endl;
			if((i % nbInfPerFit) == 15) fitParametersFile << "nCB = " << fitParameters[i] <<endl;
                        if((i % nbInfPerFit) == 16) fitParametersFile << "nCB error = " << fitParameters[i] <<endl;	


                }	
		if(fitFunction == "bwXcb")
                {
                        if((i % nbInfPerFit) == 3) fitParametersFile << "meanBW = " << fitParameters[i] <<endl;
                        if((i % nbInfPerFit) == 4) fitParametersFile << "meanBW error = " << fitParameters[i] <<endl;
                        if((i % nbInfPerFit) == 5) fitParametersFile << "widthBW = " << fitParameters[i] <<endl;
                        if((i % nbInfPerFit) == 6) fitParametersFile << "widthBW error = " << fitParameters[i] <<endl;
                        if((i % nbInfPerFit) == 7) fitParametersFile << "meanCB = " << fitParameters[i] <<endl;
                        if((i % nbInfPerFit) == 8) fitParametersFile << "meanCB error = " << fitParameters[i] <<endl;
                        if((i % nbInfPerFit) == 9) fitParametersFile << "sigmaCB = " << fitParameters[i] <<endl;
                        if((i % nbInfPerFit) == 10) fitParametersFile << "sigmaCB error = " << fitParameters[i] <<endl;
                        if((i % nbInfPerFit) == 11) fitParametersFile << "alphaCB = " << fitParameters[i] <<endl;
                        if((i % nbInfPerFit) == 12) fitParametersFile << "alphaCB error = " << fitParameters[i] <<endl;
                        if((i % nbInfPerFit) == 13) fitParametersFile << "nCB = " << fitParameters[i] <<endl;
                        if((i % nbInfPerFit) == 14) fitParametersFile << "nCB error = " << fitParameters[i] <<endl;    


                }
		if(fitFunction == "voigtianXgauss")
                {
                        if((i % nbInfPerFit) == 3) fitParametersFile << "meanV = " << fitParameters[i] <<endl;
                        if((i % nbInfPerFit) == 4) fitParametersFile << "meanV error = " << fitParameters[i] <<endl;
                        if((i % nbInfPerFit) == 5) fitParametersFile << "sigmaV = " << fitParameters[i] <<endl;
                        if((i % nbInfPerFit) == 6) fitParametersFile << "sigmaV error = " << fitParameters[i] <<endl;
                        if((i % nbInfPerFit) == 7) fitParametersFile << "widthV = " << fitParameters[i] <<endl;
                        if((i % nbInfPerFit) == 8) fitParametersFile << "widthV error = " << fitParameters[i] <<endl;
                        if((i % nbInfPerFit) == 9) fitParametersFile << "meanG = " << fitParameters[i] <<endl;
                        if((i % nbInfPerFit) == 10) fitParametersFile << "meanG error = " << fitParameters[i] <<endl;
                        if((i % nbInfPerFit) == 11) fitParametersFile << "sigmaG = " << fitParameters[i] <<endl;
                        if((i % nbInfPerFit) == 12) fitParametersFile << "sigmaG error = " << fitParameters[i] <<endl;


                }
		if((i % nbInfPerFit) == fitParameters[2] * 2 + 3) fitParametersFile << "chi2 = " << fitParameters[i] <<endl;
		if((i % nbInfPerFit) == fitParameters[2] * 2 + 4) fitParametersFile << "degrees of freedom = " << fitParameters[i] <<endl;
		if((i % nbInfPerFit) == fitParameters[2] * 2 + 5) fitParametersFile << "p-value = " << fitParameters[i] <<endl;
		if((i % nbInfPerFit) == fitParameters[2] * 2 + 6) fitParametersFile << "sigmaEff = " << fitParameters[i] <<endl;
		if((i % nbInfPerFit) == fitParameters[2] * 2 + 7) fitParametersFile << "x-value = " << fitParameters[i] <<endl;
	}


	fitParametersFile.close();

	ofstream fitParametersFileRaw(Form("%sfitsInformationsRaw.txt",directoryName.c_str()));
	for(int i = 0; i < (int) fitParameters.size(); i++)
	{
		fitParametersFileRaw << fitParameters[i] << endl;
	}
	
	fitParametersFileRaw.close();

	// --- Graph of fitVariable (%) vs cutVariable --- //

	fileName = Form("%sVS%s",fitVariable.c_str(),cutVariable.c_str());
	
	double * x = new double[nBins];
	double * xl = new double[nBins];
	double * xr = new double[nBins];
	double * y = new double[nBins];
	double * yl = new double[nBins];
	double * yr = new double[nBins];

	for(int i = 0; i < nBins; i++)
	{	
		x[i] = fitParameters[(2*fitParameters[2]+8)*i + 2*fitParameters[2]+7];
		xl[i] = lowCutVariable[i];
		xr[i] = highCutVariable[i];
		y[i] = 100 * fitParameters[(2*fitParameters[2]+8)*i + 3];
		yl[i] = 100 * fitParameters[(2*fitParameters[2]+8)*i + 4];
		yr[i] =	100 * fitParameters[(2*fitParameters[2]+8)*i + 4];

	}

	TGraphAsymmErrors fitVarVScutVar(nBins,x,y,xl,xr,yl,yr);
	fitVariableName += " (%)";
	fitVarVScutVar.GetXaxis()->SetTitle(cutVariableName.c_str());
	fitVarVScutVar.GetYaxis()->SetTitle(fitVariableName.c_str());	

	fitVarVScutVar.SetMarkerColor(4);
        fitVarVScutVar.SetLineColor(4);
	fitVarVScutVar.SetMarkerStyle(20);
        fitVarVScutVar.SetMarkerSize(0.6);

	fitVarVScutVar.Draw("AP");

	fitVarVScutVar.GetXaxis()->SetLimits(xMinCutVariable,xMaxCutVariable);
	fitVarVScutVar.GetYaxis()->SetRangeUser(yMinFitVariable,yMaxFitVariable);
	
	plotsRecording(directoryName, fileName, c1);

	c1->Clear();	

	// --- Graph of chi2 vs cutVariable --- //

	fileName = Form("chi2VS%s",cutVariable.c_str());

	yMaxChi2 = 0.0;
	for(int i = 0; i < nBins; i++)
        {
		y[i] = fitParameters[(2*fitParameters[2]+8)*i + 2*fitParameters[2]+3];
		if(y[i] > yMaxChi2) yMaxChi2 = y[i];	
	}

	TGraph chi2VScutVaraiable(nBins,x,y);
	
	chi2VScutVaraiable.GetXaxis()->SetTitle(cutVariableName.c_str());
        chi2VScutVaraiable.GetYaxis()->SetTitle("#chi^{2}");

        chi2VScutVaraiable.SetMarkerColor(4);
        chi2VScutVaraiable.SetLineColor(4);
        chi2VScutVaraiable.SetMarkerStyle(20);
        chi2VScutVaraiable.SetMarkerSize(0.6);

        chi2VScutVaraiable.Draw("AP");

        chi2VScutVaraiable.GetXaxis()->SetLimits(xMinCutVariable,xMaxCutVariable);
        chi2VScutVaraiable.GetYaxis()->SetRangeUser(0.0,yMaxChi2 + yMaxChi2 * 0.2);
     
        plotsRecording(directoryName, fileName, c1);

        c1->Clear();


	// --- Graph of p-value vs cutVariable --- //

	fileName = Form("pValueVS%s",cutVariable.c_str());

	for(int i = 0; i < nBins; i++)
        {
                y[i] = fitParameters[(2*fitParameters[2]+8)*i + 2*fitParameters[2]+5];

        }

        TGraph pValueVScutVaraiable(nBins,x,y);

	pValueVScutVaraiable.GetXaxis()->SetTitle(cutVariableName.c_str());
        pValueVScutVaraiable.GetYaxis()->SetTitle("p-value");

        pValueVScutVaraiable.SetMarkerColor(4);
        pValueVScutVaraiable.SetLineColor(4);
        pValueVScutVaraiable.SetMarkerStyle(20);
        pValueVScutVaraiable.SetMarkerSize(0.6);

        pValueVScutVaraiable.Draw("AP");

        pValueVScutVaraiable.GetXaxis()->SetLimits(xMinCutVariable,xMaxCutVariable);
        pValueVScutVaraiable.GetYaxis()->SetRangeUser(0,1);

        plotsRecording(directoryName, fileName, c1);

        c1->Clear();


	// --- Graph of sigmaEff vs cutVariable --- //

	fileName = Form("sigmaEffVS%s",cutVariable.c_str());

	yMaxSigmaEff = 0.0;
	for(int i = 0; i < nBins; i++)
        {
                y[i] = fitParameters[(2*fitParameters[2]+8)*i + 2*fitParameters[2]+6];
		if(y[i] > yMaxSigmaEff) yMaxSigmaEff = y[i];
        }

        TGraph sigmaEffVScutVaraiable(nBins,x,y);

	sigmaEffVScutVaraiable.GetXaxis()->SetTitle(cutVariableName.c_str());
        sigmaEffVScutVaraiable.GetYaxis()->SetTitle("#sigma_{eff}");

        sigmaEffVScutVaraiable.SetMarkerColor(4);
        sigmaEffVScutVaraiable.SetLineColor(4);
        sigmaEffVScutVaraiable.SetMarkerStyle(20);
        sigmaEffVScutVaraiable.SetMarkerSize(0.6);

        sigmaEffVScutVaraiable.Draw("AP");

        sigmaEffVScutVaraiable.GetXaxis()->SetLimits(xMinCutVariable,xMaxCutVariable);
        sigmaEffVScutVaraiable.GetYaxis()->SetRangeUser(0.0,yMaxSigmaEff + yMaxSigmaEff * 0.2);

        plotsRecording(directoryName, fileName, c1);

        c1->Clear();

	delete [] x;
	x = 0;
	delete [] xl;
	xl = 0;
	delete [] xr;
	xr = 0;
	delete [] y;
	y = 0;
	delete [] yl;
	yl = 0;
	delete [] yr;
	yr = 0;
	
	return 0;

}
