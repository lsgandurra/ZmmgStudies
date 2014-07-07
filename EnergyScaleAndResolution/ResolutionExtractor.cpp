#include "fitFunctions.h"
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
                cerr << "arguments should be passed : directoryName, eta, r9, dataType, xVariable, extractedValue, fitFunction, cutVariable, cutVariableValue, analysisVersion" <<endl; 
                return 1;
	}

	string directoryName = "testLog";
	string eta = "Barrel";
	string r9 = "all";
	string dataType = "data";
	string xVariable = "Mmumugamma";
	string extractedValue = "sigmaMuMuGamma2"; //sigmaMu, sigmaMuMu2, sigmaMuMuGamma2 	
	string fitFunction = "bwXcb";
	string cutVariable = "Photon_Et";
	string cutVariableValue = "0";
	string analysisVersion = "2012";

	if( argc > 1 ) directoryName = argv[1];
	if( argc > 2 ) eta = argv[2];
        if( argc > 3 ) r9 = argv[3];	
	if( argc > 4 ) dataType = argv[4];
	if( argc > 5 ) xVariable = argv[5];
	if( argc > 6 ) extractedValue = argv[6];
	if( argc > 7 ) fitFunction = argv[7];
	if( argc > 7 ) cutVariable = argv[8];
	if( argc > 8 ) cutVariableValue = argv[9];	
	if( argc > 9 ) analysisVersion = argv[10];

	string xVariableName = "";
	if(xVariable == "Mmumu") xVariableName = "M_{#mu#mu}";
	if(xVariable == "Mmumugamma") xVariableName = "M_{#mu#mu#gamma}";

	string lumi = "19.75";
        if(analysisVersion == "2011") lumi = "5.1";
        if(analysisVersion == "2012") lumi = "19.75";

	int nBins = 30;
	double xMin, xMax;
	double xMinLeg, xMaxLeg, yMinLeg, legTextSize;

	xMinLeg = 0.67;
	xMaxLeg = 0.91;	
	yMinLeg = 0.74;
	legTextSize = 0.028;


	gROOT->Reset();
        TGaxis::SetMaxDigits(3);
        setTDRStyle();
	//setEgammaStyle();

	string fileName = directoryName;

	if(extractedValue == "sigmaMu" || extractedValue == "sigmaMuMu2") directoryName += Form("/%s/Eta_%s/",extractedValue.c_str(),eta.c_str());
	else directoryName += Form("/%s/%s_cut_%s",extractedValue.c_str(),cutVariable.c_str(), cutVariableValue.c_str());
	string directoryName_2 = directoryName;

	TString cut = "";

	if(extractedValue == "sigmaMuMuGamma2")
	{
		if(eta == "Barrel") directoryName += "/Barrel_";
        	else if(eta == "Endcaps") directoryName += "/Endcaps_";
        	else if(eta == "all") directoryName += "/BarrelAndEncaps_";
		else directoryName += Form("/%s",eta.c_str());

        	if(r9 == "low") directoryName += "lowR9/";
        	if(r9 == "high") directoryName += "highR9/";
        	if(r9 == "all") directoryName += "AllR9/";
	
		cut = Form("(%s > %s && isJanLooseMMG == 1",cutVariable.c_str(), cutVariableValue.c_str());

		if(r9 == "low" && eta == "Barrel") cut += " && Photon_isEB == 1 && Photon_r9 < 0.94";
        	if(r9 == "high" && eta == "Barrel") cut += " && Photon_isEB == 1 && Photon_r9 > 0.94";
        	if(r9 == "low" && eta == "Endcaps") cut += " && Photon_isEE == 1 && Photon_r9 < 0.95";
        	if(r9 == "high" && eta == "Endcaps") cut += " && Photon_isEE == 1 && Photon_r9 > 0.95";
        	if(r9 == "all" && eta == "Barrel") cut += " && Photon_isEB == 1";
        	if(r9 == "all" && eta == "Endcaps") cut += " && Photon_isEE == 1";
		if(r9 == "all" && eta == "all") cut += " && (Photon_isEE == 1 || Photon_isEB == 1)";

		if(r9 == "all" && eta == "0_05") cut += " && Photon_isEB == 1 && abs(Photon_SC_Eta) < 0.5";
		if(r9 == "all" && eta == "05_1") cut += " && Photon_isEB == 1 && abs(Photon_SC_Eta) > 0.5 && abs(Photon_SC_Eta) < 1.0";
		if(r9 == "all" && eta == "1_15") cut += " && Photon_isEB == 1 && abs(Photon_SC_Eta) > 1.0 && abs(Photon_SC_Eta) < 1.5";
		if(r9 == "all" && eta == "15_2") cut += " && Photon_isEE == 1 && abs(Photon_SC_Eta) < 2.0";
		if(r9 == "all" && eta == "2_25") cut += " && Photon_isEE == 1 && abs(Photon_SC_Eta) > 2.0";
		if(r9 == "all" && eta == "15_25") cut += " && Photon_isEE == 1 && abs(Photon_SC_Eta) > 1.5";

		if(r9 == "all" && eta == "Barrel2") cut += " && Photon_isEB == 1 && (abs(MuonM_Eta)) < 1.5 && (abs(MuonP_Eta)) < 1.5";
		if(r9 == "all" && eta == "Endcaps2") cut += " && Photon_isEE == 1 && (abs(MuonM_Eta)) > 1.5 && (abs(MuonP_Eta)) > 1.5";

		//if(analysisVersion == "2012") cut += " && ( hltnames == \"HLT_Mu17_Mu8_v1\" || hltnames == \"HLT_Mu17_Mu8_v2\" || hltnames == \"HLT_Mu17_Mu8_v3\" || hltnames == \"HLT_Mu17_Mu8_v4\" || hltnames == \"HLT_Mu17_Mu8_v5\" || hltnames == \"HLT_Mu17_Mu8_v6\" || hltnames == \"HLT_Mu17_Mu8_v7\" || hltnames == \"HLT_Mu17_Mu8_v8\" || hltnames == \"HLT_Mu17_Mu8_v9\" || hltnames == \"HLT_Mu17_Mu8_v10\" || hltnames == \"HLT_Mu17_Mu8_v11\" || hltnames == \"HLT_Mu17_Mu8_v12\" || hltnames == \"HLT_Mu17_Mu8_v13\" || hltnames == \"HLT_Mu17_Mu8_v14\" || hltnames == \"HLT_Mu17_Mu8_v15\" || hltnames == \"HLT_Mu17_Mu8_v16\" || hltnames == \"HLT_Mu17_Mu8_v17\" || hltnames == \"HLT_Mu17_Mu8_v18\" || hltnames == \"HLT_Mu17_Mu8_v19\" || hltnames == \"HLT_Mu17_Mu8_v20\" || hltnames == \"HLT_Mu17_Mu8_v21\" || hltnames == \"HLT_Mu17_Mu8_v22\" || hltnames == \"HLT_Mu17_Mu8_v23\") ";	
		if(analysisVersion == "2012") cut += " && isGoodHLT_Mu17_Mu8 == 1";
		//if(analysisVersion == "2012") cut += " && isGoodHLT == 1";

		cut += ")*weight_pileUp*weight_Xsection"; //FIXME ^^ or vv 	

	}
	else
	{
		//if(analysisVersion == "2012") cut = "(isMM_nonFSR == 1 && ( hltnames == \"HLT_Mu17_Mu8_v1\" || hltnames == \"HLT_Mu17_Mu8_v2\" || hltnames == \"HLT_Mu17_Mu8_v3\" || hltnames == \"HLT_Mu17_Mu8_v4\" || hltnames == \"HLT_Mu17_Mu8_v5\" || hltnames == \"HLT_Mu17_Mu8_v6\" || hltnames == \"HLT_Mu17_Mu8_v7\" || hltnames == \"HLT_Mu17_Mu8_v8\" || hltnames == \"HLT_Mu17_Mu8_v9\" || hltnames == \"HLT_Mu17_Mu8_v10\" || hltnames == \"HLT_Mu17_Mu8_v11\" || hltnames == \"HLT_Mu17_Mu8_v12\" || hltnames == \"HLT_Mu17_Mu8_v13\" || hltnames == \"HLT_Mu17_Mu8_v14\" || hltnames == \"HLT_Mu17_Mu8_v15\" || hltnames == \"HLT_Mu17_Mu8_v16\" || hltnames == \"HLT_Mu17_Mu8_v17\" || hltnames == \"HLT_Mu17_Mu8_v18\" || hltnames == \"HLT_Mu17_Mu8_v19\" || hltnames == \"HLT_Mu17_Mu8_v20\" || hltnames == \"HLT_Mu17_Mu8_v21\" || hltnames == \"HLT_Mu17_Mu8_v22\" || hltnames == \"HLT_Mu17_Mu8_v23\") )*weight_pileUp*weight_Xsection";
		//if(analysisVersion == "2011" || analysisVersion == "2012") cut = "(isMM_nonFSR == 1)*weight_pileUp*weight_Xsection";
		//if(analysisVersion == "2011" || analysisVersion == "2012") cut = "(isMM_nonFSR == 1 && isGoodHLT == 1)*weight_pileUp*weight_Xsection";
		
		if(analysisVersion == "2011" || analysisVersion == "2012") cut = "(isMM_nonFSR == 1 && isGoodHLT == 1";
		if(eta == "0_05") cut += " && (abs(MuonM_Eta)) < 0.5 && (abs(MuonP_Eta)) < 0.5";
                if(eta == "05_1") cut += " && abs(MuonM_Eta) > 0.5 && abs(MuonP_Eta) > 0.5 && abs(MuonM_Eta) < 1.0 && abs(MuonP_Eta) < 1.0";
                if(eta == "1_15") cut += " && abs(MuonM_Eta) > 1.0 && abs(MuonP_Eta) > 1.0 && abs(MuonM_Eta) < 1.5 && abs(MuonP_Eta) < 1.5";
                if(eta == "15_2") cut += " && abs(MuonM_Eta) > 1.5 && abs(MuonP_Eta) > 1.5 && abs(MuonM_Eta) < 2.0 && abs(MuonP_Eta) < 2.0";
                if(eta == "2_25") cut += " && abs(MuonM_Eta) > 2.0 && abs(MuonP_Eta) > 2.0";
		if(eta == "15_25") cut += " && abs(MuonM_Eta) > 1.5 && abs(MuonP_Eta) > 1.5";
		if(eta == "Barrel2") cut += " && (abs(MuonM_Eta)) < 1.5 && (abs(MuonP_Eta)) < 1.5";
                if(eta == "Endcaps2") cut += " && (abs(MuonM_Eta)) > 1.5 && (abs(MuonP_Eta)) > 1.5";
		cut+= ")*weight_pileUp*weight_Xsection";		

		//cut = "isMM_nonFSR == 1";
		
	}

	cout<<endl<<"cut = "<<cut<<endl;

        TChain * chain = new TChain("miniTree");


	if(extractedValue == "sigmaMuMuGamma2") //FIXME
	{

		if(analysisVersion == "2011") //FIXME
		{
			if(dataType == "data")
			{
				chain->Add("");
			}
			if(dataType == "MC")
                        {
                                chain->Add("");
                        }
		}
		if(analysisVersion == "2012") //FIXME
                {
                        if(dataType == "data")
                        {
                                chain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_totouples_Run2012A_22Jan2013_v1_November2013_0_thesis_June_v1_partALL.root");
				chain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_totouples_parked_Run2012B_22Jan2013_v1_November2013_0_thesis_June_v1_partALL.root");
				chain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_totouples_parked_Run2012C_22Jan2013_v1_November2013_0_thesis_June_v1_partALL.root");
				chain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_totouples_parked_Run2012D_22Jan2013_v1_November2013_0_thesis_June_v1_partALL.root");
                        }
                        if(dataType == "MC")
                        {
                                chain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_totouple_DYToMuMu_Summer12_S10_reg5_1_thesis_June_v1_partALL.root");
                        	chain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_totouple_DYToMuMu_Summer12_S10_reg5_2_thesis_June_v1_partALL.root");
                        	chain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_totouple_TTJets_Summer12_S10_reg5_3_thesis_June_v1_partALL.root");
                        	chain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_totouples_WJetsToLNu_Summer12_S10_November2013_3_thesis_June_v1_partALL.root");
                        }
                }
	}
	else
	{

                if(analysisVersion == "2011") //FIXME
                {
                        if(dataType == "data")
                        {
                                chain->Add("");
                        }
                        if(dataType == "MC")
                        {
                                chain->Add("");
                        }
                }
                if(analysisVersion == "2012") //FIXME
                {
                        if(dataType == "data")
                        {
				chain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_muons_totouples_Run2012A_22Jan2013_v1_noSkim_v1_June_v1_reduced_partALL.root");
                        	chain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_muons_totouples_parked_Run2012B_22Jan2013_v1_noSkim_v1_part1_June_v1_reduced_partALL.root");
                        	chain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_muons_totouples_parked_Run2012B_22Jan2013_v1_noSkim_v1_part2_June_v1_reduced_partALL.root");
                        	chain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_muons_totouples_parked_Run2012C_22Jan2013_v1_noSkim_v1_part1_June_v1_reduced_partALL.root");
                        	chain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_muons_totouples_parked_Run2012C_22Jan2013_v1_noSkim_v1_part2_June_v1_reduced_partALL.root");
                        	chain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_muons_totouples_parked_Run2012C_22Jan2013_v1_noSkim_v1_part3_June_v1_reduced_partALL.root");
                        	chain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_muons_totouples_parked_Run2012D_22Jan2013_v1_noSkim_v1_part1_June_v1_reduced_partALL.root");
                        	chain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_muons_totouples_parked_Run2012D_22Jan2013_v1_noSkim_v1_part2_June_v1_reduced_partALL.root");
                        	chain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_muons_totouples_parked_Run2012D_22Jan2013_v1_noSkim_v1_part3_June_v1_reduced_partALL.root");	
			}
                        if(dataType == "MC")
                        {
				chain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_muons_totouple_DYToMuMu_Summer12_S10_reg5_noSkim_v1_June_v1_reduced_partALL.root");
                        	//chain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_muons_totouple_TTJets_Summer12_S10_reg5_noSkim_v1_June_v1_reduced_partALL.root");
                        	//chain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_muons_totouples_WJetsToLNu_Summer12_S10_noSkim_v1_part1_June_v1_reduced_partALL.root");
				//chain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_muons_totouples_WJetsToLNu_Summer12_S10_noSkim_v1_part2_June_v1_reduced_partALL.root");
				//chain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_muons_totouples_WJetsToLNu_Summer12_S10_noSkim_v1_part3_June_v1_reduced_partALL.root");
			}
                }
        }	

	TChain * reducedChain = (TChain *) chain->CopyTree(cut);

	TCanvas * c1 = new TCanvas("c1", "c1",0,0,600,600);

	vector <double> fitParameters;

	double rangeMin = 85.0;
        double rangeMax = 97.0;
	TChain * reducedChain2 = 0;

	if(extractedValue == "sigmaMu")
	{
		TString cut2 = cut;
		int maxPtMu = 225;
		int nCat = 9;	
	
		TH1D * sigmaMuTh1 = new TH1D("sigmaMuTh1","sigmaMuTh1",nCat,0,maxPtMu);
		int iter = 0;
		double mean = 0;
		for(int i = 50; i < maxPtMu; i+=25)
		{
			fileName = Form("fit_%d",i);

			fitParameters.push_back(i);

			cut2 = Form("MuonM_P > %d && MuonM_P < %d && MuonP_P > %d && MuonP_P < %d",i,i+25,i,i+25);
			if(eta == "0_05") cut2 += " && (abs(MuonM_Eta)) < 0.5 && (abs(MuonP_Eta)) < 0.5";
			if(eta == "05_1") cut2 += " && abs(MuonM_Eta) > 0.5 && abs(MuonP_Eta) > 0.5 && abs(MuonM_Eta) < 1.0 && abs(MuonP_Eta) < 1.0";
			if(eta == "1_15") cut2 += " && abs(MuonM_Eta) > 1.0 && abs(MuonP_Eta) > 1.0 && abs(MuonM_Eta) < 1.5 && abs(MuonP_Eta) < 1.5";
			if(eta == "1_15") cut2 += " && abs(MuonM_Eta) > 1.5 && abs(MuonP_Eta) > 1.5 && abs(MuonM_Eta) < 2.0 && abs(MuonP_Eta) < 2.0";
			if(eta == "1_15") cut2 += " && abs(MuonM_Eta) > 2.0 && abs(MuonP_Eta) > 2.0";			

			if(eta == "Barrel") cut2 += " && abs(MuonM_Eta) < 1.1 && abs(MuonP_Eta) < 1.1";
			if(eta == "Endcaps") cut2 += " && abs(MuonM_Eta) > 1.1 && abs(MuonP_Eta) > 1.1";

			cout<<endl<<"cut2 = "<<cut2<<endl;	
			reducedChain2 = (TChain *) reducedChain->CopyTree(cut2);
			
			TH1D * h_temp = new TH1D("h_temp","h_temp",1000,0,1000);
			reducedChain2->Draw("sqrtMpPpCos>>h_temp");
			mean = h_temp->GetMean();

			cout<<endl<<"mean = "<<mean<<endl;
			c1->Clear();

			RooRealVar variable(xVariable.c_str(), xVariableName.c_str(), rangeMin, rangeMax);
			RooRealVar weight_Xsection("weight_Xsection", "weight_Xsection", 0.0, 100);
			RooRealVar weight_pileUp("weight_pileUp", "weight_pileUp", 0.0, 100);
			RooArgSet ntplVars(variable, weight_pileUp, weight_Xsection);
			RooDataSet * dataset = new RooDataSet("dataset", "dataset", reducedChain2, ntplVars, "", "weight_Xsection * weight_pileUp");	
	
			RooPlot * fitFrame;
			fitFrame = variable.frame(rangeMin, rangeMax);
			RooBinning b;
			b.setRange(rangeMin, rangeMax);
			b.addUniform(40,85.0,97.0);
		
			RooAbsPdf * pdf_fit = NULL;
			if(fitFunction == "bwXcb") 
			{
				fitParameters.push_back(4);
				pdf_fit = bwXcb(dataset, dataset, variable, fitFrame, b, rangeMin, rangeMax, fitParameters);
				chiSquare(fitFrame,(char *)"mycurve",(char *)"myhist",fitParameters,4);
				cout << endl << "ICIIIIII ---->>>>> p-val = " << fitParameters[fitParameters.size() - 1]<<endl;
				double sigmaMuValue = 0;
				cout << endl << "fitParameters[iter * 18 + 9] = " << fitParameters[iter * 18 + 9] << endl;
				sigmaMuValue = (fitParameters[iter * 18 + 9] * 91.188) / (mean);
				cout << endl << "iter = " << iter << ", sigmaMuValue = "<<sigmaMuValue << endl;
				sigmaMuTh1->SetBinContent(iter+3,sigmaMuValue); //FIXME
				double sigmaMuValueError = 0;
				sigmaMuValueError = (fitParameters[iter * 18 + 10] * 91.188) / (mean);
				cout << endl << "iter = " << iter << ", sigmaMuValueError = "<<sigmaMuValueError << endl;
				sigmaMuTh1->SetBinError(iter+3,sigmaMuValueError); //FIXME

			}
			iter++;
			RooCurve * pdf = fitFrame->getCurve("mycurve");
			dataset->plotOn(fitFrame,Name("myhist"),Binning(b));
			fitFrame->Draw();
			pdf->Draw("SAME");	
	
			TLatex latexLabel;
                	latexLabel.SetTextFont(42);
                	latexLabel.SetTextSize(0.028);
                	latexLabel.SetNDC();

                	if(analysisVersion == "2011") latexLabel.DrawLatex(0.25, 0.96, "CMS Private, #sqrt{s} = 7 TeV");
                	if(analysisVersion == "2012") latexLabel.DrawLatex(0.25, 0.96, "CMS Private, #sqrt{s} = 8 TeV");
                	latexLabel.SetTextSize(0.030);
                	if(dataType == "MC") latexLabel.DrawLatex(0.17, 0.88, "Simulation");
                	if(dataType == "data") latexLabel.DrawLatex(0.17, 0.88, Form("Data, #int L = %s fb^{-1}",lumi.c_str()));
			if(dataType == "data") latexLabel.DrawLatex(0.17, 0.82, Form("%d < P^{#mu} < %d",i,i+25));

			plotsRecording(directoryName, fileName, c1);	
			c1->Clear();
			

			h_temp->Delete();
			h_temp = 0;
			pdf_fit->Delete();	
			pdf_fit = 0;
			dataset->Delete();
			dataset = 0;
	
			reducedChain2->Clear();	
		}

		fileName = "sigmaMu_fit";

		sigmaMuTh1->Draw("E1");//E1

		//TF1 * f1 = new TF1("f1","[0]*x+[1]",0,maxPtMu);
		TF1 * f1 = new TF1("f1","[0]*x*x+[1]*x+[2]",0,maxPtMu);
		sigmaMuTh1->Fit(f1);
		f1->Draw("SAMES");
		f1->SetLineWidth(2);
	
		double a = f1->GetParameter(0);
		double b = f1->GetParameter(1);
		double c = f1->GetParameter(2);	
	
		sigmaMuTh1->GetYaxis()->SetTitle("#sigma_{#mu} (GeV)");
        	sigmaMuTh1->GetXaxis()->SetTitle("P^{#mu} (GeV)");
        	sigmaMuTh1->GetYaxis()->SetTitleOffset(1.55);
        	sigmaMuTh1->GetXaxis()->SetTitleOffset(1.40);
		sigmaMuTh1->SetMarkerStyle(20);
		sigmaMuTh1->SetMarkerSize(0.6);	

		TLatex latexLabel;
		latexLabel.SetTextFont(42);
		latexLabel.SetTextSize(0.028);
		latexLabel.SetNDC();
		
		if(analysisVersion == "2011") latexLabel.DrawLatex(0.25, 0.96, "CMS Private, #sqrt{s} = 7 TeV");
                if(analysisVersion == "2012") latexLabel.DrawLatex(0.25, 0.96, "CMS Private, #sqrt{s} = 8 TeV");    
                latexLabel.SetTextSize(0.030);
                if(dataType == "MC") latexLabel.DrawLatex(0.17, 0.88, "Simulation");
                if(dataType == "data") latexLabel.DrawLatex(0.17, 0.88, Form("Data, #int L = %s fb^{-1}",lumi.c_str()));

		latexLabel.DrawLatex(0.17, 0.82, Form("%f * (P^{#mu})^{2} + %f * P^{#mu} + %f",a,b,c));
		
	

		plotsRecording(directoryName, fileName, c1);

		c1->Clear();
		f1->Delete();
		f1 = 0;
		sigmaMuTh1->Delete();
		sigmaMuTh1 = 0;


	}
	if(extractedValue == "sigmaMuMu2")
	{
		fileName = "sigmaMuMu2_fit";

		RooRealVar variable(xVariable.c_str(), xVariableName.c_str(), rangeMin, rangeMax);
                RooRealVar weight_Xsection("weight_Xsection", "weight_Xsection", 0.0, 100);
                RooRealVar weight_pileUp("weight_pileUp", "weight_pileUp", 0.0, 100);
                RooArgSet ntplVars(variable, weight_pileUp, weight_Xsection);
                RooDataSet * dataset = new RooDataSet("dataset", "dataset", reducedChain, ntplVars, "", "weight_Xsection * weight_pileUp");

                RooPlot * fitFrame;
                fitFrame = variable.frame(rangeMin, rangeMax);
                RooBinning b;
                b.setRange(rangeMin, rangeMax);
                b.addUniform(40,85.0,97.0);

                RooAbsPdf * pdf_fit = NULL;
                if(fitFunction == "bwXcb")
                {
                	pdf_fit = bwXcb(dataset, dataset, variable, fitFrame, b, rangeMin, rangeMax, fitParameters);
                }		

		double sigmaMuMu = fitParameters[7];

		RooCurve * pdf = fitFrame->getCurve("mycurve");
                dataset->plotOn(fitFrame,Name("myhist"),Binning(b));
                fitFrame->Draw();
                pdf->Draw("SAME");


		TLatex latexLabel;
                latexLabel.SetTextFont(42);
                latexLabel.SetTextSize(0.028);
                latexLabel.SetNDC();
        
		if(analysisVersion == "2011") latexLabel.DrawLatex(0.25, 0.96, "CMS Private, #sqrt{s} = 7 TeV");
                if(analysisVersion == "2012") latexLabel.DrawLatex(0.25, 0.96, "CMS Private, #sqrt{s} = 8 TeV");
                latexLabel.SetTextSize(0.030);
                if(dataType == "MC") latexLabel.DrawLatex(0.17, 0.88, "Simulation");
                if(dataType == "data") latexLabel.DrawLatex(0.17, 0.88, Form("Data, #int L = %s fb^{-1}",lumi.c_str()));    

	    	latexLabel.DrawLatex(0.17, 0.82, Form("#sigma_{M_{#mu#mu}} = %f",sigmaMuMu));

		plotsRecording(directoryName, fileName, c1);
                c1->Clear();

		pdf_fit->Delete();
                pdf_fit = 0;
                dataset->Delete();
                dataset = 0;

	}
	if(extractedValue == "sigmaMuMuGamma2")
	{
		fileName = "sigmaMuMuGamma2_fit";

	
		fileName = "sigmaMuMu2_fit";

                RooRealVar variable(xVariable.c_str(), xVariableName.c_str(), rangeMin, rangeMax);
                RooRealVar weight_Xsection("weight_Xsection", "weight_Xsection", 0.0, 100);
                RooRealVar weight_pileUp("weight_pileUp", "weight_pileUp", 0.0, 100);
                RooArgSet ntplVars(variable, weight_pileUp, weight_Xsection);
                RooDataSet * dataset = new RooDataSet("dataset", "dataset", reducedChain, ntplVars, "", "weight_Xsection * weight_pileUp");

                RooPlot * fitFrame;
                fitFrame = variable.frame(rangeMin, rangeMax);
                RooBinning b;
                b.setRange(rangeMin, rangeMax);
                b.addUniform(40,85.0,97.0);

                RooAbsPdf * pdf_fit = NULL;
                if(fitFunction == "bwXcb")
                {
                        pdf_fit = bwXcb(dataset, dataset, variable, fitFrame, b, rangeMin, rangeMax, fitParameters);
                }

                double sigmaMuMuGamma = fitParameters[7];
		
		RooCurve * pdf = fitFrame->getCurve("mycurve");
                dataset->plotOn(fitFrame,Name("myhist"),Binning(b));
                fitFrame->Draw();
                pdf->Draw("SAME");


                TLatex latexLabel;
                latexLabel.SetTextFont(42);
                latexLabel.SetTextSize(0.028);
                latexLabel.SetNDC();
                
		if(analysisVersion == "2011") latexLabel.DrawLatex(0.25, 0.96, "CMS Private, #sqrt{s} = 7 TeV");
                if(analysisVersion == "2012") latexLabel.DrawLatex(0.25, 0.96, "CMS Private, #sqrt{s} = 8 TeV");
                latexLabel.SetTextSize(0.030);
                if(dataType == "MC") latexLabel.DrawLatex(0.17, 0.88, "Simulation");
                if(dataType == "data") latexLabel.DrawLatex(0.17, 0.88, Form("Data, #int L = %s fb^{-1}",lumi.c_str()));

		latexLabel.DrawLatex(0.17, 0.82, Form("#sigma_{M_{#mu#mu#gamma}} = %f",sigmaMuMuGamma));

                plotsRecording(directoryName, fileName, c1);
                c1->Clear();

                pdf_fit->Delete();
                pdf_fit = 0;
                dataset->Delete();
                dataset = 0;

		TH1D * h_temp = new TH1D("h_temp","h_temp",1000000,0,1000000);
                reducedChain->Draw("resTerm_1>>h_temp",cut);
                double mean = h_temp->GetMean();
		double rms = h_temp->GetRMS();

		cout << endl << "resTerm_1 = " << mean << ", rms = " << rms << endl;

		c1->Clear();
		h_temp->Delete();
		h_temp = new TH1D("h_temp","h_temp",1000000,0,1000000);
		reducedChain->Draw("resTerm_2>>h_temp",cut);
		mean = h_temp->GetMean();
		rms = h_temp->GetRMS();
		
		cout << endl << "resTerm_2 = " << mean << ", rms = " << rms << endl;
			
		c1->Clear();
                h_temp->Delete();

	}

	
	c1->Clear();

	delete c1;
	c1 = 0;
	chain->Delete();
	chain = 0;

	return 0;


}






 
