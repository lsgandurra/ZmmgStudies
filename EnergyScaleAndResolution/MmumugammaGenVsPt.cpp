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
	string eta = "all";
	string r9 = "all";
	string dataType = "MC";
	string xVariable = "Mmumugamma_MMG_MC";
	string extractedValue = "Mean_MmumugammaGen";  	
	string fitFunction = "voigtian2";
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
	if(xVariable == "Mmumugamma_MMG_MC") xVariableName = "M_{#mu#mu#gamma}^{GEN}";

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

	if(extractedValue == "Mean_MmumugammaGen") directoryName += Form("/%s/%s_cut_%s",extractedValue.c_str(),cutVariable.c_str(), cutVariableValue.c_str());
	string directoryName_2 = directoryName;

	TString cut = "";

	if(extractedValue == "Mean_MmumugammaGen")
	{
		if(eta == "Barrel") directoryName += "/Barrel_";
        	if(eta == "Endcaps") directoryName += "/Endcaps_";
        	if(eta == "all") directoryName += "/BarrelAndEncaps_";

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

		//if(analysisVersion == "2012") cut += " && ( hltnames == \"HLT_Mu17_Mu8_v1\" || hltnames == \"HLT_Mu17_Mu8_v2\" || hltnames == \"HLT_Mu17_Mu8_v3\" || hltnames == \"HLT_Mu17_Mu8_v4\" || hltnames == \"HLT_Mu17_Mu8_v5\" || hltnames == \"HLT_Mu17_Mu8_v6\" || hltnames == \"HLT_Mu17_Mu8_v7\" || hltnames == \"HLT_Mu17_Mu8_v8\" || hltnames == \"HLT_Mu17_Mu8_v9\" || hltnames == \"HLT_Mu17_Mu8_v10\" || hltnames == \"HLT_Mu17_Mu8_v11\" || hltnames == \"HLT_Mu17_Mu8_v12\" || hltnames == \"HLT_Mu17_Mu8_v13\" || hltnames == \"HLT_Mu17_Mu8_v14\" || hltnames == \"HLT_Mu17_Mu8_v15\" || hltnames == \"HLT_Mu17_Mu8_v16\" || hltnames == \"HLT_Mu17_Mu8_v17\" || hltnames == \"HLT_Mu17_Mu8_v18\" || hltnames == \"HLT_Mu17_Mu8_v19\" || hltnames == \"HLT_Mu17_Mu8_v20\" || hltnames == \"HLT_Mu17_Mu8_v21\" || hltnames == \"HLT_Mu17_Mu8_v22\" || hltnames == \"HLT_Mu17_Mu8_v23\") ";	
		if(analysisVersion == "2012") cut += " && isGoodHLT_Mu17_Mu8 == 1";
		//if(analysisVersion == "2012") cut += " && isGoodHLT == 1";

		cut += ")*weight_pileUp*weight_Xsection"; //FIXME ^^ or vv 	

	}

	cout<<endl<<"cut = "<<cut<<endl;

        TChain * chain = new TChain("miniTree");


	if(extractedValue == "Mean_MmumugammaGen") //FIXME
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
                                chain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_totouples_Run2012A_22Jan2013_v1_November2013_0_thesis_v3f_partALL.root");
				chain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_totouples_parked_Run2012B_22Jan2013_v1_November2013_0_thesis_v3f_partALL.root");
				chain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_totouples_parked_Run2012C_22Jan2013_v1_November2013_0_thesis_v3f_partALL.root");
				chain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_totouples_parked_Run2012D_22Jan2013_v1_November2013_0_thesis_v3f_partALL.root");
                        }
                        if(dataType == "MC")
                        {
                                chain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_totouple_DYToMuMu_Summer12_S10_reg5_1_thesis_v3f_partALL.root");
                        	chain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_totouple_DYToMuMu_Summer12_S10_reg5_2_thesis_v3f_partALL.root");
                        	chain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_totouple_TTJets_Summer12_S10_reg5_3_thesis_v3f_partALL.root");
                        	chain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_totouples_WJetsToLNu_Summer12_S10_November2013_3_thesis_v3f_partALL.root");
                        }
                }
	}

	TChain * reducedChain = (TChain *) chain->CopyTree(cut);

        TCanvas * c1 = new TCanvas("c1", "c1",0,0,600,600);

        vector <double> fitParameters;

        double rangeMin = 89.0;
        double rangeMax = 93.0;
        TChain * reducedChain2 = 0;



	if(extractedValue == "Mean_MmumugammaGen")
	{
		TString cut2 = cut;
		int maxPtMu = 35;
		int nCat = 7;	
	
		TH1D * Mean_MmumugammaGen = new TH1D("Mean_MmumugammaGen","Mean_MmumugammaGen",nCat,0,maxPtMu);
		int iter = 0;
		double mean = 0;
		double xMaxHisto = 91.188;
		for(int i = 10; i < maxPtMu; i+=5)
		{
			fileName = Form("fit_%d",i);

			fitParameters.push_back(i);

			cut2 = Form("Photon_Et > %d && Photon_Et < %d && Photon_Et > %d && Photon_Et < %d",i,i+5,i,i+5);
			cout<<endl<<"cut2 = "<<cut2<<endl;	
			reducedChain2 = (TChain *) reducedChain->CopyTree(cut2);
			
			TH1D * h_temp = new TH1D("h_temp","h_temp",100000,0,1000);
			reducedChain2->Draw(Form("%s>>h_temp",xVariable.c_str()));
			mean = h_temp->GetMean();
			xMaxHisto = h_temp->GetXaxis()->GetBinCenter(h_temp->GetMaximumBin());

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
			b.addUniform(110,89.0,93.0);
		
			RooAbsPdf * pdf_fit = NULL;
			if(fitFunction == "voigtian2") 
			{
				fitParameters.push_back(4);
				pdf_fit = voigtian2(dataset, dataset, variable, fitFrame, b, rangeMin, rangeMax, xMaxHisto, fitParameters);
				chiSquare(fitFrame,(char *)"mycurve",(char *)"myhist",fitParameters,3);
				cout << endl << "ICIIIIII ---->>>>> p-val = " << fitParameters[fitParameters.size() - 1]<<endl;
				double MmmgGenFitVal = 0;
				MmmgGenFitVal = fitParameters[iter * 12 + 3];
				cout << endl << "iter = " << iter << ", MmmgGenFitVal = "<<MmmgGenFitVal << endl;
				Mean_MmumugammaGen->SetBinContent(iter+3,MmmgGenFitVal); //FIXME
				double MmmgGenFitValError = 0;
				MmmgGenFitValError = fitParameters[iter * 12 + 4];
				cout << endl << "iter = " << iter << ", MmmgGenFitValError = "<<MmmgGenFitValError << endl;
				Mean_MmumugammaGen->SetBinError(iter+3,MmmgGenFitValError); //FIXME

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
			if(dataType == "data") latexLabel.DrawLatex(0.17, 0.82, Form("%d < P_{T}^{#gamma} < %d",i,i+5));

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

		Mean_MmumugammaGen->Draw("E1");//E1
		Mean_MmumugammaGen->GetYaxis()->SetRangeUser(90.5,92.0);

		TF1 * f1 = new TF1("f1","[0]*x+[1]",0,maxPtMu);
		Mean_MmumugammaGen->Fit(f1);
		f1->Draw("SAMES");
		f1->SetLineWidth(2);
	
		double a = f1->GetParameter(0);
		double b = f1->GetParameter(1);	
	
		Mean_MmumugammaGen->GetYaxis()->SetTitle(xVariableName.c_str());
        	Mean_MmumugammaGen->GetXaxis()->SetTitle("P_{T}^{#gamma} (GeV)");
        	Mean_MmumugammaGen->GetYaxis()->SetTitleOffset(1.55);
        	Mean_MmumugammaGen->GetXaxis()->SetTitleOffset(1.40);
		Mean_MmumugammaGen->SetMarkerStyle(20);
		Mean_MmumugammaGen->SetMarkerSize(0.6);	

		TLatex latexLabel;
		latexLabel.SetTextFont(42);
		latexLabel.SetTextSize(0.028);
		latexLabel.SetNDC();
		
		if(analysisVersion == "2011") latexLabel.DrawLatex(0.25, 0.96, "CMS Private, #sqrt{s} = 7 TeV");
                if(analysisVersion == "2012") latexLabel.DrawLatex(0.25, 0.96, "CMS Private, #sqrt{s} = 8 TeV");    
                latexLabel.SetTextSize(0.030);
                if(dataType == "MC") latexLabel.DrawLatex(0.17, 0.88, "Simulation");
                if(dataType == "data") latexLabel.DrawLatex(0.17, 0.88, Form("Data, #int L = %s fb^{-1}",lumi.c_str()));

		latexLabel.DrawLatex(0.17, 0.82, Form("%f * P_{T}^{#gamma} + %f",a,b));
		
	

		plotsRecording(directoryName, fileName, c1);

		c1->Clear();
		f1->Delete();
		f1 = 0;
		Mean_MmumugammaGen->Delete();
		Mean_MmumugammaGen = 0;


	}
	
	c1->Clear();

	delete c1;
	c1 = 0;
	chain->Delete();
	chain = 0;

	return 0;


}






 
