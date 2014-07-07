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
	string cutVariable = "Photon_Et";
	string cutVariableValue = "0";
	string xVariable = "Photon_E";
	string yVariable = "HalfGreenTerm";
	int doMC = 0;

	if( argc > 1 ) directoryName = argv[1];
	if( argc > 2 ) eta = argv[2];
        if( argc > 3 ) r9 = argv[3];	
	if( argc > 4 ) cutVariable = argv[4];
	if( argc > 5 ) cutVariableValue = argv[5];	
	if( argc > 6 ) xVariable = argv[6];
	if( argc > 7 ) yVariable = argv[7];	
	if( argc > 8 ) 
        {
                std::stringstream ss ( argv[8] );
                ss >> doMC;
        }


	
	double xMin,xMax,yMin,yMax;
	int nBinsX, nBinsY;
	string yVariableName, xVariableName;

	if(xVariable == "Photon_E")
        {
		nBinsX = 300;
		xMin = 0;
		xMax = 300;
		xVariableName = "E_{#gamma}^{R}";
	}
	
	if(xVariable == "Photon_MC_E")
        {
                nBinsX = 300;
                xMin = 0;
                xMax = 300;
                xVariableName = "E_{#gamma}^{G}";
        }
	
	if(yVariable == "HalfGreenTerm")
	{
		nBinsY = 300;
		yMin = 0;
		yMax = 300;
		if(doMC == 1)
		{
			yVariableName = "(M_{#mu#mu#gamma}^{2 G} - M_{#mu#mu}^{2 G}) / 2 (E_{#mu#mu}^{G} - #vec p_{#mu#mu}^{G} #vec u_{#gamma}^{G})";
		}
		else
		{
			yVariableName = "(M_{#mu#mu#gamma}^{2 R} - M_{#mu#mu}^{2 R}) / 2 (E_{#mu#mu}^{R} - #vec p_{#mu#mu}^{R} #vec u_{#gamma}^{R})";
		}
	}

	double xMinLeg,yMinLeg,xMaxLeg,yMaxLeg,legTextSize;
	
        xMinLeg = 0.20;
       	xMaxLeg = 0.44;
	yMinLeg = 0.74;
	yMaxLeg = 0.94;
	legTextSize = 0.028;


	gROOT->Reset();
        TGaxis::SetMaxDigits(3);
        setTDRStyle();
	//setEgammaStyle();

	string fileName = directoryName;
	directoryName += Form("/%s_cut_%s",cutVariable.c_str(), cutVariableValue.c_str());
	string directoryName_2 = directoryName;

	if(eta == "Barrel") directoryName += "/Barrel_";
        if(eta == "Endcaps") directoryName += "/Endcaps_";
        if(eta == "all") directoryName += "/BarrelAndEncaps_";

        if(r9 == "low") directoryName += "lowR9/";
        if(r9 == "high") directoryName += "highR9/";
        if(r9 == "all") directoryName += "AllR9/";

	TString cut = Form("%s > %s && isAfterFSRCut1",cutVariable.c_str(), cutVariableValue.c_str());

	if(r9 == "low" && eta == "Barrel") cut += " && Photon_isEB == 1 && Photon_r9 < 0.94";
        if(r9 == "high" && eta == "Barrel") cut += " && Photon_isEB == 1 && Photon_r9 > 0.94";
        if(r9 == "low" && eta == "Endcaps") cut += " && Photon_isEE == 1 && Photon_r9 < 0.95";
        if(r9 == "high" && eta == "Endcaps") cut += " && Photon_isEE == 1 && Photon_r9 > 0.95";
        if(r9 == "all" && eta == "Barrel") cut += " && Photon_isEB == 1";
        if(r9 == "all" && eta == "Endcaps") cut += " && Photon_isEE == 1";
	if(r9 == "all" && eta == "all") cut += " && (Photon_isEE == 1 || Photon_isEB == 1)";

	//cut += " && ( hltnames == \"HLT_Mu17_TkMu8_v9\" || hltnames == \"HLT_Mu17_TkMu8_v10\" || hltnames == \"HLT_Mu17_TkMu8_v11\" || hltnames == \"HLT_Mu17_TkMu8_v12\" || hltnames == \"HLT_Mu17_TkMu8_v13\" || hltnames == \"HLT_Mu17_TkMu8_v14\" || hltnames == \"HLT_Mu17_TkMu8_v15\" || hltnames == \"HLT_Mu17_TkMu8_v16\" ) "; //FIXME	

	//cut += " && ( hltnames == \"HLT_Mu17_Mu8_v1\" || hltnames == \"HLT_Mu17_Mu8_v2\" || hltnames == \"HLT_Mu17_Mu8_v3\" || hltnames == \"HLT_Mu17_Mu8_v4\" || hltnames == \"HLT_Mu17_Mu8_v5\" || hltnames == \"HLT_Mu17_Mu8_v6\" || hltnames == \"HLT_Mu17_Mu8_v7\" || hltnames == \"HLT_Mu17_Mu8_v8\" || hltnames == \"HLT_Mu17_Mu8_v9\" || hltnames == \"HLT_Mu17_Mu8_v10\" || hltnames == \"HLT_Mu17_Mu8_v11\" || hltnames == \"HLT_Mu17_Mu8_v12\" || hltnames == \"HLT_Mu17_Mu8_v13\" || hltnames == \"HLT_Mu17_Mu8_v14\" || hltnames == \"HLT_Mu17_Mu8_v15\" || hltnames == \"HLT_Mu17_Mu8_v16\" || hltnames == \"HLT_Mu17_Mu8_v17\" || hltnames == \"HLT_Mu17_Mu8_v18\" || hltnames == \"HLT_Mu17_Mu8_v19\" || hltnames == \"HLT_Mu17_Mu8_v20\" || hltnames == \"HLT_Mu17_Mu8_v21\" || hltnames == \"HLT_Mu17_Mu8_v22\" || hltnames == \"HLT_Mu17_Mu8_v23\")";

        TChain * dataChain = new TChain("miniTree");
	TChain * dYToMuMuFSRChain = new TChain("miniTree");
	TChain * dYToMuMuNonFSRChain = new TChain("miniTree");
	TChain * ttJetsChain = new TChain("miniTree");
	TChain * wJetsChain = new TChain("miniTree");

	/*
	dataChain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_totouples_Run2012A_22Jan2013_v1_November2013_0_injRe0_v1_partALL.root");
        dataChain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_totouples_parked_Run2012B_22Jan2013_v1_November2013_0_injRe0_v1_partALL.root");
        dataChain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_totouples_parked_Run2012C_22Jan2013_v1_November2013_0_injRe0_v1_partALL.root");
        dataChain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_totouples_parked_Run2012D_22Jan2013_v1_November2013_0_injRe0_v1_partALL.root");

	dYToMuMuFSRChain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_totouples_DYToMuMu_Summer12_February2014_noskim_1_injRe0_v11_partALL.root");
	dYToMuMuNonFSRChain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_totouples_DYToMuMu_Summer12_November2013_2_injRe0_v4_partALL.root");
	ttJetsChain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_totouples_TTJets_FullLeptMGDecays_Summer12_November2013_v2_3_injRe0_v4_partALL.root");
	ttJetsChain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_totouples_TTJets_SemiLeptMGDecays_Summer12_November2013_3_injRe0_v4_partALL.root");
	ttJetsChain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_totouples_TTJets_HadronicMGDecays_Summer12_November2013_v2_temp218_3_injRe0_v4_partALL.root");
	wJetsChain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_totouples_WJetsToLNu_Summer12_S10_November2013_3_injRe0_v4_partALL.root");	
	*/

	dataChain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_totouples_Run2012A_22Jan2013_v1_November2013_0_thesis_v1_partALL.root");
        dataChain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_totouples_parked_Run2012B_22Jan2013_v1_November2013_0_thesis_v1_partALL.root");
        dataChain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_totouples_parked_Run2012C_22Jan2013_v1_November2013_0_thesis_v1_partALL.root");
        dataChain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_totouples_parked_Run2012D_22Jan2013_v1_November2013_0_thesis_v1_partALL.root");

        dYToMuMuFSRChain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_totouple_DYToMuMu_Summer12_S10_reg5_1_thesis_v1_partALL.root");
        dYToMuMuNonFSRChain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_totouple_DYToMuMu_Summer12_S10_reg5_2_thesis_v1_partALL.root");
        ttJetsChain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_totouple_TTJets_Summer12_S10_reg5_3_thesis_v1_partALL.root");
        wJetsChain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_totouples_WJetsToLNu_Summer12_S10_November2013_3_thesis_v1_partALL.root");

	
	TCanvas *c1 = new TCanvas("c1", "c1",0,0,600,600);

	//TH2D *data = new TH2D("data","data", nBinsX, xMin, xMax,nBinsY,yMin,yMax);	
	TH2D *dYToMuMuFSR = new TH2D("dYToMuMuFSR","dYToMuMuFSR", nBinsX, xMin, xMax,nBinsY,yMin,yMax);
	TH2D *dYToMuMuNonFSR = new TH2D("dYToMuMuNonFSR","dYToMuMuNonFSR", nBinsX, xMin, xMax,nBinsY,yMin,yMax);
	TH2D *ttJets = new TH2D("ttJets","ttJets", nBinsX, xMin, xMax,nBinsY,yMin,yMax);
	TH2D *wJets = new TH2D("wJets","wJets", nBinsX, xMin, xMax,nBinsY,yMin,yMax);

	/*TH1D *dYToMuMuFSR = new TH1D("dYToMuMuFSR","dYToMuMuFSR", nBins, xMin, xMax);
	TH1D *dYToMuMuNonFSR = new TH1D("dYToMuMuNonFSR","dYToMuMuNonFSR", nBins, xMin, xMax);
	TH1D *ttJets = new TH1D("ttJets","ttJets", nBins, xMin, xMax);
	TH1D *wJets = new TH1D("wJets","wJets", nBins, xMin, xMax);				
	*/
	
	string temp_string = yVariable;
	if(yVariable == "HalfGreenTerm") 
	{
		if(doMC == 1)
		{
			//temp_string = "(Mmumugamma_MMG_MC * Mmumugamma_MMG_MC - Mmumu_MMG_MC * Mmumu_MMG_MC) / (2 * ( MuonM_MC_E + MuonP_MC_E - ( (MuonM_MC_Px + MuonP_MC_Px)*Photon_MC_Px + (MuonM_MC_Py + MuonP_MC_Py)*Photon_MC_Py + (MuonM_MC_Pz + MuonP_MC_Pz)*Photon_MC_Pz ) * (  1.0 / sqrt(Photon_MC_Px*Photon_MC_Px + Photon_MC_Py*Photon_MC_Py + Photon_MC_Pz*Photon_MC_Pz) ) ) )";
			temp_string = "HalfGreenTermGEN";
		}
		else
		{
			//temp_string = "(Mmumugamma * Mmumugamma - Mmumu * Mmumu) / (2 * ( MuonM_E + MuonP_E - ( (MuonM_Px + MuonP_Px)*Photon_Px + (MuonM_Py + MuonP_Py)*Photon_Py + (MuonM_Pz + MuonP_Pz)*Photon_Pz ) * (  1.0 / sqrt(Photon_Px*Photon_Px + Photon_Py*Photon_Py + Photon_Pz*Photon_Pz) ) ) )";
			temp_string = "HalfGreenTermRECO";
		}
	}
	//dataChain->Draw("Mmumugamma:Mmumu>>data",cut);
	dYToMuMuFSRChain->Draw(Form("%s:%s>>dYToMuMuFSR",temp_string.c_str(),xVariable.c_str()),cut);
	//dYToMuMuNonFSRChain->Draw(Form("%s:%s>>dYToMuMuNonFSR",temp_string.c_str(),xVariable.c_str()),cut);	
	//ttJetsChain->Draw(Form("%s:%s>>ttJets",temp_string.c_str(),xVariable.c_str()),cut);
	//wJetsChain->Draw(Form("%s:%s>>wJets",temp_string.c_str(),xVariable.c_str()),cut);

/*dYToMuMuFSRChain->Draw(Form("%s>>dYToMuMuFSR",xVariable.c_str()),cut);
	dYToMuMuNonFSRChain->Draw(Form("%s>>dYToMuMuNonFSR",xVariable.c_str()),cut);
	ttJetsChain->Draw(Form("%s>>ttJets",xVariable.c_str()),cut);
	wJetsChain->Draw(Form("%s>>wJets",xVariable.c_str()),cut);	
*/

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


	// --- Weight by integral --- //
/*	
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
*/

	c1->Clear();

	//double meanEtData = data->GetMean();
        //double meanEtMC = dYToMuMuFSR->GetMean();
	
	//data->SetMarkerStyle(20);
	//data->SetMarkerSize(0.5);

	//dYToMuMuFSR->GetXaxis()->SetTitle(Form("%s",xVariableName.c_str()));
	//dYToMuMuFSR->GetYaxis()->SetTitle(Form("%s",yVariableName.c_str()));
	dYToMuMuFSR->GetXaxis()->SetTitle(xVariableName.c_str());
	dYToMuMuFSR->GetYaxis()->SetTitle(yVariableName.c_str());		
	dYToMuMuFSR->GetYaxis()->SetTitleOffset(1.4);

	
	dYToMuMuFSR->GetXaxis()->SetLabelFont(42);
	dYToMuMuFSR->GetXaxis()->SetTitleFont(42);
	dYToMuMuFSR->GetYaxis()->SetLabelFont(42);
	dYToMuMuFSR->GetYaxis()->SetTitleFont(42);
	dYToMuMuFSR->GetYaxis()->SetTitleOffset(1.55);
	dYToMuMuFSR->GetXaxis()->SetTitleOffset(1.40);

	
	//dYToMuMuNonFSR->Draw("");
	//ttJets->Draw("SAMES");
	//wJets->Draw("SAMES");
	dYToMuMuFSR->Draw("");
	
	
	//ttJets->Draw("");
        //wJets->Draw("SAMES");
	/*
	ttJets->GetXaxis()->SetTitle("M_{#mu#mu}");
        ttJets->GetYaxis()->SetTitle("M_{#mu#mu#gamma}");
        ttJets->GetYaxis()->SetTitleOffset(1.4);
	*/
	//data->Draw("E1SAMES");
	

	//dYToMuMuFSR->SetFillColor(2);
   	dYToMuMuFSR->SetFillColor(kGreen-7);
	dYToMuMuFSR->SetMarkerColor(kGreen-7);
	if(doMC == 1)
	{
		dYToMuMuFSR->SetFillColor(kBlue);
        	dYToMuMuFSR->SetMarkerColor(kBlue);	
	}
	else
	{
		dYToMuMuFSR->SetFillColor(kRed);
                dYToMuMuFSR->SetMarkerColor(kRed);
	}
	//dYToMuMuFSR->SetFillStyle(3001);
	//dYToMuMuNonFSR->SetFillColor(3);
	dYToMuMuNonFSR->SetFillColor(kAzure-5);
	dYToMuMuNonFSR->SetMarkerColor(kAzure-5);
        //dYToMuMuNonFSR->SetFillStyle(3001);
	//ttJets->SetFillColor(4);
	ttJets->SetFillColor(kBlue-1);
	ttJets->SetMarkerColor(kBlue-1);
        //ttJets->SetFillStyle(3001);
	//wJets->SetFillColor(5);
	wJets->SetFillColor(kCyan+2);
        wJets->SetMarkerColor(kCyan+2);
	//wJets->SetFillStyle(3001);

//	data->SetName("data");
	dYToMuMuFSR->SetName("dYToMuMuFSR");
	dYToMuMuNonFSR->SetName("dYToMuMuNonFSR");
	ttJets->SetName("ttJets");
	wJets->SetName("wJets");

	TLegend leg(xMinLeg,yMinLeg,xMaxLeg,yMaxLeg,"","brNDC");
	leg.SetTextFont(42);
        leg.SetTextSize(legTextSize);
        leg.SetFillColor(kWhite);
        leg.SetLineColor(kWhite);
        leg.SetShadowColor(kWhite);
        //leg.AddEntry(data->GetName(),"data","lep");
	//leg.AddEntry(dYToMuMuFSR->GetName(),"Z#mu#mu + #gamma FSR","f");
        //leg.AddEntry(dYToMuMuNonFSR->GetName(),"Z#mu#mu + #gamma non FSR","f");
        //leg.AddEntry(ttJets->GetName(),"t#bar{t} + jets","f");
        //leg.AddEntry(wJets->GetName(),"W + jets","f");
	//leg.Draw();

	TLatex latexLabel;
	latexLabel.SetTextFont(42);
	latexLabel.SetTextSize(0.028);
	latexLabel.SetNDC();
	//latexLabel.DrawLatex(0.25, 0.96, "CMS             2012               #sqrt{s} = 8 TeV               L = 19.6 fb^{-1}");
	//latexLabel.DrawLatex(0.25, 0.96, "CMS 2012                     #sqrt{s} = 8 TeV                       Simulation");	
        //latexLabel.DrawLatex(0.25, 0.96, "CMS Preliminary 2012                #sqrt{s} = 8 TeV                Simulation");
	latexLabel.DrawLatex(0.25, 0.96, "CMS Private 2012                 #sqrt{s} = 8 TeV                 Simulation");

	plotsRecording(directoryName, fileName, c1);


	c1->Clear();

	//data->Delete();
        //data = 0;
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






 
