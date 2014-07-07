#include "functions.h"
#include "setTDRStyle.C"


int main(int argc, char *argv[])
{
        for(int iarg = 0 ; iarg < argc; iarg++)
        {
                cout << "argv[" << iarg << "]= " << argv[iarg] << endl;
        }

        if( argc == 1 ) 
        {
                cerr << "arguments should be passed : directoryName, dataType, eta, r9, xVariable, yVariable, binFileName" <<endl; 
                return 1;

        }

	string directoryName = "SVsNVerticesProfile";
        string dataType = "MC";
	string eta = "Barrel"; 
        string r9 = "all";
	string xVariable = "Photon_Et";
	string yVariable = "GreenTerm";
	string binFileName = "Limites_GreenTerm_6bins.txt";

	if( argc > 1 ) directoryName = argv[1];
        if( argc > 2 ) dataType = argv[2];
	if( argc > 3 ) eta = argv[3];
        if( argc > 4 ) r9 = argv[4];
        if( argc > 5 ) xVariable = argv[5];
	if( argc > 6 ) yVariable = argv[6];
	if( argc > 7 ) binFileName = argv[7];

	gROOT->Reset();
        TGaxis::SetMaxDigits(4);
        setTDRStyle();
        gStyle->SetPalette(1,0);

	
	string fileName = directoryName;
	string lumi = "19.75";

	double xMin, xMax, yMin, yMax;
	xMin = xMax = yMin = yMax = 0;

	string xVariableName = "";
	string yVariableName = "";

	
	if(xVariable == "Photon_Et")
	{
		xMin = 0;
		xMax = 100;
		xVariableName = "P_{T}^{#gamma}";
	}
	
	if(yVariable == "GreenTerm")
        {
                yMin = -10.0;
                yMax = 10.0;
                yVariableName = "(E_{#mu#mu}^{G} - #vec p_{#mu#mu}^{G} #vec u_{#gamma}^{G}) / (E_{#mu#mu}^{R} - #vec p_{#mu#mu}^{R} #vec u_{#gamma}^{R}) - 1 (%)";
        }

	int nBins = rowsNumberInFile(binFileName.c_str()) - 1;

	double * xBinTab = new double[nBins+1];	
	double * xErrorTab = new double[nBins];
        double * xValueTab = new double[nBins];
	double * yErrorTab = new double[nBins];
        double * yValueTab = new double[nBins];

	ifstream binFile(binFileName.c_str());
	for(int i = 0; i <= nBins; i++)
	{
		binFile >> xBinTab[i];
	}
	for(int i = 0; i < nBins; i++)
        {
		xErrorTab[i] = (xBinTab[i+1] - xBinTab[i]) / 2.0;
                xValueTab[i] = (xBinTab[i+1] + xBinTab[i]) / 2.0;
	}	

	TChain * chain = new TChain("miniTree");

	if(dataType == "data")
        {
                chain->Add("../cvs_developpment/Selection_NewMuID/miniTree_Run2012A_13Jul2012_v1_NewMuonID_NewSelection_0_v2_partALL.root");
                chain->Add("../cvs_developpment/Selection_NewMuID/miniTree_Run2012A_recover_06Aug2012_v1_NewMuonID_NewSelection_0_v2_partALL.root");
                chain->Add("../cvs_developpment/Selection_NewMuID/miniTree_Run2012B_13Jul2012_v4_NewMuonID_NewSelection_0_v2_partALL.root");
                chain->Add("../cvs_developpment/Selection_NewMuID/miniTree_Run2012C-24Aug2012-v1_NewMuonID_NewSelection_0_v2_partALL.root");
                chain->Add("../cvs_developpment/Selection_NewMuID/miniTree_Run2012C-EcalRecover_11Dec2012-v1_NewMuonID_NewSelection_0_v2_partALL.root");
                chain->Add("../cvs_developpment/Selection_NewMuID/miniTree_Run2012C_PromptReco_v2_NewMuonID_NewSelection_0_v2_partALL.root");
                chain->Add("../cvs_developpment/Selection_NewMuID/miniTree_Run2012D_PromptReco_v1_NewMuonID_NewSelection_0_v2_partALL.root");
                //chain->Add("../cvs_developpment/Selection_NewMuID/miniTree_2011A_03Oct2011V1ReReco_toto_v2_NewSelection_0_v6_partALL.root");
                //chain->Add("../cvs_developpment/Selection_NewMuID/miniTree_2011A_05Jul2011ReReco_toto_v2_NewSelection_0_v6_partALL.root");
                //chain->Add("../cvs_developpment/Selection_NewMuID/miniTree_2011A_PromptSkimV5ReReco_toto_v2_NewSelection_0_v6_partALL.root");
                //chain->Add("../cvs_developpment/Selection_NewMuID/miniTree_2011B_PromptSkimV1ReReco_toto_v2_NewSelection_0_v6_partALL.root");
        }
        if(dataType == "MC")
        {
                //chain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_totouples_DYToMuMu_Summer12_November2013_1_injRe0_v4_partALL.root");
        	//chain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_totouple_DYToMuMu_Summer12_S10_reg5_1_thesis_v1_partALL.root");     
	  	chain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_totouple_DYToMuMu_Summer12_S10_reg5_1_thesis_v4f_partALL.root");
		//chain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_totouples_DYToMuMu_Summer12_February2014_noskim_1_injRe0_v11_partALL.root"); 
		//chain->Add("../cvs_developpment/Selection_NewMuID/miniTree_TTJets_Summer12_S7_NewMuonID_NewSelection_3_v2_partALL.root");
                //chain->Add("../cvs_developpment/Selection_NewMuID/miniTree_WJetsToLNu_Summer12_S10_NewMuonID_NewSelection_3_v2_partALL.root");
                //chain->Add("../cvs_developpment/Selection_NewMuID/miniTree_DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Fall11-PU_S6_START42_V14B-v1_September12_NewSelection_1_v6_partALL.root");
                //chain->Add("../cvs_developpment/Selection_NewMuID/miniTree_DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Fall11-PU_S6_START42_V14B-v1_September12_NewSelection_2_v6_partALL.root");
                //chain->Add("../cvs_developpment/Selection_NewMuID/miniTree_TTJets_TuneZ2_7TeV-madgraph-tauola_NewSelection_3_v1_partALL.root");
                //chain->Add("../cvs_developpment/Selection_NewMuID/miniTree_WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_NewSelection_3_v1_partALL.root");
        }

 
	TString cut = "isJanLooseMMG == 1";    

        if(r9 == "low" && eta == "Barrel") cut += " && Photon_isEB == 1 && Photon_r9 < 0.94";
        if(r9 == "high" && eta == "Barrel") cut += " && Photon_isEB == 1 && Photon_r9 > 0.94";
        if(r9 == "low" && eta == "Endcaps") cut += " && Photon_isEE == 1 && Photon_r9 < 0.95";
        if(r9 == "high" && eta == "Endcaps") cut += " && Photon_isEE == 1 && Photon_r9 > 0.95";
        if(r9 == "all" && eta == "Barrel") cut += " && Photon_isEB == 1";
        if(r9 == "all" && eta == "Endcaps") cut += " && Photon_isEE == 1";
        if(r9 == "all" && eta == "all") cut += " && (Photon_isEB == 1 || Photon_isEE == 1)";

	//cut += " && ( hltnames == \"HLT_Mu17_TkMu8_v9\" || hltnames == \"HLT_Mu17_TkMu8_v10\" || hltnames == \"HLT_Mu17_TkMu8_v11\" || hltnames == \"HLT_Mu17_TkMu8_v12\" || hltnames == \"HLT_Mu17_TkMu8_v13\" || hltnames == \"HLT_Mu17_TkMu8_v14\" || hltnames == \"HLT_Mu17_TkMu8_v15\" || hltnames == \"HLT_Mu17_TkMu8_v16\" )";

	//cut += " && ( hltnames == \"HLT_Mu17_Mu8_v1\" || hltnames == \"HLT_Mu17_Mu8_v2\" || hltnames == \"HLT_Mu17_Mu8_v3\" || hltnames == \"HLT_Mu17_Mu8_v4\" || hltnames == \"HLT_Mu17_Mu8_v5\" || hltnames == \"HLT_Mu17_Mu8_v6\" || hltnames == \"HLT_Mu17_Mu8_v7\" || hltnames == \"HLT_Mu17_Mu8_v8\" || hltnames == \"HLT_Mu17_Mu8_v9\" || hltnames == \"HLT_Mu17_Mu8_v10\" || hltnames == \"HLT_Mu17_Mu8_v11\" || hltnames == \"HLT_Mu17_Mu8_v12\" || hltnames == \"HLT_Mu17_Mu8_v13\" || hltnames == \"HLT_Mu17_Mu8_v14\" || hltnames == \"HLT_Mu17_Mu8_v15\" || hltnames == \"HLT_Mu17_Mu8_v16\" || hltnames == \"HLT_Mu17_Mu8_v17\" || hltnames == \"HLT_Mu17_Mu8_v18\" || hltnames == \"HLT_Mu17_Mu8_v19\" || hltnames == \"HLT_Mu17_Mu8_v20\" || hltnames == \"HLT_Mu17_Mu8_v21\" || hltnames == \"HLT_Mu17_Mu8_v22\" || hltnames == \"HLT_Mu17_Mu8_v23\")";

        TChain * ReducedChain = (TChain *) chain->CopyTree(cut);

	TCanvas *c1 = new TCanvas("c1", "c1",0,0,600,600);


	//TH2F* th2 = new TH2F("th2", "th2", nBins, &xBinTab[0], 10000,xMin,xMax);
	TH2F* th2 = new TH2F("th2", "th2", nBins, &xBinTab[0], 10000,yMin,yMax);
	if(yVariable == "mmg_s" || yVariable == "mmg_s_true") chain->Draw(Form("%s:%s>>th2",yVariable.c_str(),xVariable.c_str()),cut);
	else if(yVariable == "GreenTerm")
	{
		//string temp_string = "( ( MuonM_MC_E + MuonP_MC_E - ( (MuonM_MC_Px + MuonP_MC_Px)*Photon_MC_Px + (MuonM_MC_Py + MuonP_MC_Py)*Photon_MC_Py + (MuonM_MC_Pz + MuonP_MC_Pz)*Photon_MC_Pz ) * (  1.0 / sqrt(Photon_MC_Px*Photon_MC_Px + Photon_MC_Py*Photon_MC_Py + Photon_MC_Pz*Photon_MC_Pz) )  ) * 100   /   ( MuonM_E + MuonP_E - ( (MuonM_Px + MuonP_Px)*Photon_Px + (MuonM_Py + MuonP_Py)*Photon_Py + (MuonM_Pz + MuonP_Pz)*Photon_Pz ) * (  1.0 / sqrt(Photon_Px*Photon_Px + Photon_Py*Photon_Py + Photon_Pz*Photon_Pz) )  ) - 1)";

		string temp_string = "((GreenTerm - 1) * 100)";
		chain->Draw(Form("%s:%s>>th2",temp_string.c_str(),xVariable.c_str()),cut);	

	}	
	else chain->Draw(Form("%s:%s>>th2",yVariable.c_str(),xVariable.c_str()),cut);	
	TProfile* pro = (TProfile*)(th2->ProfileX());

	pro->Draw("");

	//c1->Clear();

	//pro->Draw("");
/*
	yMin = 10000;
	yMax = -10000;	
	for(int i = 1; i <= nBins; i++)
        {
        	if(yMin > pro->GetBinContent(i)) yMin = pro->GetBinContent(i);
		if(yMax < pro->GetBinContent(i)) yMax = pro->GetBinContent(i);
	}
*/


	pro->GetXaxis()->SetTitle(Form("%s",xVariableName.c_str()));
	pro->GetYaxis()->SetTitle(Form("%s",yVariableName.c_str()));
	//pro->GetYaxis()->SetTitleOffset(1.6);

	pro->GetXaxis()->SetLabelFont(42);
        pro->GetXaxis()->SetTitleFont(42);
        pro->GetYaxis()->SetLabelFont(42);
        pro->GetYaxis()->SetTitleFont(42);
        pro->GetYaxis()->SetTitleOffset(1.55);
        pro->GetXaxis()->SetTitleOffset(1.40);

        pro->SetLineColor(kGreen+1);
        //pro->SetLineColor(kBlue);
	pro->SetMarkerStyle(20);
        pro->SetMarkerSize(0.6);
        pro->SetMarkerColor(kGreen+1);
        //pro->SetMarkerColor(kBlue);
	//pro->GetYaxis()->SetRangeUser(yMin-0.07,yMax+0.07);
	pro->GetXaxis()->SetLimits(xMin,xMax);

	pro->GetYaxis()->SetRangeUser(-0.5,0.5);
	
	TLatex latexLabel;
	//latexLabel.SetTextSize(0.030);
        //latexLabel.SetNDC();
        latexLabel.SetTextFont(42);
        latexLabel.SetTextSize(0.028);
        latexLabel.SetNDC();

	//latexLabel.DrawLatex(0.13, 0.96, "CMS Private 2012, #sqrt{s} = 8 TeV");
     
        //if(dataType == "MC") latexLabel.DrawLatex(0.17, 0.88, "Simulation");
        //if(dataType == "data") latexLabel.DrawLatex(0.17, 0.88, Form("Data, #int L = %s fb^{-1}",lumi.c_str()));
        //latexLabel.DrawLatex(0.25, 0.96, "CMS 2012                      #sqrt{s} = 8 TeV                        Simulation");
	//latexLabel.DrawLatex(0.25, 0.96, "CMS Preliminary 2012                #sqrt{s} = 8 TeV                Simulation");
	latexLabel.DrawLatex(0.25, 0.96, "CMS Private 2012                 #sqrt{s} = 8 TeV                 Simulation");
	latexLabel.SetTextSize(0.030);	
	latexLabel.DrawLatex(0.17, 0.90,Form("ECAL %s",eta.c_str()));   
        if(r9 == "low" && eta == "Barrel") latexLabel.DrawLatex(0.17, 0.85,"r9 < 0,94");
        if(r9 == "low" && eta == "Endcaps") latexLabel.DrawLatex(0.17, 0.85,"r9 < 0,95");
        if(r9 == "high" && eta == "Barrel") latexLabel.DrawLatex(0.17, 0.85,"r9 > 0,94");
        if(r9 == "high" && eta == "Endcaps") latexLabel.DrawLatex(0.17, 0.85,"r9 > 0,95");
	//if(r9 == "all" && eta == "Barrel") latexLabel.DrawLatex(0.17, 0.78,"All r9");
        //if(r9 == "all" && eta == "Endcaps") latexLabel.DrawLatex(0.17, 0.78,"All r9");
	if(r9 == "all") latexLabel.DrawLatex(0.17, 0.85,"All r9");


	directoryName += Form("/%s/%s_%sR9/",dataType.c_str(),eta.c_str(),r9.c_str());
	plotsRecording(directoryName, fileName, c1);	


	c1->Clear();

	th2->Delete();
	th2 = 0;
	chain->Delete();
	chain = 0;
	delete [] xBinTab;
	xBinTab = 0;
	delete [] xErrorTab;
	xErrorTab = 0;	
	delete [] xValueTab;
	xValueTab = 0;
	delete [] yErrorTab;
        yErrorTab = 0;
        delete [] yValueTab;
        yValueTab = 0;


/*

		if(EndCaps == 0 && r9sup == 0)
                {
                        for(int a = 1; a <= 6; a++)
                        {
                                r9infEBTab[a-1] = pro->GetBinContent(a);
                                r9infEBErrorTab[a-1] = pro->GetBinError(a);
                                //cout<<endl<<"a = "<<a<<endl<<"r9infEBTab[a] = "<<r9infEBTab[a-1]<<endl<<"r9infEBErrorTab[a] = "<<r9infEBErrorTab[a-1];
                        }
                }


        TMultiGraph *mg = new TMultiGraph();

        TGraphAsymmErrors * ErecoOverEtruer9infEB = new TGraphAsymmErrors(nBins,xValueTab, r9infEBTab, xErrorTab, xErrorTab, r9infEBErrorTab, r9infEBErrorTab);

        TGraphAsymmErrors * ErecoOverEtruer9supEB = new TGraphAsymmErrors(nBins,xValueTab, r9supEBTab, xErrorTab, xErrorTab, r9supEBErrorTab, r9supEBErrorTab);

        TGraphAsymmErrors * ErecoOverEtruer9infEE = new TGraphAsymmErrors(nBins,xValueTab, r9infEETab, xErrorTab, xErrorTab, r9infEEErrorTab, r9infEEErrorTab);

        TGraphAsymmErrors * ErecoOverEtruer9supEE = new TGraphAsymmErrors(nBins,xValueTab, r9supEETab, xErrorTab, xErrorTab, r9supEEErrorTab, r9supEEErrorTab);


        c1->ToggleEventStatus();
        //ErecoOverEtruer9infEB->Draw("AP");
        //ErecoOverEtruer9supEB->Draw("SAMES");
        //ErecoOverEtruer9infEE->Draw("SAMES");
        //ErecoOverEtruer9supEE->Draw("SAMES");

        mg->Add(ErecoOverEtruer9infEB);
        mg->Add(ErecoOverEtruer9supEB);
        mg->Add(ErecoOverEtruer9infEE);
        mg->Add(ErecoOverEtruer9supEE);
        mg->Draw("AP");


        mg->GetXaxis()->SetTitle("P_{T}^{#gamma}");

        mg->GetXaxis()->SetLabelFont(42);
        mg->GetXaxis()->SetTitleFont(42);
        mg->GetXaxis()->SetLabelSize(0.03);

        mg->GetYaxis()->SetTitle("s_{RECO} = 1/k_{RECO} - 1 (%)");
	if(version == "reco") mg->GetYaxis()->SetTitle("s_{RECO} = 1/k_{RECO} - 1 (%)");
        if(version == "recoSurface") mg->GetYaxis()->SetTitle("s_{RECO Surface} = 1/k_{RECO Surface} - 1 (%)");
        if(version == "true") mg->GetYaxis()->SetTitle("s_{TRUE} = E_{RECO}/E_{TRUE} - 1 (%)");


        mg->GetYaxis()->SetLabelFont(42);
        mg->GetYaxis()->SetTitleOffset(1.24);
        mg->GetYaxis()->SetTitleFont(42);
        mg->GetYaxis()->SetLabelSize(0.03);

        ErecoOverEtruer9infEB->SetLineColor(629);
        ErecoOverEtruer9supEB->SetLineColor(597);
        ErecoOverEtruer9infEE->SetLineColor(613);
        ErecoOverEtruer9supEE->SetLineColor(429);
	//413
        //397
        ErecoOverEtruer9infEB->SetMarkerColor(629);
        ErecoOverEtruer9supEB->SetMarkerColor(597);
        ErecoOverEtruer9infEE->SetMarkerColor(613);
        ErecoOverEtruer9supEE->SetMarkerColor(429);
        //413
        //397   

        ErecoOverEtruer9infEB->SetMarkerStyle(20);
        ErecoOverEtruer9supEB->SetMarkerStyle(24);
        ErecoOverEtruer9infEE->SetMarkerStyle(21);
        ErecoOverEtruer9supEE->SetMarkerStyle(25);


        ErecoOverEtruer9infEB->SetName("ErecoOverEtruer9infEB");
        ErecoOverEtruer9supEB->SetName("ErecoOverEtruer9supEB");
        ErecoOverEtruer9infEE->SetName("ErecoOverEtruer9infEE");
        ErecoOverEtruer9supEE->SetName("ErecoOverEtruer9supEE");

        text = new TLatex();
        text->SetNDC();
        text->SetTextAlign(11);
        text->SetTextFont(42);
        text->SetTextSizePixels(17);
        text->SetTextSize(0.028);
        text->DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
        text->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV");
        text->DrawLatex(0.16, 0.80, "TProfile");
        text->DrawLatex(0.16, 0.75, "ETHZ Corrections");

        TLegend * leg = new TLegend(0.6,0.7,0.9,0.9,"","brNDC");
        //leg->SetHeader("The Legend Title");
        //leg->AddEntry(h1,"Histogram filled with random numbers","f");
        //leg->AddEntry("f1","Function abs(#frac{sin(x)}{x})","l");
        leg->SetTextSize(0.030);
        leg->SetFillColor(kWhite);
        leg->SetLineColor(kWhite);
        leg->SetShadowColor(kWhite);
        leg->AddEntry(ErecoOverEtruer9infEB->GetName(),"Barrel R_{9 #gamma} < 0.94","lep");
        leg->AddEntry(ErecoOverEtruer9supEB->GetName(),"Barrel R_{9 #gamma} > 0.94","lep");
        leg->AddEntry(ErecoOverEtruer9infEE->GetName(),"Endcaps R_{9 #gamma} < 0.95","lep");
        leg->AddEntry(ErecoOverEtruer9supEE->GetName(),"Endcaps R_{9 #gamma} > 0.95","lep");
        leg->Draw();

	gStyle->SetPadBorderMode(0);

        mg->GetYaxis()->SetRangeUser(yminErecoOverEtrue,ymaxErecoOverEtrue);
        mg->GetXaxis()->SetLimits(xminErecoOverEtrue,xmaxErecoOverEtrue);

        c1->SetTickx(1);
        c1->SetTicky(1);
        c1->SetGridx(1);
        c1->SetGridy(1);
        c1->Modified();
        c1->cd();
        c1->SetSelected(c1);
        c1->ToggleToolBar();

        //nomDossier = Form("SrecoVsPt/%s/2011/MmumuSelection%s_%s/WithCracks/%s/",variableX.c_str(),LowMmumuLim.c_str(),HightMmumuLim.c_str(),SetOfCorrections.c_str());

        nomFichier = Form("SrecoVsPt_4Cat_%s",version.c_str());

        enregistrementPlots(nomDossier, nomFichier, 2, 10000, c1);
*/	
	return 0;

}
