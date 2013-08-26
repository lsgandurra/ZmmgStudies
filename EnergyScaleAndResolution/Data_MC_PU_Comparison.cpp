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

	string directoryName = "Data_MC_PU_Comparison";
	string eta = "all";
	string r9 = "all";
	string xVariable = "nVertices";
	int log = 0;	
	string normalization = "integral"; // "integral"	
	string cutVariable = "Photon_Et";
	string cutVariableValue = "25";

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


	int nBins = 30;
	double xMin, xMax;
	double xMinLeg, xMaxLeg;

	xMinLeg = 0.67;
	xMaxLeg = 0.91;	
	string xVariableName, yVariableName;

	if(xVariable == "nVertices")
        {
                nBins = 60; 
                xMin = 0;
                xMax = 60;
                xVariableName = "nVertices";
                yVariableName = "Events / 1";
                //xMinLeg = 0.20;
                //xMaxLeg = 0.44;
        }	
	
	if(xVariable == "weight_pileUp")
        {
                nBins = 60;
                xMin = 0;
                xMax = 15;
                xVariableName = "weight_pileUp";
                yVariableName = "Events / 0.25";
                //xMinLeg = 0.20;
                //xMaxLeg = 0.44;
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
	TString cut2 = Form("(%s > %s",cutVariable.c_str(), cutVariableValue.c_str());

	if(r9 == "low" && eta == "Barrel") 
	{
		cut += " && Photon_isEB == 1 && Photon_r9 < 0.94";
        	cut2 += " && Photon_isEB == 1 && Photon_r9 < 0.94";
	}
	if(r9 == "high" && eta == "Barrel") 
	{
		cut += " && Photon_isEB == 1 && Photon_r9 > 0.94";
		cut2 += " && Photon_isEB == 1 && Photon_r9 > 0.94";
	}
        if(r9 == "low" && eta == "Endcaps") 
	{
		cut += " && Photon_isEE == 1 && Photon_r9 < 0.95";
		cut2 += " && Photon_isEE == 1 && Photon_r9 < 0.95";
	}
        if(r9 == "high" && eta == "Endcaps") 
	{
		cut += " && Photon_isEE == 1 && Photon_r9 > 0.95";
		cut2 += " && Photon_isEE == 1 && Photon_r9 < 0.95";
	}
        if(r9 == "all" && eta == "Barrel") 
	{
		cut += " && Photon_isEB == 1";
		cut2 += " && Photon_isEB == 1";
	}
        if(r9 == "all" && eta == "Endcaps") 
	{
		cut += " && Photon_isEE == 1";
		cut2 += " && Photon_isEB == 1";
	}
	if(r9 == "all" && eta == "all") 
	{
		cut += " && (Photon_isEE == 1 || Photon_isEB == 1)";
		cut2 += " && (Photon_isEE == 1 || Photon_isEB == 1)";
	}

	if(xVariable != "weight_pileUp")
	{

		if(normalization == "lumi") 
		{
			cut += ")*weight_pileUp*weight_Xsection"; 
			cut2 += ")*weight_pileUp*weight_Xsection";
		}
		if(normalization == "integral") 
		{
			cut += ")*weight_pileUp"; 
			cut2 += ")*weight_pileUp";
		}
		if(normalization == "lumi2") 
		{
			cut += ")*weight_Xsection";
			cut2 += ")*weight_Xsection";
		}
		if(normalization == "integral2") 
		{
			cut += ")";
			cut2 += ")";
		}
	
	}
	else
	{
		cut += ")";
		cut2 += ")";
	}

	//if(normalization == "lumi") cut += ")*weight_Xsection";

	cout<<endl<<"cut = "<<cut<<endl;

        TChain * dataChain = new TChain("miniTree");
	TChain * dYToMuMuFSRChain = new TChain("miniTree");
	TChain * dYToMuMuNonFSRChain = new TChain("miniTree");
	TChain * ttJetsChain = new TChain("miniTree");
	TChain * wJetsChain = new TChain("miniTree");

	dataChain->Add("/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/cvs_developpment/Selection_NewMuID/miniTree_Run2012A_13Jul2012_v1_NewMuonID_NewSelection_0_injRe0.00_v1_partALL.root");
	dataChain->Add("/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/cvs_developpment/Selection_NewMuID/miniTree_Run2012A_recover_06Aug2012_v1_NewMuonID_NewSelection_0_injRe0.00_v1_partALL.root");
	dataChain->Add("/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/cvs_developpment/Selection_NewMuID/miniTree_Run2012B_13Jul2012_v4_NewMuonID_NewSelection_0_injRe0.00_v1_partALL.root");
	dataChain->Add("/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/cvs_developpment/Selection_NewMuID/miniTree_Run2012C-24Aug2012-v1_NewMuonID_NewSelection_0_injRe0.00_v1_partALL.root");
	dataChain->Add("/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/cvs_developpment/Selection_NewMuID/miniTree_Run2012C-EcalRecover_11Dec2012-v1_NewMuonID_NewSelection_0_injRe0.00_v1_partALL.root");
	dataChain->Add("/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/cvs_developpment/Selection_NewMuID/miniTree_Run2012C_PromptReco_v2_NewMuonID_NewSelection_0_injRe0.00_v1_partALL.root");
	dataChain->Add("/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/cvs_developpment/Selection_NewMuID/miniTree_Run2012D_PromptReco_v1_NewMuonID_NewSelection_0_injRe0.00_v1_partALL.root");
	dYToMuMuFSRChain->Add("/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/cvs_developpment/Selection_NewMuID/miniTree_DYToMuMu_Summer12_NewMuonID_NewSelection_1_injRe0.00_v1_partALL.root");
	dYToMuMuNonFSRChain->Add("/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/cvs_developpment/Selection_NewMuID/miniTree_DYToMuMu_Summer12_NewMuonID_NewSelection_2_injRe0.00_v1_partALL.root");
	ttJetsChain->Add("/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/cvs_developpment/Selection_NewMuID/miniTree_TTJets_Summer12_S7_NewMuonID_NewSelection_3_injRe0.00_v1_partALL.root");
	wJetsChain->Add("/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/cvs_developpment/Selection_NewMuID/miniTree_WJetsToLNu_Summer12_S10_NewMuonID_NewSelection_3_injRe0.00_v1_partALL.root");
	//dYToMuMuFSRChain->Add("/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/cvs_developpment/Selection_NewMuID/miniTree_DYToMuMu_Summer12_NewMuonID_NewSelection_1_injRe0.00_v2_NewPU_partALL.root");
	//dYToMuMuNonFSRChain->Add("/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/cvs_developpment/Selection_NewMuID/miniTree_DYToMuMu_Summer12_NewMuonID_NewSelection_2_injRe0.00_v2_NewPU_partALL.root");
	//ttJetsChain->Add("/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/cvs_developpment/Selection_NewMuID/miniTree_TTJets_Summer12_S7_NewMuonID_NewSelection_3_injRe0.00_v2_NewPU_partALL.root");
	//wJetsChain->Add("/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/cvs_developpment/Selection_NewMuID/miniTree_WJetsToLNu_Summer12_S10_NewMuonID_NewSelection_3_injRe0.00_v2_NewPU_partALL.root");


	TCanvas *c1 = new TCanvas("c1", "c1",0,0,600,600);
	if(log == 1) c1->SetLogy();

	TH1D *data = new TH1D("data","data", nBins, xMin, xMax);	
	TH1D *dYToMuMuFSR = new TH1D("dYToMuMuFSR","dYToMuMuFSR", nBins, xMin, xMax);
	TH1D *dYToMuMuNonFSR = new TH1D("dYToMuMuNonFSR","dYToMuMuNonFSR", nBins, xMin, xMax);
	TH1D *ttJets = new TH1D("ttJets","ttJets", nBins, xMin, xMax);
	TH1D *wJets = new TH1D("wJets","wJets", nBins, xMin, xMax);				
	
	TH1D *data2 = new TH1D("data2","data2", nBins, xMin, xMax);
        TH1D *dYToMuMuFSR2 = new TH1D("dYToMuMuFSR2","dYToMuMuFSR2", nBins, xMin, xMax);
        TH1D *dYToMuMuNonFSR2 = new TH1D("dYToMuMuNonFSR2","dYToMuMuNonFSR2", nBins, xMin, xMax);
        TH1D *ttJets2 = new TH1D("ttJets2","ttJets2", nBins, xMin, xMax);
        TH1D *wJets2 = new TH1D("wJets2","wJets2", nBins, xMin, xMax);


	
	dataChain->Draw(Form("%s>>data",xVariable.c_str()),cut);
	dYToMuMuFSRChain->Draw(Form("%s>>dYToMuMuFSR",xVariable.c_str()),cut);
	dYToMuMuNonFSRChain->Draw(Form("%s>>dYToMuMuNonFSR",xVariable.c_str()),cut);
	ttJetsChain->Draw(Form("%s>>ttJets",xVariable.c_str()),cut);
	wJetsChain->Draw(Form("%s>>wJets",xVariable.c_str()),cut);	


	dataChain->Draw(Form("%s>>data2",xVariable.c_str()),cut2);
        dYToMuMuFSRChain->Draw(Form("%s>>dYToMuMuFSR2",xVariable.c_str()),cut2);
        dYToMuMuNonFSRChain->Draw(Form("%s>>dYToMuMuNonFSR2",xVariable.c_str()),cut2);
        ttJetsChain->Draw(Form("%s>>ttJets2",xVariable.c_str()),cut2);
        wJetsChain->Draw(Form("%s>>wJets2",xVariable.c_str()),cut2);

	// --- 2012 Lumi --- //
        double lumidata = 808.472 + 82.136 + 4429.0 + 495.003 + 134.242 + 6397.0 + 7274.0;
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
                /*weight = data->Integral() / dYToMuMuFSR->Integral();

                cout<<endl<<"weight = "<<weight<<endl;

                dYToMuMuFSR->Scale(weight);
                dYToMuMuNonFSR->Scale(weight);
                ttJets->Scale(weight);
                wJets->Scale(weight);
		*/


		ttJets2->Scale(lumiDY / lumiTtJets);
                wJets2->Scale(lumiDY / lumiwJets);

                dYToMuMuFSR2->Add(dYToMuMuNonFSR2);
                dYToMuMuFSR2->Add(ttJets2);
                dYToMuMuFSR2->Add(wJets2);

                dYToMuMuNonFSR2->Add(ttJets2);
                dYToMuMuNonFSR2->Add(wJets2);

                ttJets2->Add(wJets2);

                //weight = data->GetEntries() / dYToMuMuFSR->GetEntries();
                weight = dYToMuMuFSR->Integral() / dYToMuMuFSR2->Integral();

                cout<<endl<<"weight = "<<weight<<endl;

                dYToMuMuFSR2->Scale(weight);
                dYToMuMuNonFSR2->Scale(weight);
                ttJets2->Scale(weight);
                wJets2->Scale(weight);



	}
	
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
	dYToMuMuFSR2->Draw("SAMES");
	//dYToMuMuNonFSR->Draw("SAMES");
	//ttJets->Draw("SAMES");
	//wJets->Draw("SAMES");
	//data->Draw("E1SAMES");


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


	dYToMuMuFSR2->SetLineColor(kBlue);
	dYToMuMuFSR2->SetLineWidth(2);

	data->SetName("data");
	dYToMuMuFSR->SetName("dYToMuMuFSR");
	dYToMuMuNonFSR->SetName("dYToMuMuNonFSR");
	ttJets->SetName("ttJets");
	wJets->SetName("wJets");

	data2->SetName("data2");
        dYToMuMuFSR2->SetName("dYToMuMuFSR2");
        dYToMuMuNonFSR2->SetName("dYToMuMuNonFSR2");
        ttJets2->SetName("ttJets2");
        wJets2->SetName("wJets2");



	TLegend leg(xMinLeg,0.74,xMaxLeg,0.94,"","brNDC");
	leg.SetTextFont(42);
        leg.SetTextSize(0.028);
        leg.SetFillColor(kWhite);
        leg.SetLineColor(kWhite);
        leg.SetShadowColor(kWhite);
        //leg.AddEntry(data->GetName(),"data","lep");
        leg.AddEntry(dYToMuMuFSR->GetName(),"After #mu#mu#gamma cuts","f");
	leg.AddEntry(dYToMuMuFSR2->GetName(),"Before #mu#mu#gamma cuts","f");
        //leg.AddEntry(dYToMuMuNonFSR->GetName(),"Z#mu#mu + #gamma non FSR","f");
        //leg.AddEntry(ttJets->GetName(),"t#bar{t} + jets","f");
        //leg.AddEntry(wJets->GetName(),"W + jets","f");
	leg.Draw();

	TLatex latexLabel;
	latexLabel.SetTextFont(42);
	latexLabel.SetTextSize(0.028);
	latexLabel.SetNDC();
	latexLabel.DrawLatex(0.25, 0.96, "CMS Preliminary 2012               #sqrt{s} = 8 TeV               L = 19.6 fb^{-1}");
	
	std::ostringstream cutString2;
        cutString2 << setprecision (2) << fixed << dYToMuMuFSR->GetMean();
	string cutText2 = "mean_{After} = " + cutString2.str();	

	std::ostringstream cutString3;
        cutString3 << setprecision (2) << fixed << dYToMuMuFSR2->GetMean();
        string cutText3 = "mean_{Before} = " + cutString3.str();

	std::ostringstream cutString4;
        cutString4 << setprecision (2) << fixed << dYToMuMuFSR->GetRMS();
        string cutText4 = "RMS_{After} = " + cutString4.str();

	std::ostringstream cutString5;
        cutString5 << setprecision (2) << fixed << dYToMuMuFSR2->GetRMS();
        string cutText5 = "RMS_{Before} = " + cutString5.str();	


	
	latexLabel.DrawLatex(0.67, 0.70,Form("%s",cutText2.c_str()));
	latexLabel.DrawLatex(0.67, 0.65,Form("%s",cutText3.c_str()));
	latexLabel.DrawLatex(0.67, 0.60,Form("%s",cutText4.c_str()));		
	latexLabel.DrawLatex(0.67, 0.55,Form("%s",cutText5.c_str()));
	
		
	double Ymin = 0;
	double Ymax = max(dYToMuMuFSR->GetMaximum(),data->GetMaximum()) + max(dYToMuMuFSR->GetMaximum(),data->GetMaximum()) * 0.1;	
	
	if(log == 0) dYToMuMuFSR->GetYaxis()->SetRangeUser(Ymin,Ymax);
        if(log == 1) dYToMuMuFSR->SetMinimum(pow(10.0,-2));

	plotsRecording(directoryName, fileName, c1);


	c1->Clear();

	ofstream summaryFile(Form("%s/Summary.txt",directoryName_2.c_str()), ios::app);

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






 
