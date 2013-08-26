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

	string directoryName = "Data_MC_Mmumugamma_Comparison_ScaleAndSmearing/smearingFactors";
	string eta = "Barrel_1";
	string r9 = "all";
	string ptCut = "25"; // "30" "35"

	if( argc > 1 ) directoryName = argv[1];
	if( argc > 2 ) eta = argv[2];
        if( argc > 3 ) r9 = argv[3];	
	if( argc > 4 ) ptCut = argv[4];	

	


	int nBins = 50; //FIXME
	double xMin, xMax;
	xMin = 0.9; //FIXME
	xMax = 1.1; //FIXME
	
	string xVariableName, yVariableName;

	gROOT->Reset();
        TGaxis::SetMaxDigits(3);
        setTDRStyle();
	//setEgammaStyle();

	string directoryName_2 = directoryName;
	string fileName = directoryName;

	cout<<endl<<"fileName = "<<fileName<<endl;
	cout<<endl<<"directoryName = "<<directoryName<<endl;
	directoryName += Form("/%s_%sR9/",eta.c_str(), r9.c_str());
	
	TString cut = Form("(Photon_Et > %s && isJanLooseMMG == 1",ptCut.c_str());

	if(r9 == "low" && eta == "Barrel_1") cut += " && Photon_isEB == 1 && Photon_r9 < 0.94 && abs(Photon_SC_Eta) < 1";
        if(r9 == "low" && eta == "Barrel_2") cut += " && Photon_isEB == 1 && Photon_r9 < 0.94 && abs(Photon_SC_Eta) > 1";
	if(r9 == "high" && eta == "Barrel_1") cut += " && Photon_isEB == 1 && Photon_r9 > 0.94 && abs(Photon_SC_Eta) < 1";
	if(r9 == "high" && eta == "Barrel_2") cut += " && Photon_isEB == 1 && Photon_r9 > 0.94 && abs(Photon_SC_Eta) > 1";
        if(r9 == "low" && eta == "Endcaps_1") cut += " && Photon_isEE == 1 && Photon_r9 < 0.95 && abs(Photon_SC_Eta) < 2";
	if(r9 == "low" && eta == "Endcaps_2") cut += " && Photon_isEE == 1 && Photon_r9 < 0.95 && abs(Photon_SC_Eta) > 2";
        if(r9 == "high" && eta == "Endcaps_1") cut += " && Photon_isEE == 1 && Photon_r9 > 0.95 && abs(Photon_SC_Eta) < 2";
	if(r9 == "high" && eta == "Endcaps_2") cut += " && Photon_isEE == 1 && Photon_r9 > 0.95 && abs(Photon_SC_Eta) > 2";
        if(r9 == "all" && eta == "Barrel_1") cut += " && Photon_isEB == 1 && abs(Photon_SC_Eta) < 1";
	if(r9 == "all" && eta == "Barrel_2") cut += " && Photon_isEB == 1 && abs(Photon_SC_Eta) > 1";
        if(r9 == "all" && eta == "Endcaps_1") cut += " && Photon_isEE == 1 && abs(Photon_SC_Eta) < 2";
	if(r9 == "all" && eta == "Endcaps_2") cut += " && Photon_isEE == 1 && abs(Photon_SC_Eta) > 2";
	if(r9 == "all" && eta == "all") cut += " && (Photon_isEE == 1 || Photon_isEB == 1)";

	cut += ")*weight_pileUp"; 

	cout<<endl<<"cut = "<<cut<<endl;

        TChain * dYToMuMuFSRChainNew = new TChain("miniTree");
        TChain * dYToMuMuNonFSRChainNew = new TChain("miniTree");
        TChain * ttJetsChainNew = new TChain("miniTree");
        TChain * wJetsChainNew = new TChain("miniTree");

        dYToMuMuFSRChainNew->Add("/sps/cms/sgandurr/CMSSW_4_2_8_patch7/src/Zmumugamma_miniTrees_rereco_2011_lastTag/miniTree_DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Fall11-PU_S6_START42_V14B-v1_September12_NewSelection_1_scaleAndsmearing_v4_partALL.root");
        dYToMuMuNonFSRChainNew->Add("/sps/cms/sgandurr/CMSSW_4_2_8_patch7/src/Zmumugamma_miniTrees_rereco_2011_lastTag/miniTree_DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Fall11-PU_S6_START42_V14B-v1_September12_NewSelection_2_scaleAndsmearing_v4_partALL.root");
        ttJetsChainNew->Add("/sps/cms/sgandurr/CMSSW_4_2_8_patch7/src/Zmumugamma_miniTrees_rereco_2011_lastTag/miniTree_TTJets_TuneZ2_7TeV-madgraph-tauola_NewSelection_3_scaleAndsmearing_v4_partALL.root");
        wJetsChainNew->Add("/sps/cms/sgandurr/CMSSW_4_2_8_patch7/src/Zmumugamma_miniTrees_rereco_2011_lastTag/miniTree_WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_NewSelection_3_scaleAndsmearing_v4_partALL.root");
	

	TCanvas *c1 = new TCanvas("c1", "c1",0,0,600,600);

        TH1D *dYToMuMuFSRNew = new TH1D("dYToMuMuFSRNew","dYToMuMuFSRNew", nBins, xMin, xMax);
        TH1D *dYToMuMuNonFSRNew = new TH1D("dYToMuMuNonFSRNew","dYToMuMuNonFSRNew", nBins, xMin, xMax);
        TH1D *ttJetsNew = new TH1D("ttJetsNew","ttJetsNew", nBins, xMin, xMax);
        TH1D *wJetsNew = new TH1D("wJetsNew","wJetsNew", nBins, xMin, xMax);

        dYToMuMuFSRChainNew->Draw("shervinSmearing>>dYToMuMuFSRNew",cut);
        dYToMuMuNonFSRChainNew->Draw("shervinSmearing>>dYToMuMuNonFSRNew",cut);
        ttJetsChainNew->Draw("shervinSmearing>>ttJetsNew",cut);
        wJetsChainNew->Draw("shervinSmearing>>wJetsNew",cut);


	// --- 2012 Lumi --- //
        //double lumidata = 808.472 + 82.136 + 4429.0 + 495.003 + 134.242 + 6397.0 + 7274.0;
        //double lumiDY = 48819386.0 / 1914.894;
        //double lumiTtJets = 6736135.0 / 234.0;
        //double lumiwJets = 57709905.0 / 37509.25;	

	// --- 2011 Lumi --- //
	double lumidata = (0.706370 + 0.385819 + 2.741 + 1.099) * 1000;
	double lumiDY = 29743564.0 / 1665.835;
	double lumiTtJets = 3701947.0 / 165.0;
	double lumiwJets = 81345381.0 / 31314.0;	

	// --- 2011 Lumi old --- //
	//double lumiDY = 29743564.0 / 1626.0;
        //double lumiTtJets = 3701947.0 / 94.76;
        //double lumiwJets = 81345381.0 / 27770.0;

	ttJetsNew->Scale(lumiDY / lumiTtJets);
	wJetsNew->Scale(lumiDY / lumiwJets);
		
	dYToMuMuFSRNew->Add(dYToMuMuNonFSRNew);
	dYToMuMuFSRNew->Add(ttJetsNew);
	dYToMuMuFSRNew->Add(wJetsNew);

	dYToMuMuNonFSRNew->Add(ttJetsNew);
	dYToMuMuNonFSRNew->Add(wJetsNew);

	ttJetsNew->Add(wJetsNew);

	c1->Clear();


	dYToMuMuFSRNew->GetYaxis()->SetTitle("Events / 0.004"); //FIXME
       	dYToMuMuFSRNew->GetXaxis()->SetTitle("smearing factors");	

	/*
	dYToMuMuFSR->GetXaxis()->SetLabelFont(42);
	dYToMuMuFSR->GetXaxis()->SetTitleFont(42);
	dYToMuMuFSR->GetYaxis()->SetLabelFont(42);
	dYToMuMuFSR->GetYaxis()->SetTitleFont(42);
	dYToMuMuFSR->GetYaxis()->SetTitleOffset(1.65);
	*/

	dYToMuMuFSRNew->SetFillColor(kGreen-7);
        dYToMuMuNonFSRNew->SetFillColor(kAzure-5);
        ttJetsNew->SetFillColor(kBlue-1);
        wJetsNew->SetFillColor(kCyan+2);

	dYToMuMuFSRNew->Draw("");


	// --- Fit --- //

        TF1 * f1 = new TF1("f1","gaus",xMin,xMax);
	//TF1 * f1 = new TF1("f1","gaus",0.96,1.04);
	//f1->FixParameter(1,1.0); //FIXME
        //f1->FixParameter(2,1.0); //FIXME
        dYToMuMuFSRNew->Fit(f1);
        f1->SetLineColor(kBlue);
        f1->SetLineWidth(2);
        f1->Draw("SAMES");

        c1->Clear();
        dYToMuMuFSRNew->Draw("E");
        f1->Draw("SAMES");

        dYToMuMuFSRNew->SetMarkerStyle(20);
        dYToMuMuFSRNew->SetMarkerSize(0.5);
	
	TLatex latexLabel;
	latexLabel.SetTextFont(42);
	latexLabel.SetTextSize(0.028);
	latexLabel.SetNDC();
	//latexLabel.DrawLatex(0.25, 0.96, "CMS Preliminary 2011               #sqrt{s} = 7 TeV               L = 4.93 fb^{-1}");
	
	//double Ymin = 0;
	//double Ymax = max(dYToMuMuFSR->GetMaximum(),data->GetMaximum()) + max(dYToMuMuFSR->GetMaximum(),data->GetMaximum()) * 0.1;	
	//cout << "Ymax = "<<Ymax<<endl;

	//dYToMuMuFSRNew->GetYaxis()->SetRangeUser(Ymin,Ymax);

	plotsRecording(directoryName, fileName, c1);

	c1->Clear();

	f1->Delete();
	f1 = 0;
        dYToMuMuFSRNew->Delete();
        dYToMuMuFSRNew = 0;
        dYToMuMuNonFSRNew->Delete();
        dYToMuMuNonFSRNew = 0;
        ttJetsNew->Delete();
        ttJetsNew = 0;
        wJetsNew->Delete();
        wJetsNew = 0;
        dYToMuMuFSRChainNew->Delete();
        dYToMuMuFSRChainNew = 0;
        dYToMuMuNonFSRChainNew->Delete();
        dYToMuMuNonFSRChainNew = 0;
        ttJetsChainNew->Delete();
        ttJetsChainNew = 0;
        wJetsChainNew->Delete();
        wJetsChainNew = 0;


	
	delete c1;
	c1 = 0;	
	
	return 0;

}






 
