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
                cerr << "arguments should be passed : directoryName, eta, r9, log" <<endl; 
                return 1;

        }    

	string directoryName = "data_MC_Mmumugamma_Comparison";
	string eta = "Barrel";
	string r9 = "all";
	int log = 1;		

	if( argc > 1 ) directoryName = argv[1];
	if( argc > 2 ) eta = argv[2];
        if( argc > 3 ) r9 = argv[3];	
	if( argc > 4 ) 
        {
                std::stringstream ss ( argv[4] );
                ss >> log;
        }

	
	gROOT->Reset();
        TGaxis::SetMaxDigits(3);
        setTDRStyle();
	//setEgammaStyle();

	string directoryName_2 = directoryName;

	if(eta == "Barrel") directoryName += "/Barrel_";
        if(eta == "Endcaps") directoryName += "/Endcaps_";
        if(eta == "all") directoryName += "/BarrelAndEncaps_";

        if(r9 == "low") directoryName += "lowR9/";
        if(r9 == "high") directoryName += "highR9/";
        if(r9 == "all") directoryName += "AllR9/";

	string fileName = "data_MC_Mmumugamma_Comparison";
	if(log == 1) fileName += "_log";

	TString cut = "(Photon_Et > 25 && isJanLooseMMG == 1";

	if(r9 == "low" && eta == "Barrel") cut += " && Photon_isEB == 1 && Photon_r9 < 0.94";
        if(r9 == "high" && eta == "Barrel") cut += " && Photon_isEB == 1 && Photon_r9 > 0.94";
        if(r9 == "low" && eta == "Endcaps") cut += " && Photon_isEE == 1 && Photon_r9 < 0.95";
        if(r9 == "high" && eta == "Endcaps") cut += " && Photon_isEE == 1 && Photon_r9 > 0.95";
        if(r9 == "all" && eta == "Barrel") cut += " && Photon_isEB == 1";
        if(r9 == "all" && eta == "Endcaps") cut += " && Photon_isEE == 1";
	if(r9 == "all" && eta == "all") cut += " && (Photon_isEE == 1 || Photon_isEB == 1)";

	cut += ")*weight_pileUp*weight_Xsection"; // weight by luminosity
	//cut += ")*weight_pileUp"; // weight by entries

	cout<<endl<<"cut = "<<cut<<endl;

        TChain * dataChain = new TChain("miniTree");
	TChain * dYToMuMuFSRChain = new TChain("miniTree");
	TChain * dYToMuMuNonFSRChain = new TChain("miniTree");
	TChain * ttJetsChain = new TChain("miniTree");
	TChain * wJetsChain = new TChain("miniTree");

	dataChain->Add("/sps/cms/sgandurr/CMSSW_4_2_8_patch7/src/Zmumugamma_miniTrees_rereco_2011_lastTag/miniTree_2011A_03Oct2011V1ReReco_toto_v2_NewSelection_0_v6_partALL.root");
	dataChain->Add("/sps/cms/sgandurr/CMSSW_4_2_8_patch7/src/Zmumugamma_miniTrees_rereco_2011_lastTag/miniTree_2011A_05Jul2011ReReco_toto_v2_NewSelection_0_v6_partALL.root");
	dataChain->Add("/sps/cms/sgandurr/CMSSW_4_2_8_patch7/src/Zmumugamma_miniTrees_rereco_2011_lastTag/miniTree_2011A_PromptSkimV5ReReco_toto_v2_NewSelection_0_v6_partALL.root");
	dataChain->Add("/sps/cms/sgandurr/CMSSW_4_2_8_patch7/src/Zmumugamma_miniTrees_rereco_2011_lastTag/miniTree_2011B_PromptSkimV1ReReco_toto_v2_NewSelection_0_v6_partALL.root");
	dYToMuMuFSRChain->Add("/sps/cms/sgandurr/CMSSW_4_2_8_patch7/src/Zmumugamma_miniTrees_rereco_2011_lastTag/miniTree_DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Fall11-PU_S6_START42_V14B-v1_September12_NewSelection_1_v6_partALL.root");
	dYToMuMuNonFSRChain->Add("/sps/cms/sgandurr/CMSSW_4_2_8_patch7/src/Zmumugamma_miniTrees_rereco_2011_lastTag/miniTree_DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Fall11-PU_S6_START42_V14B-v1_September12_NewSelection_2_v6_partALL.root");
	ttJetsChain->Add("/sps/cms/sgandurr/CMSSW_4_2_8_patch7/src/Zmumugamma_miniTrees_rereco_2011_lastTag/miniTree_TTJets_TuneZ2_7TeV-madgraph-tauola_NewSelection_3_v6_partALL.root");
	wJetsChain->Add("/sps/cms/sgandurr/CMSSW_4_2_8_patch7/src/Zmumugamma_miniTrees_rereco_2011_lastTag/miniTree_WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_NewSelection_3_v6_partALL.root");

	

	TCanvas *c1 = new TCanvas("c1", "c1",0,0,600,600);
	if(log == 1) c1->SetLogy();

	TH1D *data = new TH1D("data","data", 60, 60, 120);	
	TH1D *dYToMuMuFSR = new TH1D("dYToMuMuFSR","dYToMuMuFSR", 60, 60, 120);
	TH1D *dYToMuMuNonFSR = new TH1D("dYToMuMuNonFSR","dYToMuMuNonFSR", 60, 60, 120);
	TH1D *ttJets = new TH1D("ttJets","ttJets", 60, 60, 120);
	TH1D *wJets = new TH1D("wJets","wJets", 60, 60, 120);				
		
	dataChain->Draw("Mmumugamma>>data",cut);
	dYToMuMuFSRChain->Draw("Mmumugamma>>dYToMuMuFSR",cut);
	dYToMuMuNonFSRChain->Draw("Mmumugamma>>dYToMuMuNonFSR",cut);
	ttJetsChain->Draw("Mmumugamma>>ttJets",cut);
	wJetsChain->Draw("Mmumugamma>>wJets",cut);	


	// --- 2011 Lumi --- //
	double lumidata = (0.706370 + 0.385819 + 2.741 + 1.099) * 1000;
	double lumiDY = 29743564.0 / 1665.835;
	double lumiTtJets = 3701947.0 / 165.0;
	double lumiwJets = 81345381.0 / 31314.0;	

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
	
	/*

	// --- Weight by entries --- //

	ttJets->Scale(lumiDY / lumiTtJets);
	wJets->Scale(lumiDY / lumiwJets);
	
	dYToMuMuFSR->Add(dYToMuMuNonFSR);
	dYToMuMuFSR->Add(ttJets);
	dYToMuMuFSR->Add(wJets);

	dYToMuMuNonFSR->Add(ttJets);
	dYToMuMuNonFSR->Add(wJets);

	ttJets->Add(wJets);

	double weight = data->GetEntries() / dYToMuMuFSR->GetEntries();	

	cout<<endl<<"weight = "<<weight<<endl;

	dYToMuMuFSR->Scale(weight);
	dYToMuMuNonFSR->Scale(weight);
	ttJets->Scale(weight);
	wJets->Scale(weight);

	*/
	
	// --- Weight by luminosity --- //

	dYToMuMuFSR->Add(dYToMuMuNonFSR);
        dYToMuMuFSR->Add(ttJets);
        dYToMuMuFSR->Add(wJets);

        dYToMuMuNonFSR->Add(ttJets);
        dYToMuMuNonFSR->Add(wJets);

        ttJets->Add(wJets);

	

	c1->Clear();

	double meanEtData = data->GetMean();
        double meanEtMC = dYToMuMuFSR->GetMean();
	
	data->SetMarkerStyle(20);
	data->SetMarkerSize(0.5);

	dYToMuMuFSR->GetXaxis()->SetTitle("M_{#mu#mu#gamma} [GeV]");
	dYToMuMuFSR->GetYaxis()->SetTitle("Events / 1 GeV");
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

	TLegend leg(0.17,0.72,0.41,0.92,"","brNDC");
	leg.SetTextFont(42);
        leg.SetTextSize(0.028);
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
	latexLabel.DrawLatex(0.25, 0.96, "CMS Preliminary 2011               #sqrt{s} = 7 TeV               L = 4.93 fb^{-1}");
	
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
	double Ymax = max(dYToMuMuFSR->GetMaximum(),data->GetMaximum()) + max(dYToMuMuFSR->GetMaximum(),data->GetMaximum()) * 0.1;	
	
	if(log == 0) dYToMuMuFSR->GetYaxis()->SetRangeUser(Ymin,Ymax);
        if(log == 1) dYToMuMuFSR->SetMinimum(pow(10.0,-2));

	plotsRecording(directoryName, fileName, c1);


	c1->Clear();

	ofstream summaryFile(Form("%s/Summary.txt",directoryName_2.c_str()), ios::app);

	if(log == 0)
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






 
