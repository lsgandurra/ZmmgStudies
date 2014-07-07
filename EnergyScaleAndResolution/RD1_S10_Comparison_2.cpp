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
	string cutVariable = "Photon_r9_2";
	string cutVariableValue = "0";
	string xVariable = "mcParticleTypeMatchedToGamma";
	string eta = "Endcaps";
	string r9 = "all";
	int log = 0;

	if( argc > 1 ) directoryName = argv[1];
	if( argc > 2 ) cutVariable = argv[2];
	if( argc > 3 ) cutVariableValue = argv[3];
	if( argc > 4 ) xVariable = argv[4];
	if( argc > 5 ) eta = argv[5];
	if( argc > 6 ) r9 = argv[6];
	if( argc > 7 ) 
        {
                std::stringstream ss ( argv[7] );
                ss >> log;
        }


	int nBins = 30;
	double xMin, xMax;
	double xMinLeg, xMaxLeg, yMinLeg, legTextSize;

	xMinLeg = 0.67;
	xMaxLeg = 0.91;	
	yMinLeg = 0.74;
	legTextSize = 0.028;
	string xVariableName, yVariableName;

	if(xVariable == "Photon_r9_2")
        {
                nBins = 50;
                xMin = 0;
                xMax = 1;
                xVariableName = "r9";
                yVariableName = "Events / 0.02";
		xMinLeg = 0.20;
		xMaxLeg = 0.44;
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
	
	if(xVariable == "Photon_SC_brem_2")
        {
                nBins = 40; 
                xMin = 0;
                xMax = 20;
                xVariableName = "Photon_SC_brem";
                yVariableName = "Events / 0.02";
                xMinLeg = 0.20;
                xMaxLeg = 0.44;
        }	


	if(xVariable == "Photon_Eta_2")
        {
                nBins = 52;
                xMin = -2.6;
                xMax = 2.6;
                xVariableName = "#eta^{#gamma}";
                yVariableName = "Events / 0.1";
                yMinLeg = 0.78;
                legTextSize = 0.026;
        }

        if(xVariable == "Photon_Phi_2")
        {
                nBins = 70;
                xMin = -3.5;
                xMax = 3.5;
                xVariableName = "#phi^{#gamma}";
                yVariableName = "Events / 0.1";
                yMinLeg = 0.78;
                legTextSize = 0.026;
        }	

	
	if(xVariable == "Photon_SC_Eta_2")
        {
                nBins = 52; 
                xMin = -2.6;
                xMax = 2.6; 
                xVariableName = "#eta{SC}_^{#gamma}";
                yVariableName = "Events / 0.1";
		yMinLeg = 0.78;
		legTextSize = 0.026;
        }

	if(xVariable == "Photon_SC_Phi_2")
        {
                nBins = 70; 
                xMin = -3.5;
                xMax = 3.5; 
                xVariableName = "#phi_{SC}^{#gamma}";
                yVariableName = "Events / 0.1";
                yMinLeg = 0.78;
                legTextSize = 0.026;
        }

	if(xVariable == "PhotonMC_Eta_2")
        {
                nBins = 52;
                xMin = -2.6;
                xMax = 2.6;
                xVariableName = "#eta{MC}_^{#gamma}";
                yVariableName = "Events / 0.1";
                yMinLeg = 0.78;
                legTextSize = 0.026;
        }

        if(xVariable == "PhotonMC_Phi_2")
        {
                nBins = 70;
                xMin = -3.5;
                xMax = 3.5;
                xVariableName = "#phi_{MC}^{#gamma}";
                yVariableName = "Events / 0.1";
                yMinLeg = 0.78;
                legTextSize = 0.026;
        }


	if(xVariable == "Photon_SC_E_2")
        {
                nBins = 50; 
                xMin = 0;
                xMax = 200;
                xVariableName = "E_{SC}^{#gamma}";
                yVariableName = "Events / 4";
                yMinLeg = 0.78;
                legTextSize = 0.026;
        }	

	if(xVariable == "Photon_SC_Eraw_2")
        {
                nBins = 50;
                xMin = 0;
                xMax = 200;
                xVariableName = "E_{SC RAW}^{#gamma}";
                yVariableName = "Events / 4";
                yMinLeg = 0.78;
                legTextSize = 0.026;
        }

	if(xVariable == "Photon_correctedE_2")
        {
                nBins = 50;
                xMin = 0;
                xMax = 200;
                xVariableName = "E_{corr}^{#gamma}";
                yVariableName = "Events / 4";
                yMinLeg = 0.78;
                legTextSize = 0.026;
        }

        if(xVariable == "Photon_E_2")
        {
                nBins = 50;
                xMin = 0;
                xMax = 200;
                xVariableName = "E_^{#gamma}";
                yVariableName = "Events / 4";
                yMinLeg = 0.78;
                legTextSize = 0.026;
        }

	if(xVariable == "PhotonMC_E_2")
        {
                nBins = 50;
                xMin = 0;
                xMax = 200;
                xVariableName = "E_{MC}^{#gamma}";
                yVariableName = "Events / 4";
                yMinLeg = 0.78;
                legTextSize = 0.026;
        }

	if(xVariable == "ParticleMC_E_2")
        {
                nBins = 50;
                xMin = 0;
                xMax = 200;
                xVariableName = "E_{MC}^{Gen part}";
                yVariableName = "Events / 4";
                yMinLeg = 0.78;
                legTextSize = 0.026;
        }

	if(xVariable == "PhotonMC_E_2-Photon_E_2")
        {
                nBins = 100; 
                xMin = -5.0;
                xMax = 5.0;
                xVariableName = "E_{MC}^{#gamma} - E_{RECO}^{#gamma}";
                yVariableName = "Events / 0.1";
                yMinLeg = 0.78;
                legTextSize = 0.026;
        }	

	if(xVariable == "ParticleMC_E_2-Photon_E_2")
        {
                nBins = 100; 
                xMin = -5.0;
                xMax = 5.0;
                xVariableName = "E_{MC}^{Gen part} - E_{RECO}^{#gamma}";
                yVariableName = "Events / 0.1";
                yMinLeg = 0.78;
                legTextSize = 0.026;
        }

	if(xVariable == "PhotonMC_Rconv_2")
        {
                nBins = 50;
                xMin = 0;
                xMax = 0.1;
                xVariableName = "#Rconv_{MC}^{#gamma}";
                yVariableName = "Events / 0.002";
                yMinLeg = 0.78;
                legTextSize = 0.026;
        }
	
	if(xVariable == "muon_Eta_2")
        {
                nBins = 52;
                xMin = -2.6;
                xMax = 2.6;
                xVariableName = "#eta^{#mu}";
                yVariableName = "Events / 0.1";
                yMinLeg = 0.78;
                legTextSize = 0.026;
        }

	if(xVariable == "muon_Pt_2")
        {
                nBins = 50;
                xMin = 0;
                xMax = 200;
                xVariableName = "P_{T}^{#mu}";
                yVariableName = "Events / 4";
                yMinLeg = 0.78;
                legTextSize = 0.026;
        }

	if(xVariable == "muon_E_2")
        {
                nBins = 50;
                xMin = 0;
                xMax = 200;
                xVariableName = "E_^{#mu}";
                yVariableName = "Events / 4";
                yMinLeg = 0.78;
                legTextSize = 0.026;
        }

        if(xVariable == "muon_Phi_2")
        {
                nBins = 70;
                xMin = -3.5;
                xMax = 3.5;
                xVariableName = "#phi^{#mu}";
                yVariableName = "Events / 0.1";
                yMinLeg = 0.78;
                legTextSize = 0.026;
        }

	if(xVariable == "muon_DeltaR_2")
        {
                nBins = 50; 
                xMin = 0;
                xMax = 5;
                xVariableName = "#Delta_{R}^{#mu_{Near},#gamma}";
                yVariableName = "Events / 0.1";
                yMinLeg = 0.78;
                legTextSize = 0.026;
        }

	if(xVariable == "mcParticleTypeMatchedToGamma")
        {
                nBins = 7000; 
                xMin = -3500;
                xMax = 3500;
                xVariableName = "Particle ID";
                yVariableName = "Events / 1";
                yMinLeg = 0.78;
                legTextSize = 0.026;
        }	

	if(xVariable == "diMuonPt_2")
        {
                nBins = 50; 
                xMin = 0;
                xMax = 200;
                xVariableName = "E_{T}^{#mu#mu}";
                yVariableName = "Events / 4";
                yMinLeg = 0.78;
                legTextSize = 0.026;
        }
	
	if(xVariable == "Photon_Et_2")
        {
                nBins = 50;
                xMin = 0;
                xMax = 200;
                xVariableName = "E_{T}^{#gamma}";
                yVariableName = "Events / 4";
                yMinLeg = 0.78;
                legTextSize = 0.026;
        }

	if(xVariable == "Photon_Et")
        {
                nBins = 50;
                xMin = 0;
                xMax = 200;
                xVariableName = "E_{T}^{#gamma}";
                yVariableName = "Events / 4";
                yMinLeg = 0.78;
                legTextSize = 0.026;
        }

	if(xVariable == "nSuperCluster")
        {
                nBins = 15;
                xMin = 0;
                xMax = 15;
                xVariableName = "nb SuperCluster";
                yVariableName = "Events / 1";
                yMinLeg = 0.78;
                legTextSize = 0.026;
        }


	
	gROOT->Reset();
        TGaxis::SetMaxDigits(3);
        setTDRStyle();
	//setEgammaStyle();

	string fileName = directoryName;
	directoryName += Form("/%s_%s/%s", cutVariable.c_str(), cutVariableValue.c_str(), xVariable.c_str());
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
	//TString cut = ""; //Form("(%s > %s && ",cutVariable.c_str(), cutVariableValue.c_str());
	//TString cut = "(Photon_SC_Eraw_2 > 70 && Photon_SC_Eraw_2 < 110";
	if(cutVariable == "NoCut") cutVariable = "Photon_r9";
	//TString cut = Form("(%s > %s && Photon_SC_E > 70 && Photon_SC_E < 110 && gammaGenMatched == 1",cutVariable.c_str(), cutVariableValue.c_str());
	TString cut = Form("(%s > %s && isJanLooseMMG == 1",cutVariable.c_str(), cutVariableValue.c_str());
	if(cutVariableValue == "runAB") cut = Form("(%s > 194000 && %s < 195000",cutVariable.c_str(), cutVariable.c_str());
	if(cutVariableValue == "runC") cut = Form("(%s > 200000 && %s < 201000",cutVariable.c_str(), cutVariable.c_str());
	if(cutVariableValue == "runD") cut = Form("(%s > 206000 && %s < 207000",cutVariable.c_str(), cutVariable.c_str());
	
	TString cut2 = "";
	if(cutVariableValue == "runAB" || cutVariableValue == "runC" || cutVariableValue == "runD") cut2 = Form("(%s > -1000000 ",cutVariable.c_str());
	else cut2 = cut;

	if(r9 == "low" && eta == "Barrel") cut += " && Photon_isEB == 1 && Photon_r9 < 0.94)";
        if(r9 == "high" && eta == "Barrel") cut += " && Photon_isEB == 1 && Photon_r9 > 0.94)";
        if(r9 == "low" && eta == "Endcaps") cut += " && Photon_isEE == 1 && Photon_r9 < 0.95)";
        if(r9 == "high" && eta == "Endcaps") cut += " && Photon_isEE == 1 && Photon_r9 > 0.95)";
        if(r9 == "all" && eta == "Barrel") cut += " && Photon_isEB == 1)";
        if(r9 == "all" && eta == "Endcaps") cut += " && Photon_isEE == 1)";
	if(r9 == "all" && eta == "all") cut += " && (Photon_isEE == 1 || Photon_isEB == 1))";

	if(r9 == "low" && eta == "Barrel") cut2 += " && Photon_isEB == 1 && Photon_r9 < 0.94)";
        if(r9 == "high" && eta == "Barrel") cut2 += " && Photon_isEB == 1 && Photon_r9 > 0.94)";
        if(r9 == "low" && eta == "Endcaps") cut2 += " && Photon_isEE == 1 && Photon_r9 < 0.95)";
        if(r9 == "high" && eta == "Endcaps") cut2 += " && Photon_isEE == 1 && Photon_r9 > 0.95)";
        if(r9 == "all" && eta == "Barrel") cut2 += " && Photon_isEB == 1)";
        if(r9 == "all" && eta == "Endcaps") cut2 += " && Photon_isEE == 1)";
        if(r9 == "all" && eta == "all") cut2 += " && (Photon_isEE == 1 || Photon_isEB == 1))";	

	cout << endl << "cut = "<<cut << endl;
	cout << endl << "cut2 = "<<cut2 << endl;

	
	// --- before muons selection --- //
	//TFile * dYToMuMuRD1File = new TFile("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniFriend2_totouples_DYToMuMu_Summer12_February2014_noskim_1_injRe0_v10_test2_partALL.root");
	//TFile * dYToMuMuPUS10File = new TFile("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniFriend2_DYToMuMu_Summer12_PU_S10_noskim_1_injRe0_v10_test2_partALL.root");

	// --- after dimuons selection --- //
	TFile * dYToMuMuRD1File = new TFile("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_totouples_DYToMuMu_Summer12_February2014_noskim_1_injRe0_v14_partALL.root");
        TFile * dYToMuMuPUS10File = new TFile("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_DYToMuMu_Summer12_PU_S10_noskim_1_injRe0_v14_partALL.root");
	
	TTree *dYToMuMuRD1Chain;
	TTree *dYToMuMuPUS10Chain;
	
	dYToMuMuRD1File->GetObject("miniTree",dYToMuMuRD1Chain);
	dYToMuMuPUS10File->GetObject("miniTree",dYToMuMuPUS10Chain);

	//dYToMuMuRD1Chain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniFriend2_totouples_DYToMuMu_Summer12_February2014_noskim_1_injRe0_v10_test2_partALL.root");
	//dYToMuMuPUS10Chain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniFriend2_DYToMuMu_Summer12_PU_S10_noskim_1_injRe0_v10_test2_partALL.root");

	TCanvas *c1 = new TCanvas("c1", "c1",0,0,600,600);
	if(log == 1) 
	{	
		c1->SetLogy();
		cout<<endl<<"c1->SetLogy() !"<<endl;
	}

	TH1D *dYToMuMuRD1 = new TH1D("dYToMuMuRD1","dYToMuMuRD1", nBins, xMin, xMax);
	TH1D *dYToMuMuPUS10 = new TH1D("dYToMuMuPUS10","dYToMuMuPUS10", nBins, xMin, xMax);
	
	
	dYToMuMuRD1Chain->Draw(Form("%s>>dYToMuMuRD1",xVariable.c_str()),cut);
	dYToMuMuPUS10Chain->Draw(Form("%s>>dYToMuMuPUS10",xVariable.c_str()),cut2);

	cout<<endl<<"coucou2";

	// --- Weight by integral --- //
	
	double weight = 1.0;
        weight = dYToMuMuRD1->Integral() / dYToMuMuPUS10->Integral();
	
	//weight = 1.0 / dYToMuMuPUS10->Integral();
	dYToMuMuPUS10->Scale(weight); //FIXME !!!!
		

	//weight = 1.0 / dYToMuMuRD1->Integral();
        //dYToMuMuRD1->Scale(weight);	

	c1->Clear();

	dYToMuMuPUS10->SetMarkerStyle(20);
	dYToMuMuPUS10->SetMarkerSize(0.5);

	dYToMuMuRD1->GetXaxis()->SetTitle(Form("%s",xVariableName.c_str()));
	dYToMuMuRD1->GetYaxis()->SetTitle(Form("%s",yVariableName.c_str()));
	/*
	dYToMuMuRD1->GetXaxis()->SetLabelFont(42);
	dYToMuMuRD1->GetXaxis()->SetTitleFont(42);
	dYToMuMuRD1->GetYaxis()->SetLabelFont(42);
	dYToMuMuRD1->GetYaxis()->SetTitleFont(42);
	dYToMuMuRD1->GetYaxis()->SetTitleOffset(1.65);
	*/
	dYToMuMuRD1->Draw("");
	dYToMuMuPUS10->Draw("E1SAMES");

	cout<<endl<<"coucou3";

	//dYToMuMuRD1->SetFillColor(2);
   	dYToMuMuRD1->SetFillColor(kGreen-7);
	//dYToMuMuRD1->SetFillStyle(3001);
	//dYToMuMuNonFSR->SetFillColor(3);
        //dYToMuMuNonFSR->SetFillStyle(3001);
	//ttJets->SetFillColor(4);
        //ttJets->SetFillStyle(3001);
	//wJets->SetFillColor(5);
        //wJets->SetFillStyle(3001);

	dYToMuMuPUS10->SetName("dYToMuMuPUS10");
	dYToMuMuRD1->SetName("dYToMuMuRD1");

	TLegend leg(xMinLeg,yMinLeg,xMaxLeg,0.94,"","brNDC");
	leg.SetTextFont(42);
        leg.SetTextSize(legTextSize);
        leg.SetFillColor(kWhite);
        leg.SetLineColor(kWhite);
        leg.SetShadowColor(kWhite);
        leg.AddEntry(dYToMuMuPUS10->GetName(),"DYToMuMu PU_S10","lep");
        leg.AddEntry(dYToMuMuRD1->GetName(),"DYToMuMu RD1","f");
	leg.Draw();

	TLatex latexLabel;
	latexLabel.SetTextFont(42);
	latexLabel.SetTextSize(0.028);
	latexLabel.SetNDC();
	latexLabel.DrawLatex(0.25, 0.96, "CMS Preliminary 2012               #sqrt{s} = 8 TeV               L = 19.6 fb^{-1}");
/*	
	std::ostringstream cutString2;
        cutString2 << setprecision (2) << fixed << nMCw;
	string cutText = "N_{MC} = " + cutString2.str();	
*/
/*
	if(xVariable == "Photon_Et")
        {
                latexLabel.DrawLatex(0.40, 0.83,Form("%s",cutText3.c_str()));
                latexLabel.DrawLatex(0.40, 0.78,Form("%s",cutText4.c_str()));
                latexLabel.DrawLatex(0.40, 0.73,Form("%s",cutText3b.c_str()));    
                latexLabel.DrawLatex(0.40, 0.68,Form("%s",cutText4b.c_str()));
        }
*/
	//if(r9sup == 2 && EndCaps == 2) latexLabel.DrawLatex(0.47, 0.63,Form("%s",cutText5.c_str()));
	
	cout<<endl<<"coucou4";
		
	double Ymin = 0;
	double Ymax = max(dYToMuMuRD1->GetMaximum(),dYToMuMuPUS10->GetMaximum()) + max(dYToMuMuRD1->GetMaximum(),dYToMuMuPUS10->GetMaximum()) * 0.1;	
	
	if(log == 0) dYToMuMuRD1->GetYaxis()->SetRangeUser(Ymin,Ymax);
	if(log == 0 && (xVariable == "MuonM_Eta" || xVariable == "MuonP_Eta" || xVariable == "Photon_SC_Eta_2")) dYToMuMuRD1->GetYaxis()->SetRangeUser(Ymin,Ymax + Ymax * 0.2); 
        if(log == 1) dYToMuMuRD1->SetMinimum(pow(10.0,-2));
	if(log == 1 && (xVariable == "MuonM_Eta" || xVariable == "MuonP_Eta" || xVariable == "Photon_SC_Eta_2" || xVariable == "Mmumugamma" || xVariable == "Mmumu")) dYToMuMuRD1->SetMaximum(30 * Ymax);

	//dYToMuMuRD1->GetXaxis()->SetLimits(-200,200);

	plotsRecording(directoryName, fileName, c1);


	c1->Clear();



	dYToMuMuPUS10->Delete();
	dYToMuMuPUS10 = 0;
	dYToMuMuRD1->Delete();
        dYToMuMuRD1 = 0;
	//dYToMuMuRD1Chain->Delete();
        //dYToMuMuRD1Chain = 0;
        //dYToMuMuPUS10Chain->Delete();
        //dYToMuMuPUS10Chain = 0;

	dYToMuMuRD1File->Delete();
	dYToMuMuPUS10File->Delete();
	
	cout<<endl<<"coucou5";
	delete c1;
	c1 = 0;	
	
	return 0;


}






 
