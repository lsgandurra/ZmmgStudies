// -------------------------------------------------------
// program implemented by Louis Sgandurra (September 2012)
// -------------------------------------------------------

#include "TCanvas.h"
#include "TAxis.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TH2F.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TBranch.h"
#include "TChain.h"
#include "TFile.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TString.h"
#include "TBits.h"
#include "TSystem.h"
#include "TPaveStats.h"
#include "TPaveText.h"
#include "TGraphAsymmErrors.h"
#include "TLatex.h"
#include "TGraphErrors.h"
#include "TPad.h"
#include "TGraph.h"
#include "TText.h"
#include "TLine.h"
#include "TGaxis.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstring>
#include <iomanip>
#include "TROOT.h"
#include "TRint.h"
#include "TMultiGraph.h"
#include "setTDRStyle.C"
#include <cmath>
#include <sys/types.h>
#include <sys/stat.h>
#include <algorithm>
#include <vector>
#include <set>
#include <list>
#include <time.h>
#include <stdio.h>
#pragma optimize 0

using namespace std;

//int main()
int StatCategories()
{
	double  lumiData = 706.370 + 385.819 + 2741 + 1099;
/*
        double  xSectionDY = 1914.894;
        double  xSectionTtJets = 234;
        double  xSectionWJets = 37509.25;
*/
	double  xSectionDY = 1665.835;
	double  xSectionTtJets = 165.0;
	double  xSectionWJets = 31314.0;

	double InitialNumberDYToMuMu = 29743564.0;
	double InitialNumberTTJets = 3701947.0;
	double InitialNumberWJetsToLNu = 81345381.0;

	double lumiDY = InitialNumberDYToMuMu / xSectionDY;
	double lumiTtJets = InitialNumberTTJets / xSectionTtJets;
	double lumiWJets = InitialNumberWJetsToLNu / xSectionWJets;	

	double purity = 0;

        TChain * dataChain;
	TChain * dataChain1;
	TChain * dataChain2;
	TChain * dataChain3;
	TChain * dataChain4;
	//TChain * dataChain5;
	//TChain * dataChain6;
	//TChain * dataChain7;	
	TChain * dyToMuMuFSRChain;
	TChain * dyToMuMuNonFSRChain;
	TChain * ttJetsChain;
	TChain * wJetsToLNuChain;	
		
	TString cuts = "isJanLooseMMG == 1";

	TChain * reducedChain;
	TChain * reducedChain1;
	TChain * reducedChain2;
	TChain * reducedChain3;
	TChain * reducedChain4;
	//TChain * reducedChain5;
	//TChain * reducedChain6;
	//TChain * reducedChain7;
	TChain * reducedChainDyToMuMuFSR;
	TChain * reducedChainDyToMuMuNonFSR;
	TChain * reducedChainTTJets;
	TChain * reducedChainWJetsToLNu;


	cout<<endl<<"**** Data (4.932 fb-1) ****"<<endl;
	for(int i = 0; i < 3; i++){

	if(i == 0) cout<<endl<<"--- No Pt cuts ---";
	if(i == 1) cout<<endl<<"--- Pt > 20 GeV ---";
	if(i == 2) cout<<endl<<"--- Pt > 25 GeV ---";
	for(int EndCaps = 0; EndCaps <= 2; EndCaps++)
	{
		for(int r9sup = 0; r9sup <= 2; r9sup++)
		{
			dataChain = new TChain("miniTree");
		        dataChain1 = new TChain("miniTree");
		        dataChain2 = new TChain("miniTree");
		        dataChain3 = new TChain("miniTree");
		        dataChain4 = new TChain("miniTree");
		        //dataChain5 = new TChain("miniTree");
		        //dataChain6 = new TChain("miniTree");
		        //dataChain7 = new TChain("miniTree");   
		        //dyToMuMuFSRChain = new TChain("miniTree");
		        //dyToMuMuNonFSRChain = new TChain("miniTree");
		        //ttJetsChain = new TChain("miniTree");
		        //wJetsToLNuChain = new TChain("miniTree");    
		     
		        dataChain->Add("../Zmumugamma_miniTrees_rereco_2011_lastTag/miniTree_2011A_03Oct2011V1ReReco_toto_v2_NewSelection_0_v4_partALL.root");
		        dataChain->Add("../Zmumugamma_miniTrees_rereco_2011_lastTag/miniTree_2011A_05Jul2011ReReco_toto_v2_NewSelection_0_v4_partALL.root");
		        dataChain->Add("../Zmumugamma_miniTrees_rereco_2011_lastTag/miniTree_2011A_PromptSkimV5ReReco_toto_v2_NewSelection_0_v4_partALL.root");
		        dataChain->Add("../Zmumugamma_miniTrees_rereco_2011_lastTag/miniTree_2011B_PromptSkimV1ReReco_toto_v2_NewSelection_0_v4_partALL.root");
		        //dataChain->Add("../Zmumugamma_miniTrees_rereco_2011_lastTag/miniTree_Run2012C-EcalRecover_11Dec2012-v1_NewSelection_0_v1_partALL.root");
		        //dataChain->Add("../Zmumugamma_miniTrees_rereco_2011_lastTag/miniTree_Run2012C_PromptReco_v2_v2_NewSelection_0_v1_partALL.root");
		        //dataChain->Add("../Zmumugamma_miniTrees_rereco_2011_lastTag/miniTree_Run2012D_PromptReco_v1_v2_NewSelection_0_v1_partALL.root");
		        dataChain1->Add("../Zmumugamma_miniTrees_rereco_2011_lastTag/miniTree_2011A_03Oct2011V1ReReco_toto_v2_NewSelection_0_v4_partALL.root");
		        dataChain2->Add("../Zmumugamma_miniTrees_rereco_2011_lastTag/miniTree_2011A_05Jul2011ReReco_toto_v2_NewSelection_0_v4_partALL.root");
		        dataChain3->Add("../Zmumugamma_miniTrees_rereco_2011_lastTag/miniTree_2011A_PromptSkimV5ReReco_toto_v2_NewSelection_0_v4_partALL.root");
		        dataChain4->Add("../Zmumugamma_miniTrees_rereco_2011_lastTag/miniTree_2011B_PromptSkimV1ReReco_toto_v2_NewSelection_0_v4_partALL.root");
		        //dataChain5->Add("../Zmumugamma_miniTrees_rereco_2011_lastTag/miniTree_Run2012C-EcalRecover_11Dec2012-v1_NewSelection_0_v1_partALL.root");
		        //dataChain6->Add("../Zmumugamma_miniTrees_rereco_2011_lastTag/miniTree_Run2012C_PromptReco_v2_v2_NewSelection_0_v1_partALL.root");
		        //dataChain7->Add("../Zmumugamma_miniTrees_rereco_2011_lastTag/miniTree_Run2012D_PromptReco_v1_v2_NewSelection_0_v1_partALL.root");
		

			if((EndCaps == 2 && r9sup == 0) || (EndCaps == 2 && r9sup == 1)) continue;

			if(i == 0) cuts = "weight_Xsection*(isJanLooseMMG == 1";
			if(i == 1) cuts = "weight_Xsection*(isJanLooseMMG == 1 && Photon_Et > 20";
			if(i == 2) cuts = "weight_Xsection*(isJanLooseMMG == 1 && Photon_Et > 25";

	       		if(EndCaps == 0 && r9sup == 1) cuts += " && Photon_isEB == 1 && Photon_r9 > 0.94)";
	       		if(EndCaps == 0 && r9sup == 0) cuts += " && Photon_isEB == 1 && Photon_r9 < 0.94)";
	       		if(EndCaps == 1 && r9sup == 1) cuts += " && Photon_isEE == 1 && Photon_r9 > 0.95)";
	       		if(EndCaps == 1 && r9sup == 0) cuts += " && Photon_isEE == 1 && Photon_r9 < 0.95)";
	       		if(EndCaps == 0 && r9sup == 2) cuts += " && Photon_isEB == 1)";
	       		if(EndCaps == 1 && r9sup == 2) cuts += " && Photon_isEE == 1)";
			if(EndCaps == 2 && r9sup == 2) cuts += ")";
			//if(EndCaps == 2 && r9sup == 2) cuts += " && (Photon_isEE == 1 || Photon_isEB == 1))";
			//if(EndCaps == 2 && r9sup == 2) cuts += " && (abs(Photon_SC_Eta) < 1.4442 || abs(Photon_SC_Eta) > 1.566))";		
	
			//cout<<endl<<"cuts = "<<cuts<<endl;	
			reducedChain = (TChain *) dataChain->CopyTree(cuts);
			reducedChain1 = (TChain *) dataChain1->CopyTree(cuts);
			reducedChain2 = (TChain *) dataChain2->CopyTree(cuts);
			reducedChain3 = (TChain *) dataChain3->CopyTree(cuts);
			reducedChain4 = (TChain *) dataChain4->CopyTree(cuts);
			//reducedChain5 = (TChain *) dataChain5->CopyTree(cuts);
			//reducedChain6 = (TChain *) dataChain6->CopyTree(cuts);
			//reducedChain7 = (TChain *) dataChain7->CopyTree(cuts);
	
			if(EndCaps == 0 && r9sup == 2) cout<<endl<<">>EB All r9 : "<<reducedChain->GetEntries();
			if(EndCaps == 0 && r9sup == 1) cout<<endl<<">>EB high r9 : "<<reducedChain->GetEntries();
			if(EndCaps == 0 && r9sup == 0) cout<<endl<<">>EB low r9 : "<<reducedChain->GetEntries();
				
			if(EndCaps == 1 && r9sup == 2) cout<<endl<<">>EE All r9 : "<<reducedChain->GetEntries();
	                if(EndCaps == 1 && r9sup == 1) cout<<endl<<">>EE high r9 : "<<reducedChain->GetEntries();
	                if(EndCaps == 1 && r9sup == 0) cout<<endl<<">>EE low r9 : "<<reducedChain->GetEntries();
			if(EndCaps == 2 && r9sup == 2) cout<<endl<<">>EE + EB All r9 (with the transition region) : "<<reducedChain->GetEntries();
			cout<<endl<<"2011A_03Oct2011V1ReReco : "<<reducedChain1->GetEntries();
			cout<<endl<<"2011A_05Jul2011ReReco : "<<reducedChain2->GetEntries();
			cout<<endl<<"2011A_PromptSkimV5ReReco : "<<reducedChain3->GetEntries();
			cout<<endl<<"2011B_PromptSkimV1ReReco : "<<reducedChain4->GetEntries();
			//cout<<endl<<"Run2012C-EcalRecover_11Dec2012-v1 : "<<reducedChain5->GetEntries();
			//cout<<endl<<"Run2012C_PromptReco_v2 : "<<reducedChain6->GetEntries();
			//cout<<endl<<"Run2012D_PromptReco_v1 : "<<reducedChain7->GetEntries();
	
			dataChain->Delete();
		        dataChain = 0;
		        dataChain1->Delete();
		        dataChain1 = 0;
		        dataChain2->Delete();
		        dataChain2 = 0;
		        dataChain3->Delete();
		        dataChain3 = 0;
		        dataChain4->Delete();
		        dataChain4 = 0;
		        //dataChain5->Delete();
		        //dataChain5 = 0;
		        //dataChain6->Delete();
		        //dataChain6 = 0;
		        //dataChain7->Delete();
		        //dataChain7 = 0;
		
			reducedChain->Delete();	
			reducedChain = 0;
			reducedChain1->Delete();
			reducedChain1 = 0;
			reducedChain2->Delete();
			reducedChain2 = 0;
			reducedChain3->Delete();
			reducedChain3 = 0;
			reducedChain4->Delete();
			reducedChain4 = 0;
			//reducedChain5->Delete();
			//reducedChain5 = 0;
			//reducedChain6->Delete();
			//reducedChain6 = 0;
			//reducedChain7->Delete();
			//reducedChain7 = 0;
			cout<<endl;
		}
		cout<<endl;
	}
	cout<<endl;
	}

	cout<<endl;

	cout<<endl<<"**** MC : DYtoMuMu (FSR + nonFSR) + ttjets + Wjets ****"<<endl;

	for(int i = 0; i < 3; i++){

        if(i == 0) cout<<endl<<"--- No Pt cuts ---";
        if(i == 1) cout<<endl<<"--- Pt > 20 GeV ---";
        if(i == 2) cout<<endl<<"--- Pt > 25 GeV ---";
        for(int EndCaps = 0; EndCaps <= 2; EndCaps++)
        {
                for(int r9sup = 0; r9sup <= 2; r9sup++)
                {
			dyToMuMuFSRChain = new TChain("miniTree");
		        dyToMuMuNonFSRChain = new TChain("miniTree");
		        ttJetsChain = new TChain("miniTree");
		        wJetsToLNuChain = new TChain("miniTree");
		
		        dyToMuMuFSRChain->Add("../Zmumugamma_miniTrees_rereco_2011_lastTag/miniTree_DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Fall11-PU_S6_START42_V14B-v1_September12_NewSelection_1_v4_partALL.root");
		        dyToMuMuNonFSRChain->Add("../Zmumugamma_miniTrees_rereco_2011_lastTag/miniTree_DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Fall11-PU_S6_START42_V14B-v1_September12_NewSelection_2_v4_partALL.root");
		        ttJetsChain->Add("../Zmumugamma_miniTrees_rereco_2011_lastTag/miniTree_TTJets_TuneZ2_7TeV-madgraph-tauola_NewSelection_3_v4_partALL.root");
		        wJetsToLNuChain->Add("../Zmumugamma_miniTrees_rereco_2011_lastTag/miniTree_WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_NewSelection_3_v4_partALL.root");
		
                        if((EndCaps == 2 && r9sup == 0) || (EndCaps == 2 && r9sup == 1)) continue;

                        if(i == 0) cuts = "weight_Xsection*(isJanLooseMMG == 1";
                        if(i == 1) cuts = "weight_Xsection*(isJanLooseMMG == 1 && Photon_Et > 20";
                        if(i == 2) cuts = "weight_Xsection*(isJanLooseMMG == 1 && Photon_Et > 25";

                        if(EndCaps == 0 && r9sup == 1) cuts += " && Photon_isEB == 1 && Photon_r9 > 0.94)";
                        if(EndCaps == 0 && r9sup == 0) cuts += " && Photon_isEB == 1 && Photon_r9 < 0.94)";
                        if(EndCaps == 1 && r9sup == 1) cuts += " && Photon_isEE == 1 && Photon_r9 > 0.95)";
                        if(EndCaps == 1 && r9sup == 0) cuts += " && Photon_isEE == 1 && Photon_r9 < 0.95)";
                        if(EndCaps == 0 && r9sup == 2) cuts += " && Photon_isEB == 1)";
                        if(EndCaps == 1 && r9sup == 2) cuts += " && Photon_isEE == 1)";
                        if(EndCaps == 2 && r9sup == 2) cuts += ")";
                        //if(EndCaps == 2 && r9sup == 2) cuts += " && (Photon_isEE == 1 || Photon_isEB == 1))";
                        //if(EndCaps == 2 && r9sup == 2) cuts += " && (abs(Photon_SC_Eta) < 1.4442 || abs(Photon_SC_Eta) > 1.566))";            

                        //cout<<endl<<"cuts = "<<cuts<<endl;   
			reducedChainDyToMuMuFSR = (TChain *) dyToMuMuFSRChain->CopyTree(cuts);
			reducedChainDyToMuMuNonFSR = (TChain *) dyToMuMuNonFSRChain->CopyTree(cuts);
        		reducedChainTTJets = (TChain *) ttJetsChain->CopyTree(cuts);
			reducedChainWJetsToLNu = (TChain *) wJetsToLNuChain->CopyTree(cuts);		
			

                        if(EndCaps == 0 && r9sup == 2) cout<<endl<<">>EB All r9 : ";
                        if(EndCaps == 0 && r9sup == 1) cout<<endl<<">>EB high r9 : ";
                        if(EndCaps == 0 && r9sup == 0) cout<<endl<<">>EB low r9 : ";

                        if(EndCaps == 1 && r9sup == 2) cout<<endl<<">>EE All r9 : ";
                        if(EndCaps == 1 && r9sup == 1) cout<<endl<<">>EE high r9 : ";
                        if(EndCaps == 1 && r9sup == 0) cout<<endl<<">>EE low r9 : ";
                        if(EndCaps == 2 && r9sup == 2) cout<<endl<<">>EE + EB All r9 (with the transition region) : ";
			cout<< reducedChainDyToMuMuFSR->GetEntries() * (lumiData / lumiDY) + reducedChainDyToMuMuNonFSR->GetEntries() * (lumiData / lumiDY) + reducedChainTTJets->GetEntries() * (lumiData / lumiTtJets) + reducedChainWJetsToLNu->GetEntries() * (lumiData / lumiWJets);
			cout<<endl<<"DYToMuMuFSR : "<<reducedChainDyToMuMuFSR->GetEntries()<<", "<<reducedChainDyToMuMuFSR->GetEntries() * (lumiData / lumiDY); 
			cout<<endl<<"DYToMuMuNonFSR : "<<reducedChainDyToMuMuNonFSR->GetEntries()<<", "<<reducedChainDyToMuMuNonFSR->GetEntries() * (lumiData / lumiDY);
			cout<<endl<<"TTJets : "<<reducedChainTTJets->GetEntries()<<", "<<reducedChainTTJets->GetEntries() * (lumiData / lumiTtJets);
			cout<<endl<<"WJetsToLNu : "<<reducedChainWJetsToLNu->GetEntries()<<", "<<reducedChainWJetsToLNu->GetEntries() * (lumiData / lumiWJets);
			purity = ( reducedChainDyToMuMuFSR->GetEntries() * (lumiDY / lumiDY) ) / ( reducedChainDyToMuMuFSR->GetEntries() * (lumiDY / lumiDY) + reducedChainDyToMuMuNonFSR->GetEntries() * (lumiDY / lumiDY) + reducedChainTTJets->GetEntries() * (lumiDY / lumiTtJets) + reducedChainWJetsToLNu->GetEntries() * (lumiDY / lumiWJets) );
			cout<<endl<<"Purity : "<< 100 * purity<<" %";

			dyToMuMuFSRChain->Delete();
		        dyToMuMuFSRChain = 0;
		
		        dyToMuMuNonFSRChain->Delete();
		        dyToMuMuNonFSRChain = 0;
		
		        ttJetsChain->Delete();
		        ttJetsChain = 0;
		
		        wJetsToLNuChain->Delete();
		        wJetsToLNuChain = 0;
	
			reducedChainDyToMuMuFSR->Delete();	
			reducedChainDyToMuMuFSR = 0;
			reducedChainDyToMuMuNonFSR->Delete();
			reducedChainDyToMuMuNonFSR = 0;
			reducedChainTTJets->Delete();
			reducedChainTTJets = 0;
			reducedChainWJetsToLNu->Delete();
			reducedChainWJetsToLNu = 0;
                	cout<<endl;
		}
		cout<<endl;
        }
	cout<<endl;
        }

/*
	dataChain->Delete();
	dataChain = 0;
	dataChain1->Delete();
        dataChain1 = 0;
	dataChain2->Delete();
        dataChain2 = 0;
	dataChain3->Delete();
        dataChain3 = 0;
	dataChain4->Delete();
        dataChain4 = 0;
	dataChain5->Delete();
        dataChain5 = 0;
	dataChain6->Delete();
        dataChain6 = 0;
	dataChain7->Delete();
        dataChain7 = 0;

	dyToMuMuFSRChain->Delete();
	dyToMuMuFSRChain = 0;

	dyToMuMuNonFSRChain->Delete();
	dyToMuMuNonFSRChain = 0;
	
	ttJetsChain->Delete();
	ttJetsChain = 0;

	wJetsToLNuChain->Delete();
	wJetsToLNuChain = 0;
*/
	
	return 0;

}






 
