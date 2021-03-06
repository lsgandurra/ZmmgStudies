#include "functions.h"
#include "setTDRStyle.C"


int main(int argc, char *argv[])
{

	system("echo blahhhhh");

	for(int iarg = 0 ; iarg < argc; iarg++)
        {
                cout << "argv[" << iarg << "]= " << argv[iarg] << endl;
        }

        if( argc == 1 ) 
        {
                cerr << "arguments should be passed : directoryName, dataType, surfaceBinning, nbJobs, option, fitFunction" <<endl; 
                return 1;

        }	

	string directoryName = "Surface_generation_v1";
	string dataType = "MC";
	int surfaceBinning = 5;
	int nbJobs = 200;
	string option = "plot2D_pro";//"plot2D_pro" "plot3D_pro" "plot2D_fit" "plot3D_fit" "MmumuPro" "MmumuFit";
	string fitFunction = "voigtian";
	string MmumuOption = "MmumugammaGEN"; // "MmumuRECO" "MmumuGEN" "MmumugammaRECO" "MmumugammaGEN"

	if( argc > 1 ) directoryName = argv[1];
	if( argc > 2 ) dataType = argv[2];
	if( argc > 3 )
        {
                std::stringstream ss ( argv[3] );
                ss >> surfaceBinning;
        }	
	if( argc > 4 )
        {
                std::stringstream ss ( argv[4] );
                ss >> nbJobs;
        }
	if( argc > 5 ) option = argv[5];
	if( argc > 6 ) fitFunction = argv[6];
	if( argc > 7 ) MmumuOption = argv[7];

	string zAxisTitle = "";
	string yAxisTitle = "";
	string xAxisTitle = "";

	if(MmumuOption == "MmumuRECO") 
	{
		zAxisTitle = "M_{#mu#mu} (GeV)";
		yAxisTitle = "P_{T}^{#mu leading}";
        	xAxisTitle = "P_{T}^{#mu trailing}";
	}
	if(MmumuOption == "MmumuGEN")
        {
                zAxisTitle = "M_{#mu#mu GEN} (GeV)";
                yAxisTitle = "P_{T}^{#mu leading GEN}";
                xAxisTitle = "P_{T}^{#mu trailing GEN}";
        }
	if(MmumuOption == "MmumugammaRECO")
        {
                zAxisTitle = "M_{#mu#mu#gamma} (GeV)";
                yAxisTitle = "P_{T}^{#mu(+#gamma) leading}";
                xAxisTitle = "P_{T}^{#mu(+#gamma) trailing}";
        }
        if(MmumuOption == "MmumugammaGEN")
        {
                zAxisTitle = "M_{#mu#mu#gamma GEN} (GeV)";
                yAxisTitle = "P_{T}^{#mu(+#gamma) leading GEN}";
                xAxisTitle = "P_{T}^{#mu(+#gamma) trailing GEN}";
        }
	


	TChain * chain = new TChain("miniTree");
	TChain * reducedChain = 0;
	if(dataType == "data")
        {
		if(MmumuOption == "MmumuRECO")
		{
			chain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_muons_totouples_Run2012A_22Jan2013_v1_noSkim_v1_March_v2_reduced_partALL.root");
			chain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_muons_totouples_parked_Run2012B_22Jan2013_v1_noSkim_v1_part1_March_v2_reduced_partALL.root");
			chain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_muons_totouples_parked_Run2012B_22Jan2013_v1_noSkim_v1_part2_March_v2_reduced_partALL.root");
			chain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_muons_totouples_parked_Run2012C_22Jan2013_v1_noSkim_v1_part1_March_v2_reduced_partALL.root");
			chain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_muons_totouples_parked_Run2012C_22Jan2013_v1_noSkim_v1_part2_March_v2_reduced_partALL.root");
			chain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_muons_totouples_parked_Run2012C_22Jan2013_v1_noSkim_v1_part3_March_v2_reduced_partALL.root");
			chain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_muons_totouples_parked_Run2012D_22Jan2013_v1_noSkim_v1_part1_March_v2_reduced_partALL.root");
			chain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_muons_totouples_parked_Run2012D_22Jan2013_v1_noSkim_v1_part2_March_v2_reduced_partALL.root");
			chain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_muons_totouples_parked_Run2012D_22Jan2013_v1_noSkim_v1_part3_March_v2_reduced_partALL.root");
		}
		if(MmumuOption == "MmumugammaRECO")//FIXME
		{
			chain->Add("/sps/cms/sgandurr/CMSSW_5_3_6_RECO_5_3_3_v4/src/cvs_developpment/Selection_NewMuID/miniTree_Run2012A_13Jul2012_v1_NewMuonID_NewSelection_0_injRe0_v6_partALL.root");


		}
	}
	if(dataType == "MC")
        {
		
		if(MmumuOption == "MmumuRECO" || MmumuOption == "MmumuGEN")
		{
			chain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_muons_totouple_DYToMuMu_Summer12_S10_reg5_noSkim_v1_March_v2_reduced_partALL.root");
			chain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_muons_totouple_TTJets_Summer12_S10_reg5_noSkim_v1_March_v2_reduced_partALL.root");
			chain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_muons_totouples_WJetsToLNu_Summer12_S10_noSkim_v1_part1_March_v2_reduced_partALL.root");
                        chain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_muons_totouples_WJetsToLNu_Summer12_S10_noSkim_v1_part2_March_v2_reduced_partALL.root");
                        chain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_muons_totouples_WJetsToLNu_Summer12_S10_noSkim_v1_part3_March_v2_reduced_partALL.root");
		}
		if(MmumuOption == "MmumugammaRECO" || MmumuOption == "MmumugammaGEN")
                {
			chain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_totouple_DYToMuMu_Summer12_S10_reg5_1_thesis_v1f_recoEnergy_s10_partALL.root");
			chain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_totouple_DYToMuMu_Summer12_S10_reg5_2_thesis_v1f_recoEnergy_s10_partALL.root");
			chain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_totouple_TTJets_Summer12_S10_reg5_3_thesis_v1f_recoEnergy_s10_partALL.root");
			chain->Add("/sps/cms/sgandurr/CMSSW_5_3_11_patch3_RECO_5_3_11_v1_hggpaperV5/src/ZmmgStudies/Selection/miniTree_totouples_WJetsToLNu_Summer12_S10_November2013_3_thesis_v1f_recoEnergy_s10_partALL.root");					
		}
	}
	
	gROOT->Reset();
        TGaxis::SetMaxDigits(3);
        setTDRStyle();
	gStyle->SetPalette(1,0);
	//gStyle->SetPalette(1);
	TCanvas *c1 = new TCanvas("c1", "c1",0,0,600,600);

	directoryName += Form("/%s/Binning_%dGeV/%s/%s/",dataType.c_str(),surfaceBinning,option.c_str(),MmumuOption.c_str());
	
	string fileName = Form("%s",option.c_str()); //FIXME

	TString cut = "isMM_nonFSR == 1 && MuonL_Pt < 200 && MuonS_Pt < 200";

	int xMax = 100;
	int yMax = 100;

	int nbBins = xMax / surfaceBinning;

	int i_begin = 0;
	int i_end = 0;
	int j_begin = 0;
        int j_end = 0;
	int nbFits = 0;
	int nbFitsMax = ( nbBins * nbBins / 2 ) / nbJobs; 
	int iterSub = 0;
	string chainSub = "";

	system(Form("rm binsToFit_%s_%s_%s*.txt",dataType.c_str(),fitFunction.c_str(),MmumuOption.c_str())); //suppression of old bins to fit.
	if(option == "MmumuFit") 
	{
		//system(Form("rm fileTosubmit_%s.sh",dataType.c_str())); //suppression of old submission file.
		system(Form("rm fileTosubmit_%s_%s_%s.sh",dataType.c_str(),fitFunction.c_str(),MmumuOption.c_str()));
		//ofstream temp_file(Form("fileTosubmit_%s.sh",dataType.c_str()));
		ofstream temp_file(Form("fileTosubmit_%s_%s_%s.sh",dataType.c_str(),fitFunction.c_str(),MmumuOption.c_str()));
		temp_file << "#! /usr/local/bin/bash -l" <<endl;
		temp_file.close();
	}
	cout<<endl<<"nbFitsMax = "<<nbFitsMax<<", nbBins = "<<nbBins<<endl;
	if(option == "MmumuFit")
	{
		if(MmumuOption == "MmumuRECO" || MmumuOption == "MmumuGEN")
		{
			cut = "(isMM_nonFSR == 1"; //FIXME == 1 ou == 2 !
			//cut += " && ( hltnames == \"HLT_Mu17_TkMu8_v9\" || hltnames == \"HLT_Mu17_TkMu8_v10\" || hltnames == \"HLT_Mu17_TkMu8_v11\" || hltnames == \"HLT_Mu17_TkMu8_v12\" || hltnames == \"HLT_Mu17_TkMu8_v13\" || hltnames == \"HLT_Mu17_TkMu8_v14\" || hltnames == \"HLT_Mu17_TkMu8_v15\" || hltnames == \"HLT_Mu17_TkMu8_v16\" )";
			//cut += " && ( hltnames == \"HLT_Mu17_Mu8_v1\" || hltnames == \"HLT_Mu17_Mu8_v2\" || hltnames == \"HLT_Mu17_Mu8_v3\" || hltnames == \"HLT_Mu17_Mu8_v4\" || hltnames == \"HLT_Mu17_Mu8_v5\" || hltnames == \"HLT_Mu17_Mu8_v6\" || hltnames == \"HLT_Mu17_Mu8_v7\" || hltnames == \"HLT_Mu17_Mu8_v8\" || hltnames == \"HLT_Mu17_Mu8_v9\" || hltnames == \"HLT_Mu17_Mu8_v10\" || hltnames == \"HLT_Mu17_Mu8_v11\" || hltnames == \"HLT_Mu17_Mu8_v12\" || hltnames == \"HLT_Mu17_Mu8_v13\" || hltnames == \"HLT_Mu17_Mu8_v14\" || hltnames == \"HLT_Mu17_Mu8_v15\" || hltnames == \"HLT_Mu17_Mu8_v16\" || hltnames == \"HLT_Mu17_Mu8_v17\" || hltnames == \"HLT_Mu17_Mu8_v18\" || hltnames == \"HLT_Mu17_Mu8_v19\" || hltnames == \"HLT_Mu17_Mu8_v20\" || hltnames == \"HLT_Mu17_Mu8_v21\" || hltnames == \"HLT_Mu17_Mu8_v22\" || hltnames == \"HLT_Mu17_Mu8_v23\"))*weight_pileUp*weight_Xsection";
			cut += " && isGoodHLT == 1)*weight_pileUp*weight_Xsection";
		}	
		if(MmumuOption == "MmumugammaRECO" || MmumuOption == "MmumugammaGEN") 
                {
                	cut = "(isJanLooseMMG == 1";
                        //cut += " && ( hltnames == \"HLT_Mu17_TkMu8_v9\" || hltnames == \"HLT_Mu17_TkMu8_v10\" || hltnames == \"HLT_Mu17_TkMu8_v11\" || hltnames == \"HLT_Mu17_TkMu8_v12\" || hltnames == \"HLT_Mu17_TkMu8_v13\" || hltnames == \"HLT_Mu17_TkMu8_v14\" || hltnames == \"HLT_Mu17_TkMu8_v15\" || hltnames == \"HLT_Mu17_TkMu8_v16\" )";
                	//cut += " && ( hltnames == \"HLT_Mu17_Mu8_v1\" || hltnames == \"HLT_Mu17_Mu8_v2\" || hltnames == \"HLT_Mu17_Mu8_v3\" || hltnames == \"HLT_Mu17_Mu8_v4\" || hltnames == \"HLT_Mu17_Mu8_v5\" || hltnames == \"HLT_Mu17_Mu8_v6\" || hltnames == \"HLT_Mu17_Mu8_v7\" || hltnames == \"HLT_Mu17_Mu8_v8\" || hltnames == \"HLT_Mu17_Mu8_v9\" || hltnames == \"HLT_Mu17_Mu8_v10\" || hltnames == \"HLT_Mu17_Mu8_v11\" || hltnames == \"HLT_Mu17_Mu8_v12\" || hltnames == \"HLT_Mu17_Mu8_v13\" || hltnames == \"HLT_Mu17_Mu8_v14\" || hltnames == \"HLT_Mu17_Mu8_v15\" || hltnames == \"HLT_Mu17_Mu8_v16\" || hltnames == \"HLT_Mu17_Mu8_v17\" || hltnames == \"HLT_Mu17_Mu8_v18\" || hltnames == \"HLT_Mu17_Mu8_v19\" || hltnames == \"HLT_Mu17_Mu8_v20\" || hltnames == \"HLT_Mu17_Mu8_v21\" || hltnames == \"HLT_Mu17_Mu8_v22\" || hltnames == \"HLT_Mu17_Mu8_v23\"))*weight_pileUp*weight_Xsection";
			cut += " && isGoodHLT_Mu17_Mu8 == 1)*weight_pileUp*weight_Xsection";
		}	

		reducedChain = (TChain *) chain->CopyTree(cut);
		for(int i = 0; i < nbBins; i++)
		{
			for(int j = 0; j < nbBins; j++)
                	{
				if(MmumuOption == "MmumuRECO" || MmumuOption == "MmumuGEN") 
				{	
					cut = Form("(Mmumu > 70.0 && Mmumu < 110.0 && MuonL_Pt > %d && MuonL_Pt <= %d && MuonS_Pt > %d && MuonS_Pt <= %d)*weight_pileUp*weight_Xsection",i * surfaceBinning,(i+1) * surfaceBinning,j * surfaceBinning,(j+1) * surfaceBinning); //Mmumu or Mmumu_Muons_MC should be approximately similar
					 
				}
				if(MmumuOption == "MmumugammaRECO" || MmumuOption == "MmumugammaGEN") 
                                {
                                        cut = Form("(Mmumu > 70.0 && Mmumu < 110.0 && ptMuGammaL > %d && ptMuGammaL <= %d && ptMuGammaS > %d && ptMuGammaS <= %d)*weight_pileUp*weight_Xsection",i * surfaceBinning,(i+1) * surfaceBinning,j * surfaceBinning,(j+1) * surfaceBinning);
				}
			
				//if(reducedChain->GetEntries(cut) != 0) nbFits++;
				if(reducedChain->GetEntries(cut) > 10) 
				{
					nbFits++;
					//ofstream txtFile(Form("binsToFit_%s_%d.txt",dataType.c_str(),iterSub),ios::app);//FIXME
					ofstream txtFile(Form("binsToFit_%s_%s_%s_%d.txt",dataType.c_str(),fitFunction.c_str(),MmumuOption.c_str(),iterSub),ios::app);
					txtFile << i << " " << j <<endl;
					txtFile.close();

				}
				if(nbFits == nbFitsMax || (i == nbBins-1 && j == nbBins-1))
				{
					nbFits = 0;
					i_end = i;
					j_end = j;
				
					directoryName = argv[1];	
					cout<<endl<<"job n : "<<iterSub + 1;
				
					ofstream fileTosubmit(Form("fileTosubmit_%s_%s_%s.sh",dataType.c_str(),fitFunction.c_str(),MmumuOption.c_str()),ios::app);
					
					chainSub = Form("qsub batchJob_surface.sh %s %s %d %s %d %d %d %d %d %s binsToFit_%s_%s_%s_%d.txt",directoryName.c_str(),dataType.c_str(),surfaceBinning,fitFunction.c_str(), i_begin, i_end, j_begin, j_end, iterSub, MmumuOption.c_str(),dataType.c_str(),fitFunction.c_str(),MmumuOption.c_str(),iterSub) ;			
					fileTosubmit << chainSub << endl;	

					//system(Form("qsub batchJob_surface.sh %s %s %d %s %d %d %d %d %d %s",directoryName.c_str(),dataType.c_str(),surfaceBinning,fitFunction.c_str(), i_begin, i_end, j_begin, j_end, iterSub, MmumuOption.c_str()));	
					cout<<endl<<chainSub.c_str()<<endl;
					fileTosubmit.close();
				
					iterSub++;

					if(j != (nbBins - 1)) 
					{
						i_begin = i; 
						j_begin = j + 1; 
					}
					else
					{
						i_begin = i + 1;
						j_begin = 0;
					}

				}

			}
		}

		return 0;
	}	
	
	
	if(option == "MmumuPro")
        {

	}

	TH2D *h1 = NULL;
	TH2D *h2 = NULL;
	TH2D *h3 = NULL;
	TH2D *h4 = NULL;

	h1 = new TH2D("h1", "h1", nbBins, 0, 100, nbBins, 0, 200); //FIXME 100 ?
	h2 = new TH2D("h2", "h2", nbBins, 0, 100, nbBins, 0, 100); //FIXME 100 ?

	if(option == "plot2D_pro" || option == "plot3D_pro")
        {
		if(MmumuOption == "MmumuRECO")
		{	
			cut = "Mmumu*(isMM_nonFSR == 1)"; //FIXME HLT
			//cut = "Mmumu*(isMM == 1)";
			chain->Draw("MuonL_Pt:MuonS_Pt>>h1", cut, "LEGO2"); //FIXME Gen
			h3 = (TH2D*)gDirectory->Get("h1");
			cut = "isMM_nonFSR == 1"; //FIXME HLT
			//cut = "isMM == 1";
			chain->Draw("MuonL_Pt:MuonS_Pt>>h2", cut, "LEGO2");
		}
		if(MmumuOption == "MmumuGEN")
                {
                        cut = "Mmumu_Muons_MC*(isMM_nonFSR == 1)"; //FIXME HLT
                        //cut = "Mmumu*(isMM == 1)";
                        chain->Draw("MuonL_Pt:MuonS_Pt>>h1", cut, "LEGO2"); //FIXME Gen
                        h3 = (TH2D*)gDirectory->Get("h1");
                        cut = "isMM_nonFSR == 1"; //FIXME HLT
                        //cut = "isMM == 1";
                        chain->Draw("MuonL_Pt:MuonS_Pt>>h2", cut, "LEGO2");
                }	
		if(MmumuOption == "MmumugammaRECO")
                {
			cut = "Mmumugamma * ( isJanLooseMMG )";
			chain->Draw("((sqrt(pow(MuonN_Px + Photon_Px,2) + pow(MuonN_Py + Photon_Py,2))>MuonF_Pt)?sqrt(pow(MuonN_Px + Photon_Px,2) + pow(MuonN_Py + Photon_Py,2)):MuonF_Pt):((sqrt(pow(MuonN_Px + Photon_Px,2) + pow(MuonN_Py + Photon_Py,2))<MuonF_Pt)?sqrt(pow(MuonN_Px + Photon_Px,2) + pow(MuonN_Py + Photon_Py,2)):MuonF_Pt) >>h1", cut, "LEGO2");
			h3 = (TH2D*)gDirectory->Get("h1");
			cut = "isJanLooseMMG";
			chain->Draw("((sqrt(pow(MuonN_Px + Photon_Px,2) + pow(MuonN_Py + Photon_Py,2))>MuonF_Pt)?sqrt(pow(MuonN_Px + Photon_Px,2) + pow(MuonN_Py + Photon_Py,2)):MuonF_Pt):((sqrt(pow(MuonN_Px + Photon_Px,2) + pow(MuonN_Py + Photon_Py,2))<MuonF_Pt)?sqrt(pow(MuonN_Px + Photon_Px,2) + pow(MuonN_Py + Photon_Py,2)):MuonF_Pt) >>h2", cut, "LEGO2");

		}
		if(MmumuOption == "MmumugammaGEN")
                {
                        //cut = "Mmumugamma_MMG_MC * ( isJanLooseMMG )";
                        //chain->Draw("((sqrt(pow(MuonN_Px + Photon_Px,2) + pow(MuonN_Py + Photon_Py,2))>MuonF_Pt)?sqrt(pow(MuonN_Px + Photon_Px,2) + pow(MuonN_Py + Photon_Py,2)):MuonF_Pt):((sqrt(pow(MuonN_Px + Photon_Px,2) + pow(MuonN_Py + Photon_Py,2))<MuonF_Pt)?sqrt(pow(MuonN_Px + Photon_Px,2) + pow(MuonN_Py + Photon_Py,2)):MuonF_Pt) >>h1", cut, "LEGO2");
                        //h3 = (TH2D*)gDirectory->Get("h1");
                        cut = "isJanLooseMMG";
                        chain->Draw("((sqrt(pow(MuonN_Px + Photon_Px,2) + pow(MuonN_Py + Photon_Py,2))>MuonF_Pt)?sqrt(pow(MuonN_Px + Photon_Px,2) + pow(MuonN_Py + Photon_Py,2)):MuonF_Pt):((sqrt(pow(MuonN_Px + Photon_Px,2) + pow(MuonN_Py + Photon_Py,2))<MuonF_Pt)?sqrt(pow(MuonN_Px + Photon_Px,2) + pow(MuonN_Py + Photon_Py,2)):MuonF_Pt) >>h2", cut, "LEGO2");
		
			cut = "(isJanLooseMMG";
                        cut += " && isGoodHLT_Mu17_Mu8 == 1)*weight_pileUp*weight_Xsection";
                        chain->Draw("((sqrt(pow(MuonN_Px + Photon_Px,2) + pow(MuonN_Py + Photon_Py,2))>MuonF_Pt)?sqrt(pow(MuonN_Px + Photon_Px,2) + pow(MuonN_Py + Photon_Py,2)):MuonF_Pt):((sqrt(pow(MuonN_Px + Photon_Px,2) + pow(MuonN_Py + Photon_Py,2))<MuonF_Pt)?sqrt(pow(MuonN_Px + Photon_Px,2) + pow(MuonN_Py + Photon_Py,2)):MuonF_Pt) >>h1",cut);
			h3 = (TH2D*)gDirectory->Get("h1");                

		}
		
		h4 = (TH2D*)gDirectory->Get("h2");
		//h3->Divide(h4); //FIXME !!
		c1->Clear();
        	
	}

	string tempChain = "";
	double temp_number = 0;
	int nbJobsReal = 0;
	if(option == "plot2D_fit" || option == "plot3D_fit")
        {
		if(MmumuOption == "MmumuRECO" || MmumuOption == "MmumuGEN")
		{
			cut = "(isMM_nonFSR == 1"; //FIXME == 1 ou == 2;
                        //cut += " && ( hltnames == \"HLT_Mu17_TkMu8_v9\" || hltnames == \"HLT_Mu17_TkMu8_v10\" || hltnames == \"HLT_Mu17_TkMu8_v11\" || hltnames == \"HLT_Mu17_TkMu8_v12\" || hltnames == \"HLT_Mu17_TkMu8_v13\" || hltnames == \"HLT_Mu17_TkMu8_v14\" || hltnames == \"HLT_Mu17_TkMu8_v15\" || hltnames == \"HLT_Mu17_TkMu8_v16\" )";
			//cut += " && ( hltnames == \"HLT_Mu17_Mu8_v1\" || hltnames == \"HLT_Mu17_Mu8_v2\" || hltnames == \"HLT_Mu17_Mu8_v3\" || hltnames == \"HLT_Mu17_Mu8_v4\" || hltnames == \"HLT_Mu17_Mu8_v5\" || hltnames == \"HLT_Mu17_Mu8_v6\" || hltnames == \"HLT_Mu17_Mu8_v7\" || hltnames == \"HLT_Mu17_Mu8_v8\" || hltnames == \"HLT_Mu17_Mu8_v9\" || hltnames == \"HLT_Mu17_Mu8_v10\" || hltnames == \"HLT_Mu17_Mu8_v11\" || hltnames == \"HLT_Mu17_Mu8_v12\" || hltnames == \"HLT_Mu17_Mu8_v13\" || hltnames == \"HLT_Mu17_Mu8_v14\" || hltnames == \"HLT_Mu17_Mu8_v15\" || hltnames == \"HLT_Mu17_Mu8_v16\" || hltnames == \"HLT_Mu17_Mu8_v17\" || hltnames == \"HLT_Mu17_Mu8_v18\" || hltnames == \"HLT_Mu17_Mu8_v19\" || hltnames == \"HLT_Mu17_Mu8_v20\" || hltnames == \"HLT_Mu17_Mu8_v21\" || hltnames == \"HLT_Mu17_Mu8_v22\" || hltnames == \"HLT_Mu17_Mu8_v23\")";
			cut += " && isGoodHLT == 1)*weight_pileUp*weight_Xsection";
			chain->Draw("MuonL_Pt:MuonS_Pt>>h1",cut);

		}
		if(MmumuOption == "MmumugammaGEN")
                {
			cut = "(isJanLooseMMG";
			cut += " && isGoodHLT_Mu17_Mu8 == 1)*weight_pileUp*weight_Xsection";
			chain->Draw("((sqrt(pow(MuonN_Px + Photon_Px,2) + pow(MuonN_Py + Photon_Py,2))>MuonF_Pt)?sqrt(pow(MuonN_Px + Photon_Px,2) + pow(MuonN_Py + Photon_Py,2)):MuonF_Pt):((sqrt(pow(MuonN_Px + Photon_Px,2) + pow(MuonN_Py + Photon_Py,2))<MuonF_Pt)?sqrt(pow(MuonN_Px + Photon_Px,2) + pow(MuonN_Py + Photon_Py,2)):MuonF_Pt) >>h1",cut);
		}
		
		h4 = (TH2D*)gDirectory->Get("h1");

		directoryName = argv[1];
		//directoryName += Form("/%s/Binning_%dGeV/MmumuFit/",dataType.c_str(),surfaceBinning);
		directoryName += Form("/%s/Binning_%dGeV/%s_Fit/",dataType.c_str(),surfaceBinning,MmumuOption.c_str());

		if(FileExists(Form("%sMmumu.txt",directoryName.c_str()))) //Mmumu.txt creation
		{
			cout << "Mmumu.txt already exists" <<endl;
			cout << "Regeneration..."<<endl;
			system(Form("rm %sMmumu.txt",directoryName.c_str()));

		}
		else
		{
			cout << "Mmumu.txt doesn't exists" <<endl;
			cout << "Generation..."<<endl;
		}

		system("rm temp_file");
		system(Form("ls %sMmumu*.txt >> temp_file",directoryName.c_str()));
		nbJobsReal = rowsNumberInFile("temp_file");
        	
		for(int i = 0; i < nbJobsReal; i++)
        	{
			//tempChain = Form("cat %sMmumu_%d.txt >> %sMmumu.txt",directoryName.c_str(),i,directoryName.c_str());	
			//cout<<endl<<tempChain;
                	system(Form("cat %sMmumu_%d.txt >> %sMmumu.txt",directoryName.c_str(),i,directoryName.c_str()));

        	}    

		h3 = new TH2D("h4", "h4", nbBins, 0, 200, nbBins, 0, 200);
		ifstream fitFile(Form("%sMmumu.txt",directoryName.c_str()));
		for(int i = 1; i <= nbBins; i++) //ATTENTION 0 or 1 to begin ??
                {
                        for(int j = 1; j <= nbBins; j++) //ATTENTION 0 or 1 to begin ??
                        {
				//cout<<endl<<"fit number : "<<(i-1)*40 + j;
                                fitFile>>temp_number;
                                if(h4->GetBinContent(j,i) != 0)
                                {
                                        h3->SetBinContent(j,i,temp_number);
                                }
                        }
                }    
                fitFile.close();
        }


	if(option == "plot2D_pro" || option == "plot2D_fit") h3->Draw("COLZ");
        if(option == "plot3D_pro" || option == "plot3D_fit") h3->Draw("LEGO2");
	gPad->Update(); 
	
	h3->GetYaxis()->SetTitle(yAxisTitle.c_str());
	h3->GetXaxis()->SetTitle(xAxisTitle.c_str());
	h3->GetZaxis()->SetTitle(zAxisTitle.c_str());
	h3->GetYaxis()->SetTitleOffset(1.55);
	h3->GetXaxis()->SetTitleOffset(1.40);
	if(option == "plot2D_fit") 
	{
		h3->GetXaxis()->SetRangeUser(0,100);
		h3->GetYaxis()->SetRangeUser(0,100);
		h3->GetZaxis()->SetRangeUser(86,93);
	}

	TPaletteAxis *palette = NULL;
	if(option == "plot2D_pro" || option == "plot2D_fit")
        {
		//c1->Update(); 
                //TPaletteAxis palette(96,0,100,100,h3);
		palette = (TPaletteAxis*)h3->GetListOfFunctions()->FindObject("palette");
		palette->SetX1NDC(0.92);
		palette->SetX2NDC(0.95);
		palette->SetLabelSize(0.025);
		palette->SetLabelOffset(0.007);
		palette->SetTitleOffset(0.007);
		//palette->SetY2NDC(0.7);
		//palette = new TPaletteAxis(96,0,100,100,h3);
		/*palette->SetTitleSize(0.02);
		palette->SetX1NDC(90);
		palette->SetX2NDC(100);
		palette->SetY1NDC(10);
		palette->SetY2NDC(100);
		palette->SetLabelColor(1);
                palette->SetLabelFont(42);
                palette->SetLabelOffset(0.007);
                palette->SetLabelSize(0.02);
                palette->SetTitleOffset(1);
                palette->SetTitleSize(0.02);
                palette->SetFillColor(100);
                palette->SetFillStyle(1001);
                */
		//h3->GetListOfFunctions()->Add(palette,"br");
		
        }

	plotsRecording(directoryName, fileName, c1);
	c1->Clear();
	h3->Delete();
	h3 = 0;
	h4->Delete();
       	h4 = 0;


	return 0;

}
