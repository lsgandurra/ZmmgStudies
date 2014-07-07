#include "fitFunctions.h"
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
                cerr << "arguments should be passed : directoryName, dataType, surfaceBinning, fitFunction, i_begin, i_end, j_begin, j_end, iterSub, MmumuOption, binsToFit" <<endl; 
                return 1;

        }	

	string directoryName = "Surface_generation_v1";
	string dataType = "data";
	int surfaceBinning = 5;
	string fitFunction = "voigtian";
	int i_begin = 35;
	int i_end = 39;
        int j_begin = 4;
        int j_end = 39;
	int iterSub = 69;
	string MmumuOption = "MmumuRECO";
	string binsToFit = "binsToFit_data_69.txt";


	if( argc > 1 ) directoryName = argv[1];
	if( argc > 2 ) dataType = argv[2];
	if( argc > 3 )
        {
                std::stringstream ss ( argv[3] );
                ss >> surfaceBinning;
        }	
	if( argc > 4 ) fitFunction = argv[4];
	if( argc > 5 )
        {
                std::stringstream ss ( argv[5] );
                ss >> i_begin;
        }
	if( argc > 6 )
        {
                std::stringstream ss ( argv[6] );
                ss >> i_end;
        }
	if( argc > 7 )
        {
                std::stringstream ss ( argv[7] );
                ss >> j_begin;
        }
	if( argc > 8 )
        {
                std::stringstream ss ( argv[8] );
                ss >> j_end;
        }
	if( argc > 9 )
        {
                std::stringstream ss ( argv[9] );
                ss >> iterSub;
        }
	if( argc > 10 ) MmumuOption = argv[10];
        if( argc > 11 ) binsToFit = argv[11];

	//directoryName += Form("/%s/Binning_%dGeV/MmumuFit/",dataType.c_str(),surfaceBinning);
	cout<<endl<<"directoryName = "<<directoryName<<endl;
	
	system(Form("mkdir -p %s",directoryName.c_str()));

	// --- style of plots --- //
        gROOT->Reset();
        TGaxis::SetMaxDigits(3);
        setTDRStyle();
	gStyle->SetPalette(1,0);
        TCanvas * c1 = new TCanvas("c1", "c1",0,0,600,600);
	TLatex latexLabel;

	TChain * chain = new TChain("miniTree");
	//TChain * reducedChain_temp = 0;
	TChain * reducedChain = 0;
	
	if(dataType == "data")
        {
                if(MmumuOption == "MmumuRECO")
                {
			chain->Add("miniTree_muons_totouples_Run2012A_22Jan2013_v1_noSkim_v1_March_v2_reduced_partALL.root");
			chain->Add("miniTree_muons_totouples_parked_Run2012B_22Jan2013_v1_noSkim_v1_part1_March_v2_reduced_partALL.root");
			chain->Add("miniTree_muons_totouples_parked_Run2012B_22Jan2013_v1_noSkim_v1_part2_March_v2_reduced_partALL.root");
			chain->Add("miniTree_muons_totouples_parked_Run2012C_22Jan2013_v1_noSkim_v1_part1_March_v2_reduced_partALL.root");
			chain->Add("miniTree_muons_totouples_parked_Run2012C_22Jan2013_v1_noSkim_v1_part2_March_v2_reduced_partALL.root");
			chain->Add("miniTree_muons_totouples_parked_Run2012C_22Jan2013_v1_noSkim_v1_part3_March_v2_reduced_partALL.root");
			chain->Add("miniTree_muons_totouples_parked_Run2012D_22Jan2013_v1_noSkim_v1_part1_March_v2_reduced_partALL.root");
			chain->Add("miniTree_muons_totouples_parked_Run2012D_22Jan2013_v1_noSkim_v1_part2_March_v2_reduced_partALL.root");
			chain->Add("miniTree_muons_totouples_parked_Run2012D_22Jan2013_v1_noSkim_v1_part3_March_v2_reduced_partALL.root");
	        }
                if(MmumuOption == "MmumugammaRECO") //FIXME
                {
                        chain->Add("miniTree_Run2012D_PromptReco_v1_NewMuonID_NewSelection_0_injRe0_v6_partALL.root");


                }
        }
	if(dataType == "MC")
        {

                if(MmumuOption == "MmumuRECO" || MmumuOption == "MmumuGEN")
                {
			chain->Add("miniTree_muons_totouple_DYToMuMu_Summer12_S10_reg5_noSkim_v1_March_v2_reduced_partALL.root");
			chain->Add("miniTree_muons_totouple_TTJets_Summer12_S10_reg5_noSkim_v1_March_v2_reduced_partALL.root");					   chain->Add("miniTree_muons_totouples_WJetsToLNu_Summer12_S10_noSkim_v1_part1_March_v2_reduced_partALL.root");
			chain->Add("miniTree_muons_totouples_WJetsToLNu_Summer12_S10_noSkim_v1_part2_March_v2_reduced_partALL.root");
			chain->Add("miniTree_muons_totouples_WJetsToLNu_Summer12_S10_noSkim_v1_part3_March_v2_reduced_partALL.root");	
		}

		if(MmumuOption == "MmumugammaRECO" || MmumuOption == "MmumugammaGEN")
                {
                        chain->Add("miniTree_totouple_DYToMuMu_Summer12_S10_reg5_1_thesis_v1f_recoEnergy_s10_partALL.root");
                        chain->Add("miniTree_totouple_DYToMuMu_Summer12_S10_reg5_2_thesis_v1f_recoEnergy_s10_partALL.root");
                        chain->Add("miniTree_totouple_TTJets_Summer12_S10_reg5_3_thesis_v1f_recoEnergy_s10_partALL.root");
                        chain->Add("miniTree_totouples_WJetsToLNu_Summer12_S10_November2013_3_thesis_v1f_recoEnergy_s10_partALL.root");        
                }

        }

	string fileName = ""; 

	TString cut = "isMM_nonFSR == 1 && MuonL_Pt < 200 && MuonS_Pt < 200"; //FIXME == 1 ou == 2 !

	string fitVariable = "Mmumu";
	string fitVariableName = "M_{#mu#mu}";
	if(MmumuOption == "MmumuRECO")
	{
		fitVariable = "Mmumu";
		fitVariableName = "M_{#mu#mu}";
	}
	if(MmumuOption == "MmumuGEN")
        {
                fitVariable = "Mmumu_Muons_MC";
                fitVariableName = "M_{#mu#mu GEN}";
        }
	if(MmumuOption == "MmumugammaRECO")
        {
                fitVariable = "Mmumugamma";
                fitVariableName = "M_{#mu#mu#gamma}";
        }
	if(MmumuOption == "MmumugammaGEN")
        {
                fitVariable = "Mmumugamma_MMG_MC";
                fitVariableName = "M_{#mu#mu#gamma GEN}";
        }	

	directoryName += Form("/%s/Binning_%dGeV/%s_Fit/",dataType.c_str(),surfaceBinning,MmumuOption.c_str());

	double fitPercentage = 50;
	vector<double> mumuVector;	
	double rangeMin = 0;
	double rangeMax = 0;
	float variableF = 0;	
	vector <double> fitParameters;
	string lumi = "19.6";
	double xMinFitVariable = 70.0;
	double xMaxFitVariable = 110.0;
	int accessNumber = 0;
	double pValue = 0;
	double mean = 0;
	double meanError = 0;
	int xMax = 200; // Change if you want larger surface
        int yMax = 200; // Change if you want larger surface
        int nbBins = xMax / surfaceBinning;
	double fitPercentageF = 1.0;
	int iterDo = 0;
	int iterLoop = 0;
	RooHist* residuals = NULL;
        RooHist* pulls = NULL;
	RooDataSet * dataset = NULL;
	RooDataSet * dataset2 = NULL;

	int nbRows = 0;
	nbRows = rowsNumberInFile(binsToFit.c_str());	
	cout<<endl<<"nbRows = "<<nbRows<<endl;
	pair<double,double> * pairBinsToFit;
	if(nbRows != 0)
	{
		pairBinsToFit = new pair<double,double>[nbRows];
		ifstream txtFile(binsToFit.c_str());
		for(int k = 0; k < nbRows; k++)
		{
			txtFile >> pairBinsToFit[k].first;
			txtFile >> pairBinsToFit[k].second;
			cout<<endl<<"pairBinsToFit[k].first = "<<pairBinsToFit[k].first;
			cout<<endl<<"pairBinsToFit[k].second = "<<pairBinsToFit[k].second;
		}
		txtFile.close();
	}
	

	int inTheLoop = 0;

	TH1D *histo = NULL;
	double xMaxHisto = 91;

	for(int i = i_begin; i <= i_end; i++)
	{
		cout <<endl<< "i = "<<i;
                for(int j = 0; j < nbBins; j++)
		{

			cout <<endl<< "j = "<<j;
			if(i == i_begin && iterLoop == 0) j = j_begin; 
			iterLoop++;


			inTheLoop = 0;
			for(int k = 0; k < nbRows; k++)
			{
				if(pairBinsToFit[k].first == i && pairBinsToFit[k].second == j) 
				{	
					
					inTheLoop = 1;
					if(MmumuOption == "MmumuRECO" || MmumuOption == "MmumuGEN") 
		                        {
						cut = "(isMM_nonFSR == 1"; //FIXME !
						//cut += " && ( hltnames == \"HLT_Mu17_TkMu8_v9\" || hltnames == \"HLT_Mu17_TkMu8_v10\" || hltnames == \"HLT_Mu17_TkMu8_v11\" || hltnames == \"HLT_Mu17_TkMu8_v12\" || hltnames == \"HLT_Mu17_TkMu8_v13\" || hltnames == \"HLT_Mu17_TkMu8_v14\" || hltnames == \"HLT_Mu17_TkMu8_v15\" || hltnames == \"HLT_Mu17_TkMu8_v16\" )";	
		                                //cut += " && ( hltnames == \"HLT_Mu17_Mu8_v1\" || hltnames == \"HLT_Mu17_Mu8_v2\" || hltnames == \"HLT_Mu17_Mu8_v3\" || hltnames == \"HLT_Mu17_Mu8_v4\" || hltnames == \"HLT_Mu17_Mu8_v5\" || hltnames == \"HLT_Mu17_Mu8_v6\" || hltnames == \"HLT_Mu17_Mu8_v7\" || hltnames == \"HLT_Mu17_Mu8_v8\" || hltnames == \"HLT_Mu17_Mu8_v9\" || hltnames == \"HLT_Mu17_Mu8_v10\" || hltnames == \"HLT_Mu17_Mu8_v11\" || hltnames == \"HLT_Mu17_Mu8_v12\" || hltnames == \"HLT_Mu17_Mu8_v13\" || hltnames == \"HLT_Mu17_Mu8_v14\" || hltnames == \"HLT_Mu17_Mu8_v15\" || hltnames == \"HLT_Mu17_Mu8_v16\" || hltnames == \"HLT_Mu17_Mu8_v17\" || hltnames == \"HLT_Mu17_Mu8_v18\" || hltnames == \"HLT_Mu17_Mu8_v19\" || hltnames == \"HLT_Mu17_Mu8_v20\" || hltnames == \"HLT_Mu17_Mu8_v21\" || hltnames == \"HLT_Mu17_Mu8_v22\" || hltnames == \"HLT_Mu17_Mu8_v23\")";
						cut += " && isGoodHLT == 1)*weight_pileUp*weight_Xsection";
						cut += Form(" && MuonL_Pt > %d && MuonL_Pt <= %d && MuonS_Pt > %d && MuonS_Pt <= %d",i * surfaceBinning,(i+1) * surfaceBinning,j * surfaceBinning,(j+1) * surfaceBinning);
		                        }
		                        if(MmumuOption == "MmumugammaRECO" || MmumuOption == "MmumugammaGEN")
		                        {
						cut = "(isJanLooseMMG == 1";
						//cut += " && ( hltnames == \"HLT_Mu17_TkMu8_v9\" || hltnames == \"HLT_Mu17_TkMu8_v10\" || hltnames == \"HLT_Mu17_TkMu8_v11\" || hltnames == \"HLT_Mu17_TkMu8_v12\" || hltnames == \"HLT_Mu17_TkMu8_v13\" || hltnames == \"HLT_Mu17_TkMu8_v14\" || hltnames == \"HLT_Mu17_TkMu8_v15\" || hltnames == \"HLT_Mu17_TkMu8_v16\" )";
						//cut += " && ( hltnames == \"HLT_Mu17_Mu8_v1\" || hltnames == \"HLT_Mu17_Mu8_v2\" || hltnames == \"HLT_Mu17_Mu8_v3\" || hltnames == \"HLT_Mu17_Mu8_v4\" || hltnames == \"HLT_Mu17_Mu8_v5\" || hltnames == \"HLT_Mu17_Mu8_v6\" || hltnames == \"HLT_Mu17_Mu8_v7\" || hltnames == \"HLT_Mu17_Mu8_v8\" || hltnames == \"HLT_Mu17_Mu8_v9\" || hltnames == \"HLT_Mu17_Mu8_v10\" || hltnames == \"HLT_Mu17_Mu8_v11\" || hltnames == \"HLT_Mu17_Mu8_v12\" || hltnames == \"HLT_Mu17_Mu8_v13\" || hltnames == \"HLT_Mu17_Mu8_v14\" || hltnames == \"HLT_Mu17_Mu8_v15\" || hltnames == \"HLT_Mu17_Mu8_v16\" || hltnames == \"HLT_Mu17_Mu8_v17\" || hltnames == \"HLT_Mu17_Mu8_v18\" || hltnames == \"HLT_Mu17_Mu8_v19\" || hltnames == \"HLT_Mu17_Mu8_v20\" || hltnames == \"HLT_Mu17_Mu8_v21\" || hltnames == \"HLT_Mu17_Mu8_v22\" || hltnames == \"HLT_Mu17_Mu8_v23\")";
		                                cut += " && isGoodHLT_Mu17_Mu8 == 1";
						cut += Form(" && ptMuGammaL > %d && ptMuGammaL <= %d && ptMuGammaS > %d && ptMuGammaS <= %d)*weight_pileUp*weight_Xsection",i * surfaceBinning,(i+1) * surfaceBinning,j * surfaceBinning,(j+1) * surfaceBinning);
					}
		
					cout<<endl<<"cut = "<<cut;
					rangeMin = 85.0; //FIXME
					rangeMax = 96.0; //FIXME
					reducedChain = (TChain *) chain->CopyTree(cut);
					//if(chain->GetEntries(cut) != 0) 
					if(reducedChain->GetEntries("Mmumu > 70.0 && Mmumu < 110.0") > 10) //FIXME select the lowest number of entries 
					{
						cout<<endl<<"reducedChain->GetEntries() = "<<reducedChain->GetEntries()<<endl;
						iterDo = 0;
						histo = new TH1D("histo","histo", 80, 70.0, 110.0);
						reducedChain->Draw(Form("%s>>histo",fitVariable.c_str()));
						xMaxHisto = histo->GetXaxis()->GetBinCenter(histo->GetMaximumBin());
						if(xMaxHisto > 100) xMaxHisto = 91.188;
						cout<<endl<<"xMaxHisto = "<<xMaxHisto<<endl;
						do{
							fitPercentageF = fitPercentage / 100.0;
							//rangeEstimator(fitPercentageF, chain, cut, rangeMin, rangeMax, fitVariable, variableF, fitParameters); 
							//fitParameters.push_back(chain->GetEntries(cut));
							fitParameters.push_back(reducedChain->GetEntries());
							if(iterDo == 0)
							{
								rangeMin = xMaxHisto - 6.0;
								rangeMax = xMaxHisto + 6.0;

							}
							//if(iterDo > 0 && rangeMin < 89.0 && rangeMax > 92.0)
							if(iterDo > 0 && rangeMin < (xMaxHisto - 1.5) && rangeMax > (xMaxHisto + 1.5))
							{
								rangeMin += 0.5; //FIXME
		                                        	rangeMax -= 0.5; //FIXME
							}
							cout<<"rangeMin = "<<rangeMin<<", rangeMax = "<<rangeMax<<endl;
							RooRealVar variable(fitVariable.c_str(), fitVariableName.c_str(), rangeMin, rangeMax);	
							RooRealVar weight_pileUp("weight_pileUp", "weight_pileUp", 0.0, 100);
							RooRealVar weight_Xsection("weight_Xsection", "weight_Xsection", 0.0, 100);
							//RooRealVar MuonL_Pt("MuonL_Pt","MuonL_Pt",0,200,"");
							//RooRealVar MuonS_Pt("MuonS_Pt","MuonS_Pt",0,200,"");	
							//RooRealVar isMM_nonFSR("isMM_nonFSR","isMM_nonFSR",0.0,1.0,"");			
			
							RooArgSet ntplVars(variable,weight_pileUp,weight_Xsection);
							dataset = new RooDataSet("dataset", "dataset", reducedChain, ntplVars, "", "weight_Xsection * weight_pileUp");
							// --- Second RooDataSet (only for cosmetic considerations) --- //
							RooRealVar variable2(fitVariable.c_str(), fitVariableName.c_str(), xMinFitVariable, xMaxFitVariable);
							RooArgSet ntplVars2(variable2,weight_pileUp,weight_Xsection);
							dataset2 = new RooDataSet("dataset2", "dataset2", reducedChain, ntplVars2, "", "weight_Xsection * weight_pileUp");	
			
							RooPlot * fitFrame; 
			                		fitFrame = variable.frame(xMinFitVariable,xMaxFitVariable);
							RooBinning b;
							b.setRange(70.0,110.0);
			                        	b.addUniform(80,70.0,110.0);
		
							cout<<endl<<"coucou avant voigtian"<<endl;	
							if(fitFunction == "voigtian") voigtian_surface(dataset, dataset2, variable, fitFrame, b, rangeMin, rangeMax, xMaxHisto, fitParameters);
							cout<<endl<<"coucou apres voigtian"<<endl;
							// --- Compute the chi2 and the associated p-value and save informations in fitParameters --- //
							chiSquare(fitFrame,(char *)"mycurve",(char *)"myhist",fitParameters, 3);
							cout<<endl<<"coucou avant chi2"<<endl;
			
							// --- Draw a nice plot --- //
			                		RooCurve * pdf = fitFrame->getCurve("mycurve");
			                		fitFrame = variable2.frame();
			                		dataset2->plotOn(fitFrame,Name("myhist"),Binning(b));
			                		fitFrame->Draw();
			                		pdf->Draw("SAME");
			
							latexLabel.SetTextSize(0.030);
			                		latexLabel.SetNDC();
			               			latexLabel.DrawLatex(0.13, 0.96, "CMS Preliminary 2012, #sqrt{s} = 8 TeV");
			
			                		if(dataType == "MC") latexLabel.DrawLatex(0.17, 0.88, "Simulation");
			                		if(dataType == "data") latexLabel.DrawLatex(0.17, 0.88, Form("Data, #int L = %s fb^{-1}",lumi.c_str()));
							accessNumber = 2;
							
							if(fitFunction == "voigtian")
							{
								latexLabel.DrawLatex(0.65, 0.88, Form("#color[1]{mean = %f #pm %f}",fitParameters[accessNumber], fitParameters[accessNumber + 1]));	
								latexLabel.DrawLatex(0.65, 0.83, Form("#color[1]{#sigma = %f #pm %f}",fitParameters[accessNumber + 2], fitParameters[accessNumber + 3]));
								latexLabel.DrawLatex(0.65, 0.78, Form("#color[1]{width = %f #pm %f}",fitParameters[accessNumber + 4], fitParameters[accessNumber + 5]));
							}
		
							
							fitPercentage -= 1;
							iterDo++;	
			
							pValue = fitParameters[accessNumber + 8];
							mean = fitParameters[accessNumber];
							meanError = fitParameters[accessNumber + 1];
			
							fitParameters.erase(fitParameters.begin(),fitParameters.end());
				
							cout<<endl<<"*****HERE*****";	
							cout<<endl<<"pValue = "<<pValue;
							cout<<endl<<"mean = "<<mean;
							cout<<endl<<"meanError = "<<meanError;
							fileName = Form("fit_L_%d_%d_S_%d_%d",i * surfaceBinning,(i+1) * surfaceBinning,j * surfaceBinning,(j+1) * surfaceBinning);
		                                	cout<<endl<<"directoryNameBeforeRecord = "<<directoryName<<endl;
							plotsRecording(directoryName, fileName, c1);
			
		
							pdf->Delete();
				
							dataset->Delete();
		                			dataset = 0;
		     
		                			dataset2->Delete();
		                			dataset2 = 0;   
		     
		                			fitFrame->Delete();
		                			fitFrame = 0;
							
							//if(rangeMin == 89.0 && rangeMax == 92.0) break;
							if(rangeMin == (xMaxHisto - 1.5) && rangeMax == (xMaxHisto + 1.5)) break;	
	
						}while(pValue < 0.01 || meanError > (mean / 2.0));
		
						//fileName = Form("fit_L_%d_%d_S_%d_%d",i * surfaceBinning,(i+1) * surfaceBinning,j * surfaceBinning,(j+1) * surfaceBinning);
						//plotsRecording(directoryName, fileName, c1);
		
						c1->Clear();
						
						mumuVector.push_back(mean);
		
					}
					else
	                        	{
	                                	mumuVector.push_back(91.1876);
	                                	cout<<endl<<"else bizarre !!!!!"<<endl;
	                        	}
					histo->Delete();
					histo = 0;
	
				}
			}
			if(inTheLoop == 0)
			{
				mumuVector.push_back(91.1876);
				cout<<endl<<"mumuVector push back"<<endl;
			}
		
			if((j == j_end) && (i == i_end)) 
			{
				cout<<endl<<"coucou dans le  break 1"<<endl;
				cout<<endl<<"iterLoop = "<<iterLoop<<endl;
				break; // We stop the loop after the last fit
			}
		}
	}

	if(nbRows != 0) 
	{
		delete pairBinsToFit;
		pairBinsToFit = 0;
	}


	system(Form("mkdir -p %s",directoryName.c_str()));
	ofstream mumuFile(Form("%sMmumu_%d.txt",directoryName.c_str(),iterSub));

	for(unsigned int i = 0; i < mumuVector.size(); i++)
	{
		mumuFile << mumuVector[i] << endl;
	}

	cout<<endl<<"closing file"<<endl;
	mumuFile.close();

	system(Form("ls -ltrh %s",directoryName.c_str()));

	system(Form("touch Mmumu_%d.done",iterSub));
	return 0;

}
