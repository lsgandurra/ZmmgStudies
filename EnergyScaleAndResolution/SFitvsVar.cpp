//Do TGraphs of the selected s values versus a given variable 
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
                cerr << "arguments should be passed : directoryName, dataType, fitVariable, cutVariable, injectedResolution" <<endl; 
                return 1;

        }    

        string directoryName = "Results_v8_surface";
        string dataType = "MC"; //data, MC, both, bothPlusTrue
        string fitVariable = "mmg_s";
	string cutVariable = "Photon_Et";
	string binFileName = "LimitesAllPtOneBin.txt";
	string fitFunction = "voigtian";
	string injectedResolution = "0";

	if( argc > 1 ) directoryName = argv[1];
	if( argc > 2 ) dataType = argv[2];
	if( argc > 3 ) fitVariable = argv[3];
	if( argc > 4 ) cutVariable = argv[4];
	if( argc > 5 ) binFileName = argv[5];
	if( argc > 6 ) fitFunction = argv[6];
	if( argc > 7 ) injectedResolution = argv[7];


	double xMinCutVariable, xMaxCutVariable, xMin, xMax;
	string cutVariableName = "";
	if(cutVariable == "Photon_Et") 
        {	
		cutVariableName = "E^{#gamma} (GeV)";
                xMinCutVariable = 0.0; 
                xMaxCutVariable = 250.0;
		xMin = 0;
		xMax = 100;
	}

	string yTitle = "";
	if(fitVariable == "mmg_s") yTitle = "s (%)";
	if(fitVariable == "mmg_s_MZ_Surface") yTitle = "s_{Surface} (%)";

	int nBins = rowsNumberInFile(binFileName.c_str()) - 1;
	int nbRows = 0;
	int accessNumber = 0;
        double temp_number = 0;

	double * binsTab = new double[nBins+1];
	
	ifstream binFile(binFileName.c_str());
	for(int i = 0; i < nBins+1; i++)
	{
		binFile >> binsTab[i];	
	}
	binFile.close();	

	string fitVariable2 = "";
	string eta = "";
	string r9 = "";
	string fileName = "";
	string directoryNameSave = "";
	string fileNameSave = "";

	double * x = new double[nBins];
        double * xl = new double[nBins];
        double * xr = new double[nBins];
        double * y = new double[nBins];
        double * yl = new double[nBins];
        double * yr = new double[nBins];

	double * x_mc = new double[nBins];
        double * xl_mc = new double[nBins];
        double * xr_mc = new double[nBins];
        double * y_mc = new double[nBins];
        double * yl_mc = new double[nBins];
        double * yr_mc = new double[nBins];

	double * x_true = new double[nBins];
        double * xl_true = new double[nBins];
        double * xr_true = new double[nBins];
        double * y_true = new double[nBins];
        double * yl_true = new double[nBins];
        double * yr_true = new double[nBins];

	gROOT->Reset();
        TGaxis::SetMaxDigits(3);
        setTDRStyle();

	TCanvas *c1 = new TCanvas("c1", "c1",0,0,600,600);
	TLegend * leg2 = NULL;
	TLegend * leg3 = NULL; 

	TMultiGraph * mg2 = NULL;
	TMultiGraph * mg3 = NULL;
	TGraphAsymmErrors * dataGraph = NULL;
	TGraphAsymmErrors * mcGraph = NULL; 
	TGraphAsymmErrors * trueGraph = NULL;


	vector <double> fitParameters;



		for(int i = 0; i < 2; i++)
		{
			if(i == 0) eta = "Barrel";
			if(i == 1) eta = "Endcaps";
	
			for(int j = 0; j < 3; j++)
			{
				if(j == 0) r9 = "low";
				if(j == 1) r9 = "high";
				if(j == 2) r9 = "all";
					
				for(int h = 0; h < 3; h++)
        			{
			                if(h == 0 && fitVariable == "mmg_s")
			                {
			                        dataType = "data";
			                        fitVariable2 = "mmg_s";
						fitFunction = "voigtian";
			                }
			                if(h == 1 && fitVariable == "mmg_s")
			                {
			                        dataType = "MC";
			                        fitVariable2 = "mmg_s";
						fitFunction = "voigtian";
			                }
			                if(h == 2 && fitVariable == "mmg_s")
			                {
			                        dataType = "MC";
			                        fitVariable2 = "mmg_s_true";
						fitFunction = "cruijff";
			                }
			                if(h == 0 && fitVariable == "mmg_s_MZ_Surface")
			                {
			                        dataType = "data";
			                        fitVariable2 = "mmg_s_MZ_Surface";
						fitFunction = "voigtian";
			                }
			                if(h == 1 && fitVariable == "mmg_s_MZ_Surface")
			                {
			                        dataType = "MC";
			                        fitVariable2 = "mmg_s_MZ_Surface";
						fitFunction = "voigtian";
			                }
			                if(h == 2 && fitVariable == "mmg_s_MZ_Surface")
			                {
			                        dataType = "MC";
			                        fitVariable2 = "mmg_s_true";
						fitFunction = "cruijff";
			                }
				

	
					for(int k = 0; k < nBins; k++)
					{
		
						fileName = Form("%s/InjectedResolution_%sPercent/%s/Selected_Fits/Bin_%d/%s_%s_%sR9_%s/fitsInformationsRaw.txt",directoryName.c_str(),injectedResolution.c_str(),dataType.c_str(),k,fitVariable2.c_str(),eta.c_str(),r9.c_str(),fitFunction.c_str());
				
						nbRows = rowsNumberInFile(fileName.c_str());
		
						ifstream fitParametersFile(fileName.c_str());
		
		                		if(fitParametersFile)    
		                		{
		                        		cout << "OK : The file is open." << endl;
		                		}
		                		else
		                		{
		                        		cout << "ERROR: Impossible to open the file." << endl;
		                		}
		     
		                		for(int l = 0; l < nbRows; l++)
		                		{
		                        		fitParametersFile >> temp_number;
		                        		fitParameters.push_back(temp_number);
		                		}
		
						fitParametersFile.close();
		
						accessNumber = k * ( 2 + 6 + fitParameters[2] * 2 ) + 1 + 2; 
						if(h == 0)
						{
							y[k] = 100 * fitParameters[accessNumber];
							yl[k] = yr[k] = 100 * fitParameters[accessNumber + 1];
						}
						if(h == 1)
	                                        {
	                                                y_mc[k] = 100 * fitParameters[accessNumber];
	                                                yl_mc[k] = yr_mc[k] = 100 * fitParameters[accessNumber + 1];
	                                        }
						if(h == 2)
	                                        {
	                                                y_true[k] = 100 * fitParameters[accessNumber];
	                                                yl_true[k] = yr_true[k] = 100 * fitParameters[accessNumber + 1];
	                                        }
						accessNumber = k * ( 2 + 6 + fitParameters[2] * 2 ) + fitParameters[2] * 2 + 7;				
						if(h == 0) x[k] = fitParameters[accessNumber];
						if(h == 1) x_mc[k] = fitParameters[accessNumber];
						if(h == 2) x_true[k] = fitParameters[accessNumber];	
						if(h == 0)
						{
							xl[k] = x[k] - binsTab[k]; 		
							xr[k] = binsTab[k + 1] - x[k];
						}
						if(h == 1)
	                                        {
	                                        	xl_mc[k] = x_mc[k] - binsTab[k];
	                                                xr_mc[k] = binsTab[k + 1] - x_mc[k];
	                                        }
						if(h == 2)
	                                        {
	                                        	xl_true[k] = x_true[k] - binsTab[k];
	                                                xr_true[k] = binsTab[k + 1] - x_true[k];
	                                        }
					
						fitParameters.erase(fitParameters.begin(),fitParameters.end());
					
					}
	
			
				}
					
				directoryNameSave = Form("%s/InjectedResolution_%sPercent/CombinedGraphs/%s_%sR9/",directoryName.c_str(),injectedResolution.c_str(),eta.c_str(),r9.c_str());			
				system(Form("mkdir -p %s",directoryNameSave.c_str()));
				//Continue the code from here 

                                mg2 = new TMultiGraph();
				mg3 = new TMultiGraph();
				dataGraph = new TGraphAsymmErrors(nBins,x,y,xl,xr,yl,yr);
				mcGraph = new TGraphAsymmErrors(nBins,x_mc,y_mc,xl_mc,xr_mc,yl_mc,yr_mc);
				trueGraph = new TGraphAsymmErrors(nBins,x_true,y_true,xl_true,xr_true,yl_true,yr_true);
				
				dataGraph->SetMarkerColor(kViolet-3);
				mcGraph->SetMarkerColor(kAzure+7);
				trueGraph->SetMarkerColor(kBlue);
				dataGraph->SetLineColor(kViolet-3);
                                mcGraph->SetLineColor(kAzure+7);
                                trueGraph->SetLineColor(kBlue);
				dataGraph->SetMarkerStyle(20);
				dataGraph->SetMarkerSize(0.6);
				mcGraph->SetMarkerStyle(20);
                                mcGraph->SetMarkerSize(0.6);
				trueGraph->SetMarkerStyle(20);
                                trueGraph->SetMarkerSize(0.6);
				dataGraph->SetName("dataGraph");
				mcGraph->SetName("mcGraph");
				trueGraph->SetName("trueGraph");


				dataGraph->Draw("AP");			
				dataGraph->GetXaxis()->SetTitle(cutVariableName.c_str());
				dataGraph->GetYaxis()->SetTitle(yTitle.c_str());
				fileNameSave = Form("data");
				dataGraph->GetXaxis()->SetLimits(xMin,xMax);
				plotsRecording(directoryNameSave, fileNameSave, c1);
				
				c1->Clear();
	
				mcGraph->Draw("AP");
                                mcGraph->GetXaxis()->SetTitle(cutVariableName.c_str());
                                mcGraph->GetYaxis()->SetTitle(yTitle.c_str());
                                fileNameSave = Form("MC");
				mcGraph->GetXaxis()->SetLimits(xMin,xMax);
                                plotsRecording(directoryNameSave, fileNameSave, c1);
				
				c1->Clear();

                                trueGraph->Draw("AP");
                                trueGraph->GetXaxis()->SetTitle(cutVariableName.c_str());
                                trueGraph->GetYaxis()->SetTitle(yTitle.c_str());
                                fileNameSave = Form("true");
				trueGraph->GetXaxis()->SetLimits(xMin,xMax);
                                plotsRecording(directoryNameSave, fileNameSave, c1);
				
				c1->Clear();
			
				mg2->Add(dataGraph);
				mg2->Add(mcGraph);
				mg2->Draw("AP");
				mg2->GetXaxis()->SetTitle(cutVariableName.c_str());
                                mg2->GetYaxis()->SetTitle(yTitle.c_str());	
				fileNameSave = Form("dataMC");
				mg2->GetXaxis()->SetLimits(xMin,xMax);
        			leg2 = new TLegend(0.6,0.7,0.9,0.9,"","brNDC");                        
				leg2->SetTextSize(0.030);
				leg2->SetLineColor(kWhite);
				leg2->SetShadowColor(kWhite);
				leg2->AddEntry(dataGraph->GetName(),"data scale","lep");
				leg2->AddEntry(mcGraph->GetName(),"MC scale","lep");
				leg2->Draw();
				plotsRecording(directoryNameSave, fileNameSave, c1);	
					
				c1->Clear();

				mg3->Add(dataGraph);
                                mg3->Add(mcGraph);
				mg3->Add(trueGraph);
                                mg3->Draw("AP");
				mg3->GetXaxis()->SetTitle(cutVariableName.c_str());
                                mg3->GetYaxis()->SetTitle(yTitle.c_str());
				fileNameSave = Form("dataMCtrue");
				mg3->GetXaxis()->SetLimits(xMin,xMax);
                                plotsRecording(directoryNameSave, fileNameSave, c1);

                                c1->Clear();

				leg2->Delete();
				leg2 = 0;
				dataGraph->Delete();
				dataGraph = 0;
				mcGraph->Delete();
				mcGraph = 0;
				trueGraph->Delete();
				trueGraph = 0;
				/*mg2->Delete();
				mg2 = 0;
				mg3->Delete();
                                mg3 = 0;
				*/
			}
		}
	


	delete [] x;
	x = 0;
	delete [] xl;
        xl = 0;
	delete [] xr;
        xr = 0;
	delete [] y;
        y = 0;
	delete [] yl;
        yl = 0;
	delete [] yr;
        yr = 0;

	delete [] x_mc;
        x_mc = 0;
        delete [] xl_mc;
        xl_mc = 0;
        delete [] xr_mc;
        xr_mc = 0;
        delete [] y_mc;
        y_mc = 0;
        delete [] yl_mc;
        yl_mc = 0;
        delete [] yr_mc;
        yr_mc = 0;

	delete [] x_true;
        x_true = 0;
        delete [] xl_true;
        xl_true = 0;
        delete [] xr_true;
        xr_true = 0;
        delete [] y_true;
        y_true = 0;
        delete [] yl_true;
        yl_true = 0;
        delete [] yr_true;
        yr_true = 0;

	return 0;

}
