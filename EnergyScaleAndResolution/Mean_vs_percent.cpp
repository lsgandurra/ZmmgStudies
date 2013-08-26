//Script by Louis Sgandurra (March 2012)
//combine P-Value for different ranges of fit, in one graph.
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
                cerr << "arguments should be passed : directoryName, dataType, fitVariable, cutVariable, eta, r9, fitFunction, lowFitRange, highFitRange" <<endl; 
                return 1;

        }    

        string directoryName = "Results_v2";
        string dataType = "data";
        string fitVariable = "mmg_s";
       	string cutVariable = "Photon_Et";
	string eta = "Barrel"; 
        string r9 = "all";
        string fitFunction = "voigtian";
        int lowFitRange = 60; 
        int highFitRange = 100;

	if( argc > 1 ) directoryName = argv[1];
        if( argc > 2 ) dataType = argv[2];
        if( argc > 3 ) fitVariable = argv[3];
        if( argc > 4 ) cutVariable = argv[4];
	if( argc > 5 ) eta = argv[5];
        if( argc > 6 ) r9 = argv[6];
        if( argc > 7 ) fitFunction = argv[7];
        if( argc > 8 ) 
        {
                std::stringstream ss ( argv[8] );
                ss >> lowFitRange;
        }
        if( argc > 9 ) 
        {
                std::stringstream ss ( argv[9] );
                ss >> highFitRange;
        }


	gROOT->Reset();
	setTDRStyle();
	TGaxis::SetMaxDigits(3);

	string fileName = Form("%s/%s/%dPercents/%s_%s_%sR9_%s/fitsInformationsRaw.txt",directoryName.c_str(),dataType.c_str(),lowFitRange,fitVariable.c_str(),eta.c_str(),r9.c_str(),fitFunction.c_str());
	system(Form("mkdir -p %s/%s/CombinedGraphs/",directoryName.c_str(),dataType.c_str()));
	string recordingDirectory = Form("%s/%s/CombinedGraphs/",directoryName.c_str(),dataType.c_str());
	string recordingFile = Form("MuVSPercentage_%s_%s_%s_%sR9_%s",fitVariable.c_str(),dataType.c_str(),eta.c_str(),r9.c_str(),fitFunction.c_str());

	int nbFits = 0;
	int nbRows = rowsNumberInFile(fileName.c_str());	
	int nbFitParams = 0;

	ifstream tempFile(fileName.c_str());
	tempFile >> nbFitParams;
	tempFile >> nbFitParams;	
	tempFile >> nbFitParams;
	tempFile.close();

	nbFits = nbRows / (2 + 6 + 2 * nbFitParams);

	double temp_number = 0;
	int temp_iter = 0;
		
	vector <double> fitParameters;

	double yminMean, ymaxMean;
	double xminMean, xmaxMean;
	
	yminMean = -2.0;
	ymaxMean = 2.0;
	xminMean = lowFitRange - 1;
        xmaxMean = highFitRange + 1;

	double * meanTab_temp = new double[(highFitRange - lowFitRange + 1) * nbFits];
	double * meanTabError_temp = new double[(highFitRange - lowFitRange + 1) * nbFits];
	double * meanTab = new double[highFitRange - lowFitRange + 1];
	double * meanTabError = new double[highFitRange - lowFitRange + 1];
	double * xValueTab = new double[highFitRange - lowFitRange + 1];
	double * xValueTabError = new double[highFitRange - lowFitRange + 1];
	
	

	for(int i = lowFitRange; i <= highFitRange; i++)
	{
		fileName = Form("%s/%s/%dPercents/%s_%s_%sR9_%s/fitsInformationsRaw.txt",directoryName.c_str(),dataType.c_str(),i,fitVariable.c_str(),eta.c_str(),r9.c_str(),fitFunction.c_str());
		ifstream fitParametersFile(fileName.c_str());
		
		if(fitParametersFile)    
                {
                        cout << "OK : The file is open." << endl;
                }
                else
                {
                        cout << "ERROR: Impossible to open the file." << endl;
                }
	
		for(int j = 0; j < nbRows; j++)
                {
                        fitParametersFile >> temp_number;
                        fitParameters.push_back(temp_number);
                }

		fitParametersFile.close();

	
		for(int j = 0; j < nbFits; j++)
		{ 
			temp_iter = j * ( 2 + 6 + fitParameters[2] * 2 ) + 3;
			meanTab_temp[(i-lowFitRange) * nbFits + j] = fitParameters[temp_iter];
			meanTabError_temp[(i-lowFitRange) * nbFits + j] = fitParameters[temp_iter + 1];	
		}

		xValueTab[i-lowFitRange] = i;
		xValueTabError[i-lowFitRange] = 0;

		fitParameters.erase(fitParameters.begin(),fitParameters.end());
	}

	TCanvas * c1 = new TCanvas("c1", "c1",0,0,600,600);	

	c1->SetTickx(1);
        c1->SetTicky(1);
        c1->SetGridx(1);
        c1->SetGridy(1);

	double averageMean = 0;

	TGraphErrors * meanGraph = 0;

	for(int i = 0; i < nbFits; i++)
	{
		averageMean = 0;

		for(int j = 0; j <= (highFitRange - lowFitRange); j++)
		{
			meanTab[j] = meanTab_temp[j*nbFits + i] * 100;
			meanTabError[j] = meanTabError_temp[j*nbFits + i] * 100;
			averageMean += meanTab[j];
		}
	
		averageMean = averageMean / (highFitRange - lowFitRange + 1);

		recordingFile = Form("MuVSPercentage_%s_%s_%s_%sR9_%s",fitVariable.c_str(),dataType.c_str(),eta.c_str(),r9.c_str(),fitFunction.c_str());
		recordingFile += Form("_%d",i+1);

		meanGraph = new TGraphErrors((highFitRange - lowFitRange) + 1,xValueTab, meanTab, xValueTabError, meanTabError);

		meanGraph->Draw("ACP");

		meanGraph->GetXaxis()->SetTitle("Fit ranges (%)");
		meanGraph->GetYaxis()->SetTitle("#mu(E_{#gamma})");	
		meanGraph->GetXaxis()->SetLabelFont(42);
        	meanGraph->GetXaxis()->SetTitleFont(42);
        	meanGraph->GetXaxis()->SetLabelSize(0.03);
		meanGraph->GetYaxis()->SetLabelFont(42);
        	meanGraph->GetYaxis()->SetTitleOffset(1.24);
        	meanGraph->GetYaxis()->SetTitleFont(42);
       	 	meanGraph->GetYaxis()->SetLabelSize(0.03);	
	
		meanGraph->SetMarkerColor(1);
        	meanGraph->SetLineColor(4);
		meanGraph->SetLineWidth(2);
        	meanGraph->SetMarkerStyle(21);
        	meanGraph->SetMarkerSize(0.8);
	
		yminMean = averageMean - 3.5;
		ymaxMean = averageMean + 3.5;
	
		meanGraph->GetYaxis()->SetRangeUser(yminMean,ymaxMean);
        	meanGraph->GetXaxis()->SetLimits(xminMean,xmaxMean);

		plotsRecording(recordingDirectory, recordingFile, c1);

		meanGraph->Delete();
		meanGraph = 0;
		
		c1->Clear();
	}

	delete [] meanTab_temp;
	meanTab_temp = 0;
	delete [] meanTabError_temp;
        meanTabError_temp = 0;
	delete [] meanTab;
        meanTab = 0;
	delete [] meanTabError;
        meanTabError = 0;
	delete [] xValueTab;
	xValueTab = 0;
	delete [] xValueTabError;
        xValueTabError = 0;	

	return 0;
}
