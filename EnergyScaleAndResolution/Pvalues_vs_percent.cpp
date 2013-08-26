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
	string recordingFile = Form("PvalueVSPercentage_%s_%s_%s_%sR9_%s",fitVariable.c_str(),dataType.c_str(),eta.c_str(),r9.c_str(),fitFunction.c_str());

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

	double yminPValue, ymaxPValue;
	double xminPValue, xmaxPValue;
	
	yminPValue = 0.0;
	ymaxPValue = 1.0;
	xminPValue = lowFitRange - 1;
        xmaxPValue = highFitRange + 1;

	double * pValueTab_temp = new double[(highFitRange - lowFitRange + 1) * nbFits];
	double * pValueTab = new double[highFitRange - lowFitRange + 1];
	//double * pValueTabError = new double[highFitRange - lowFitRange];
	double * xValueTab = new double[highFitRange - lowFitRange + 1];
	
	

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
			temp_iter = j * ( 2 + 6 + fitParameters[2] * 2 ) + fitParameters[2] * 2 + 5; 
			pValueTab_temp[(i-lowFitRange) * nbFits + j] = fitParameters[temp_iter];
			cout<<endl<<"(i-lowFitRange) * nbFits + j = "<<(i-lowFitRange) * nbFits + j<<", fitParameters[temp_iter] = "<<fitParameters[temp_iter];
		}

		xValueTab[i-lowFitRange] = i;

		fitParameters.erase(fitParameters.begin(),fitParameters.end());
	}

	TCanvas * c1 = new TCanvas("c1", "c1",0,0,600,600);	

	c1->SetTickx(1);
        c1->SetTicky(1);
        c1->SetGridx(1);
        c1->SetGridy(1);

	TGraph * pValueGraph = 0;

	for(int i = 0; i < nbFits; i++)
	{

		for(int j = 0; j < (highFitRange - lowFitRange); j++)
		{
			pValueTab[j] = pValueTab_temp[j*nbFits + i];
		}

		recordingFile = Form("PvalueVSPercentage_%s_%s_%s_%sR9_%s",fitVariable.c_str(),dataType.c_str(),eta.c_str(),r9.c_str(),fitFunction.c_str());
		recordingFile += Form("_%d",i+1);

		pValueGraph = new TGraph((highFitRange - lowFitRange) + 1,xValueTab, pValueTab);

		pValueGraph->Draw("ACP");

		pValueGraph->GetXaxis()->SetTitle("Fit ranges (%)");
		pValueGraph->GetYaxis()->SetTitle("p-values");	
		pValueGraph->GetXaxis()->SetLabelFont(42);
        	pValueGraph->GetXaxis()->SetTitleFont(42);
        	pValueGraph->GetXaxis()->SetLabelSize(0.03);
		pValueGraph->GetYaxis()->SetLabelFont(42);
        	pValueGraph->GetYaxis()->SetTitleOffset(1.24);
        	pValueGraph->GetYaxis()->SetTitleFont(42);
       	 	pValueGraph->GetYaxis()->SetLabelSize(0.03);	
	
		pValueGraph->SetMarkerColor(1);
        	pValueGraph->SetLineColor(4);
		pValueGraph->SetLineWidth(2);
        	pValueGraph->SetMarkerStyle(21);
        	pValueGraph->SetMarkerSize(0.8);
		
		pValueGraph->GetYaxis()->SetRangeUser(yminPValue,ymaxPValue);
        	pValueGraph->GetXaxis()->SetLimits(xminPValue,xmaxPValue);

		plotsRecording(recordingDirectory, recordingFile, c1);

		pValueGraph->Delete();
		pValueGraph = 0;
		
		c1->Clear();
	}

	delete [] pValueTab_temp;
	pValueTab_temp = 0;
	delete [] pValueTab;
        pValueTab = 0;
	delete [] xValueTab;
	xValueTab = 0;

	return 0;
}
