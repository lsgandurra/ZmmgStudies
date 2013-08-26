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
                cerr << "arguments should be passed : directoryName, dataType, fitVariable, eta, r9, fitFunction" <<endl; 
                return 1;

        }    

        //string directoryName = "Results_v2/InjectedResolution_1Percent";
        string directoryName = "Results_v4";
	string dataType = "data";
        string fitVariable = "mmg_s";
	string eta = "Barrel"; 
        string r9 = "all";
        string fitFunction = "voigtian";

	if( argc > 1 ) directoryName = argv[1];
        if( argc > 2 ) dataType = argv[2];
        if( argc > 3 ) fitVariable = argv[3];
	if( argc > 5 ) eta = argv[4];
        if( argc > 6 ) r9 = argv[5];
        if( argc > 7 ) fitFunction = argv[6];

	gROOT->Reset();
	setTDRStyle();
	TGaxis::SetMaxDigits(3);

	string fileName = Form("%s/InjectedResolution_0Percent/%s/Selected_Fits/%s_%s_%sR9_%s/fitsInformationsRaw.txt",directoryName.c_str(),dataType.c_str(),fitVariable.c_str(),eta.c_str(),r9.c_str(),fitFunction.c_str());
	system(Form("mkdir -p %s/CombinedGraphs/",directoryName.c_str()));
	string recordingDirectory = Form("%s/CombinedGraphs/",directoryName.c_str());
	string recordingFile = Form("MuVSInjectedResolution_%s_%s_%s_%sR9_%s",fitVariable.c_str(),dataType.c_str(),eta.c_str(),r9.c_str(),fitFunction.c_str());

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
	xminMean = -0.5;
        xmaxMean = 4.5;

	double upperLim = 2.5;
        if(eta == "Barrel") upperLim = 1.5;
	if(eta == "Endcaps") upperLim = 2.5;
	if(eta == "Barrel" && r9 == "high" && dataType == "data") upperLim = 1.5;
        if(eta == "Barrel" && r9 == "all" && dataType == "data") upperLim = 1;
        if(eta == "Endcaps" && r9 == "all" && dataType == "data") upperLim = 2.5;
        if(eta == "Barrel" && r9 == "high" && dataType == "MC") upperLim = 1.5;
        if(eta == "Barrel" && r9 == "all" && dataType == "MC") upperLim = 1;
        if(eta == "Endcaps" && r9 == "all" && dataType == "MC") upperLim = 2;	

	int points = (upperLim / 0.5) + 1;

	double * meanTab_temp = new double[points * nbFits];
	double * meanTabError_temp = new double[points * nbFits];
	double * meanTab = new double[points];
	double * meanTabError = new double[points];
	double * xValueTab = new double[points];
	double * xValueTabError = new double[points];

	int iter1 = 0;

	//for(int i = 0; i <= 5; i++)
	for(double i = 0; i <= upperLim; i+= 0.5)
	{
		cout<<endl<<"i = "<<i<<endl;
		std::ostringstream iString;
        	if( (iter1 + 2) % 2 == 0) iString << setprecision (0) << fixed << i;
		else iString << setprecision (1) << fixed << i;
		string iText = iString.str();

		fileName = Form("%s/InjectedResolution_%sPercent/%s/Selected_Fits/%s_%s_%sR9_%s/fitsInformationsRaw.txt",directoryName.c_str(),iText.c_str(),dataType.c_str(),fitVariable.c_str(),eta.c_str(),r9.c_str(),fitFunction.c_str());
		cout<<endl<<"fileName = "<<fileName<<endl;
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
			meanTab_temp[iter1 * nbFits + j] = fitParameters[temp_iter];
			meanTabError_temp[iter1 * nbFits + j] = fitParameters[temp_iter + 1];	
		}

		xValueTab[iter1] = i;
		xValueTabError[iter1] = 0.25;

		fitParameters.erase(fitParameters.begin(),fitParameters.end());
	
		iter1 += 1;
	}

	iter1 = 0;

	TCanvas * c1 = new TCanvas("c1", "c1",0,0,600,600);	

	c1->SetTickx(1);
        c1->SetTicky(1);
        c1->SetGridx(1);
        c1->SetGridy(1);

	double averageMean = 0;
	double averageMeanError = 0;

	TGraphErrors * meanGraph = 0;

	double systematic = 0;

	for(int i = 0; i < nbFits; i++)
	{
		systematic = 0;
		averageMean = 0;

		//for(int j = 0; j <= 5; j++)
		for(double j = 0; j <= upperLim; j+= 0.5)
		{
			meanTab[iter1] = meanTab_temp[iter1*nbFits + i] * 100;
			meanTabError[iter1] = meanTabError_temp[iter1*nbFits + i] * 100;
			averageMean += meanTab[iter1];
			averageMeanError += meanTabError[iter1];
			cout<<endl<<"meanTab[iter1] = "<<meanTab[iter1]<<endl;
			cout<<endl<<"meanTabError[iter1] = "<<meanTabError[iter1]<<endl;
			if(abs(meanTab[0] - meanTab[iter1]) > systematic ) systematic = abs(meanTab[0] - meanTab[iter1]);
			cout<<endl<<"systematic = "<<systematic;
			iter1 += 1; 
		}
	
		averageMean = averageMean / points;
		averageMeanError = averageMeanError / points;

		recordingFile = Form("MuVSInjectedResolution_%s_%s_%s_%sR9_%s",fitVariable.c_str(),dataType.c_str(),eta.c_str(),r9.c_str(),fitFunction.c_str());
		recordingFile += Form("_%d",i+1);

		meanGraph = new TGraphErrors(points,xValueTab, meanTab, xValueTabError, meanTabError);

		meanGraph->Draw("AP");

		meanGraph->GetXaxis()->SetTitle("Injected extra resolution (%)");
		meanGraph->GetYaxis()->SetTitle("#mu(E_{#gamma}) (%)");	
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
	
		//yminMean = averageMean - 1.0;
		//ymaxMean = averageMean + 1.0;
		yminMean = averageMean - 2 * averageMeanError;
		ymaxMean = averageMean + 2 * averageMeanError; 	

		double tempNumber = 0.2;
		if(eta == "Barrel") tempNumber = 0.125;
		if(eta == "Endcaps") tempNumber = 0.225;

		yminMean = averageMean - 2 * tempNumber;
                ymaxMean = averageMean + 2 * tempNumber;	

		meanGraph->GetYaxis()->SetRangeUser(yminMean,ymaxMean);
        	meanGraph->GetXaxis()->SetLimits(xminMean,xmaxMean);

		plotsRecording(recordingDirectory, recordingFile, c1);


		meanGraph->Delete();
		meanGraph = 0;
		
		c1->Clear();
		
		ofstream summaryFile(Form("%sSummary_%d.txt",recordingDirectory.c_str(), i+1), ios::app);
        	summaryFile << fitVariable << " " << dataType << " >> " << r9 << " r9 " << eta << ", " << fitFunction << ", systematic (resolution) = " << systematic << " (%)"<< endl;

        	summaryFile.close();


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
