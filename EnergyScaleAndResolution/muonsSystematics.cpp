#include "functions.h"
#include "setTDRStyle.C"
using namespace std;

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
	string directoryName = "Results_v6_RecoEnergy_MuSys";
        string dataType = "MC";
        string fitVariable = "mmg_s";
        string eta = "Endcaps"; 
        string r9 = "low";
        string fitFunction = "voigtian";

	if( argc > 1 ) directoryName = argv[1];
	if( argc > 2 ) dataType = argv[2];
        if( argc > 3 ) fitVariable = argv[3];	
	if( argc > 4 ) eta = argv[4];	
	if( argc > 5 ) r9 = argv[5];
	if( argc > 6 ) fitFunction = argv[6];


	int nBins = 150; //FIXME
	double xMin, xMax;
	xMin = -0.03; //FIXME
	xMax = 0.03; //FIXME
	
	string xVariableName, yVariableName;

	xVariableName = "Fitted Voigtian mean";
	yVariableName = "Events / 0.0001"; //FIXME

	string fileName = "fitsInformationsRaw.txt"; 

	string directoryNameOpen = directoryName + Form("/%s/toy_1/%s_%s_%sR9_%s/",dataType.c_str(),fitVariable.c_str(),eta.c_str(),r9.c_str(),fitFunction.c_str());

	vector <double> fitParameters;
        vector <double> reducedFitParameters;
        double temp_number = 0;
        int temp_iter = 0;
        int nbFits = 0;
        int nbRows = rowsNumberInFile(Form("%s%s",directoryNameOpen.c_str(),fileName.c_str()));
        int nbFitParams = 0;

        ifstream tempFile(Form("%s%s",directoryNameOpen.c_str(),fileName.c_str()));
        tempFile >> nbFitParams;
        tempFile >> nbFitParams;
        tempFile >> nbFitParams;
        tempFile.close();

        nbFits = nbRows / (3 + 5 + 2 * nbFitParams);

	
	for(int toy = 1; toy <= 100; toy++)
	{	
		directoryNameOpen = directoryName + Form("/%s/toy_%d/%s_%s_%sR9_%s/",dataType.c_str(),toy,fitVariable.c_str(),eta.c_str(),r9.c_str(),fitFunction.c_str());

		cout<<endl<<"directoryNameOpen = "<<directoryNameOpen;

        	ifstream fitParametersFile(Form("%s%s",directoryNameOpen.c_str(),fileName.c_str()));
		
		cout<<endl<<"fileName = "<<fileName;
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

		for(int i = 0; i < nbFits; i++)
        	{
			temp_iter = i * ( 2 + 6 + fitParameters[2] * 2 ) + 3;
			reducedFitParameters.push_back(fitParameters[temp_iter]);
		}


		fitParameters.erase(fitParameters.begin(),fitParameters.end());
		fitParametersFile.close();
	}

	gROOT->Reset();
        TGaxis::SetMaxDigits(5);
        setTDRStyle();

	TCanvas *c1 = new TCanvas("c1", "c1",0,0,600,600);

	TH1D * h1 = 0;
	TF1 * f1 = 0;  

	directoryName += Form("/%s/MuonSystematics/",dataType.c_str()); 

	double muonSystematics = 0;
	double mean = 0;

	for(int i = 0; i < nbFits; i++)
	{
		xMin = -0.06;  //FIXME
        	xMax = 0.06;  //FIXME
		fileName = Form("%s_%s_%sR9_%s_%d",fitVariable.c_str(),eta.c_str(),r9.c_str(),fitFunction.c_str(),i+1);;		

		h1 = new TH1D("h1","h1", nBins,xMin,xMax);

		for(int j = 1; j <= 100; j++)
		{
			h1->Fill(reducedFitParameters[ (j - 1) * nbFits + i ]);		
		}
		h1->Draw("");

		// --- Fit --- //

        	f1 = new TF1("f1","gaus",xMin,xMax);
        	h1->Fit(f1);
        	f1->Draw("SAMES");
		
		mean = f1->GetParameter(1);
		cout<<endl<<"mean = "<<mean;
		//for(double x = -0.03; x < 0.03; x += 0.0008)
		for(double x = -0.06; x < 0.06; x += 0.0004) //FIXME
		{

			if((mean > x) && (mean < ( x + 0.0004)))
			{
				xMin = x - 0.005 + 0.0004;
				xMax = x + 0.005 + 0.0004;
			}
		}
/*
		h1->Delete();
                h1 = 0;
		f1->Delete();
                f1 = 0;
	
		c1->Clear();

		h1 = new TH1D("h1","h1", nBins,xMin,xMax);
		for(int j = 1; j <= 100; j++)
                {
                        h1->Fill(reducedFitParameters[ (j - 1) * nbFits + i ]);    
                }
		f1 = new TF1("f1","gaus",xMin,xMax);
		h1->Draw("");
		h1->Fit(f1);
*/
		c1->Clear();
  
                h1->Draw("E1");
                f1->Draw("SAMES");
		
		h1->GetXaxis()->SetRangeUser(xMin,xMax);
		h1->SetMarkerStyle(20);
		h1->SetMarkerSize(0.8);
		f1->SetLineColor(kBlue);
		f1->SetLineWidth(3);
		h1->GetYaxis()->SetTitle(yVariableName.c_str()); 
                h1->GetXaxis()->SetTitle(xVariableName.c_str());
		h1->GetXaxis()->SetNdivisions(509);

		cout<<endl<<"xMin = "<<xMin<<", xMax = "<<xMax;
		muonSystematics = f1->GetParameter(2); //FIXME
		cout<<endl<<"muonSystematics = "<<muonSystematics;


        	plotsRecording(directoryName, fileName, c1);

        	c1->Clear();

		f1->Delete();
		f1 = 0;

		h1->Delete();
		h1 = 0;

		ofstream summaryFile(Form("%sSummary_MuonSystematics.txt",directoryName.c_str()), ios::app);

                summaryFile << fitVariable << " " << dataType << " >> " << r9 << " r9 " << eta << ", fitFunction = " << fitFunction <<", systematics : " << muonSystematics * 100 << " %" << endl;;

                summaryFile.close();
		

	
	}

	delete c1;
	c1 = 0;	

	
	return 0;

}






 
