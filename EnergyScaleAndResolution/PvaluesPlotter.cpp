//Script by Louis Sgandurra (January 2012)
//combine P-Value vs Ptgamma graphs, for different ranges of fit, in one graph.
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

        string directoryName = "Results_v1";
        string dataType = "data";
        string fitVariable = "mmg_s";
       	string cutVariable = "Photon_Et";
	string eta = "Barrel"; 
        string r9 = "all";
        string fitFunction = "voigtian";
        int lowFitRange = 70; 
        int highFitRange = 80;

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
	string recordingFile = Form("PvalueVS%s_%s_%s_%s_%sR9_%s_%d_%d",cutVariable.c_str(),fitVariable.c_str(),dataType.c_str(),eta.c_str(),r9.c_str(),fitFunction.c_str(),lowFitRange,highFitRange-1);

	int nbFits = 0;
	int nbRows = rowsNumberInFile(fileName.c_str());	
	int nbFitParams = 0;

	ifstream tempFile(fileName.c_str());
	tempFile >> nbFitParams;
	tempFile >> nbFitParams;	
	tempFile >> nbFitParams;
	tempFile.close();

	nbFits = nbRows / (2 + 6 + 2 * nbFitParams);

	double * xValueTab = new double[nbFits];
	double * xErrorLTab = new double[nbFits];
	double * xErrorRTab = new double[nbFits];

	double * pValueTab1 = new double[nbFits];
	double * pValueTab2 = new double[nbFits];
	double * pValueTab3 = new double[nbFits];
	double * pValueTab4 = new double[nbFits];
	double * pValueTab5 = new double[nbFits];
	double * pValueTab6 = new double[nbFits];
        double * pValueTab7 = new double[nbFits];
        double * pValueTab8 = new double[nbFits];
        double * pValueTab9 = new double[nbFits];
        double * pValueTab10 = new double[nbFits];
	double * pValueTab11 = new double[nbFits];

	double meanPvalueTab[11] = {0};
	double * PValueErrorTab = new double[nbFits];

	double temp_number = 0;
	int temp_iter = 0;
		
	vector <double> fitParameters;

	double yminPValue, ymaxPValue;
	double xminPValue, xmaxPValue;
	
	yminPValue = 0.0;
	ymaxPValue = 1.0;
	xminPValue = 0.0;
        xmaxPValue = 3.0;

	if(cutVariable == "Photon_SC_Eta") xminPValue = -3.0;
        if(cutVariable == "Photon_SC_rawEt") xmaxPValue = 100.0;
        if(cutVariable == "Photon_Et") xmaxPValue = 100.0;
        if(cutVariable == "Photon_E") xmaxPValue = 200.0;
        if(cutVariable == "Photon_SC_Eta") xmaxPValue = 3.0;
        if(cutVariable == "Photon_SC_brem") xmaxPValue = 15.0;


	for(int i = lowFitRange; i < highFitRange; i++)
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
			if(i == lowFitRange) pValueTab1[j] = fitParameters[temp_iter];
			if(i == (lowFitRange+1)) pValueTab2[j] = fitParameters[temp_iter];
			if(i == (lowFitRange+2)) pValueTab3[j] = fitParameters[temp_iter];
			if(i == (lowFitRange+3)) pValueTab4[j] = fitParameters[temp_iter];
			if(i == (lowFitRange+4)) pValueTab5[j] = fitParameters[temp_iter];
			if(i == (lowFitRange+5)) pValueTab6[j] = fitParameters[temp_iter];
			if(i == (lowFitRange+6)) pValueTab7[j] = fitParameters[temp_iter];
			if(i == (lowFitRange+7)) pValueTab8[j] = fitParameters[temp_iter];
			if(i == (lowFitRange+8)) pValueTab9[j] = fitParameters[temp_iter];
			if(i == (lowFitRange+9)) pValueTab10[j] = fitParameters[temp_iter];
			if(highFitRange == 101 && i == (lowFitRange+10)) pValueTab11[j] = fitParameters[temp_iter];
		}

		if(i == lowFitRange) 
		{
			for(int j = 0; j < nbFits; j++) 
			{
				temp_iter = j * ( 2 + 6 + fitParameters[2] * 2 ) + fitParameters[2] * 2 + 7; 
				xValueTab[j] = fitParameters[temp_iter];
			}
		}

		fitParameters.erase(fitParameters.begin(),fitParameters.end());

	}

	for(int i = 0; i < nbFits; i++)
	{
		meanPvalueTab[0] += pValueTab1[i];
		meanPvalueTab[1] += pValueTab2[i];
                meanPvalueTab[2] += pValueTab3[i];
                meanPvalueTab[3] += pValueTab4[i];
                meanPvalueTab[4] += pValueTab5[i];
                meanPvalueTab[5] += pValueTab6[i];
		meanPvalueTab[6] += pValueTab7[i];
                meanPvalueTab[7] += pValueTab8[i];
                meanPvalueTab[8] += pValueTab9[i];
                meanPvalueTab[9] += pValueTab10[i];
		if(highFitRange == 101) meanPvalueTab[10] += pValueTab11[i];

	}

	for(int i = 0; i < 11; i++)
        {
		meanPvalueTab[i] *= 1.0 / nbFits;
	}

	TCanvas * c1 = new TCanvas("c1", "c1",0,0,600,600);

	TMultiGraph * mg = new TMultiGraph();

	TGraph * pValueGraph1 = new TGraph(nbFits,xValueTab, pValueTab1);
	TGraph * pValueGraph2 = new TGraph(nbFits,xValueTab, pValueTab2);
	TGraph * pValueGraph3 = new TGraph(nbFits,xValueTab, pValueTab3);
	TGraph * pValueGraph4 = new TGraph(nbFits,xValueTab, pValueTab4);
	TGraph * pValueGraph5 = new TGraph(nbFits,xValueTab, pValueTab5);
	TGraph * pValueGraph6 = new TGraph(nbFits,xValueTab, pValueTab6);
	TGraph * pValueGraph7 = new TGraph(nbFits,xValueTab, pValueTab7);
	TGraph * pValueGraph8 = new TGraph(nbFits,xValueTab, pValueTab8);
	TGraph * pValueGraph9 = new TGraph(nbFits,xValueTab, pValueTab9);
	TGraph * pValueGraph10 = new TGraph(nbFits,xValueTab, pValueTab10);
	TGraph * pValueGraph11 = new TGraph(nbFits,xValueTab, pValueTab11);

	mg->Add(pValueGraph1);
	mg->Add(pValueGraph2);
	mg->Add(pValueGraph3);
	mg->Add(pValueGraph4);
	mg->Add(pValueGraph5);
	mg->Add(pValueGraph6);
	mg->Add(pValueGraph7);
	mg->Add(pValueGraph8);
	mg->Add(pValueGraph9);
	mg->Add(pValueGraph10);
	if(highFitRange == 101) mg->Add(pValueGraph11);
	
	mg->Draw("AP");

        if(cutVariable == "Photon_SC_rawEt") mg->GetXaxis()->SetTitle("P_{T RAW}^{#gamma}");
        if(cutVariable == "Photon_Et") mg->GetXaxis()->SetTitle("P_{T}^{#gamma}");
        if(cutVariable == "Photon_E") mg->GetXaxis()->SetTitle("E^{#gamma}");
        if(cutVariable == "Photon_SC_Eta") mg->GetXaxis()->SetTitle("#eta");
        if(cutVariable == "Photon_SC_brem") mg->GetXaxis()->SetTitle("#sigma_{#phi}/#sigma_{#eta}");
	
	mg->GetXaxis()->SetLabelFont(42);
        mg->GetXaxis()->SetTitleFont(42);
        mg->GetXaxis()->SetLabelSize(0.03);
        mg->GetYaxis()->SetTitle("P-Value");
        mg->GetYaxis()->SetLabelFont(42);
        mg->GetYaxis()->SetTitleOffset(1.24);
        mg->GetYaxis()->SetTitleFont(42);
        mg->GetYaxis()->SetLabelSize(0.03);
	//mg->SetTitle("");
        //mg->SetMarkerColor(4);
        //mg->SetMarkerStyle(21);
        //mg->SetMarkerSize(0.6);

	pValueGraph1->SetLineColor(629);
	pValueGraph1->SetMarkerColor(629);
	pValueGraph1->SetMarkerStyle(20);
	pValueGraph1->SetMarkerSize(1.15);
	pValueGraph1->SetName("pValueGraph1");
        pValueGraph2->SetLineColor(629);
        pValueGraph2->SetMarkerColor(629);
	pValueGraph2->SetMarkerStyle(24);
	pValueGraph2->SetMarkerSize(1.15);
	pValueGraph2->SetName("pValueGraph2");
	pValueGraph3->SetLineColor(397);
	pValueGraph3->SetMarkerColor(397);
	pValueGraph3->SetMarkerStyle(21);
	pValueGraph3->SetMarkerSize(1.15);
	pValueGraph3->SetName("pValueGraph3");
        pValueGraph4->SetLineColor(397);
	pValueGraph4->SetMarkerColor(397);
	pValueGraph4->SetMarkerStyle(25);
	pValueGraph4->SetMarkerSize(1.15);
	pValueGraph4->SetName("pValueGraph4");
        pValueGraph5->SetLineColor(413);
	pValueGraph5->SetMarkerColor(413);
	pValueGraph5->SetMarkerStyle(22);
	pValueGraph5->SetMarkerSize(1.15);
	pValueGraph5->SetName("pValueGraph5");
	pValueGraph6->SetLineColor(413);
        pValueGraph6->SetMarkerColor(413);
        pValueGraph6->SetMarkerStyle(26);
        pValueGraph6->SetMarkerSize(1.15);
        pValueGraph6->SetName("pValueGraph6");
	pValueGraph7->SetLineColor(597);
	pValueGraph7->SetMarkerColor(597);
	pValueGraph7->SetMarkerStyle(23);
	pValueGraph7->SetMarkerSize(1.15);
	pValueGraph7->SetName("pValueGraph7");
	pValueGraph8->SetLineColor(597);
	pValueGraph8->SetMarkerColor(597);
	pValueGraph8->SetMarkerStyle(32);
	pValueGraph8->SetMarkerSize(1.15);
	pValueGraph8->SetName("pValueGraph8");
	pValueGraph9->SetLineColor(613);
	pValueGraph9->SetMarkerColor(613);
	pValueGraph9->SetMarkerStyle(34);
	pValueGraph9->SetMarkerSize(1.15);
	pValueGraph9->SetName("pValueGraph9");
	pValueGraph10->SetLineColor(613);
	pValueGraph10->SetMarkerColor(613);
	pValueGraph10->SetMarkerStyle(28);
	pValueGraph10->SetMarkerSize(1.15);
	pValueGraph10->SetName("pValueGraph10");
	if(highFitRange == 101)
	{
		pValueGraph11->SetLineColor(1);
		pValueGraph11->SetMarkerColor(1);
		pValueGraph11->SetMarkerStyle(29);
		pValueGraph11->SetMarkerSize(1.15);
		pValueGraph11->SetName("pValueGraph11");
	}


	TLegend leg(0.65,0.3,0.9,0.9,"","brNDC");
        leg.SetTextSize(0.023);
        leg.SetFillColor(kWhite);
        leg.SetLineColor(kWhite);
        leg.SetShadowColor(kWhite);
	string temp = "";
	temp = Form("%d",lowFitRange); temp += "%, m = "; temp += Form("%f",meanPvalueTab[0]);
        leg.AddEntry(pValueGraph1->GetName(),temp.c_str(),"p");
	temp = Form("%d",lowFitRange+1); temp += "%, m = "; temp += Form("%f",meanPvalueTab[1]);
        leg.AddEntry(pValueGraph2->GetName(),temp.c_str(),"p");
	temp = Form("%d",lowFitRange+2); temp += "%, m = "; temp += Form("%f",meanPvalueTab[2]);
        leg.AddEntry(pValueGraph3->GetName(),temp.c_str(),"p");
	temp = Form("%d",lowFitRange+3); temp += "%, m = "; temp += Form("%f",meanPvalueTab[3]);
        leg.AddEntry(pValueGraph4->GetName(),temp.c_str(),"p");
	temp = Form("%d",lowFitRange+4); temp += "%, m = "; temp += Form("%f",meanPvalueTab[4]);
        leg.AddEntry(pValueGraph5->GetName(),temp.c_str(),"p");
	temp = Form("%d",lowFitRange+5); temp += "%, m = "; temp += Form("%f",meanPvalueTab[5]);
        leg.AddEntry(pValueGraph6->GetName(),temp.c_str(),"p");	
	temp = Form("%d",lowFitRange+6); temp += "%, m = "; temp += Form("%f",meanPvalueTab[6]);
        leg.AddEntry(pValueGraph7->GetName(),temp.c_str(),"p");
	temp = Form("%d",lowFitRange+7); temp += "%, m = "; temp += Form("%f",meanPvalueTab[7]);
        leg.AddEntry(pValueGraph8->GetName(),temp.c_str(),"p");
	temp = Form("%d",lowFitRange+8); temp += "%, m = "; temp += Form("%f",meanPvalueTab[8]);
        leg.AddEntry(pValueGraph9->GetName(),temp.c_str(),"p");
	temp = Form("%d",lowFitRange+9); temp += "%, m = "; temp += Form("%f",meanPvalueTab[9]);
        leg.AddEntry(pValueGraph10->GetName(),temp.c_str(),"p");
	if(highFitRange == 101)
        {
        	temp = Form("%d",lowFitRange+10); temp += "%, m = "; temp += Form("%f",meanPvalueTab[10]);
                leg.AddEntry(pValueGraph11->GetName(),temp.c_str(),"p");
        }
	leg.Draw();


	TLatex latexLabel;
        latexLabel.SetNDC();
        latexLabel.SetTextAlign(11);
        latexLabel.SetTextFont(42);
	latexLabel.SetTextSize(0.030);
        latexLabel.SetNDC();
        latexLabel.DrawLatex(0.13, 0.96, "CMS Preliminary 2012, #sqrt{s} = 8 TeV");     

	//if(dataType == "MC") latexLabel.DrawLatex(0.17, 0.88, "Simulation");
        //if(dataType == "data") latexLabel.DrawLatex(0.17, 0.88, Form("Data, #int L = %s fb^{-1}",lumi.c_str()));
        //latexLabel.DrawLatex(0.17, 0.83,Form("ECAL %s",eta.c_str()));   
        //if(r9 == "low") latexLabel.DrawLatex(0.17, 0.78,"E_{T}^{#gamma} > 25 GeV, r9 < 0,94");
        //if(r9 == "high") latexLabel.DrawLatex(0.17, 0.78,"E_{T}^{#gamma} > 25 GeV, r9 < 0,95");
        //if(r9 == "all") latexLabel.DrawLatex(0.17, 0.78,"E_{T}^{#gamma} > 25 GeV, All r9");
 
	gStyle->SetPadBorderMode(0);

        mg->GetYaxis()->SetRangeUser(yminPValue,ymaxPValue);
        mg->GetXaxis()->SetLimits(xminPValue,xmaxPValue);

        c1->SetTickx(1);
        c1->SetTicky(1);
        c1->SetGridx(1);
        c1->SetGridy(1);
        c1->Modified();
        c1->cd();
        c1->SetSelected(c1);
        c1->ToggleToolBar();

	plotsRecording(recordingDirectory, recordingFile, c1);

	//mg->Delete();
	//mg = 0;
	pValueGraph1->Delete();
	pValueGraph1 = 0;
	pValueGraph2->Delete();
        pValueGraph2 = 0;
	pValueGraph3->Delete();
        pValueGraph3 = 0;
	pValueGraph4->Delete();
        pValueGraph4 = 0;
	pValueGraph5->Delete();
        pValueGraph5 = 0;
	pValueGraph6->Delete();
        pValueGraph6 = 0;
	pValueGraph7->Delete();
        pValueGraph7 = 0;
	pValueGraph8->Delete();
        pValueGraph8 = 0;
	pValueGraph9->Delete();
        pValueGraph9 = 0;
	pValueGraph10->Delete();
        pValueGraph10 = 0;
	pValueGraph11->Delete();
        pValueGraph11 = 0;

	delete [] xValueTab;
	xValueTab = 0;
	delete [] xErrorLTab;
	xErrorLTab = 0;
	delete [] xErrorRTab;
	xErrorRTab = 0;
	delete [] pValueTab1;
	pValueTab1 = 0;	
	delete [] pValueTab2;
	pValueTab2 = 0;
	delete [] pValueTab3;
	pValueTab3 = 0;
	delete [] pValueTab4;
	pValueTab4 = 0;
	delete [] pValueTab5;
	pValueTab5 = 0;	
	delete [] pValueTab6;
	pValueTab6 = 0;
	delete [] pValueTab7;
	pValueTab7 = 0;
	delete [] pValueTab8;
	pValueTab8 = 0;
	delete [] pValueTab9;
	pValueTab9 = 0;
	delete [] pValueTab10;
	pValueTab10 = 0;
	delete [] pValueTab11;
	pValueTab11 = 0;
	delete [] PValueErrorTab;
	PValueErrorTab = 0;

	return 0;
}
