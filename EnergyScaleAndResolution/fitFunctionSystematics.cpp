//#include "fitFunctions.h"
#include "RooCruijff.hh"
#include "functions.h"
#include "setTDRStyle.C"

int main(int argc, char *argv[])
{

        // --- Initialization --- //

        for(int iarg = 0 ; iarg < argc; iarg++)
        {
                cout << "argv[" << iarg << "]= " << argv[iarg] << endl;
        }

        if( argc == 1 ) 
        {
                cerr << "arguments should be passed : directoryName, dataType, fitVariable, eta, r9, trueModel, testModel" <<endl; 
                return 1;

        }    

        string directoryName = "Results_v3";
        string dataType = "data";
        string fitVariable = "mmg_s";
        string eta = "Barrel"; 
        string r9 = "low";
        string trueModel = "voigtian";
	string testModel = "cruijff"; 

        if( argc > 1 ) directoryName = argv[1];
        if( argc > 2 ) dataType = argv[2];
        if( argc > 3 ) fitVariable = argv[3];
        if( argc > 4 ) eta = argv[4];
        if( argc > 5 ) r9 = argv[5];
	if( argc > 6 ) trueModel = argv[6];
	if( argc > 7 ) testModel = argv[7];


	string directoryNameOpen = directoryName + Form("/%s/Selected_Fits/%s_%s_%sR9_%s/",dataType.c_str(),fitVariable.c_str(),eta.c_str(),r9.c_str(),trueModel.c_str());
	string fileName = "fitsInformationsRaw.txt";

	cout<<endl<<"directoryNameOpen = "<<directoryNameOpen<<endl;
	cout<<endl<<"fileName = "<<fileName<<endl;

	double xMinFitVariable, xMaxFitVariable, xMin, xMax;
	double nBins = 40; 
        string fitVariableName = "";    
	if(fitVariable == "mmg_s")
        {
		fitVariableName = "s";
                xMinFitVariable = -0.5;
                xMaxFitVariable = 0.5;
		xMin = -0.06;
		xMax = 0.06;
		nBins = 60;
        }
        if(fitVariable == "mmg_s_true")
        {	
		fitVariableName = "s_{TRUE}";
                xMinFitVariable = -0.2;
                xMaxFitVariable = 0.2;
        	xMin = -0.04;
                xMax = 0.04;	
	}
        if(fitVariable == "Mmumugamma")
        {
		fitVariableName = "M_{#mu#mu#gamma} (GeV)";
                xMinFitVariable = 70.0;
                xMaxFitVariable = 110.0;
		xMin = 70;
                xMax = 110;
        }

	double rangeMin, rangeMax;
	rangeMin = rangeMax = 0;
	double fitFunctionSystematics = 0;
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

	
	ifstream fitParametersFile(Form("%s%s",directoryNameOpen.c_str(),fileName.c_str()));
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

	directoryName += Form("/%s/Selected_Fits/FitFunctionSystematics/",dataType.c_str()); //FIXME
	fileName = Form("%s_%s_%sR9_%s",fitVariable.c_str(),eta.c_str(),r9.c_str(),trueModel.c_str());

	RooRealVar variable(fitVariable.c_str(), fitVariableName.c_str(), xMinFitVariable, xMaxFitVariable);
	RooAbsData * absData = 0;
	RooPlot* frame = 0;
	double chi2 = 0;
		

	gROOT->Reset();
        TGaxis::SetMaxDigits(4);
        setTDRStyle();
	TCanvas * c1 = new TCanvas("c1", "c1",0,0,600,600);
	TLatex latexLabel;
	for(int i = 0; i < nbFits; i++)
	{
		fileName += Form("_%d",i);
		if(trueModel == "voigtian" && testModel == "cruijff")
		{
			temp_iter = i * ( 2 + 6 + fitParameters[2] * 2 ) + 3;	
			reducedFitParameters.push_back(fitParameters[temp_iter - 2]); //entries	
			reducedFitParameters.push_back(fitParameters[temp_iter]); //mean
			reducedFitParameters.push_back(fitParameters[temp_iter + 2]); //sigma
			reducedFitParameters.push_back(fitParameters[temp_iter + 4]); //width	
	

			cout<<endl<<"reducedFitParameters[0] = "<<reducedFitParameters[0];
			cout<<endl<<"reducedFitParameters[1] = "<<reducedFitParameters[1];
			cout<<endl<<"reducedFitParameters[2] = "<<reducedFitParameters[2];
			cout<<endl<<"reducedFitParameters[3] = "<<reducedFitParameters[3];
		

			// --- Voigtian --- //
						
			RooRealVar mean("mean","mean",reducedFitParameters[1],"GeV");
        		//mean.setConstant(kTRUE);
        		RooRealVar sigma("sigma","sigma",reducedFitParameters[2],"GeV");
        		//sigma.setConstant(kTRUE);
        		RooRealVar width("width","width",reducedFitParameters[3],"GeV");
        		//width.setConstant(kTRUE);
			
			RooVoigtian genPdf("genPdf","Voigtian",variable,mean,sigma,width);		

			// --- Cruijff --- //

			//RooRealVar m0("m0","Fitted Cruijff mean",0.0,-0.1,0.1);
        		//RooRealVar m0("m0","Fitted Cruijff mean",0.0,-0.04,0.04);
			RooRealVar m0("m0","Fitted Cruijff mean",0.0,xMin,xMax);
			RooRealVar sigmaL("sigmaL","sigmaL",0.5,0.0,1.0);
        		RooRealVar sigmaR("sigmaR","sigmaR",0.5,0.0,1.0);
        		RooRealVar alphaL("alphaL","alphaL",0.5,0.0,1.0);
        		RooRealVar alphaR("alphaR","alphaR",0.5,0.0,1.0);
				
			RooCruijff fitPdf("fitPdf","Cruijff",variable,m0,sigmaL,sigmaR,alphaL,alphaR);


			//RooRealVar meanG("meanG","mean of gaussian",0.0,-0.1,0.1) ;
			//RooRealVar sigmaG("sigmaG","width of gaussian",0.01,0.0,0.1) ;
			//RooGaussian gauss("gauss","gaussian PDF",m0,meanG,sigmaG) ;


			//RooMCStudy mcs(gauss,gauss,x,"","mhv");
			RooMCStudy mcs (genPdf,variable,FitModel(fitPdf),Silence()) ;
			mcs.generateAndFit(1000,(int) reducedFitParameters[0]); 

			//mcs.fitParDataSet().Print("v");
			RooDataSet data = mcs.fitParDataSet();
			absData = data.reduce(RooArgSet(m0));
			//absData->Print("v");	


			RooRealVar meanG("meanG","mean of gaussian",0.0,xMin,xMax) ;
                        RooRealVar sigmaG("sigmaG","width of gaussian",0.01,0.0,0.1) ;
                        RooGaussian gauss("gauss","gaussian PDF",m0,meanG,sigmaG) ;


			/*
			frame = mcs.plotParam(m0,Name("myhist")) ;
			frame->Draw();
			
			hist = (RooDataHist*) frame->findObject("myhist", RooHist::Class());
			frame->Clear();	

			cout << endl << "Th1 creation "<<endl;	

			th1 = (TH1*) hist->createHistogram("th1 test",m0,50);			

			cout << endl << "bdata creation "<<endl;

			RooDataHist bdata("bdata","bdata",m0,th1);	
			*/
			RooBinning b;
			b.setRange(xMin,xMax);
			b.addUniform(nBins,xMin,xMax);
			frame = mcs.plotParam(m0,Name("myhist"),Binning(b),Range(xMin,xMax));
			//frame = mcs.plotParam(m0,Name("myhist"));
			//absData->plotOn(frame,Binning(b));
			
			gauss.fitTo(*absData);
			if(dataType == "data") 
			{	
				rangeMin = meanG.getVal() - 3 * sigmaG.getVal();
				rangeMax = meanG.getVal() + 3 * sigmaG.getVal();	
			}
			if(dataType == "MC") 
                        {           
                                rangeMin = meanG.getVal() - 2.5 * sigmaG.getVal();
                                rangeMax = meanG.getVal() + 2.5 * sigmaG.getVal();    
                        }	
			gauss.fitTo(*absData,Range(rangeMin, rangeMax));
			gauss.plotOn(frame,Name("mycurve"));
			frame->Draw();	
			chi2 = frame->chiSquare(2);	

			latexLabel.SetTextSize(0.030);
                	latexLabel.SetNDC();
			latexLabel.DrawLatex(0.17, 0.88,Form("ECAL %s, %s r9",eta.c_str(),r9.c_str()));
			latexLabel.DrawLatex(0.17, 0.83,Form("#color[4]{m_{0} = %f #pm %f}",meanG.getVal(),meanG.getError()));
			latexLabel.DrawLatex(0.17, 0.78,Form("#color[4]{#sigma = %f #pm %f}",sigmaG.getVal(),sigmaG.getError()));
			latexLabel.DrawLatex(0.17, 0.73,Form("#color[4]{#chi^{2} / ndf =  %f}",chi2));
				

			
			plotsRecording(directoryName, fileName, c1);
			
			fitFunctionSystematics = fabs(reducedFitParameters[1] - meanG.getVal());
				
			system(Form("mkdir -p %s",directoryName.c_str()));


        		ofstream summaryFile(Form("%sSummary_fitFunctionSystematics.txt",directoryName.c_str()), ios::app);

        		summaryFile << fitVariable << " " << dataType << " >> " << r9 << " r9 " << eta << ", trueModel = " << trueModel <<", testModel = " << testModel << ", systematics : " << fitFunctionSystematics * 100 << " %" << endl;;

        		summaryFile.close();


		}
		c1->Clear();
		reducedFitParameters.erase(reducedFitParameters.begin(),reducedFitParameters.end());
	}


	return 0;
}
 
