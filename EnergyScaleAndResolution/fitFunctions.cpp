#include "fitFunctions.h"


void voigtian(RooDataSet *dataset, RooDataSet *dataset2, RooRealVar &variable, RooPlot *fitFrame, RooBinning b, double rangeMin, double rangeMax, vector <double> &fitParameters)
{
	double median = (rangeMax + rangeMin) /2.0;
	double sigH = (rangeMax-rangeMin) /2.0;
        cout<<endl<<"median = "<<median;
	cout<<endl<<"sigH = "<<sigH;
	//RooRealVar mean("mean","mean",0.0,-0.1,0.1); //FIXME
        RooRealVar mean("mean","mean",median,rangeMin,rangeMax);
	//RooRealVar sigma("sigma","sigma",0.5,0.0,1.0); //FIXME
        RooRealVar sigma("sigma","sigma",0.5,0.0,sigH);
	//RooRealVar width("width","width",0.5,0.0,1.0); //FIXME
	RooRealVar width("width","width",0.5,0.0,sigH);	

	cout<<endl<<"coucou avant pdf"<<endl;

        RooVoigtian * pdf = new RooVoigtian("pdf","Voigtian",variable,mean,sigma,width);

	cout<<endl<<"coucou apres pdf"<<endl;
        dataset->plotOn(fitFrame,Name("myhist"),Binning(b),DataError(RooAbsData::SumW2));

        RooFitResult * res = pdf->fitTo(*dataset, Range(rangeMin, rangeMax),Save(),SumW2Error(kTRUE));
	res->Print();
	cout<<endl<<"coucou apres fitTo"<<endl;

	//minNll = res->minNll();
        
	pdf->plotOn(fitFrame,Name("mycurve"));


	fitParameters.push_back(3); // nb of fit parameters

	fitParameters.push_back(mean.getVal());
	fitParameters.push_back(mean.getError());

	fitParameters.push_back(sigma.getVal());
        fitParameters.push_back(sigma.getError());
	
	fitParameters.push_back(width.getVal());
        fitParameters.push_back(width.getError());		
	

}

void cruijff(RooDataSet *dataset, RooDataSet *dataset2, RooRealVar &variable, RooPlot *fitFrame, RooBinning b, double rangeMin, double rangeMax, vector <double> &fitParameters)
{

        RooRealVar m0("m0","m0",0.0,-0.1,0.1);
        RooRealVar sigmaL("sigmaL","sigmaL",0.5,0.0,1.0);
        RooRealVar sigmaR("sigmaR","sigmaR",0.5,0.0,1.0);
        RooRealVar alphaL("alphaL","alphaL",0.5,0.0,1.0);
        RooRealVar alphaR("alphaR","alphaR",0.5,0.0,1.0);


        RooCruijff * pdf = new RooCruijff("pdf","Cruijff",variable,m0,sigmaL,sigmaR,alphaL,alphaR);

        dataset->plotOn(fitFrame,Name("myhist"),Binning(b),DataError(RooAbsData::SumW2));

        RooFitResult * res = pdf->fitTo(*dataset, Range(rangeMin, rangeMax),Save(),SumW2Error(kTRUE));
        res->Print();

        //minNll = res->minNll();

        pdf->plotOn(fitFrame,Name("mycurve"));


        fitParameters.push_back(5); // nb of fit parameters

        fitParameters.push_back(m0.getVal());
        fitParameters.push_back(m0.getError());

        fitParameters.push_back(sigmaL.getVal());
        fitParameters.push_back(sigmaL.getError());

        fitParameters.push_back(sigmaR.getVal());
        fitParameters.push_back(sigmaR.getError());

        fitParameters.push_back(alphaL.getVal());
        fitParameters.push_back(alphaL.getError());

        fitParameters.push_back(alphaR.getVal());
        fitParameters.push_back(alphaR.getError());


}

void voigtianXcb(RooDataSet *dataset, RooDataSet *dataset2, RooRealVar &variable, RooPlot *fitFrame, RooBinning b, double rangeMin, double rangeMax, vector <double> &fitParameters)
{

	// --- voigtian --- //

	RooRealVar meanV("meanV","meanV",91.1876,"GeV");
        meanV.setConstant(kTRUE);
	RooRealVar sigmaV("sigmaV","sigmaV",0.5,0.0,5.0);
        //RooRealVar sigmaV("sigmaV","sigmaV",1.21,"GeV");
	//sigmaV.setConstant(kTRUE);
	RooRealVar widthV("widthV","widthV",2.4952,"GeV");
	widthV.setConstant(kTRUE);

        RooVoigtian * voigtian = new RooVoigtian("voigtian","voigtian",variable,meanV,sigmaV,widthV);

	// --- CB --- //

	RooRealVar meanCB("meanCB","meanCB",0.43,0.0,3.0);
        RooRealVar sigmaCB("sigmaCB","sigmaCB",1.655,0.0,3.5);
        RooRealVar alphaCB("alphaCB","alphaCB",2.18,-2.5,2.5);
        RooRealVar nCB ("nCB","nCB",8.51,0.0,20.0);

        RooCBShape * cb = new RooCBShape("cb","cb",variable,meanCB,sigmaCB,alphaCB,nCB) ;
        
	// --- Convolution --- //

	RooFFTConvPdf * pdf  = new RooFFTConvPdf("pdf","pdf",variable,*voigtian,*cb); 


        dataset->plotOn(fitFrame,Name("myhist"),Binning(b),DataError(RooAbsData::SumW2));

        RooFitResult * res = pdf->fitTo(*dataset, Range(rangeMin, rangeMax),Save(),SumW2Error(kTRUE));
        res->Print();

        //minNll = res->minNll();
     
        pdf->plotOn(fitFrame,Name("mycurve"));


        fitParameters.push_back(7); // nb of fit parameters

        fitParameters.push_back(meanV.getVal());
        fitParameters.push_back(meanV.getError());

        fitParameters.push_back(sigmaV.getVal());
        fitParameters.push_back(sigmaV.getError());
     
        fitParameters.push_back(widthV.getVal());
        fitParameters.push_back(widthV.getError()); 

	fitParameters.push_back(meanCB.getVal());
        fitParameters.push_back(meanCB.getError());

        fitParameters.push_back(sigmaCB.getVal());
        fitParameters.push_back(sigmaCB.getError());   

	fitParameters.push_back(alphaCB.getVal());
        fitParameters.push_back(alphaCB.getError());

        fitParameters.push_back(nCB.getVal());
        fitParameters.push_back(nCB.getError());
     
}

void bwXcb(RooDataSet *dataset, RooDataSet *dataset2, RooRealVar &variable, RooPlot *fitFrame, RooBinning b, double rangeMin, double rangeMax, vector <double> &fitParameters)
{

	// --- BW --- //

	//RooRealVar meanBW("meanBW","meanBW",91.1876,90,92);
	RooRealVar meanBW("meanBW","meanBW",91.1876,"GeV");
        meanBW.setConstant(kTRUE);
        RooRealVar widthBW("widthBW","widthBW",2.4952,"GeV");
	widthBW.setConstant(kTRUE);

	//RooRealVar meanBW("meanBW","meanBW",91.1876,89,93);
	//RooRealVar widthBW("widthBW","widthBW",2.4952,0.0,5.0);

	RooBreitWigner * bw = new RooBreitWigner("bw","bw",variable,meanBW,widthBW);

	// --- CB --- //

	RooRealVar meanCB("meanCB","meanCB",0.43,0.0,3.0);
        RooRealVar sigmaCB("sigmaCB","sigmaCB",1.655,0.0,3.5);
        RooRealVar alphaCB("alphaCB","alphaCB",2.18,-2.5,2.5);
        RooRealVar nCB ("nCB","nCB",8.51,0.0,20.0);

        RooCBShape * cb = new RooCBShape("cb","cb",variable,meanCB,sigmaCB,alphaCB,nCB) ;
        
	// --- Convolution --- //

	RooFFTConvPdf * pdf  = new RooFFTConvPdf("pdf","pdf",variable,*bw,*cb); 


        dataset->plotOn(fitFrame,Name("myhist"),Binning(b),DataError(RooAbsData::SumW2));

        RooFitResult * res = pdf->fitTo(*dataset, Range(rangeMin, rangeMax),Save(),SumW2Error(kTRUE));
        res->Print();

        //minNll = res->minNll();
     
        pdf->plotOn(fitFrame,Name("mycurve"));


        fitParameters.push_back(6); // nb of fit parameters

        fitParameters.push_back(meanBW.getVal());
        fitParameters.push_back(meanBW.getError());
     
        fitParameters.push_back(widthBW.getVal());
        fitParameters.push_back(widthBW.getError()); 

	fitParameters.push_back(meanCB.getVal());
        fitParameters.push_back(meanCB.getError());

        fitParameters.push_back(sigmaCB.getVal());
        fitParameters.push_back(sigmaCB.getError());   

	fitParameters.push_back(alphaCB.getVal());
        fitParameters.push_back(alphaCB.getError());

        fitParameters.push_back(nCB.getVal());
        fitParameters.push_back(nCB.getError());
     
}

void voigtianXgauss(RooDataSet *dataset, RooDataSet *dataset2, RooRealVar &variable, RooPlot *fitFrame, RooBinning b, double rangeMin, double rangeMax, vector <double> &fitParameters)
{

	// --- voigtian --- //

	RooRealVar meanV("meanV","meanV",91.1876,"GeV");
        meanV.setConstant(kTRUE);
	RooRealVar sigmaV("sigmaV","sigmaV",0.5,0.0,5.0);
        //RooRealVar sigmaV("sigmaV","sigmaV",1.21,"GeV");
	//sigmaV.setConstant(kTRUE);
	RooRealVar widthV("widthV","widthV",2.4952,"GeV");
	widthV.setConstant(kTRUE);

        RooVoigtian * voigtian = new RooVoigtian("voigtian","voigtian",variable,meanV,sigmaV,widthV);

	// --- Gaussian --- //

	RooRealVar meanG("meanG","meanG",0.0,0.0,2.0);
        RooRealVar sigmaG("sigmaG","sigmaG",0.5,0.0,2.0);

        RooGaussian * gauss = new RooGaussian("gauss","gauss",variable,meanG,sigmaG);
        
	// --- Convolution --- //

	RooFFTConvPdf * pdf  = new RooFFTConvPdf("pdf","pdf",variable,*voigtian,*gauss); 


        dataset->plotOn(fitFrame,Name("myhist"),Binning(b),DataError(RooAbsData::SumW2));

        RooFitResult * res = pdf->fitTo(*dataset, Range(rangeMin, rangeMax),Save(),SumW2Error(kTRUE));
        res->Print();

        //minNll = res->minNll();
     
        pdf->plotOn(fitFrame,Name("mycurve"));


        fitParameters.push_back(5); // nb of fit parameters

        fitParameters.push_back(meanV.getVal());
        fitParameters.push_back(meanV.getError());

        fitParameters.push_back(sigmaV.getVal());
        fitParameters.push_back(sigmaV.getError());
     
        fitParameters.push_back(widthV.getVal());
        fitParameters.push_back(widthV.getError()); 

	fitParameters.push_back(meanG.getVal());
        fitParameters.push_back(meanG.getError());

        fitParameters.push_back(sigmaG.getVal());
        fitParameters.push_back(sigmaG.getError());   

}






