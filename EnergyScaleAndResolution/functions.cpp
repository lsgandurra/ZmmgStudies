#include "functions.h"

// --- Test the existence of a file --- //
bool FileExists(const char* FileName)
{
	FILE* fp = NULL;

    	//will not work if you do not have read permissions

    	//to the file, but if you don't have read, it

    	//may as well not exist to begin with.

    	fp = fopen( FileName, "rb" );
    	if( fp != NULL )
    	{
        	fclose( fp );
        	return true;
    	}

    	return false;
}


// --- Transform a 'double' into a 'string' --- //
string doubleToString(double x)
{

        std::string s;
        {
                std::ostringstream oss;
                oss << x;
                s = oss.str();
        }
        std::cout << "x = " << x << " s = " << s << std::endl;

        return s;
}


// --- Return the number of rows in a file --- //
int rowsNumberInFile(string filename)
{
        ifstream in(filename.c_str()); 

        string row = "";
        int nbRows = 0;

        while(std::getline(in, row)) nbRows++;

        in.close(); 

        return nbRows;
}

// --- To save the plots --- //
//void plotsRecording(string directoryName, string fileName, string eta, int iteration, TCanvas * c1)
void plotsRecording(string directoryName, string fileName, TCanvas * c1)
{    
     
	//directoryName += eta;
	//mkdir(directoryName.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	system(Form("mkdir -p %s", directoryName.c_str()));
	//if(iteration != 10000) fileName += Form("%d",iteration);
        c1->Print(Form("%s%s.root",directoryName.c_str(),fileName.c_str()));
        c1->Print(Form("%s%s.C",directoryName.c_str(),fileName.c_str()));
        c1->Print(Form("%s%s.pdf",directoryName.c_str(),fileName.c_str()));
        //c1->Print(Form("%s%s.ps",directoryName.c_str(),fileName.c_str()));
        c1->Print(Form("%s%s.png",directoryName.c_str(),fileName.c_str()));
        return;
}


// --- Compute the effective sigma of a given histogram --- //
Double_t effSigma(TH1 * hist)
{

  TAxis *xaxis = hist->GetXaxis();
  Int_t nb = xaxis->GetNbins();
  if(nb < 10) {
    cout << "effsigma: Not a valid histo. nbins = " << nb << endl;
    return 0.;
  }
  
  Double_t bwid = xaxis->GetBinWidth(1);
  if(bwid == 0) {
    cout << "effsigma: Not a valid histo. bwid = " << bwid << endl;
    return 0.;
  }
  Double_t xmax = xaxis->GetXmax();
  Double_t xmin = xaxis->GetXmin();
  Double_t ave = hist->GetMean();
  Double_t rms = hist->GetRMS();

  Double_t total=0.;
  for(Int_t i=0; i<nb+2; i++) {
    total+=hist->GetBinContent(i);
  }
  if(total < 100.) {
    cout << "effsigma: Too few entries " << total << endl;
    return 0.;
  }
  Int_t ierr=0;
  Int_t ismin=999;
  
  Double_t rlim=0.683*total;
  Int_t nrms=rms/(bwid);    // Set scan size to +/- rms
  if(nrms > nb/10) nrms=nb/10; // Could be tuned...

  Double_t widmin=9999999.;
  for(Int_t iscan=-nrms;iscan<nrms+1;iscan++) { // Scan window centre
    Int_t ibm=(ave-xmin)/bwid+1+iscan;
    Double_t x=(ibm-0.5)*bwid+xmin;
    Double_t xj=x;
    Double_t xk=x;
    Int_t jbm=ibm;
    Int_t kbm=ibm;
    Double_t bin=hist->GetBinContent(ibm);
    total=bin;
    for(Int_t j=1;j<nb;j++){
      if(jbm < nb) {
        jbm++;
        xj+=bwid;
        bin=hist->GetBinContent(jbm);
        total+=bin;
        if(total > rlim) break;
      }
      else ierr=1;
      if(kbm > 0) {
        kbm--;
        xk-=bwid;
        bin=hist->GetBinContent(kbm);
        total+=bin;
        if(total > rlim) break;
      }
      else ierr=1;
    }
    Double_t dxf=(total-rlim)*bwid/bin;
    Double_t wid=(xj-xk+bwid-dxf)*0.5;
    if(wid < widmin) {
      widmin=wid;
      ismin=iscan;
    }   
  }
  if(ismin == nrms || ismin == -nrms) ierr=3;
  if(ierr != 0) cout << "effsigma: Error of type " << ierr << endl;
  
  return widmin;
  
}

// --- Give the range limit corresponding to a given percentage --- //
void rangeEstimator(double percentage, TChain * chain, TString cut, double &minRange, double &maxRange, string variableName, float variable, vector <double> &fitParameters)
{

	cout << endl << "coucou 3 in 1" << endl;

        //TChain * reducedChain = (TChain *) chain->CopyTree(cut);

        TChain * reducedChain = chain;
	reducedChain->SetBranchAddress(variableName.c_str(),&variable);

        vector <float> variableVector;

	fitParameters.push_back(reducedChain->GetEntries());
        for (int ievt = 0 ; ievt < reducedChain->GetEntries() ; ievt++)
        {
                reducedChain->GetEntry(ievt);
                variableVector.push_back(variable);

        }
	cout << endl << "coucou 3 in 2" << endl;

        sort(variableVector.begin(), variableVector.end());
        size_t interval_entries = TMath::Ceil(percentage * variableVector.size());
	
        vector<float>::iterator lower = variableVector.begin();
        vector<float>::iterator upper = variableVector.begin() + interval_entries - 1; 
        double dx = *upper - *lower;
        for(vector<float>::iterator first = lower, last = upper; last < variableVector.end(); first++, last++)
        {
                if((*last - *first) < dx)
                {

                        lower = first;
                        upper = last;
                        dx = *upper - *lower;
                }
                      
        }
        minRange = *lower;
        maxRange = *upper;
	cout << endl << "coucou 3 in 3" << endl;
	//reducedChain->Clear();
}


// --- chiSquare estimation between an histogram and a fit function --- //
Double_t chiSquare(RooPlot* plot, char* pdfname, char* histname, vector <double> &fitParameters, int nbPar) //int* fewBins 
{
  // Calculate the chi^2/NDOF of this curve with respect to the histogram
  // 'hist' accounting nFitParam floating parameters in case the curve
  // was the result of a fit

  // Find curve object
  RooCurve* curve = (RooCurve*) plot->findObject(pdfname, RooCurve::Class());
  //RooCurve* curve = plot->getCurve(pdfname);  
  //curve->Print();

  if (!curve) {
    cout<<endl << "cit::RooChi2Calculator(plotname=" << plot->GetName()
         << ")::chiSquare(..) cannot find curve" << endl ;
    return 0 ;
  }

  // Find histogram object
  RooHist* hist = (RooHist*) plot->findObject(histname, RooHist::Class()) ;
  //RooHist* hist = plot->getHist(histname);

  if (!hist) {
    cout<<endl << "cit::RooChi2Calculator(plotname=" << plot->GetName()
         << ")::chiSquare(..) cannot find histogram" << endl ;
    return 0 ;
  }


  Int_t i,np = hist->GetN() ;
  Double_t x,y,/*eyl,eyh,*/ xl,xh ;

  // Find starting and ending bin of histogram based on range of RooCurve
  Double_t xstart,xstop ;

#if ROOT_VERSION_CODE >= ROOT_VERSION(4,0,1)
  curve->GetPoint(0,xstart,y) ;
  curve->GetPoint(curve->GetN()-1,xstop,y) ;
#else
  const_cast<RooCurve*>(curve)->GetPoint(0,xstart,y) ;
  const_cast<RooCurve*>(curve)->GetPoint(curve->GetN() - 1,xstop,y) ;
#endif

  Int_t nbin(0) ;

  Double_t chisq(0) ;
  for (i=0 ; i<np ; i++) {   

    // Retrieve histogram contents
    hist->GetPoint(i,x,y) ;
    xl = x - hist->GetEXlow()[i] ;
    xh = x + hist->GetEXhigh()[i] ;
    // eyl = hist->GetEYlow()[i] ;
    // eyh = hist->GetEYhigh()[i] ;

    // Check if the whole bin is in range of curve
    if (xl < xstart || xstop < xh) continue ;

    //if(y != 0 && y < 35.0)
    //if(y == 0)
    if(y < 5.0)
    {
    	cout<<endl<<"Too few entries : "<<y<<" in the bin : "<<i<<"  >> Need to reduce the binning for the p-value calculation!"<<endl;
	//*fewBins = 1;
	//break;
	continue;	
    
    }
    //else *fewBins = 0;

    nbin++ ;

    // Integrate function over this bin.
    // Start a hack to work around a bug in RooCurve::interpolate
    // that sometimes gives a wrong result.
    Double_t avg = curve->average(xl, xh);
    Double_t avg2 = 0.5 * (curve->average(xl, x) + curve->average(x, xh));
    if (avg + avg2 > 0 &&
	(avg2 - avg) / (avg2 + avg) > 0.1) {
      avg = curve->interpolate(x);
    }
    // End of hack around the bug in RooCurve::interpolate

    // JV: Adjust observed and expected number of events for bin width to represent
    // number of events.
    Double_t norm = (xh - xl) / plot->getFitRangeBinW();
    y *= norm;
    avg *= norm;

    if (avg < 5.) {
      cout << "cit::RooChi2Calculator(plotname=" << plot->GetName()
			    << ")::chiSquare(..) expectation in bin "
			    << i << " is " << avg << " < 5!" << endl ;
    }

    // JV: Use the expected number of events for the y uncertainty,
    // See (33.34) of http://pdg.lbl.gov/2011/reviews/rpp2011-rev-statistics.pdf

    // Add pull^2 to chisq
    if (avg != 0) {      
      Double_t resid = y - avg;
      chisq += (resid * resid / avg) ;
    }
  }

  //int nFitParam = fitParameters.size() / 2; //FIXME !!
  //int nFitParam = fitParameters[2];
  int nFitParam = nbPar;
  fitParameters.push_back(chisq / (nbin - nFitParam));
  fitParameters.push_back(nbin - nFitParam);
  fitParameters.push_back(TMath::Prob(chisq, nbin - nFitParam));

	cout<<endl<<"chisq = "<<chisq<<", nbin = "<<nbin<<", nFitParam = "<<nFitParam<<endl;

  return chisq / (nbin - nFitParam) ;
}

// --- RooHist of chi2 residuals --- //  
RooHist* residHist(RooPlot* plot, char *histname, char* curvename, bool normalize, string recordingDirectory, int iteration)
{
  if(normalize == true)
  {
	if(FileExists(Form("%sPullsX_%d.txt",recordingDirectory.c_str(), iteration)))
	{	
		system(Form("rm %sPullsX_%d.txt",recordingDirectory.c_str(), iteration));
	}
	system(Form("touch %sPullsX_%d.txt",recordingDirectory.c_str(), iteration));

	if(FileExists(Form("%sPullsErrorX_%d.txt",recordingDirectory.c_str(), iteration)))
        {           
                system(Form("rm %sPullsErrorX_%d.txt",recordingDirectory.c_str(), iteration));
        }   
        system(Form("touch %sPullsErrorX_%d.txt",recordingDirectory.c_str(), iteration));

	if(FileExists(Form("%sPullsY_%d.txt",recordingDirectory.c_str(), iteration)))
        {           
                system(Form("rm %sPullsY_%d.txt",recordingDirectory.c_str(), iteration));
        }   
        system(Form("touch %sPullsY_%d.txt",recordingDirectory.c_str(), iteration));

	if(FileExists(Form("%sPullsErrorY_%d.txt",recordingDirectory.c_str(), iteration)))
        {           
                system(Form("rm %sPullsErrorY_%d.txt",recordingDirectory.c_str(), iteration));
        }   
        system(Form("touch %sPullsErrorY_%d.txt",recordingDirectory.c_str(), iteration));	
  }

  // Create and return RooHist containing  residuals w.r.t to given curve->
  // If normalize is true, the residuals are normalized by the histogram
  // errors creating a RooHist with pull values

  // Find curve object
  RooCurve* curve = (RooCurve*) plot->findObject(curvename, RooCurve::Class());
  if (!curve) {
    cout << "cit::RooChi2Calculator(plotname=" << plot->GetName()
         << ")::residHist(..) cannot find curve" << endl ;
    return 0 ;
  }

  // Find histogram object
  RooHist* hist = (RooHist*) plot->findObject(histname, RooHist::Class()) ;
  if (!hist) {
    cout << "cit::RooChi2Calculator(plotname=" << plot->GetName()
         << ")::residHist(..) cannot find histogram" << endl ;
    return 0 ;
  }

  // Copy all non-content properties from hist
  RooHist* ret = new RooHist(plot->getFitRangeBinW()) ;
  if (normalize) {
    ret->SetName(Form("pull_%s_%s", hist->GetName(), curve->GetName())) ;
    ret->SetTitle(Form("Pull of %s and %s", hist->GetTitle(), curve->GetTitle())) ;
  } else {
    ret->SetName(Form("resid_%s_%s", hist->GetName(), curve->GetName())) ;
    ret->SetTitle(Form("Residual of %s and %s", hist->GetTitle(), curve->GetTitle())) ;
  }

  // Determine range of curve
  Double_t xstart, xstop, y ;
#if ROOT_VERSION_CODE >= ROOT_VERSION(4,0,1)
  curve->GetPoint(0, xstart, y) ;
  curve->GetPoint(curve->GetN()-1,xstop,y) ;
#else
  const_cast<RooCurve&>(curve)->GetPoint(0, xstart, y) ;
  const_cast<RooCurve&>(curve)->GetPoint(curve->GetN()-1, xstop, y) ;
#endif
  // cout << "cit::RooChi2Calculator::residHist dumping curve:\n";
  // for (int i=0; i<curve->GetN(); ++i){
  //   Double_t xi, yi;
  //   curve->GetPoint(i, xi, yi);
  //   printf("i=%d x,y: %.3g, %.3g\n", i, xi, yi);
  // }

  // cout << "cit::RooChi2Calculator::residHist  adding bins with error:\n";

  // Add histograms, calculate Poisson confidence interval on sum value
  for(Int_t i=0 ; i < hist->GetN() ; i++) {
    Double_t x, point;
#if ROOT_VERSION_CODE >= ROOT_VERSION(4,0,1)
    hist->GetPoint(i,x,point) ;
#else
    const_cast<RooHist&>(hist)->GetPoint(i,x,point) ;
#endif
    Double_t xl = x - hist->GetErrorXlow(i);
    Double_t xh = x + hist->GetErrorXhigh(i);

    // Only calculate pull for bins inside curve range
    if (xl < xstart || xstop < xh) continue ;

    Double_t norm = (xh - xl) / plot->getFitRangeBinW();
    point *= norm;

    // Start a hack to work around a bug in RooCurve::interpolate
    // that sometimes gives a wrong result.
    Double_t avg = curve->average(xl, xh);
    Double_t avg2 = 0.5 * (curve->average(xl, x) + curve->average(x, xh));
    Double_t yexpected;
    if (avg + avg2 > 0 && (avg2 - avg) / (avg2 + avg) > 0.1) {
      yexpected = curve->interpolate(x);
    } else {
      yexpected = avg;
    }
    // End of hack around the bug in RooCurve::interpolate

    // Correct the expected number of events in this bin for the non-uniform
    // bin width.
    yexpected *= norm;

    Double_t yy = point - yexpected;
    // Normalize to the number of events per bin taking into account
    // variable bin width.
    Double_t dy = TMath::Sqrt(yexpected);
    if (normalize) {
	if (dy==0.) {
	  cout << "cit::RooChi2Calculator::residHist(histname ="
               << hist->GetName() << ", ...) WARNING: point "
               << i << " has zero error, setting residual to zero"
               << endl ;
	  yy=0 ;
	  dy=0 ;
	} else {
	  yy /= dy;
	  dy = 1.;
	}
    }
    // printf("bin=%3d n=%5.3g nu=%5.3g x=%5.3g .. %5.3g y=%5.3g +/- %5.3g "
    //	   "norm=%5.3g\n", i, point, yexpected, xl, xh, yy, dy, norm);
    ret->addBinWithError(x,yy,dy,dy);
  
    if(normalize == true)
    {
        ofstream filePullsX(Form("%sPullsX_%d.txt",recordingDirectory.c_str(), iteration), ios::app);
        filePullsX << yy <<endl;
        filePullsX.close();

        ofstream filePullsErrorX(Form("%sPullsErrorX_%d.txt",recordingDirectory.c_str(), iteration), ios::app);
        filePullsErrorX << dy <<endl;
        filePullsErrorX.close();

        ofstream filePullsY(Form("%sPullsY_%d.txt",recordingDirectory.c_str(), iteration), ios::app);
        filePullsY << x <<endl;
        filePullsY.close();

        ofstream filePullsErrorY(Form("%sPullsErrorY_%d.txt",recordingDirectory.c_str(), iteration), ios::app);
        filePullsErrorY << dy <<endl;
        filePullsErrorY.close();

    }

  }
  return ret;
}

// --- RooHist of chi2 pulls --- //  
RooHist* pullHist(RooPlot* plot_, char* histname, char* pdfname, string dossierSauvegardePull, int iteration)
{  
        return residHist(plot_, histname, pdfname, true, dossierSauvegardePull, iteration); 
}



/*
void RangeEstimator(double percentage, double centralValue, TChain * chain, int Endcaps, double * MinRange, double * MaxRange, double Varmin, double Varmax)
{

	float mmg_s;
        float Photon_SC_Eta;
        float Photon_r9;
        float isJanLooseMMG;
        float isMultipleCandidate;
        float Photon_Et;
        chain->SetBranchAddress("mmg_s",&mmg_s);
        chain->SetBranchAddress("Photon_SC_Eta",&Photon_SC_Eta);
        chain->SetBranchAddress("Photon_r9",&Photon_r9);
        chain->SetBranchAddress("isJanLooseMMG",&isJanLooseMMG);
        chain->SetBranchAddress("isMultipleCandidate",&isMultipleCandidate);
        chain->SetBranchAddress("Photon_Et",&Photon_Et);
        vector <float> mmg_sVector;

	for (int ievt = 0 ; ievt < chain->GetEntries() ; ievt++)
        {
                chain->GetEntry(ievt);

                if(Endcaps == 0)
                {

                        if(( (fabs(Photon_SC_Eta)) > 0.018 || ( (fabs(Photon_SC_Eta)) < 0.423 && (fabs(Photon_SC_Eta)) > 0.461 ) || ( (fabs(Photon_SC_Eta)) < 0.770 && (fabs(Photon_SC_Eta)) > 0.806 ) || ( (fabs(Photon_SC_Eta)) < 1.127 && (fabs(Photon_SC_Eta)) > 1.163 ) ) && (fabs(Photon_SC_Eta)) < 1.4442 && isMultipleCandidate == 0 && isJanLooseMMG != 0 && Photon_r9 > 0.94 && Photon_Et > Varmin && Photon_Et <= Varmax)
                        {

                                mmg_sVector.push_back(mmg_s);

                        }


                }

                if(Endcaps == 1)
                {

        //cout<<endl<<BremVector[1000]<<endl;
                        if((fabs(Photon_SC_Eta)) > 1.566 && isJanLooseMMG != 0 && isMultipleCandidate == 0 && Photon_r9 > 0.95 && Photon_Et > Varmin && Photon_Et <= Varmax)  mmg_sVector.push_back(mmg_s);

                }


        }
        sort(mmg_sVector.begin(), mmg_sVector.end());

	cout<<endl<<"mmg_sVector.size() = "<<mmg_sVector.size()<<endl;

	int meanVector = 0;
	for(int h = 0; h < mmg_sVector.size(); h++)
	{
		if(mmg_sVector[h] >= centralValue)
		{
			meanVector = h;
			h = mmg_sVector.size();
			break;
		}


	}

        
	cout<<endl<<"meanVector = "<<meanVector<<endl;

        double stopVectorMoins = meanVector - mmg_sVector.size() * percentage / 2.0;
        double stopVectorPlus = meanVector + mmg_sVector.size() * percentage / 2.0;
        cout<<endl<<"stopVectorPlus = "<<stopVectorPlus<<endl;
        int intStopVectorMoins = (int) stopVectorMoins;
        int intStopVectorPlus = (int) stopVectorPlus;
        cout<<endl<<"intStopVectorMoins = "<<intStopVectorMoins<<endl;
        cout<<endl<<"intStopVectorPlus = "<<intStopVectorPlus<<endl;

        cout<<endl<<"mmg_sVector[meanVector] = "<<mmg_sVector[meanVector]<<endl;

        cout<<endl<<"Le range est compris entre : "<<mmg_sVector[intStopVectorMoins]<<" et "<<mmg_sVector[intStopVectorPlus]<<endl;

	*MinRange = mmg_sVector[intStopVectorMoins];
	*MaxRange = mmg_sVector[intStopVectorPlus]; 


}
*/
/*
void RangeEstimator2(double percentage, TChain * chain, int Endcaps, double * MinRange, double * MaxRange, double Varmin, double Varmax)
{

	float mmg_s;
        float Photon_SC_Eta;
        float Photon_r9;
        float isJanLooseMMG;
        float isMultipleCandidate;
        float Photon_Et;
        chain->SetBranchAddress("mmg_s",&mmg_s);
        chain->SetBranchAddress("Photon_SC_Eta",&Photon_SC_Eta);
        chain->SetBranchAddress("Photon_r9",&Photon_r9);
        chain->SetBranchAddress("isJanLooseMMG",&isJanLooseMMG);
        chain->SetBranchAddress("isMultipleCandidate",&isMultipleCandidate);
        chain->SetBranchAddress("Photon_Et",&Photon_Et);
        vector <float> mmg_sVector;

        for (int ievt = 0 ; ievt < chain->GetEntries() ; ievt++)
        {
                chain->GetEntry(ievt);

                if(Endcaps == 0)
                {

                        if(( (fabs(Photon_SC_Eta)) > 0.018 || ( (fabs(Photon_SC_Eta)) < 0.423 && (fabs(Photon_SC_Eta)) > 0.461 ) || ( (fabs(Photon_SC_Eta)) < 0.770 && (fabs(Photon_SC_Eta)) > 0.806 ) || ( (fabs(Photon_SC_Eta)) < 1.127 && (fabs(Photon_SC_Eta)) > 1.163 ) ) && (fabs(Photon_SC_Eta)) < 1.4442 && isMultipleCandidate == 0 && isJanLooseMMG != 0 && Photon_r9 > 0.94 && Photon_Et > Varmin && Photon_Et <= Varmax)
                        {

                                mmg_sVector.push_back(mmg_s);

                        }


                }

                if(Endcaps == 1)
                {

        //cout<<endl<<BremVector[1000]<<endl;
                        if((fabs(Photon_SC_Eta)) > 1.566 && isJanLooseMMG != 0 && isMultipleCandidate == 0 && Photon_r9 > 0.95 && Photon_Et > Varmin && Photon_Et <= Varmax)  mmg_sVector.push_back(mmg_s);

                }


        }
        sort(mmg_sVector.begin(), mmg_sVector.end());

        double Min = mmg_sVector.size() * (percentage / 4.0);
        int MinInt = (int) Min;
        double Max = mmg_sVector.size() * (3.0 * percentage / 4.0);
        int MaxInt = (int) Max;

        *MinRange = mmg_sVector[MinInt];
        *MaxRange = mmg_sVector[MaxInt];

}
*/

/*
void SymetricRangeEstimator(TChain * chain, double centralValue, double * MinRange, double * MaxRange, double Entries, double percentage, TString temp)
{

	double min = 0.0;
	double max = 0.0;
	string minresult;
	string maxresult;
	std::ostringstream oss;
	std::ostringstream oss2;	
	TString temp2;
	

	for(double sigma = 0.001; sigma < 1.0; sigma += 0.001)
	{
		temp2.Clear();
		oss.str("");
		oss2.str("");
		minresult.clear();
		maxresult.clear();

		min = centralValue - sigma; 
		max = centralValue + sigma;

        	oss << min;
        	minresult = oss.str();
		
                oss2 << max; 
                maxresult = oss2.str();	
		
		TH1D *histo = new TH1D("histo","histo", 200, 0.0, 2.0);
		
		temp2 = temp;
		temp2 += " && mmg_s > ";
		temp2 += minresult;
		temp2 += " && mmg_s < ";
		temp2 += maxresult;

		chain->Draw("mmg_s>>histo", temp2);
		
		if(histo->GetEntries() >= (Entries * percentage))
		{
			*MinRange = centralValue - sigma;
			*MaxRange = centralValue + sigma;
			sigma = 1.0;
		}

		histo->Delete();	
	}



}
*/

/*
void SymetricRangeEstimator2(TChain * chain, double lastX, double * MinRange, double * MaxRange, double Entries, double percentage, TString temp)
{

        double min = 0.0;
        double max = 0.0;
        int itermin = 0;
        int itermax = 0;
	int loop1 = 0;
	int loop2 = 0;
        string minresult;
        string maxresult;
        std::ostringstream oss;
        std::ostringstream oss2;
        TString temp2;

	double limite = ((1.0 - percentage) / 2.0) * Entries;

        for(double sigma = 0.0; sigma < 2.0; sigma += 0.01)
        {
                temp2.Clear();
                oss.str("");
                oss2.str("");
                minresult.clear();
                maxresult.clear();

                min = sigma;
                max = lastX - sigma;

                oss << min;
                minresult = oss.str();

                oss2 << max;
                maxresult = oss2.str();

                TH1D *histo = new TH1D("histo","histo", 2000, 0.0, 2.0);
                TH1D *histo2 = new TH1D("histo2","histo2", 2000, 0.0, 2.0);

		temp2 = temp;
                temp2 += " && mmg_s < ";
                temp2 += minresult;

                chain->Draw("mmg_s>>histo", temp2);


                if(histo->GetEntries() >= limite && itermin == 0)
		{
                        *MinRange = sigma;
                        itermin = 1;
		}

                temp2.Clear();

                temp2 = temp;
                temp2 += " && mmg_s > ";
                temp2 += maxresult;

		chain->Draw("mmg_s>>histo2", temp2);


                if(histo2->GetEntries() >= limite && itermax == 0)
		{
		        *MaxRange =  lastX - sigma;
                        itermax = 1;
                }

		if (itermin == 1 && itermax == 1) sigma = 2.0;

		

		if((histo->GetEntries() < (limite * 0.1)) && (histo2->GetEntries() < (limite * 0.1)) && loop1 == 0)
		{
			 sigma +=0.05;
		}
		else loop1 = 1;
		

                histo->Delete();
                histo2->Delete();
        }


}

*/


/*
void SymetricRangeEstimator3(TChain * chain, double centralValue, double * MinRange, double * MaxRange, double Entries, double percentage, TString temp)
{
	

        double min = 0.0;
        double max = 0.0;
	int iterMin = 0;
	int iterMax = 0;
        string minresult;
        string maxresult;
        std::ostringstream oss;
        std::ostringstream oss2;
        TString temp2;


        for(double sigma = 0.001; sigma < 1.0; sigma += 0.001)
        //for(double sigma = 0.01; sigma < 1.0; sigma += 0.01)
	{
		temp2.Clear();
                oss.str("");
                oss2.str("");
                minresult.clear();
                maxresult.clear();

                min = centralValue - sigma; 
                max = centralValue + sigma;

                oss << min; 
                minresult = oss.str();
      
                oss2 << max; 
                maxresult = oss2.str(); 
      
		if(iterMin == 0)
		{
                	TH1D *histo = new TH1D("histo","histo", 200, 0.0, 2.0);
      
                	temp2 = temp;
                	temp2 += " && mmg_s > ";
                	temp2 += minresult;
                	temp2 += " && mmg_s < ";
                	temp2 += centralValue;

                	chain->Draw("mmg_s>>histo", temp2);
			if(histo->GetEntries() >= (Entries * (percentage / 2.0)))
                	{
                        	*MinRange = centralValue - sigma;
                        	iterMin = 1;
                	}
			histo->Delete();
		}

		if(iterMax == 0)
                {		

			TH1D *histo2 = new TH1D("histo2","histo2", 200, 0.0, 2.0);
                	
			temp2 = temp;
                	temp2 += " && mmg_s > ";
                	temp2 += centralValue;
                	temp2 += " && mmg_s < ";
                	temp2 += maxresult;
                
			chain->Draw("mmg_s>>histo2", temp2);            
	
			if(histo2->GetEntries() >= (Entries * (percentage / 2.0)))
                	{
                        	*MaxRange = centralValue + sigma;
                        	iterMax = 1;
                	}

                	histo2->Delete();
		}
		
		if(iterMin == 1 && iterMax == 1) sigma = 1.0;

        }
		if(iterMin == 0) *MinRange = 0.0;
		if(iterMax == 0) *MaxRange = 1.5;		

}

*/


/*

void chi2Chain(char * buffer, double chi2)
{
	
	string chain("#chi^{2} / ndf = ");
        string finalChain;
	
	finalChain = chain + doubleToString(chi2);

        //size_t size = chaineFinale.size() + 1; 
        //char * buffer = new char[ size ];
        //strncpy( buffer, chaineFinale.c_str(), size );
	strncpy( buffer, finalChain.c_str(), 25 );

}

double sigmaR(TF1* function, double xmin, double xmax)
{

	double y = function->GetMaximum(xmin, xmax) * 1.0/exp(1.0);
	double maxX = function->GetMaximumX(xmin, xmax);
	double sigma = function->GetX(y, maxX, xmax) - maxX;

	return sigma;
}

double sigmaL(TF1* function, double xmin, double xmax)
{

        double y = function->GetMaximum(xmin, xmax) * 1.0/exp(1.0);
        double maxX = function->GetMaximumX(xmin, xmax);
	double sigma = maxX - function->GetX(y, xmin, maxX);

        return sigma;
}
*/

