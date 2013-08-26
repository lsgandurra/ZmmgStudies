//Script by Louis Sgandurra (February 2013)
//Select the best fit for all categories and different fit ranges.
//Compute fit range systematics
#include "functions.h"

int main(int argc, char *argv[])
{
	for(int iarg = 0 ; iarg < argc; iarg++)
        {
                cout << "argv[" << iarg << "]= " << argv[iarg] << endl;
        }

        if( argc == 1 ) 
        {
                cerr << "arguments should be passed : directoryName, dataType, fitVariable, eta, r9, fitFunction, lowFitRange, highFitRange" <<endl; 
                return 1;

        }    

        string directoryName = "Results_v1";
        string dataType = "MC";
        string fitVariable = "mmg_s";
        string eta = "Barrel"; 
        string r9 = "all";
        string fitFunction = "voigtian";
        int lowFitRange = 60;
	int highFitRange = 100;

	if( argc > 1 ) directoryName = argv[1];
	if( argc > 2 ) dataType = argv[2];
	if( argc > 3 ) fitVariable = argv[3];
	if( argc > 4 ) eta = argv[4];
	if( argc > 5 ) r9 = argv[5];
	if( argc > 6 ) fitFunction = argv[6];
	if( argc > 7 ) 
	{
		std::stringstream ss ( argv[7] );
		ss >> lowFitRange;
	}
	if( argc > 8 ) 
	{
		std::stringstream ss ( argv[8] );
		ss >> highFitRange;
	}

	double temp_number = 0;

	int nbFits = 0;
        int nbRows = rowsNumberInFile(fileName.c_str());
        int nbFitParams = 0;

        ifstream tempFile(fileName.c_str());
        tempFile >> nbFitParams;
        tempFile >> nbFitParams;
        tempFile >> nbFitParams;
        tempFile.close();

        nbFits = nbRows / (2 + 6 + 2 * nbFitParams);


	string fileName = "";

	double * average_pvalue = new double[nbFits];
	double average_mean_error = 0;
	double ** median_mean_error_Tab = new double*[nbFits];
	double *  median_mean_error = new double[nbFits];
	double temp_var = 0;
	int best_fit = 0;
	double *min_error = new double[nbFits];
	
	int fit_number = (highFitRange - lowFitRange + 1) * nbFits;	

	double * meanTab = new double[fit_number];
	double * mean_errorTab = new double[fit_number];
	double * pvalue_tab = new double[fit_number];

	vector <double> fitParameters;
	vector <double> mean_errorVec;

	for(int i = lowFitRange; i <= highFitRange; i++)
	{
		fileName = Form("%s/%s/%dPercents/%s_%s_%sR9_%s/fitsInformationsRaw.txt",directoryName.c_str(),dataType.c_str(),i,fitVariable.c_str(),eta.c_str(),r9.c_str(),fitFunction.c_str());		

		//nbRows = rowsNumberInFile(fileName);

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


		for(int j = 0; j < nbFits; j++)
                {
                        temp_iter = j * ( 2 + 6 + fitParameters[2] * 2 ) + fitParameters[2] * 2 + 5;
                        pvalue_tab[ (i-lowFitRange) * nbFits + j ] = fitParameters[temp_iter];
                        cout<<endl<<"(i-lowFitRange) * nbFits + j = "<<(i-lowFitRange) * nbFits + j<<", fitParameters[temp_iter] = "<<fitParameters[temp_iter];
                	average_pvalue[j] += pvalue_tab[ (i-lowFitRange) * nbFits + j  ];

			temp_iter = j * ( 2 + 6 + fitParameters[2] * 2 ) + 1 + 2;
			meanTab[ (i-lowFitRange) * nbFits + j ] = fitParameters[temp_iter];
			mean_errorTab[ (i-lowFitRange) * nbFits + j ] = fitParameters[temp_iter+1];
			mean_errorVec.push_back(fitParameters[temp_iter+1]);

			if((i-lowFitRange) == 0) median_mean_error_Tab[j] = new double [fit_number];
			median_mean_error_Tab[j][i] = fitParameters[temp_iter+1];
		}
	
		fitParameters.erase(fitParameters.begin(),fitParameters.end());
	
		fitParametersFile.close();

	}


	for(int j = 0; j < nbFits; j++)
        {
		sort(&median_mean_error_Tab[j][0], &median_mean_error_Tab[j][fit_number]);	
		median_mean_error[j] = median_mean_error_Tab[j][(highFitRange - lowFitRange) / 2];	
        	average_pvalue[j] = average_pvalue[j] / fit_number;
		min_error[j] = median_mean_error[j] * 1.5;

		cout<< "median_mean_error[j] = "<<median_mean_error[j]<<endl;
		cout<< "min_error[j] = "<<min_error[j]<<endl;
		//cout<< "lowFitRange = "<<lowFitRange<<", highFitRange = "<<highFitRange<<endl;
	}


////// Modify after this line to take into account multiple fits /////
	
	///// Fit range selection /////
	
	for(int i = highFitRange; i >= lowFitRange; i--)
	{
		//if(i == 100 || i == 89 || i = 79 || i = 69) min_error = average_mean_error;

		if(pvalue_tab[i-lowFitRange] > 0.001)
		{
			if(i >= 90 && mean_errorTab[i-lowFitRange] < min_error && mean_errorTab[i-lowFitRange] > (median_mean_error * 0.1))
			{
				min_error = mean_errorTab[i-lowFitRange];
                                best_fit = i;
			}
			if(i >= 80 && i < 90 && best_fit == 0 && mean_errorTab[i-lowFitRange] < min_error && mean_errorTab[i-lowFitRange] > (median_mean_error * 0.1))
			{
                                min_error = mean_errorTab[i-lowFitRange];
                                best_fit = i;
                        }
			if(i >= 70 && i < 80 && best_fit == 0 && mean_errorTab[i-lowFitRange] < min_error && mean_errorTab[i-lowFitRange] > (median_mean_error * 0.1))
                        {
                                min_error = mean_errorTab[i-lowFitRange];
                                best_fit = i;
                        }
			if(i >= 60 && i < 70 && best_fit == 0 && mean_errorTab[i-lowFitRange] < min_error && mean_errorTab[i-lowFitRange] > (median_mean_error * 0.1))
                        {
                                min_error = mean_errorTab[i-lowFitRange];
                                best_fit = i;
                        }
		}
	}
	if(best_fit == 0) 
	{
		for(int i = highFitRange; i >= lowFitRange; i--)
        	{
			if(pvalue_tab[i-lowFitRange] > 0.0001)
	                {
		               	if(i >= 90 && mean_errorTab[i-lowFitRange] < min_error && mean_errorTab[i-lowFitRange] > (median_mean_error * 0.1)) 
	                        {
	                                min_error = mean_errorTab[i-lowFitRange];
	                                best_fit = i;
	                        }
	                        if(i >= 80 && i < 90 && best_fit == 0 && mean_errorTab[i-lowFitRange] < min_error && mean_errorTab[i-lowFitRange] > (median_mean_error * 0.1))
	                        {
	                                min_error = mean_errorTab[i-lowFitRange];
	                                best_fit = i;
	                        }
	                        if(i >= 70 && i < 80 && best_fit == 0 && mean_errorTab[i-lowFitRange] < min_error && mean_errorTab[i-lowFitRange] > (median_mean_error * 0.1))
	                        {
	                                min_error = mean_errorTab[i-lowFitRange];
	                                best_fit = i;
	                        }
	                        if(i >= 60 && i < 70 && best_fit == 0 && mean_errorTab[i-lowFitRange] < min_error && mean_errorTab[i-lowFitRange] > (median_mean_error * 0.1))
	                        {
	                                min_error = mean_errorTab[i-lowFitRange];
	                                best_fit = i;
	                        } 
			}
	
		}
	}

	if(best_fit == 0) 
	{
	
		for(int i = highFitRange; i >= lowFitRange; i--)
                {
                        if(pvalue_tab[i-lowFitRange] > average_pvalue)
                        {
	                        if(i >= 90 && mean_errorTab[i-lowFitRange] < min_error && mean_errorTab[i-lowFitRange] > (median_mean_error * 0.1)) 
	                        {
	                                min_error = mean_errorTab[i-lowFitRange];
	                                best_fit = i;
	                        }
	                        if(i >= 80 && i < 90 && best_fit == 0 && mean_errorTab[i-lowFitRange] < min_error && mean_errorTab[i-lowFitRange] > (median_mean_error * 0.1))
	                        {
	                                min_error = mean_errorTab[i-lowFitRange];
	                                best_fit = i;
	                        }
	                        if(i >= 70 && i < 80 && best_fit == 0 && mean_errorTab[i-lowFitRange] < min_error && mean_errorTab[i-lowFitRange] > (median_mean_error * 0.1))
	                        {
	                                min_error = mean_errorTab[i-lowFitRange];
	                                best_fit = i;
	                        }
	                        if(i >= 60 && i < 70 && best_fit == 0 && mean_errorTab[i-lowFitRange] < min_error && mean_errorTab[i-lowFitRange] > (median_mean_error * 0.1))
	                        {
	                                min_error = mean_errorTab[i-lowFitRange];
	                                best_fit = i;
	                        }
			}
     
                }	
	}

	if(best_fit == 0)
        {
		cout<<endl<<"Unable to find a good fit"<<endl;
                best_fit = 200;
	}


	cout<<endl<<"Best fit range = "<<best_fit<<" %"<<endl;
	cout<<endl<<"mean_errorTab[best_fit] = "<<mean_errorTab[best_fit-lowFitRange]<<", pvalue_tab[best_fit] = "<<pvalue_tab[best_fit-lowFitRange]<<endl;
	//FIXME v (Results_v6)
	//if(eta == "Barrel" && r9 == "low" && dataType == "MC" && fitVariable == "mmg_s" && fitFunction == "voigtian") best_fit = 78;
	//if(eta == "Barrel" && r9 == "all" && dataType == "MC" && fitVariable == "mmg_s" && fitFunction == "voigtian") best_fit = 76;

	cout<<endl<<"Best fit range = "<<best_fit<<" %"<<endl;
        cout<<endl<<"mean_errorTab[best_fit] = "<<mean_errorTab[best_fit-lowFitRange]<<", pvalue_tab[best_fit] = "<<pvalue_tab[best_fit-lowFitRange]<<endl;	

	///// Fit range systematics /////

	double systematics = 0;
	double temp_error = 0;

	//int lowBoundary = ((best_fit - 10) > lowFitRange) ? best_fit - 10 : lowFitRange;
	//int highBoundary = ((best_fit + 10) < highFitRange) ? best_fit + 10 : highFitRange;

	int lowBoundary = best_fit - 10;
	int highBoundary = best_fit + 10;
	if(lowBoundary < lowFitRange)
	{
		lowBoundary = lowFitRange;
		highBoundary = lowFitRange + 20;
	}
	if(highBoundary > highFitRange)
        {   
                highBoundary = highFitRange;
                lowBoundary = highFitRange - 20;
        }	

	for(int i = lowBoundary; i <= highBoundary; i++)
	{
		if(mean_errorTab[i-lowFitRange] > fabs(10 * mean_errorTab[best_fit-lowFitRange]) || pvalue_tab[i-lowFitRange] < 0.000001)
		{
			continue;
		}	
		
		if(meanTab[i-lowFitRange] < meanTab[best_fit-lowFitRange])
		{
			temp_error = fabs(meanTab[best_fit-lowFitRange] - meanTab[i-lowFitRange]);
		}
		else
		{
			temp_error = fabs(meanTab[i-lowFitRange] - meanTab[best_fit-lowFitRange]);	
		}
		if(temp_error > systematics) systematics = temp_error;
		//cout<<endl<<"temp_error = "<<temp_error;
	}

	cout<<endl<<"mean = "<<meanTab[best_fit-lowFitRange]<<" +- "<<mean_errorTab[best_fit-lowFitRange]<<" (stat) +- "<<systematics<<" (range) %"<<endl;

	system(Form("mkdir -p %s/%s/Selected_Fits/",directoryName.c_str(),dataType.c_str()));

	//system(Form("source create_directory.sh %sSelected_Fits/",directoryName.c_str()));

	system(Form("cp -r %s/%s/%dPercents/%s_%s_%sR9_%s/ %s/%s/Selected_Fits/%s_%s_%sR9_%s/",directoryName.c_str(),dataType.c_str(),best_fit,fitVariable.c_str(),eta.c_str(),r9.c_str(),fitFunction.c_str(),directoryName.c_str(),dataType.c_str(),fitVariable.c_str(),eta.c_str(),r9.c_str(),fitFunction.c_str()));

	ofstream summaryFile(Form("%s/%s/Selected_Fits/Summary_%s.txt",directoryName.c_str(),dataType.c_str(),fitVariable.c_str()), ios::app);

        summaryFile << fitVariable << " " << dataType << " >> " << r9 << " r9 " << eta << ", " << fitFunction <<" : mu = " << 100 * meanTab[best_fit-lowFitRange] << " +- " << 100 * mean_errorTab[best_fit-lowFitRange] << " (stat) +- " << 100 * systematics << " (range) %" << ", fit range = "<< best_fit <<"%, p-value = " << pvalue_tab[best_fit-lowFitRange] << ", scale factor = " << 1.0 / ( meanTab[best_fit-lowFitRange] + 1 ) << endl;

        summaryFile.close();
	
	delete [] meanTab;
	meanTab = 0;
	delete [] mean_errorTab;
	mean_errorTab = 0;
	delete [] pvalue_tab;
	pvalue_tab = 0;

	return 0;
}
