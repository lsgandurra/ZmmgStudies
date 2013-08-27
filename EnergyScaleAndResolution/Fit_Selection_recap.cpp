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
	int nbRows = 0;

	string fileName = "";

	double average_pvalue = 0;
	double average_mean_error = 0;
	double median_mean_error = 0;
	double temp_var = 0;
	int best_fit = 0;
	double min_error = 100;
	
	int fit_number = highFitRange - lowFitRange + 1;	

	double * meanTab = new double[fit_number];
	double * mean_errorTab = new double[fit_number];
	double * pvalue_tab = new double[fit_number];

	vector <double> fitParameters;
	vector <double> mean_errorVec;

	for(int i = lowFitRange; i <= highFitRange; i++)
	{
		fileName = Form("%s/%s/%dPercents/%s_%s_%sR9_%s/fitsInformationsRaw.txt",directoryName.c_str(),dataType.c_str(),i,fitVariable.c_str(),eta.c_str(),r9.c_str(),fitFunction.c_str());		

		nbRows = rowsNumberInFile(fileName);

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

		pvalue_tab[i-lowFitRange] = fitParameters[fitParameters[2] * 2 + 5]; 
		average_pvalue += pvalue_tab[i-lowFitRange];

		meanTab[i-lowFitRange] = fitParameters[3];
		mean_errorTab[i-lowFitRange] = fitParameters[4];
		mean_errorVec.push_back(fitParameters[4]);	
	
		fitParameters.erase(fitParameters.begin(),fitParameters.end());
	
		fitParametersFile.close();

	}

	//mean_errorVec.sort();
	sort(mean_errorVec.begin(), mean_errorVec.end());
	median_mean_error = mean_errorVec[(highFitRange - lowFitRange) / 2];
	average_pvalue = average_pvalue / fit_number;
	//average_mean_error = average_mean_error / fit_number;
		
	min_error = median_mean_error * 1.5;

	cout<< "median_mean_error = "<<median_mean_error<<endl;
	cout<< "min_error = "<<min_error<<endl;
	cout<< "lowFitRange = "<<lowFitRange<<", highFitRange = "<<highFitRange<<endl;
	
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
