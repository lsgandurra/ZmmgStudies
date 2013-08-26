#include <fstream>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstring>
#include <iomanip>
#include <TMath.h>

#pragma optimize 0
using namespace std;

double GlobalPvalue (string data_file_adress)
{

		//1 Initialization

	int gamma_categories_int = 0, parameters_int =0, i = 0 /*incrementer*/;
	double chi_square = 0.0;

		//2 Notice

/*
The data file which contains the values of theta_true, sigma_true, theta_fit and sigma_fit must follow this patern:

	Number of different categories of photons (ndcp)
	Number of parameters
	theta_true[0]
	...
	theta_true[ndcp-1]
	sigma_true[0]
	...
	sigma_true[ndcp-1]
	theta_fit[0]
	...
	theta_fit[ndcp-1]
	sigma_fit[0]
	...
	sigma_fit[ndcp-1]
*/

		//3 Reading of the number of different categories of photons and the number of parameters

	ifstream myStream (data_file_adress.c_str());

	if (myStream)
	{

		myStream >> gamma_categories_int;
		myStream >> parameters_int;
	}

	else
        {

        cout << "ERROR: Cannot open the file in reading mode" << endl;

        }

		//4 Reading of theta_true, sigma_true, theat_fit and sigma_fit

	double *theta_true = new double [gamma_categories_int];
        double *sigma_true = new double [gamma_categories_int];
        double *theta_fit = new double [gamma_categories_int];
        double *sigma_fit = new double [gamma_categories_int];

	for (i = 0 ; i <= 4*gamma_categories_int-1 ; i++)
	{
		if (i <= gamma_categories_int-1)
		{
			myStream >> theta_true[i];
		}
		
		else
		{
			if (i >= gamma_categories_int && i <= 2*gamma_categories_int-1)
			{
				myStream >> sigma_true[i-gamma_categories_int];
			}

			else
			{
				if (i >= 2*gamma_categories_int && i <= 3*gamma_categories_int-1)
				{
					myStream >> theta_fit[i-2*gamma_categories_int];
				}

				else
				{
					myStream >> sigma_fit[i-3*gamma_categories_int];
				}
			}
		}
	}

	

		//5 Chi_square test & p-value

	for (i = 0; i <= gamma_categories_int-1; i++)
	{
		chi_square += (theta_fit[i] - theta_true[i])*(theta_fit[i] - theta_true[i])/(sigma_fit[i]*sigma_fit[i] + sigma_true[i]*sigma_true[i]);
	}

	delete [] theta_true; delete [] sigma_true; delete [] theta_fit; delete [] sigma_fit;
	theta_true = 0; sigma_true = 0; theta_fit = 0; sigma_fit = 0;
//	cout<<endl<<"chi_square = "<<chi_square<<endl<<"gamma_categories_int = "<<gamma_categories_int<<endl<<"parameters_int = "<<parameters_int<<endl;

	return TMath::Prob(chi_square, gamma_categories_int - parameters_int);

}


int main(int argc, char *argv[])
{
	for(int iarg = 0 ; iarg < argc; iarg++)
        {
                cout << "argv[" << iarg << "]= " << argv[iarg] << endl;
        }

        if( argc == 1 ) 
        {
                cerr << "arguments should be passed : txtFile" <<endl; 
                return 1;

        }    

        string txtFile = "Global.txt";

        if( argc > 1 ) txtFile = argv[1];


	cout<<endl<<endl<<"GlobalPvalue = "<<GlobalPvalue(txtFile)<<endl<<endl;

	return 0;

}

