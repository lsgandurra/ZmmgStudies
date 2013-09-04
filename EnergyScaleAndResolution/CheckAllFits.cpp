//Check if all fit of s are ok, if not, produce a txt file with fits to resubmit.
#include "functions.h"

int main(int argc, char *argv[])
{
	for(int iarg = 0 ; iarg < argc; iarg++)
        {
                cout << "argv[" << iarg << "]= " << argv[iarg] << endl;
        }

        if( argc == 1 ) 
        {
                cerr << "arguments should be passed : directoryName, dataType, fitVariable, lowFitRange, highFitRange" <<endl; 
                return 1;

        }    

        string directoryName = "Results_v8_surface";
        string dataType = "MC";
        string fitVariable = "mmg_s_MZ_Surface";
	string cutVariable = "Photon_Et";
	string binFileName = "LimitesAllPt.txt";
        int lowFitRange = 60;
	int highFitRange = 100;
	string injectedResolution = "0";

	if( argc > 1 ) directoryName = argv[1];
	if( argc > 2 ) dataType = argv[2];
	if( argc > 3 ) fitVariable = argv[3];
	if( argc > 4 ) cutVariable = argv[4];
	if( argc > 5 ) binFileName = argv[5];
	if( argc > 6 ) 
	{
		std::stringstream ss ( argv[6] );
		ss >> lowFitRange;
	}
	if( argc > 7 ) 
	{
		std::stringstream ss ( argv[7] );
		ss >> highFitRange;
	}
	if( argc > 8 ) injectedResolution = argv[8];

	string jobsToResubmitString = Form("jobsToResubmit_%s_%s_%s.sh",directoryName.c_str(),dataType.c_str(),fitVariable.c_str());	

	if( FileExists( jobsToResubmitString.c_str() ) )
	{
		system(Form("rm %s",jobsToResubmitString.c_str()));
	}

	int nBins = rowsNumberInFile(binFileName.c_str()) - 1;

	string eta = "";
	string r9 = "";
	string fitFunction = "";
	string fileName = "";
	string resub = "";
	int good = 1;

	for(int i = lowFitRange; i <= highFitRange; i++)
	{
		good = 1;
		for(int k = 0; k < 12; k++) //FIXME if you want more fonctions
		{ 

			if(k == 0)
			{
				eta = "Barrel";
				r9 = "low";
				fitFunction = "voigtian";
			}				
			if(k == 1)
                        {
        	                eta = "Barrel";
                                r9 = "low";
                                fitFunction = "cruijff";
                        }
			if(k == 2)
                        {
                                eta = "Barrel";
                                r9 = "high";
                                fitFunction = "voigtian";
                        }
			if(k == 3)
                        {
                                eta = "Barrel";
                                r9 = "high";
                                fitFunction = "cruijff";
                        }
			if(k == 4)
                        {
                                eta = "Barrel";
                                r9 = "all";
                                fitFunction = "voigtian";
                        }
			if(k == 5)
                        {
                                eta = "Barrel";
                                r9 = "all";
                                fitFunction = "cruijff";
                        }
			if(k == 6)
                        {
                                eta = "Endcaps";
                                r9 = "low";
                                fitFunction = "cruijff";
                        }
			if(k == 7)
                        {
                                eta = "Endcaps";
                                r9 = "low";
                                fitFunction = "cruijff";
                        }
			if(k == 8)
                        {
                                eta = "Endcaps";
                                r9 = "high";
                                fitFunction = "voigtian";
                        }
			if(k == 9)
                        {
                                eta = "Endcaps";
                                r9 = "high";
                                fitFunction = "cruijff";
                        }
			if(k == 10)
                        {
                                eta = "Endcaps";
                                r9 = "all";
                                fitFunction = "voigtian";
                        }
			if(k == 11)
                        {
                                eta = "Endcaps";
                                r9 = "all";
                                fitFunction = "cruijff";
                        }

			for(int j = 1; j <= nBins; j++)
                	{
				fileName = Form("%s/InjectedResolution_%sPercent/%s/%dPercents/%s_%s_%sR9_%s/fit_%d.png",directoryName.c_str(),injectedResolution.c_str(),dataType.c_str(),i,fitVariable.c_str(),eta.c_str(),r9.c_str(),fitFunction.c_str(),j);
			
				if( FileExists( fileName.c_str() ) ) cout << endl << "OK >> " << fileName << " exists.";
				else
				{	
					cout << endl << "!! Warning, need to resubmit >> " << fileName;
					good = 0;
				}
			}
		}

		if(good == 0)
		{
			ofstream filesToResubmit(jobsToResubmitString.c_str(), ios::app);
			resub = Form("qsub batchJob.sh input %s %s %s %s %s %d %s",directoryName.c_str(),dataType.c_str(),fitVariable.c_str(),cutVariable.c_str(),binFileName.c_str(),i,injectedResolution.c_str());
			filesToResubmit << resub << endl;
			filesToResubmit.close();
		}
	}


	return 0;

}
