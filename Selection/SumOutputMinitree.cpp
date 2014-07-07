#include <iostream>
#include <fstream>
#include <cstdlib>
#include "TFile.h"

using namespace std;



int main(int argc, char *argv[])
{
	for(int iarg = 0 ; iarg < argc; iarg++)
        {
                cout << "argv[" << iarg << "]= " << argv[iarg] << endl;
        }

        if( argc == 1 ) 
        {
                cerr << "arguments should be passed : logfile" <<endl; 
                return 1;

        }    

        string logfile = "";

	if( argc > 1 ) logfile = argv[1];
	
	string systemChain = Form("ls stored_output/%s*.out",logfile.c_str());
	systemChain += " >> temp_file";
	cout << endl << systemChain;
	cout << endl << "coucou";
	system(systemChain.c_str());

	string fileChain = "";
	string file2Chain = "";
	int file2_int = 0;

	fstream file;
	file.open("temp_file",ios::in);

	vector<int> nObjectsAfterCuts;
	vector<int> nEventsAfterCuts;

	
	while(! file.eof())
	{
		file >> fileChain;

		fstream file2;
		file2.open(fileChain.c_str(),ios::in);
		
		while(file2Chain != "TOTALnbMuonsAfterID[0]")
		{
			file2 >> file2Chain;
		}
		file2 >> file2_int;
		cout << endl << "file2_int = " << file2_int <<endl;
		while(file2Chain != "TOTALnbMuonsAfterID[0]")
                {
                        file2 >> file2Chain;
                }
                file2 >> file2_int;
		cout << endl << "file2_int = " << file2_int <<endl;

		file2.close();
	}  
	file.close();

	return 0;
}
