#ifndef _CORRECTIONS_H
#define _CORRECTIONS_H

// C++ headers
#include <sstream>
#include <iostream>
#include <fstream>
#include <utility>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

// Root headers
#include "TMath.h"

using namespace std;

void parameters_Ceta(vector<double> &param, string correctionSet);
void parameters_fbrem(vector<double> &param, string correctionSet, bool isEB);
void parameters_feteta(vector<double> &param, string correctionSet, bool isEB);
double fEta(vector<double> param, double eta); 
double BremCor(vector<double> param, double brem);
double EtEtaCor(vector<double> param, double et, double eta, bool isEB); 
float ETHZ_fBremEta(float sigmaPhiSigmaEta, float eta, int algorithm);
float ETHZ_fEt(float ET, int algorithm); 
float ETHZ_fEnergy(float E, int algorithm); 

#endif


