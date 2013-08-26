#include "TMVA/Reader.h"
#include "TMVA/Tools.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TBranch.h"
#include "TChain.h"
#include "TFile.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TString.h"
#include "TBits.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TTimeStamp.h"
#include "TSystem.h"
#include "TProof.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <utility>
#include <algorithm>

//#pragma optimize 0

#include "interface/TRootBardak.h"
#include "interface/TRootBeamSpot.h"
#include "interface/TRootCluster.h"
#include "interface/TRootDummyEvent.h"
#include "interface/TRootEcalRecHit.h"
#include "interface/TRootElectron.h"
#include "interface/TRootEvent.h"
#include "interface/TRootJet.h"
#include "interface/TRootMCParticle.h"
#include "interface/TRootMCPhoton.h"
#include "interface/TRootMET.h"
#include "interface/TRootMuon.h"
#include "interface/TRootParticle.h"
#include "interface/TRootPhoton.h"
#include "interface/TRootRun.h"
#include "interface/TRootSignalEvent.h"
#include "interface/TRootSuperCluster.h"
#include "interface/TRootTopTop.h"
#include "interface/TRootTrack.h"
#include "interface/TRootVertex.h"


using namespace std;



int factorial(int number) {
	int temp;

	if(number <= 1) return 1;

	temp = number * factorial(number - 1);
	return temp;
}


double DeltaR( double eta1, double phi1, double eta2, double phi2)
{
	double DeltaEta = fabs( eta1-eta2 );
	double DeltaPhi = fabs( phi1-phi2 );
	// Returning DeltaPhi in the correct range (0, 2pi)
	while (DeltaPhi >	 TMath::Pi()) DeltaPhi -= 2*TMath::Pi();
	while (DeltaPhi <= -TMath::Pi()) DeltaPhi += 2*TMath::Pi();
	return sqrt(DeltaEta*DeltaEta + DeltaPhi*DeltaPhi);
}

void doGenInfo(TRootParticle* myparticle, TClonesArray* mcParticles, float* particule_trueE, float* particule_truePx, float* particule_truePy, float* particule_truePz, float* particule_trueEta, float* particule_truePhi, int particle_pdgId = 0)
{
  TRootMCParticle* mygenparticle;
  int NbMCpartInCone=0;
  double bestPtdiff=500.0;
  int igpsl=-1;
  for (int igp=0; igp<mcParticles->GetEntriesFast(); igp++) {
    mygenparticle = (TRootMCParticle*) mcParticles->At(igp);
    //cout<<endl<<"deltaR = "<<DeltaR(mygenparticle->Eta(), mygenparticle->Phi(), myparticle->Eta(), myparticle->Phi())<<endl;
    //cout<<"mygenparticle->Mag() = "<<mygenparticle->Mag()<<endl;
    if (DeltaR(mygenparticle->Eta(), mygenparticle->Phi(), myparticle->Eta(), myparticle->Phi())<0.3){
      if ( (mygenparticle->status()==1) && ( (particle_pdgId==0)?true:((mygenparticle->type())==particle_pdgId) ) ){
        //HistoMCpartStatus1InConeId->Fill(mygenparticle->type());
        NbMCpartInCone++;
        if (fabs( (mygenparticle->Pt()) - (myparticle->Pt()) )<bestPtdiff){
          bestPtdiff=fabs(mygenparticle->Pt()-myparticle->Pt());
          igpsl=igp;
        }
      }
    }
  }
  if (igpsl!=-1){
    mygenparticle = (TRootMCParticle*) mcParticles->At(igpsl);
    *particule_trueE = mygenparticle->Energy();
    *particule_truePx = mygenparticle->Px();
    *particule_truePy = mygenparticle->Py();
    *particule_truePz = mygenparticle->Pz();
    *particule_truePhi = mygenparticle->Phi();
    *particule_trueEta = mygenparticle->Eta();

  }
  return;
}

double applySidra( double _pt, double charge, double eta, double phi, TRandom3* generator)
{
	double pt = _pt;
// Correct MC
	double a = 0.0650687e-3;
	double b = 0.212987e-3;
	double c = 1.53414;
	pt = (double)1.0/(double)pt;
	pt -= a - b * charge * sin( phi + c );
// Apply Corrections
	double A = 0.143812;
	double B = 0.0404834;
	double Ap = 0.0995898;
	double Bp = -0.0692569;
	double Cp = 0.0952174;
	double phi0 = -1.08881;
	pt += ( A + B * eta * eta ) * (generator->Gaus(0,1)) /1000. + (Ap + Cp * charge * sin( phi + phi0 )+ Bp * charge * eta )/1000.;
	pt = (double)1.0/(double)pt;
	return pt;
}

double applyMuScleFit(double _pt, double charge, double eta, double phi)
{
	double b = -5.03313e-6;
	double c = -4.41463e-5;
	double d0 = -0.000148871;
	double e0 = 1.59501;
	double d1 = 7.95495e-5;
	double e1 = -0.364823;
	double d2 = 0.000152032;
	double e2 = 0.410195;
	double d = eta > .9 ? d1 : (eta < -.9 ? d2 : d0);
	double e = eta > .9 ? e1 : (eta < -.9 ? e2 : e0);
	double pt = _pt;
	double sgn_eta = eta >= 0.0 ? 1.0 : -1.0;
	pt = pt * (1 + b * pt  + c * charge * pt * sgn_eta * eta * eta + charge * d * pt * sin( phi + e ));
	return pt;

}

TLorentzVector getCorrectedMuon();

int main(int argc, char *argv[]);




