#include <TLegend.h>
#include <TStyle.h>
#include <TFormula.h>
#include <TF1.h>
#include <TSystem.h>
#include <TClonesArray.h>
#include <TLeaf.h>
#include <TChain.h>
#include <TObject.h>
#include <string.h>
#include <algorithm>
#include <TBranch.h>
#include <TString.h>
#include <TBits.h>
#include <TMath.h>
#include "TROOT.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <TLatex.h>

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


#pragma optimize 0

using namespace std;

void DrawDataDataMCMCplot(TTree *Data_miniTree, TTree *DYToMuMu_miniTree, TTree *OLD_Data_miniTree, TTree *OLD_DYToMuMu_miniTree, string var, string pic, string limits, string cut, string name, string Title, bool inlog, bool drawUnderOverFsubleading, TCanvas *c1);

//void DrawDataMCplot(TTree *Data_miniTree, TTree *FSR_DYToMuMu_miniTree, TTree *nonFSR_DYToMuMu_miniTree, TTree *QCDMu_miniTree, TTree* TTJets_miniTree, TTree *WJetsToLNu_miniTree, string var, string pic, string limits, string cut, string name, string Title, bool inlog, bool drawUnderOverFsubleading, TCanvas *c1, bool doFit = false);
//void DrawDataMCplot(TTree *Data_miniTree, TTree *FSR_DYToMuMu_miniTree, TTree *nonFSR_DYToMuMu_miniTree, TTree *TTJets_miniTree, TTree *WJetsToLNu_miniTree, string var, string pic, string limits, string cut, string name, string Title, bool inlog, bool drawUnderOverFsubleading, TCanvas *c1, bool doFit = false);

//void DrawDataMCplot_TH1I(TTree *Data_miniTree, TTree *FSR_DYToMuMu_miniTree, TTree *nonFSR_DYToMuMu_miniTree, TTree *TTJets_miniTree, TTree *WJetsToLNu_miniTree, string var, string pic, string limits, string cut, string name, string Title, bool inlog, bool drawUnderOverFsubleading, TCanvas *c1, bool doFit = false);

//void DrawDataMCplot_TH2F(TTree *Data_miniTree, TTree *FSR_DYToMuMu_miniTree, TTree *nonFSR_DYToMuMu_miniTree, TTree *TTJets_miniTree, TTree *WJetsToLNu_miniTree, string var1, string var2, string limits, string cut, string name, string Title_var1, string Title_var2, string pic, bool inlog_var1, bool inlog_var2, TCanvas *c1, bool doFit = false);

